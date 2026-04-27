/*

  VSEARCH: a versatile open source tool for metagenomics

  Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
  All rights reserved.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway

  This software is dual-licensed and available under a choice
  of one of two licenses, either under the terms of the GNU
  General Public License version 3 or the BSD 2-Clause License.

  GNU General Public License version 3 / BSD 2-Clause License:
  See LICENSE_GNU_GPL3.txt or LICENSE.txt.

*/

/*
  Daemon and submit-query commands.

  Daemon design
  =============
  - Reference database is loaded ONCE at startup via search_load_db().
  - The watch directory is monitored with kqueue (macOS) via dirwatch.
  - For each arriving FASTA file (or batch of files):
      1. A timestamp is generated (YYYYMMDDTHHMMSS).
      2. Queued files are grouped into sub-batches by max_batch_sequences.
      3. Single-file batches: existing fork/search/move path (unchanged).
      4. Multi-file batches:
           a. All files claimed atomically to .processing staging paths.
           b. Sequences merged into one FASTA with VSDF#### tags in headers.
           c. One child process searches the merged file.
           d. Parent splits output files by tag into per-file outputs.
           e. Staging files moved to output_dir or errors_dir.

  Batch merging rationale
  =======================
  With many small query files queued, thread workers would be underutilised
  because each individual search runs too few sequences to saturate all
  threads.  By merging, a single search saturates all opt_threads workers
  continuously. Fixed per-search overhead (fork, output-file open/close,
  thread setup) is also amortized across the batch.

  VSDF tag format
  ===============
  Each query header in the merged file is prefixed with "VSDF####_" where
  #### is a zero-padded 4-digit file index (0000–9999).  After searching,
  output files are scanned for this prefix to route each result line to the
  correct per-file output, and the prefix is stripped before writing.

  Why fork()?
  ===========
  vsearch's fatal() calls std::exit().  Without fork(), a single bad
  query file would kill the entire daemon and lose the loaded database.
  The child inherits the in-memory database via copy-on-write; since the
  search workers only read the DB, no physical copies are made.

  Atomic rename (submit_query)
  ============================
  Writers call --submit_query to drop files into the watch directory.
  The file is first written to a .tmp path, then renamed to the final
  .fasta path.  The daemon's dirwatch skips .tmp files, so it only ever
  sees complete files.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "vsearch.h"
#include "config_yaml.h"
#include "daemon.h"
#include "dirwatch.h"
#include "search.h"
#include "utils/fatal.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

#include <fcntl.h>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>


/* ------------------------------------------------------------------ */
/* Shutdown flag and per-daemon state                                   */
/* ------------------------------------------------------------------ */

static volatile sig_atomic_t g_shutdown = 0;

/* cmdline is stored so forked children can pass it to
   search_open_output_files() without extra arguments. */
static char * g_cmdline = nullptr;

static auto handle_signal(int /*sig*/) -> void {
  g_shutdown = 1;
}


/* ------------------------------------------------------------------ */
/* Utility: generate ISO 8601 compact timestamp                         */
/* ------------------------------------------------------------------ */

static auto make_timestamp(char * buf, std::size_t bufsize) -> void {
  std::time_t now = std::time(nullptr);
  struct tm tm_info;
  localtime_r(&now, &tm_info);
  std::strftime(buf, bufsize, "%Y%m%dT%H%M%S", &tm_info);
}


/* ------------------------------------------------------------------ */
/* Utility: extract the stem (filename without last extension)          */
/* ------------------------------------------------------------------ */

static auto stem_of(const std::string & path) -> std::string {
  /* basename part */
  auto slash = path.rfind('/');
  std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
  /* strip last extension */
  auto dot = name.rfind('.');
  if (dot != std::string::npos) {
    name = name.substr(0, dot);
  }
  return name;
}


/* ------------------------------------------------------------------ */
/* Utility: atomic file copy then rename                                */
/* ------------------------------------------------------------------ */

static auto copy_file_atomic(const char * src,
                              const char * dst_final) -> bool
{
  /* Write to a .tmp alongside dst_final, then rename. */
  std::string tmp = std::string(dst_final) + ".tmp";

  int in_fd = open(src, O_RDONLY);
  if (in_fd < 0) { return false; }

  int out_fd = open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (out_fd < 0) { close(in_fd); return false; }

  char buf[65536];
  ssize_t n = 0;
  bool ok = true;
  while ((n = read(in_fd, buf, sizeof(buf))) > 0) {
    if (write(out_fd, buf, static_cast<std::size_t>(n)) != n) {
      ok = false;
      break;
    }
  }
  if (n < 0) { ok = false; }

  close(in_fd);
  close(out_fd);

  if (!ok) {
    unlink(tmp.c_str());
    return false;
  }

  if (rename(tmp.c_str(), dst_final) != 0) {
    unlink(tmp.c_str());
    return false;
  }
  return true;
}


/* ------------------------------------------------------------------ */
/* Utility: mkdir if it doesn't exist                                   */
/* ------------------------------------------------------------------ */

static auto ensure_dir(const std::string & path) -> void {
  struct stat st;
  if (stat(path.c_str(), &st) == 0) {
    if (S_ISDIR(st.st_mode)) { return; }
    fatal("Path exists but is not a directory", path.c_str());
  }
  if (mkdir(path.c_str(), 0755) != 0 && errno != EEXIST) {
    fatal("Unable to create directory", path.c_str());
  }
}


/* ------------------------------------------------------------------ */
/* Utility: set output file opt_* globals for one query file            */
/* ------------------------------------------------------------------ */

static auto set_output_paths(const daemon_config_t & cfg,
                              const std::string & stem_ts,
                              std::vector<std::string> & allocated_strings) -> void
{
  /* Helper lambda: allocate a path string and point the opt_* global at it. */
  auto set_opt = [&](char *& opt, const char * ext) {
    std::string p = cfg.output_dir + "/" + stem_ts + ext;
    char * s = static_cast<char *>(std::malloc(p.size() + 1));
    if (s == nullptr) { fatal("daemon: out of memory"); }
    std::memcpy(s, p.c_str(), p.size() + 1);
    allocated_strings.push_back(p);  /* keep track for potential cleanup */
    opt = s;
  };

  if (cfg.enable_alnout)            { set_opt(opt_alnout,            ".alnout");   }
  if (cfg.enable_blast6out)         { set_opt(opt_blast6out,         ".blast6out"); }
  if (cfg.enable_uc)                { set_opt(opt_uc,                ".uc");       }
  if (cfg.enable_userout)           { set_opt(opt_userout,           ".userout");  }
  if (cfg.enable_samout)            { set_opt(opt_samout,            ".sam");      }
  if (cfg.enable_matched)           { set_opt(opt_matched,           "_matched.fasta"); }
  if (cfg.enable_notmatched)        { set_opt(opt_notmatched,        "_notmatched.fasta"); }
  if (cfg.enable_fastapairs)        { set_opt(opt_fastapairs,        ".fastapairs"); }
  if (cfg.enable_lcaout)            { set_opt(opt_lcaout,            ".lcaout");   }
  if (cfg.enable_qsegout)           { set_opt(opt_qsegout,           ".qsegout");  }
  if (cfg.enable_tsegout)           { set_opt(opt_tsegout,           ".tsegout");  }
  if (cfg.enable_otutabout)         { set_opt(opt_otutabout,         ".otutab");   }
  if (cfg.enable_mothur_shared_out) { set_opt(opt_mothur_shared_out, ".shared");   }
  if (cfg.enable_biomout)           { set_opt(opt_biomout,           ".biom");     }
}


/* ------------------------------------------------------------------ */
/* Utility: clear output opt_* globals (reset between files in child)  */
/* ------------------------------------------------------------------ */

static auto clear_output_opts() -> void {
  opt_alnout = nullptr;
  opt_blast6out = nullptr;
  opt_uc = nullptr;
  opt_userout = nullptr;
  opt_samout = nullptr;
  opt_matched = nullptr;
  opt_notmatched = nullptr;
  opt_fastapairs = nullptr;
  opt_lcaout = nullptr;
  opt_qsegout = nullptr;
  opt_tsegout = nullptr;
  opt_otutabout = nullptr;
  opt_mothur_shared_out = nullptr;
  opt_biomout = nullptr;
}


/* ------------------------------------------------------------------ */
/* Process one file: run in a forked child                             */
/* ------------------------------------------------------------------ */

/* Returns the child exit status (0 = success). */
static auto process_file_forked(const daemon_config_t & cfg,
                                 const std::string & input_path,
                                 const std::string & stem_ts) -> int
{
  pid_t pid = fork();
  if (pid < 0) {
    std::fprintf(stderr, "daemon: fork() failed: %s\n", strerror(errno));
    return 1;
  }

  if (pid == 0) {
    /* ---- child ---- */
    clear_output_opts();

    std::vector<std::string> alloc;
    set_output_paths(cfg, stem_ts, alloc);

    /* Open per-file output file handles, run search, close them. */
    search_open_output_files(g_cmdline, g_cmdline);
    int rc = search_process_query_file(input_path.c_str());
    search_close_output_files();
    std::exit(rc);
  }

  /* ---- parent: wait for child ---- */
  int status = 0;
  while (waitpid(pid, &status, 0) < 0) {
    if (errno != EINTR) {
      std::fprintf(stderr, "daemon: waitpid() error: %s\n", strerror(errno));
      return 1;
    }
  }

  if (WIFEXITED(status)) {
    return WEXITSTATUS(status);
  }
  /* Child killed by signal */
  return 1;
}


/* ================================================================== */
/* Batch merging support                                               */
/* ================================================================== */

/*
  VSDF tag: "VSDF####_" prepended to each query label during merging.
  #### = zero-padded 4-digit file index within the batch (0000–9999).
  Total tag length: 9 characters.
*/
static constexpr std::size_t VSDF_IDX_DIGITS = 4;
static constexpr std::size_t VSDF_TAG_LEN    = 4 + VSDF_IDX_DIGITS + 1; /* "VSDF####_" */


/* ------------------------------------------------------------------ */
/* Count sequences in a FASTA file (fast buffered scan for '>')        */
/* ------------------------------------------------------------------ */

static auto count_fasta_sequences(const std::string & path) -> std::size_t {
  std::FILE * fp = std::fopen(path.c_str(), "r");
  if (fp == nullptr) { return 0; }
  std::size_t count = 0;
  char buf[65536];
  std::size_t n = 0;
  bool at_sol = true;   /* at start of line */
  while ((n = std::fread(buf, 1, sizeof(buf), fp)) > 0) {
    for (std::size_t i = 0; i < n; ++i) {
      if (at_sol && buf[i] == '>') { ++count; }
      at_sol = (buf[i] == '\n');
    }
  }
  std::fclose(fp);
  return count;
}


/* ------------------------------------------------------------------ */
/* Merge multiple FASTA files into one with VSDF-tagged headers        */
/* ------------------------------------------------------------------ */

static auto merge_fasta_files(const std::vector<std::string> & paths,
                               const std::string & merged_path) -> bool {
  std::FILE * out = std::fopen(merged_path.c_str(), "w");
  if (out == nullptr) { return false; }

  char buf[65536];
  bool ok = true;
  for (std::size_t idx = 0; idx < paths.size() && ok; ++idx) {
    std::FILE * fp = std::fopen(paths[idx].c_str(), "r");
    if (fp == nullptr) { ok = false; break; }
    while (std::fgets(buf, static_cast<int>(sizeof(buf)), fp) != nullptr) {
      if (buf[0] == '>') {
        /* Rewrite ">label..." as ">VSDF####_label..." */
        if (std::fprintf(out, ">VSDF%0*zu_%s",
                         static_cast<int>(VSDF_IDX_DIGITS), idx, buf + 1) < 0) {
          ok = false; break;
        }
      } else {
        if (std::fputs(buf, out) == EOF) { ok = false; break; }
      }
    }
    std::fclose(fp);
  }
  std::fclose(out);
  if (!ok) { std::remove(merged_path.c_str()); }
  return ok;
}


/* ------------------------------------------------------------------ */
/* Parse VSDF file index from a label string                           */
/* ------------------------------------------------------------------ */

/* Returns the 0-based file index, or SIZE_MAX on parse failure or
   if idx >= n_files. */
static auto parse_vsdf_idx(const char * label,
                            std::size_t n_files) -> std::size_t {
  if (std::strncmp(label, "VSDF", 4) != 0) { return SIZE_MAX; }
  char tmp[VSDF_IDX_DIGITS + 1] = {};
  std::memcpy(tmp, label + 4, VSDF_IDX_DIGITS);
  char * end = nullptr;
  std::size_t idx = static_cast<std::size_t>(std::strtoul(tmp, &end, 10));
  if (end != tmp + VSDF_IDX_DIGITS) { return SIZE_MAX; }
  if (idx >= n_files) { return SIZE_MAX; }
  return idx;
}


/* ------------------------------------------------------------------ */
/* Split a tab-delimited output file (blast6out / uc) by VSDF tag      */
/* ------------------------------------------------------------------ */

/* query_field: 0-based index of the tab-delimited field holding the
   query label.  blast6out = 0, uc = 8. */
static auto split_tabular_output(const std::string & merged_path,
                                  const std::vector<std::string> & per_file_paths,
                                  int query_field) -> void {
  std::FILE * in = std::fopen(merged_path.c_str(), "r");
  if (in == nullptr) { return; }   /* no output file = nothing to split */

  std::size_t n_files = per_file_paths.size();
  std::vector<std::FILE *> fps(n_files, nullptr);
  char line[65536];

  while (std::fgets(line, static_cast<int>(sizeof(line)), in) != nullptr) {
    /* Strip trailing CR/LF */
    std::size_t llen = std::strlen(line);
    while (llen > 0 && (line[llen - 1] == '\n' || line[llen - 1] == '\r')) {
      line[--llen] = '\0';
    }
    if (llen == 0) { continue; }

    /* Walk to the query field by counting tab separators */
    char * field_start = line;
    bool ok = true;
    for (int f = 0; f < query_field; ++f) {
      char * tab = std::strchr(field_start, '\t');
      if (tab == nullptr) { ok = false; break; }
      field_start = tab + 1;
    }
    if (!ok) { continue; }

    /* field_start now points to the query label */
    char * field_end = std::strchr(field_start, '\t');  /* nullptr if last field */

    std::size_t idx = parse_vsdf_idx(field_start, n_files);
    if (idx == SIZE_MAX) { continue; }

    if (fps[idx] == nullptr) {
      fps[idx] = std::fopen(per_file_paths[idx].c_str(), "a");
      if (fps[idx] == nullptr) { continue; }
    }

    /* Write: everything before the query field */
    if (field_start > line) {
      std::fwrite(line, 1, static_cast<std::size_t>(field_start - line), fps[idx]);
    }
    /* Write: original label (skip "VSDF####_" prefix) */
    const char * orig = field_start + VSDF_TAG_LEN;
    std::size_t orig_len = (field_end != nullptr)
                           ? static_cast<std::size_t>(field_end - orig)
                           : std::strlen(orig);
    std::fwrite(orig, 1, orig_len, fps[idx]);
    /* Write: tab and remaining fields (if any), then newline */
    if (field_end != nullptr) {
      std::fputs(field_end, fps[idx]);
    }
    std::fputc('\n', fps[idx]);
  }

  for (auto * fp : fps) { if (fp != nullptr) { std::fclose(fp); } }
  std::fclose(in);
}


/* ------------------------------------------------------------------ */
/* Split a FASTA output file (matched / notmatched) by VSDF tag        */
/* ------------------------------------------------------------------ */

static auto split_fasta_output(const std::string & merged_path,
                                const std::vector<std::string> & per_file_paths) -> void {
  std::FILE * in = std::fopen(merged_path.c_str(), "r");
  if (in == nullptr) { return; }

  std::size_t n_files = per_file_paths.size();
  std::vector<std::FILE *> fps(n_files, nullptr);
  std::size_t cur = SIZE_MAX;
  char line[65536];

  while (std::fgets(line, static_cast<int>(sizeof(line)), in) != nullptr) {
    if (line[0] == '>') {
      cur = parse_vsdf_idx(line + 1, n_files);
      if (cur == SIZE_MAX) { continue; }
      if (fps[cur] == nullptr) {
        fps[cur] = std::fopen(per_file_paths[cur].c_str(), "a");
      }
      if (fps[cur] == nullptr) { cur = SIZE_MAX; continue; }
      /* Write header with original label (skip "VSDF####_") */
      std::fprintf(fps[cur], ">%s", line + 1 + VSDF_TAG_LEN);
    } else if (cur != SIZE_MAX && fps[cur] != nullptr) {
      std::fputs(line, fps[cur]);
    }
  }

  for (auto * fp : fps) { if (fp != nullptr) { std::fclose(fp); } }
  std::fclose(in);
}


/* ------------------------------------------------------------------ */
/* Split an alnout file (multi-line) by VSDF tag                       */
/* ------------------------------------------------------------------ */

/* alnout query blocks start with "Query >label  len=N".
   All subsequent lines belong to that query until the next "Query >" line. */
static auto split_alnout_output(const std::string & merged_path,
                                 const std::vector<std::string> & per_file_paths) -> void {
  std::FILE * in = std::fopen(merged_path.c_str(), "r");
  if (in == nullptr) { return; }

  std::size_t n_files = per_file_paths.size();
  std::vector<std::FILE *> fps(n_files, nullptr);
  std::size_t cur = SIZE_MAX;
  char line[65536];

  while (std::fgets(line, static_cast<int>(sizeof(line)), in) != nullptr) {
    if (std::strncmp(line, "Query >", 7) == 0) {
      /* "Query >VSDF####_original_label  len=N\n" */
      const char * label = line + 7;   /* points to "VSDF####_..." */
      cur = parse_vsdf_idx(label, n_files);
      if (cur == SIZE_MAX) { continue; }
      if (fps[cur] == nullptr) {
        fps[cur] = std::fopen(per_file_paths[cur].c_str(), "a");
      }
      if (fps[cur] == nullptr) { cur = SIZE_MAX; continue; }
      /* Write "Query >original_label  len=N\n" */
      std::fprintf(fps[cur], "Query >%s", label + VSDF_TAG_LEN);
    } else if (cur != SIZE_MAX && fps[cur] != nullptr) {
      std::fputs(line, fps[cur]);
    }
  }

  for (auto * fp : fps) { if (fp != nullptr) { std::fclose(fp); } }
  std::fclose(in);
}


/* ------------------------------------------------------------------ */
/* Split all enabled output formats after a batch search               */
/* ------------------------------------------------------------------ */

static auto split_batch_outputs(const daemon_config_t & cfg,
                                  const std::string & batch_stem,
                                  const std::vector<std::string> & stems_ts) -> void {
  std::size_t n = stems_ts.size();

  /* Helper: build per-file path list for a given extension. */
  auto per_file = [&](const char * ext) -> std::vector<std::string> {
    std::vector<std::string> paths(n);
    for (std::size_t i = 0; i < n; ++i) {
      paths[i] = cfg.output_dir + "/" + stems_ts[i] + ext;
    }
    return paths;
  };

  /* Helper: merged output path for a given extension. */
  auto merged = [&](const char * ext) -> std::string {
    return cfg.output_dir + "/" + batch_stem + ext;
  };

  if (cfg.enable_blast6out) {
    split_tabular_output(merged(".blast6out"), per_file(".blast6out"), 0);
    std::remove(merged(".blast6out").c_str());
  }
  if (cfg.enable_uc) {
    split_tabular_output(merged(".uc"), per_file(".uc"), 8);
    std::remove(merged(".uc").c_str());
  }
  if (cfg.enable_matched) {
    split_fasta_output(merged("_matched.fasta"), per_file("_matched.fasta"));
    std::remove(merged("_matched.fasta").c_str());
  }
  if (cfg.enable_notmatched) {
    split_fasta_output(merged("_notmatched.fasta"), per_file("_notmatched.fasta"));
    std::remove(merged("_notmatched.fasta").c_str());
  }
  if (cfg.enable_alnout) {
    split_alnout_output(merged(".alnout"), per_file(".alnout"));
    std::remove(merged(".alnout").c_str());
  }
  /* Formats without per-file splitting (userout, samout, lcaout, etc.) produce
     a single merged output file in output_dir named with batch_stem. */
}


/* ------------------------------------------------------------------ */
/* Process a merged batch file: run in a forked child                  */
/* ------------------------------------------------------------------ */

static auto process_batch_forked(const daemon_config_t & cfg,
                                   const std::string & merged_input,
                                   const std::string & batch_stem) -> int {
  pid_t pid = fork();
  if (pid < 0) {
    std::fprintf(stderr, "daemon: fork() failed: %s\n", strerror(errno));
    return 1;
  }

  if (pid == 0) {
    /* ---- child ---- */
    clear_output_opts();
    std::vector<std::string> alloc;
    set_output_paths(cfg, batch_stem, alloc);
    search_open_output_files(g_cmdline, g_cmdline);
    int rc = search_process_query_file(merged_input.c_str());
    search_close_output_files();
    std::exit(rc);
  }

  /* ---- parent: wait for child ---- */
  int status = 0;
  while (waitpid(pid, &status, 0) < 0) {
    if (errno != EINTR) {
      std::fprintf(stderr, "daemon: waitpid() error: %s\n", strerror(errno));
      return 1;
    }
  }
  return WIFEXITED(status) ? WEXITSTATUS(status) : 1;
}


/* ------------------------------------------------------------------ */
/* cmd_usearch_global_daemon                                            */
/* ------------------------------------------------------------------ */

auto cmd_usearch_global_daemon(struct Parameters const & parameters,
                                char * cmdline,
                                char * /*progheader*/) -> void
{
  if (parameters.opt_config == nullptr) {
    fatal("--usearch_global_daemon requires --config <yaml>");
  }

  /* Load config (sets opt_* globals for usearch_global parameters). */
  daemon_config_t cfg;
  config_yaml_load(parameters.opt_config, cfg);

  /* CLI --max_batch_sequences overrides config value. */
  if (parameters.opt_max_batch_sequences > 0) {
    cfg.max_batch_sequences = parameters.opt_max_batch_sequences;
  }

  /* Re-apply weak_id clamp: the pre-dispatch fixup in vsearch.cc ran
     before config_yaml_load(), so opt_weak_id was already clamped to
     opt_id's unset sentinel (-1.0).  After config loading opt_id now
     holds the real threshold (e.g. 0.97) but opt_weak_id is still -1.0.
     Any negative value here is the sentinel; reset it to opt_id so that
     only hits meeting the configured id threshold are written to output. */
  if (opt_id >= 0.0 && opt_weak_id < 0.0) {
    opt_weak_id = opt_id;
  }

  /* Validate daemon paths. */
  if (cfg.watch_dir.empty())  { fatal("config: watch_dir is required");  }
  if (cfg.output_dir.empty()) { fatal("config: output_dir is required"); }
  if (cfg.errors_dir.empty()) { fatal("config: errors_dir is required"); }

  ensure_dir(cfg.watch_dir);
  ensure_dir(cfg.output_dir);
  ensure_dir(cfg.errors_dir);

  if (opt_db == nullptr) {
    fatal("--usearch_global_daemon requires --db or db: in config");
  }

  /* Store cmdline for use by forked children. */
  g_cmdline = cmdline;

  /* Load the reference database ONCE. */
  std::fprintf(stderr, "vsearch daemon: loading database %s\n", opt_db);
  search_load_db(cmdline);
  std::fprintf(stderr, "vsearch daemon: database loaded, watching %s\n",
               cfg.watch_dir.c_str());
  if (cfg.max_batch_sequences > 0) {
    std::fprintf(stderr,
                 "vsearch daemon: batch merging enabled (max_batch_sequences: %" PRId64 ")\n",
                 cfg.max_batch_sequences);
  } else {
    std::fprintf(stderr, "vsearch daemon: batch merging enabled (max_batch_sequences: unlimited)\n");
  }

  /* Install signal handlers for clean shutdown. */
  struct sigaction sa;
  std::memset(&sa, 0, sizeof(sa));
  sa.sa_handler = handle_signal;
  sigaction(SIGINT,  &sa, nullptr);
  sigaction(SIGTERM, &sa, nullptr);

  /* Start directory watcher. */
  int kq = dirwatch_init(cfg.watch_dir.c_str());

  /* ---------------------------------------------------------------- */
  /* Main loop                                                          */
  /* ---------------------------------------------------------------- */

  while (g_shutdown == 0) {
    std::vector<std::string> ready;
    dirwatch_poll(kq, cfg.watch_dir.c_str(), ready);

    if (ready.empty()) { continue; }

    /* ---- Count sequences in each ready file ---- */
    std::vector<std::size_t> seq_counts;
    seq_counts.reserve(ready.size());
    for (const auto & f : ready) {
      seq_counts.push_back(count_fasta_sequences(f));
    }

    /* ---- Group files into sub-batches ---- */
    /* Each sub-batch has total sequences <= max_batch_sequences
       (or unlimited when max_batch_sequences == 0). */
    std::vector<std::vector<std::size_t>> batches;  /* indices into ready[] */
    {
      std::vector<std::size_t> cur_batch;
      std::size_t cur_total = 0;
      for (std::size_t i = 0; i < ready.size(); ++i) {
        bool would_overflow = (cfg.max_batch_sequences > 0) &&
                              !cur_batch.empty() &&
                              (cur_total + seq_counts[i] >
                               static_cast<std::size_t>(cfg.max_batch_sequences));
        if (would_overflow) {
          batches.push_back(std::move(cur_batch));
          cur_batch.clear();
          cur_total = 0;
        }
        cur_batch.push_back(i);
        cur_total += seq_counts[i];
      }
      if (!cur_batch.empty()) {
        batches.push_back(std::move(cur_batch));
      }
    }

    /* ---- Process each sub-batch ---- */
    for (const auto & batch_indices : batches) {
      if (g_shutdown != 0) { break; }

      /* Generate a shared timestamp for this batch. */
      char ts[32];
      make_timestamp(ts, sizeof(ts));

      /* ---- Single-file path (no merging overhead) ---- */
      if (batch_indices.size() == 1) {
        const auto & input_path = ready[batch_indices[0]];

        std::string stem    = stem_of(input_path);
        std::string stem_ts = stem + "_" + ts;

        /* Atomically claim the file. */
        std::string staging = input_path + ".processing";
        if (rename(input_path.c_str(), staging.c_str()) != 0) {
          std::fprintf(stderr,
                       "daemon: warning: could not claim %s (skipping): %s\n",
                       input_path.c_str(), strerror(errno));
          continue;
        }

        std::fprintf(stderr, "vsearch daemon: processing %s\n", input_path.c_str());

        int rc = process_file_forked(cfg, staging, stem_ts);

        if (rc == 0) {
          std::string dest = cfg.output_dir + "/" + stem_ts + ".fasta";
          if (rename(staging.c_str(), dest.c_str()) != 0) {
            std::fprintf(stderr,
                         "daemon: warning: could not move %s to %s: %s\n",
                         staging.c_str(), dest.c_str(), strerror(errno));
          }
          std::fprintf(stderr, "vsearch daemon: done \342\206\222 %s\n", dest.c_str());
        } else {
          std::string err_fasta = cfg.errors_dir + "/" + stem_ts + ".fasta";
          std::string err_file  = cfg.errors_dir + "/" + stem_ts + ".err";

          if (rename(staging.c_str(), err_fasta.c_str()) != 0) {
            std::fprintf(stderr,
                         "daemon: warning: could not move %s to %s: %s\n",
                         staging.c_str(), err_fasta.c_str(), strerror(errno));
          }

          std::FILE * ef = std::fopen(err_file.c_str(), "w");
          if (ef != nullptr) {
            std::fprintf(ef, "error processing: %s\nchild exit code: %d\n",
                         input_path.c_str(), rc);
            std::fclose(ef);
          }

          std::fprintf(stderr,
                       "vsearch daemon: error processing %s (exit %d)\n",
                       input_path.c_str(), rc);
        }
        continue;
      }

      /* ---- Multi-file batch path ---- */
      /* Claim all files in the batch first (atomic rename to .processing). */
      std::vector<std::string> staging_paths;
      std::vector<std::string> input_labels;   /* original paths for logging */
      std::vector<std::string> stems_ts;

      for (std::size_t bi : batch_indices) {
        const auto & input_path = ready[bi];
        std::string staging = input_path + ".processing";
        if (rename(input_path.c_str(), staging.c_str()) != 0) {
          std::fprintf(stderr,
                       "daemon: warning: could not claim %s (skipping): %s\n",
                       input_path.c_str(), strerror(errno));
          continue;
        }
        staging_paths.push_back(staging);
        input_labels.push_back(input_path);
        stems_ts.push_back(stem_of(input_path) + "_" + ts);
      }

      if (staging_paths.empty()) { continue; }

      /* If only one file was claimed (others failed), use single-file path. */
      if (staging_paths.size() == 1) {
        std::fprintf(stderr, "vsearch daemon: processing %s\n",
                     input_labels[0].c_str());
        int rc = process_file_forked(cfg, staging_paths[0], stems_ts[0]);
        if (rc == 0) {
          std::string dest = cfg.output_dir + "/" + stems_ts[0] + ".fasta";
          rename(staging_paths[0].c_str(), dest.c_str());
          std::fprintf(stderr, "vsearch daemon: done \342\206\222 %s\n", dest.c_str());
        } else {
          std::string err_fasta = cfg.errors_dir + "/" + stems_ts[0] + ".fasta";
          std::string err_file  = cfg.errors_dir + "/" + stems_ts[0] + ".err";
          rename(staging_paths[0].c_str(), err_fasta.c_str());
          std::FILE * ef = std::fopen(err_file.c_str(), "w");
          if (ef != nullptr) {
            std::fprintf(ef, "error processing: %s\nchild exit code: %d\n",
                         input_labels[0].c_str(), rc);
            std::fclose(ef);
          }
          std::fprintf(stderr, "vsearch daemon: error processing %s (exit %d)\n",
                       input_labels[0].c_str(), rc);
        }
        continue;
      }

      std::fprintf(stderr,
                   "vsearch daemon: batch processing %zu files (%zu–%zu seqs)\n",
                   staging_paths.size(),
                   seq_counts[batch_indices.front()],
                   seq_counts[batch_indices.back()]);

      /* Build merged FASTA in output_dir (dot-prefixed, removed after split). */
      std::string batch_stem   = std::string("._batch_") + ts;
      std::string merged_input = cfg.output_dir + "/" + batch_stem + ".fasta";

      if (!merge_fasta_files(staging_paths, merged_input)) {
        std::fprintf(stderr, "daemon: error: could not create merged batch file\n");
        /* Move all staged files to errors_dir. */
        for (std::size_t i = 0; i < staging_paths.size(); ++i) {
          std::string err_fasta = cfg.errors_dir + "/" + stems_ts[i] + ".fasta";
          std::string err_file  = cfg.errors_dir + "/" + stems_ts[i] + ".err";
          rename(staging_paths[i].c_str(), err_fasta.c_str());
          std::FILE * ef = std::fopen(err_file.c_str(), "w");
          if (ef != nullptr) {
            std::fprintf(ef, "error: batch merge failed for: %s\n",
                         input_labels[i].c_str());
            std::fclose(ef);
          }
        }
        continue;
      }

      /* Search the merged file. */
      int rc = process_batch_forked(cfg, merged_input, batch_stem);
      std::remove(merged_input.c_str());

      if (rc == 0) {
        /* Split merged outputs into per-file outputs. */
        split_batch_outputs(cfg, batch_stem, stems_ts);

        /* Move staged inputs to output_dir. */
        for (std::size_t i = 0; i < staging_paths.size(); ++i) {
          std::string dest = cfg.output_dir + "/" + stems_ts[i] + ".fasta";
          if (rename(staging_paths[i].c_str(), dest.c_str()) != 0) {
            std::fprintf(stderr,
                         "daemon: warning: could not move %s to %s: %s\n",
                         staging_paths[i].c_str(), dest.c_str(), strerror(errno));
          }
          std::fprintf(stderr, "vsearch daemon: done \342\206\222 %s\n", dest.c_str());
        }
      } else {
        /* Batch search failed: move all staged files to errors_dir. */
        for (std::size_t i = 0; i < staging_paths.size(); ++i) {
          std::string err_fasta = cfg.errors_dir + "/" + stems_ts[i] + ".fasta";
          std::string err_file  = cfg.errors_dir + "/" + stems_ts[i] + ".err";
          rename(staging_paths[i].c_str(), err_fasta.c_str());
          std::FILE * ef = std::fopen(err_file.c_str(), "w");
          if (ef != nullptr) {
            std::fprintf(ef, "error processing batch: %s\nchild exit code: %d\n",
                         input_labels[i].c_str(), rc);
            std::fclose(ef);
          }
        }
        std::fprintf(stderr,
                     "vsearch daemon: error processing batch of %zu files (exit %d)\n",
                     staging_paths.size(), rc);
      }
    }  /* for each sub-batch */
  }  /* while !g_shutdown */

  std::fprintf(stderr, "vsearch daemon: shutting down\n");
  dirwatch_close(kq);
  search_unload_db();
}


/* ------------------------------------------------------------------ */
/* cmd_submit_query                                                     */
/* ------------------------------------------------------------------ */

auto cmd_submit_query(struct Parameters const & parameters) -> void
{
  if (parameters.opt_submit_query == nullptr) {
    fatal("--submit_query requires an input file path");
  }
  if (parameters.opt_config == nullptr) {
    fatal("--submit_query requires --config <yaml>");
  }

  daemon_config_t cfg;
  config_yaml_load(parameters.opt_config, cfg);

  if (cfg.watch_dir.empty()) {
    fatal("config: watch_dir is required for --submit_query");
  }

  const char * src = parameters.opt_submit_query;

  /* Verify source file exists. */
  struct stat st;
  if (stat(src, &st) != 0) {
    fatal("submit_query: input file not found", src);
  }

  char ts[32];
  make_timestamp(ts, sizeof(ts));

  std::string stem    = stem_of(std::string(src));
  std::string final_name = stem + "_" + ts + ".fasta";
  std::string final_path = cfg.watch_dir + "/" + final_name;

  if (!copy_file_atomic(src, final_path.c_str())) {
    fatal("submit_query: failed to copy file to watch directory");
  }

  std::fprintf(stderr, "submitted: %s\n", final_path.c_str());
}
