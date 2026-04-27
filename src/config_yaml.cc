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
  Flat key:value YAML config parser for vsearch daemon mode.

  Supported syntax:
    key: value        (leading/trailing whitespace stripped from both)
    # comment         (entire line ignored)
    blank lines       (ignored)

  No nested structures, lists, or multi-line values are supported —
  all usearch_global parameters are scalars so a flat parser suffices.
*/

#include "vsearch.h"
#include "config_yaml.h"
#include "mask.h"
#include "utils/fatal.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

/* Build a fatal message and call fatal(). */
static auto config_fatal(const char * fmt, const char * a, const char * b) -> void
{
  char msg[512];
  std::snprintf(msg, sizeof(msg), fmt, a, b);
  fatal(msg);
}


/* ------------------------------------------------------------------ */
/* Internal helpers                                                     */
/* ------------------------------------------------------------------ */

static auto trim(std::string & s) -> void {
  const char * ws = " \t\r\n";
  auto first = s.find_first_not_of(ws);
  if (first == std::string::npos) {
    s.clear();
    return;
  }
  auto last = s.find_last_not_of(ws);
  s = s.substr(first, last - first + 1);
}

/* Parse a boolean value: true/yes/1 → true, false/no/0 → false.
   Calls fatal() on unrecognised value. */
static auto parse_bool(const std::string & key, const std::string & val) -> bool {
  if (val == "true" || val == "yes" || val == "1") {
    return true;
  }
  if (val == "false" || val == "no" || val == "0") {
    return false;
  }
  config_fatal("config: key '%s' expects true/false/yes/no/0/1, got: %s",
               key.c_str(), val.c_str());
  return false; // unreachable
}

/* Parse an int64_t value.  Calls fatal() on parse error. */
static auto parse_int64(const std::string & key,
                        const std::string & val) -> int64_t {
  char * end = nullptr;
  errno = 0;
  long long v = std::strtoll(val.c_str(), &end, 10);
  if (errno != 0 || end == val.c_str() || *end != '\0') {
    config_fatal("config: key '%s' expects an integer, got: %s",
                 key.c_str(), val.c_str());
  }
  return static_cast<int64_t>(v);
}

/* Parse a double value.  Calls fatal() on parse error. */
static auto parse_double(const std::string & key,
                         const std::string & val) -> double {
  char * end = nullptr;
  errno = 0;
  double v = std::strtod(val.c_str(), &end);
  if (errno != 0 || end == val.c_str() || *end != '\0') {
    config_fatal("config: key '%s' expects a number, got: %s",
                 key.c_str(), val.c_str());
  }
  return v;
}

/* Parse a mask value: none/0 → MASK_NONE, dust/1 → MASK_DUST,
   soft/2 → MASK_SOFT.  Calls fatal() on unrecognised value. */
static auto parse_mask(const std::string & key,
                       const std::string & val) -> int64_t {
  if (val == "none" || val == "0") { return MASK_NONE; }
  if (val == "dust" || val == "1") { return MASK_DUST; }
  if (val == "soft" || val == "2") { return MASK_SOFT; }
  config_fatal("config: key '%s' expects none/dust/soft, got: %s",
               key.c_str(), val.c_str());
  return MASK_NONE; // unreachable
}

/* Duplicate a std::string to a heap char* (like strdup).
   The caller is responsible for the memory; in practice vsearch keeps
   these alive for the entire run. */
static auto dup_string(const std::string & s) -> char * {
  char * p = static_cast<char *>(std::malloc(s.size() + 1));
  if (p == nullptr) {
    fatal("config: out of memory");
  }
  std::memcpy(p, s.c_str(), s.size() + 1);
  return p;
}


/* ------------------------------------------------------------------ */
/* Public API                                                           */
/* ------------------------------------------------------------------ */

auto config_yaml_load(const char * path, daemon_config_t & daemon_cfg) -> void
{
  std::FILE * fp = std::fopen(path, "r");
  if (fp == nullptr) {
    fatal("Unable to open config file", path);
  }

  char linebuf[4096];
  int lineno = 0;

  while (std::fgets(linebuf, static_cast<int>(sizeof(linebuf)), fp) != nullptr)
    {
      ++lineno;
      std::string line(linebuf);
      trim(line);

      /* skip blanks and comments */
      if (line.empty() || line[0] == '#') {
        continue;
      }

      /* find the first ':' separator */
      auto colon = line.find(':');
      if (colon == std::string::npos) {
        std::fprintf(stderr,
                     "config warning: no ':' on line %d, skipping: %s\n",
                     lineno, linebuf);
        continue;
      }

      std::string key = line.substr(0, colon);
      std::string val = line.substr(colon + 1);
      trim(key);
      trim(val);

      if (key.empty()) {
        std::fprintf(stderr,
                     "config warning: empty key on line %d, skipping\n",
                     lineno);
        continue;
      }

      /* ------------------------------------------------------------ */
      /* Daemon-specific settings                                       */
      /* ------------------------------------------------------------ */

      if (key == "watch_dir")            { daemon_cfg.watch_dir             = val; continue; }
      if (key == "output_dir")           { daemon_cfg.output_dir            = val; continue; }
      if (key == "errors_dir")           { daemon_cfg.errors_dir            = val; continue; }
      if (key == "max_batch_sequences")  { daemon_cfg.max_batch_sequences   = parse_int64(key, val); continue; }

      /* Output format enable flags */
      if (key == "output_alnout")            { daemon_cfg.enable_alnout            = parse_bool(key, val); continue; }
      if (key == "output_blast6out")         { daemon_cfg.enable_blast6out         = parse_bool(key, val); continue; }
      if (key == "output_uc")                { daemon_cfg.enable_uc                = parse_bool(key, val); continue; }
      if (key == "output_userout")           { daemon_cfg.enable_userout           = parse_bool(key, val); continue; }
      if (key == "output_samout")            { daemon_cfg.enable_samout            = parse_bool(key, val); continue; }
      if (key == "output_matched")           { daemon_cfg.enable_matched           = parse_bool(key, val); continue; }
      if (key == "output_notmatched")        { daemon_cfg.enable_notmatched        = parse_bool(key, val); continue; }
      if (key == "output_fastapairs")        { daemon_cfg.enable_fastapairs        = parse_bool(key, val); continue; }
      if (key == "output_lcaout")            { daemon_cfg.enable_lcaout            = parse_bool(key, val); continue; }
      if (key == "output_qsegout")           { daemon_cfg.enable_qsegout           = parse_bool(key, val); continue; }
      if (key == "output_tsegout")           { daemon_cfg.enable_tsegout           = parse_bool(key, val); continue; }
      if (key == "output_otutabout")         { daemon_cfg.enable_otutabout         = parse_bool(key, val); continue; }
      if (key == "output_mothur_shared_out") { daemon_cfg.enable_mothur_shared_out = parse_bool(key, val); continue; }
      if (key == "output_biomout")           { daemon_cfg.enable_biomout           = parse_bool(key, val); continue; }

      /* ------------------------------------------------------------ */
      /* Standard usearch_global parameters → global opt_* variables  */
      /* ------------------------------------------------------------ */

      /* string options */
      if (key == "db")               { opt_db               = dup_string(val); continue; }
      if (key == "log")              { opt_log              = dup_string(val); continue; }
      if (key == "alnout")           { opt_alnout           = dup_string(val); continue; }
      if (key == "blast6out")        { opt_blast6out        = dup_string(val); continue; }
      if (key == "uc")               { opt_uc               = dup_string(val); continue; }
      if (key == "userout")          { opt_userout          = dup_string(val); continue; }
      if (key == "samout")           { opt_samout           = dup_string(val); continue; }
      if (key == "matched")          { opt_matched          = dup_string(val); continue; }
      if (key == "notmatched")       { opt_notmatched       = dup_string(val); continue; }
      if (key == "dbmatched")        { opt_dbmatched        = dup_string(val); continue; }
      if (key == "dbnotmatched")     { opt_dbnotmatched     = dup_string(val); continue; }
      if (key == "fastapairs")       { opt_fastapairs       = dup_string(val); continue; }
      if (key == "otutabout")        { opt_otutabout        = dup_string(val); continue; }
      if (key == "mothur_shared_out"){ opt_mothur_shared_out= dup_string(val); continue; }
      if (key == "biomout")          { opt_biomout          = dup_string(val); continue; }
      if (key == "lcaout")           { opt_lcaout           = dup_string(val); continue; }
      if (key == "qsegout")          { opt_qsegout          = dup_string(val); continue; }
      if (key == "tsegout")          { opt_tsegout          = dup_string(val); continue; }

      /* double options */
      if (key == "id")              { opt_id              = parse_double(key, val); continue; }
      if (key == "maxid")           { opt_maxid           = parse_double(key, val); continue; }
      if (key == "weak_id")         { opt_weak_id         = parse_double(key, val); continue; }
      if (key == "query_cov")       { opt_query_cov       = parse_double(key, val); continue; }
      if (key == "target_cov")      { opt_target_cov      = parse_double(key, val); continue; }
      if (key == "mid")             { opt_mid             = parse_double(key, val); continue; }
      if (key == "maxqt")           { opt_maxqt           = parse_double(key, val); continue; }
      if (key == "minqt")           { opt_minqt           = parse_double(key, val); continue; }
      if (key == "maxsl")           { opt_maxsl           = parse_double(key, val); continue; }
      if (key == "minsl")           { opt_minsl           = parse_double(key, val); continue; }

      /* int64 options */
      if (key == "maxaccepts")      { opt_maxaccepts      = parse_int64(key, val); continue; }
      if (key == "maxrejects")      { opt_maxrejects      = parse_int64(key, val); continue; }
      if (key == "maxhits")         { opt_maxhits         = parse_int64(key, val); continue; }
      if (key == "strand")          { opt_strand          = parse_int64(key, val); continue; }
      if (key == "threads")         { opt_threads         = parse_int64(key, val); continue; }
      if (key == "wordlength")      { opt_wordlength      = parse_int64(key, val); continue; }
      if (key == "iddef")           { opt_iddef           = parse_int64(key, val); continue; }
      if (key == "mincols")         { opt_mincols         = parse_int64(key, val); continue; }
      if (key == "minseqlength")    { opt_minseqlength    = parse_int64(key, val); continue; }
      if (key == "maxseqlength")    { opt_maxseqlength    = parse_int64(key, val); continue; }
      if (key == "maxdiffs")        { opt_maxdiffs        = parse_int64(key, val); continue; }
      if (key == "maxgaps")         { opt_maxgaps         = parse_int64(key, val); continue; }
      if (key == "maxsubs")         { opt_maxsubs         = parse_int64(key, val); continue; }
      if (key == "hardmask")        { opt_hardmask        = parse_int64(key, val); continue; }
      if (key == "fulldp")          { opt_fulldp          = parse_int64(key, val); continue; }
      if (key == "self")            { opt_self            = parse_int64(key, val); continue; }
      if (key == "selfid")          { opt_selfid          = parse_int64(key, val); continue; }
      if (key == "output_no_hits")  { opt_output_no_hits  = parse_int64(key, val); continue; }
      if (key == "top_hits_only")   { opt_top_hits_only   = parse_int64(key, val); continue; }
      if (key == "uc_allhits")      { opt_uc_allhits      = parse_int64(key, val); continue; }
      if (key == "minwordmatches")  { opt_minwordmatches  = parse_int64(key, val); continue; }
      if (key == "notrunclabels")   { opt_notrunclabels   = parse_int64(key, val); continue; }
      if (key == "randseed")        { opt_randseed        = parse_int64(key, val); continue; }
      if (key == "rowlen")          { opt_rowlen          = parse_int64(key, val); continue; }
      if (key == "alignwidth")      { opt_alignwidth      = parse_int64(key, val); continue; }

      /* mask options (string → int constant) */
      if (key == "dbmask") { opt_dbmask = parse_mask(key, val); continue; }
      if (key == "qmask")  { opt_qmask  = parse_mask(key, val); continue; }

      /* bool options */
      if (key == "sizein")         { opt_sizein         = parse_bool(key, val); continue; }
      if (key == "sizeout")        { opt_sizeout        = parse_bool(key, val); continue; }
      if (key == "quiet")          { opt_quiet          = parse_bool(key, val); continue; }
      if (key == "no_progress")    { opt_no_progress    = parse_bool(key, val); continue; }
      if (key == "samheader")      { opt_samheader      = parse_bool(key, val); continue; }
      if (key == "lengthout")      { opt_lengthout      = parse_bool(key, val); continue; }
      if (key == "sizeorder")      { opt_sizeorder      = parse_bool(key, val); continue; }
      if (key == "xsize")          { opt_xsize          = parse_bool(key, val); continue; }
      if (key == "eeout")          { opt_eeout          = parse_bool(key, val); continue; }

      /* gap parameters (int) */
      if (key == "gap_open_query_interior")    { opt_gap_open_query_interior    = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_open_query_left")        { opt_gap_open_query_left        = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_open_query_right")       { opt_gap_open_query_right       = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_open_target_interior")   { opt_gap_open_target_interior   = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_open_target_left")       { opt_gap_open_target_left       = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_open_target_right")      { opt_gap_open_target_right      = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_extension_query_interior")  { opt_gap_extension_query_interior  = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_extension_query_left")      { opt_gap_extension_query_left      = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_extension_query_right")     { opt_gap_extension_query_right     = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_extension_target_interior") { opt_gap_extension_target_interior = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_extension_target_left")     { opt_gap_extension_target_left     = static_cast<int>(parse_int64(key, val)); continue; }
      if (key == "gap_extension_target_right")    { opt_gap_extension_target_right    = static_cast<int>(parse_int64(key, val)); continue; }

      /* unrecognised key */
      std::fprintf(stderr,
                   "config warning: unknown key '%s' on line %d, ignoring\n",
                   key.c_str(), lineno);
    }

  std::fclose(fp);
}
