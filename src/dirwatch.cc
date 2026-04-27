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
  kqueue-based directory watcher for macOS.

  Principle:
  - We open the directory with open() and attach EVFILT_VNODE + NOTE_WRITE.
  - kevent() wakes when anything inside the directory changes (rename, write,
    create, delete).
  - On wake (or on timeout), we scan the directory with opendir()/readdir(),
    collecting all .fasta/.fa/.fna files that are NOT .tmp.
  - We sort collected entries by mtime ascending (oldest first) so the daemon
    processes files in arrival order.

  The 1-second kevent() timeout means the daemon's shutdown flag is checked
  at most once per second.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "vsearch.h"
#include "dirwatch.h"
#include "utils/fatal.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef __APPLE__
#include <sys/event.h>
#include <sys/time.h>
#endif


/* ------------------------------------------------------------------ */
/* Helpers                                                              */
/* ------------------------------------------------------------------ */

/* Returns true if the filename has a recognised FASTA extension. */
static auto is_fasta_extension(const char * name) -> bool {
  const char * dot = std::strrchr(name, '.');
  if (dot == nullptr) { return false; }
  return (std::strcmp(dot, ".fasta") == 0 ||
          std::strcmp(dot, ".fa")    == 0 ||
          std::strcmp(dot, ".fna")   == 0);
}

/* Returns true if the filename ends in ".tmp" or ".processing".
   Both are transient states: .tmp = writer not finished,
   .processing = daemon claimed the file but child hasn't finished yet. */
static auto is_tmp(const char * name) -> bool {
  const char * dot = std::strrchr(name, '.');
  if (dot == nullptr) { return false; }
  return (std::strcmp(dot, ".tmp") == 0 ||
          std::strcmp(dot, ".processing") == 0);
}

struct dir_entry_t {
  std::string path;
  time_t      mtime;
};


/* ------------------------------------------------------------------ */
/* kq → dir_fd mapping so dirwatch_close can clean up properly          */
/* ------------------------------------------------------------------ */

#ifdef __APPLE__
static std::unordered_map<int, int> kq_to_dir_fd;
#endif


/* ------------------------------------------------------------------ */
/* Public API                                                           */
/* ------------------------------------------------------------------ */

auto dirwatch_init(const char * dir_path) -> int
{
#ifdef __APPLE__
  int dir_fd = open(dir_path, O_RDONLY);
  if (dir_fd < 0) {
    fatal("dirwatch: cannot open directory", dir_path);
  }

  int kq = kqueue();
  if (kq < 0) {
    close(dir_fd);
    fatal("dirwatch: kqueue() failed");
  }

  struct kevent ev;
  EV_SET(&ev, dir_fd, EVFILT_VNODE,
         EV_ADD | EV_CLEAR,
         NOTE_WRITE | NOTE_EXTEND | NOTE_RENAME | NOTE_DELETE,
         0, nullptr);

  if (kevent(kq, &ev, 1, nullptr, 0, nullptr) < 0) {
    close(dir_fd);
    close(kq);
    fatal("dirwatch: kevent registration failed");
  }

  kq_to_dir_fd[kq] = dir_fd;
  return kq;

#else
  (void)dir_path;
  fatal("dirwatch: directory watching is not supported on this platform");
  return -1;
#endif
}


auto dirwatch_poll(int kq_fd,
                   const char * dir_path,
                   std::vector<std::string> & ready_files) -> int
{
#ifdef __APPLE__
  /* Wait up to 1 second for a directory change event. */
  struct timespec timeout;
  timeout.tv_sec  = 1;
  timeout.tv_nsec = 0;

  struct kevent triggered;
  /* Ignore the return value — even on timeout (ret==0) we scan,
     because files may have been dropped before we started watching. */
  kevent(kq_fd, nullptr, 0, &triggered, 1, &timeout);

  /* Scan the directory for ready files. */
  DIR * dir = opendir(dir_path);
  if (dir == nullptr) {
    return 0;
  }

  std::vector<dir_entry_t> entries;
  struct dirent * de = nullptr;
  while ((de = readdir(dir)) != nullptr)
    {
      const char * name = de->d_name;
      if (name[0] == '.') { continue; }        /* skip . .. and hidden */
      if (is_tmp(name))   { continue; }        /* skip in-flight writes */
      if (!is_fasta_extension(name)) { continue; }

      std::string full = std::string(dir_path) + "/" + name;
      struct stat st;
      if (stat(full.c_str(), &st) != 0) { continue; }

      entries.push_back({full, st.st_mtime});
    }
  closedir(dir);

  /* Sort oldest-first. */
  std::sort(entries.begin(), entries.end(),
            [](const dir_entry_t & a, const dir_entry_t & b) {
              return a.mtime < b.mtime;
            });

  for (const auto & e : entries) {
    ready_files.push_back(e.path);
  }

  return static_cast<int>(entries.size());

#else
  (void)kq_fd;
  (void)dir_path;
  (void)ready_files;
  fatal("dirwatch: directory watching is not supported on this platform");
  return -1;
#endif
}


auto dirwatch_close(int kq_fd) -> void
{
#ifdef __APPLE__
  auto it = kq_to_dir_fd.find(kq_fd);
  if (it != kq_to_dir_fd.end()) {
    close(it->second);  /* close the watched directory fd */
    kq_to_dir_fd.erase(it);
  }
  close(kq_fd);
#else
  (void)kq_fd;
#endif
}
