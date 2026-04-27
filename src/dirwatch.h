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

#pragma once

/*
  Directory watcher — macOS implementation using kqueue.

  The watcher detects when new files appear in a directory.
  Files with a ".tmp" suffix are ignored (they are in-flight writes).
  Recognised FASTA extensions: .fasta  .fa  .fna

  Usage:
    int kq = dirwatch_init(watch_dir);
    while (!done) {
        std::vector<std::string> ready;
        dirwatch_poll(kq, watch_dir, ready);   // blocks up to ~1 s
        for (const auto & path : ready) {
            process(path);
        }
    }
    dirwatch_close(kq);

  Future: a Linux implementation using inotify can be substituted
  behind the same interface by providing dirwatch_linux.cc compiled
  when HAVE_SYS_INOTIFY_H is defined.
*/

#include <string>
#include <vector>

/* Open a kqueue and attach a EVFILT_VNODE watch on dir_path.
   Returns the kqueue file descriptor.  Calls fatal() on error. */
auto dirwatch_init(const char * dir_path) -> int;

/* Block until the directory changes (up to ~1 second), then scan it.
   Appends to ready_files the full paths of all .fasta/.fa/.fna files
   that are NOT .tmp, sorted oldest-first by mtime.
   Returns the number of files appended (may be 0 on timeout). */
auto dirwatch_poll(int kq_fd,
                   const char * dir_path,
                   std::vector<std::string> & ready_files) -> int;

/* Close the kqueue file descriptor. */
auto dirwatch_close(int kq_fd) -> void;
