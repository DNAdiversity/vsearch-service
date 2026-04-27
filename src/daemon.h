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

struct Parameters;

/* Daemon command: --usearch_global_daemon <db> --config <yaml>
   Loads the reference database once, then watches the configured
   input directory and processes arriving FASTA files one at a time.
   Runs in the foreground; logs to stderr.
   Exits cleanly on SIGINT or SIGTERM. */
auto cmd_usearch_global_daemon(struct Parameters const & parameters,
                                char * cmdline,
                                char * progheader) -> void;

/* Submit command: --submit_query <input.fasta> --config <yaml>
   Copies the input file into the configured watch directory using an
   atomic rename, adding a timestamp to avoid filename collisions.
   The daemon will detect the file via the directory watcher. */
auto cmd_submit_query(struct Parameters const & parameters) -> void;
