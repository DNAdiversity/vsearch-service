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

#include <cstdint>
#include <string>

/* Daemon-specific configuration populated by config_yaml_load().
   Standard usearch_global parameters are written directly into the
   global opt_* variables. */
struct daemon_config_t {
  /* Required daemon paths */
  std::string watch_dir;
  std::string output_dir;
  std::string errors_dir;

  /* Batch merging: maximum total sequences merged into one search.
     Files are never split — if adding the next file would exceed this limit,
     the current batch is flushed and the new file starts a fresh batch.
     A single file larger than this limit is processed alone.
     0 = unlimited (merge all queued files into one batch). */
  int64_t max_batch_sequences = 1000;

  /* Output format flags: daemon generates one file per input query
     for each enabled format, named: stem_TIMESTAMP.ext */
  bool enable_alnout          = false;
  bool enable_blast6out       = false;
  bool enable_uc              = false;
  bool enable_userout         = false;
  bool enable_samout          = false;
  bool enable_matched         = false;
  bool enable_notmatched      = false;
  bool enable_fastapairs      = false;
  bool enable_lcaout          = false;
  bool enable_qsegout         = false;
  bool enable_tsegout         = false;
  bool enable_otutabout       = false;
  bool enable_mothur_shared_out = false;
  bool enable_biomout         = false;
};

/* Load a flat key:value YAML config file.
   - Standard usearch_global parameters are applied to the global opt_*
     variables immediately (they can be overridden by subsequent CLI parsing).
   - Daemon-specific settings and output-format flags are returned via
     the daemon_config_t struct.
   - Unknown keys produce a warning to stderr.
   - Calls fatal() on file-open error or malformed value. */
auto config_yaml_load(const char * path, daemon_config_t & daemon_cfg) -> void;
