# vsearch Daemon Mode

## Overview

Loading a large reference database into memory and building the k-mer index can take 30–60 seconds or more for typical metagenomic reference sets. For workflows that search many query files against the same reference, this cost is paid on every invocation — wasted time that scales linearly with the number of searches.

Daemon mode solves this by loading the reference database **once** and keeping it resident in memory, then watching a directory for incoming FASTA files and searching them as they arrive. Combined with **batch merging**, multiple queued files are merged into a single search, fully saturating all worker threads without repeated per-search setup overhead.

---

## Architecture

### Core components

| File | Role |
|------|------|
| `src/daemon.cc` | Main daemon loop, batch merging, output splitting, submit command |
| `src/dirwatch.cc` | kqueue-based directory watcher (macOS) |
| `src/config_yaml.cc` | Flat YAML config parser |
| `src/search.cc` | Exposes `search_load_db()`, `search_unload_db()`, `search_process_query_file()` |

### Why fork()?

vsearch's `fatal()` calls `std::exit()` — it cannot be caught. A single malformed query file would kill the daemon and lose the loaded database. Each file (or batch) is processed in a forked child process. The child inherits the in-memory database via copy-on-write: because search workers only read the database, no physical copies are made and memory usage stays at ~1×.

### Atomic file handoff (submit_query)

Writers use `--submit_query` to drop files into the watch directory:

1. File is copied to `watch_dir/stem_TIMESTAMP.fasta.tmp`
2. Atomically renamed to `watch_dir/stem_TIMESTAMP.fasta`

The daemon's directory watcher skips `.tmp` files, so it only ever sees complete files. When a file is claimed for processing, it is immediately renamed to `stem_TIMESTAMP.fasta.processing`, preventing same-second resubmits from overwriting it while the child is running.

### Batch merging

When multiple files queue up before the previous batch finishes, they are merged into a single FASTA before searching:

1. Each query header is prefixed with a `VSDF####_` tag (e.g., `>VSDF0002_original_label`) where `####` is the zero-padded file index within the batch.
2. The merged FASTA is searched in one child process, saturating all worker threads continuously.
3. After the child exits, output files are scanned for the `VSDF####_` prefix and each result line is routed to the correct per-file output, stripping the prefix before writing.

**Output formats with per-file splitting:**
- `blast6out` — query label in field 0
- `uc` — query label in field 8
- `matched` / `notmatched` — FASTA format, query label on `>` line
- `alnout` — multi-line format, routed by `Query >VSDF####_` lines

**Other formats** (`userout`, `samout`, `lcaout`, etc.) produce a single merged output file named `._batch_TIMESTAMP.ext` in `output_dir`. Submit one file at a time if per-file split output is required for these formats.

---

## Configuration

All options can be set via a YAML config file and optionally overridden on the command line.

### Example config file

```yaml
# Reference database
db: /path/to/reference.fasta

# Daemon directories
watch_dir:  /path/to/watch
output_dir: /path/to/output
errors_dir: /path/to/errors

# Batch merging cap (default: 1000; 0 = unlimited)
max_batch_sequences: 1000

# usearch_global search parameters
id:         0.97
maxaccepts: 1
maxrejects: 32
threads:    8
wordlength: 8
dbmask:     dust
strand:     1        # 1 = plus only, 3 = both strands

# Output formats to produce (true/false)
output_uc:       true
output_blast6out: false
output_matched:  false
```

### Config keys

| Key | Type | Description |
|-----|------|-------------|
| `db` | string | Path to reference FASTA database |
| `watch_dir` | string | Directory to watch for incoming query files |
| `output_dir` | string | Directory for completed results |
| `errors_dir` | string | Directory for failed queries |
| `max_batch_sequences` | integer | Max total sequences per merged batch (default: 1000; 0 = unlimited). Files are never split — if the next file would exceed the limit, the current batch flushes and a new one starts. A single file larger than this limit is processed alone. |
| `id` | float | Minimum identity threshold (e.g. `0.97`) |
| `maxaccepts` | integer | Stop per-query search after N accepted hits |
| `maxrejects` | integer | Stop per-query search after N rejected candidates |
| `threads` | integer | Worker threads (0 = auto-detect) |
| `strand` | integer | 1 = plus only, 3 = both strands |
| `wordlength` | integer | k-mer word length for pre-filter |
| `dbmask` | string | `none`, `dust`, or `soft` |
| `output_uc` | bool | Enable UC output |
| `output_blast6out` | bool | Enable BLAST tabular output |
| `output_matched` | bool | Enable matched sequences FASTA |
| `output_notmatched` | bool | Enable unmatched sequences FASTA |
| `output_alnout` | bool | Enable alignment output |

---

## Usage

### 1. Start the daemon

```bash
vsearch --usearch_global_daemon /dev/null \
        --config /path/to/config.yaml
```

The `--usearch_global_daemon` argument is vestigial (the database path comes from `db:` in the config). The daemon logs to stderr and runs until it receives `SIGINT` or `SIGTERM`.

**Override config values on the command line:**

```bash
vsearch --usearch_global_daemon /dev/null \
        --config config.yaml \
        --threads 16 \
        --max_batch_sequences 200
```

Command-line flags take precedence over config file values.

### 2. Submit a query file

```bash
vsearch --submit_query /path/to/query.fasta \
        --config /path/to/config.yaml
```

This copies `query.fasta` into `watch_dir` as `query_TIMESTAMP.fasta` using an atomic `.tmp` → rename handoff. The daemon picks it up within ~1 second.

### 3. File lifecycle

```
Submitted:   watch_dir/query_20260426T120000.fasta
Claimed:     watch_dir/query_20260426T120000.fasta.processing
On success:  output_dir/query_20260426T120000.fasta
             output_dir/query_20260426T120000.uc
             output_dir/query_20260426T120000.blast6out  (etc.)
On failure:  errors_dir/query_20260426T120000.fasta
             errors_dir/query_20260426T120000.err
```

### 4. Shutdown

```bash
kill -SIGTERM <daemon_pid>
```

The daemon finishes any in-flight batch, then exits cleanly.

---

## Benchmark Results

Measured on macOS with the MIDORI2 LONGEST CO1 reference database:
- **237,000 sequences, 226 MB**
- DB load time: ~39.5 s
- Query: 8 sequences submitted 20 times
- Thread count: 4 (of 8 logical cores available)

### Raw numbers

| Mode | Total (20 runs) | Per-run average |
|------|----------------|-----------------|
| `--usearch_global` (cold, DB reloaded each run) | 792 s (13.2 min) | 39.6 s |
| Daemon warm search | 10.9 s | **0.55 s** |
| Daemon total (load + 20 searches) | 50.5 s | — |

**Overall speedup (daemon vs loop, 20 runs): 15.7×**
**Warm search speedup (excluding DB load): 72.5×**

### Speedup by number of runs (batch size = 8 sequences, 4 threads)

| Runs | Loop total | Daemon total | Speedup | Time saved |
|------|-----------|--------------|---------|------------|
| 20 | 13.2 min | 50.5 s | 15.7× | 12.4 min |
| 100 | 1.10 hr | 94.1 s | 42.1× | 1.07 hr |
| 200 | 2.20 hr | 2.5 min | 53.3× | 2.16 hr |
| 400 | 4.40 hr | 4.3 min | 61.4× | 4.33 hr |
| 800 | 8.80 hr | 7.9 min | 66.5× | 8.67 hr |
| 1600 | 17.6 hr | 15.2 min | 69.4× | 17.35 hr |

Asymptote (N→∞): **72.5×** (= 39.6 s per cold run ÷ 0.55 s per warm search)

### Speedup by batch size (8 → 1600 sequences per file, 20 runs, 4 threads)

As batch size grows, per-search time dominates over DB-load time, narrowing the daemon's advantage:

| Sequences/file | Loop total | Daemon total | Speedup |
|----------------|-----------|--------------|---------|
| 8 | 13.4 min | 50.5 s | 15.9× |
| 100 | 15.5 min | 2.93 min | 5.3× |
| 200 | 17.7 min | 5.21 min | 3.4× |
| 400 | 22.3 min | 9.76 min | 2.3× |
| 800 | 31.4 min | 18.9 min | 1.66× |
| 1600 | 49.6 min | 37.1 min | 1.34× |

**The daemon is most valuable for high-frequency submission of small query files**, which is the typical interactive metagenomics use case.

### Effect of thread count

The speedup ceiling for N runs is N× (achieved when warm search time → 0). More threads compress warm search time toward zero, pushing closer to this ceiling:

| Threads | Warm search (8 seqs) | Daemon total (20 runs) | Speedup |
|---------|---------------------|----------------------|---------|
| 4 | 0.55 s | 50.5 s | 15.9× |
| 8 | ~0.27 s | 44.9 s | 17.7× |
| ∞ | ~0 s | 39.5 s (load only) | 20.0× |

More threads benefit the daemon disproportionately: loop mode's per-run savings are offset by the repeated DB load, while the daemon's one-time load is unaffected.

---

## Changes Made

### New files

| File | Description |
|------|-------------|
| `src/daemon.cc` | Daemon loop, submit command, batch merge + split |
| `src/daemon.h` | Public API declarations |
| `src/dirwatch.cc` | kqueue directory watcher (macOS) |
| `src/dirwatch.h` | Public API |
| `src/config_yaml.cc` | Flat YAML config parser |
| `src/config_yaml.h` | `daemon_config_t` struct and parser declaration |

### Modified files

| File | Change |
|------|--------|
| `src/search.cc` | Extracted `search_load_db()`, `search_unload_db()`, `search_open_output_files()`, `search_close_output_files()`, `search_process_query_file()` |
| `src/search.h` | Declarations for extracted functions |
| `src/vsearch.cc` | Added `--usearch_global_daemon`, `--submit_query`, `--config`, `--max_batch_sequences` options; `number_of_options` bumped to 251 |
| `src/vsearch.h` | Added `opt_usearch_global_daemon`, `opt_submit_query`, `opt_config`, `opt_max_batch_sequences` to Parameters struct |
| `src/Makefile.am` | Added new source files |
| `configure.ac` | Added `AC_CHECK_HEADERS([sys/event.h])` |

### Key design decisions

1. **`.processing` staging**: Input files are renamed to `.processing` before forking. This prevents same-second resubmits (same filename) from overwriting the file while a child is running. `dirwatch.cc` skips `.processing` files during directory scans.

2. **`opt_weak_id` sentinel fix**: vsearch's pre-dispatch fixup sets `opt_weak_id = opt_id` before `opt_id` is loaded from config. In daemon mode this left `opt_weak_id = -1.0`, causing all sub-threshold hits to be written as "weak" hits. Fixed by re-clamping after `config_yaml_load()`:
   ```cpp
   if (opt_id >= 0.0 && opt_weak_id < 0.0) { opt_weak_id = opt_id; }
   ```

3. **`strand` in YAML**: The YAML parser uses `parse_int64()` for `strand`, so use `strand: 1` (plus only) or `strand: 3` (both) rather than `strand: plus`.

4. **Batch temp files**: Merged FASTA and output files use a dot-prefix name (`._batch_TIMESTAMP.*`) in `output_dir`, making them invisible to most directory listings. They are deleted immediately after splitting.
