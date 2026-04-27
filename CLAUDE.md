# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

VSEARCH is a high-performance C++ bioinformatics tool for metagenomic sequence analysis — an open-source alternative to USEARCH. It supports chimera detection, clustering, dereplication, global sequence search, paired-end merging, sequence masking, and taxonomic classification (SINTAX).

It can be used both as a CLI tool and as an embeddable static library (`libvsearch.a`). See `LIBRARY_API.md` and `src/vsearch_api.h` for the library API.

## Build Commands

```bash
# First-time setup (from a fresh git clone)
./autogen.sh

# Configure (optimization flags recommended)
./configure CFLAGS="-O2" CXXFLAGS="-O2"

# Build the CLI tool
make ARFLAGS="cr" -j$(nproc)

# Build only the static library (for embedding)
make -C src libvsearch.a
```

**Optional configure flags:**
- `--enable-debug` — enable debug build with extra compiler warnings
- `--enable-profiling` — enable profiling build
- `--disable-zlib` / `--disable-bzip2` — disable compression support
- `--prefix=DIR` — customize installation prefix

## Testing

Tests live in a **separate repository** (`vsearch-tests`):

```bash
git clone https://github.com/frederic-mahe/vsearch-tests
cd vsearch-tests
bash ./run_all_tests.sh
```

Tests are bash scripts that run vsearch commands and diff output against ground-truth files. To run a single test, execute its individual bash script directly. CI uses `.github/workflows/build-and-test.yml`.

There is no separate lint command — compiler warnings (`-Wall -Wextra -Wpedantic`) are enforced at build time.

## Architecture

**Entry point:** `src/vsearch.cc` — parses ~200 `opt_*` global config variables from CLI args, then dispatches to the appropriate operation module.

**Key architectural layers:**

| Layer | Files | Role |
|-------|-------|------|
| CLI / dispatch | `vsearch.cc` | Argument parsing, global `opt_*` vars, operation dispatch |
| Database | `db.cc`, `dbhash.cc`, `dbindex.cc` | Load sequences into memory, build k-mer indices |
| Alignment core | `align_simd.cc` | SIMD-vectorized global alignment (1 query vs 8 DB seqs simultaneously) |
| Search | `search.cc`, `searchcore.cc`, `search_exact.cc` | Needle-Wunsch global alignment search, exact matching |
| Operations | `cluster.cc`, `chimera.cc`, `derep.cc`, `mergepairs.cc`, `sintax.cc`, etc. | One file per major vsearch operation |
| File I/O | `fasta.cc`, `fastq.cc`, `fastx.cc` | FASTA/FASTQ parsing/writing with gzip/bzip2 support |
| Utilities | `utils/` | CIGAR ops, threading wrappers, sequence comparison, maps |

**SIMD strategy:** Native implementations for SSE2/SSSE3 (x86_64), AltiVec/VSX (ppc64le), and Neon (ARMv8). Falls back to the SIMDe library for other architectures (riscv64, mips64el). Architecture detection happens in `configure.ac`.

**Threading model:** Global database state is read-only after initialization. Per-thread state handles computation. The library API uses mutex-protected session initialization.

**Build system:** GNU Autotools (`configure.ac` → `Makefile.am` → `src/Makefile.am`). `config.h.in` is generated at configure time and defines `HAVE_*` feature macros used throughout the source.
