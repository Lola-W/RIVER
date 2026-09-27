# Synthetic RIVER demo

Run these commands from the repository root after following the [installation instructions](../README.md#installation):

```bash
bash examples/run_demo.sh
```

This single command generates inputs, runs RIVER with one preparation thread, and compares the results with [expected/](expected/). No study data, external reference download, network access, or random seed is needed after installation.

The default output is `examples/demo/` (ignored by Git). Existing directories are never overwritten. To rerun, choose a new directory:

```bash
bash examples/run_demo.sh /path/to/new-demo
```

Use `RIVER_PYTHON` and `RIVER_SAMTOOLS` for executables outside the active environment:

```bash
RIVER_PYTHON=/path/to/python RIVER_SAMTOOLS=/path/to/samtools \
  bash examples/run_demo.sh /path/to/new-demo
```

On a Slurm cluster, request an interactive compute allocation using your site's normal `salloc` or `srun` options. Inside the allocation, activate the `river` environment, change to the repository root, and run `bash examples/run_demo.sh`. The selected Python and SAMtools must be available inside that allocation. Queue wait is excluded from the reported runtime.

## Generated inputs

[generate_synthetic.py](generate_synthetic.py) follows the small-alignment approach in `tests/test_river.py`, but runs independently of the test suite. It creates:

- `data/reference.fa` and its FASTA index: **500 artificial A bases** on a toy contig named `chr1`. These are not GRCh38 sequence or real genomic loci.
- Four coordinate-sorted BAMs and their BAI indexes in `data/bams/`: bulk, raw ppmSeq, and two PTA cells; **45 alignments total**.
- `data/samples.tsv` and `data/variants.tsv`, copied from the checked-in example tables. BAM paths are relative to the generated sample table; use this copy when running RIVER.

Every alignment is 100 bases long, has mapping quality 60 and base quality 40, and covers either position 50 or 200. The bulk and retained ppmSeq reads have eight reference and two alternate observations at each covered position. Both cells have an alternate read at position 50 and a reference read at position 200. Position 400 has no coverage.

The ppmSeq BAM also contains an extra alternate read at position 50 with `st=MIXED, et=PLUS`. Preparation must discard it, retaining **20 of 21** raw reads. The variant table records the three alleles as bulk-origin calls and records position 50 again as a ppmSeq-origin call.

## Individual commands

The runner is a convenience wrapper around these commands. Use fresh input and output directories:

```bash
python examples/generate_synthetic.py --outdir examples/demo/data

river run \
  --variants examples/demo/data/variants.tsv \
  --samples examples/demo/data/samples.tsv \
  --reference examples/demo/data/reference.fa \
  --outdir examples/demo/results \
  --threads 1

python examples/check_demo.py examples/demo/results
```

The generator uses pysam to write and index its inputs. RIVER uses the external SAMtools executable for preparation and pileup.

## Expected results

The full reference tables are committed in [expected/](expected/):

| File | Data rows | What to check |
|---|---:|---|
| [read_support.tsv](expected/read_support.tsv) | 12 | Three unique alleles × four queried BAMs |
| [variant_support.tsv](expected/variant_support.tsv) | 4 | Three bulk origins plus the repeated ppmSeq origin |
| [summary.tsv](expected/summary.tsv) | 2 | One row per source sample/method/callset |

Expected read counts below are **per dataset**; both cells are queried separately:

| Position | Dataset(s) | REF | ALT | Depth | VAF | Qualifying support |
|---|---|---:|---:|---:|---:|---|
| 50 | Bulk and ppmSeq, each | 8 | 2 | 10 | 0.2 | 1 |
| 50 | Cell 01 and cell 02, each | 0 | 1 | 1 | 1.0 | 1 |
| 200 | Bulk and ppmSeq, each | 8 | 2 | 10 | 0.2 | 1 |
| 200 | Cell 01 and cell 02, each | 1 | 0 | 1 | 0.0 | 0 |
| 400 | All four datasets, each | 0 | 0 | 0 | NA | 0 |

The zero-depth rows have `coverage_status=no_ref_alt_depth`, `query_status=no_pileup_observations`, and undefined VAF/confidence bounds (`NA`).

| Position | Origin | Other supporting methods | Number | Class |
|---|---|---|---:|---|
| 50 | `bulk_library` (IlluminaSeq) | PTA-Seq, ppmSeq | 2 | >=2 |
| 200 | `bulk_library` (IlluminaSeq) | ppmSeq | 1 | 1 |
| 400 | `bulk_library` (IlluminaSeq) | none | 0 | 0 |
| 50 | `duplex_library` (ppmSeq) | IlluminaSeq, PTA-Seq | 2 | >=2 |

Two supporting cells count as **one PTA-Seq method**. The entire originating method is excluded from each row's other-method count. The uncovered position remains in the bulk denominator.

| Source | Candidates | Class 0 | Class 1 | Class >=2 | Proportions (0, 1, >=2) |
|---|---:|---:|---:|---:|---|
| `bulk_library` | 3 | 1 | 1 | 1 | 1/3, 1/3, 1/3 |
| `duplex_library` | 1 | 0 | 0 | 1 | 0, 0, 1 |

[check_demo.py](check_demo.py) compares all columns and row order with the expected TSVs, allowing `1e-12` absolute/relative tolerance for floating-point fields. It also checks `run.json` for completion and ppmSeq preparation counts. A successful check prints:

```text
Demo matches expected tables: 12 read rows, 4 variant rows, 2 summary rows.
```

`run.json` contains paths, timestamps, versions, and preparation details specific to each run, so it is not compared byte for byte. The output directory also contains intermediate candidate positions and the selected ppmSeq BAM.

## Runtime and resources

Allow **under 10 seconds** after installation for generation, RIVER, and comparison on one CPU. The demo needs no GPU; a 1 GB RAM allocation is sufficient, and generated inputs/results occupy less than 1 MB. Environment installation and scheduler queue time are separate. These estimates apply only to this tiny synthetic dataset, not full-depth study BAMs.

Measured on **2026-09-25** on a TSCC compute node: Rocky Linux 9.7 x86-64, Intel Xeon Gold 6240 at 2.60 GHz, one allocated CPU, Python 3.12.12, pysam 0.23.3, and external SAMtools 1.21 / HTSlib 1.24. The expected TSVs were generated with this environment and reviewed against the counts above.

| Operation | Wall time | Peak resident memory |
|---|---:|---:|
| Generate and index inputs | 0.16 s | 28 MiB |
| Run RIVER | 0.26 s | 32 MiB |
| Complete runner, including verification | 0.49 s | 32 MiB |

Times were measured with `/usr/bin/time`, exclude Slurm launch/queue time, and are single-run observations. Inputs and results used about **0.5 MB** of allocated disk on the shared filesystem; the installed environment used about **362 MB**, excluding the package cache. The complete 21-test suite passed in about **3 seconds**. The demo also matched the same expected tables using SAMtools 1.18 / HTSlib 1.18 in the existing Python 3.12.12 / pysam 0.23.3 environment.
