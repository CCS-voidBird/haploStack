# haploStack I/O reference — v0.7.7

**A genotype file is read as text.**

This release changes one thing in the package: how a file's bytes become
values. Nothing about the `GenoMeta` container, its fields or its contract
moved. But the read step was making a decision it had no business making, so
this document sets out the whole path — **from the file on disk to the
`metadata` list inside the container** — and says exactly where each field
comes from.

Companion documents: [`geno_meta_structure.md`](geno_meta_structure.md) for the
container's anatomy, and
[`haploStack_io_reference_0.6.0.md`](haploStack_io_reference_0.6.0.md) for the
missing-allele and centring contracts, which are unchanged.

---

## 1. What was wrong

`utils::read.table()` runs `type.convert()` on every column it reads.
`type.convert()` tries, in order, logical → integer → double → character, and
it accepts `T` and `F` as logical literals. So a column whose every value is
`"T"` is returned as `TRUE`.

In a nucleotide genotype file, that column is an allele.

```
PED on disk           read.table()          what the package then saw
-----------           ------------          -------------------------
T T A T               logical  TRUE          "TRUE" "TRUE" ...
T T T T      ──────►  logical  TRUE   ──────► mode(x) <- "character"
T T A A               character "A"/"T"       ... and the read aborts
```

Three readers were affected, in three different ways:

| Reader | What happened | How it showed up |
|---|---|---|
| `read_plink_ped()` | `"TRUE"` is four characters; `.reshape_wide_alleles()` requires one | **Hard failure.** *"cells contain multi-character values. Wide format requires exactly one character per cell."* The file could not be read at all. |
| `read_effect_file()` | the `Allele` column became `TRUE`, which matches no genotype | **Wrong numbers.** `.compute_effect_size()` scores an absent allele as 0, so every marker contributed nothing and GEBVs came out flat. Nothing warns at read time. |
| `read_dart_genotype()` | `as.numeric("TRUE")` is `NA`, so the column fails the 0/1/2 sample-column test | **Mis-parse.** `first_sample` shifts and a sample column is reclassified as marker metadata. |

### It is the allele *column*, not the marker

This is the detail that kept the defect hidden. A PED stores `ploidy` columns
per marker. For an A/T SNP whose **first copy** happens to be fixed `T` across
the panel, column `2k−1` is uniformly `"T"` while the marker itself still
segregates perfectly well on column `2k`.

On the GRDC barley AGG1 panel (`data/GRDC/hordeum_vulgare/GRDC_barley_AGG1/
plink/geno.ped`, 1107 lines × 9856 markers):

| | |
|---|---|
| PED allele columns | 19 712 |
| Columns `type.convert()` returns as logical | **511** |
| Distinct markers involved | 511 of 9 856 |
| Markers that are monomorphic-T | **0** |

And the count is a property of the **panel**, not of the file format: a single
missing call (`0`) anywhere in a column is enough to keep it character. Subset
the lines and the number moves in either direction. That is why every fixture
in the test suite happened to escape it, and why a file that reads today can
stop reading tomorrow when you add a line.

---

## 2. The rule now

> Every column of a genotype file arrives as `character`. Each reader then
> coerces, by name and on purpose, only the columns it knows to be numeric.

The shared primitive is `.read_table_chr()` (internal, `R/format_helper.R`):

```r
.read_table_chr(file, header = FALSE, sep = "auto", strict = FALSE)
```

- `data.table::fread(colClasses = "character")` — no inference, and far faster
  on a wide file than `read.table()`.
- Explicit errors for a missing file and for an empty one.
- `""` is normalised to `NA` (no allele vocabulary contains the empty string).
- **Ragged input is checked here, not delegated to `fread`.** See below.

### Why raggedness is checked in R rather than left to `fread`

`fread(fill = FALSE)` does **not** reliably error on a ragged file. Given a PED
whose first line is wider than the rest, it settles on the narrower width and
**silently drops the wider row** — returning a table with one fewer sample, no
warning, and a container that looks perfectly valid:

```
F1 S1 0 0 0 -9 A T C G      <- 10 fields, DISCARDED
F2 S2 0 0 0 -9 A T          <-  8 fields, kept
```

read as one sample, two markers, **ploidy 1**. Silently losing a line is worse
than the type-inference defect this release exists to fix, so the fread call
always passes `fill = TRUE` (short rows padded, no row ever discarded) and the
reader checks the shape itself:

| Check | Applies | Behaviour |
|---|---|---|
| every data line appears as a row | always | error naming both counts |
| the table is rectangular (no padded cell in the last column) | `strict = TRUE` | error naming the first short line |

`strict = TRUE` is used for the **PED and the MAP**, whose lines are
fixed-width by construction: a short one is a truncated file, and padding it
into missing calls would hide the damage. `strict = FALSE` keeps the lenient
`read.table(fill = TRUE)` behaviour that the effect-file and DArT readers
always had.

`.read_block_file()` deliberately still uses `read.table()`. A block file holds
identifiers and coordinates, never allele codes; `Start`/`End` are *wanted* as
integers, and those columns cannot be uniformly T/F.

---

## 3. From input format to metadata — `read_plink_ped()`

This is the full path. Every arrow is one step in the function.

```
geno.map                              geno.ped
(chr, snp_id, cM, bp)                 (FID IID PAT MAT SEX PHENO | a1 a2 a1 a2 …)
   │                                     │
   │ .read_table_chr(header = FALSE)     │ .read_table_chr(header = FALSE)
   │   every column character            │   every column character  ← THE FIX
   ▼                                     ▼
 map_df                               ped
 Chrom  = as.character(col 1)         n_allele_cols = ncol(ped) - 6
 Marker = as.character(col 2)         ploidy        = n_allele_cols %/% n_markers
 Pos    = as.numeric(col 4 or 3)  ←── explicit coercion, the only one
   │                                     │
   │ .process_map(map_preset = "PL",     │ allele_mat = as.matrix(ped[, 7:ncol])
   │              pos_unit)              │ rownames   = ped[[2]]  (IID, verbatim)
   ▼                                     ▼
 map_info                             .reshape_wide_alleles(mat, n_markers, ploidy)
 Chrom | Marker | Pos (bp)              "T" + "A"  →  "TA"      one cell per marker
   │                                     │
   │                                     ▼
   │                                   split_haplotypes(split_char = "", ploidy)
   │                                     │
   │                                     ▼
   └──────────────────────────────►  geno_arrays   [n_samples × n_markers × ploidy]
                                         │
                                         ▼
                                   new_geno_meta(..., extra_na = "0")
                                         │
                       normalise_missing_alleles() runs HERE, at the one portal
                       — PLINK's "0" becomes NA, per allele (0.6.0 contract)
                                         ▼
                                     GenoMeta
```

### Where each `metadata` field comes from

| `metadata` field | Derived from | Note |
|---|---|---|
| `ploidy` | `(ncol(ped) − 6) %/% nrow(map_df)` | detected, never assumed; errors if not divisible |
| `n_samples` | `nrow(ped)` | |
| `n_markers` | `nrow(map_df)` — the **MAP** decides | the PED must agree or ploidy detection fails |
| `char_dict` | the distinct observed alleles, after `"0"` and `NA` are excluded | on the barley panel: `{T, A}` |
| `encoding` | `"nucleotide"` if any observed allele matches `[ACGTacgt]`, else `"selectiontools"` | **declared, not detected** — this is why a 1/2-coded PED is not mistaken for dosage, where `"0"` would be a real genotype |
| `source_format` | the literal `"plink_ped"` | records *how it got here*, distinct from `encoding` |

And beside `metadata`:

| Container field | From |
|---|---|
| `geno_arrays` | the allele columns, reshaped and split |
| `map_info` | MAP columns 1, 2 and 4 (or 3), `Pos` in bp |
| `sample_info` | PED columns 1–5 → `Family, Sample, Parent1, Parent2, Sex` |
| `block_info` / `block_marker_map_info` | `block_file`, when supplied |

### What changed in this table

Only one row, and it is `char_dict` / the array contents: alleles are now the
letters the file spells. Everything else is computed exactly as before.

Plus one deliberate behaviour change that is worth knowing about:

> **Identifiers are no longer type-converted.** A sample id written `007` stays
> `"007"`. It used to arrive as the integer `7` and become `"7"`, which silently
> broke the join back to the source file. An id that reads as a boolean stays
> `"TRUE"`.

If you have code that joined haploStack sample ids against a zero-padded id
column by accident-of-coercion, that join now works — and if you compensated
for the old behaviour by stripping zeros yourself, remove the compensation.

---

## 4. The other two readers

### `read_effect_file(path, sep = " ", header = FALSE)`

```
file: Marker  Allele  Effect  [extra numeric columns …]
        │       │       │
        │       │       └── coerced with as.numeric(), explicitly, cols 3+
        │       └────────── stays character  ← THE FIX
        └────────────────── stays character (so "0100" keeps its spelling)
                                  │
                                  ▼
                          transform_effects()  → Marker | Allele | Effect
```

The allele column is the one that must never be inferred. Column 1 staying
character is a bonus with the same motive as the sample ids above.

### `read_dart_genotype(file, …)`

Read as text with `fill = TRUE` and `sep` chosen from the extension. Nothing
downstream had to change, because the function already re-derived every number
it needed with `as.numeric(as.character(.))` — the dosages, the positions, and
the metadata/sample column split. The fix simply removes the case where that
re-derivation was handed `"TRUE"` and returned `NA`.

---

## 5. Migration checklist

- [ ] **Nothing to do for a container you already saved.** This is a read-path
      change only; `.rds` files on disk are unaffected and no repair is needed.
- [ ] **If a PED previously refused to load** with *"cells contain
      multi-character values"*, it loads now. No workaround required —
      `ready_to_use/workshop/_common.R`'s `read_plink_ped_safe()` has been
      withdrawn, and module 01c calls `read_plink_ped()` directly.
- [ ] **If you use a single-allele effect file**, recheck any result built from
      it under ≤ 0.7.6. Symptom: GEBVs that are all zero, or a
      `.validate_effect_encoding()` warning naming `TRUE` as an unknown allele.
      Re-read the file and recompute.
- [ ] **If you match sample or marker ids across systems**, check for
      zero-padded ids. They keep their leading zeros now.
- [ ] **If you rely on ragged PED lines being padded**, they now error instead.
      That is intentional; fix the file. (Under ≤ 0.7.6 such a line was not
      padded either — depending on where it sat in the file it could be dropped
      outright, taking the whole sample with it.)

---

*haploStack 0.7.7, 2026-09-20. Suite: 2 584 pass, 0 fail, 0 error, 12
environmental skips. Verified on the GRDC barley AGG1 PED: the container
`read_plink_ped()` returns is `identical()` — `geno_arrays` and `map_info`
both — to the one the withdrawn workaround produced.*
