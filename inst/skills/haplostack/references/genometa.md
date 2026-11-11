# The GenoMeta container

## Elements

| Element | Shape | Notes |
|---|---|---|
| `geno_arrays` | 3-D character array `[samples x markers x ploidy]`, dimnames = (sample IDs, marker IDs, copy labels) | allele labels (`"A"`, or digits such as PLINK's `"1"`/`"2"`); missing = `NA`. Not a list. |
| `map_info` | `Chrom, Marker, Pos` | in the array's marker order |
| `block_marker_map_info` | `Chrom, Marker, Block` or `NULL` | needed by catalogue and simulation |
| `block_info` | `Chrom, Block, Start, End, Pos, Num_SNPs, Block_length` or `NULL` | `Region_Start` / `Region_End` after a region map |
| `sample_info` | `Family, Sample, Parent1, Parent2, Sex` (PLINK FAM) | `Sample` = the genotype sample IDs; `"0"` = unknown / founder |
| `marker_effects` | named list, one `(Marker, Allele, Effect)` table per trait, or `NULL` | |
| `phenotype` | `Sample` + numeric trait columns, or `NULL` | |
| `metadata` | `ploidy, n_samples, n_markers, source_format, encoding, missing_tokens, char_dict, blocks, schema_version, built_with` (+ `family_history`, `signal_selection`, `annotation`) | counters are derived from the array |

`is_geno_meta(meta, verbose = TRUE)` checks the shape and says what is wrong;
`strict = TRUE` adds a full scan. After a hand edit (`meta$marker_effects <-
...`) run it again. It fails when the counters disagree with the array, when
`map_info` is not in array order, or when a marker sits in two blocks.

## What each transform keeps

| Transform | `sample_info` | `marker_effects` | `phenotype` |
|---|---|---|---|
| `filter_geno_meta()` | subset | pruned | pruned |
| `attach_block_info()` / LD resolvers | kept | kept | kept |
| `recode_geno_meta()` | kept | kept | kept |
| `attach_phenotype()` | kept | kept | set |
| `attach_family_history()` | `Family` written | kept | kept |

Re-blocking clears only what depended on the old blocks (signal selection,
gene annotation, block provenance). Offspring containers from the simulation
drop the founders' phenotype and family history.

## Building one

```r
meta <- read_plink_ped(ped, map, block_file = blocks)           # PLINK text PED/MAP
meta <- read_genotype(geno_file, map_file, preset = "PL", block_file = blocks)
meta <- build_geno_meta(geno_2d, map, fam = NULL, blocks = block_df)
```

- `build_geno_meta()` wants `geno_2d` with rows = samples, columns = markers,
  both named (no automatic transpose); cells like `"A/T"` (`split_char = "/"`).
- `read_genotype()` presets: `GENO_FORMAT_PRESETS` (`"PL"`, `"ST"`, `"GS"`,
  `"default"`); `resolve_geno_format()` shows what a preset means.
- Every text reader reads all columns as text first, so IDs keep their
  spelling and an all-`"T"` allele column stays `"T"`. Blank lines are
  skipped. A line of separators (`,,`, common at the end of spreadsheet CSVs)
  reads as a row of `NA`: drop it.
- The DArT and VCF readers convert dosage to `A/B` labels and filter on the
  way in; prefer PLINK or a delimited matrix.

## Saving and loading

`save_geno_meta(meta, path)` / `load_geno_meta(path)`. Load only files you
trust (an `.rds` can hold any R object). A container saved by an older
haploStack (schema 1, or SelectionTools allele codes) is upgraded in memory
with one message; `upgrade_geno_meta(x)` does the same for an object. A
container from a newer haploStack is refused (`hs_error_schema`).

## Missing data

Inside a container a missing allele is `NA`, normalised once when it is built.
Missing tokens: `""`, `"NA"`, `"N"`, `"n"`, `"."`, `"-"`, `"0"` - except in
`0/1/2` dosage, where `"0"` is a genotype. A half call such as `"A0"` keeps its
`A`. Declare an extra token per read with `missing_chars =`, or for the session
with `use_missing_chars("#", mode = "nucleotide")`. The container records the
vocabulary it was built with; a later registration does not change it.

## Encodings

`"nucleotide"` (A/C/G/T, or digit labels used as alleles) and `"numeric"`
(`0/1/2` dosage). `detect_genotype_encoding()` returns one of the two.
SelectionTools' integer codes are refused unless declared; use `st_convert()`.

## Centring (local GEBVs, catalogue, favourable blocks)

| `center` | contribution of allele a | use for |
|---|---|---|
| `"none"` (default) | `e_a` | plain local GEBVs |
| `"locus_mean"` | `e_a - sum(p * e)` per marker | selection; comparison with HapSelect `mean_adjust` |
| `"allele_freq"` | `e_a * (1 - p_a)` | reproducing catalogues from <= 0.5.5 |

`level = "genotype"` sums all copies; `"haplotype"` is one copy in reader
order, which is **not phase** on unphased data.
