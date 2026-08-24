# haploStack I/O reference — v0.6.0

A user guide to the **two contracts** that changed in **haploStack 0.6.0**, both
of which move numbers you may already have published:

- **What a missing allele is.** Inside a `GenoMeta` it is `NA` and nothing else,
  normalised at one place on the way in.
- **What "centring" means.** One vocabulary with three named modes, replacing a
  logical, and the catalogue's default is no longer to centre.

New exports: `repair_geno_meta()`, `use_missing_chars()`, `missing_chars()`,
`reset_missing_chars()`.

Companion sheets:
[`haploStack_io_reference_0.5.5.md`](haploStack_io_reference_0.5.5.md),
[`haploStack_io_reference_0.5.3.md`](haploStack_io_reference_0.5.3.md),
[`haploStack_io_reference_0.5.2.md`](haploStack_io_reference_0.5.2.md),
[`haploStack_io_reference_0.4.7.md`](haploStack_io_reference_0.4.7.md).
Container spec: [`geno_meta_structure.md`](geno_meta_structure.md).

---

## 1. Missing alleles

### The rule

> **Inside a `GenoMeta`, a missing allele is `NA`. Nothing else — regardless of
> the source encoding.**

It is applied **per allele**, once, at `new_geno_meta()` — the one constructor
every route reaches. `build_geno_meta()` goes through it, and so do
`read_plink_ped()`, `read_dart_genotype()` and `read_vcf_genotype()`, which do
**not** go through `build_geno_meta()`. There is no way in that skips it.

| source format | characters read as missing |
|---|---|
| `nucleotide` | `""`, `"NA"`, `"N"`, `"n"`, `"."`, **`"-"`**, **`"0"`** |
| `selectiontools` | the same set — `"0"` is new in 0.6.0 |
| `numeric` / dosage | the same **minus `"0"`**, which is a real genotype there |

Because the rule is per allele, a diploid cell written `"00"`, `"--"`, `"0/0"`
or `"-/-"` becomes a missing **genotype** — the reader splits the cell first, so
both copies land on the rule. A half-called cell keeps what was observed:

```
"AT"  ->  c("A", "T")
"A0"  ->  c("A", NA)     # the A still counts toward allele frequencies
"00"  ->  c(NA,  NA)
```

### What moves in your output

On any panel with missing calls, read through any reader:

| Output | Before 0.6.0 | Now |
|---|---|---|
| marker effects (`fit_genomic_prediction()`) | the missing code could rank as an allele and take the effect | effects sit on real alleles only |
| local GEBVs at `center = "allele_freq"` | frequencies computed over a denominator that counted missing | correct denominator |
| block signatures | `A0A` | `A-A`, and responds to `options(haploStack.fill_na)` |
| `>2 alleles` warnings | fired on biallelic panels with missing calls | gone |

**Unchanged:** MAF, `Major`/`Minor`, call rate, block definitions, LD, GEBVs at
`center = "none"`, and every dosage container.

### Containers you already saved

A `.rds` written before 0.6.0 may still hold a source token as an allele.
Loading one tells you so and names the fix; it does **not** rewrite your file.

```r
meta <- load_geno_meta("panel.rds")
#> load_geno_meta: 'panel.rds' holds 1462492 allele cell(s) across 9759 marker(s)
#>   as a missing token ("-") rather than NA ... The file has NOT been modified.

meta <- repair_geno_meta(meta)
#> repair_geno_meta: normalised 1462492 missing allele cell(s) to NA.

save_geno_meta(meta, "panel.rds")     # your call, not the loader's
```

`is_geno_meta()` warns on such a container but still returns `TRUE` — it is
structurally valid, its content merely predates the invariant. For a pipeline
that wants a hard gate:

```r
options(haploStack.strict_missing = TRUE)   # the warning becomes an error
```

### When your file uses a character the package does not know

Register it for the session. `mode` is required — a character that means
"absent" in one vocabulary can be a genotype in another.

```r
use_missing_chars("#", mode = "nucleotide")
#> use_missing_chars: "#" now read as missing for 'nucleotide' (session only).
#>   Effective set: "" "NA" "N" "n" "." "-" "0" "#".

missing_chars("nucleotide")                 # the effective set
missing_chars("nucleotide", chars_only = TRUE)   # just what you registered
reset_missing_chars()                       # clear
```

One character only, and never a real allele of that mode — A/C/G/T for
`nucleotide`, the codes 1–4 for `selectiontools`, and **any digit** for
`numeric`. A dosage panel that genuinely writes `9` for missing declares it on
the read instead, where you own the risk:

```r
meta <- build_geno_meta(geno, map, encoding = "numeric", extra_na = "9")
```

`extra_na =` is on `build_geno_meta()`, `read_genotype()`, `read_st_genotype()`,
`read_gs_genotype()` and `new_geno_meta()`.

### It only goes one way

Once a container exists, the character the **source file** wrote is gone. You
can render missing back out for export, but that is a new rendering, not a
recovery — nothing can tell you the file said `"0"`.

```r
out <- recode_geno_meta(meta, to_type = "selectiontools", na_as = "-1")
#> recode_geno_meta: rendered 396750 missing allele(s) as "-1". The result is no
#>   longer a normalised container ...
```

Such a container is deliberately outside the invariant; `is_geno_meta()` says
so, and `repair_geno_meta()` returns it to `NA`.

`recode_geno_meta()` also declares both ends now: `from_type` defaults to the
recorded provenance but is **validated against the data**, and a dosage source
with an allele target is refused rather than pushed through the nucleotide
dictionary.

---

## 2. Centring

### The vocabulary

`center_by_freq` (logical) is replaced by `center` (three named modes). With
`e_a` an allele's effect, `p_a` its frequency, and
`mean_i = sum_a p_a * e_a` the locus mean of marker `i`:

| `center =` | contribution | allele contrast | block variance | use it for |
|---|---|---|---|---|
| `"none"` | `e_a` | `e_A - e_B` | reference | the raw quantity; **the new default** |
| `"locus_mean"` | `e_a - mean_i` | **preserved** | **unchanged** | selection, and any cross-package comparison |
| `"allele_freq"` | `e_a * (1 - p_a)` | becomes `q*e_A - p*e_B` | changed | reproducing output published with ≤ 0.5.5 |

`"locus_mean"` subtracts **one** value per marker from every allele at it — a
pure location shift. It is the only mode that is invariant to adding a constant
to both alleles' effects, and it reproduces HapSelect's `mean_adjust = TRUE`
**value for value**, because
`d*e_A + (ploidy - d)*e_B - ploidy*mean_i` is exactly
`(d - ploidy*p_A) * (e_A - e_B)` — which is HapSelect's `(d - mean_d) * beta`.

`"allele_freq"` centres each allele by **its own** frequency, so the two alleles
of a marker are weighted differently. That is a reweighting, not a shift: block
variances change and rankings can move with them.

### The default flipped

`run_haplotype_catalogue()` centred by allele frequency by default;
`compute_profile_localgebv()` did not. The two entry points to the same engine
now agree, at `"none"`.

```r
# 0.5.5 numbers, explicitly:
run_haplotype_catalogue(meta, effect_list = eff, center = "allele_freq")

# 0.6.0 default:
run_haplotype_catalogue(meta, effect_list = eff)          # center = "none"
```

**This changes published catalogue block variances** — a median factor of about
3–4 in measurement, varying block by block. If you have a table from an earlier
release, name its setting in the file rather than assuming it matches a new one.

### The old argument still works

`center_by_freq` is deprecated for one release, with a one-time message:

```r
center_by_freq = TRUE   ->  center = "allele_freq"
center_by_freq = FALSE  ->  center = "none"
```

Supplying both is an error. Every existing call keeps working and keeps its
numbers; only the unspecified default moved.

`center =` is available on `compute_profile_localgebv()`,
`compute_block_catalogue()`, `run_haplotype_catalogue()`,
`count_block_subtypes()`, `count_favourable_subtypes()`,
`count_favourable_blocks()`, `analyze_block_catalogue()` and `predict_gebv()`.
Only `run_haplotype_catalogue()`'s **default** changed; the selection functions
still default to `"allele_freq"`.

---

## 3. Two things worth knowing before you compare across packages

**Levels.** `level = "genotype"` sums over every ploidy copy;
`level = "haplotype"` is a single copy. Across levels the factor is between
`ploidy` and `ploidy^2` — a population property you cannot divide out. Compare a
genotype-level figure from another package against `level = "genotype"`.

**`level = "haplotype"` is a slot, not phase.** It reports the genotype matrix's
copies in the order the reader wrote them. An unphased source writes every
heterozygote in one fixed order, so slot 1 can be a whole-panel artefact rather
than a chromosome. Nothing in the package checks this and no warning fires.

---

## 4. Migration checklist

1. **Re-read, don't repair, where you can.** Rebuilding a container from its
   source file is cleaner than repairing a saved one — the reader normalises.
2. **Repair what you cannot rebuild**, with `repair_geno_meta()`, and
   `save_geno_meta()` it yourself.
3. **Name the centring in every output file** you keep. `"allele_freq"`
   reproduces your ≤ 0.5.5 catalogue numbers; `"locus_mean"` is what to use for
   a comparison against another package.
4. **Expect marker effects and centred GEBVs to move** on any panel with missing
   calls. That is the fix, not a regression — the old numbers had an effect
   estimated for missingness in them.
5. **Drop any `"00" -> NA` pre-processing of your own.** The package does it now,
   and doing it twice is harmless but no longer needed.
