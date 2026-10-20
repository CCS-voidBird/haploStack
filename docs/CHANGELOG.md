# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com), and the project adheres to
semantic versioning.

> **Reorganised 2026-07-21 (documentation-management pass, two edits).**
> *First edit:* a `## [0.4.7]` section was added below for the marker-stats /
> block-allele-sheet release (previously undocumented here), and the prior single
> `## [Unreleased]` block was relabelled as development history. *Second edit:*
> explicit `## [0.4.5]` (2026-07-07) and `## [0.4.6]` (2026-07-13) sections were
> **reconstructed from evidence** — diffing the archived
> `ready_to_use/haploStack_0.4.5.tar.gz` and `_0.4.6.tar.gz` builds (source dates,
> `NAMESPACE` export deltas, and file presence). The two 2026-07-07-dated blocks
> (encoding/fitness fixes, block-catalogue revision) were confirmed present in the
> 0.4.5 build and moved under `[0.4.5]`; the 0.4.6 export delta (simulation-data
> generators + helpers) became `[0.4.6]`. Existing entry **text** was not altered,
> only regrouped under version headings; remaining older content sits under
> `## [≤ 0.4.4]`. Boundaries **below 0.4.5 are not recoverable** from the sources,
> so that block stays coarse. Pre-edit original: verbatim at
> [`docs/CHANGELOG.md.orig`](CHANGELOG.md.orig); version lineage table in
> [`docs/project_record/INDEX.md`](project_record/INDEX.md) (Release history) and
> [`docs/project_record/threads/ready-to-use-workflow.md`](project_record/threads/ready-to-use-workflow.md).

> **Updated 2026-07-22.** Added the `[0.4.9]` alignment-plot-performance section
> (Rcpp `collapse_genotype_sorted_cpp`), and **reconstructed a `[0.4.8]` section
> from the tarball diff** (`ready_to_use/haploStack_0.4.7.tar.gz` vs `_0.4.8.tar.gz`):
> 0.4.8 was shipped without a changelog entry — no export change, an internal
> catalogue-output column-ordering pass. Pre-consolidation baseline still at
> [`docs/CHANGELOG.md.orig`](CHANGELOG.md.orig).

> **Updated 2026-07-28.** Added the `[0.4.11]` signal-marker-selection section,
> and **reconstructed a `[0.4.10]` section from the code** — 0.4.10 was in
> `DESCRIPTION` but had no changelog entry, thread, or tarball (the newest build
> is still `ready_to_use/haploStack_0.4.9.tar.gz`). Evidence for the
> reconstruction is cited inside that section. Pre-consolidation baseline still
> at [`docs/CHANGELOG.md.orig`](CHANGELOG.md.orig).

> **Updated 2026-07-29.** Added the `[0.5.0]` alignment-chromosome-sub-panels
> section. Note that 0.5.0 was previously pencilled in for the stacking-simulation
> rebuild (#013) and the genomicSimulation persistent session (GSPS); those are
> still in progress and are **re-targeted to 0.5.x** — see INDEX discrepancy 13.

> **Updated 2026-07-29 (second pass).** Added `[0.5.1]` — the encoding-vocabulary
> unification. It closes two long-standing INDEX discrepancies (12: the broken
> `@examples` that had been failing `R CMD check` since ≤ 0.4.9; 15: both items).

> **Updated 2026-08-13.** Second documentation-only pass over
> `docs/geno_meta_structure.md`, re-verifying it against the 0.5.5 source
> (`[Unreleased]` → *Documentation (2026-08-13)*). Three new boundary findings
> open as INDEX discrepancy **20**; no code changed and the version stays 0.5.5.
> The same pass added the source-only `ready_to_use/genoMetaBuild/` deliverable
> (`[Unreleased]` → *Deliverable (2026-08-13)*), which is verbatim copies of the
> container-layer sources — again no package code change.

> **Updated 2026-08-17.** Added `[Unreleased]` → *Distribution (2026-08-17)* —
> the GitHub release assets for 0.5.5 in `ready_to_use/release_v0.5.5/`, the
> first release to ship a Windows binary beside the portable source tarball.
> Appended to thread [#026](project_record/threads/026-distribution-build.md).
> No code changed; the version stays 0.5.5. The dev-tree vs published-repo
> `DESCRIPTION` identity split found during the build is open as INDEX
> discrepancy **21**.

> **Updated 2026-08-17 (second pass).** Added `[Unreleased]` → *Metadata
> (2026-08-17)* — the canonical `DESCRIPTION` identity corrected to match the
> published repo (author **Chensong Chen**, `chensong.chen@uq.edu.au`,
> `CCS-voidBird/haploStack`), with `man/haploStack-package.Rd` regenerated. This
> **closes the identity half of discrepancy 21**. No code changed; the version
> stays 0.5.5.

> **Updated 2026-08-19.** Added `[Unreleased]` → *Prototype (2026-08-19)* — the
> **HapSelect bridge** in `dev/hapselect_bridge/`, a GenoMeta entry point for
> the package that can replace the SelectionTools blocking route. New dev-tree
> files only; **no package code, no export, no dependency** (HapSelect is
> deliberately absent from `DESCRIPTION`), version stays 0.5.5. Thread
> [#028](project_record/threads/028-hapselect-bridge.md); opens INDEX
> discrepancy **22**.

> **Updated 2026-08-21.** Added `[Unreleased]` → *Analysis (2026-08-21)* — the
> first **real-panel** exercise of the HapSelect bridge (GRDC barley, 1107 x 9856),
> which found that the PLINK missing code can capture a marker's trait effect.
> No package code changed; version stays 0.5.5. Thread
> [#028](project_record/threads/028-hapselect-bridge.md); opens INDEX
> discrepancy **23**.

> **Updated 2026-08-21 (second pass).** Raised **MAN** —
> `dev/PRD_missing_allele_normalisation.md` + its upgrade plan, targeting
> **0.6.0**: one representation of a missing allele inside a `GenoMeta`.
> Scoped, not implemented; decisions open as **Q21**. Registered under
> `[Unreleased]` -> *Proposed*.

> **Updated 2026-08-24.** **MAN and CEN shipped as `[0.6.0]`** — all five MAN
> stages and the whole CEN track, in one day. The two `[Unreleased]` sections
> that tracked the work in progress (*Proposed*, 2026-08-21, and *In progress*,
> 2026-08-24) are struck through and kept as the spec-level and staged records;
> the release-facing account, with both migration notes, is under `[0.6.0]`.
> Thread [#029](project_record/threads/029-missing-allele-normalisation.md);
> tarball `ready_to_use/haploStack_0.6.0.tar.gz`, verified by installing into a
> throwaway library.

> **Updated 2026-09-15.** Added `[0.7.0]` — **GBS: Kemper's two modes as
> written, and a bounded penalty scale** (thread
> [#030](project_record/threads/030-gbs-kemper-modes.md), plan
> `dev/UPGRADE_PLAN_opt_kemper.md` §4). `gbs_method = "LOG"` becomes Kemper's
> `τ·log p`, the absent-state penalty becomes `gbs_penalty · s_t` with a bounded
> risk × hope scale, and the ideal genotype is founder-fixed by default —
> **changes published GBS numbers**; TS / OHS / OCS are byte-identical. Tarball
> `ready_to_use/haploStack_0.7.0.tar.gz`, verified by installing into a
> throwaway library. The OPT scheme and the KCT rename from the same plan are
> **not** in this release.

> **Updated 2026-09-15 (second pass).** Added `[0.7.1]` — **direct pair
> selection: the placeholder and its budget** (thread
> [#031](project_record/threads/031-direct-pair-selection.md), Q23 i / Q24).
> `selection_unit = c("parents", "pairs")` and `pair_budget = c("families",
> "parents")` are reserved on `run_stacking_simulation()`, the budget that ties
> a direct-pairing cycle to the parents-first traversal is resolved and
> reported, and the exported `pair_selection_budget()` gives the numbers. No
> selector is implemented — `"pairs"` stops at the first cycle. Default runs are
> byte-identical to 0.7.0. Tarball `ready_to_use/haploStack_0.7.1.tar.gz`.

> **Updated 2026-09-20.** Added `[0.7.7]` — **a genotype file is read as text**.
> `read_plink_ped()`, `read_effect_file()` and `read_dart_genotype()` route
> their file read through the new internal `.read_table_chr()`, so R's
> `type.convert()` can no longer turn an allele column into a logical. This
> **closes INDEX discrepancy 30**, which the same day's `[Unreleased]` →
> *Distribution and deliverables (2026-09-20)* note had recorded as found but
> not fixed; the workshop's `read_plink_ped_safe()` workaround is withdrawn.
> The read path from file to `metadata` is written up in
> [`docs/haploStack_io_reference_0.7.7.md`](haploStack_io_reference_0.7.7.md).

> **Updated 2026-09-30.** Added `[0.7.8]` — **an exported effect file reads
> back with the defaults**: `read_effect_file()` defaults become
> `sep = "auto", header = "auto"`, so the tab-separated, headed file
> `export_marker_effects()` writes round-trips as its roxygen always said
> (closes INDEX discrepancy 31, intake 01). Same release: the workshop's
> HapSelect route now blocks at the declared `LD_TOLERANCE` (discrepancy 32,
> intake 03), and a false "ragged file" error on blank lines is recorded as
> found, not fixed (discrepancy 33, intake 08).

> **Updated 2026-09-30 (second pass).** Added `[0.7.9]` — **SelectionTools LD
> blocks straight from a GenoMeta** (`define_haplotype_blocks_meta()`, thread
> #034, discrepancy 34, intake 04; `define_haplotype_blocks()`'s chromosome-name
> join, rule check and working-directory side files fixed) and **a blank line
> is not a record** (`.read_table_chr()`, discrepancy 33 closed, intake 08).
> Third pass, same day: the first 0.7.9 build's blank-line fix refused a file
> with `,,` lines — which 0.7.8 read — as ragged; reworked and the assets
> rebuilt before release (see `[0.7.9]` *Fixed*). Fourth pass: the second build
> failed on a script re-run after `rm(list = ls())` (SelectionTools' global
> `st.input.dir` / `st.output.dir` gone); now supplied for the call, and the
> assets rebuilt again (see `[0.7.9]` *Added*).

> **Updated 2026-09-30 (evening).** Added `[0.7.10]` — **a block map first; LD
> blocking is a named side route** (thread #035): `hapS_ld_blocks()` (HapSelect,
> promoted from the dev bridge, the default resolver), `st_ld_blocks()` /
> `st_ld_blocks_files()` (SelectionTools, experimental); the 0.7.9 names
> deprecated.

## [Unreleased]

_The remaining roadmap items (stacking-simulation rebuild #013, genomicSimulation
persistent session GSPS) are re-targeted to **v0.5.x** and tracked in
[`docs/project_record/INDEX.md`](project_record/INDEX.md). An opening PRD for the
remaining input formats (allelic DArT, VCF) is at
`dev/PRD_input_formats_dart_vcf.md` — questions only, nothing scoped._

### Distribution and documentation (2026-09-20, second pass) — no package code changed

Version stays **0.7.7**; nothing in `R/`, `src/` or `NAMESPACE` moved.

**Published.** The GitHub mirror `ready_to_use/haploStack` now carries one
commit and one annotated tag per version — `0a6a2f5` / `haploStack-0.7.6` and
`4dda05a` / `haploStack-0.7.7` — pushed to `CCS-voidBird/haploStack`, and both
versions are published as GitHub releases from the `ready_to_use/release_v0.7.x/`
bundles (source tarball, Windows binary, `SHA256SUMS.txt`, release notes);
0.7.7 is marked *latest*. Tags 0.7.3–0.7.5, previously local only, went up with
the same push.

**Workshop.** `ready_to_use/workshop/_common.R` now stops at start-up when
`packageVersion("haploStack") < "0.7.7"`, with the reason, so the PLINK route's
requirement is enforced rather than merely stated. Header bumped to v0.7.7.

**Record fixes.** `test-ped-text-read.R` assertion count corrected to 47 in
INDEX discrepancy 30 (was 40); `release_v0.7.6/RELEASE_NOTES.md` known-issue
section rewritten in the past tense now that 0.7.7 exists; the
`block_annotation` bundle and `block_annotation_0.3.0.zip` carry the 0.7.7-pinned
bridge pointer. Details in thread
[#032](project_record/threads/032-genotype-text-read.md).

### Distribution and deliverables (2026-09-20) — no package code changed

Version stays **0.7.6**; nothing in `R/`, `src/` or `NAMESPACE` moved.

**Released.** `ready_to_use/release_v0.7.6/` — the portable source tarball, a
Windows binary built from it (`R CMD INSTALL --build`, R 4.4.3 ucrt + Rtools44),
`SHA256SUMS.txt` and `RELEASE_NOTES.md`. Verified by installing the binary into
the build machine's library and running the package's own suite **against the
installed build**: 2 537 pass, 0 fail, 0 error, 12 environmental skips (9 ×
`read_vcf_genotype` with `vcfR` absent, 1 × the monitor's `future.apply` path,
2 guards that skip because `circlize` and `genomicSimulation` *are* installed).
The GitHub mirror `ready_to_use/haploStack/` is synced to 0.7.6.

**Workshop extended to annotation.** `ready_to_use/workshop/` gains a single
path block, `_paths.R` (every file path, subset control and toolchain location,
each with an `HS_*` environment override, plus `check_paths()`), and two
modules: **`01c_build_geno_meta_PLINK.R`** (raw PLINK `.ped`/`.map` → GenoMeta →
HapSelect LD blocks → rrBLUP effects → block coordinate files) and
**`07_block_annotation.R`** (the three handover files → `block_annotation
0.3.0` offline against a local GFF3 → joined back on `Block`).
`00_run_all.R` takes `BUILD = "SIM" | "REAL" | "PLINK"`. Run end to end on the
GRDC barley AGG1 panel (1H + 2H, 300 lines, 2 677 markers): 1 774 blocks, 228
with at least one overlapping gene, 3 111 gene overlaps, 86 marker–gene
co-localisations, ~8.5 min. Artefacts kept in
`ready_to_use/workshop_products_0.7.7/`.

**Bridge rewritten.** `docs/haplostack_annotation_bridge.md` is pinned to
haploStack 0.7.6 / block_annotation 0.3.0 and now documents the raw-genotype
entry point, the offline local-GFF3 route (including GO / InterPro through an
xref dump), and the asymmetry that made the offline route silently mismatch:
`annotate_blocks_ensembl()` falls back to `SPECIES_CONFIG[species]["chr_map"]`
while `annotate_blocks_local()` does **not**, so barley blocks must be
normalised to the plain integer and the species `chr_map` passed explicitly.

**Found, not fixed.** `read_plink_ped()` cannot read a nucleotide PED whose
panel contains a monomorphic-`T` allele column: `utils::read.table()` returns
that column as **logical**, and `.reshape_wide_alleles()` rejects the resulting
`"TRUE"`. 511 of 19 712 allele columns of the GRDC barley PED are affected.
Proposed fix for 0.7.7 is `colClasses = "character"` on that call; the workshop
ships `read_plink_ped_safe()` as a documented workaround. INDEX discrepancy
**30**.

### ~~In progress (2026-08-24)~~ — **SHIPPED as [0.6.0]** the same day

> Kept as the staged record of how the release was built: stages 1 and 2 landed
> first, additive and numerically inert, and the inventory below is what they
> produced. Stages 3-5 and the CEN track followed on the same day — see
> **[0.6.0]** for the release-facing account, including the migration notes.

#### MAN stages 1-2, additive, no numeric change

First code for **MAN** (`dev/PRD_missing_allele_normalisation.md`, plan stages 1
and 2). Both stages are additive by design, so **no numeric output moves yet** —
the behaviour change is stage 3. Version stays **0.5.5**; the bump to 0.6.0 is
stage 5. Thread [#029](project_record/threads/029-missing-allele-normalisation.md).

**Added**

- `normalise_missing_alleles(arr, source_format, extra_na = NULL)` (internal,
  `R/format_helper.R`) — **FR1**, the primitive that maps a source format's
  missing tokens to `NA_character_` and returns the array plus a count of cells
  changed. The argument is the format being *parsed*, supplied by the caller,
  not a value read back out of a container (**Q21.d**). Nothing calls it yet.
- `repair_geno_meta(meta, verbose = TRUE, source_format = NULL)` — **FR5**,
  **exported**. The migration path for containers already on disk. Explicit and
  opt-in; it never touches the file. `source_format` overrides the recorded
  provenance for a container whose `metadata$encoding` is absent or wrong.
- `options(haploStack.strict_missing = TRUE)` — **FR4/Q21.b**, hidden option
  that turns the new warning into an error for pipelines wanting a hard gate.
- `use_missing_chars()`, `missing_chars()`, `reset_missing_chars()` (**FR2b**,
  all **exported**, `R/alphabet.R`) — see the amendment below.
- `tests/testthat/test-missing-normalisation.R` — 127 assertions covering **V2**
  (a dosage source keeps its `0`), **V8** (a `selectiontools` source loses it),
  **V5**, **V9**, **V10**, **V18-V20** (FR2b), and **`read_plink_ped()`'s first
  tests ever** (dimensions, ploidy detection, FAM columns, MAP parsing, error
  paths).

**Changed**

- `.missing_tokens()` gains `"0"` for `selectiontools` — **FR2/Q21.a**. A
  1/2/3/4-coded PLINK file read as ST wrote `"0"` for a missing call and kept it
  as a real allele code. `"0"` is still **never** missing for `numeric`/dosage
  data, where it is a real genotype; the primitive enforces that even against an
  explicit `extra_na`, because that is the one failure worse than the bug.
- `is_geno_meta()` now checks one **content** invariant beside the structural
  ones — a missing allele is `NA` and nothing else (**FR4**). A container
  carrying a source token stays *valid* (still returns `TRUE`) but warns, naming
  the token, the affected marker count and `repair_geno_meta()`. The warning
  carries the condition class `haploStack_missing_token`.
- `load_geno_meta()` **mentions, never repairs** (**FR6/Q21.c**): loading a
  pre-0.6.0 container reports what it holds, names `repair_geno_meta()`, and
  leaves the file on disk byte-identical. Auto-repair would silently change what
  the user saved.

**Migration inventory** (stage 2 gate) —
`dev/missing_normalisation_check/inventory.tsv` + `README.md`, produced by
`inventory_scan.R`. 26 `.rds` files scanned, 12 hold a container, **4 need
repair** — and the finding is that **the token in the wild is `"-"`, not
`"0"`**: `GRDC_project/AGG_CV_geno_meta.rds` (1,462,492 cells, 9,759 of 10,128
markers), `output/geno_meta.rds` and `output/geno_meta_0.5_t0.rds` (1,365,084
cells, 11,312 of 11,386 markers each), `catalogue_output/geno_meta_0.5_t0.rds`
(49,780 cells). The four barley/HapSelect containers are clean **because stage
01's workaround already removes PLINK's `"0"`** — which is the evidence V6 will
retire. `repair_geno_meta()` handles `"-"` with no source-format argument.

**Amended the same day — FR2b, the missing vocabulary stated and made
extensible** (PRD decisions **j-n**). Still additive; suite still identical to
baseline.

- **The default rule is now written down.** For a `nucleotide` source **both**
  `"0"` (PLINK's code) and `"-"` are missing, and the rule is **per allele** —
  so a cell written `"00"`, `"--"`, `"0/0"` or `"-/-"` is a missing *genotype*
  by consequence, while a half-called cell (`"A0"`) keeps its observed allele
  and marks only the absent copy. Neither token is new; the statement is.
  Stated in `?normalise_missing_alleles`, `?use_missing_chars` and beside
  `.missing_tokens()`.
- **`use_missing_chars(chars, mode)`** (**exported**) registers an extra
  missing character for the session — `use_missing_chars("#", mode =
  "nucleotide")` — with `missing_chars()` to inspect the effective set and
  `reset_missing_chars()` to clear. It is read inside `.missing_tokens()`, so
  one registration reaches the readers, the `is_geno_meta()` check,
  `repair_geno_meta()` and the legacy frequency/IBS/clustering paths at once.
  Named for `use_sample_order()`, the package's existing "validate, then set a
  global option" precedent; the option is `haploStack.missing_chars`.
- **`mode` is required, and a registration can never be a real allele of it.**
  There is no package-wide default because a character that means "absent" in
  one vocabulary is a genotype in another. Refused: anything longer than one
  character (missing is decided per allele), blank/whitespace (already missing
  everywhere), A/C/G/T for `nucleotide`, the codes 1-4 for `selectiontools`, and
  **any digit** for `numeric` — a dosage panel that genuinely writes `9` for
  missing declares it on the read (`extra_na = "9"`), where the caller owns the
  risk. A malformed option cannot change what is missing.
- **`metadata$missing_tokens`** — a container now records the vocabulary it was
  normalised against, so it stays self-describing in a session where nothing is
  registered. Written by `repair_geno_meta()`; the constructors and readers will
  write it too once they normalise (stage 3). Detection and repair take the
  **union** of the stamp and the current set — more sensitive is the safe
  direction for a detector — and the dosage exclusion survives the union.

Gates **V18-V20**; `tests/testthat/test-missing-normalisation.R` is now 127
assertions.

**Verification.** Full suite green with the skip and warning list **identical to
the pre-change baseline**; `run_examples()` exit 0. No existing test fixture
carries a missing token, so the stage-3 no-op claim has a clean starting point.

### Investigation (2026-08-21, no code change)

**Why haploStack and HapSelect report different block effects** —
[`docs/block_effect_semantics_haploStack_vs_hapselect.md`](block_effect_semantics_haploStack_vs_hapselect.md),
raised as a follow-up to **Q21.e**. Neither package is wrong, and the recorded
suspicion (genotype- vs haplotype-based computation) is **not** the cause — both
figures are genotype-level.

- The ~2.5x gap is **`center_by_freq`**: the catalogue centres by population
  allele frequency (its default), the HapSelect run did not. Ratio 2.5197 against
  the centred figure, 0.9807 against the uncentred one.
- The residual ~2% is the **missing-data policy** (HapSelect mean-imputes,
  haploStack scores a missing allele as effect 0). Proven rather than argued: on
  the **737 blocks with no missing data the ratio is 1.000000 exactly**, and
  across all 7604 blocks it correlates -0.383 with the block's missing rate.
- Both tables are internally exact against their own definitions (max abs
  difference 0 and 4.34e-19), so this is a reporting problem — the two published
  tables were never like-for-like and neither names its setting.
- **Side finding (upstream, HapSelect):** `Num_Uniq_Hap` counts missing as a
  haplotype state — `Num_Uniq_Hap == dosage classes + 1 when the block has
  missing` holds for **6181 of 6181** single-marker blocks, so 4630 biallelic
  SNPs report four haplotypes. It feeds `Unique_Haplo_Block_Var`, which
  `select_top_blocks()` can rank on; `Block_Var` is unaffected.

Explicitly **not** a MAN requirement: normalisation makes *what* is missing
unambiguous, but the *policy* for missing stays a per-function choice on both
sides.

**Amended the same day** with three additions:

- **Level semantics.** `.assemble_localgebv_wide()` makes `level = "genotype"`
  the sum over **every ploidy copy** (2 strands diploid, 4 tetraploid) and
  `level = "haplotype"` a **single strand whatever the ploidy**. Verified:
  `max |genotype - (hap1 + hap2)| = 0` over 500 blocks x 1107 samples. So
  HapSelect's `Block_Var` (genotype-level) may only be compared against
  `level = "genotype"`; across levels the factor is between x ploidy and
  x ploidy^2, a population property that cannot be divided out.
- **`mean_adjust` is a location shift, not a mirror of `center_by_freq`.**
  Both of its branches mean-impute, so they differ by the constant
  `sum(mean * beta)` — identical `var_haplo_effects` to 12 significant figures,
  difference sd 1.3e-16. **HapSelect has no setting that reaches
  `center_by_freq`**, so alignment has exactly one direction:
  `center_by_freq = FALSE` on the haploStack side. Adopted for the re-run
  (plan stage 5, V6).
- **Caveat on unphased input:** the barley PED writes every heterozygote as `AT`
  and never `TA` (7,584,722 vs 0), so strand 1 is the A-strand panel-wide and
  `level = "haplotype"` reports field order rather than phase. Nothing in
  haploStack checks this and no warning fires.

### ~~Proposed (2026-08-21)~~ — MAN + CEN, **shipped as [0.6.0]** on 2026-08-24

> Kept as the spec-level record of what was proposed and why, including the
> solution selection and the four measured consequences. Everything below was
> implemented; the release-facing account is under **[0.6.0]**.

**One representation of missing.** `dev/PRD_missing_allele_normalisation.md`
(plan: `dev/UPGRADE_PLAN_missing_allele_normalisation.md`; decisions: **Q21**)
proposes that inside a `GenoMeta` a missing allele is **`NA` and nothing else**,
normalised at construction where the encoding is known, validated by a
`is_geno_meta()` warning, with an exported `repair_geno_meta()` for containers
already on disk. Closes INDEX discrepancy **23** and the effect half of **18**.

The container currently permits three representations — `NA`, a reader token
such as PLINK's `"0"`, and the `fill_na` placeholder `"-"` — and three different
idioms decide which is which (`na_vals`, `strip_blank`,
`!is.na & nzchar(trimws())`). Four consequences were measured, not inferred:
marker effects estimated for the missing code (8751 of 9856 markers mis-coded on
the barley panel, 34 outright); **catalogue local GEBVs shifted under
`center_by_freq = TRUE`, the default** (max abs difference 0.318 on a fixture;
exactly 0 with `center_by_freq = FALSE`, isolating the frequency denominator as
the mechanism); signature strings labelled `A0A` instead of `A-A`, ignoring
`options(haploStack.fill_na)`; and a container that reports call rate 0.70 while
listing `0` as an allele.

**The solution was selected, not assumed.** Patching the four offending call
sites is rejected as insufficient — it leaves two legal representations of
missing and does not reach the catalogue filter or the signature layer — and a
central read seam is rejected because `geno_arrays` is documented public, so
user code reading it directly would still see the token. Normalisation must be
**encoding-aware**, because `0` is missing for nucleotide but a **real genotype**
for a dosage container, which is what fixes its home at the constructor and the
readers. **0.6.0 rather than 0.5.6** because published marker effects and
default-setting catalogue GEBVs change; LD, blocks, MAF, call rate and
`center_by_freq = FALSE` GEBVs are verified unchanged.

**Second workstream in the same release — centring semantics (CEN).** Decisions:
**Q22**; same PRD, *Second workstream*. haploStack's only centring
(`center_by_freq = TRUE`, the catalogue default) weights each allele by **its
own** frequency, which moves the A-vs-B contrast from `e_A - e_B` to
`q*e_A - p*e_B` and changes block variance — median factor **3.03** on the barley
panel. HapSelect's `mean_adjust` is the other kind, a pure location shift that
cannot change variance at all, so neither package could reach the other's
setting.

- **The catalogue's centring default flips to `"none"`**, matching
  `compute_profile_localgebv()`, which already defaulted that way. **Changes
  published catalogue block variances**; `center = "allele_freq"` reproduces the
  old numbers.
- **A new `center = "locus_mean"` mode**: every allele at a locus is centred by
  the *shared* locus mean `sum_a AlleleFreq(a) * Effect(a)`. It preserves the
  allele contrast, leaves block variance identical to uncentred (a per-block
  constant `ploidy * sum_i mean_i`), and reproduces HapSelect's
  `mean_adjust = TRUE` **value for value** — verified **6.94e-18** max absolute
  difference across 300 barley blocks, variance difference **2.71e-20**. It is
  also invariant to adding a constant to both alleles' effects, which is the
  property 0.5.3 wanted from `center_by_freq` and did not get.
- `center_by_freq` is retained as a deprecated logical
  (`TRUE -> "allele_freq"`, `FALSE -> "none"`).
- **Statements, not behaviour changes:** `level = "haplotype"` is a slot
  assignment in reader order, not phase — correct behaviour for an unordered
  genotype matrix, now to be documented (INDEX discrepancy 24 **demoted** to a
  docs requirement). Dosage sources get told they carry no strand resolution, and
  pointed at HapSelect for LD/blocking.
- **Found while scoping that:** INDEX discrepancy **25** — a dosage container is
  built ploidy-2 with an entirely `NA` second layer, so `compute_marker_stats()`
  reports `Call_rate` **0.5** on every marker and `detect_genotype_encoding()`
  errors on the collapsed matrix. Documented, **not** fixed in 0.6.0: the fix
  changes `Call_rate` on every dosage panel and needs its own spec.

### Analysis (2026-08-21, no code change)

**GRDC barley run — the HapSelect bridge on a real panel.** Story:
[`docs/stories/2026-08-21_barley_hapselect_run.md`](stories/2026-08-21_barley_hapselect_run.md);
thread [#028](project_record/threads/028-hapselect-bridge.md). Pipeline and
outputs live outside the package, in `analysis/barley_hapselect/`, and are
delivered to `GRDC_project/hapSelect_/`.

- The bridge carried 1107 samples x 9856 markers x 2 through
  `HapSelect::pairwise_ld()` (7,180,212 within-chromosome pairs, 13.9 min on 15
  workers), `def_blocks()` and back into the container. The catalogue then ran
  **unchanged** on HapSelect-defined blocks — the claim #028 was built to test,
  now demonstrated on real data rather than a simulated panel.
- **Found: the missing code can capture a marker's trait effect** — INDEX
  discrepancy **23**, documented not fixed. `.gp_allele_dosage()`
  (`R/genomic_prediction.R:83`) counts alleles without the missing-token set that
  `compute_marker_stats()` passes (`R/marker_stats.R:231`), so PLINK's `"0"` ranks
  as an allele. On this panel it warned for **8751 of 9856** markers and put the
  estimated effect on `Allele == "0"` for **34** of them. `.meta_to_dosage_matrix()`
  has the same gap (`R/format_helper.R:2143`). NA-ing the missing code at build
  time removes the warning and moves every effect onto a real allele, leaving LD
  and blocks bit-identical. The reach is **panel-wide, not just those 34**: the
  run fits the same rrBLUP model through both packages and compares them, and
  cross-package effect agreement goes from **r = 0.936 / 0.928 before the fix to
  r = 1.000000000000000 after** (max |diff| ~1e-9) — so every one of the 8751
  warned markers carried a slightly wrong effect.
- Two `hapselect_blocks_to_meta(preserve_fields = TRUE)` restorations fired on the
  real container (`phenotype`, `sample_info`), confirming discrepancy 20 bites in
  practice and that the bridge's guard is load-bearing rather than defensive.

### Prototype (2026-08-19, no code change)

**HapSelect bridge — a GenoMeta entry point for the SelectionTools
replacement.** Story:
[`docs/stories/2026-08-19_hapselect_bridge.md`](stories/2026-08-19_hapselect_bridge.md);
thread [#028](project_record/threads/028-hapselect-bridge.md).

- `dev/hapselect_bridge/hapselect_bridge.R` — **new**, converts between the
  GenoMeta container and [HapSelect](https://github.com/wrshf7/HapSelect) 0.1.0
  (LD → haploblocks → localGEBV → GA parent selection). Every function is
  prefixed `hapselect_`; `hapselect_inputs(meta, trait)` is the entry point and
  `hapselect_blocks_to_meta()` the return path. The join is HapSelect's **first**
  stage: it needs a dosage frame + map, haploStack needs back the
  `(Chrom, Marker, Block)` table `attach_block_info()` already takes — the same
  shape `define_haplotype_blocks()` returns — so **nothing downstream of
  `attach_block_info()` changes**, which is what makes it a drop-in for the
  SelectionTools route (roadmap #002, never installable here).
- The conversion that matters: haploStack effects are **per allele**, HapSelect
  effects are **per marker on a dosage scale**, so the bridge emits the **allele
  substitution effect** `eff(counted) − eff(other)`. The two models then differ
  by a per-block constant only — block variances and sample spread, which is
  everything selection reads, are identical. Verified rather than asserted:
  difference-sd ≤ 3.22e-16, per-block correlation 1.000000000000, block-variance
  difference ≤ 4.44e-16 over 49 blocks × 60 samples.
- `dev/hapselect_bridge/run_hapselect_bridge_demo.R` — **new**, 22 checks over
  the whole route on a **simulated** founder-haplotype panel (a local simulator
  was needed: `sim_geno_matrix()` draws alleles independently, so a
  `sim_dataset()` container has no LD and `def_blocks()` returns only
  singletons). **22/22 pass**, including `run_haplotype_catalogue()` running
  unchanged on HapSelect-defined blocks. Exits non-zero on failure.
- **No package code changed.** Nothing in `R/`, `NAMESPACE`, `DESCRIPTION` or
  any existing file; **HapSelect is deliberately not a declared dependency** (not
  Imports, not Suggests) — the bridge is `source()`d and guards itself with
  `hapselect_require()`. Promotion to `R/` would mean adding it to `Suggests`.
- Two findings, both **documented not fixed** — INDEX discrepancy **22** for the
  first, discrepancy 20 for the second: `compute_marker_stats()` reports the
  **same allele as `Major` and `Minor` at exactly MAF 0.5** (`which.max` and
  `which.min` both take the first tied index), which silently zeroes a
  substitution effect — 3 markers in 90 here, dropping cross-package correlation
  to 0.956 before the bridge stopped relying on it; and `attach_block_info()`
  dropping `marker_effects`/`phenotype`, which the bridge restores by default.

### Documentation (2026-07-29, no code change)

Maintenance pass — story:
[`docs/stories/2026-07-29_geno_meta_metadata_doc.md`](stories/2026-07-29_geno_meta_metadata_doc.md);
the four findings below are open as INDEX discrepancy **16**.

- `docs/geno_meta_structure.md` — new **`metadata` keys** section documenting
  every key the container can carry (`ploidy`, `char_dict`, `n_samples`,
  `n_markers`, `source_format`, `encoding`, `family_history`,
  `signal_selection`), each with its type, its writer, what later updates it, and
  a `str()` dump of real content. Verified against the code, not the prose: the
  eight char_dict shapes the builders actually emit (list vs named integer
  vector; `A,C,G,T,-` vs `A,T,C,G` ordering), the `"selectiontools"` /
  `"nucleotide"` / `"numeric"` values that reach `metadata$encoding` versus the
  `HS_ENCODINGS` set and the tokens removed in 0.5.1, the `source_format` labels
  each reader writes, `family_history`'s `FID`/`IID`/`Label` content and its
  `gen_1.F_1` write-back into `sample_info$Family`, and `signal_selection`'s four
  fields. Records four behaviours worth knowing at the boundary: the counters
  are written and never re-derived (so they go stale on a hand edit),
  `signal_selection$n_blocks` is not refreshed by `filter_geno_meta()`, a
  named-vector `char_dict` errors inside `.detect_minor_alleles()`
  (`char_dict[["na_char"]]` → subscript out of bounds) where the list form
  returns `NULL`, and `recode_geno_meta()` trusts the stored `encoding` rather
  than the array — on a `"numeric"` container it silently returns
  `1 → "A"`, `2 → "C"`, `0 → NA`. Documented, not fixed; no code changed.
- Same file: `block_info` row now lists the `(Pos, Num_SNPs, Block_length)`
  coordinate columns added by `expand_block_info()` /
  `build_block_location_info()`.

### Documentation (2026-08-13, no code change)

Second maintenance pass over the container manual, re-verifying
`docs/geno_meta_structure.md` against the **0.5.5** source (it was last checked at
0.5.1, before four releases of new field consumers). Story:
[`docs/stories/2026-08-13_geno_meta_carryover_doc.md`](stories/2026-08-13_geno_meta_carryover_doc.md);
the three findings below are open as INDEX discrepancy **20**.

- `docs/geno_meta_structure.md` — new **"What each transform carries forward"**
  section: a per-transform matrix of which optional fields (`sample_info`,
  `marker_effects`, `phenotype`, `metadata`) survive each rebuild. Every
  transform except `attach_family_history()` rebuilds through `new_geno_meta()`,
  and a field the rebuild call does not name comes back at the constructor
  default — `NULL`, or a dummy FAM — while the result still passes
  `is_geno_meta()`, since none of those fields is required. Three findings,
  confirmed by running each call on a live container (`sim_geno_meta()` +
  effects + phenotype + traced pedigree + a signal selection):
  - **`attach_block_info()` is destructive beyond blocks.**
    `.replace_block_in_meta()` (`R/format_helper.R:2289`) forwards only
    `geno_arrays`, `map_info`, the two block tables and `metadata`, so
    re-blocking also **drops `marker_effects` and `phenotype` and resets the
    whole FAM** (`Family`, `Parent1`, `Parent2`, `Sex` → `"0"`) — while
    `metadata$family_history` rides on unchanged, leaving the pedigree table and
    the tracing history out of sync in the one direction
    `attach_family_history()` exists to close. Only the dropped `block_info`
    signal columns are announced.
  - **`recode_geno_meta()` drops `phenotype`.** It forwards `sample_info` and
    `marker_effects` but not `phenotype` (`R/catalogue_profile.R:222`).
  - **`metadata$signal_selection` outlives the columns it describes.**
    `attach_block_info()` drops the `block_info` `Signal_*` columns (with a
    message) but does not clear the selection record, so a container can carry a
    selection record with nothing selected.
  The section closes with the safe build order — read → recode → blocks →
  effects → phenotype → pedigree → signal markers — and notes that
  `filter_geno_meta()` is the only transform lossless across all four fields.
- Same file: coverage of the 0.5.2–0.5.5 functions that read or write container
  fields — `export_marker_effects()` / `plot_marker_effects()` as read-only
  `marker_effects` consumers, **`predict_gebv()`** as a `phenotype` writer
  (appends `paste0(trait, "_", gebv_suffix)`, default `"<trait>_GEBV"`, via
  `attach_phenotype()` so existing traits survive) alongside
  `fit_genomic_prediction()`, and the 0.5.5 row-order contract
  (`use_sample_order()`, the `haploStack.sample_order` option, `sample_order=`,
  `plot_meta_alignment(check_order=)`) with the note that none of it touches the
  container.
- Same file: new **"Reading `geno_arrays`"** section — tokens are stored as
  character, `geno_meta_views()` / `.meta_to_allele_matrix()` are the two read
  seams, and the container's **NA-any rule** (a cell is uncalled when *any* copy
  is missing) is stated where consumers meet it, with the 0.5.2 `fill_na`
  placeholder as its consequence. Also: the offspring builders
  (`.cross_block()`, `update_geno_meta_from_matrix()`, the driver's advance step)
  named as `geno_arrays` creators that carry the parents' `metadata`, and a
  version stamp at the top of the file recording what it was verified against.
- No code, tests, `NAMESPACE` or `DESCRIPTION` changed; the version stays
  **0.5.5**. The package was loaded with `pkgload::load_all()` only to exercise
  the transforms.

### Deliverable (2026-08-13) — `ready_to_use/genoMetaBuild/`, source-only container layer

A second delivery shape beside the tarball and templates, for users whose job is
only **getting data in**: the container layer of 0.5.5, `source()`d rather than
installed. No package code changed — the fork's `R/` files are verbatim copies.

- **`load_genometa.R`** sources 14 files into the session, checks `data.table` /
  `dplyr`, and sets the two `haploStack.*` options that `.onLoad()` would.
- **In scope:** build (`read_plink_ped()`, `build_geno_meta()`, `new_geno_meta()`,
  `sim_geno_meta()`), validate (`is_geno_meta()`, `geno_meta_views()`), filter,
  I/O (`save`/`load_geno_meta()`, `write_geno_meta()`, `read_effect_file()`,
  `export_marker_effects()`, `read_phenotype()`) and the container updates
  (`attach_block_info()`, `attach_phenotype()`, `attach_family_history()`,
  `attach_signal_markers()`, `recode_geno_meta()`, `build_block_location_info()`).
  **Out:** catalogue, local GEBV, plotting, genomic prediction, simulation.
- **`GenoMeta_build_guide.Rmd`** — a runnable new-user guide restricted to
  **PLINK PED/MAP and nucleotide tables**, on simulated demo data
  (`example_data/`, regenerable via `make_example_data.R`).
- **One non-verbatim file, `R/00_engine_shim.R`**: `.HS_ID_COLS` and
  `.lead_identifier_cols` lifted verbatim from `catalogue_engine.R` so the fork
  need not carry the catalogue layer, plus a **fork-specific `.hs_use_rcpp()`
  that always returns FALSE**. The package's version would return TRUE whenever
  an installed haploStack is found — and then fail, because the kernels are
  namespace-private and a *sourced* `split_haplotypes()` cannot see
  `split_genotype_cpp` (reproduced on R 4.4.3 with 0.5.5 installed). Forcing the
  pure-R fallbacks also makes the fork behave the same on every machine
  regardless of what version happens to be installed; the two paths are
  parity-tested in-package, so results are identical.
- Verified on **R 4.4.3 and R 4.2.0**, identical on both: **109** functions
  defined, `codetools::checkUsageEnv()` clean apart from the three now-unreachable
  Rcpp kernel names, all 13 copied files byte-identical (SHA-256) to their
  parents, and the guide knits end to end. Recorded in thread
  [RTU](project_record/threads/ready-to-use-workflow.md) (caution 6);
  `ready_to_use/README.md` points at it.

### Distribution (2026-08-17) — GitHub release assets for 0.5.5

The first **tagged GitHub release** of the package, cut from the published repo
`CCS-voidBird/haploStack` (main @ `0073898`). No package code changed and the
version stays **0.5.5** — this is a packaging event, recorded in thread
[#026](project_record/threads/026-distribution-build.md).

- **`ready_to_use/release_v0.5.5/`** (new folder) — `haploStack_0.5.5.tar.gz`
  (423 KB, 249 files, portable source), `haploStack_0.5.5.zip` (1.47 MB, Windows
  binary for R 4.4.x x86_64) and `RELEASE_NOTES.md` (the release description,
  with SHA-256 sums for both assets).
- **A Windows binary is now shipped for the first time**, alongside — never
  instead of — the source tarball. It exists so Windows users without Rtools can
  install; the source tarball remains the portable, all-platform artefact and the
  standing convention. The binary was built **from the tarball**
  (`R CMD INSTALL --build`), so the two assets are provably the same code.
- The flat lineage archive `ready_to_use/haploStack_*.tar.gz` is **untouched** —
  the 2026-07-31 `haploStack_0.5.5.tar.gz` built from the dev tree is still there
  as the record's evidence baseline. The new tarball is a **different build of
  the same code**: it carries the published repo's `DESCRIPTION` identity
  (Chensong Chen, `CCS-voidBird/haploStack`), its GitHub `README.md` and its
  `quick_start.Rmd`, which is why it does not overwrite that slot. This
  identity split between the dev tree and the published repo is open as INDEX
  discrepancy **21**.
- Verified as a recipient would: the binary was installed into a throwaway
  library and the package's **own test suite run against the installed build** —
  **1641 pass, 0 fail, 0 error**, 12 skips (optional deps absent: `vcfR`,
  `future.apply`, and two guards that skip *because* `circlize` /
  `genomicSimulation` are present). All six Rcpp kernels register.
- Still open from #026 and unchanged: no Linux/macOS install has been observed,
  and `devtools::check()` remains un-run (0.5.2–0.5.5) for want of pandoc.

### Metadata (2026-08-17) — canonical package identity corrected

The dev tree's `DESCRIPTION` named the author as
`person("University of Queensland", "Bioinformatics", …)` with a personal Gmail
maintainer address, and pointed `URL` / `BugReports` at
`github.com/uq-bioinformatics/haploStack` — an org that was never confirmed to
exist and does **not** host the package. It is now the same identity the
published repo ships:

| Field | Was | Now |
|---|---|---|
| `Authors@R` | University of Queensland Bioinformatics | **Chensong Chen** |
| maintainer email | `ccsconnet@gmail.com` | **`chensong.chen@uq.edu.au`** |
| `URL` | `uq-bioinformatics/haploStack` | **`CCS-voidBird/haploStack`** |
| `BugReports` | `uq-bioinformatics/haploStack/issues` | **`CCS-voidBird/haploStack/issues`** |

- `man/haploStack-package.Rd` regenerated (`devtools::document()`, roxygen2
  7.3.3 — matching `RoxygenNote`). **It was the only `.Rd` of 150 to change** and
  `NAMESPACE` is byte-identical, so the regeneration introduced nothing else.
  `Authors@R` re-parses and derives `Author: Chensong Chen [aut, cre]` /
  `Maintainer: Chensong Chen <chensong.chen@uq.edu.au>`, matching the published
  copy exactly.
- **No code, tests, exports or `Version` changed** — the version stays **0.5.5**,
  and the 2026-08-17 release assets need no rebuild: they were cut from the
  published clone, which already carried this identity. This is the dev tree
  catching up to them.
- Closes the identity half of INDEX discrepancy **21** and settles the
  personal-Gmail item left open under `dev/questions_sheet.md` **Q20.f**. Still
  open there: whether `BugReports` should name `/issues` at all (Q20.f.1 decided
  Issues stay **off**, so that URL may not resolve), the `LICENSE` two-line stub
  form, and the missing `CITATION.cff`. `LICENSE` is untouched — its copyright
  holder is still "University of Queensland", which is a rights statement, not an
  authorship one, and it already matches the published copy byte for byte.

## [0.8.3] - 2026-10-20 — A ring can name categories

Needed for the owner's real-data circos request of 2026-10-20 (VP1
manuscript, MorexV3): "a ring that recorded both all blocks and the selected
51 blocks (use different color to specify)". Before this, `circos_track()`
coerced every value to numeric, so a two-set ring could only be a 0–1
colourbar. Thread [#033](project_record/threads/033-block-annotation-native.md).
150 exports (unchanged).

### Added

- **Categorical rect rings.** In `circos_track(type = "rect", value_col =)`,
  a character or factor column makes the ring categorical:
  - one colour per category, from a palette named by category (it must name
    every category, else the call stops naming the missing ones), or from an
    unnamed palette (default Okabe–Ito) recycled in level order;
  - level order follows the factor's levels, else the named palette's order;
  - rows are drawn in level order, so later levels sit on top (a highlighted
    subset stays visible over thousands of grey blocks);
  - the legend shows one square chip per category instead of a colourbar.
  The track carries `$categories`. Numeric and logical columns stay
  continuous, as before, and non-rect tracks are unaffected.

### Tests

- `test-circos-meta.R` 145 → 160 assertions:
  - level order and colours from a named palette, a factor and an unnamed
    palette;
  - the error for a missing colour;
  - the default palette;
  - logical and point columns stay continuous;
  - a categorical ring renders warning-free with the centre and bottom
    legends;
  - the square chip draws.
- Suite **3 046 pass / 0 fail / 12 skips** against the dev tree (= 0.8.2's
  3 031 + 15).

### Verified

- Used by the VP1 MorexV3 circos run
  (`GRDC_project/.../Archive_VP1_manuscript/haploStack_MorexV3_circos/`): all
  4 223 blocks in grey, the 51 selected blocks in orange, on the installed
  0.8.3 build.

## [0.8.2] - 2026-10-20 — The legend can leave the circle

Owner request 2026-10-20, a stage-6 follow-up on thread
[#033](project_record/threads/033-block-annotation-native.md): the 0.8.1 SIM
circos (`data/SIM/.../haplostack_annotation/circos/sim_meta_circos_183mm.png`)
has its legend in the centre — "add the option that can place it to bottom
(right/left) or somewhere else to make clean main plot". 150 exports
(unchanged).

### Added

- **`plot_meta_circos(legend_position = )`** — `"center"` (the default and
  0.8.1's placement; `"centre"` accepted), `"bottomright"`, `"bottomleft"`,
  `"topright"`, `"topleft"` (one column tucked into that corner, outside the
  circle), `"bottom"` (the entries in rows under the circle, wrapped at the
  page width, each row centred, a chip on its row's colourbar centre line), or
  `"none"`. It is the last argument, so positional calls are unaffected;
  `legend = FALSE` still turns the legend off.
- **An outside legend keeps the circle at its declared size.** The canvas
  grows down (bottom positions) or up (top positions), never sideways, so
  the declared journal width holds. It grows only as far as the legend box
  needs to clear the chr bar's tick labels (5 mm reach + 1.5 mm air); a
  corner box uses the free corner first. On the SIM six-ring figure the corners add
  nothing at 183 mm and make the 89 mm figure 112.3 mm tall; `"bottom"`
  gives 183.8 mm and 99.8 mm (two rows).
- **`hs_meta_circos$height_mm(width_mm)`** — the figure height at a width.
  `save_journal_figure(height = NULL)` takes the page height from it; an
  object saved by ≤ 0.8.1 has no such field and stays square.
- The legend is measured in mm before any device exists, from Helvetica's
  advance widths (R's own `grDevices/afm/Helvetica.afm.gz`; Arial, the cairo /
  ragg sans, is metric-compatible) — internal `.hs_text_mm()`,
  `.hs_legend_layout()`, `.hs_circos_frame()`, `.hs_legend_anchor()`,
  `.hs_draw_legend()`.

### Fixed

- **Colourbar titles overprinted their bars** in the centre legend (the
  0.8.1 SIM 183 mm PNG: "genes / block", "GY effect", "NFNB effect"). Two
  causes. The title sat 0.7 mm under its row top with the bar top only
  1.4 mm down, so descenders reached the bar. And the legend's mm went
  through `.mm_frac()`, whose mm come out ~8 % short on the page (circlize
  plots with R's 4 % axis padding per side), so the geometry shrank while the
  pt text did not. The legend is now laid out in true mm read off the live
  device (`par("usr")` / `par("pin")`): title centre 0.9 mm and bar top
  2.3 mm under the entry top. Ring heights still go through `.mm_frac()` and
  are unchanged — INDEX discrepancy **36**.
- **The centre-legend fit check measures the box.** The box's half-diagonal
  plus 1.5 mm must fit inside the hole, less whatever an innermost label
  ring's text reaches into it. The old rule was hole radius ≥ 18 mm. The SIM
  figure comes out the same: the legend fits at 183 mm and is skipped at
  89 mm. The skip message now names `legend_position`.

### Tests

- `test-circos-meta.R` 66 → 145 assertions:
  - the AFM table matches R's file;
  - the stack and flow layout arithmetic;
  - for every outside position at 89 and 183 mm, the anchored box clears the
    circle by the tick-label reach plus air, stays on the page, and the page
    aspect equals `height_mm`;
  - every position renders without skipping;
  - the centre-skip message names `legend_position`;
  - the saver's page height (pdftools-gated).
- Suite **3 031 pass / 0 fail / 12 skips** against the dev tree (= 0.8.1's 2 952 + 79).

### Verified

- The SIM six-ring figure was re-rendered at every position × 89 / 183 mm and
  checked by eye: the corner legend sits clear of the chr 3 / chr 4 tick labels,
  the bottom strip is one row at 183 mm and two at 89 mm, and the titles clear
  their bars. Renders are in
  `data/SIM/.../haplostack_annotation/circos/legend_position_0.8.2/`.

## [0.8.1] - 2026-10-19 — The catalogue draws the annotation: the circos meta plot and the journal standard

BAN stage 6, the viz phase (thread
[#033](project_record/threads/033-block-annotation-native.md); Q25 n/o
accepted 2026-10-19; figure quality per `dev/figure_standards_nature.md`).
All of it lives in the **catalogue layer** and consumes the standard
annotation / catalogue table formats from any source — never the annotation
module (the 0.8.0 layer boundary).

### Added

- **The circos meta plot** (`R/catalogue_circos_meta.R`; 139 → 150 exports
  across the release) —
  the polar analogue of `plot_meta_alignment()` (Q25.n):
  - `plot_meta_circos(tracks, ...)` — the chr–pos bar outline (one sector per
    chromosome in natural order, Okabe–Ito colours, labelled Mb ticks with
    unlabelled minors, the chromosome name inside the bar) around a
    \strong{ranked track list}: `tracks[[1]]` is the outermost data ring, the
    last the innermost; `rel_heights` scales per ring. A centre legend names
    each ring, with a mini colourbar for rect rings. Returns an
    `hs_meta_circos` whose `$render(width_mm)` re-draws the figure at any
    declared width.
  - `circos_track()` + five constructors: `track_annotation()` (gene / QTL
    counts from the annotation blocks table, container, or bare df),
    `track_block_effect()` (one trait ring from `analyze_block_catalogue()`,
    diverging anchored at 0), `track_block_variance()` (line or rect;
    `Start_Pos`/`End_Pos` spellings accepted), `track_marker_stat()` (points;
    df or GenoMeta + trait = |effect|), `track_marker_labels()` (the Q25.o df
    through the standardizer; leader lines + greedy angular deconfliction,
    the Python circos' rule).
  - **Journal sizing is structural**: ring heights are physical mm and text
    is specified in pt *at the declared width* — the engine computes `cex`
    per device, so the 5 pt floor holds at 89 mm and 183 mm alike. The
    radial budget is checked: too many rings at single column ⇒ a message
    recommending double column, and the centre legend is skipped (with a
    message) when the hole cannot hold it. Axis numerals and chromosome
    names render screen-upright (`facing = "downward"`) because tangential
    text at the 3 / 9 o'clock flip boundary can render a "600" as "009" and
    a chromosome "6" as "9" — found by eyeballing the SIM renders.
- **The journal figure standard as code**:
  - `theme_hs_journal(base_size = 7)` — the `figure_standards_nature.md`
    contract post hoc: 5–7 pt sans, black axis text, 0.5 pt lines, outward
    ticks, no grid, 2 pt margins, borderless bottom legend, titles nulled.
  - `save_journal_figure(plot, file, width = "single"|"onehalf"|"double"|mm)`
    — cairo_pdf (fonts embedded) or ragg PNG at dpi; warns when a ggplot
    still carries a title; **re-renders** an `hs_meta_circos` at the saved
    width; **refuses** a plain `recordedplot` (rescaling breaks the pt
    contract).
- **The Q25.o overlay normaliser** (`R/catalogue_annotation_overlay.R`):
  - `standardize_marker_annotations(x, meta, match_bp = 10 kb)` — the single
    owner of the marker-level annotation / QTL df contract: `Marker` route or
    `Chrom`+`Pos` route (nearest marker within `match_bp`), optional
    `Type`/`Trait`, unplaced rows kept with one warning; returns `Marker,
    Chrom, Pos, Block, Name, Type, Trait, matched_dist_bp`.
  - `attach_marker_annotations(meta, x)` — stores the resolved table as
    `meta$annotation$custom` + its own provenance; independent of the gene
    tables (attachable to a never-annotated container, in either order);
    `filter_geno_meta()` prunes rows whose marker left (position-only rows
    stay), `attach_block_info()` drops the whole annotation as before,
    `is_geno_meta()` validates `custom` markers and now accepts a
    custom-only annotation.
- **Tests** (`test-circos-meta.R`, +66 assertions; suite 2 886 → 2 952 / 0
  fail / 12 skips): track-spec validation, all five constructors from their standard
  sources, render + re-render-by-width to files, the deconfliction rule,
  theme contract, saver page size (pdftools-gated) / title warning /
  recorded-plot refusal / circos re-render, both standardizer routes +
  unplaced warning, custom storage / coexistence / carry-forward /
  validation.

### Verified

- SIM test drive (`data/SIM/.../haplostack_annotation/circos/`): six rings —
  genes/block, GY block variance (the dataset's HapSelect sheet), GY and
  NFNB block effects (`analyze_block_catalogue` on the true effects),
  |true GY effect| points, 12 ground-truth major-QTL labels — on the real
  MorexV3 pseudomolecule extents; saved single + double column PDF + 600 dpi
  PNG through `save_journal_figure()`. Rendered output eyeballed at both
  widths across four iterations (tick collisions, legend collisions, the
  6/9 digit flips) before acceptance; the classic
  `plot_block_effect_circos()` output on the same analysis saved alongside.

### Deferred (still #033)

- `annotations =` overlay arguments on `plot_block_effect_manhattan()` /
  `plot_block_effect_genome()` / `plot_marker_signal()`.
- Rebuilding `plot_block_effect_circos()` on the meta engine: its acceptance
  ("visually unchanged") conflicts with the engine's upgraded frame (ticks,
  upright digits, legend); an owner call on which wins — until then the
  classic implementation stands untouched.
- Stage-4 remainder as before (workshop toggle, bridge §0, io reference).

## [0.8.0] - 2026-10-19 — Which genes sit in each block, without leaving R

BAN stages 1–3 + 5 (thread
[#033](project_record/threads/033-block-annotation-native.md);
`dev/PRD_block_annotation_r_native.md`, Q25 a–p accepted 2026-10-19). The
offline annotation route — a local GFF3, optionally an Ensembl xref dump,
optionally a QTL table, against the container's own block spans — no longer
crosses the Python boundary: it is three `data.table::foverlaps()` passes on
tables the `GenoMeta` already holds. The Ensembl REST route and the Python
figures stay in `block_annotation` (the R plotting phase is 0.8.1, Q25.n/o/p).

### Added

- **Five exports (134 → 139)**, all table producers — no plotting in the
  annotation module (owner clarification 2026-10-19):
  - `read_gff3_genes(path, feature_type, chr_map)` — gene rows of a (gzipped)
    GFF3; `##FASTA` tail ended, `gene_biotype=`/`gene_id=`/`gene_name=` NCBI
    synonyms accepted, IDs kept as the file spells them (`gene:` prefix
    included), no attribute required but the coordinates (`R/annotation_read.R`).
  - `read_gene_xrefs(path)` — an Ensembl Plants xref dump collapsed to one row
    per gene, replicating `load_gene_xrefs()` exactly: GO and
    INTERPRO/INTERPROSCAN as in-order de-duplicated `ID:description` joins,
    `"LOC" + xref_id` symbol from NCBI GENE, description priority Swiss-Prot >
    TrEMBL > NCBI GENE > RFAM > tRNAscan-SE, no-InterPro warning (Q25.g).
  - `compute_block_annotation(meta, genes, xrefs, qtl, chr_map, flank_bp,
    strip_prefix, browser_base)` — object in, object out: per-block gene (and
    QTL) overlap, the (block, gene) record table, and the (marker, gene) pair
    table with the **marker-in-its-own-block gate** (`R/annotation_compute.R`).
    1-based inclusive overlap on both sides, identical to the Python route.
    `chr_map` is the only chromosome mapping — explicit, no species table — and
    a panel-wide zero-overlap prints both vocabularies side by side (Q25.e).
    Gene ids stored as spelled; the `gene:` prefix is stripped **only at the
    xref join** (Q25.f/g). `browser_base` builds the Ensembl browser URLs the
    Python species config built (explicit, same format).
  - `attach_block_annotation(meta, annotation)` — stores the three tables as
    the optional container element `meta$annotation` (plain data.frames) +
    `metadata$annotation` provenance (inputs, counts, `chr_map`, `flank_bp`,
    date, version); refuses an annotation whose blocks or markers are not the
    container's (Q25.a).
  - `export_block_annotation(x, output_dir, prefix)` — the Python step files
    (`step1_block_coordinates`, `step1_block_marker_map`,
    `step3_local_annotation`, `step3_local_annotation_records`,
    `step4_deep_annotation`), **line-identical to block_annotation 0.3.0 at
    `flank_bp = 0`** (Q25.c).
- **`flank_bp` matching threshold, default 10 000 bp** (Q25.d-iv, owner
  insertion 2026-10-19): pads every block for the gene overlap and the marker
  gate; the reported spans are never padded. The default **deliberately
  diverges from the Python route** (no padding there): parity is defined at
  `flank_bp = 0`, which the golden test pins, and exporting a non-zero flank
  messages the difference. The 10 kb constant (`.HS_ANNOT_DEFAULT_BP`) is
  shared with the 0.8.1 overlay's `match_bp`.
- **Carry-forward** (Q25.i): `filter_geno_meta()` prunes the stored annotation
  to the surviving blocks and markers; `attach_block_info()` **drops it with a
  message** (new blocks invalidate every overlap);
  `attach_phenotype()` / `attach_signal_markers()` / `recode_geno_meta()` and
  the simulation's offspring builders (`.cross_block()`,
  `.stk_rename_generation()`) keep it. `is_geno_meta()` validates a stored
  annotation's keys against the container when one is present.
- **Tests** (+77 assertions, suite 2 809 → 2 886 / 0 fail / 12 skips):
  `test-block-annotation.R` (readers incl. gz + `##FASTA` + NCBI synonyms,
  xref collapse rules, flank monotonicity, the marker–block gate on a
  straddling gene, vocabulary warning + `chr_map`, QTL counts, xref join on
  the stripped id, pos-unit refusal, attach/validate, export schema +
  flank message, all carry-forward rows, `.rds` round-trip `identical()`);
  `test-block-annotation-parity.R` — the **golden parity gate**: at
  `flank_bp = 0` the five exported files are `identical()` line-for-line to
  the Python reference run (`ready_to_use/workshop_products_0.7.7/
  07_annotation/`: 1 774 blocks / 228 with ≥ 1 gene / 3 111 records / 86
  pairs / 70 unique genes), browser URLs included. The reference inputs live
  outside the package, so the test skip-gates on their presence (a deviation
  from the PRD's carved-fixture plan, recorded in thread #033).

### Verified

- Full suite 2 886 pass / 0 fail / 12 skips against the dev tree; parity
  counts + byte-level step3/step4 re-verified from the installed
  `ready_to_use/haploStack_0.8.0.tar.gz` in a throwaway library.
- Test drive on the simulated panel `data/SIM/.../SIM_barley_MorexV3_500x10K`
  (400 DH × 10 000 markers, 3 856 HapSelect blocks, real MorexV3 coordinates)
  against the real MorexV3 gene models + the simulation's 352 ground-truth QTL
  points: flank 0 → 1 373 blocks with ≥ 1 gene (35.6 %), 14 678 records,
  423 pairs, 321 QTL-hit blocks; flank 10 kb → 1 915 (49.7 %), 15 763,
  440, 328 — counts monotone, 21 s end-to-end including the 35 826-gene GFF3
  read. Products in `data/SIM/.../haplostack_annotation/`.

### Deferred (still #033)

- Stage 4's workshop `ANNOTATION_ENGINE` toggle and bridge-doc §0; an
  `haploStack_io_reference_0.8.0.md`. The `geno_meta_structure.md` element +
  carry-forward rows shipped in this pass.
- Stage 6 (0.8.1): the catalogue-side viz — `standardize_marker_annotations()`
  / `meta$annotation$custom`, `annotations =` overlays, `plot_meta_circos()`
  + track constructors, `theme_hs_journal()` + `save_journal_figure()` per
  `dev/figure_standards_nature.md`.

## [0.7.11] - 2026-10-12 — A blackout is an emergency; a missing understudy is not

Stage 1 of the GBS multi-subtype core upgrade (thread
[#036](project_record/threads/036-gbs-multi-subtype.md);
`dev/UPGRADE_PLAN_gbs_multi_subtype.md`, PRD decision D5, Q26 a–c accepted
2026-10-09). GBS's metric treated every target state as an independent
coordinate, so a block whose rank-1 subtype was fixed still fired a full
absence penalty for its rank-2 twin — picks were spent re-acquiring a +0.03
understudy while true blackouts waited ("phantom emergencies"). The absence
penalty now carries a within-block **rank multiplier**. Default off:
`penalty_rank_scale = 1` is **bit-identical to 0.7.10** (T1 gate across
method × objective × engine × penalty mode × preset, fixed seeds, plus a
fixed-seed driver-run equality). No export change — **134 exports**.

### Added
- **`penalty_rank_scale = 1`** on `gbs_selection()`, `gbs_pair_selection()` and
  `run_stacking_simulation()` (threaded through the driver's `opts` like
  `penalty_mode`; **not** preset-backed — `"kemper"` and `"haploStack"` keep
  their published meaning). A numeric vector of per-rank multipliers of the
  absent-state penalty `beta_t`:
  - **Rank** (`R/stacking_gbs.R` `.gbs_rank_scale()`): within every block,
    states rank by `TargetFreq` descending — effect-free, so a bare sheet
    works, and it coincides with effect order whenever
    `construct_ideal_genotype()` made the sheet. Ties: `CFreq` descending,
    then state key (`Effect`) ascending (Q26.a — deterministic, never
    run-order-dependent). Ranks beyond the vector recycle the last element
    (Q26.b): `scale[min(r, length(scale))]`.
  - **β only.** Applied in `.gbs_need_score()`, `.gbs_set_objective()`,
    `.gbs_marginal_scores()` and `.pair_scorer_gbs()` (where β enters the
    loss-weighted `L_t·β_t`, the hard penalty and the candidate-objective D
    alike). Deviations are untouched: a held minor state at its target has
    deviation ≈ 0 anyway, and the **unscaled** deviation is what polices
    overshoot — a rank-2 subtype drifting high is eating rank-1's slots and
    is pushed back at full strength (PRD D5).
  - **Dial reading** (documented on `gbs_selection()`): `c(1, 1, 1)` ≡ the
    current full-penalty regime — an implicit diversity-hoarding policy;
    `c(1, 0, 0)` ≈ argmax-flavoured; `c(1, 0.01, 0.001)` keeps a rank-1
    blackout at −10…−100 while a rank-2 absence (−0.1…−1) competes with
    ordinary frequency work and rank-3 sinks to the noise floor. Composed
    with `penalty_mode = "dynamic"`, the D4 pathology is neutralised by
    construction: even full risk×hope escalation ×0.01 is small, and the
    lost-cause damping (hope = 0 → s_t = 1) survives underneath.
  - `penalty_rank_scale[1] != 1` is allowed but **messaged** — the rank-1
    (blackout) penalty itself is being rescaled (Q26.c; `.gbs_check_rank_scale()`).
  - `tune_gbs_penalty()` documents the **tuning trap**: its unweighted D2Ideal
    score charges full price for every absent minor state, so the all-ones
    vector wins that grid by definition — the rank vector is benched on
    fitness (plan §5), never tuned there (Q26.g).
- `.gbs_hapn_matrix()` now carries the ideal table's `CFreq` through to
  `mm$targets` when present — solely the Q26.a tie-break input; pool
  frequencies are still computed from the supplied pool, never read from
  `CFreq`.
- Tests: `tests/testthat/test-gbs-rank-scale.R` (70 assertions, plan gates
  T1–T7) — T1 default parity across the full option grid and a fixed-seed
  driver-run equality (the formula oracles of `test-stacking-gbs-modes.R`
  remain the 0.7.10 behaviour pin); T2 single-target blocks invariant under
  any vector; T3 `c(1,1,1)` ≡ `c(1)` exactly, selections and scores; T4 rank
  assignment with ties and recycling, with and without `CFreq`; T5 dynamic-
  mode composition — a pool-extinct rank-2 state takes `−|gbs_penalty|·1·
  scale[2]`; T6 pairs parity at both `gbs_pair_level`s; T7 the phantom-
  emergency fix asserted on need scores directly (rank-1-held/rank-2-absent
  block at −1 vs blackout at −100 under `c(1, 0.01)`, equal under the old
  regime); plus Q26.c message and argument validation, including the
  driver's up-front check.

### Verified
- Suite **2 809 pass / 0 fail / 12 skips** — the 0.7.10 baseline's 2 739
  plus exactly the 70 new assertions; the 12 skips unchanged; the 29 warnings
  are the long-standing ggplot `inherit.aes` set from the catalogue-viz /
  cluster / plot-order / format files, none from the stacking code.
- `devtools::run_examples()` exit 0.
- `ready_to_use/haploStack_0.7.11.tar.gz` (277 files, 545 KB, sources only,
  SHA-256 `f6b6a8da…`) installed from source into a throwaway library in a
  fresh R process: version 0.7.11, **134 exports**, Rcpp kernels live;
  `gbs_selection()` default vs `penalty_rank_scale = 1` identical from the
  installed build, the `c(1, 0.01)` vector and a `gbs_pair_selection()` run
  with `c(1, 0.5)` both live.
- Against the installed build (the analysis skill's loop): suite
  **2 809 / 0 / 12** again (T1); contract smoke **29 PASS / 3 KNOWN / 0
  FAIL**, the same set as 0.7.10 (T2); both `ready_to_use` pipelines run end
  to end on SIM data (T3). The user's default R library was updated to
  0.7.11 (the 2026-09-30 keep-at-latest standing request).

### Notes
- Next: the BarleyBench rank-vector sweep (`[1,1,1]` control, `[1,0.1,0.01]`,
  `[1,0.01,0.001]`, `[1,0,0]`; G2 subset, 3 000 cycles,
  `ref_arm = OHS_pairs_block10`) decides the extinction-guard question Q26.h
  from data. Stages 2–4 (effect-free target sheet v2 + `state_id` switch,
  slots, dominance/genotype pair level) target 0.9.0–0.9.2; 0.8.0 stays
  reserved for BAN (#033).
- GitHub mirror not synced, nothing pushed (as for 0.7.10 — the outward step
  waits for the owner).

## [0.7.10] - 2026-09-30 — A block map first; LD blocking is a side route, named after its backend

Owner direction of the same day: haploStack's recommended use is **a block map
the user already has**, given to `read_plink_ped(block_file =)`,
`build_geno_meta(blocks =)` or `attach_block_info()`, with no LD parameter to
configure. The LD engines are built-in resolvers for users who really want to
build blocks locally. HapSelect is the default, SelectionTools the alternative
("restored, but not stable"). Each gets a function name that cannot be mistaken
for the main route or for the other engine: `hapS_` / `st_`. Thread
[#035](project_record/threads/035-ld-block-resolvers.md). 131 → **134
exports**.

### Added
- **`hapS_ld_blocks(geno_meta, ld_threshold = 0.5, tolerance = 3L, tol_reset =
  TRUE, method = "flanking", start = "LD", attach = TRUE, output_dir = NULL,
  prefix = "hapS_blocks", parallelize = FALSE, verbose = TRUE)`**
  (`R/blocks_hapselect.R`) — the default local resolver.
  - **Origin.** It is the LD half of the dev-tree HapSelect bridge (thread
    #028), promoted into the package. It keeps one chromosome key, the
    minor-allele dosage, and markers in genomic order (HapSelect reads row order
    as genomic order), then runs `pairwise_ld()` + `def_blocks()` +
    `block_obj_to_df()`.
  - **Same answer as the bridge.** On the workshop subset (GRDC barley 1H + 2H,
    300 lines, 2 677 markers) its block map, block IDs included, is
    **identical** to `ready_to_use/workshop_products_0.7.8/`, which the bridge
    made: 1 425 blocks in 7.8 min.
  - **The rule is passed.** HapSelect's own defaults are threshold 0.7 /
    tolerance 1; the resolver defaults to the repo standard, r² ≥ 0.5 /
    tolerance 3 / `tol_reset`.
  - **Unplaced markers** are left out, and become singletons on attach.
  - **Constant-dosage markers.** The "standard deviation is zero" warnings are
    muffled and the markers counted: 4 518 warnings from 4 markers on that
    subset.
  - **Fields kept.** Effects, phenotype and FAM survive the attach.
  - **Side effects contained.** HapSelect replaces the user's `future` plan
    when `parallelize = TRUE`, and leaves progressr's handlers option set; both
    are restored.
  - HapSelect goes into **Suggests** (`Remotes: wrshf7/HapSelect`); haploStack
    installs and runs without it.
- **`st_ld_blocks()`** — `define_haplotype_blocks_meta()` of 0.7.9, renamed
  and unchanged (defaults 0.5 / tolerance 0), documented as **experimental**.
- **`st_ld_blocks_files()`** — `define_haplotype_blocks()` of ≤ 0.7.9 (the
  SelectionTools file route), renamed, moved to `R/blocks_selectiontools.R`,
  and unchanged.
- `R/blocks_ld_common.R`: what the resolvers share — rule check, placed
  markers, the attach / table tail.
- Tests:
  - `tests/testthat/test-hapselect-blocks.R`, 45 assertions (H1–H8): planted LD
    groups at tolerance 0 and 3; fields and chromosome names on attach;
    tolerance; unplaced markers; rule validation; constant-dosage markers;
    side effects; agreement with SelectionTools on planted structure;
  - `test-selectiontools-blocks.R` moved to the new names, plus B9 (4
    assertions: the old names warn and return identical results);
  - the toy panels moved to `helper-ld-toy.R`.

### Deprecated
- `define_haplotype_blocks()` → `st_ld_blocks_files()`, and
  `define_haplotype_blocks_meta()` → `st_ld_blocks()`. Both still work: each
  raises a `deprecatedWarning` naming the replacement, then forwards its
  arguments. The results are identical (B9).

### Documentation
- `attach_block_info()` says the block map is the recommended route and names
  the resolvers. Its example no longer calls the SelectionTools file route.
- README pipeline stage 2: "bring your own block map (recommended)", with the
  resolvers as the side route (was "Block Definition (requires
  SelectionTools)"). Also a blocks paragraph, the file map, and dependencies.
- DESCRIPTION `Description`: blocks "from a block map the user supplies
  (optionally built locally from LD with HapSelect or SelectionTools)" (was
  "via SelectionTools").
- `quick_start` vignette: the block-map route first, then a "no block map"
  section. `ibd_ibs_signal` vignette updated too.
- `run_haplotype_catalogue()` / `build_block_location_info()` see-also point at
  the resolvers.

### Deliverables (`ready_to_use/` — not package code)
- **Workshop 01c.**
  - `BLOCK_FILE` set → `BLOCK_METHOD` defaults to `"file"`, the recommended
    route.
  - Otherwise `"hapselect"` calls `hapS_ld_blocks()`, which replaces the
    sourced dev bridge; the partition is the same.
  - `"selectiontools"` calls `st_ld_blocks()`.
  - An unknown method now stops.
  - `HAPSELECT_BRIDGE` is gone from `_paths.R`.
  - The workshop README gained a "Blocks" section.
- The template scaffold of the user-level `haplostack` skill does the same.

### Not changed
- `run_haplotype_catalogue(ld_threshold = 0.5, tolerance = 0)` still stamps
  those labels on its output file names whatever made the blocks — including an
  external map. Raised with the owner (INDEX discrepancy 35); not changed here.

### Verified
- **Suite**, against the **installed** Windows binary built from the tarball:
  2 752 expectations — **2 739 pass, 0 fail, 0 error**, 12 environmental
  skips. That is 2 690 from 0.7.9 plus 45 + 4 new.
- **Fresh install.** The source tarball was installed from source in a fresh
  R process: 134 exports, kernels live. Both resolvers recover planted LD
  groups, and `define_haplotype_blocks_meta()` raises `deprecatedWarning`.
- **Contract smoke test:** 29 PASS / 3 KNOWN (new C13c, C14, C14b), in a
  throwaway library and in the user's default library, where 0.7.10 is
  installed.
- **`hapS_ld_blocks()` against the dev bridge**, on the workshop subset: the
  block map is identical to `workshop_products_0.7.8/`, IDs included — 1 425
  blocks, 7.8 min. After attaching, effects, phenotype and FAM are identical to
  the input.
- **Workshop, PLINK route**, end to end with the new module 01c: 11.4 min.
  The counts equal the shipped products: 1 425 blocks, 193 with a gene, 4 155
  genes, 86 marker–gene co-localisations. `block_genes.tsv` is byte-identical
  to the shipped file.
- **Workshop and templates on SIM data:** both complete.
- **The user-level scaffold.** The SIM demo passes. The `"hapselect"` branch on
  barley 1H × 100 lines gave 495 blocks in 2.3 min, with 15 all-heterozygous
  markers flagged and 15 221 warnings muffled.
- **Assets:** tarball SHA-256 `8e002e99…3391d` (550 516 bytes, 276 files);
  Windows zip `09bdfb9d…7c860`.

## [0.7.9] - 2026-09-30 — SelectionTools LD blocks straight from a GenoMeta; a blank line is not a record

Two owner requests of the same day. First: check every bridge between a
`GenoMeta` and SelectionTools, and give LD blocking a GenoMeta-in route that
takes the r² threshold and tolerance and handles SelectionTools' files itself
(thread [#034](project_record/threads/034-selectiontools-bridge.md), INDEX
discrepancy **34**, intake 04). Second: fix the blank-line defect of the shared
text reader (INDEX discrepancy **33**, intake 08). 130 → **131 exports**.

### Added
- **`define_haplotype_blocks_meta(geno_meta, ld_threshold = 0.5, tolerance =
  0, attach = TRUE, output_dir = NULL, prefix = "st_blocks", keep_files =
  FALSE, verbose = TRUE)`** (`R/blocks_selectiontools.R`). Writes the container
  into SelectionTools' own formats — matrix `m` (tab-separated, `m/i` header)
  and map `mcp` (with header) — in a **private temporary directory**, under
  aliases (`m1…`, `i1…`, integer chromosome keys; each marker's alleles recoded
  `1, 2, …`, missing `-`), runs `st.calc.ld()` + `st.def.hblocks(ld.criterion =
  "flanking")` **inside it** on the dedicated data set `"hs_ld_blocks"`, maps
  the blocks back to the container's marker IDs and chromosome **names**,
  deletes the directory (`keep_files = TRUE` keeps it), and by default attaches
  the blocks with `marker_effects`, `phenotype` and `sample_info` carried across
  (`attach_block_info()` alone drops them — discrepancy 20). Unplaced markers
  (no finite position > 0, or chromosome `0`/`UN`/`unknown`) are not blocked
  and become singleton blocks on attach; markers with > 5 allele states are
  left out with a warning. `ld_threshold` must be in (0, 1]: **0 ends the R
  session inside SelectionTools**, so it is refused. `attach = FALSE` returns
  `block_marker_map` + `block_info`; `output_dir` writes the canonical files.
  Whole GRDC barley panel (1 107 × 9 856): ~72 s — HapSelect's pure-R LD on the
  same panel is ~14 min.
  **SelectionTools' directory variables.** SelectionTools resolves file names
  against `st.input.dir` / `st.output.dir`. It holds neither itself: its
  `.onLoad` creates both, as `""`, in the **global environment**, and its
  functions find them there.
  - Set to a directory, they are refused, and the error names the remedy.
  - Gone — after `rm(list = ls())` once SelectionTools has loaded, i.e. a
    script re-run in the same session — every `st.*` reader stops with
    *"object 'st.input.dir' not found"*. Both routes now supply `""` for the
    call and remove it again on exit, even on error.

  The second 0.7.9 build treated "gone" as fine and failed there. Caught before
  release by probing the lookup; the assets were rebuilt.
- `tests/testthat/test-selectiontools-blocks.R`, 46 assertions (B1–B8): planted
  LD groups recovered exactly, chromosome names and carried fields on attach,
  tolerance semantics, unplaced markers, rule validation, no files in the
  working directory / temp handling / `output_dir`, the file route agreeing
  with the GenoMeta route on a named-chromosome map, and SelectionTools'
  directory variables (supplied when gone and removed again; refused when set).
- `tests/testthat/test-ped-text-read.R` T8–T8d (33 assertions): blank lines;
  a line of separators kept as in 0.7.8; a file parsed from its lines equal to
  the same file without its blank line.

### Fixed
- **`.read_table_chr()`: a blank line is not a record** (discrepancy 33, intake
  08). The ragged guard counted blank lines out while `fread` returned some of
  them as rows, so two or more trailing blank lines, or any blank line between
  records, made a well-formed file — PED, MAP, effect file, DArT report — fail
  as *"ragged"*. A blank line is now one that is empty or holds only spaces
  and tabs. `fread(blank.lines.skip = TRUE)` skips an empty line; a line of
  spaces or tabs it can return as a row of empty fields (`"\t\t"` under
  `sep = "\t"`), so a file holding one is parsed from its lines with the blank
  ones taken out (`fread(text = )`), every other file from disk as before. A
  line of separators such as `,,` (a spreadsheet's CSV export often ends in
  them) is **not** blank and reads, as up to 0.7.8, as a row of missing values.
  Data lines are counted by file line number, so the strict-mode error now
  names the short record's real line. T7's guarantees (no row ever lost, a
  truncated PED line refused) unchanged. Parity check over 1 157 layouts
  (blank-line kinds × positions × separators × header × LF/CRLF, separator-only
  lines, a byte order mark, ragged files, the effect / PED / MAP / DArT-style
  readers): **every file 0.7.8 read is read identically**; 948 it refused now
  read, and each equals 0.7.8's reading of the same file with its blank lines
  removed.
  *The first 0.7.9 build* (earlier the same day; installed only in throwaway
  libraries, never published) dropped every all-empty row instead. Its guard
  still counted a `,,` line as a record, so a file that 0.7.8 read was refused
  as ragged (18 of the parity cases), as was a CSV holding a `" \t "` line.
  Caught before release by that parity check; the assets were rebuilt.
- **`define_haplotype_blocks()`** (the file route):
  - SelectionTools keeps only the leading integer of a chromosome (`1H` → 1);
    the function joined its output to the map on `(Marker, Chrom)`, so a named
    map lost every position (`Start`/`End` `NA`, `Chrom` `"1"`). It now joins
    on `Marker` and takes `Chrom` from the map. `Marker`/`Chrom` are read as
    text.
  - `ld_threshold` / `tolerance` are validated (0 used to abort R).
  - SelectionTools runs in a private temporary directory holding its own copy
    of the map: its side files (`ld`, `blocks.hmap`, `blocks.hall`,
    `<map>.srt`) no longer appear in the working directory or next to the input
    map. Data are read straight into data set `"ld.blocks"`; the user's
    `"default"` data set is no longer replaced.
  - `st.calc.ld(auxfiles = FALSE)`: the pairwise LD table was written to disk
    and read back into R only to be discarded.
  - **Outputs unchanged otherwise** — on an integer-chromosome map, 0.7.8 and
    0.7.9 return identical tables at 4 threshold × tolerance settings and write
    16 of 16 byte-identical files.

### Documentation
- `write_geno_meta()`: its `"ST"` output is a CSV carrying SelectionTools'
  allele *encoding*; SelectionTools' reader cannot read it (tab / blank / `;`
  separators only) — use `define_haplotype_blocks_meta()`.
- `define_haplotype_blocks()`: the input formats it actually needs (matrix `m`,
  not comma-separated; map `name chrom pos` with a header) and the chromosome
  caveat. The page also had three errors:
  - it listed two output files where the function writes four (`<tag>` =
    `<geno>_<thr>_t<tol>`: `_block_info.csv`, `_block_marker_map.csv`, `.ld.txt`,
    `.haplotype_map.txt`);
  - it named the column `Chr`, which is `Chrom`;
  - its example passed `geno.csv`.

### Deliverables (`ready_to_use/` — not package code)
- Workshop 01c `BLOCK_METHOD = "selectiontools"` calls
  `define_haplotype_blocks_meta()` (it passed arguments `geno_meta_to_matrix()`
  does not have — intake 04) and checks for ≥ 0.7.9. Run end to end on the
  PLINK route (1H + 2H, 300 lines, 0.5 / 3) through module 07: **294 blocks** in
  4.6 s. HapSelect gives 1 425 on the same subset — the two engines do not
  partition alike (thread #034 has the comparison); compare block sets only
  within one engine.

### Verified
- Suite against the **installed** Windows binary built from the tarball: 2 703
  expectations — **2 690 pass, 0 fail, 0 error**, 12 environmental skips
  (2 611 from 0.7.8 + 46 + 33 new). Source tarball installed from source in a
  fresh R process: 131 exports, kernels live, the new function runs.
- Contract smoke test: 26 PASS / 3 KNOWN (new C12 blank lines, C12b a CSV
  ending in `,,` lines, C13 / C13b the SelectionTools route), in a throwaway
  library and in the user's default library, where 0.7.9 is installed.
- Both `ready_to_use` pipelines on SIM data (workshop, templates) complete on
  0.7.9.
- Assets: tarball SHA-256 `096f13c3…6e741ac` (537 468 bytes, 269 files),
  Windows zip `fc8ba907…a44be6`. Two superseded builds of the same day were
  never published: `106243c9…398ab` / `67d7cb2d…17e` (blank-line rework) and
  `9d247ddc…92363` / `360d9ce2…8f8d9` (SelectionTools directory variables).

### Found, not reproducible
- `raw/Barley_filtered_genos_blocks.txt`, the external barley block set, was
  made by SelectionTools (its `b000000` IDs) but **on a different marker set**
  than `raw/Barley_filtered_genos.*`: 399 of its multi-marker members and 514 of
  its singleton positions are not in that 9 856-marker map. It cannot validate
  a blocking of this panel, and its stated settings (0.5 / tolerance 3) cannot
  be checked from the files here. Thread #034.

## [0.7.8] - 2026-09-30 — An exported effect file reads back with the defaults

`export_marker_effects()` writes one **tab-separated file with a header row**
(`Marker\tAllele\tEffect`) per trait, and its roxygen has promised since 0.5.2
that *"a written file round-trips through `read_effect_file()` unchanged"*. It
never did with the reader's defaults (`sep = " ", header = FALSE`): the whole
line landed in one column and the read stopped with *"transform_effects: need
at least 3 columns"*. Found by the `haplostack-analysis` skill's contract smoke
test (intake 01); closes INDEX discrepancy **31**. Released with the workshop
blocking fix (intake 03, discrepancy **32**, see *Deliverables*).

### Changed
- **`read_effect_file(sep = "auto", header = "auto")`** — new defaults (were
  `" "` and `FALSE`).
  - `sep = "auto"` takes a tab if the first non-blank line holds one, else a
    comma, else a space — decided on that line only (internal
    `.effect_file_sep()`), **not** by `fread`'s own candidate scoring: a
    space-separated file of DArT allele IDs (`100013|F|0-12:A>G A 0.5`) splits
    into three fields on `|` as consistently as on a space. A space-separated
    file therefore reaches `.read_table_chr()` with `sep = " "`, exactly as
    before.
  - `header = "auto"`: the first line is a header when its third field is
    present, is not a missing-value token (`NA`, `NaN`, `.`) and does not parse
    as a number (internal `.is_effect_header()`). The file is parsed once with
    no header and the row dropped after the check, so the decision and the
    table come from the same parse. A first data row whose effect is `NA`
    stays data.
  - Explicit arguments keep their meaning. `header` must be `TRUE`, `FALSE` or
    `"auto"`; a file holding only a header line now errors by name.
- **Effect on existing files — verified against 0.7.7 on 17 cases**, two of
  them real files: every file that read under 0.7.7 reads **`identical()`**
  now (space-separated, no header: LF and CRLF, 4 columns, aligned double
  spaces, all-`"T"` alleles, zero-padded IDs, DArT-style IDs, a first row with
  an `NA` effect; and every explicit-argument call). Newly read, each equal to
  its explicit-argument read: tab- or comma-separated, with or without a header
  — including `analysis/barley_hapselect/outputs/04_marker_effects_haploStack_NFNB.tsv`
  (19 712 rows, the export format). **One intended change:** a space-separated
  file with a header line used to read that line as a data row
  (`Marker Allele NA`, silently); it is now recognised and dropped.
- `export_marker_effects()` roxygen states the round trip holds with the
  reader's defaults since 0.7.8.

### Added
- `tests/testthat/test-effect-file-read.R`, 27 assertions (E1–E7): the export
  round trip with the defaults; parity with the old defaults for the
  space / no-header layout; tab × comma × header variants give one table; an
  `NA` first effect is data; DArT-style IDs; explicit arguments; the header-only
  and bad-`header` errors.

### Build
- `.Rbuildignore` excludes `Rplots.pdf` anywhere (`(^|/)Rplots\.pdf$`): up to
  0.7.7 the pattern was anchored to the root, so every tarball shipped
  `tests/testthat/Rplots.pdf`. Tarball: 266 files, 522 KB.

### Deliverables (`ready_to_use/` — not package code)
- **Workshop 01c blocks at the rule it declares** (intake 03, INDEX
  discrepancy 32). The HapSelect route called `def_blocks()` without
  `tolerance`, so it blocked at HapSelect's default tolerance **1** whatever
  `_paths.R`'s `LD_TOLERANCE` (3) said — on 1H × 300 lines, 610 blocks instead
  of 511 from the same LD table. It now passes `tolerance = LD_TOLERANCE,
  tol_reset = TRUE` (the setting of `analysis/dataset_ld_segmentation/`) and
  prints both. Module 03 stamps that rule into the catalogue file names on the
  PLINK LD route (`…_0.5_t3_…`; other routes keep `…_0.5_t0_…`).
  `workshop_products_0.7.7/` is kept as built (tolerance 1, a dated note says
  so); **`workshop_products_0.7.8/`** is the fixed route's run.
- `release_v0.7.8/` staged: source tarball, Windows binary built from it,
  `SHA256SUMS.txt`, `RELEASE_NOTES.md`. Not pushed to GitHub.

### Tooling (agent skills — not package code)
- `.agents/skills/haplostack-analysis/` rebuilt against 0.7.7 → 0.7.8 and
  linked into the repo's `.claude/skills/`; new testing and issue-raising
  sub-skills; a user-level `haplostack` skill for other projects. Findings go
  to the new intake `.scratch/intake/` (`needs-triage`); items 01 and 03 are
  closed by this release.

### Verified
- Suite against the **installed** Windows binary (built from the tarball):
  2 624 expectations — **2 611 pass, 0 fail, 0 error**, 12 environmental skips
  (2 584 from 0.7.7 + 27 new). Source tarball installed into a throwaway
  library in a fresh R process: kernels live, the round trip holds.
- Contract smoke test (`haplostack-analysis/templates/smoke_test.R`): 22 PASS,
  3 KNOWN — check C5 (the round trip) moved from KNOWN to PASS.
- No export change (130). No other function's defaults changed.

### Found, not fixed
- `.read_table_chr()`'s ragged-file guard counts blank lines out but `fread`
  returns some of them as rows, so a file with **two or more trailing blank
  lines, or any blank line between data lines**, is refused as *"ragged"* —
  every text reader, `read_plink_ped()` included (one trailing blank line
  reads, as T7 pins). Same on 0.7.7. Intake 08, INDEX discrepancy 33.

## [0.7.7] - 2026-09-20 — A genotype file is read as text

Three readers parsed their file with `utils::read.table()`, which runs
`type.convert()` on every column. `type.convert()` maps a column whose values
are all `"T"` (or all `"F"`) to **logical** — and in a nucleotide genotype file
that is an ordinary allele column. Nothing in the package could tell the
difference afterwards, because by then the cells read `"TRUE"`.

The rule now is explicit and shared: **a genotype file is read as text.** Every
column arrives as `character` through the new internal `.read_table_chr()`, and
each reader coerces, by name and on purpose, only the columns it knows to be
numeric. Closes INDEX discrepancy **30**.

### ⚠ Fixed — one silent wrong answer, one hard failure

1. **`read_effect_file()` — a wrong answer, the serious one.** An effect table
   whose `Allele` column is uniformly `"T"` (a single-allele export, a
   substitution-effect table, a filtered subset) came back with the allele
   `TRUE`. `TRUE` matches no genotype, and `.compute_effect_size()` scores an
   absent allele as 0 (`R/catalogue_gebv.R`: `long_dt[is.na(Effect), Effect :=
   0]`) — so **every marker contributed nothing and the GEBVs came out flat**.
   Demonstrated on a 3-sample × 3-marker container: local GEBVs
   `0, 0, 0` from the type-converted table against `0.6, 1.7, 0` from the same
   file read as text.
   Nothing warns at **read** time — `read_effect_file()` and
   `transform_effects()` both accept the table without comment. The one guard
   is downstream: `compute_profile_localgebv()` warns through
   `.validate_effect_encoding()` (*"effect_df contains allele(s) not in the
   'nucleotide' encoding dict: TRUE. These rows will not contribute to local
   GEBV"*) and then continues with zeros. So the failure is loud enough to
   catch if you read your warnings and fatal if you do not, and it is invisible
   at the point where the mistake is actually made.
   Effect files whose allele column carries more than one distinct letter — the
   normal genome-wide case — were never affected.
2. **`read_plink_ped()` — a hard failure.** After
   `mode(allele_mat) <- "character"` a logical column reads `"TRUE"`, and
   `.reshape_wide_alleles()` refuses it: *"cells contain multi-character values.
   Wide format requires exactly one character per cell."* The file could not be
   read at all. It is the **allele column** that must be uniform, not the
   marker: for an A/T SNP whose first copy is fixed T across the panel, column
   `2k−1` is all `"T"` while the marker itself still segregates. On
   `data/GRDC/hordeum_vulgare/GRDC_barley_AGG1/plink/geno.ped` (1107 × 9856)
   **511 of 19 712 allele columns** are logical on a full-file read, from 511
   distinct markers, and no marker is monomorphic-T. The count is a property of
   the panel, not the file: one missing call anywhere in a column is enough to
   keep it character, so sub-setting lines moves the number in either
   direction. That is why the defect stayed invisible — the fixture panels in
   the test suite happen never to produce a uniform column.
3. **`read_dart_genotype()` — a mis-parse.** The metadata/sample column split
   is decided by `all(as.numeric(as.character(col)) %in% c(0, 1, 2))`.
   `as.numeric("TRUE")` is `NA`, so a logical column fails the test and
   `first_sample` moves, silently reclassifying a sample column as marker
   metadata. Plausible on any DArT report carrying a T/F flag column.

### Added
- `.read_table_chr(file, header, sep, strict)` (internal,
  `R/format_helper.R`) — the shared text reader.
  `data.table::fread(colClasses = "character")`, an explicit not-found /
  empty-file error, and `""` normalised to `NA`.
  **Ragged input is handled here, not delegated to `fread`.** During this work
  `fread(fill = FALSE)` was found *not* to error reliably on a ragged file:
  given a PED whose first line is wider than the rest it settles on the
  narrower width and **silently drops the wider row**, returning a table with
  one fewer sample and no warning — strictly worse than the defect being
  fixed. The fread call therefore always passes `fill = TRUE` (short rows
  padded, no row ever discarded) and the reader checks the shape itself: every
  data line in the file must appear as a row, and under `strict = TRUE` (the
  PED and MAP) the table must be rectangular. A PED line is fixed-width by
  construction, so a short one is a truncated file and is refused rather than
  padded into missing calls. `strict = FALSE` keeps the lenient behaviour the
  effect-file and DArT readers always had.
- `tests/testthat/test-ped-text-read.R`, 47 assertions: the monomorphic-T
  column that used to abort the read (T1), the all-F and mixed T/F columns
  (T2), identifiers preserved verbatim (T3), a clean panel unchanged and a
  1/2/3/4 PED still declared `selectiontools` (T4), the primitive itself (T5),
  and **`read_effect_file()`'s first tests ever** (T6) — the all-T allele
  column, the zeroed-GEBV consequence, and a mixed-allele file reading as it
  always did. `read_effect_file()` had no test coverage at all before this
  release, which is how the silent case stayed silent. T7 pins the ragged-input
  guards: a wider-first-line file never loses a row, a truncated PED line is
  refused instead of padded, and a trailing blank line still reads.

### Changed
- **Identifiers are no longer type-converted.** A sample or marker id is
  whatever the file spells: `007` stays `"007"` (it used to arrive as integer 7
  and become `"7"`), and an id that reads as a boolean stays `"TRUE"`. This is
  a behaviour change on files with zero-padded or boolean-looking ids, and it
  is the intended one — an identifier that silently changes cannot be joined
  back to its source.
- `read_plink_ped()` is also **substantially faster**: the 43 MB / 19 718-column
  GRDC PED reads in ~52 s where `read.table()` needed minutes.
- No export change (130). No argument added to any exported function.

### Verified
- Full suite under `devtools::load_all()`: **2 584 pass, 0 failed, 0 errors**
  (12 pre-existing environmental skips) — 2 537 from 0.7.6 plus the 47 new.
- On the real GRDC barley PED: `read_plink_ped()` now returns a valid container
  — 1107 × 9856 × 2, ploidy 2, `encoding = "nucleotide"`,
  `source_format = "plink_ped"`, `char_dict` = {T, A}, 396 750 `NA` cells, no
  `"TRUE"` anywhere — and it is **`identical()` to the container the workshop's
  `read_plink_ped_safe()` workaround produced**, for `geno_arrays` and
  `map_info` both. That workaround is withdrawn from
  `ready_to_use/workshop/_common.R`.
- The read path from file to `metadata` is documented end to end in
  `docs/haploStack_io_reference_0.7.7.md`.

### Deliberately unchanged
- `.read_block_file()` still uses `read.table()`. A block file holds
  identifiers and coordinates, never allele codes, and `Start` / `End` are
  wanted as integers — type inference is correct there and the columns cannot
  be uniformly T/F.

## [0.7.6] - 2026-09-19 — Views keep the observed allele of a half-called cell

### Fixed
- `geno_meta_views()` builds the per-copy allele table straight from the
  ploidy layers of `geno_arrays` (new internal `.build_alleles_layers()`), no
  longer by re-splitting `.meta_to_allele_matrix()`'s collapsed genotype
  string. That string is `NA` whenever *any* copy is missing, so a half-called
  cell `(A, NA)` came back as `(NA, NA)` and its observed allele was dropped
  from every consumer of the views: `compute_profile_localgebv()` at both
  levels — and with it every simulation offspring `mean_F1mean`, `max_F1mean`,
  `D2Ideal`, `ideal_hap_gebv` and the haplotype table OHS / GBS / TS select on
  from cycle 1 — `compute_block_catalogue()`, block signatures, clustering,
  block IBD / signal scoring and the marker-stats views. Found 2026-09-18 by
  the barley bench review (manuscript folder `06` §3): founders read from PED
  are never half-called (missing is `0 0`, both copies), but block- and
  genomicSimulation-backend offspring carry 1.6–3.6 % such cells, so the
  founders' stacking ceiling rose from 4.8747 to 4.8842 at cycle 1 and
  cycle-10 population means were 0.3–1.5 % off (individuals by up to 0.06;
  rank agreement 0.98–0.99). The 0.6.0 contract already said "a half-called
  cell keeps what was observed" — the reader honoured it, the views did not.
  `.build_alleles_st()`, `.build_alleles_nuc()` and `.bind_haplotypes()` are
  removed. `.meta_to_allele_matrix()` is unchanged and still serves encoding
  detection and `geno_meta_to_matrix()` export, where a half-called cell is
  written as missing (PED cannot carry one); `.meta_to_dosage_matrix()` still
  counts it as `NA`. Dosage containers are untouched. **Numbers change only on
  containers with half-called cells**; a container without any is
  byte-identical to 0.7.5 (tested). INDEX discrepancy **29**; thread #029 note.

### Verified
- `tests/testthat/test-half-called-views.R`, 412 assertions: the views keep
  the observed copy (nucleotide and SelectionTools encodings); byte-identical
  to the 0.7.5 route when no cell is half-called; local GEBVs equal a hand
  computation from the array at both levels on the Rcpp kernel and the pure-R
  path; block signatures keep the observed copy and fill only the absent one;
  a 3-cycle block-backend run's `ideal_hap_gebv` never rises above the
  founders'.
- `dev/half_called_check/verify_half_called_views.R` on the barley bench
  subset (120 founders × 300 NFNB blocks): one TS cycle yields 2 351
  half-called offspring cells (2.5 %); engine = hand to 1e-17 at both levels
  on both engines; ceiling 4.8747 → 4.7745 → 4.4109 → 4.3618 over three cycles
  (never above cycle 0); the review's single-cell puncture reproduces the
  expected values. Results in `results.tsv` beside it. Every bench table
  produced under ≤ 0.7.5 (`dev/bench_barley_2026-09-17/`, incl.
  `hpc_grid01/`) is pre-fix: rankings stand, decimals are superseded by the
  re-run.
- Suite: 2 537 pass / 0 fail / 12 skips (the same optional-package skips as
  0.7.5: vcfR, future.apply, and the genomicSimulation-absent guard).

## [0.7.5] - 2026-09-18 — `max_F1mean` (top elite) in the simulation summary; genomicSimulation read-back keeps missing alleles NA; benchmark on both backends

### Added
- `run_stacking_simulation()$summary` gains **`max_F1mean`** — the breeding
  value of the cycle's best single offspring (the top elite), beside the
  population mean `mean_F1mean`. It is the realised counterpart of a set-level
  "ultimate genotype" objective (HapSelect's OHS GA), so the two OHS families
  can be compared on the field each one optimises. The founder row carries the
  best founder. Placed after `var_F1mean`; every other column is unchanged.
- `dev/bench_barley_2026-09-17/`: `bench_arm.R` takes `key=value` options
  (`reps=`, `seed=`, `out=`, `metrics=all|nocoan`, `tag=`, `backend=block|gs`)
  and four Kemper-preset GBS arms (`GBSk_{LOG,ARCSIN}_{parents,pairs}`,
  `gbs_preset = "kemper"`); `bench_aggregate.R` reports the top elite;
  `bench_cap_aggregate.R` (single-arm cap runs) and `bench_backend_compare.R`
  (block vs genomicSimulation side by side) are new. **Grid run:** 15 arms ×
  2 backends × 100 cycles × 3 reps (`cycles_100_block/`, `cycles_100_gs/`,
  `backend_compare/`). The genomicSimulation backend (real linkage, 1 cM/Mb)
  reproduces the block ranking — GBS direct pairs at the founder ceiling on
  both (4.81 / 4.84), every GBS / TS endpoint within ±0.07 — and costs the
  OHS family: haploStack OHS ends 0.02–0.14 lower and needs 12–14 more cycles
  to 90 % of its gain, HapSelect's OHS GA 0.23 lower and still segregating.
  Linkage lifts the GEBV-driven arms early (+0.2–0.4 at cycle 10) and the gap
  closes by cycle 20. Kemper's fixed penalty and the haploStack dynamic
  penalty are indistinguishable (≤ 0.07, no consistent sign; LOG pairs fix on
  the same stack). The top elite runs 0.05–0.5 ahead of the mean and reaches
  its plateau one to two cycles earlier; HapSelect's realised elite is behind
  every haploStack arm at every horizon ≤ 100 cycles. Results: thread #013.

### Fixed
- genomicSimulation backend: a missing founder allele travels to
  genomicSimulation as the nucleotide missing token `"N"` and came back as a
  literal allele in every offspring container, so `is_geno_meta()` warned on
  each cycle of a gS-backend run ("1402 allele cell(s) across 83 marker(s)
  hold a missing token" on the barley subset). `.gs_genes_to_array()` now
  restores the nucleotide missing tokens to `NA` (the 0.6.0 container
  contract). Breeding values were not affected on the barley subset (the
  affected cells carry no effect); the warning is gone. Test added.

### Verified
- Suite: 0 failures, 12 skips (optional packages; genomicSimulation installed
  so its bridge and session tests ran).

## [0.7.4] - 2026-09-17 — OHS: `GA_coverage` removed; barley benchmark of the selection arms

### Removed
- **`ohs_method = "GA_coverage"`** (owner's decision, 2026-09-17: "dispose the
  GA_coverage mode as its poor performance"). It scored a parent set by the
  fraction of blocks with at least one positive-effect carrier (genotype level,
  binary) — not an OHS objective — and was the worst set on every metric in
  `dev/review_ohs_gbs_2026-09-17.md` §1.3 (mid-parent 5.4 vs 20.8). `ohs_method`
  is now `c("GA_pairwise", "Greedy")` in `ohs_selection()`, `ohs_pair_selection()`
  and `run_stacking_simulation()`; the deprecated `select_method = "GA"` maps to
  `OHS(GA_pairwise)`. `genetic_algorithm_selection()` itself stays exported (it
  is the stage-4 legacy engine on genotype-level tables) and `.collapse_hap_to_geno()`
  stays (the driver's truncation scorer uses it). Tests, the selection grid,
  the example vignette and the skill sheet updated; the vignette's replicate
  note now names OCS only (the OHS GA restores the RNG since 0.7.3).

### Added
- `dev/bench_barley_2026-09-17/` — benchmark of eleven selection arms on the
  GRDC barley AGG1 panel (HapSelect blocks, haploStack rrBLUP per-allele
  effects, trait NFNB), reduced to 120 random AGG parents × the 300
  highest-variance blocks (838 markers): `bench_prep.R` (the subset),
  `bench_arm.R <arm>` (one arm through `run_stacking_simulation()`, block
  backend, 8 parents, ~112 offspring per cycle, 10 cycles, 3 replicates),
  `bench_aggregate.R` (tables + `trajectories.png`). Arms: OHS `GA_pairwise` at
  `elite_ratio` 1 / 0.1 / one family; HapSelect's `haplotype_parent_selection(
  strategy = "OHS")` plugged in as the parent selector; GBS ARCSIN / LOG,
  parents-first and direct pairing; TS with all / the top 10 % / the top 1
  families by mid-parent kept (offspring per family scaled to the same
  population size). Results in the thread #013 note and the endpoint table.
  Same day, the same arms to **50 cycles** with the per-cycle step timing
  recorded (`bench_arm.R <arm> [cycles]` → `cycles_50/`, `timing_<arm>.tsv`,
  `timing_table.tsv`, `milestones.tsv`): every TS arm fixes by cycle 12–30
  and is overtaken by the OHS / GBS arms at cycles 14–27; GBS direct pairing
  reaches the founders' stacking ceiling, OHS stalls 0.5 below it, HapSelect's
  set-level OHS GA is last but the only arm still segregating. Method cost per
  cycle: TS 0.04 s, GBS 0.06 s (parents) / 0.5 s (pairs), OHS 2.2–3.0 s,
  HapSelect GA 21.7 s. Finding: the `metrics = "all"` block-IBS coancestry
  (`compute_parent_coancestry()`, an R double loop over offspring pairs ×
  blocks) costs ~47 s of every ~50 s cycle-step — recorded under Q18, not
  changed. Sub-test `hapselect_cap/` (`bench_arm.R OHS_hapselect 200 reps=1
  seed=<s> metrics=nocoan`, `bench_cap_aggregate.R`): HapSelect's set-level
  OHS GA alone to 200 cycles reaches 97 % of the founder ceiling (gain 4.69 of
  4.83), 90 % of its final gain at cycle 120 and 99 % at 181, passing the
  fixed haploStack OHS endpoints at cycles 126–131 — its cap is the ceiling,
  reached ~4× slower.

### Verified
- Full suite under `devtools::load_all()`: **2114 tests, 0 failed, 0 errors**
  (12 pre-existing skips; 4 expectations fewer than 0.7.3 = the removed
  `GA_coverage` tests).
- `ready_to_use/haploStack_0.7.4.tar.gz` (264 files, 521 KB, no `.o`/`.dll`)
  installed from source into a throwaway library: 130 exports, kernels live,
  `Greedy`, `GA_pairwise(elite_level = "block")`, `ohs_pair_selection()` and
  the driver with `ohs_elite_ratio` run.

## [0.7.3] - 2026-09-17 — OHS: elite subsetting of the families, the greedy's first pick, the GA's RNG

The OHS objective gains an optional **elite ratio**: instead of every family a
parent set implies, the search may value only the best `q = ceiling(ratio · P)`
of its `P = k (k − 1) / 2` families (owner's request, 2026-09-17, after the
review `dev/review_ohs_gbs_2026-09-17.md`). At ratio 1 nothing changes; the
crossing design has the same size at every ratio — the ratio changes what the
search values, not which crosses are made. Thread
[#013](project_record/threads/013-stacking-simulation.md) (OHS issue 01) and
[#031](project_record/threads/031-direct-pair-selection.md) (the pairs side).

### ⚠ Changes published numbers — OHS runs only

1. **`ohs_method = "Greedy"` starts from the best row sum of `M`.** The
   0.6.0–0.7.2 greedy started from the *first sample column* whatever the data
   (the diagonal of `M` is 0, so the first marginal-gain vector was all zeros
   and `which.max()` returned index 1 — the origin `greedy_selection()` had
   the same defect). Reversing the column order changed its answer; on the
   review's 80 × 60 panel the fix lifts the greedy from 4607.9 to 4693.0 =
   the GA optimum. Pinned by test (first pick, column-order invariance, never
   worse than the old start).
2. **`GA_pairwise` and `GA_coverage` (`genetic_algorithm_selection()`) restore
   the global RNG stream on exit** (the `.gbs_ga_set()` pattern). The GAs' own
   draws are unchanged, so their picks are byte-identical to 0.7.2 at ratio 1
   (parity script, 3 random fixtures), but the driver's crossing step no longer
   runs on a stream reseeded by the selection call: OHS replicates in one
   multi-rep call now differ, as they should. **Closes discrepancy 28 for OHS**
   (`ocs_selection()` still reseeds; Q24 o).

### Added
- `ohs_selection(elite_ratio = 1, elite_level = c("pair", "block"))`:
  - `"pair"` (default): the objective is the sum of the `q` largest family
    ceilings `M[i, j]` within the set. As the ratio shrinks it tends to the
    **single best cross**; a tie-break of `1e-6 ×` the mean family ceiling keeps
    the non-elite parents ordered (below ratio `(k − 1) / P` some parents can
    sit in no elite family and are fillers as far as the elite objective goes).
  - `"block"`: per block, the `q` largest per-block family ceilings
    `M_b[i, j]`, summed over blocks — the elite families may differ from block
    to block, so the small-ratio limit is the **best two haplotypes of every
    block from any two parents**: Kemper et al. (2012)'s "ultimate genotype",
    HapSelect's `Haploid_OHS` (pinned by test against a transcription of that
    fitness; HapSelect's `OHS` differs only by forcing the two haplotypes onto
    different parents). At ratio 1 both levels are the full objective
    `sum(M[S, S])`.
  - Both engines take the objective: `GA_pairwise` through its fitness,
    `Greedy` by re-evaluating the objective for every candidate (the closed-form
    marginal gain is kept at ratio 1). `GA_coverage` has no family reading and
    refuses a ratio below 1.
  - Internals: `.ohs_objective()` (the set objective for a ratio / level),
    `.ohs_pair_block_values()` (per-block family ceilings from the sorted
    haplotype array, vectorised for ploidy 2), `.ohs_col_topq_sum()`,
    `.ohs_elite_q()`, `.ohs_check_elite()`;
    `.pairwise_block_value_matrix(return_array = TRUE)` returns the array too.
- `ohs_pair_selection(elite_ratio, elite_level, shortlist_families)` (direct
  pairing): under the `"parents"` budget the knobs go to the set engine and its
  traversal is returned as before. Under the family greedy, `"pair"` leaves the
  ranking unchanged (the families taken best-first *are* the pool's elite
  subset, whatever the ratio); `"block"` with a ratio below 1 scores the
  **family set** by `sum_b` of the `q = ceiling(ratio · n_families)` largest
  per-block ceilings over its families, re-scored after every pick — the
  marginal gain of a family is `sum_b max(0, v_fb − τ_b)`, `τ_b` the `q`-th
  largest chosen value at block `b` (`.pair_scorer_ohs_block()`, non-additive)
  — over the `shortlist_families` best families of the pool by `M` (default the
  larger of `4 · n_families` and 100). At ratio 1 it ranks exactly as the pair
  level (pinned). The budget, and so the number of families made, is the same at
  every ratio.
- `run_stacking_simulation(ohs_elite_ratio = 1, ohs_elite_level = "pair",
  ohs_shortlist_families = NULL)`, passed to both units; `GA_coverage` with a
  ratio below 1 is refused up front.
- Tests: 8 in `test-stacking-ohs.R` (ratio 1 = the pre-0.7.3 picks for both
  levels; the pair and block objectives against independent transcriptions;
  the block level at ratio → 0 = HapSelect's `Haploid_OHS`; a small pair ratio
  keeps the best single cross; the greedy's first pick and column-order
  invariance; the RNG; validation) and 4 in `test-stacking-pairs-schemes.R`
  (pass-through under `"parents"`; block level at ratio 1 = pair level and the
  pair level unchanged by the ratio; the hand greedy of the first two families
  at a small block ratio; shortlist and budgets; the driver under both units).
- `dev/review_ohs_gbs_2026-09-17.md` (+ `dev/review_ohs_gbs_2026-09-17/`
  scripts and outputs): the review that motivated this release — haploStack's
  `GA_pairwise` (sum of every family's ceiling) vs HapSelect's built-in OHS GA
  (Kemper's set-level ultimate genotype, Yadav et al. 2025 TAG `Σ a_ij = 2`),
  cross-scored on one panel; GBS presets × modes with a 10-cycle grid.

### Changed
- `dev/UPGRADE_PLAN_opt_kemper.md` §2: the claim that `GA_pairwise` is "the
  group's port of that GA (Tong et al. 2025)" is corrected — it is the port of
  the origin's pairwise-sum GA; Yadav / Tong 2025 and HapSelect implement
  Kemper's set-level objective, which the block-level elite limit now reaches.
- No export change (130).

### Verified
- Full suite under `devtools::load_all()` on the final code: **2118 tests, 0
  failed, 0 errors** (12 pre-existing skips for absent optional dependencies).
- Parity with the installed 0.7.2 (separate process, 3 random fixtures):
  `GA_pairwise` picks identical; `Greedy` objective 4741.8 → 4838.5,
  9201.2 → 9225.3, 1912.8 → 2096.0.
- `ready_to_use/haploStack_0.7.3.tar.gz` (264 files, 521 KB, no `.o`/`.dll`)
  installed from source into a throwaway library: 130 exports, kernels live,
  `Greedy`, `GA_pairwise(elite_level = "block")`, `ohs_pair_selection(elite_level
  = "block")` and the driver with `ohs_elite_ratio` run; two replicates of an
  OHS-GA run now differ at cycle 2.

## [0.7.2] - 2026-09-15 — Direct pair selection: the selector

The 0.7.1 placeholder becomes a working unit. Under `selection_unit = "pairs"`
every scheme scores **families** — the expected F1 of one cross — and one shared
greedy adds families until the budget of `pair_selection_budget()` stops it
(plan `dev/UPGRADE_PLAN_pair_selection.md`, Phases A–B; thread #031). The idea:
**a scheme is a scorer, the unit is the engine** — the scoring currency of each
scheme is unchanged, only the object it is applied to changes.

### Added
- `R/stacking_pairs.R`: the shared engine (`.pair_candidates()`,
  `.pairs_greedy()` with the two budgets, the overshoot rule under `"parents"`
  — `"skip"` continues among the parents in use up to the traversal cap,
  `"stop"` ends the search — a per-parent usage cap `max_use`, and the selfing
  rule) and one scorer per scheme:
  - **TS** — the mid-parent value `(g_i + g_j)/2` (the exact expected F1 mean),
    or `ts_pair_score = "usefulness"` = mid-parent `+ usefulness_k · sd(F1)`
    with the F1 standard deviation exact under the block model (`.pair_f1_sd()`).
  - **OHS** — `ohs_pair_selection()` (exported): the family's stacking
    ceiling `M[i, j]`. Under `"families"` the top-`F` of the triangle (exact);
    under `"parents"` the optimum of an additive family value with a name limit
    is the best parent *set* and every family among it — `sum(M[S, S])`, the
    parents-first objective — so that budget runs `ohs_selection()` and returns
    its traversal ranked by `M`. `"GA_coverage"` has no family reading and is
    refused. A self's ceiling is the parent's haplotypes doubled, top-`ploidy`
    per block (`.pair_self_ceiling()`).
  - **GBS** — `gbs_pair_selection()` (exported): the ideal-genotype objective
    on the family's expected offspring `q_t(i, j) = (HapN[i,t] + HapN[j,t]) /
    (2m)`. Two levels (`gbs_pair_level`): `"population"` (default) scores the
    marginal gain of the set objective `F(P)` at the offspring-pool frequency
    `p_t(P) = mean_f q_t(f)`, which **reduces exactly to the parents-first
    `p_t(S)` when `P` is the traversal of `S`** (pinned by test); `"family"`
    scores the family's own distance to the ideal. `gbs_objective = "candidate"`
    is accepted and is separable over a family (`= s_i + s_j`; pinned by test),
    so the pairs default is `"set"` for both presets. `pair_freq = "projected"`
    (Kemper's projection to the final family count; default under
    `"families"`) or `"running"` (default under `"parents"`).
    `pair_penalty = "loss"` (default) replaces the hard absence test by the
    exact loss probability `L_t = prod_f a_f(t)^M` that none of the
    `num_offspring = M` progeny carries state `t` — the parents-first risk term
    made exact for a chosen design (`"hard"` keeps the 0.7.0 rule; the two
    coincide as `M → ∞`, pinned by test). The penalty scale is sized with
    `top_founders`, so `β_t` equals the matched parents-first run's.
    `pair_pool = "shortlist"` (default; the `shortlist_k = 4 · top_founders`
    parents whose single gamete helps most at each step) or `"all"`.
  - Output class `"pair_selection"` (`Parent1, Parent2, score, rank`;
    attributes `n_families`, `n_parents_used`, `pair_budget`, `scheme`) with a
    `print()` method.
- `run_stacking_simulation()`: `.stk_select_pairs()` dispatches TS / OHS / GBS
  to the scorers and validates the family table (names in the pool, no
  duplicate in either orientation, selfs only when allowed; a scheme returning
  no family **errors with its name** — no truncation fallback for this unit).
  New arguments `gbs_pair_level`, `pair_freq`, `pair_penalty`, `pair_pool`,
  `shortlist_k`, `max_use`, `overshoot`, `ts_pair_score`, `usefulness_k`.
  `"OCS"` under `"pairs"` is refused until the KCT rename (pair plan Phase C);
  `gbs_engine = "GA"` under `"pairs"` is refused (Phase D).
- `summary` gains **`n_families`** and **`n_parents_used`** (filled for both
  units — parents-first reports `choose(n, 2)` and `n`); `plans` gains
  **`pair_score`** (`NA` under parents-first). Appended columns only.
- 31 new tests: `test-stacking-pairs-schemes.R` (hand-computed `M` and F1
  variance on a 3 × 2 toy; the reduction identity; separability; the first set
  pick maximises `F` over every single family; loss → hard as `M` grows; the
  two levels agree with one family; shortlist vs all; both projections × both
  links × both presets × selfing; driver under both budgets and both backends)
  and the engine's stop rules in `test-stacking-pairs.R` (top-`F` = `order()`,
  the name limit with both overshoot rules, `max_use`, selfing).
- `dev/pair_selection_toy/run_toy.R`: the owner's toy comparison — one family
  × 16 offspring × 20 cycles × 5 reps on 60 simulated founders, ten arms
  (parents-first TS / OHS / OCS / GBS vs pairs TS / TS-usefulness / OHS / GBS /
  GBS-candidate / GBS-kemper); outputs beside it.

### Changed
- **`avoid_selfing` default flips to `TRUE`** (Q24 c iii). No numeric change:
  the parents-first `all_pairs` design never contained a self-cross
  (discrepancy 27, now closed). Under `"pairs"` the flag is honoured by the
  candidate pool, the budget (`allow_selfing = !avoid_selfing`) and the design
  alike, so `avoid_selfing = FALSE` offers selfs and counts them
  (`choose(n, 2) + n` families).
- `gbs_objective = NULL` resolves to `"set"` under `"pairs"` (preset value
  `"candidate"` under parents-first, unchanged).
- Parents-first path otherwise untouched: default runs identical to 0.7.1
  apart from the appended columns.

### Docs, templates, release (same day, later)
- New `vignettes/selection_grid.md` — the matrix of scheme × selection unit ×
  budget × options, with the caveats (OCS not recommended; the multi-rep RNG
  re-seed of OHS-GA / OCS). Copied to `ready_to_use/SELECTION_GRID.md`.
- `vignettes/example_stacking_simulation.Rmd`: scheme table with the OCS
  caveat, the comparison narrowed to `TS` / `OHS` / `GBS`, a new §3 *Direct
  pair selection* (both budgets, per-scheme family scores, the counts and
  scores in the records); `skill_stacking_simulation.md`, `quick_start.Rmd`
  and `simulation_stacking.Rmd` carry the same OCS note and pairs snippet.
- `ready_to_use/templates/06_stacking_simulation.R` and
  `workshop/06_stacking_simulation.R`: `SELECTION_UNIT`, `PAIR_BUDGET`,
  `OHS_METHOD` parameters, the OCS warning, family / parent counts in the
  final table. `ready_to_use/HaploSim/run_method_comparison_v2.R`:
  `MMC_SELECTION_UNIT` / `MMC_PAIR_BUDGET` route knobs (recorded in
  `metadata.dcf` and every task log; preflight refuses pairs under < 0.7.2,
  OCS or `GA_coverage` under pairs), default `MMC_SCHEMES=TS,OHS,GBS`, default
  `MMC_AVOID_SELFING=1` (the package default; under pairs `0` would offer
  selfs). Smoke-tested end to end in quick mode under pairs.
- Release bundle `ready_to_use/release_v0.7.2/`: source tarball (rebuilt with
  the docs above, 264 files, 513 KB), Windows binary built from it,
  `SHA256SUMS.txt`, `RELEASE_NOTES.md`; the test suite passes against the
  installed binary (2080 / 0 / 0).

### Found, not fixed (discrepancy 28)
- `.ga_pairwise()`, `genetic_algorithm_selection()` and `ocs_selection()` call
  `set.seed(seed)` with the driver's fixed `seed` on every cycle, so the
  crossing draws that follow are identical in every replication: under
  parents-first, **OHS (GA engines) and OCS replicates are not independent**
  (the toy run shows five identical cycle-20 values). Only `.gbs_ga_set()`
  restores the RNG (0.7.0). Campaigns that run one replicate per call with
  its own seed (the HaploSim driver) are unaffected. Fixing it (restore `.Random.seed` around the
  scheme call) changes every published parents-first OHS / OCS trajectory —
  owner's call (Q24 o).

## [0.7.1] - 2026-09-15 — Direct pair selection: the placeholder and its budget

The stacking driver has always run the **parents-first protocol**: choose a
subgroup of `top_founders`, then cross every pair of them (`make_cross_design(
"all_pairs")`, `choose(top_founders, 2)` families per cycle). The owner's next
step is **direct pairing** — choose the crosses, not the parents (plan #030
§2.4, Q23 i). This release reserves the seam for it and settles the one design
question that has to be answered before any selector is written: how a
direct-pairing cycle is matched to the parents-first run. Thread
[#031](project_record/threads/031-direct-pair-selection.md). **Nothing is
selected by pairs yet** — see *Placeholder* below. Default runs are
byte-identical to 0.7.0 (`test-stacking-pairs.R` compares an explicit
`selection_unit = "parents"` run with a default run; the full suite passes
unchanged).

### Added

- **`selection_unit = c("parents", "pairs")`** on `run_stacking_simulation()`.
  `"parents"` (default) is the implemented protocol. `"pairs"` is the direct
  pairing placeholder.
- **`pair_budget = c("families", "parents")`**, required when `selection_unit =
  "pairs"` and ignored (with a warning) otherwise. The two ways of tying a
  direct-pairing cycle to the parents-first run **conflict**, so there is no
  default — the caller has to say which one holds:
  - `"families"` — the same offspring population size: as many families as the
    traversal of `top_founders` parents makes (`choose(top_founders, 2)`); the
    number of distinct parents is then free (roughly `sqrt(2F)` to `2F`).
  - `"parents"` — the same number of parents: pairs are added until
    `top_founders` distinct parent names are in use; the number of families is
    then free (`ceiling(top_founders / 2)` to the full traversal).
- **`pair_selection_budget(top_founders, pair_budget, allow_selfing = FALSE)`**
  (exported, 128 exports): the "quick calculation" of what a given
  `top_founders` implies under each budget — `n_families` (fixed under
  `"families"`), `n_parents_max` (fixed under `"parents"`), and the
  `n_families_min` / `n_families_max` range. `allow_selfing = TRUE` adds the
  `top_founders` self-crosses to the traversal count (16 parents: 120 families,
  136 with selfing; `"parents"`: 8 to 120 families).
- Internals: `.stk_n_families()`, `.stk_selection_budget()` (resolves the
  budget for either unit; `"pairs"` without a `pair_budget` fails **before**
  anything runs), `.stk_budget_text()` (one-line statement, printed in the
  `verbose` header — `Selection unit - parents-first: 16 parents -> 120 families
  per cycle (...)`), and `.stk_select_pairs()` — the **placeholder**: the seam
  in `.stk_advance()` where a scheme would return `data.table(Parent1,
  Parent2)` for `make_cross_design(method = "custom")`. Its contract is written
  in the function comment (score candidate crosses against the ideal, stop by
  the budget; the GBS internals are already pair-agnostic).
- `tests/testthat/test-stacking-pairs.R` (12 tests): family counts, the two
  budgets and their ranges, up-front validation, the stop message carries the
  resolved budget, the warning under `"parents"`, the verbose header, and the
  default path unchanged.

### Placeholder — what `selection_unit = "pairs"` does now

Validates `pair_budget`, resolves the budget, prints it under `verbose`, then
the first cycle stops: `selection_unit = "pairs" (direct pair selection) is a
placeholder and is not implemented yet (thread #031). Budget resolved for this
run — direct pairing, pair_budget = "families": 120 families per cycle (= the
traversal of 16 parents; self-crosses not in the all-pairs design), distinct
parents free.` The compact family-ID width (`C<cycle>F<family>_<child>`) is
sized from the budget's `n_families_max` — identical to `choose(top_founders,
2)` under the default unit, so IDs do not change.

### Recorded, not fixed

- **`avoid_selfing` does not reach the all-pairs design** (INDEX discrepancy
  27, Q24 c). `make_cross_design(method = "all_pairs")` enumerates `combn(
  parents, 2)`, so a parents-first cycle never contains a self-cross whatever
  the driver's `avoid_selfing` (default `FALSE`, documented as "exclude
  self-crosses") says; the flag only acts on `"sequential"` (odd-n tail) and
  `"custom"` tables. Every run to date therefore made `choose(top_founders, 2)`
  families with no selfs. The budget is computed from what the design actually
  produces (`allow_selfing = FALSE` from the driver), so the reported family
  count is right; the `allow_selfing` term in `pair_selection_budget()` is there
  for the day the design honours the flag — a numeric change to every default
  run, hence the owner's call. `test-stacking-pairs.R` pins the current
  behaviour; the `avoid_selfing` roxygen now says so.
- Under `pair_budget = "parents"`, what happens when the next-best pair would
  overshoot the parent limit (skip it and continue, or stop) is open — Q24 b.

## [0.7.0] - 2026-09-15 — GBS: Kemper's two modes as written, and a bounded penalty scale

The genotype-building scheme (GBS) now runs Kemper et al. (2012) as published
under `gbs_preset = "kemper"` — ARCSIN or LOG deviation, sequential greedy,
fixed −100 = base −10 × scale 10 — and haploStack's adjustments are options
layered on top. Plan: `dev/UPGRADE_PLAN_opt_kemper.md` (#030) §4, Phase 1 and
its driver hooks. Thread
[#030](project_record/threads/030-gbs-kemper-modes.md). Decisions:
`dev/questions_sheet.md` Q23 b, c, d, j (owner, 2026-09-14/15). **Not in this
release:** the OPT scheme (plan §3, Phases 2–3) and the KCT rename (Phase 0) —
both wait on Q23 a / e / f.

### ⚠ Changes published numbers — GBS runs only

TS, OHS and OCS selection are byte-identical to 0.6.0 (the full suite passes
unchanged and the parity script below confirms the ARCSIN path). Three GBS
defaults move:

1. **`gbs_method = "LOG"` is Kemper's form** `d_t = τ_t · log p_t`, ranked
   ascending (a need score, like ARCSIN). 0.6.0's `"LOG"` was the 2026-06-19
   rebuild's own signed log-ratio — `log p_t − log τ_t` in the candidate score
   and `(log(p_t + 1e-10) − log τ_t)²` in D2Ideal — introduced to cure the
   origin's inverted (descending) rank. Both are **dropped, not renamed**: the
   ascending rank removes the reason for them, they have no oracle (neither
   Kemper's nor the origin's), and no reported number used them (one
   `expect_length` test). The formulas are quoted here so a 0.6.0 LOG run can be
   recomputed by hand if ever needed. Under LOG the D2Ideal metric is now
   `KL(τ ‖ p) = Σ_t w_t τ_t (log τ_t − log p_t)` with an absent state floored at
   half a copy, `1 / (2·|S|·ploidy)`.
2. **The absent-state ("diversity") penalty is `β_t = gbs_penalty · s_t`**, a
   vector over target states (owner's design). Base `gbs_penalty = −10` (the
   origin's constant). `penalty_mode = "fixed"`: `s_t = penalty_scale` for every
   state — 10 gives Kemper's −100, 1 the origin's −10. `penalty_mode =
   "dynamic"` (default): `s_t = 1 + (penalty_scale − 1) · risk_t · hope_t`,
   bounded in `[1, penalty_scale]`, with `risk_t = (1 − pop_f_t)^(n_parents ·
   ploidy)` (probability a random parental set carries no copy of the state;
   `pop_f_t` its frequency in the whole candidate pool, computed from the pool
   supplied, never read from the ideal table) and `hope_t = min(1,
   copies_in_pool_t / (τ_t · n_parents · ploidy))` (can the pool supply the
   target?). The scale rises from the origin's −10 for common states to near
   Kemper's −100 for states with a handful of copies, then declines for the
   near-hopeless, and the same pool frequency is penalised more with few parents
   than with many (`docs/haploStack_TAG_manuscript/audit/penalty_curves_2026-09-15.{R,png}`).
   0.6.0's rule — one scale-relative constant `−penalty_scale · max_t |d_t(p_min)|`
   with `penalty_scale = 2`, `p_min = 1/(n_parents·ploidy)` — is dropped;
   feeding that constant back as `gbs_penalty` with `penalty_mode = "fixed",
   penalty_scale = 1` reproduces every 0.6.0 ARCSIN pick exactly
   (`audit/gbs_parity_2026-09-15.R`: 3 fixtures × {default, effect weight 0.5,
   −10}, all identical). **`penalty_scale` changes meaning** — it multiplied
   `max|d|`, it now multiplies `gbs_penalty`.
3. **`ideal_scope = "founder"`** is the new driver default: the ideal genotype
   is built once from the founders and kept for every cycle — Kemper's design
   ("upon losing segments, ARCSIN simply tries to obtain the remaining segments
   *at the target frequency*"), and the only scope under which the D2Ideal
   trajectory is comparable across cycles. GBS selection is then also pinned to
   the founders' `select_top_blocks` so a target block cannot drop out of a later
   cycle's top-N; TS/OHS/OCS keep the per-cycle top-N and the metric is measured
   on the full offspring table, so neither is affected. `ideal_scope = "cycle"`
   reproduces 0.6.0 (rebuild every cycle; a lost state drops out and the
   remaining targets are re-normalised).

### Added

- `gbs_selection()` — `gbs_preset = c("haploStack", "kemper")`,
  `gbs_objective = c("candidate", "set")` (marginal gain of `F(S) = Σ_t w_t
  f_t(p_t(S))`, `f_t = −d_t²` for ARCSIN, `f_t = d_t` — the cross-entropy summand
  — for LOG), `gbs_engine = c("greedy", "GA")` (GA over subsets on `F(S)`,
  seeded with the greedy set; requires the set objective; restores the global
  RNG stream on exit), `penalty_mode`, `ga_control`, `seed`. A preset fills
  every knob not passed explicitly; an explicit argument always wins.
- `run_stacking_simulation()` — `gbs_preset`, `gbs_objective`, `gbs_engine`,
  `penalty_mode`, `penalty_scale`, `ideal_scope`; `gbs_penalty` and
  `gbs_effect_weight` now default to `NULL` = the preset value.
- `tune_gbs_penalty()` — `penalty_mode`, `gbs_penalty`; grids `penalty_scale`
  (the fixed scale or the dynamic cap; default grid `c(1, 2, 5, 10, 20)`) with
  the base held fixed; `best_penalty = gbs_penalty × best_penalty_scale`.
- Internals (`R/stacking_gbs.R`): `.gbs_deviation()` — **the only place the two
  formulas live**; every other function (need score, set term, D2Ideal, penalty
  scale) calls it — `.gbs_set_term()`, `.gbs_need_score()`,
  `.gbs_set_objective()`, `.gbs_marginal_scores()`, `.gbs_greedy()`,
  `.gbs_ga_set()`, `.gbs_penalty_scale()`, `.gbs_context()`,
  `.gbs_apply_preset()` + `.GBS_PRESETS`, `.gbs_check_penalty_args()`;
  `.gbs_hapn_matrix()` also returns `pop_f` and `n_pool`. The internals are
  written pair-agnostic (plan §2.4) for a later pair-level scheme.
- Tests: `tests/testthat/test-stacking-gbs-modes.R` (17 tests — plan Phase 1
  gates (a)–(h): hand-computed scores and `F(S)`; **the `kemper` preset
  reproduces an independent transcription of Kemper's sequential greedy for both
  modes** on a 30 × 25 fixture, with and without the effect weight; fixed = base
  × scale; the dynamic curve's table values, bounds and monotonicity in
  `n_parents`; absent target recruited first; the set greedy equals the
  brute-force argmax at every step; GA ≥ greedy `F(S)`; RNG restored; the
  2 × 2 × 2 × 2 smoke; D2Ideal(LOG) zero at target and ≥ 0; `ideal_scope`
  founder vs cycle and the Kemper / GA presets through the driver).
  `test-stacking-gbs.R` updated for the new penalty signature.

### Changed

- `.resolve_gbs_penalty(gbs_penalty, penalty_mode, penalty_scale, pop_f,
  target_freq, n_parents, ploidy, n_pool)` — new signature, returns a vector.
- `.stk_advance()` takes the founder ideal from `opts$founder_ideal` under the
  founder scope and pins GBS selection to `opts$founder_blocks`.
- `gbs_selection()` validates `gbs_penalty <= 0`, `penalty_scale >= 1`, and
  refuses `gbs_engine = "GA"` without `gbs_objective = "set"`.
- No export change (127).

### Verified

- Full suite under `devtools::load_all()`: **583 tests, 0 failed, 0 errors**
  (12 pre-existing skips for absent optional dependencies; the one warning is
  `test-plot-order.R`'s expected-warning test).
- `ready_to_use/haploStack_0.7.0.tar.gz` (481 KB, no `.o`/`.dll`) installed
  into a throwaway library from source: kernels live, `gbs_preset = "kemper",
  gbs_method = "LOG"` and the dynamic scale run from the installed build.

### Recorded, not fixed

- The origin GRDC `OCS_fitness_function` is **not** a sequential greedy: it
  picks all `select_top_n` at once on the initial score (every target absent →
  rank by target-copy count) and then makes one position-wise replacement
  sweep (`Origin_workflow_example/_Haplotype_Stacking/_haplotype_combin_function_GA.R`
  L1318–1416). The owner confirmed Kemper's engine is the sequential greedy,
  which 0.6.0 already implemented and 0.7.0 keeps; the origin is therefore an
  oracle for Kemper's *formula* (used as such by the tests) but not for the
  engine. INDEX discrepancy **26**; plan checkpoint K4b.

## [0.6.0] - 2026-08-24 — One representation of missing, one centring vocabulary

Two workstreams, one release, for the same reason: both change numeric output
that has already been published, and both were found by the same investigation
(`docs/block_effect_semantics_haploStack_vs_hapselect.md`). Spec:
`dev/PRD_missing_allele_normalisation.md`; plan:
`dev/UPGRADE_PLAN_missing_allele_normalisation.md`; thread
[#029](project_record/threads/029-missing-allele-normalisation.md).

### ⚠ Migration — MAN: a missing allele is `NA`, and normalisation is one-way

> **Inside a `GenoMeta`, a missing allele is `NA`. Nothing else — regardless of
> the source encoding.**

Enforced at **one portal**, `new_geno_meta()`. That is a property of the call
graph rather than a convention: `read_plink_ped()`, `read_dart_genotype()` and
`read_vcf_genotype()` all bypass `build_geno_meta()` and construct directly, so
the constructor is the only function every route passes through. Normalisation
runs after the encoding is resolved and before anything is returned.

**What changes for you.** On any panel with missing calls, read through any
reader:

- **marker effects move** — the missing code no longer ranks as an allele, so
  `fit_genomic_prediction()` stops estimating an effect for missingness (on the
  barley panel that warning fired for 8751 of 9856 markers, and 34 had their
  effect keyed to `Allele == "0"`);
- **local GEBVs at `center_by_freq = TRUE` move** — a retained token inflated
  `MarkerTotal` and deflated every real allele's frequency;
- **block signatures move**, `A0A` → `A-A`, and now respond to
  `options(haploStack.fill_na)`.

**Unchanged, and verified so on a paired fixture:** MAF, `Major`/`Minor`, call
rate, block definitions, and local GEBVs at `center = "none"`. A dosage
(`"numeric"`) container is untouched — `0` is a real genotype there and the
token map never carries it.

**Containers saved before 0.6.0** may still hold a source token.
`load_geno_meta()` says so and names `repair_geno_meta()`; it never rewrites
your file. `is_geno_meta()` warns (the container stays *valid*), and
`options(haploStack.strict_missing = TRUE)` turns that into an error.

```r
meta <- load_geno_meta("panel.rds")   # warns if it predates 0.6.0
meta <- repair_geno_meta(meta)
save_geno_meta(meta, "panel.rds")     # your call, not the loader's
```

**Normalisation is one-way.** Once a container exists, the character the source
file wrote for a missing allele is gone. `recode_geno_meta(na_as = "-")` renders
missing back out for export, but that is a new rendering, not a recovery.

### ⚠ Migration — CEN: `run_haplotype_catalogue()` no longer centres by default

`center_by_freq = TRUE` was its default and `compute_profile_localgebv()`
already defaulted to no centring, so the two entry points to the same engine
disagreed about what the default answer was. They now agree, at
`center = "none"`.

**This changes published catalogue block variances** — median factor ~3 on the
barley panel, 3.9 (range 1.8–8.8) on the release fixture
(`dev/missing_normalisation_check/cen_parity.tsv`).
**`center = "allele_freq"` reproduces the 0.5.5 numbers exactly.**

### Added

- **`use_missing_chars(chars, mode)`**, `missing_chars()`,
  `reset_missing_chars()` — register an extra missing character for the session
  when a panel marks absent calls with something the package does not know
  (`use_missing_chars("#", mode = "nucleotide")`). Read inside
  `.missing_tokens()`, so one registration reaches the readers, the
  `is_geno_meta()` check, `repair_geno_meta()` and the frequency/IBS/clustering
  paths at once. `mode` is required; one character only; never a real allele of
  that mode (A/C/G/T; the codes 1–4; **any digit** for `numeric`). Option:
  `haploStack.missing_chars`.
- **`repair_geno_meta(meta, verbose, source_format)`** — the migration path for
  containers already on disk.
- **`center = c("none", "locus_mean", "allele_freq")`** on
  `compute_profile_localgebv()`, `compute_block_catalogue()`,
  `run_haplotype_catalogue()`, `count_block_subtypes()`,
  `count_favourable_subtypes()`, `count_favourable_blocks()`,
  `analyze_block_catalogue()` and `predict_gebv()`.
  **`"locus_mean"` is new**: every allele at a locus is centred by the *shared*
  locus mean `sum_a p_a * e_a`. It preserves the allele contrast, leaves block
  variance identical to uncentred, is invariant to adding a constant to both
  alleles' effects, and reproduces HapSelect's `mean_adjust = TRUE`
  **value for value** — the identity
  `d*e_A + (ploidy-d)*e_B - ploidy*mean_i == (d - ploidy*p_A)*(e_A - e_B)`,
  checked to 1e-12 in `test-centring-modes.R`. It is the mode to use for
  selection and for any cross-package comparison.
- `metadata$missing_tokens` — the vocabulary a container was normalised
  against, so it stays self-describing in a later session. Detection and repair
  take the union of it and the session's set.
- `extra_na=` on `build_geno_meta()`, `read_genotype()`, `read_st_genotype()`,
  `read_gs_genotype()` and `new_geno_meta()` — the per-read escape hatch for a
  source whose absence code is not in the vocabulary, including a dosage panel
  that writes a digit for missing.
- `normalise_missing_alleles()` (internal) — the primitive.
- `tests/testthat/test-missing-normalisation.R` (165 assertions, V1–V5, V8–V10,
  V18–V20, and **`read_plink_ped()`'s first tests ever**) and
  `tests/testthat/test-centring-modes.R` (22, V12–V16 + C2).

### Changed

- `.missing_tokens()` gains `"0"` for `selectiontools` — a 1/2/3/4-coded file
  writes it for an absent call, never as a real allele code. Never for
  `numeric`.
- **`detect_genotype_encoding()` classifies on the calls, not the absences.**
  Cells whose every allele is a missing token are set aside before the
  vocabulary is decided. This fixes a latent misdetection: a heavily-missing PED
  used to read as `"numeric"`, and `read_plink_ped()` now declares its format
  outright so a 1/2-coded PED cannot be taken for dosage (which would have kept
  its `"0"` as a genotype).
- **`recode_geno_meta()` declares both ends.** New `from_type` defaults to the
  recorded provenance but is **validated against the data** rather than trusted,
  and a dosage source with an allele target is refused. This closes the defect
  where a container recorded as `"numeric"` had its dosages pushed through the
  nucleotide dictionary to `1 -> "A"`, `2 -> "C"`, `0 -> NA` with no warning
  (discrepancy 16, item 4). New `na_as` renders missing back out for export.
- **`metadata$encoding` is provenance, not an interpreter switch** — a record of
  the source format, not an instruction consumers read to decide what an allele
  means. A stale value is now harmless where it used to be silently wrong.
- `center_by_freq` is **deprecated** (one release) with a one-time message;
  `TRUE` → `"allele_freq"`, `FALSE` → `"none"`. Supplying both it and `center`
  is an error. Every existing call keeps working and keeps its numbers.
- Defence-in-depth (no-ops on any 0.6.0 container, correctness on a legacy one):
  `.gp_allele_dosage()` and `.meta_to_dosage_matrix()` pass the token set, and
  `.catalogue_core()`'s filter uses the shared set instead of a third
  `nzchar(trimws())` idiom. The three disagreeing missing-detection idioms are
  now one.

### Fixed

- INDEX discrepancy **23** (the missing code capturing a marker's trait effect)
  and the effect half of **18**; item 4 of **16** via `recode_geno_meta()`.

### Distribution

- `ready_to_use/haploStack_0.6.0.tar.gz` — portable source, 468 KB, no `.o` or
  `.dll`; verified by installing into a throwaway library (kernels live, 127
  exports, the portal and all three centring modes exercised from the installed
  copy).
- **Published to `CCS-voidBird/haploStack` (2026-08-25).** The mirror
  `ready_to_use/haploStack/` is committed (`8a73229` on `main`) and tagged
  **`haploStack-0.6.0`**, both pushed. The mirror is an extracted release
  tarball plus a curated `docs/` and its own README, so it was re-seeded from
  the tarball — package content is byte-identical to what shipped, except
  `README.md` and `vignettes/quick_start.Rmd`, which are deliberately
  repo-facing. `docs/CHANGELOG.md` and `docs/geno_meta_structure.md` refreshed,
  the new I/O sheet added, and the README given an upgrade warning naming both
  migration settings. The install instructions in the mirror's README and
  vignette were corrected from `uq-bioinformatics/haploStack` to
  `CCS-voidBird/haploStack`, where the release actually lives — closing the
  remaining half of INDEX discrepancy 21.
- `ready_to_use/release_v0.6.0/` — release assets: the source tarball, a Windows
  binary built *from* it (verified by running the package suite against the
  installed build: **1836 pass, 0 fail, 0 error**), and `RELEASE_NOTES.md` with
  both migration settings, SHA-256s and build provenance. **The GitHub Release
  page is not yet created** — `gh` is absent on this machine; the tag is live
  and the assets are ready to attach.
- `docs/haploStack_io_reference_0.6.0.md` — the per-release I/O sheet: what a
  missing allele is under each source mode, the new registry and repair path,
  why normalisation is one-way, what each centring mode preserves, and a
  migration checklist.

### Verification

- Full suite green with the skip and warning lists **identical to the
  pre-change baseline** at every stage; `run_examples()` exit 0.
- **Stage-3 parity table** (`dev/missing_normalisation_check/parity/`): on
  paired token/NA/dosage fixtures, every output the plan says must move did, and
  every output that must not did not — including `center = "none"` GEBVs, MAF
  and block definitions byte-identical, and a dosage container untouched
  throughout. Stage 4's four edits are byte-identical no-ops on top of stage 3.
- **V6 on the real GRDC barley panel** (1107 × 9856, 198,375 missing cells):
  with the stage-01 workaround **removed**, the container built from the raw PED
  is byte-identical to the workaround build *and* to the stored 2026-08-21
  container; 396,750 NA alleles, no token surviving
  (`dev/missing_normalisation_check/v6_barley_workaround.R`). The workaround has
  been deleted from `analysis/barley_hapselect/scripts/01_build_geno_meta.R`.
  **Not run here:** V6's cross-package half (LD, 7604 blocks, effect agreement
  r = 1.0) — HapSelect is not installed on this machine.

## [0.5.5] - 2026-07-31 — One row order for every figure

Before this release a sample order could reach exactly **one** figure:
`plot_haplotype_alignment(sample_order=)`. Every other sample-axis figure invented
its own row order, so two panels of the same population put the same line on
different rows — including through `plot_meta_alignment()`, whose purpose is to
stack panels that share an axis. 0.5.5 makes the order a **package-wide
contract** (spec `dev/PRD_meta_plot_clustering_order.md`, thread #027, decisions
`dev/questions_sheet.md` Q17).

**No default output changes.** Every figure drawn without a `sample_order` and
without a broadcast produces exactly what 0.5.4 produced; regression tests pin
the y-axis factor levels of both alignment figures.

### Added

- **`use_sample_order(x)`** — broadcast one order to every figure with a sample
  axis, by setting the new **`haploStack.sample_order`** option. Takes an
  `hs_sample_clustering` from `cluster_samples()`, a character vector of sample
  IDs, or `NULL` to clear; returns the previous value invisibly, so
  `on.exit(use_sample_order(old))` restores it. Validated **on set**, so a bad
  order fails at the call the user can see rather than inside the next figure.
- **`sample_order =` on four more figures**: `plot_localgebv_alignment()`,
  `plot_favourable_blocks()`, `plot_phenotype_bars()` and
  `plot_top_block_heatmap()`. Each accepts a character vector, an
  `hs_sample_clustering`, `TRUE` (require the broadcast; error when unset) or
  `FALSE` (ignore the broadcast for this call).
- **`plot_meta_alignment(check_order = TRUE)`** — warns when composed panels were
  drawn in different sample orders. The reference panel's order is imposed on
  every other panel, so without this the disagreement is resolved silently, which
  is the exact failure this contract exists to prevent. Panels carrying no order
  record are skipped rather than assumed to agree.
- **Provenance on the figure.** Plots carry `attr(p, "sample_order")` and
  `attr(p, "sample_order_source")` (`"default"` / `"argument"` / `"option"`), and
  a clustering-derived order adds `· rows: ward.D2 / hamming, k=4` to the
  subtitle. A reader must be able to tell a clustered panel from an alphabetical
  one without re-running the code.
- **`R/plot_order.R`** — the one resolver (`.hs_resolve_sample_order()`) plus
  `.hs_order_within()`, `.hs_warn_order_mismatch()` and the provenance helpers.

### Changed

- **`plot_localgebv_alignment()` orders rows *within* each family group**, never
  across it. The offspring/founder split and its boundary line are not the
  order's to move: they are the panel's main readable feature. Samples missing
  from an order follow the ordered ones, and haplotype copies stay adjacent and
  in index order.
- **`plot_clustered_alignment()` uses a broadcast clustering when
  `clustering = NULL`**, instead of computing a fresh one. Otherwise this panel
  would cluster independently of every other figure in the session and quietly
  disagree with them. A broadcast that is only a character order cannot drive the
  dendrogram, so that case still computes one.
- **`.resolve_sample_order()`** is now a thin wrapper over the shared contract.
  `plot_haplotype_alignment()`'s strict validation and alphabetical default are
  unchanged and pinned by tests.
- **An argument is strict; a broadcast is best-effort.** An order passed *to a
  figure* must be a permutation of that panel's samples — a mismatch is an error,
  because you named those samples for that figure. A *broadcast* order is
  intersected with whatever the panel holds, with unmatched panel samples
  appended. Otherwise setting one order would make every filtered panel in the
  session an error at a distance.
- **`plot_phenotype_bars(sample_levels=)` still outranks `sample_order`** — it is
  the explicit lock used to pin a panel to a companion figure.

### Notes

- Precedence, everywhere: per-call argument > `haploStack.sample_order` > the
  figure's own default. Same shape as `fill_na` (0.5.2), so there is one rule.
- The option is **unset** by default and `.onLoad()` deliberately does not set
  it: there is no sensible package-wide row order, only per-figure defaults.
- Not included (out of scope, `dev/PRD_meta_plot_clustering_order.md` §5):
  `plot_block_counts()` (its axis is a count, and it is produced inside
  `run_haplotype_catalogue()`), block/column-axis ordering (blocks stay in
  genomic order — that is what makes the 0.5.0 chromosome sub-panels readable),
  and tabular output, which keeps its current row order.
- New suite `tests/testthat/test-plot-order.R`.

---

## [0.5.4] - 2026-07-30 — Cross-trait analysis is raw by default

**Breaking (defaults).** Cross-trait analysis rescaled its inputs unless told not
to: `analyze_block_catalogue()` z-scored every trait's effects, and
`create_fitness_index()` divided the trait weights by `sum(abs(weights))`. Both
now use the values **as given**, and rescaling happens only when asked for
(thread #025).

### Changed

- **`analyze_block_catalogue(normalize = )` now defaults to `"none"`** (was
  `"zscore"`). Nothing is rescaled unless requested, so the effects reaching the
  cross-trait selection sum, the reported `S`/`N` tables and the plots are the
  local GEBVs as computed, in each trait's own units and comparable with
  `compute_block_catalogue()`'s `Effect` column. Under the default the `N` table
  is the identity transform of `S`.
  - **`value_source = "normalized"` counts as asking**, so it promotes
    `normalize = "zscore"` (the historical transform) with a message rather than
    returning raw numbers under a normalized label. Callers that were already
    passing `value_source = "normalized"` therefore see **no change**.
  - Naming **both** `value_source = "normalized"` and an explicit
    `normalize = "none"` is a contradiction — `N` would be `S` — and now errors
    naming the fix instead of silently mislabelling.
  - The result carries **`$normalize`**, the transform actually applied
    (`"none"`, a method name, or `"custom"` for a function).
- **`create_fitness_index(normalise_weights = )` now defaults to `FALSE`** (was
  `TRUE`), and likewise on `update_fitness_index()` and
  `run_fitness_pipeline()`. The composite is the plain `sum(weight * effect)`, so
  the index — and any GEBV computed from it — stays in raw units. Previously unit
  weights `c(A = 1, B = 1)` silently **halved** the composite. Rescaling is still
  one argument away, and because it is a positive constant factor it never
  changed any ranking or selection decision — only the magnitudes.
- **`plot_block_effect_genome()` no longer claims a standard score it did not
  compute.** Its subtitle reads "white = average effect (within-trait standard
  score)" only when a transform actually ran; a raw-default result re-plotted with
  `value_source = "normalized"` now reads "white = zero effect (raw)". Results
  from before 0.5.4 carry no `$normalize` field and are treated as rescaled, which
  is what they were.
- **`ready_to_use/` cross-trait templates** (`templates/` and `workshop/`
  `04_cross_trait_analysis.R`) switched from `zscore`/`normalized` to
  `none`/`raw`, with a comment naming when to switch back. They passed the old
  values explicitly, so the default change alone would not have reached them.

### Notes

- **What is *not* touched:** `center_by_freq` (default `TRUE`) still mean-centres
  allele effects by population frequency inside the catalogue. That is a
  within-trait part of the local-GEBV definition — it is what makes the sign of a
  block effect mean "above / below the panel mean" rather than "which allele the
  effect file called the reference" — and it is applied identically to every
  trait, so it is not a cross-trait rescaling. Changing it would move every
  catalogue number in the package away from the origin GRDC workflow.
- Anyone comparing **absolute** cross-trait or fitness-index magnitudes against
  output from 0.5.3 or earlier will see different numbers. Rankings, selections
  and block picks are unchanged except where a genuinely different scale changes
  the summed argmax — which is the point of the flag.

### Packaging — first distribution build since 0.4.9 (#026)

- **`.Rbuildignore` added** (package root). The package had **none**, so
  `R CMD build` shipped everything it found: 0.4.7 and 0.4.8 went out at **5.5 MB
  carrying `dev/` run outputs, `docs/`, the `.scratch/` issue tracker, and the
  compiled Windows objects `src/*.o` + `src/haploStack.dll`**. Shipping those
  makes a "source" tarball non-portable, which is the real defect; the size was
  just the symptom. The ignore list covers compiled artefacts, `dev/`, `docs/`,
  `sim/`, `.scratch/`, `CLAUDE.md`, `.agents/`, and scratch output
  (`Rplots.pdf`, `test_load.txt`, `vignettes/figure`).
- **`ready_to_use/haploStack_0.5.4.tar.gz`** — a portable **source** tarball
  (500 KB), sources only, so one file installs on Windows, Linux and macOS by
  compiling `src/*.cpp` on the target machine. Built with `vignettes = FALSE`
  because pandoc is absent here; the vignette *sources* ship, only rendered HTML
  is missing. Verified by installing it into a throwaway library — the C++
  recompiled, the Rcpp kernels loaded, and 0.5.3/0.5.4 behaviour was checked
  through the installed copy.
- **This closes INDEX discrepancy 14**, open since 2026-07-28: the distributed
  workflow had fallen **seven releases** behind source (newest tarball 0.4.9).
  `ready_to_use/README.md`'s install line now names 0.5.4.

### Verification

`devtools::document()` clean (`NAMESPACE` unchanged at 122) · `devtools::test()`
**1578 pass / 0 fail / 12 skip** (was 1556/0/12 at 0.5.3) ·
`devtools::run_examples()` **exit 0** · the 0.5.4 tarball installs from source
into a clean library and runs (`count_favourable_blocks()`,
`analyze_block_catalogue()` raw-default, `create_fitness_index()` unhalved).

---

## [0.5.3] - 2026-07-30 — Stacking evaluation: favourable blocks by local-GEBV sign

`run_haplotype_catalogue()` has always computed per-sample positive / zero /
negative block counts, but only as a by-product of a **full catalogue run** — it
needs an `output_dir`, writes seven files per trait, and returns the counts in
long `(Sample, Label, Count)` form at **genotype** level. There was no way to ask
the question on its own, at haplotype level, against a pedigree. This release is
that question asked directly (`docs/Update_notes.txt` item 1, thread #024, closing
`dev/qa_suggestions.md` **E4**).

### Added

- **`count_favourable_blocks()`** (`R/stacking_evaluate.R`) — one row per sample
  with the requested schema `GID, FID, IID, Label, nBlock, GEBV, nPositive,
  nNegative, nPreferred`, ordered by generation, family, then descending
  `nPreferred`. Standalone: it needs only a `GenoMeta` and marker effects,
  resolved exactly as elsewhere in the catalogue layer (explicit `effect_df`
  wins, else `marker_effects` picked by `trait`).
  - **Counting is always at haplotype level**, never genotype level, so each
    block contributes `ploidy` copies per sample and
    `nPositive + nNegative + nZero == nBlock * ploidy`. `nZero` is that
    identity's remainder rather than a column. Genotype-level counts remain
    `run_haplotype_catalogue()`'s `block_counts`; there is deliberately no
    `count_by` argument, so the two cannot drift.
  - **`prefer = c("positive", "negative")`** sets which sign is favourable, and
    therefore `nPreferred` — `nPositive` and `nNegative` are reported either
    way, so one call describes both directions.
  - **`center_by_freq = TRUE`** by default, matching `run_haplotype_catalogue()`
    rather than `compute_profile_localgebv()`'s bare `FALSE`. Uncentred, the sign
    of a local GEBV largely reflects which allele the effect file called the
    reference; centred, it means "above / below the panel mean for this block",
    which is the only reading under which counting positive blocks is a breeding
    statement.
  - **`nBlock` is per sample**, counted over cells that hold a value. Under the
    default `fill_na = 0` that is every block; `fill_na = NA_real_` keeps "no
    effect information" distinct from "zero effect", and `nBlock` then varies per
    sample while `GEBV` and the sign counts are unchanged.
  - Returns a list with fixed element names — `table`, `plots`, `dummy`,
    `prefer`, `level` — and writes `table` to `output_file` when given.
- **`plot_favourable_blocks()`** (`R/stacking_evaluate.R`) — the 1/0 dummy
  alignment: a preferred-sign copy is `1`, everything else `0`. It **calls
  `plot_localgebv_alignment()`** on an internally built dummy matrix rather than
  reimplementing a heatmap, so the chromosome separators, per-family panels,
  offspring/founder boundary and genomic block ordering are the same figure
  machinery as the continuous alignment. Accepts either a
  `count_favourable_blocks()` result or a bare dummy table. A `GenoMeta` with no
  family history **warns and returns `NULL`** instead of erroring — counting
  works without a pedigree, only the per-family figure needs one, and `plot =
  TRUE` is the default.

### Changed

- **`plot_localgebv_alignment()` gains `fill_name` (default `"Local GEBV"`)** —
  the one edit to an existing function, purely additive. A legend reading "Local
  GEBV" over a 1/0 dummy field is wrong; the default reproduces previous output
  and a regression test pins the default title.
- `utils::globalVariables()` gains the bare symbols the new file uses (11 column
  names + 3 update-join `i.*` symbols). "No visible binding" findings stay at
  **83**, with **0** from the new code.

### Notes

- **Pedigree columns are resolved to one row per sample.** A sample usually
  appears in several generations (a parent in `gen_2` is an offspring in
  `gen_1`), so the new internal `.sample_generation_map()` follows
  `.sample_family_map()`'s existing convention: the **oldest generation in which
  the sample is an offspring wins** (its birth family), falling back to the
  oldest generation in which it appears as a parent. Samples absent from the
  history — and every sample when there is no history — get
  `GID = FID = Label = NA`.
- **This is a different metric from 0.5.2's `count_favourable_subtypes()`**,
  which counts *identity matches* against catalogue-chosen subtypes. This one
  counts *local-GEBV signs*. The two tables join on `IID` / `Sample`.
- **The +/- classification is `run_haplotype_catalogue()`'s own `fcase()`
  idiom**, reused verbatim rather than refactored into a shared helper: the
  catalogue driver needs a long factor-ordered `(Sample, Label, Count)` table for
  `plot_block_counts()` while this needs wide per-sample counts, and the shared
  part is one three-branch expression. The duplication is deliberate and
  recorded, not accidental.

### Verification

`devtools::document()` clean (`NAMESPACE` 120 → 122) · `devtools::test()`
**1556 pass / 0 fail / 12 skip** (was 1448/0/12 at 0.5.2 — +108 assertions, no
regressions) · `devtools::run_examples()` **exit 0** · "no visible binding"
findings **83 → 83**, none from the new file.

---

## [0.5.2] - 2026-07-29 — Marker-effect export, subtype search, positional signatures

Two read-only ways to interrogate a GenoMeta that were missing. A fit attaches
`marker_effects`, but there was **no way to get them back out to a file** and no
marker-level effect figure (`plot_manhattan_effect()` hardcodes `"Block Effect"`
in its title). And while `compute_block_catalogue(group_by = "signature")`
already counts subtypes, it *requires* effects and has no block filter, so there
was no cheap way to ask "which haplotypes exist in this block, and how common are
they?". Both new files only read a GenoMeta; no existing function changed
behaviour.

### Added

- **`export_marker_effects()`** (`R/marker_effects.R`) — returns
  `geno_meta$marker_effects` as one stacked table and optionally writes **one
  tab-separated file per trait** (`output_dir` + `prefix`, trait name in the file
  name) or a single named file (`output_file`, one trait only; the two arguments
  are mutually exclusive). Written files carry exactly the stored columns
  `Marker, Allele, Effect`, so they round-trip through `read_effect_file()`; the
  trait lives in the file name, not an extra column. `annotate = TRUE` adds
  `Chrom`/`Pos`/`Block` to the *returned* table only. Paths written come back as
  `attr(x, "files")`. Trait names are sanitised for the file name (a trait called
  `"Yield (t/ha)"` no longer produces an invalid path) with a message naming the
  substitution; the returned `Trait` column keeps the original string.
- **`plot_marker_effects()`** (`R/marker_effects.R`) — the marker-level
  counterpart of the catalogue's block-effect scatter, sharing the same internal
  plotting body (`.plot_manhattan()`) so axis, palette and theme match
  `plot_manhattan_effect()` exactly. Because effect tables hold one row per
  `(Marker, Allele)` and most fits store the reference allele at `0`, the default
  `collapse = "max_abs"` keeps one point per marker — the largest effect in
  absolute value, **sign preserved** — instead of drawing a dense meaningless
  line of zeros beneath the signal; `collapse = "none"` plots every stored row.
  Markers absent from `map_info` are dropped with a counted warning. The plotted
  table comes back as `attr(p, "effect_table")`.
- **`count_block_subtypes()`** (`R/block_subtypes.R`) — per-block subtype
  (allele-combination) counts for one, several, or all blocks, with the fixed
  six-column schema `Trait, Block, Subtype, Count, Freq, Effect` (always all
  present, so results of several calls bind cleanly). `count_by = "haplotype"`
  (default) counts haplotype copies; `"genotype"` counts distinct samples
  carrying the subtype on at least one copy. `trait = NULL` (default) leaves
  `Trait`/`Effect` as `NA` and does **no** effect computation at all; naming a
  trait joins the block local GEBV from `compute_block_catalogue()`. Block IDs
  not present warn and are skipped rather than erroring — "no hits" is a valid
  answer for a search helper.

- **A three-layer subtype search** (`R/block_subtypes.R`), for asking whether a
  panel — usually a one-sample sub-GenoMeta — carries particular haplotypes:
  - **`count_subtype_occurrence()`** — the primitive: one block, one signature,
    how many times it occurs. A subtype that does not segregate returns `0`
    rather than erroring; an unknown block warns and returns `0`.
  - **`count_subtype_carriage()`** — the same question for a set of
    `(Block, Subtype)` targets, given as a table (any of `Subtype` / `signature`
    / `Signature`, so `count_block_subtypes()`, `write_block_allele_sheet()` and
    `compute_block_catalogue()` output all pass straight in) or a named vector.
    Returns a list: `n_target`, `n_samples`, `n_carried`, `total`, a per-sample
    `by_sample` table and a full `detail` grid in which uncarried targets remain
    as explicit zeros. One row per sample, so a one-sample sub-GenoMeta gives the
    headline count directly and a panel is ranked without a wrapper loop.
  - **`count_favourable_subtypes()`** — the pipeline: take a catalogue, pick the
    single most favourable subtype per block, count its carriage. The
    winner-per-block rule is `analyze_block_catalogue()`'s own
    `.select_block_subtype()`, so `fitness_fn` accepts the same values
    (`"max_effect"` default, `"frequent_best"`, `"expected"`, `"min_effect"`, or a
    custom function) and *favourable* means one thing across the package.
    `catalogue` accepts a `compute_block_catalogue()` table, a per-trait
    `run_haplotype_catalogue()`-shaped list, or an `analyze_block_catalogue()`
    result (whose `selected` picks are then honoured verbatim); `NULL` computes
    one from `geno_meta` + `trait`. `top_n` restricts to the strongest blocks.

### Changed

- **`.compute_block_signatures()` gained `blocks=`**, restricting the melt to
  those blocks' member markers. A single-block signature build on a
  200-sample x 2000-marker panel drops from 1.25 s to 0.42 s, and
  `count_block_subtypes()` uses it whenever no `trait` is given (with a `trait`
  the full table is still built, because `compute_block_catalogue()` labels every
  block). The absent-block warning now checks the GenoMeta's block map rather
  than the signature table, so it still fires when *every* requested block is
  missing and the restricted table comes back empty.

### Fixed

- **Block haplotype signatures are now positional — distinct allele combinations
  no longer collide.** An unobserved allele is written as a placeholder
  (`fill_na`, default `"-"`, the canonical missing character of
  `DEFAULT_ENCODING_DICT`) instead of being dropped, so every signature is
  exactly as wide as its block and each character keeps its marker's position.
  Previously, over markers `(M1, M2, M3)`, the copies `(A, T, -)` and
  `(A, -, T)` both pasted to `"AT"` and were merged into one subtype.
  - `.compute_block_signatures()` fills rather than drops, and **no longer drops
    the copy itself**: a haplotype copy with nothing observed in a block now
    appears as an all-`fill_na` subtype, so `sum(Count)` within a block is again
    `n_samples * ploidy`.
  - `write_block_allele_sheet()` uses the same placeholder in `Signature` and
    `Alleles`. It had the **same string collision** — `1;1;;4` also pastes to
    `"114"` — although its *counts* were always correct, because it groups by the
    per-marker columns rather than by the pasted string. Both producers now emit
    identical strings, so `Signature` and `count_block_subtypes()`'s `Subtype`
    join directly.
  - **Affects existing output:** signature/`Subtype` strings from
    `compute_block_catalogue(group_by = "signature")`,
    `analyze_block_catalogue()` and `write_block_allele_sheet()` change wherever
    a block carries missing data (clean data is unchanged), and blocks with
    missing data gain the previously-dropped copies. Note the container's NA-any
    rule: a genotype cell is uncalled when *any* copy is missing, so on a diploid
    one missing copy shows as `fill_na` on both copies at that marker.
  - Configurable globally with `options(haploStack.fill_na = )` or per call via
    the new `fill_na` argument on `write_block_allele_sheet()` and
    `count_block_subtypes()`. An empty placeholder is refused — it would
    reintroduce the collision.
  - New export: **`HS_FILL_NA_DEFAULT`** (`"-"`).
- **`R CMD check` "no visible binding" findings: 106 → 83** (measured with
  `tools:::.check_code_usage_in_package()` against the 0.5.1 tree and this one).
  Thirteen symbols added to `utils::globalVariables()` (`R/zzz.R`):
  - `"Freq"` and `"Subtype"` retire **5 pre-existing** findings —
    `write_block_allele_sheet()`, `.select_block_subtype()`,
    `analyze_block_catalogue()` (×2) and `plot_block_effect_manhattan()`.
    `"signature"` is registered alongside them for consistency; it produced no
    finding on its own (it is always `:=`-assigned before it is read).
  - `"i.Block"`, `"i.Chrom"`, `"i.Effect"`, `"i.Pos"` — the data.table
    update-join RHS symbols (`x[i, on =, col := i.col]`), supplied by the join
    and never lexically bound. Registering them retires **22** findings across
    `.hap_effects_rcpp()`, `.marker_signal_table()` and the catalogue joins.
  - `"S"` and `"Vbase"` (the block-subtype selection value columns) retire **4**
    — two pre-existing in `analyze_block_catalogue()` / `.select_block_subtype()`
    and the two `.resolve_top_subtypes()` would otherwise have added.
  - `"Carried"`, `"Copies"`, `"Frac"`, `"Target"`, `"Total"`, `"nCarried"`,
    `"nTarget"` for the new carriage tables.
  - Net effect: the new code contributes **0** findings of its own, and the
    pre-existing backlog shrinks by 23.

### Notes

- `count_block_subtypes()` derives both its counts and its effect-join key from a
  **single** `.compute_block_signatures()` call, which is then injected into
  `compute_block_catalogue(.signatures = )`, so the join key is identical by
  construction rather than by convention. This mattered more before the
  positional-signature fix above, when the package's two signature producers
  disagreed on missing data; they now agree, but the single-source construction
  is kept because it is the property that guarantees the join.
- `count_block_subtypes(trait = )` runs a **genome-wide** local-GEBV pass even for
  a single block, because `compute_block_catalogue()` has no block filter.
  Mitigations are named in the docs (`filter_geno_meta()` first, or reuse
  `.signatures=`). Left as-is rather than adding a filter to an existing
  hot-path function.
- The signature collision was first recorded here as INDEX discrepancy 17 and
  deliberately left unfixed; it was then **fixed in the same release** on a
  follow-up request (see *Fixed*, above). Discrepancy 17 is closed.
- `count_subtype_carriage()` answers for a whole target set in one keyed join
  rather than by looping `count_subtype_occurrence()`, because the irreducible
  per-call cost is `geno_meta_views()`. On a 200-sample x 2000-marker panel the
  full 400-block pipeline runs in **0.94 s** against roughly **2.7 min** for the
  loop. The primitive still defines the contract, and a test pins the two to
  agree target by target.
- Exports 113 → **120**. Verification: `document()` clean, `test()`
  **1448 pass / 0 fail / 12 skip** (from 1223/0/12 — +225 assertions, no
  regressions), `run_examples()` exit 0. `devtools::check()` cannot complete on
  this machine — pandoc is absent, so the vignette rebuild fails during the
  source-tarball step, before checking starts. The globals figures above were
  therefore taken with `tools:::.check_code_usage_in_package()`, which is the
  same pass that emits the NOTE.

## [0.5.1] - 2026-07-29 — One encoding vocabulary

The format module had **two names for one concept**. `"char"` / `"number"` was
the live vocabulary of the storage layer (`translate_alleles()`,
`.resolve_preset()`, `transform_effects()`); `"nucleotide"` / `"selectiontools"`
was the read-side and metadata vocabulary (`detect_genotype_encoding()`,
`GENO_FORMAT_PRESETS`, `metadata$encoding`). They met in `process_genotype()`,
which overwrites the preset's `encode` with the *detected* one — so an `ST`
preset produced `"selectiontools"` and handed it to `translate_alleles()`, which
only knew `"number"`. One vocabulary now.

### Removed

- **The `"number"` and `"char"` encoding tokens.** They are removed, not
  aliased: passing either now raises an error naming its replacement
  (`encoding "number" was removed in 0.5.1 — use "selectiontools" instead`).
  Silent normalisation was rejected as an option — it would have left the
  vocabulary split in fact while claiming it was unified. **Breaking:**
  `translate_alleles(arr, "char", "number")` and
  `build_geno_meta(encoding = "number")` no longer work.
- **`.resolve_geno_spec()`** — dead code (nothing called it) carrying a third
  copy of the retired vocabulary. Live counterpart: `resolve_geno_format()` /
  `GENO_FORMAT_PRESETS`.

### Added

- **`HS_ENCODINGS`** (`R/format_detect.R`) — the canonical vocabulary
  (`"selectiontools"`, `"nucleotide"`, `"allele"`), defined in one place.
- **`.encode_arg()`** — validates every encoding argument and always refuses a
  removed token by name, whatever the caller's `allowed` set. `strict = FALSE`
  serves `build_geno_meta()`, which legitimately accepts detector output it has
  no dictionary for.

### Fixed

- **`transform_genotype()`'s round-trip example, broken since ≤ 0.4.9**
  (project-record discrepancy 12). `devtools::run_examples()` now completes
  **exit 0** across the package; the ST → GS → ST round-trip is a true inverse.
  The discrepancy note had guessed this was "a real encoding-resolution defect
  behind the example" — it was exactly that.
- **`build_geno_meta()` could build a GenoMeta nothing could use**
  (discrepancy 15, F1): `"number"` was accepted for `char_dict` selection, then
  stored verbatim in `metadata$encoding`, which `geno_meta_views()` rejects. The
  object built cleanly and failed on first use.
- **`dev/alignment_check/run_alignment_check.R`** (discrepancy 15, F2) — stale
  and non-running: it passed a markers × samples matrix where `build_geno_meta()`
  fixes samples × markers, and used the retired tokens. Both corrected; its
  two-encoding rendering coverage is reachable again.

### Tests

- `tests/testthat/test-format.R` — three `translate_alleles` tests renamed to
  the canonical vocabulary; two added, covering the migration errors at every
  layer that used to accept the old tokens, and the contract that
  `build_geno_meta()` only ever stores an encoding `geno_meta_views()` accepts.

## [0.5.0] - 2026-07-29 — Chromosomes as alignment sub-panels

The haplotype alignment showed chromosomes as text and nothing else. A
whole-genome `plot_haplotype_alignment()` was one undifferentiated strip of
tiles with `Chr 1H,2H,…` in the subtitle, so nothing in the figure said where 2H
began — while blocks, the finer grouping, had both boundary lines and labels.
The two column-level groupings are now drawn at their own scales: chromosomes as
sub-panels, blocks annotated inside them.

### Changed

- **`plot_haplotype_alignment()` (`R/catalogue_viz.R`) draws one sub-panel per
  chromosome** — `facet_grid(cols = vars(chromosome), scales = "free_x",
  space = "free_x")`, so every panel is as wide as its share of the markers and
  is named by the strip above it. `x_pos` remains the global marker index, so
  `p$data` keeps its 0.4.9 contract (plus a new `chromosome` column) and the
  parity tests are untouched. New `chrom_panels` argument (default `TRUE`);
  `chrom_panels = FALSE` restores the single continuous marker axis and the
  `Chr …` subtitle, and a single-chromosome panel is never faceted regardless.
- **Block labels moved inside the panel, into reserved headroom.** They kept
  their vertical style but now anchor at `hjust = 1`, so the text grows
  *downward* from the panel top instead of upward past it, and
  `scale_y_discrete()` reserves `max(2, ceiling(0.1 * rows))` blank row-units for
  them. This both (a) keeps them clear of the chromosome strips and (b) fixes a
  standing defect: at `hjust = 0` the labels grew into a 10 pt `plot.margin` and
  were clipped by the device (`plot_marker_signal()` reserves 46 pt for the same
  construct).
- **Block boundaries and labels are layer data, not `annotate()`/bare
  `xintercept`.** `.alignment_block_boundaries()` and
  `.alignment_block_labels()` now return the owning `chromosome` (and, for
  boundaries, a `cross` flag) so ggplot routes each annotation to its own panel —
  `annotate()` replicates into *every* panel under faceting. Boundaries that
  straddle a chromosome break are dropped in panel mode, where the panel gap
  already states them, and kept in single-panel mode, which therefore keeps its
  exact pre-0.5.0 boundary set.
- **Natural chromosome order** via the new internal `.natural_chrom_levels()`:
  chromosomes whose names start with digits sort numerically first
  (`1H < 2H < 10H`), the rest alphabetically after. The previous
  `order(chromosome, position)` put `"10H"` ahead of `"2H"`, disagreeing with
  `.create_genomic_positions()` and the Manhattan plots. Identical to the old
  order for `1H`–`7H` and `1`–`9`.
- New arguments, all appended so positional calls are unaffected:
  `chrom_panels` (`TRUE`), `show_chrom_labels` (`TRUE`), `block_label_size`
  (`2.5`). The X scale expands additively by half a tile so edge tiles are no
  longer cut in half.

### Fixed

- **`plot_clustered_alignment()` dendrogram pairing.** The reserved block-label
  headroom shifts the alignment rows, which would have slid the dendrogram
  leaves out of register with the rows they index. `plot_haplotype_alignment()`
  now publishes its additive Y expansion as `attr(p, "row_expand")`,
  `plot_cluster_dendrogram()` gains a `leaf_expand` argument, and
  `plot_clustered_alignment()` mirrors one onto the other.

### Tests

- `tests/testthat/test-catalogue-viz.R` — six new tests: natural chromosome
  ordering; one panel per chromosome in natural order with every tile carrying
  its chromosome; `chrom_panels = FALSE` back to one panel + `Chr …` subtitle;
  single-chromosome never faceted; boundary/label layers routed per panel with
  the three cross-chromosome boundaries dropped; headroom reserved only when
  there are labels to place.
- `tests/testthat/test-catalogue-cluster.R` — the dendrogram leaf axis matches
  the alignment's expanded row range, and the unmirrored default does not.

## [0.4.11] - 2026-07-28 — Signal-marker selection

Scoring now ends in a decision. `score_marker_signal()` ranked every marker but
never picked one, so "which marker tags this block?" was re-derived by hand at
every use site — and the marker-signal plot invented its own top-quantile
candidate set at draw time, highlighting a different set than any downstream
analysis. The selection is now a first-class result that is recorded on the
`GenoMeta` and simply rendered by the plot.

### Added

- `select_signal_markers()` (`catalogue_similarity.R`): ranks the markers of
  each haplotype block by a signal metric (default `final_mcc`) and keeps the
  top `n` (default 1) as that block's **signal markers**. Ties break
  deterministically (`ties = "position"` by default, or `"maf"` / `"first"`);
  `min_score` leaves a block unselected rather than giving it a weak
  representative. Returns `(Block, Chrom, Marker, Pos, score_val, signal_rank)`.
- `attach_signal_markers()`: records a selection on `block_info` as
  `Signal_Marker` (comma-separated in rank order), `Signal_Score` (the rank-1
  marker's score) and `N_Signal`, plus the selection parameters in
  `metadata$signal_selection`. Derives `block_info` from
  `block_marker_map_info` when the container has none.
- `get_signal_markers()`: the single supported parser of the packed column —
  a tidy `(Block, Marker, signal_rank, Signal_Score)` table, or the bare marker
  vector with `long = FALSE`.
- `save_marker_signal()`: writes a marker-signal plot and its marker × block
  metrics table under one stem (`marker_signal.png` + `marker_signal.csv`).
- `score_marker_signal()` gained `n_signal` (default 1) / `select_score` and a
  fourth return element, `$signal`, carrying the selection.
- `plot_marker_signal()` gained `highlight`, `signal_markers`,
  `show_marker_labels`, `show_title` and `table_file`; every plot now carries
  its marker × block table as `attr(p, "signal_table")` (marker metrics joined
  to block, with `is_signal`, `block_score`, `block_n_markers`,
  `block_highlighted`).
- `tests/testthat/test-similarity.R`: 20 new cases covering selection, the
  container round-trip, the lifecycle guards below, and the plot defaults.

### Changed

- **`plot_marker_signal()` default output changed.** It now highlights the
  selected signal markers (`highlight = "auto"`: explicit `signal_markers` →
  the selection stored on `geno_meta` → the old quantile rule with a message),
  draws **no marker labels** (`show_marker_labels = FALSE`) and **no title**
  (`show_title = FALSE`) so the vertical block labels stay legible. Block
  labels, the block strip, the quantile line and the theoretical-maximum line
  are unchanged. *Migration:* `show_marker_labels = TRUE, show_title = TRUE`
  restores the 0.4.10 look; `highlight = "quantile"` restores the old
  highlighting rule. Supplying `title` still implies a title.
- `filter_geno_meta()` prunes signal markers that leave the panel, recomputing
  `N_Signal` and clearing a now-meaningless `Signal_Score` (messages once).
- `attach_block_info()` reports that it dropped a stored selection — redefining
  blocks invalidates a selection made against the old ones.
- `build_block_location_info()` joins the signal columns back on `Block` after
  recomputing coordinates, and writes them to `<prefix>_block_info.csv`;
  `.read_block_file()` reads them back when present.
- `docs/geno_meta_structure.md`: `block_info` / `metadata` rows updated and a
  new "Signal markers" section documents the columns and their lifecycle.
- `vignettes/ibd_ibs_signal.Rmd` and `vignettes/example_family_tracing.Rmd` show
  the score → select → attach → plot flow; `dev/r2_functional_test/run_functional_test.R`
  replaces its hand-rolled top-5% candidate CSV with the selection API.

### Fixed

- `combine_haplotypes()`'s `@examples` built a 1 × 2 matrix with only one column
  name (`dimnames = list("Mkr1", "S1")`), so `R CMD check` failed at
  *"length of 'dimnames' [2] not equal to array extent"*. Now
  `list("Mkr1", c("S1","S2"))`. Pre-existing (identical in the 0.4.9 tarball),
  found while checking this release. **`R CMD check` still errors on the next
  example in the chain** — `transform_genotype()`'s round-trip example hits
  *"unsupported encoding 'selectiontools' → 'char'"*, also pre-existing and
  unrelated to this release; see INDEX discrepancy 12.

Suite after the change: 22 files, 379 tests, 0 failures, 12 skips.

## [0.4.10] - 2026-07-27 — Sample clustering

> **Reconstructed from the code on 2026-07-28** — 0.4.10 was set in
> `DESCRIPTION` but shipped with no changelog entry and no tarball (newest build
> remains `haploStack_0.4.9.tar.gz`). Evidence: `R/catalogue_cluster.R` (new,
> mtime 2026-07-27) and the `catalogue_viz.R` edit of the same date; the
> `NAMESPACE` delta against the 0.4.9 tarball — `+export(cluster_samples)`,
> `+export(plot_cluster_dendrogram)`, `+export(plot_clustered_alignment)`,
> `+S3method(print, hs_sample_clustering)`; and `tests/testthat/test-catalogue-cluster.R`
> (11 tests) with a new `_snaps/` directory.

### Added

- `R/catalogue_cluster.R` — sample clustering from haplotype-level allele calls.
  Distances are computed on the haplotype copies but the clustering **leaves are
  whole samples**: the ploidy copies are collapsed to one genotype-level
  distance (`copy_merge = "average" | "pool" | "best"`). Marker-level Hamming
  distance with pairwise-missing exclusion, computed via a single BLAS
  `tcrossprod` over one-hot allele indicators rather than an R pair loop, so it
  scales to ~10K markers with no new dependencies.
  - `cluster_samples()` → `hs_sample_clustering` (with `samples`/`markers`/
    `blocks`/`chromosomes` pre-filtering, `linkage`, optional `k` cut).
  - `plot_cluster_dendrogram()` — dendrogram panel, horizontal or vertical,
    optionally coloured by cluster.
  - `plot_clustered_alignment()` — clustering → row-ordered alignment heatmap
    (+ dendrogram).
  - `print.hs_sample_clustering()`.

### Changed

- `plot_haplotype_alignment()` gained `sample_order`, accepting either a
  character permutation of the panel's samples or an `hs_sample_clustering`
  (its leaf order). Validated by `.resolve_sample_order()`, which errors on
  unknown/missing/duplicated IDs so a stale order cannot silently drop rows. In
  `"haplotype"` mode the ploidy copies stay grouped within each sample.

## [0.4.9] - 2026-07-22 — Alignment plot performance

### Changed

- `plot_haplotype_alignment()` (`catalogue_viz.R`): rebuilt the long-format
  data-prep so it scales to large panels. The `method = "genotype"` path no
  longer melts the full 3-D `geno_arrays` and then runs an N×M
  `by = .(Sample, marker)` group-by of `sort()` + `paste0()`; it now collapses
  the two haplotype copies to a sorted, separator-less genotype matrix in a
  single pass. Both methods replace `as.data.frame.table()` with a vectorised
  `rep()`/`as.character()` build. Plot output is **byte-identical** to 0.4.8.
  Benchmark (800 samples × 12 000 markers): genotype collapse **155.7 s → 2.9 s
  (~53×)**; haplotype long-format build ~2.4× faster.

### Added

- `collapse_genotype_sorted_cpp` (new Rcpp kernel in `src/format_cells.cpp`,
  Kernel A): sorts each cell's non-missing ploidy alleles and concatenates them
  with no separator — the compiled twin of `paste0(sort(allele), collapse = "")`.
  Wired through the existing `.hs_use_rcpp()` toggle with a vectorised pure-R
  fallback (`.collapse_genotype_sorted()` / `.collapse_genotype_sorted_r()` in
  `R/format_split.R`), matching the established kernel-dispatch pattern.
- `tests/testthat/test-catalogue-viz.R` (new): plot-data parity of both row
  methods against the original melt+group-by logic, under both engines.
  `test-format-cells.R` extended with sorted-collapse kernel-vs-pure-R parity
  (NA-drop, both-NA → `""`, multi-character sort order) and input-validation
  checks.
- Source tarball `haploStack_0.4.9.tar.gz` in `ready_to_use/`.

## [0.4.8] - 2026-07-21 — Catalogue output column order

> Reconstructed 2026-07-22 from the `ready_to_use/haploStack_0.4.7.tar.gz` vs
> `_0.4.8.tar.gz` source diff (build 2026-07-21 16:41). No `NAMESPACE` export
> change (106 → 106); the release was internal presentation only.

### Changed

- **Identification columns now lead every catalogue output table.** A new
  internal helper `.lead_identifier_cols()` (`catalogue_engine.R`, priority list
  `.HS_ID_COLS`) reorders the ID keys — `Trait, Chrom, Pos, Block, Sample,
  Haplotype, signature/Subtype …` — to the front of each returned/written
  `data.table`. Applied across `analyze_block_catalogue()`
  (`catalogue_analysis.R`: `subtypes`, `selected`, `S`, `N`, `top_blocks`),
  `compute_block_catalogue()` (`catalogue_gebv.R`), and the per-block effect
  summary in `run_haplotype_catalogue()` (`catalogue_main.R`). Column *values*
  are unchanged — only their order — so downstream keyed joins are unaffected.
- `tests/testthat/test-catalogue-analysis.R` / `test-catalogue.R` updated for the
  new column order; `man/compute_block_catalogue.Rd` /
  `man/run_haplotype_catalogue.Rd` regenerated.
- `DESCRIPTION`: version **0.4.7 → 0.4.8**.
- Source tarball `haploStack_0.4.8.tar.gz` in `ready_to_use/`.

## [0.4.7] - 2026-07-21 — Panel descriptors

Spec: [`dev/PRD_marker_stats_block_allele_sheet.md`](../dev/PRD_marker_stats_block_allele_sheet.md).
I/O guide: [`docs/haploStack_io_reference_0.4.7.md`](haploStack_io_reference_0.4.7.md)
(+ `_slide.html`).

### Added

- `marker_stats.R` (new): **`compute_marker_stats()`** — per-marker QC summary
  keyed on a `GenoMeta`, one tidy `data.table` row per marker in genomic order
  (`Marker, Chrom, Pos, Block, N_obs, Missing_rate, Call_rate, MAF, Major, Minor,
  Het`). All statistics come from a single integer-encoding pass over
  `geno_arrays` (promotes the internal `.encode_marker_codes` / `.marker_maf` /
  `.marker_allele_counts` helpers to a user-facing entry point). Optional TSV
  write via `data.table::fwrite` (default `marker_stats.txt`); the table is
  returned either way.
- `catalogue_gebv.R`: **`write_block_allele_sheet()`** — per haplotype block, the
  distinct marker-allele combinations (haplotype subtypes) present in the panel
  and their frequency, as one tidy `data.table`
  (`Block, Chrom, Signature, Markers, Alleles, Count, Freq`), ordered by block
  then descending `Count`. Ragged member-marker sets live inside `;`-joined
  `Markers`/`Alleles` cells (no per-marker columns), so the sheet scales to whole
  genomes. `count_by = c("haplotype","genotype")` switches copy vs distinct-carrier
  counting; missing alleles coded `""` so per-block `Count` sums exactly. Mirrors
  the internal `group_by = "signature"` catalogue path (`.compute_block_signatures`).
- `tests/testthat/test-marker-stats.R` (new): coverage for both functions.
- `docs/haploStack_io_reference_0.4.7.md` + `docs/haploStack_io_reference_0.4.7_slide.html`
  (new): user-facing I/O reference and a 16:9 analysis slide.
- Source tarball `haploStack_0.4.7.tar.gz` at the package root.

### Changed

- `DESCRIPTION`: version **0.4.6 → 0.4.7**.

## [0.4.6] - 2026-07-13

> Reconstructed from the `NAMESPACE` export delta (95 → 104) and file dates of
> `ready_to_use/haploStack_0.4.6.tar.gz` vs the 0.4.5 build.

### Added

- `simulate_dataset.R` (new): **synthetic-data generators** — `sim_dataset()`
  plus `sim_geno_meta()`, `sim_geno_matrix()`, `sim_map_data()`, `sim_block_map()`,
  and `sim_effect_df()`. Produce a self-contained `GenoMeta` (genotypes, map,
  blocks, effects) for demos, tests, and the ready-to-use templates.
- `update_geno_meta_from_matrix()` — rebuild a `GenoMeta`'s genotype arrays from a
  2-D matrix in place (round-trip helper).
- genomicSimulation session helpers **`gs_close()`** and
  **`gs_sync_group_names()`** exported (crossing-backend session management;
  groundwork for the deferred persistent-session work, `dev/PRD_gs_persistent_session.md`).

### Changed

- `DESCRIPTION`: version **0.4.5 → 0.4.6**.

## [0.4.5] - 2026-07-07

> The two blocks below are dated 2026-07-07 and were confirmed present in
> `ready_to_use/haploStack_0.4.5.tar.gz` (build 2026-07-07): `catalogue_analysis.R`
> carries `value_source`, and `fitness_functions.R` normalises by `sum(abs())`.
> This build also first bundles the `src/format_cells.cpp` + `src/marker_stats.cpp`
> Rcpp kernels and the `metrics`-gated `run_stacking_simulation()` driver (see the
> Rcpp / stacking-metrics notes at the end of this section).

### Fixed — Encoding & fitness (ready-to-use workflow review, 2026-07-07)

- **SelectionTools numeric codes are now alphabetical: `A=1, C=2, G=3, T=4,
  -=-1`** (`DEFAULT_ENCODING_DICT`). Previously `C` and `G` were swapped
  (`G=2, C=3`), so `recode_geno_meta(to_type = "selectiontools")` mis-coded C/G.
  `REVERSED_ENCODING_DICT` and the SelectionTools colour maps
  (`BASE_COLOR_MAP` / `COLOR_MAP`) were realigned to mirror the nucleotide
  palette under the new mapping. Nucleotide `A, C, G, T, -` remains the primary
  dictionary (workflow-review D2).
- **`create_fitness_index()` normalises trait weights by `sum(abs(weights))`**
  (was `sum(weights)`), so **signed** weights stay balanced and cannot blow up
  when they sum to ~0. A **negative** weight now correctly flips a
  lower-is-better trait; documented in `?create_fitness_index` (workflow-review
  D4).

### Changed — Block-catalogue analysis (BREAKING, 2026-07-07 revision)

- **`analyze_block_catalogue()` now selects ONE representative subtype per block,
  shared across all traits** (previously one subtype per *(block, trait)*). The
  block's subtype is `argmax_subtype Σ_t w_t · value(subtype, t)`, where `value`
  is the raw effect `S` or the within-trait normalized `N`, chosen by a new
  **`value_source = c("raw","normalized")`** argument. `value_source` drives both
  the summed selection and the plot fill, and is **stored on the result** so the
  plot layer honours it without the caller re-specifying (D2).
- **`normalize = "zscore"` is now a true standard score `(x − mean)/sd`** (centred
  and scaled) computed within each trait, replacing the previous divide-by-SD-only
  (zero-anchored) transform (D1). `minmax` / `rank` remain zero-anchored.
  Normalization now runs even for a single-trait analysis.
- **`trait_weights` is now signed direction/scale only** — applied to each
  subtype's value *before* the cross-trait selection sum (a negative weight flips a
  lower-is-better trait). It is no longer a weighted composite; the default is unit
  weights.
- **Removed the cross-trait adjustments:** the weighted composite `fitness`, the
  `consensus` table, the `corrected` column (`N − mean_t N`), and the
  `optimizing_potential` column / `.optimizing_potential()` helper. `fitness_fn`
  now only filters candidates / defines the per-subtype value (it no longer does
  the cross-trait combine).
- **Return shape changed (breaking):** `analyze_block_catalogue()` now returns
  `subtypes, selected` (one row per block), per-`(block, trait)` raw `S` and
  normalized `N` tables (each with a per-trait `is_top` flag), `top_blocks`,
  `gebv`, `plots`, and `value_source`. Gated writes drop `block_metric.txt` /
  `block_consensus.txt` and add `block_S.txt` / `block_N.txt`.
- **Plot layer (`catalogue_analysis_viz.R`):** `plot_block_effect_genome()` and
  `plot_block_effect_circos()` replace `metric = c("N","S","corrected")` with
  `value_source = c("raw","normalized")` (the `corrected` fill is gone);
  `plot_block_effect_manhattan()` / `plot_block_effect_density()` rename their
  `value` argument to `value_source`. Every plotter defaults `value_source` to the
  source stored on the result and accepts an explicit override for re-plotting.
  The diverging fill's white midpoint now means *average* effect under
  `normalized` (standard score) versus *zero* effect under `raw`.

### First bundled in this build (present in 0.4.5, added to the changelog 2026-07-21)

> These features carry no dated entry above but are **present in the 0.4.5
> tarball**; their exact introduction version (0.4.5 or an earlier 0.4.x) is not
> recoverable, so they are recorded against the earliest build that proves them.
> See `docs/project_record/threads/016-rcpp-optimization.md` and `017-stacking-metrics.md`.

- **Rcpp round 2 kernels** (`src/format_cells.cpp`, `src/marker_stats.cpp`):
  `split_genotype_cpp()` / `combine_haplotype_cpp()` (per-cell genotype split /
  combine) and `marker_allele_counts_cpp()` / `dosage_from_codes_cpp()` (per-marker
  tally / `{0,1,2}` dosage), extending the Kernel A dispatch pattern. Pure-R
  fallback retained.
- **Optional / optimised stacking metrics**: `run_stacking_simulation(metrics = )`
  (default `"gebv"`) gates the per-cycle bundle — `d2ideal`, `gene_diversity`,
  `coancestry`, and the cheap `ideal_hap_gebv` ceiling metric are computed only on
  request (measurement-only; selection and offspring are seed-identical regardless).

## [≤ 0.4.4] — earlier development history (not individually versioned)

> Everything below was accumulated across the pre-0.4.5 0.4.x line under a single
> `[Unreleased]` heading; the sources do not support splitting it into individual
> release sections. Dates on entries are their development dates.

### Performance (catalogue compute)

- `catalogue_engine.R` (new): the haplotype-catalogue compute path was
  restructured around a single **trait-invariant compute core**
  (`.catalogue_core()`). The long-format allele table, the population
  allele-frequency table, and the block/sample orderings are now built **once
  per `run_haplotype_catalogue()` call** and shared across every trait, instead
  of being rebuilt twice per trait (once at genotype level inside
  `compute_profile_localgebv()`, again at haplotype level inside
  `compute_block_catalogue()`). `compute_profile_localgebv()` and
  `compute_block_catalogue()` gained an internal `.core` argument so the run-level
  orchestrator threads the shared core into every trait (and into `future`
  workers). `compute_block_catalogue()` also gained an internal `.signatures`
  argument, and `analyze_block_catalogue()` (`.resolve_subtype_tables()`) now
  builds the core **and** the trait-invariant per-copy block-signature table once
  and shares them across all traits (previously rebuilt per trait on the
  signature + genotype-GEBV paths). Outputs are **byte-identical** to the previous
  pure-R pipeline (parity-gated by `dev/benchmark_catalogue.R`; cross-engine
  parity in `test-catalogue-engine.R`).
- `src/local_gebv.cpp` (new): **Kernel A** — `local_gebv_block_sum_cpp()`, an
  Rcpp scatter loop that sums per-block local GEBV directly from an integer-coded
  genotype matrix in one `O(n_copies × n_markers)` pass, with no long-format
  table materialised. It is **encoding-agnostic by construction**: all allele
  strings / encodings are resolved to integer codes and a numeric contribution
  table in R (`.catalogue_core()` / `.hap_effects_rcpp()`), so the kernel only
  ever indexes doubles by integer codes. Selected via the new global option
  `options(haploStack.use_rcpp = ...)` (default `TRUE`); falls back automatically
  to the pure-R `data.table` path when the option is `FALSE` **or** the kernel was
  not compiled (no compiler at install time). On `center_by_freq = TRUE` (the
  catalogue default) the per-trait block summation is ~5–6× faster than pure-R;
  results match pure-R to ~1e-13 (continuous GEBV), with all discrete structure
  (block ids, sign classification, `round(Effect, 10)` grouping, representative
  carrier) identical. Resolves dev Issue #001 (Rcpp not wired) for the catalogue
  path.
- `Rcpp` moved from `Suggests` to `Imports` + `LinkingTo`; `useDynLib` /
  `importFrom(Rcpp, sourceCpp)` wired via roxygen in `haploStack-package.R`.
- `tests/testthat/test-catalogue-engine.R` (new): core invariants, kernel input
  asserts, the `use_rcpp` toggle, and pure-R↔Rcpp parity across
  `compute_profile_localgebv()`, `compute_block_catalogue()`, and
  `run_haplotype_catalogue()`.
- `dev/benchmark_catalogue.R` (new): scalable golden-output + timing harness and
  the structure-exact / numeric-tolerant parity gate. Pre-refactor
  implementations archived verbatim in `dev/codebank/2026-06-26_catalogue_compute.R`.

### Added

- `catalogue_main.R`: `run_haplotype_catalogue()` now emits a per-block **effect
  summary** for every trait — one row per block with `BlockVariance`, `BlockMax`,
  `BlockMin` and `BlockMean` of its across-sample local GEBV, plus block location
  (`Chrom`, `Pos`). Returned as the new `block_summary` element and written to
  `Block_Summary_<ld>_t<tol>_<trait>.txt` (one file per trait). It also feeds the
  per-trait variance scatter (`plot_manhattan_variance()`), so the plotted
  variance and the on-disk table are guaranteed consistent.

### Changed

- `catalogue_viz.R`: `plot_manhattan_effect()` and `plot_manhattan_variance()`
  now share a single internal genomic-scatter body (`.plot_manhattan()`); the two
  figures differ only in the value column drawn (`Effect` vs `BlockVariance`) and
  the effect plot's dashed zero reference line. Public signatures are unchanged.

- `genomic_prediction.R`: `fit_genomic_prediction()` — bridge from a `GenoMeta`
  to the standard genomic-prediction engines. Transforms the haplo-wise
  `geno_arrays` into a numeric alt-allele dosage matrix and fits `"rrblup"`
  (`rrBLUP::mixed.solve`), `"gblup"` (`rrBLUP::A.mat` + `kin.blup`, with marker
  effects back-solved via the SNP-BLUP equivalence) or `"bayesc"`
  (`BGLR::BGLR`). Returns the estimates already shaped for re-use:
  `marker_effects` as a **named list with one `(Marker, Allele, Effect)` table
  per trait** (the effect_list shape consumed by `run_haplotype_catalogue()` /
  `compute_profile_localgebv()`), a GEBV `phenotype` table (`Sample, GEBV_<trait>`,
  drop-in for `attach_phenotype()`), the fitted `model`, and a `geno_meta` copy
  with both attached and revalidated (fitted effects merged into any already
  stored, one table per trait). `rrBLUP` / `BGLR` added to `Suggests`
  (guarded with `requireNamespace()`).
- `tests/testthat/test-genomic-prediction.R`: coverage for the dosage transform
  (`{0,1,2}` values, shape, verbatim marker names), all three model fits
  (skipped when the backend is absent), the GenoMeta round-trip through
  `compute_profile_localgebv()`, and error paths.
- `phenotype.R`: new `phenotype` element on `GenoMeta` — an optional
  `Sample` + numeric-trait data.frame, keyed on `Sample` the way `marker_effects`
  is keyed on `Marker`. Threaded through `new_geno_meta()`
  (`.standardise_phenotype()` aligns it to the genotype sample set: drop unknown
  samples, NA-fill missing rows, reorder), `is_geno_meta()` (optional validation),
  and `filter_geno_meta()` (sample-pruned). Round-trips via `save_geno_meta()` /
  `load_geno_meta()`.
- `phenotype.R`: `read_phenotype()` (delimited `.pheno`/CSV/TSV reader),
  `attach_phenotype()` (filter to shared samples + store), and `get_phenotype()`
  (accessor). Ports the data-loading half of the standalone trait-analysis script
  into the GenoMeta pipeline.
- `phenotype_viz.R`: `plot_phenotype_distribution()`, `plot_phenotype_correlation()`,
  and `phenotype_summary()` — portable (accept a `GenoMeta` or any data.frame with
  a `Sample` column) ggplot2 reimplementations of the distribution/correlation
  tasks, with no `corrplot`/`GGally`/`ggridges` dependency.
- `phenotype_viz.R`: `plot_phenotype_bars()` — per-sample dodged bar chart of up to
  two traits, with a `"vertical"` (samples on X, direct) and `"horizontal"`
  (samples on Y, merge-ready) orientation.
- `phenotype_viz.R`: `plot_meta_alignment()` — merges ggplots that share a discrete
  sample axis (e.g. a genotype-level `plot_localgebv_alignment()` panel + horizontal
  phenotype bars) into one figure: reorders every panel's sample axis to the
  reference order for identical ticks and pools the legends into a shared bottom
  row (uses `cowplot`, already in Suggests).
- `tests/testthat/test-phenotype.R`: coverage for sample intersection / NA-fill /
  warnings, validation, filter pruning, save/load round-trip, and plot outputs.

- `io/format_helper.R`: `GENO_FORMAT_PRESETS` — exported registry of genotype
  file-format *parameter sets* (`default`, `PL`, `ST`, `GS`), each describing
  on-disk orientation, header, embedded FAM columns, map column order, allele
  delimiter, and encoding. [FORMAT: PLINK/ST/GS]
- `io/format_helper.R`: `resolve_geno_format()` — dictionary resolver that maps
  a preset name (or alias, e.g. `"plink_ped"`), or a user-supplied named list,
  to a validated parameter set, applying per-field overrides. Shared by the
  reader and the writers so each format is described once. [FORMAT: PLINK/ST/GS]
- `io/format_helper.R`: `.read_genotype_files()` and `.read_map_file()` —
  internal file readers that, given a resolved parameter set, return the
  genotype (re-oriented to the fixed `sample × marker` layout), map
  (`Chrom, Marker, Pos`), and fam tables **without modifying their contents**.
- `tests/testthat/test-format.R`: tests for `resolve_geno_format()`, a PL-preset
  read with an embedded FAM block, and a preset-driven `write_geno_meta()`
  round-trip.

### Changed

- `catalogue_analysis_viz.R`: `plot_block_effect_genome()` now draws each block
  as a `geom_rect()` rectangle spanning its genomic extent
  (`Block_start`–`Block_end`, sourced from the selected-subtype table / GenoMeta
  block locations) instead of a single-point `geom_tile()`. Traits map to integer
  y-bands; point/singleton blocks get a small visible width. `geom_rect` added to
  the `ggplot2` imports.
- **Marker-effects format unified to a named list of per-trait tables.** The
  canonical shape for `geno_meta$marker_effects`, the `effect_list` of
  `run_haplotype_catalogue()`, and the `marker_effects` returned by
  `fit_genomic_prediction()` is now a **named list** (one `(Marker, Allele,
  Effect)` data.frame per trait, the list name giving the trait).
  `new_geno_meta()` / `filter_geno_meta()` / `is_geno_meta()` accept and
  normalise this list; a bare data.frame (optionally with a `Trait` column to
  split on) is still accepted for back-compatibility. Marker effects are
  **never required** for a valid `GenoMeta`.
- `catalogue_main.R`: `run_haplotype_catalogue()` — `effect_list` is now
  optional (default `NULL`) and falls back to `geno_meta$marker_effects`, so
  effects estimated by `fit_genomic_prediction()` are picked up automatically.
  `trait_names` defaults to the resolved effect traits.
- `catalogue_gebv.R`: `compute_profile_localgebv()` and
  `compute_block_catalogue()` gain a `trait` argument that selects one table
  from a stored per-trait `marker_effects` list when `effect_df` is `NULL` (a
  lone stored trait is picked automatically).
- `genomic_prediction.R`: `fit_genomic_prediction()` returns `marker_effects`
  as a one-element named list (no inline `Trait` column) and **merges** the
  fitted table into any effects already on the `GenoMeta` rather than replacing
  the whole field.
- `io/format_helper.R`: `build_geno_meta()` now enforces a **fixed input
  contract** — `geno_2d` must be `(n_samples × n_markers)` with named rows
  (samples) and columns (markers); `map` columns must be in the order
  `Chrom, Marker, Pos`; `fam` columns in the order
  `Family, Sample, Parent1, Parent2, Sex`. Orientation auto-detection and
  wide-format reshaping were removed (now handled by the reader). The function
  exposes only the changeable format parameters `split_char`, `pos_unit`, and
  `encoding`.
- `io/format_helper.R`: read-side encoding vocabulary renamed from
  `"number"`/`"char"` to `"selectiontools"`/`"nucleotide"`, matching
  `detect_genotype_encoding()` and the catalogue pipeline. The resolved encoding
  is now recorded in `metadata$encoding`. Legacy `"number"`/`"char"` remain
  accepted as silent aliases for back-compatibility.
- `io/format_helper.R`: `read_genotype()` reworked to delegate file reading to
  `.read_genotype_files()` and assembly to `build_geno_meta()`; it accepts a
  `preset`/parameter-set and the `presets` registry instead of the previous
  `row_axis`/`map_preset`/`encode` arguments.
- `io/format_helper.R`: `geno_meta_to_matrix()` and `write_geno_meta()` are now
  preset-aware (same parameter-set mechanism as the reader). They accept
  `preset`/`presets` and retain the legacy `output_format` argument
  (`"plink_ped"`/`"st"`/`"gs"`) as a deprecated alias so the Shiny GUI continues
  to work. Output orientation, allele delimiter, and map column order are taken
  from the resolved spec; the genomicSimulation write orientation is now
  marker-as-row, consistent with the GS read layout.

### Fixed

- `catalogue_profile.R`: `geno_meta_views()` no longer mangles non-syntactic
  marker IDs. The internal allele-table builder (`.bind_haplotypes()`) now calls
  `as.data.frame(..., check.names = FALSE)`, so IDs such as `"1H-123"` survive
  verbatim instead of becoming `"X1H.123"`. Previously the mangled names broke the
  `Marker` join in `compute_profile_localgebv()` (silently zeroing local GEBVs)
  and propagated mangled names into any estimated `marker_effects`.

### Migration notes

- Callers of `build_geno_meta()` must pass `geno_2d` already in
  `sample × marker` orientation. In-memory data previously relying on automatic
  orientation detection should be transposed beforehand.
- `read_genotype(..., row_axis=, map_preset=, encoding="number"/"char")` should
  migrate to `read_genotype(..., preset=)` (e.g. `preset = "ST"`), or pass a
  custom parameter-set list / `encoding = "selectiontools"|"nucleotide"`.
