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

## [Unreleased]

_The remaining roadmap items (stacking-simulation rebuild #013, genomicSimulation
persistent session GSPS) are re-targeted to **v0.5.x** and tracked in
[`docs/project_record/INDEX.md`](project_record/INDEX.md). An opening PRD for the
remaining input formats (allelic DArT, VCF) is at
`dev/PRD_input_formats_dart_vcf.md` — questions only, nothing scoped._

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
