# Selection grid — what can be combined with what (haploStack 0.7.2)

The stacking driver `run_stacking_simulation()` is a matrix of **scheme ×
selection unit × budget × engine/options × crossing backend**. This sheet lists
every valid cell, what it optimises, and its status. Focus on **TS, OHS and
GBS**; `OCS` has known issues (see §5).

Legend: ✔ available · ⚠ available with a caveat · ✖ refused (error) · ◻ planned

---

## 1  Scheme × selection unit

| Scheme | `selection_unit = "parents"` (default) — picks a parent subgroup of `top_founders`, crosses every pair | `selection_unit = "pairs"` — picks the families directly under `pair_budget` | Status |
|---|---|---|---|
| **TS** | top-`n` by total local GEBV `g_i` | families ranked by the **mid-parent** `(g_i + g_j)/2` (= expected F1 mean); or **usefulness** `+ usefulness_k · sd(F1)`, sd exact under the block model | ✔ / ✔ |
| **OHS** | maximise `Σ_{i<j∈S} M[i, j]` — `M` = the family's stacking ceiling (`ohs_method`: `GA_pairwise`, `Greedy`, `GA_coverage`) | families ranked by `M[i, j]`; under `"families"` the top-`F` of the triangle (exact); under `"parents"` = the best set's traversal, i.e. `ohs_selection()` is run | ✔ / ✔ (`GA_coverage` ✖ under pairs) |
| **GBS** | sequential greedy / GA on the ideal-genotype objective at the selected-set frequency `p_t(S)` | the same objective on the family's expected offspring `q_t(i,j) = (HapN_i + HapN_j)/(2m)`; `gbs_pair_level` = `"population"` (marginal gain of `F(P)` at the offspring-pool frequency, which equals `p_t(S)` for a traversal) or `"family"` (the family's own distance) | ✔ / ✔ (`gbs_engine = "GA"` ✖ under pairs) |
| **OCS** | greedy on `g_i − asin√(mean K)`, cap on the set's mean kinship; **falls back to truncation** when the cap stops it early | — | ⚠ / ✖ (no family form until the KCT rename) |
| OPT (Kemper objective 2, Meuwissen contributions) | ◻ plan #030 | ◻ contributions → cross table | ◻ |
| KCT (audited rename of OCS) | ◻ plan #030 Phase 0 | ◻ per-family kinship cap (pair plan Phase C) | ◻ |

---

## 2  The pair budget (required under `"pairs"`)

`pair_selection_budget(top_founders, pair_budget, allow_selfing)` gives the numbers.

| `pair_budget` | Held equal to the parents-first run | Free quantity | Range (16 founders) | Stop rule |
|---|---|---|---|---|
| `"families"` | offspring population size — `choose(top_founders, 2)` families (`+ top_founders` selfs when `avoid_selfing = FALSE`) | distinct parents | 120 families (136 with selfs); parents ≈ 16 … 240 | exactly `n_families` |
| `"parents"` | number of distinct parents — `top_founders` | families | 16 parents; 8 … 120 families | add families best-first until `top_founders` names are in use; then `overshoot = "skip"` (default) continues among those names up to the traversal, `"stop"` ends there |

Shared knobs under `"pairs"`: `max_use` (per-parent usage cap, default `Inf`),
`avoid_selfing` (default `TRUE`; `FALSE` offers selfs and counts them in the
budget). Under `"parents"` the flag has no effect (the all-pairs design never
selfs).

---

## 3  Per-scheme options

### TS
| Argument | Values | Unit | Note |
|---|---|---|---|
| `select_priority` | `"high"` (default) / `"low"` | both | direction |
| `ts_pair_score` | `"midparent"` (default) / `"usefulness"` | pairs | usefulness = mid-parent + `usefulness_k` × F1 sd (`c(M)` ≈ 1.03 at 4 offspring, 1.77 at 16) |
| `usefulness_k` | numeric, default `1` | pairs | |

### OHS
| Argument | Values | Unit | Note |
|---|---|---|---|
| `ohs_method` | `"GA_pairwise"` (default) / `"Greedy"` / `"GA_coverage"` | parents | GA engines re-seed the RNG each cycle — replicates in one multi-rep call are identical (§5) |
| `ohs_method` | `"GA_pairwise"` / `"Greedy"` | pairs, `"parents"` budget | runs the set engine, returns its traversal |
| `ohs_method` | irrelevant | pairs, `"families"` budget | the sort on `M` is exact |
| `ohs_method = "GA_coverage"` | — | pairs | ✖ refused (set-coverage objective, no family reading) |

### GBS
| Argument | Values | Unit | Note |
|---|---|---|---|
| `gbs_method` | `"ARCSIN"` (default) / `"LOG"` | both | Kemper's two deviation links |
| `gbs_preset` | `"haploStack"` (default) / `"kemper"` | both | kemper = fixed −100 penalty (`penalty_mode = "fixed"`, `penalty_scale = 10`), founder-fixed ideal |
| `gbs_objective` | `"candidate"` (parents default) / `"set"` (pairs default) | both | the candidate score is separable over a family (`s_i + s_j`) |
| `gbs_engine` | `"greedy"` (default) / `"GA"` (needs `"set"`) | parents | `"GA"` ✖ under pairs (planned) |
| `penalty_mode`, `gbs_penalty`, `penalty_scale` | `"dynamic"`/`"fixed"`, ≤ 0, ≥ 1 | both | `β_t = gbs_penalty · s_t`; dynamic = bounded risk × hope scale of the pool frequency |
| `gbs_effect_weight` | `[0, 1]` | both | |
| `ideal_scope` | `"founder"` (default) / `"cycle"` | both | founder-fixed ideal (Kemper) vs rebuilt every cycle |
| `construct_method`, `k_penalty` | `"custom"`/`"argmax"`/`"ratio"`, numeric | both | the ideal genotype |
| `gbs_pair_level` | `"population"` (default) / `"family"` | pairs | the two levels coincide when one family is made per cycle |
| `pair_freq` | `NULL` (by budget) / `"projected"` / `"running"` | pairs | projected = Kemper's projection to the final family count (default under `"families"`); running = mean over the families chosen so far (default under `"parents"`) |
| `pair_penalty` | `"loss"` (default) / `"hard"` | pairs | loss = exact probability `Π_f a_f(t)^M` that none of the `num_offspring` progeny carries the state; hard = the 0.7.0 absence test; loss → hard as `M → ∞` |
| `pair_pool`, `shortlist_k` | `"shortlist"` (default) / `"all"`; `NULL` = `4 · top_founders` | pairs | shortlist = families among the `shortlist_k` parents whose single gamete helps most, re-picked every step; `"all"` is exact but `~ choose(N, 2) × targets` memory per step |

### OCS (⚠ see §5)
| Argument | Values | Unit | Note |
|---|---|---|---|
| `max_kinship` | numeric, default `0.125` | parents | on the VanRaden GRM; a cap that stops the greedy early triggers a silent truncation fall-back |
| `freq_penalty_weight` | numeric | parents | no-op inside the driver (no ideal table is passed) |

---

## 4  Shared axes (every scheme, both units)

| Axis | Values | Note |
|---|---|---|
| `cross_backend` | `"block"` (default) / `"genomicSimulation"` | whole-block inheritance (no intra-block recombination) vs marker-level linkage; both units run on both |
| `top_founders`, `num_offspring`, `select_cycles`, `replications` | integers | under pairs `top_founders` is the reference the budget is derived from |
| `selfing_cycles`, `ga_start_by` | integers | intervening truncation generations; delay of the chosen scheme |
| `selected_blocks`, `select_top_blocks` | character / integer | what selection sees; crossing transmits every block |
| `metrics` | `"gebv"` (+ `"d2ideal"`, `"coancestry"`, `"gene_diversity"`, `"ideal_hap_gebv"`, or `"all"`) | measurement only; never affects selection |
| `monitor`, `n_cores`, `output_dir` | | live record; parallel reps × schemes |

Per-cycle records under both units: `summary$n_families`, `summary$n_parents_used`;
per-family `plans$pair_score` (`NA` under parents-first).

---

## 5  Caveats and known issues

- **OCS is not recommended.** It is a kinship-penalised truncation, not
  Meuwissen-1997 optimum-contribution selection nor Kemper's OPT (no
  contribution vector, a unit-dependent penalty, silent truncation fall-back
  when the cap binds — audit `docs/haploStack_TAG_manuscript/05_ocs_audit.md`).
  Its rename to `KCT` with the audit fixes, and a true `OPT` scheme, are
  planned (`dev/UPGRADE_PLAN_opt_kemper.md`). It has no family form.
- **Replicates in one multi-rep call — OHS GA engines and OCS.**
  `.ga_pairwise()`, `genetic_algorithm_selection()` and `ocs_selection()` call
  `set.seed(seed)` with the driver's fixed `seed` every cycle, so in a single
  `run_stacking_simulation(replications = k)` call their `k` replicates are
  identical (the crossing draws that follow are re-seeded). `TS`, `GBS`
  (greedy; the GBS GA restores the RNG) and `ohs_method = "Greedy"` are
  unaffected, and so are campaigns that run one replicate per call with its
  own seed (the HaploSim driver does). Recorded as project-record discrepancy 28.
- **Pairs with one family per cycle** (`top_founders = 2`, `"families"`):
  the `"population"` and `"family"` GBS levels coincide, and TS-pairs equals
  parents-first TS (the top pair is the two best parents). Population-level
  complementarity needs `top_founders ≥ 4`.
- **`"parents"` budget + additive scores** (TS, OHS): with families free the
  greedy fills the traversal of the chosen names, so the family count reaches
  `choose(top_founders, 2)` unless `max_use` caps it; for GBS (set objective)
  it stops when no family improves the objective.

---

## 6  Recommended comparison arms

| Arm | Call |
|---|---|
| parents-first baseline | `scheme = c("TS", "OHS", "GBS"), ohs_method = "Greedy"` |
| pairs, same population size | `... selection_unit = "pairs", pair_budget = "families"` |
| pairs, same parent count | `... selection_unit = "pairs", pair_budget = "parents"` |
| TS usefulness vs mid-parent | pairs + `ts_pair_score = "usefulness"` |
| GBS family vs population level | pairs + `gbs_pair_level = "family"` / `"population"` |
| GBS Kemper original vs adjusted | `gbs_preset = "kemper"` vs default, either unit |
| GBS hard vs loss penalty | pairs + `pair_penalty = "hard"` / `"loss"` |

Toy result (one family × 16 offspring × 20 cycles, `dev/pair_selection_toy/`):
OHS-pairs gave the highest gain and retained the most diversity; GBS arms gave
the lowest distance to the ideal genotype; TS-usefulness beat TS mid-parent.
