# Troubleshooting

First: `packageVersion("haploStack")`, `is_geno_meta(meta, verbose = TRUE)`,
`dim(meta$geno_arrays)`. Reproduce on a small clip or `sim_dataset()` before
reading anything big. Errors carry classes (`class(e)` in a `tryCatch`):
`hs_error_input`, `hs_error_block_map`, `hs_error_format`, `hs_error_schema`,
`hs_error_dependency` (install command in `e$hint`).

| Symptom | Cause -> fix |
|---|---|
| `... is not a valid GenoMeta: <reason>` | read the reason: counters vs array, map order, a marker in two blocks, newer schema. Rebuild through the readers rather than editing fields by hand |
| `hs_error_block_map` on a block file | header not `Chrom, Marker, Block` / `Chrom, Block, Start, End`, or no header -> rename columns; a block on two chromosomes; overlapping blocks -> `overlap = "first"`; chromosome names differ -> `chrom_map =` |
| warning: fewer than half the markers matched the block map | marker ID schemes differ (`1_123` vs `chr1:123`) |
| far fewer markers than the input | genotype and map are intersected on `Marker`; ID schemes differ |
| GEBVs all zero / flat | effect `Marker` IDs or `Allele` labels do not match the panel: `length(intersect(eff$Marker, meta$map_info$Marker))`, `table(meta$geno_arrays)` vs `unique(eff$Allele)` |
| effect file has fewer than 3 columns | `read_effect_file()` needs `Marker, Allele, Effect` (separator and header are detected) |
| `The file is ragged` | lines really differ in field count: `table(count.fields(f))` finds the odd one |
| rows that are all `NA` | separator-only lines (`,,`) at the end of a spreadsheet CSV; drop rows whose ID is `NA` |
| `must have row names ... and column names` | `build_geno_meta()` wants rows = samples, columns = markers, both named |
| integer codes refused / `hs_error_format` naming `st_convert()` | SelectionTools codes (`1..4`, `-1`) -> `use_st_formats(TRUE)` and `st_convert()`, or declare the encoding |
| encoding `"numeric"` | `0/1/2` dosage; convert to alleles first (ask which allele was counted) |
| `SelectionTools formats are off` | `use_st_formats(TRUE)`, or prefer `hapS_ld_blocks()` / your own block map |
| `hs_error_dependency` (e.g. HapSelect, rrBLUP, genomicSimulation, circlize) | install what `e$hint` says, or take a route that does not need it |
| thousands of `the standard deviation is zero` warnings | `HapSelect::pairwise_ld()` called directly on invariant markers; `hapS_ld_blocks()` muffles and counts them |
| SelectionTools and HapSelect blocks disagree | expected - different engines; compare within one engine |
| every block a singleton on simulated data | `sim_geno_meta()` simulates independent markers; LD blocking finds nothing |
| message on load about an upgrade | container from an older haploStack, upgraded in memory; re-save to keep the new form |
| `hs_error_schema` on load | container written by a newer haploStack; upgrade the package |
| catalogue numbers differ from an old run | centring mode or package version differs; recompute rather than compare |
| output files named `_0.5_t0_` though blocks came from elsewhere | `ld_threshold` / `tolerance` on the catalogue are labels; pass the real values |
| `no family-history individual matches the genotype sample set` | pedigree IDs differ from sample IDs, or the FAM has no parents; family figures will be empty |
| `plot_manhattan_effect: need columns Chrom, Pos, Effect` | it takes a marker table; for a container use `plot_marker_effects(meta, trait =)` |
| `plot_simulation_metrics: 'summary' needs columns ...` | pass `sim$summary`, not `sim` |
| `"pairs" needs 'pair_budget'` / OCS has no pair form | set `pair_budget`; run OCS under `"parents"` |
| OCS ran though nobody asked | an old call named it, or a pre-0.9.0 build (the default is now TS / OHS / GBS) |
| `the effects hold N traits ... pass trait =` | `effect_df =` one table or `trait = "<name>"` (0.9.0 no longer takes the first trait silently) |
| `select_priority = "low" is not supported by scheme` | TS honours it; GBS only with `construct_method = "argmax"`; OHS / OCS never |
| `top_founders (n) exceeds the population` | lower `top_founders` (0.9.0 no longer caps it silently) |
| no files after the run | `output_dir` defaults to `NULL` since 0.9.0; name a directory or `write_outputs(sim, dir)` |
| `genomicSimulation holds one session per R process` | `gs_close()` the open handle, or `replace = TRUE`; one gS job per process |
| warning `<code>: ...` at the end of a run | a noted fallback; `hs_notes(sim)` has the rows (job, first generation, count) |
| simulation trajectory flat or falling | `select_priority` inverted, or `effect_df = NULL` took another trait |
| `summarise_simulation: no metric column with a value` | the trajectories table has no `MeanGEBV` / `SDev` values (pass `sim$trajectories`) |
| `deprecated` warning | the message names the replacement; `?"haploStack-deprecated"` lists every retired name |
| stray `Rplots.pdf` in the working directory | `options(device = function(...) grDevices::pdf(NULL))` in scripts |
| progress messages flood the log | wrap the call in `hs_quiet()` |
| `... have no gene on that seqid ... A chr_map that would match: ...` | the container and the GFF3 spell chromosomes differently (`1` vs `1H`); pass the suggested `chr_map =` (container -> annotation). A QTL table warns only when no chromosome matches at all |
| `... positions this small look like Mb or cM` on annotation | `map_info$Pos` is not bp (`metadata$pos_unit` unset); rebuild the container with `pos_unit = "bp"` positions |
| `defunct` error (`hs_error_defunct`) | the name is past its warning release: `define_haplotype_blocks*()` -> `st_ld_blocks()` / `hapS_ld_blocks()`; `center_by_freq = TRUE/FALSE` -> `center = "allele_freq"` / `"none"` |
| `greedy_selection()` / `genetic_algorithm_selection()` warn | deprecated 0.9.0: `ohs_selection(method = "Greedy")` / `ohs_selection()`, or `run_stacking_simulation(scheme = "OHS")` |
| `output_format =` warns | the old spelling of `preset =` (`"PL"`, `"ST"`, `"GS"`) on `write_geno_meta()` / `geno_meta_to_matrix()` |
| BayesC left files in the working directory (<= 0.8.8) | since 0.9.0 the chain files are private and removed; pass `saveAt =` to keep them |
| SelectionTools: `st.input.dir` changed by a call | since 0.9.0 the bridge restores your values and never changes the working directory; if you set them yourself, they are back after the call |
| a reused R worker behaves differently from a fresh one | `hs_reset_session()` between jobs (clears the per-job settings, closes a genomicSimulation session, removes the package's scratch); the contract is `?haploStack_server` |

