haploStack annotation parity fixture (decision-form row 36 b; design D5-10)

Synthetic markers laid on the public MorexV3 (Ensembl Plants) gene rows of
chromosomes 1H and 2H. No project genotype data is included.

  genes_1H2H.gff3.gz         the 1H + 2H gene rows of the MorexV3 GFF3, verbatim
  snp_map.txt.gz             2 400 synthetic markers (name, chr, pos), 1 200 per
                             chromosome, uniform over the gene span, seed 20261111
  block_marker_map.csv.gz    1 195 blocks of 1-5 consecutive markers (Chrom, Marker, Block)
  block_info.csv.gz          their spans (Chrom, Block, Start, End)
  step*.tsv.gz               the Python block_annotation 0.3.0 local route on these
                             inputs (run_pipeline, use_ensembl = False, chr_map
                             1 -> 1H, 2 -> 2H, deep annotation, no xref), re-cut once
                             on 2026-11-11; python_run.log is its log

tests/testthat/test-block-annotation-parity.R reproduces the five step files
line for line with compute_block_annotation(flank_bp = 0) + export_block_annotation().
The generator and the re-cut driver are in .scratch/maintenance-0.9.0/parity/
(make_parity_fixture.R, recut_reference.py).
