# Four original Leptidea assemblies

This pipeline uses the recovered original male assemblies LsinapisSpaM, LsinapisSweM,
LjuvernicaM and LrealiM from DOI 10.1007/s10577-023-09713-z. It reuses completed
Liftoff mappings and genome downloads. No substitute species enter the comparison.

One valid longest translated isoform is selected for each quality-screened unique
transferred donor gene. Protein sequence searches, rather than equal donor IDs,
define discovery matches. Six independent reciprocal BLASTP/MCScanX pairwise runs
use five distinct best subjects per query and MCScanX maximum gene gap 10.
Three adjacent pairs are scored; all six pairs provide direct reference anchors.

Gene coordinates are converted from zero-based half-open to MCScanX one-based
inclusive, then back to zero-based half-open for SynTangle. Whole chromosomes
retain their full recovered lengths; chromosomes with no adjacent scored blocks
are listed explicitly in provenance. This is newly inferred evidence from
transferred annotations, not recovery of the original MAKER annotations or exact
published ribbons. Input order is the original numbered assembly order, not a
verified reconstruction of the published panel's order. Colours use Spanish
L. sinapis as reference through direct block overlap.

Submit `bash validation/real_data/submit_leptidea_paper_comparison.sh` from
SynTangle_pipeline_test. Defaults use recovered source_recovery_CRtNNj/leptidea,
the already installed annotated_lep_rjv857/tools environment, and gffread from
rhynchospora_native_mNeWPw/tools. Override SCT_LEPTIDEA_TOOLS and
SCT_LEPTIDEA_GFFREAD if those cached tools have moved. Submission checks required
executables before creating jobs. Setup translates proteins on one CPU; all six
pair jobs run concurrently with four CPUs each. The figure job uses one CPU,
a 300-second SynTangle budget, and the existing GENESPACE R environment.

Logs: st_lep_protein, st_lep_paper_pair and st_lep_paper_fig. The latest archive
is local_results/SynTangle_Leptidea_paper_comparisons.zip. It is only promoted
when the complete comparison returns successfully. A complete comparison does
not imply proof of optimality: consult result.json for bounds and status.
