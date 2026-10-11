# Rhynchospora OrthoFinder compatibility repair

The original comparison job completed but OrthoFinder 2.5.5 failed during DendroBLAST tree distance calculations with `TypeError: only 0-dimensional arrays can be converted to Python scalars`. Its generic worker failure handler described this as running out of RAM. The traceback does not establish an OOM. In the upstream 2.5.5 source, `mins` and `maxes_inv` are (n,1) NumPy arrays, so indexing by one integer produces a one-dimensional value passed to a scalar multiprocessing buffer assignment. Modern NumPy rejects that conversion.

GENESPACE 1.3.1 can proceed with preexisting ordinary `Orthogroups.tsv` and search files, and disables useHOGs when N0.tsv is absent. Thus the completed original figure is a useful ordinary-orthogroup comparison; the upstream phylogenetic completion claim needed correction.

The repair creates a separate Python 3.11 / NumPy 1.26.4 / SciPy 1.13.1 / OrthoFinder 2.5.5 environment. It copies cached search files into an isolated continuation directory and uses upstream `-b` to avoid rerunning the nine DIAMOND searches. No `-o` is passed: upstream 2.5.5 only permits that flag with new `-f` runs. It then exposes the completed result to GENESPACE with local copied search files, since GENESPACE expects them in the immediate WorkingDirectory.

The subprocess must exit zero. All input peptides must occur exactly once across assigned and unassigned ordinary orthogroups, with correct species membership. Rooted species-tree and nonempty hierarchical-group results must exist and contain the expected species and peptide IDs. The previous ordinary group membership is compared by sorted gene sets, ignoring arbitrary OG labels. Search/input hashes, versions and commands are saved.

GENESPACE is rebuilt from the same prepared peptide/BED inputs and successfully completed OrthoFinder results, explicitly using useHOGs=FALSE to keep the ordinary orthogroup choice of the first figure. Newly available orthologue information can still change downstream results; identical drawings are not assumed. The corrected comparison uses the same physical chromosome lengths, reference chromosome palette, regions and whole-chromosome moves. Original outputs remain available unchanged. The new archive is `local_results/SynTangle_rhynchospora_native_repaired.zip` and includes the repair audit and subprocess log.

Run `bash validation/real_data/submit_rhynchospora_native_repair.sh`. Environment setup is a small nano job; analysis uses four cores and 32G, matching the earlier allocation rather than increasing RAM. No partial archive is promoted if repair/comparison fails.

Sources:
- https://github.com/davidemms/OrthoFinder/blob/2.5.5/scripts_of/orthologues.py
- https://github.com/davidemms/OrthoFinder/blob/2.5.5/scripts_of/__main__.py
- https://github.com/davidemms/OrthoFinder#using-pre-computed-blast-results
- GENESPACE 1.3.1 run_genespace.R and set_syntenyParams.R.
