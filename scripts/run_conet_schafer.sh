#!/bin/sh
set -eu

root_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cache_dir="$root_dir/analysis_cache/schafer_phyllosphere_2022/conet"
jar="$root_dir/external_data/conet/conet-1.1.1.beta.jar"
main=be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser
methods=correl_pearson/correl_spearman/dist_bray/dist_kullbackleibler/sim_mutInfo

mkdir -p "$cache_dir"
python3 "$root_dir/scripts/prepare_conet_input.py" \
  "$root_dir/cleaned_data/schafer_phyllosphere_2022_abundance.csv" \
  "$cache_dir/abundance_minocc20.tsv" --min-occurrence 20

java -Xmx8g -cp "$jar" "$main" \
  --input "$cache_dir/abundance_minocc20.tsv" --matrixtype count \
  --method ensemble --ensemblemethods "$methods" \
  --thresholdguessing edgeNumber --guessingparam 500 --topbottom \
  --output "$cache_dir/thresholds_minocc20.txt" --verbosity fatal

java -Xmx8g -cp "$jar" "$main" \
  --input "$cache_dir/abundance_minocc20.tsv" --matrixtype count \
  --method ensemble --ensemblemethods "$methods" \
  --ensembleparamfile "$cache_dir/thresholds_minocc20.txt" \
  --networkmergestrategy union --multigraph --filter rand \
  --randroutine edgeScores --resamplemethod shuffle_rows --renorm \
  --iterations 100 --edgethreshold 0.05 \
  --randscorefile "$cache_dir/permutation_minocc20_100.txt" --scoreexport \
  --format gml --output "$cache_dir/permutation_minocc20_100.gml" --verbosity fatal

java -Xmx8g -cp "$jar" "$main" \
  --input "$cache_dir/abundance_minocc20.tsv" --matrixtype count \
  --method ensemble --ensemblemethods "$methods" \
  --ensembleparamfile "$cache_dir/thresholds_minocc20.txt" \
  --networkmergestrategy union --multigraph --filter rand \
  --randroutine edgeScores --resamplemethod bootstrap --iterations 100 \
  --edgethreshold 0.05 --pvaluemerge brown --multicorr benjaminihochberg \
  --nulldistribfile "$cache_dir/permutation_minocc20_100.txt" \
  --randscorefile "$cache_dir/bootstrap_minocc20_100.txt" --scoreexport \
  --format gml/tab_table --output "$cache_dir/conet_minocc20_100" --verbosity fatal

/Users/newuser/opt/miniconda3/bin/python3 "$root_dir/scripts/evaluate_conet_schafer.py"
