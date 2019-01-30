id=builder
threshold=150
mappers=$1
timestamp=$(date +%s)

mkdir $timestamp
echo "Results will end up in directory $timestamp"

echo "correct mq score length unaligned known.nodes known.bp novel.nodes novel.bp is_on_rare_variant aligner" > results.tsv
for mapper in $(echo $mappers| tr "," "\n")
    do
    echo $mapper
    cat $mapper.compare | awk -v name="${mapper}" 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, name}' >> results.tsv
done
../plot-roc.R results.tsv

( cat results.tsv | head -1; cat results.tsv | tail -n+2 | awk '{ if ($8 == 0) print }' ) | gzip >results-known.tsv
( cat results.tsv | head -1; cat results.tsv | tail -n+2 | awk '{ if ($8 > 0) print }' ) | gzip >results-novel.tsv

../plot-roc.R results.tsv $timestamp/roc-$id.png
echo rendering known ROC
../plot-roc.R results-known.tsv $timestamp/roc-known-$id.png
echo rendering novel ROC
../plot-roc.R results-novel.tsv $timestamp/roc-novel-$id.png

exit

( cat vg-pan-se.compare | awk 'BEGIN { OFS="\t"; print "correct", "mq", "score", "length", "unaligned", "known.nodes", "known.bp", "novel.nodes", "novel.bp", "is_on_rare_variant", "aligner"; } { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "vg" }' ;
  #cat two_step_graph_mapper.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "two-step-using-linear" }' ;
  #cat two_step_graph_mapper_traversemapped.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "two-step-traversemapper" }' ;
#  cat two_step_graph_mapper_vg.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "two-step-vg" }' ;
#  cat vg-pan-se-mitty.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "vg-mitty" }' ;
  cat seven_bridges.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "Seven-Bridges" }' ;
#  cat seven_bridges_mitty.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "Seven-Bridges-Mitty" }' ;
#  cat vg-pan-se-1pc.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "vg-pruned-graph" }' ;
#  cat bwa-untuned.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11,"BWA-MEM-untuned" }'
  cat bwa.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11,"BWA-MEM-tuned" }'
  ) | gzip >results-$id.tsv.gz

( zcat results-$id.tsv.gz | head -1; zcat results-$id.tsv.gz | tail -n+2 | awk '{ if ($8 == 0) print }' ) | gzip >results-known-$id.tsv.gz
( zcat results-$id.tsv.gz | head -1; zcat results-$id.tsv.gz | tail -n+2 | awk '{ if ($8 > 0) print }' ) | gzip >results-novel-$id.tsv.gz
#( zcat results-$id.tsv.gz | head -1; zcat results-$id.tsv.gz | tail -n+2 | awk '{ if ($10 > 0) print }' ) | gzip >results-rare-$id.tsv.gz

# This can then be rendered using scripts in the vg repo
echo rendering ROC
../plot-roc.R results-$id.tsv.gz roc-$id.png
echo rendering known ROC
../plot-roc.R results-known-$id.tsv.gz roc-known-$id.png
echo rendering novel ROC
../plot-roc.R results-novel-$id.tsv.gz roc-novel-$id.png
#../plot-roc.R results-rare-$id.tsv.gz roc-rare-$id.png
#echo rendering QQ
#../plot-qq.R results-$id.tsv.gz qq-$id.pdf


