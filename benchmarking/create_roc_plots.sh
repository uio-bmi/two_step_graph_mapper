id=builder
threshold=150


echo combining results
( cat bwa.compare | awk 'BEGIN { OFS="\t"; print "correct", "mq", "score", "length", "unaligned", "known.nodes", "known.bp", "novel.nodes", "novel.bp", "is_on_rare_variant", "aligner"; } { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "BWA-MEM" }' ;
  cat two_step_graph_mapper.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "two-step-graph-mapper-linear" }' ;
  cat two_step_graph_mapper_traversemapped.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, "two-step-graph-mapper-traversemapper" }' ;
  cat vg-pan-se.compare | awk 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11,"vgpan" }'
  ) | gzip >results-$id.tsv.gz

( zcat results-$id.tsv.gz | head -1; zcat results-$id.tsv.gz | tail -n+2 | awk '{ if ($8 == 0) print }' ) | gzip >results-known-$id.tsv.gz
( zcat results-$id.tsv.gz | head -1; zcat results-$id.tsv.gz | tail -n+2 | awk '{ if ($8 > 0) print }' ) | gzip >results-novel-$id.tsv.gz
( zcat results-$id.tsv.gz | head -1; zcat results-$id.tsv.gz | tail -n+2 | awk '{ if ($10 > 0) print }' ) | gzip >results-rare-$id.tsv.gz

# This can then be rendered using scripts in the vg repo
echo rendering ROC
../plot-roc.R results-$id.tsv.gz roc-$id.png
echo rendering known ROC
../plot-roc.R results-known-$id.tsv.gz roc-known-$id.png
echo rendering novel ROC
../plot-roc.R results-novel-$id.tsv.gz roc-novel-$id.png
../plot-roc.R results-rare-$id.tsv.gz roc-rare-$id.png
#echo rendering QQ
#../plot-qq.R results-$id.tsv.gz qq-$id.pdf


