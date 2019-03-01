id=builder
threshold=150
mappers=$1
timestamp=$(date +%s)

mkdir $timestamp
echo "Results will end up in directory $timestamp"

echo "correct mq score length unaligned known.nodes known.bp novel.nodes novel.bp is_on_rare_variant aligner" > $timestamp/results.tsv
for mapper in $(echo $mappers| tr "," "\n")
    do
    echo $mapper
    cat $mapper.compare | awk -v name="${mapper}" 'BEGIN { OFS="\t"} { print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, name}' >> $timestamp/results.tsv
done

( cat $timestamp/results.tsv | head -1; cat $timestamp/results.tsv | tail -n+2 | awk '{ if ($8 == 0) print }' ) | gzip >$timestamp/results-known.tsv
( cat $timestamp/results.tsv | head -1; cat $timestamp/results.tsv | tail -n+2 | awk '{ if ($8 > 0) print }' ) | gzip >$timestamp/results-novel.tsv

../plot-roc.R $timestamp/results.tsv $timestamp/roc-$id.png
echo rendering known ROC
../plot-roc.R $timestamp/results-known.tsv $timestamp/roc-known-$id.png
echo rendering novel ROC
../plot-roc.R $timestamp/results-novel.tsv $timestamp/roc-novel-$id.png

cd $timestamp
python3 ../../make_html_report.py $mappers > report.html
cd ..

echo "Report created at $timestamp/report.html"


