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

cd $timestamp
python3 ../../make_html_report.py $mappers > report.html
cd ..

echo "Report created at $timestamp/report.html"


