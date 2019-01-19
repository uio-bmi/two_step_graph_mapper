# Two-step graph mapper
## About

## Install
BWA mem and Minimap2 are requirements. Make sure they are installed first.
Then install by cloning:
```bash
git clone 
cd two_step_graph_mapper 
pip3 install .
```

## How to use 
First, you need some graphalignments, e.g. from running vg or rough_graph_mapper:
```bash
rough_graph_mapper map_linear_to_graph -r linear_reference.fa -f reads.fa -d graphs_dir/ --chromosomes 1,2,3 > some.graphalignments
```

The first step is to predict a path through the graph:
```bash
two_step_graph_mapper predict_path -d graphdir/ -a some.graphalignments -c 1,2,3 -o my_predicted_path
```
NB: The above command also runs bwa index, which will take an hour or two extra. You can skip that step by adding `-s`.

Map reads to that path:
```bash
two_step_graph_mapper map_to_path -r my_predicted_path.fa -f my_reads.fa -o aligned_reads.sam
```
You now have a `aligned_reads.sam` which is a sam-file relative to your predicted path. 
We can convert the coordinates to linear reference genome coordinates by running:
```bash
two_step_graph_mapper convert_to_reference_positions -s aligned_reads.sam -d graph_dir/ -l my_predicted_path -c 1,2,3 -o converted.sam
```

## Benchmarking and comparison with vg
### Quick and rough benchmarking
This test acts as an integration test, to see that everything is working and to get a quick overview of the
performance. The test takes about 30 minutes to run on a laptop, and will produce plots similar to those below:

All necessary graphs for running this test are included in this repository. 
The only thing you need is [vg](https://github.com/vgteam/vg) installed on your system.
Simply run the run_benchmarking.sh script like this, standing in the mhc_benchmark folder:
```bash
cd benchmarking/mhc_benchmark
../run_benchmark.sh None ../mhc_graph_data/linear_ref.fa None ../mhc_graph_data/wg ../mhc_graph_data/giab_chr6_haplotype0 ../mhc_graph_data/giab_chr6_haplotype1 \
       ../mhc_graph_data/giab ../mhc_graph_data/giab_reference 75 "--forward-only -n 250000 -e 0.01 -i 0.002 -l 150" 2358792 150 "" ../mhc_graph_data/ 6  ../mhc_graph_data/giab_chr6.nobg ../mhc_graph_data/giab_reference_path_6.intervalcollection.indexed ../mhc_graph_data/1000genomes_variants.vcf ../mhc_graph_data/giab_variants.vcf
```
Three plots will be produces in the folder you are running from: roc-builder.png, roc-novel-builder.png and roc-known-builder.png.
If the vg roc-plot scripts fail with some xml error, you might need to install lib-xml (on ubuntu: `sudo apt-get install libxml2-dev`)

### Thorough benchmarking on whole genome 1000 genomes graphs
This benchmark takes about 10-15 hours to run. First download all the graphs and graph indices (these are not included in this repo). 
* [Full 1000 genomes graph](http://), containing ~80m variants
* [Pruned 1000 genomes graph](http://), containing ~14m variants (only those with allele frequency >= 1%)

Download one of these and extract the contents somewhere. This is your `graph_dir`.

Create a folder inside the tests folder, e.g. `mkdir whole_genome_benchmarking`

Then run this command to run the benchmark (insert the path to your graph_dir):
```bash
graph_dir=your_graph_dir/
../run_benchmarks.sh .......
```



