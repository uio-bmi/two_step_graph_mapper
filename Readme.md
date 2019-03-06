# Two-step graph mapper

## Install
BWA-MEM and Minimap2 are requirements. Make sure they are installed first.

For running the experiments described further down on this page, also [vg](http://github.com/vgteam/vg) and 
[Mitty](https://github.com/sbg/Mitty) (Seven Bridges' tool for simulating reads) are needed.

Assuming you have BWA-MEM and Minimap 2 installed, you can install Two-step Graph Mapper by cloning and using pip:
```bash
git clone https://github.com/uio-bmi/two_step_graph_mapper.git 
cd two_step_graph_mapper 
pip3 install .
```

## How to use 
First, you need some graphalignments, e.g. from running vg or rough_graph_mapper. 
If using the rough_graph_mapper, you will need to specify the chromosomes that you are running on after `--chromosomes`.
The `graphs_dir/` is a directory containing graphs. You can test with the directory `benchmarking/mhc_graph_data/` in this repository, which contains
one test-graph for chromosome 6. You will also find a link to a full human graph under the Benchmarking section further down.

```bash
rough_graph_mapper map_linear_to_graph -r linear_reference.fa -f reads.fa -d graphs_dir/ --chromosomes 1,2,3 > some.graphalignments
```


After having some rough alignments, the first step of the two-step approach is to predict a path through the graph:
```bash
two_step_graph_mapper predict_path -d graphs_dir/ -a some.graphalignments -c 1,2,3 -o my_predicted_path
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

# Benchmarking
The following explaind how to run the benchmarks as presented in the manuscript. 

**NOTE**: Seven Bridges is not directly assessed with these scripts, since a license is required in order to run the Seven Bridges mapper. 
For including Seven Bridges in these tests, [get a license](http://sevenbridges.com/graph-genome-academic-release) and run Seven Bridges 
using the sim.fq and reads_mitty.fq files produced by first running this benchmark. Then, simply put the sam files produces by the Seven Bridges mapper in this directory (they must be called `seven_bridges.sam` and `seven_bridges_mitty.sam`), and they will be included in the report.

## Initial setup
There are many software dependencies required for comparing all the mapping tools, and install all of them is a tedious process.
Thus, we have created a Docker image that contains everything that is needed. To get started, install this image by cloning the
repository containing the Docker file that will build the image, and build the Docker image:
```bash
git clone https://github.com/uio-bmi/graph_mapping_benchmarking.git
cd graph_mapping_benchmarking
docker build -t graph_mapping_benchmarking .
```
Building the image should take 15-20 minutes.

When the image is built, you can enter the container, and you will then be ready to run the experiments (see the two different experiments below):
```bash
docker run -it graph_mapping_benchmarking
cd two_step_graph_mapper   # Go inside this repository inside the container, and you are ready to run the benchmarks
```

### Quick and rough benchmarking
This test acts as an integration test, to see that everything is working and to get a quick overview of the
performance. The test takes about 30 minutes to run on a laptop.

Assuming you are inside the Docker container (see guide above), simply download the MHC graphs and run the run_benchmarking.sh script like this, standing in the mhc_benchmark folder:
```bash
# Be positoned in the benchmarking directory
wget http://158.39.75.109/mhc_graph_data.tar.gz && tar -xzf mhc_graph_data.tar.gz
# Go to the benchmark directory and run the benchmarks
mkdir mhc_benchmark && cd mhc_benchmark
# Copy and paste everything below:
../run_benchmark.sh None ../mhc_graph_data/linear_ref.fa None ../mhc_graph_data/wg ../mhc_graph_data/giab_chr6_haplotype0 \
       ../mhc_graph_data/giab_chr6_haplotype1 ../mhc_graph_data/giab ../mhc_graph_data/giab_reference 75 \
       "--forward-only -n 250000 -e 0.01 -i 0.002 -l 150" 2358792 150 "" ../mhc_graph_data/ 6  \
       ../mhc_graph_data/giab_chr6.nobg ../mhc_graph_data/giab_reference_path_6.intervalcollection.indexed \
       ../mhc_graph_data/1000genomes_variants.vcf ../mhc_graph_data/giab_variants.vcf.gz 6 4970557 \
       ../mhc_graph_data/hisat2_index ../mhc_graph_data/graph_minimap_index ../mhc_graph_data/numpy_graph.npz
```

If successfully run, you will end up with a lot of `.compare` files, one for each read mapper. 
You can now generate ROC plots for those you want. See the "Generating ROC plots" section below.

### Reproducing the results from the manuscript: Thorough benchmarking on whole genome 1000 genomes graph
This benchmark takes about 10-15 hours to run. First download all the graph data and graph indices (these are not included in this repo). 
* [Pruned 1000 genomes graph](http://158.39.75.109/human_pruned_1pc.tar), containing ~14m variants (only those with allele frequency >= 1%)
* [Graphs representing the GIAB sample that we simulate reads from](http://158.39.75.109/simulation_data.tar)

Download both of these and extract the contents somewhere. The first one is your `graph_dir`, the second one is your `simulation_data_dir`

Create a folder inside the benchmarking folder of this repository, e.g. `mkdir whole_genome_benchmarking`

Then run this command to run the benchmark (insert the paths to the two directories):
```bash
graph_dir=your_graph_dir/
simulation_data_dir=your_simulation_data_dir/
threads=32  # Number of threads you want to use

# Run the benchmarks
../run_benchmark.sh None $simulation_data_dir/hg19_chr1-Y.fa None $graph_dir/wg $simulation_data_dir/haplotype0_only_chr20_no_paths $simulation_data_dir/haplotype1_only_chr20_no_paths \
    $simulation_data_dir/giab $simulation_data_dir/giab_only_reference $threads "--forward-only -n 5000000 -e 0.01 -i 0.002 -l 150" 2358792 150 "" \
    $graph_dir 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X None None None None 20 63025520 
```

If successfully run, you will end up with a lot of `.compare` files, one for each read mapper. 
You can now generate ROC plots for those you want. See the "Generating ROC plots" section below. 

## Generating ROC plots
These are the commands used to generate the figures in the manuscript.
```bash
# You should be positioned in the directory that you ran the benchmarking from, e.g. "mhc_benchmark" for the mhc benchmarking
# Figure 1
../create_roc_plots.sh vg,vg_mitty,seven_bridges,seven_bridges_mitty

# Figure 2
../create_roc_plots.sh vg,bwa,bwa_untuned

# Figure 3
../create_roc_plots.sh vg,two_step_graph_mapper_vg,bwa

# Figure 4
../create_roc_plots.sh vg,bwa,two_step_graph_mapper_vg

# Figure 5
../create_roc_plots.sh vg,bwa,seven_bridges,two_step_graph_mapper_traversemapped,two_step_graph_mapper_linearmapped
```





