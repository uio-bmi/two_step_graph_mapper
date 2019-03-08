# Two-step graph mapper

## Install
BWA-MEM and Minimap2 are requirements. Make sure they are installed first.

If you want to run experiments comparing Two Step Graph Mapper with other graph-based read mappers, see the "Benchmarking" section further down.

Assuming you have BWA-MEM and Minimap 2 installed, you can install Two-step Graph Mapper by cloning and using pip:
```bash
git clone https://github.com/uio-bmi/two_step_graph_mapper.git 
cd two_step_graph_mapper 
pip3 install .
```

## How to use 
First, you need some graphalignments, e.g. from running vg or rough_graph_mapper. 
If using the rough_graph_mapper, you will need to specify the chromosomes that you are running on after `--chromosomes`.
The `graphs_dir/` is a directory containing graphs. You can test with any of the two graph directories linked to under the Benchmarking
 section below. 

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
The following explains how to run the benchmarks presented in the manuscript. 

**NOTE**: Seven Bridges is not directly assessed with these scripts, since a license is required in order to run the Seven Bridges mapper. 
For including Seven Bridges in these tests, [get a license](http://sevenbridges.com/graph-genome-academic-release) and run Seven Bridges 
using the sim.fq and reads_mitty.fq files produced by first running this benchmark. Then, simply put the sam files produces
 by the Seven Bridges mapper in the same directory that you run the bencmarks from, and they will be included
  in the comparison. The files must be called `seven_bridges.sam` and `seven_bridges_mitty.sam` to automatically
  be included.

## Initial setup
There are many software dependencies required for comparing all the mapping tools, and installing all of them is a tedious process.
Thus, we have created a Docker image that contains everything that is needed. To get started, run the following commands
to install the Docker image (this should take about 30 minutes):
```bash
git clone https://github.com/uio-bmi/graph_mapping_benchmarking.git
cd graph_mapping_benchmarking
docker build -t graph_mapping_benchmarking .
```

When the image is built, run the image interactively and change directory to this repository
(which automatically is included in the container):
```bash
docker run -it graph_mapping_benchmarking
cd two_step_graph_mapper/benchmarking  # Go inside this repository inside the container, and you are ready to run the benchmarks
```

You are now ready to run the experiments. There are two experiments available. The first ("Quick and rough benchmarking") is quick to run (approximately 30 minutes),
and can be run on a normal computer. The second is a full benchmark which requires about 50 GB of memory and takes many hours to run. 

### Quick and rough benchmarking
This acts as an integration test to see that everything is working and to get a quick overview of the
performance of the different mapping methods. The test takes about 30 minutes to run on a laptop.

Assuming you are inside the Docker container (see guide above), simply download the MHC graphs and run the run_benchmarking.sh script like this
 (make sure you are positioned in the two_step_graph_mapper/benchmarking directory before running this):
```bash
# Be positoned in the benchmarking directory
wget -O mhc_graph_data.tar.gz https://zenodo.org/record/2586090/files/mhc_graph_data.tar.gz?download=1 && tar -xzf mhc_graph_data.tar.gz
# We now create a dedicated directory for running the mhc benchmarks:
mkdir mhc_benchmark && cd mhc_benchmark
# Now, the benchmarks are run with this single command. Copy and paste everything below:
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
This benchmark takes several hours to run. First download all the graph data and graph indices. 
* [Pruned 1000 genomes graph](https://zenodo.org/record/2586090/files/human_pruned_1pc.tar.gz?download=1), containing ~14m variants (only those with allele frequency >= 1%)
* [Graphs representing the GIAB sample that we simulate reads from](https://zenodo.org/record/2586090/files/simulation_data.tar.gz?download=1)

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
    $graph_dir 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X None None None None 20 63025520 $graph_dir/hisat2_index None None
```

If successfully run, you will end up with a lot of `.compare` files, one for each read mapper. 
You can now generate ROC plots for those you want. See the "Generating ROC plots" section below. 

## Generating ROC plots
These are the commands used to generate the figures in the manuscript.
```bash
# You should be positioned in the directory that you ran the benchmarking from, e.g. "mhc_benchmark" for the mhc benchmarking
# Figure 1
../create_roc_plots.sh vg,vg_mitty,seven_bridges,seven_bridges_mitty,hisat,hisat_mitty

# Figure 2
../create_roc_plots.sh vg,bwa,bwa_untuned

# Figure 3
../create_roc_plots.sh vg,seven_bridges,hisat,bwa

# Figure 4
../create_roc_plots.sh vg,two_step_graph_mapper_vg,bwa

# Figure 5
../create_roc_plots.sh vg,bwa,seven_bridges,hisat,two_step_graph_mapper_linearmapped
```
When running one of those commands, the final ouput should be something like `Report created at 1551899104/report.html`.
This means you now have a html report with figures, which you can open in a web browser. 
The number is a timestamp and will be unique each time. If you ran the bencmarking inside a Docker repo, you 
typically cannot open the html report inside the container. In that case, you can copy the report out:

```bash
exit  # Exit the docker container
docker cp [container ID]:/two_step_graph_mapper/benchmarking/mhc_benchmark/[report directory id] .
```
Replace [container ID] with the ID of the container (run `docker ps -a` to see all containers with ID). Replace [report directory ID]
with the directory of the report, e.g. 1551899104. After copying, you should be able to open the report in a web browser.




