# Two-step graph mapper

## Install
BWA-MEM2 and Minimap2 are requirements. Make sure they are installed first.

If you want to run experiments comparing Two Step Graph Mapper with other graph-based read mappers, see the "Benchmarking" section further down.

Assuming you have BWA-MEM2 (not the "standard" BWA-MEM1) and Minimap 2 installed, you can install Two-step Graph Mapper by cloning and using pip:
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
# Map linear first
bwa-mem2 mem linear_reference.fa reads.fa > mapped.sam
# Fit linear alignments to graph
rough_graph_mapper mdz_align_bam -b mapped.sam -d graphs_dir/ -o rough_graph_alignments -c 1,2,3
```


After having some rough alignments, the first step of the two-step approach is to predict a path through the graph:
```bash
two_step_graph_mapper predict_path -d graphs_dir/ -a rough_graph_alignments -c 1,2,3 -o my_predicted_path --input-is-edgecounts True
```
NB: We send the rough_graph_alignments as input, but the above command can also take a vg json file as input. 
In that case, `--input-is-edgecounts` should be set to false, since the input is full alignments.
NB: The above command also runs bwa index, which will take an hour or two extra. You can skip that step by adding `-s`.

You now have a linear reference path through the graph, and can map reads to that with BWA-MEM 2:
```bash
bwa-mem2 mem my_predicted_path.fa my_reads.fa > aligned_reads.sam
```

You now have a `aligned_reads.sam` which is a sam-file relative to your predicted path. 
We can convert the coordinates to linear reference genome coordinates by running:

```bash
cat aligned_reads.sam | two_step_graph_mapper convert_to_reference_positions -d graph_dir/ -l my_predicted_path -c 1,2,3 > converted.sam
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

You are now ready to run the experiments. 

### Reproducing the results from the manuscript: Thorough benchmarking on whole genome 1000 genomes graph
This benchmark takes several hours to run. First download all the graph data and graph indices. 
* [Pruned 1000 genomes graph](https://zenodo.org/record/2586090/files/human_pruned_1pc.tar.gz?download=1), containing ~14m variants (only those with allele frequency >= 1%)
* [Graphs representing the GIAB sample that we simulate reads from](https://zenodo.org/record/2586090/files/simulation_data.tar.gz?download=1)

Download both of these and extract the contents somewhere. The first one is your `graph_dir`, the second one is your `simulation_data_dir`

You need to make a BWA-MEM2 index of the linear reference genome (this index is not included in the Zenodo repository due to file size limitations, but is straight-forward to create):
```bash
cd graph_dir/
bwa-mem2 index hg19_chr1-Y.fa
```

Init some variables:
```bash
graph_dir=/path/to/downloaded/graph/dir/
simulation_data_dir=/path/to/downloaded/simulation_data/dir/
chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X


```

There are four different experiments, each with a different folder inside the benchmarking folder of this repository.

### Experiment 1 (figure 1, 2 and 3): Comparing vg, Hisat2 and Seven Bridges
Note: You will have to run Seven Bridges manually and make the seven_bridges.npz and seven_bridges_low_error_rate.npz if you want to include Seven Bridges in the figures.
To skip Seven Bridges, simply modify the commands for creating the figures below to not include Seven Bridges.
```bash
cd benchmarking/all_mappers
./run.sh /simulation_data_dir/variants.vcf.gz HG002 $chromosomes 23 $graph_dir 60 $graph_dir/wg $graph_dir/hisat2_index $graph_dir/hg19_chr1-Y.bwa-mem2.fa $graph_dir

# You can now produce an html report which includes figure 1:
numpy_alignments make_report -f figure1 --names="vg (high read error rate),vg (low read error rate), Seven Bridges (high read error rate), Seven Bridges (low read error rate), Hisat2 (high read error rate), Hisat2 (low read error rate)"  truth vg,vg_low_error,seven_bridges,seven_bridges_low_error,hisat2,hisat2_low_error "#3355FF,#79B4FF,#BC35C2,#D797DA,#148046,#5BC38B"

# Similar for figure 2 and 3:
numpy_alignments make_report -f figure2 --names="vg, BWA-MEM untuned, Linear mapper"  truth vg,bwa_untuned,linear "#3355FF,#D68B8B,#B81B1B"
numpy_alignments make_report -f figure3 --names="vg, Seven Bridges, Hisat2, Linear mapper"  truth vg,seven_bridges,hisat2,linear "#3355FF,#BC35C2,#D68B8B,#B81B1B"
```

### Experiment 2 (figure 4):
```bash
cd benchmarking/two_step_approach_vg
./run.sh /simulation_data_dir/variants.vcf.gz HG002 $chromosomes 23 $graph_dir 60 $graph_dir/wg $graph_dir/hg19_chr1-Y.bwa-mem2.fa $graph_dir/

# Make figure 4
numpy_alignments make_report -f figure4 --names="vg, Linear mapper, Two-step approach using vg-alignments" truth vg,linear,two_step_approach_vg "#3355FF,#B81B1B,#de9000"
```

### Experiment 3 (figure 6) 
Note: This experiment does a whole genome 30x read simulation. Running vg alone on this data takes about 28 hours using 24 threads. The whole experiment can take a couple of days to run 
and will need several hundre gigabytes of free disk space.
```bash
cd benchmarking/two_step_approach
./run.sh /simulation_data_dir/variants.vcf.gz HG002 $chromosomes 23 $graph_dir 60 $graph_dir/wg $graph_dir/hg19_chr1-Y.bwa-mem2.fa $graph_dir/

# Make figure 6
numpy_alignments make_report -f figure6 --names="vg, Linear mapper, Two-step approach" truth vg,linear,two_step_approach "#3355FF,#B81B1B,#89D2D9"
```

### Experiment 4 (figure 7)
In this experiment, we simulate with two different read depths, and copy the results back to the two_step_approach directory to compare with those simulations.
```bash
# 15x
cd benchmarking/two_step_approach_15x/
./run.sh /simulation_data_dir/variants.vcf.gz HG002 $chromosomes 23 $graph_dir 60 $graph_dir/wg $graph_dir/hg19_chr1-Y.bwa-mem2.fa $graph_dir/
cp two_step_approach_15x ../two_step_approach/.

# 7.5x
cd benchmarking/two_step_approach_7.5x/
./run.sh /simulation_data_dir/variants.vcf.gz HG002 $chromosomes 23 $graph_dir 60 $graph_dir/wg $graph_dir/hg19_chr1-Y.bwa-mem2.fa $graph_dir/
cp two_step_approach_7.5x ../two_step_approach/.

# Now we can make figure 7:
numpy_alignments make_report -f figure7 --names="Two-step approach 30x, Two-step approach 15x, Two-step approach 7.5x, vg" truth two_step_approach,two_step_approach_15x,two_step_approach7.5x,vg "#6E2C00,#BA4A00,#EB984E,#3355FF"
```


### Experiment 5 (variant calling analysis)
We have made a separate github repository to describe this analysis. Follow the instructions there:
https://github.com/uio-bmi/variant_calling_benchmarks_two_step_mapper


### After running the experiments
You are now inside the docker container, and will probably want to export the html reports. You can do this with `docker cp`:
```bash
exit  # Exit the docker container
docker cp [container ID]:/two_step_graph_mapper/benchmarking/[experiment]/[report directory id] .
```
Replace [container ID] with the ID of the container (run `docker ps -a` to see all containers with ID). Replace [report directory ID]
with the directory of the report, e.g. `figure3`. After copying, you should be able to open the report in a web browser and see the figures.




