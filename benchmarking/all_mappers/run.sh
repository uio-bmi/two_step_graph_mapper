#!/usr/bin/env bash
population_vcf=$1
individual_vcf=$2
individual_ID=$3
chromosomes=$4
n_chromosomes=$5
linear_ref_fasta=$6
n_threads=$7
vg_graph_base=$8
hisat2_index=$9
bwa_index=${10}

# Prepare the simulation data
graph_read_simulator_prepare_data $population_vcf $individual_vcf $linear_ref_fasta $chromosomes $n_chromosomes $individual_ID

# Simulate reads with two different error rates
coverage=0.5

# "Normal" error rate
cat haplotypes.txt | parallel --line-buffer -j 44 "graph_read_simulator simulate_reads -s 0.01 -i 0.001 -d 0.001 {} $coverage" | graph_read_simulator assign_ids positions.tsv simulated_reads.fa

# Low error rate 0.0026
cat haplotypes.txt | parallel --line-buffer -j 44 "graph_read_simulator simulate_reads -s 0.0026 -i 0.0005 -d 0.0005 {} $coverage" | graph_read_simulator assign_ids positions_low_error.tsv simulated_reads_low_error.fa

# Get number of reads
n_reads=$(wc -l positions.tsv | awk '{print $1}')
echo "Total reads: $n_reads"

# Run hisat2
echo "Running hisat2"
hisat2 -p $n_threds -f simulated_reads.fa --no-spliced-alignment -x $hisat2_index | numpy_alignments store sam hisat2 $n_reads
hisat2 -p $n_threds -f simulated_reads_low_error.fa --no-spliced-alignment -x $hisat2_index | numpy_alignments store sam hisat2_low_error $n_reads

# Run bwa untuned
bwa-mem2 mem -t $n_threads $bwa_index simulated_reads.fa | numpy_alignments store sam bwa_untuned $n_reads

# Run bwa tuned
map_linear $n_threads $bwa_index $linear_ref_fasta simulated_reads.fa | numpy_alignments store sam linear $n_reads

# Run vg
echo "Running vg"
map_vg $n_threads $vg_graph_base simulated_reads.fa && cat vg.pos | numpy_alignments store pos vg $n_reads
map_vg $n_threads $vg_graph_base simulated_reads_low_error.fa && cat vg.pos | numpy_alignments store pos vg_low_error $n_reads


