#!/usr/bin/env bash
population_vcf=$1
individual_vcf=$2
individual_ID=$3
chromosomes=$4
n_chromosomes=$5
linear_ref_fasta=$6
n_threads=$7
vg_graph_base=$8
bwa_index=$9
obg_graph_dir=${10}

# Prepare the simulation data
graph_read_simulator_prepare_data $population_vcf $individual_vcf $linear_ref_fasta $chromosomes $n_chromosomes $individual_ID

# Simulate reads with two different error rates
coverage=3.75

cat haplotypes.txt | parallel --line-buffer -j 44 "graph_read_simulator simulate_reads -s 0.01 -i 0.001 -d 0.001 {} $coverage" | graph_read_simulator assign_ids positions.tsv simulated_reads.fa

# Get number of reads
n_reads=$(wc -l positions.tsv | awk '{print $1}')
echo "Total reads: $n_reads"

# Store truth positions as numpy alignments
cat positions.tsv | numpy_alignments store truth truth $n_reads

# Run two step
two_step_map $chromosomes $n_threads /home/ivar/data/human_1pc/ simulated_reads.fa $n_reads 2 $bwa_index $linear_ref_fasta
cat two_step_final.sam | numpy_alignments store sam two_step_approach_7.5x $n_reads
numpy_alignments set_correctness truth two_step_approach_7.5x


wait

echo "Done"
