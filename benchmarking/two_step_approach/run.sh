#!/usr/bin/env bash
individual_vcf=$1
individual_ID=$2
chromosomes=$3
n_chromosomes=$4
linear_ref_fasta=$5
n_threads=$6
vg_graph_base=$9
bwa_index=$10
obg_graph_dir=${11}

# Prepare the simulation data
graph_read_simulator_prepare_data $population_vcf $individual_vcf $linear_ref_fasta $chromosomes $n_chromosomes $individual_ID

# Simulate reads with two different error rates
coverage=15

cat haplotypes.txt | parallel --line-buffer -j 44 "graph_read_simulator simulate_reads -s 0.01 -i 0.001 -d 0.001 {} $coverage" | graph_read_simulator assign_ids positions.tsv simulated_reads.fa

# Get number of reads
n_reads=$(wc -l positions.tsv | awk '{print $1}')
echo "Total reads: $n_reads"

# Store truth positions as numpy alignments
cat positions.tsv | numpy_alignments store truth truth $n_reads

# Run bwa tuned
map_linear $n_threads $bwa_index $linear_ref_fasta simulated_reads.fa | numpy_alignments store sam linear $n_reads

# Run vg
echo "Running vg"
map_vg $n_threads $vg_graph_base simulated_reads.fa && cat vg.pos | numpy_alignments store pos vg $n_reads

# Run two step
two_step_map $chromosomes $n_threads /home/ivar/data/human_1pc/ simulated_reads.fa $n_reads 5 $bwa_index $linear_ref_fasta
cat two_step_final.sam | numpy_alignments store sam two_step_approach $n_reads

#Set correctness
numpy_alignments set_correctness truth linear &
numpy_alignments set_correctness truth vg &
numpy_alignments set_correctness truth two_step_approach &

wait

echo "Done"
