#!/bin/bash
set -e
# Partly copy from vg's benchmark script (https://github.com/vgteam/vg/blob/master/scripts/map-sim)

#sudo docker pull quay.io/vgteam/vg:v1.12.1
#alias dockervg="sudo docker run quay.io/vgteam/vg:v1.12.1 vg"
    #sudo docker run quay.io/vgteam/vg:v1.12.1 vg "$@"
#}

if [ $# -ne 24 ];
then
    echo "Invalid input parameters. Check Readme for examples on how to run."
    exit
fi

output=$1
fasta=$2
ref=$3
pan=$4
haps1=$5
haps2=$6
sim=$7
sim_ref=$8
threads=$9
read_spec=${10}
seed=${11}
threshold=${12}
vg_map_opts=${13}
obg_graph_dir=${14}
chromosomes=${15}
simulation_ob_graph=${16}
simulation_ob_reference_path=${17}
graph_vcf_file=${18}
simulation_vcf_file=${19}
simulation_chromosome=${20}
simulation_chromosome_size=${21}
hisat2_index=${22}
graph_minimap_index=${23}
ob_numpy_graphs=${24}

pan_xg=$pan.xg
pan_gcsa=$pan.gcsa
ref_xg=$ref.xg
ref_gcsa=$ref.gcsa
haps1_xg=$haps1.xg
haps2_xg=$haps2.xg
sim_xg=$sim.xg
sim_ref_xg=$sim_ref.xg
vg_scripts_dir=/home/ivar/dev/vg/scripts/

echo $sim_xg $ref_xg $pan_xg

mkdir -p $output

# Get the vg id
id=$(vg version | cut -f 3 -d- | tail -c 8 | head -c 7)

# generate the simulated reads if we haven't already
if [ ! -e sim.gam.truth.tsv ];
then

    # Simualate reads using Mitty (Seven Bridges way of doing it)
    echo "$simulation_chromosome 0 $simulation_chromosome_size . +" > region.bed
    echo "Filter variants"
    mitty filter-variants $simulation_vcf_file HG002 region.bed mitty_filtered_variants.vcf
    echo "Gzip"
    bgzip -f mitty_filtered_variants.vcf
    tabix -p vcf mitty_filtered_variants.vcf.gz
    echo "Simulate mitty reads"
    mitty -v4 generate-reads --threads $threads --unpair $fasta mitty_filtered_variants.vcf.gz HG002 region.bed hiseq-2500-v1-pcr-free.pkl 24 1 reads_mitty.fq reads_mitty.fq.txt
    gzip -c reads_mitty.fq > reads_mitty.fq.gz
    # Add sequencing errors to the reads
    mitty corrupt-reads --threads 40 hiseq-2500-v1-pcr-free.pkl reads_mitty.fq.gz reads_mitty_with_errors.fq reads_mitty.fq.txt reads-corrupt-1q.txt 1
    # Get correct alignment position for each reads (seems we have to extract this from the fq name). We are only interested in reads that are novel, we add 1 novel node to all these (vg roc script only cares about this column)
    grep '^@' reads_mitty_with_errors.fq | awk -F ":" '/148=/{print $0 " $simulation_chromosome " $3 " 148 0 0 148 0 0 0"} !/148=/{print $0, " $simulation_chromosome " $3 " 148 0 0 148 1 1 0"}' | sed -e 's/@//g' | sort > sim.gam.truth.mitty.tsv

    # Simulate using vg
    echo generating simulated reads
    vg sim $read_spec -s $seed -x $haps1_xg -a >sim1.gam &
    vg sim $read_spec -s $(echo "$seed + 1" | bc) -x $haps2_xg -a >sim2.gam &
    wait
    cat sim1.gam sim2.gam >sim.gam

    rm -f sim1.gam sim2.gam
    vg annotate -p -x $sim_xg -a sim.gam | vg view -a - | jq -c -r '[ .name, .refpos[0].name, .refpos[0].offset ] | @tsv' | sort >truth.tsv
    vg annotate -n -x $sim_ref_xg -a sim.gam | tail -n+2 | sort >novelty.tsv
    # Find out which reads cover rare variants
    #python3 ../find_rare_variants.py $simulation_ob_graph sim.json $simulation_ob_reference_path $graph_vcf_file $simulation_vcf_file | sort > rare_variants.pos
    #join truth.tsv novelty.tsv | join - rare_variants.pos > sim.gam.truth.tsv
    # Add zeros instead of rare variant information:
    join truth.tsv novelty.tsv | awk '{print $0, " 0"}' > sim.gam.truth.tsv

    # split the file into the mates
    vg view -X sim.gam | gzip >sim.fq.gz &
    vg view -a sim.gam | jq -cr 'select(.name | test("_1$"))' | vg view -JaG - | vg view -X - | sed s/_1$// | gzip >sim_1.fq.gz &
    vg view -a sim.gam | jq -cr 'select(.name | test("_2$"))' | vg view -JaG - | vg view -X - | sed s/_2$// | gzip >sim_2.fq.gz &
    wait

    # Convert fq to fasta, neede for hybrid aligner
    gunzip -c sim.fq.gz > sim.fq
    #fastq_to_fasta -i sim.fq -o sim.fa
    paste - - - - < sim.fq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > sim.fa

fi


# Map these simulated reads in three different ways:

# Method 2: Graph minimap
rough_graph_mapper remove_reads_from_fasta -f sim.fa -a linear_to_graph_mapped.graphalignments > sim_without_linear_to_graph_mapped.fa
graph_minimap sim_without_linear_to_graph_mapped.fa $graph_minimap_index $ob_numpy_graphs $threads graph_minimap.graphalignments.tmp
cat linear_to_graph_mapped.graphalignments graph_minimap.graphalignments.tmp > graph_minimap.graphalignments

# 1) Two step mapper
# Run the two different methods for initial rough mapping
# Method 1:
rough_graph_mapper map_linear_to_graph -t $threads -r $fasta -f sim.fa -d $obg_graph_dir -c $chromosomes -o linear_to_graph_mapped.graphalignments
two_step_graph_mapper predict_path -t $threads -d $obg_graph_dir -a linear_to_graph_mapped.graphalignments -c $chromosomes -o predicted_path
two_step_graph_mapper map_to_path -t $threads -r predicted_path.fa -f sim.fa -o two_step_graph_mapper.sam
# convert to linear pos
two_step_graph_mapper convert_to_reference_positions -s two_step_graph_mapper.sam -d $obg_graph_dir/ -l predicted_path -c $chromosomes -o two_step_graph_mapper_on_reference.sam
awk '$2!=2048 && $2 != 2064' two_step_graph_mapper_on_reference.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort >two_step_graph_mapper.pos
join two_step_graph_mapper.pos sim.gam.truth.tsv | ../vg_sim_pos_compare.py $threshold >two_step_graph_mapper_linearmapped.compare


# Hisat 2
hisat2 -p 30 -q sim.fq --no-spliced-alignment -x $hisat2_index > hisat.sam
awk '$2 < 256' hisat.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort >hisat.pos
join hisat.pos sim.gam.truth.tsv | ../vg_sim_pos_compare.py $threshold  > hisat.compare

# Hisat 2 using Mitty reads
hisat2 -p 30 -q reads_mitty_with_errors.fq --no-spliced-alignment -x $hisat2_index > hisat_mitty.sam
awk '$2 < 256' hisat_mitty.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort >hisat_mitty.pos
join hisat_mitty.pos sim.gam.truth.mitty.tsv | ../vg_sim_pos_compare.py $threshold  > hisat_mitty.compare

# Method 2, traversemapper
#rough_graph_mapper traversemapper -t 70 -r $fasta -f sim.fa -d $obg_graph_dir -c $chromosomes -o traversemapped.graphalignments
#rough_graph_mapper filter --min-mapq 50 -a traversemapped.graphalignments > traversemapped.filtered.graphalignments
#two_step_graph_mapper predict_path -t $threads -d $obg_graph_dir -a traversemapped.filtered.graphalignments -c $chromosomes --linear-ref-bonus 1 -o predicted_path_traversemapped
#two_step_graph_mapper map_to_path -t $threads -r predicted_path_traversemapped.fa -f sim.fa -o two_step_graph_mapper_traversemapped.sam
#two_step_graph_mapper convert_to_reference_positions -s two_step_graph_mapper_traversemapped.sam -d $obg_graph_dir/ -l predicted_path_traversemapped -c $chromosomes -o two_step_graph_mapper_on_reference_traversemapped.sam
#awk '$2!=2048 && $2 != 2064' two_step_graph_mapper_on_reference_traversemapped.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort >two_step_graph_mapper_traversemapped.pos
#join two_step_graph_mapper_traversemapped.pos sim.gam.truth.tsv | ../vg_sim_pos_compare.py $threshold >two_step_graph_mapper_traversemapped.compare

# 2) Normal linear mapping
two_step_graph_mapper map_to_path -t $threads -r $fasta -f sim.fa -o bwa.sam
awk '$2!=2048 && $2 != 2064' bwa.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort >bwa.pos
join bwa.pos sim.gam.truth.tsv | ../vg_sim_pos_compare.py $threshold >bwa.compare

# 3) bwa mem with no tuning
bwa mem -t threads $fasta sim.fa > bwa-untuned.sam
awk '$2!=2048 && $2 != 2064' bwa-untuned.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort >bwa-untuned.pos
join bwa-untuned.pos sim.gam.truth.tsv | ../vg_sim_pos_compare.py $threshold > bwa_untuned.compare


# 4) vg
echo "vg pan single mappping"
time vg map $vg_map_opts -G sim.gam -x $pan_xg -g $pan_gcsa -t $threads --refpos-table | sort  > vg.pos
join vg-pan-se.pos sim.gam.truth.tsv | ../vg_sim_pos_compare.py $threshold >vg.compare

# Vg using mitty reads
time vg map -f reads_mitty_with_errors.fq -x $pan_xg -g $pan_gcsa -t $threads --refpos-table | sort  > vg-pan-se-mitty.pos
join vg-pan-se-mitty.pos sim.gam.truth.mitty.tsv | ../vg_sim_pos_compare.py $threshold > vg_mitty.compare

# Seven bridges
# This will only be run if seven bridges have been run seperately and sam files are present
if [ -e  seven_bridges.sam ];
then
    awk '$2!=2048 && $2 != 2064' seven_bridges.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort > seven_bridges.pos
    join seven_bridges.pos sim.gam.truth.tsv | ../vg_sim_pos_compare.py $threshold >seven_bridges.compare

    # Seven bridges using mitty
    awk '$2!=2048 && $2 != 2064' seven_bridges_mitty.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort > seven_bridges_mitty.pos
    join  seven_bridges_mitty.pos sim.gam.truth.mitty.tsv | ../vg_sim_pos_compare.py $threshold > seven_bridges_mitty.compare

else
    echo " ====================
    NOT RUNNING SEVEN BRIDGES. MUST BE RUN SEPERATELY
============ "
fi


# 5) Two-step using vg alignments
vg map -G sim.gam -x $pan_xg -g $pan_gcsa -t $threads --output-json  > vg.json
two_step_graph_mapper predict_path -t $threads -d $obg_graph_dir -a vg.json -c $chromosomes --linear-ref-bonus 1 -o predicted_path_vg
two_step_graph_mapper map_to_path -t $threads -r predicted_path_vg.fa -f sim.fa -o two_step_graph_mapper_vg.sam
two_step_graph_mapper convert_to_reference_positions -s two_step_graph_mapper_vg.sam -d $obg_graph_dir/ -l predicted_path_vg -c $chromosomes -o two_step_graph_mapper_on_reference_vg.sam
awk '$2!=2048 && $2 != 2064' two_step_graph_mapper_on_reference_vg.sam | grep -v ^@ | awk -v OFS="\t" '{$4=($4 + 0); print}' | cut -f 1,3,4,5,14 | sed s/AS:i:// | sort >two_step_graph_mapper_vg.pos
join two_step_graph_mapper_vg.pos sim.gam.truth.tsv | ../vg_sim_pos_compare.py $threshold >two_step_graph_mapper_vg.compare

#../create_roc_plots.sh

