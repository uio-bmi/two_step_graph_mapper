# Make reference graph
vg construct -C -R 6 -r linear_ref.fa -v 1000genomes_variants.vcf.gz -t 20 -m 32 > 6.vg
vg index -x wg.xg 6.vg
vg prune -r 6.vg > 6.pruned.vg
vg index -g wg.gcsa 6.pruned.vg

vg view -Vj 6.vg > 6.json
graph_peak_caller create_ob_graph 6.json
graph_peak_caller find_linear_path -g 6.nobg 6.json 6 6_linear_pathv2.interval


# Make giab graph (our individual)
vg construct -C -R 6 -r linear_ref.fa -v giab_variants.vcf.gz -t 20 -m 1000 --alt-paths > giab_chr6.vg

# Index the graph (make a gbwt index that will include all phased variants as paths. We will use this gbwt index later to pull out paths containing phased variants)
vg index -x giab.xg -G giab.gbwt -v giab_variants.vcf.gz giab_chr6.vg

# Get all haplotype paths in gam for each chromosome and each haplotype
# Also, get haplotype paths in vg format
# (need the vg paths for later removing all nodes not part of paths to make graphs for single haplotypes, which will be used for read simulation)
vg paths --gbwt giab.gbwt --extract-vg -x giab.xg -q _thread_HG002_6_0 > haplotype_6_0.paths &
vg paths --gbwt giab.gbwt --extract-vg -x giab.xg -q _thread_HG002_6_1 > haplotype_6_1.paths &

vg paths --gbwt giab.gbwt --extract-gam -x giab.xg -q _thread_HG002_6_0 > haplotype_6_0.gam &
vg paths --gbwt giab.gbwt --extract-gam -x giab.xg -q _thread_HG002_6_1 > haplotype_6_1.gam &
wait


# Extract a linear reference path through each chromosome
vg paths -X -x giab.xg > giab_reference_paths.gam
vg view -aj giab_reference_paths.gam > giab_reference_paths.json
python3 ../parse_paths.py giab_reference_paths.json giab_reference_path.intervalcollection



# Convert haplotype paths and vg graphs to json
vg view -aj haplotype_6_0.gam > haplotype_6_0.json &
vg view -aj haplotype_6_1.gam > haplotype_6_1.json &
vg view -Vj giab_chr6.vg > giab_chr6.json &

# Create ob graph for each chromosome
graph_peak_caller create_ob_graph giab_chr6.json 

# Note: The haplotype paths are not complete. They may have gaps in them
# Thus, we want to traverse the giab graph, folloing the haplotype paths, but fill in with reference paths where possible
python3 ../make_linear_reference_and_interval.py giab_chr6.nobg giab_reference_path_6.intervalcollection haplotype_6_0.json haplotype_6_1.json haplotype_6_ 6 


cat haplotype_*__0.fasta >> haplotype0.fasta
cat haplotype_*__1.fasta >> haplotype1.fasta


# Index reference and haplotype intervals (we need them later to convert coordinates from haplotype interval offset to linear ref offset)
graph_peak_caller index_interval -g giab_chr6.nobg giab_reference_path_6.intervalcollection &
graph_peak_caller index_interval -g giab_chr6.nobg haplotype_6__0.intervalcollection &
graph_peak_caller index_interval -g giab_chr6.nobg haplotype_6__1.intervalcollection &

# Make modified graphs only containing the haplotypes, will be used to simulate reads from
vg mod -D giab_chr6.vg > giab_chr6_haplotype0.tmp && \
cat haplotype_6_0.paths >> giab_chr6_haplotype0.tmp && \
vg mod -N giab_chr6_haplotype0.tmp | vg mod -D - > giab_chr6_haplotype0.vg && \
vg index -x giab_chr6_haplotype0.xg giab_chr6_haplotype0.vg &

vg mod -D giab_chr6.vg > giab_chr6_haplotype1.tmp && \
cat haplotype_6_1.paths >> giab_chr6_haplotype1.tmp && \
vg mod -N giab_chr6_haplotype1.tmp | vg mod -D - > giab_chr6_haplotype1.vg && \
vg index -x giab_chr6_haplotype1.xg giab_chr6_haplotype1.vg &


# Create a whole genome graph with only reference, used by vg to annotate positions on
cat giab_chr?.vg giab_chr??.vg > giab_all_chromosomes.vg
# Extract reference paths in vg format for each chromosome (bit hackish since vg only allows fetching by prefix, not by actual name)
> reference_paths.vg
chromosomes="6"
for chromosome in $(echo $chromosomes | tr "," "\n")
	do
	vg paths -Q $chromosome -V -x giab.xg >> reference_paths.vg
done

# Verify by vg paths -L -v reference_paths.vg
# Remove all paths from giab_all_chromosomes.vg and add these paths:
vg mod -D giab_all_chromosomes.vg > tmp.vg
cat reference_paths.vg >> tmp.vg
vg mod -N tmp.vg > giab_reference.vg
vg index -x giab_reference.xg giab_reference.vg





