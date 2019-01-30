# Analyse number of variants in genome vs number of reads mapped to that part of the genome by dividing the genome into bins
import numpy as np
import logging
import matplotlib.pyplot as plt
import tkinter
logging.basicConfig(level=logging.INFO)

bin_size = 200
variants_file_name = "filtered_20.vcf"
alignments_position_file_name = "vg-multimapped-with-correct-status.pos"

chromosome_size = 63025520
n_bins = chromosome_size // bin_size

variant_counts = np.zeros(n_bins)
alignment_counts = np.zeros(n_bins)
total_bins_with_n_variants = np.zeros(30)

logging.info("Processing variants")
with open(variants_file_name) as f:
    for line in f:
        if not line.startswith("@"):
            pos = int(line.split()[1])

            bin = pos // bin_size
            variant_counts[bin] += 1

for bin in range(len(variant_counts)):
    total_bins_with_n_variants[int(variant_counts[bin])] += 1

total_alignments_in_bins_with_n_variants = np.zeros(30)


logging.info("Processing alignments")
with open(alignments_position_file_name) as f:
    for line in f:
        l = line.split()
        try:
            if False and l[5] == "1":
                continue
        except IndexError:
            continue

        if l[1] == "20":
            pos = int(l[2]) + 75
            bin = pos // bin_size

            alignment_counts[bin] += 1

            n_variants_here = variant_counts[bin]
            total_alignments_in_bins_with_n_variants[int(n_variants_here)] += 1


n_alignments_for_n_variants = [[] for i in range(int(max(variant_counts)+1))]
for bin in range(len(alignment_counts)):
    n_alignments = alignment_counts[bin]
    n_variants = int(variant_counts[bin])
    n_alignments_for_n_variants[n_variants].append(n_alignments)

with open("variant_counts.txt", "w") as f:
    f.write(','.join([str(v) for v in variant_counts]))

with open("alignment_counts.txt", "w") as f:
    f.write(",".join([str(a) for a in alignment_counts]))

print(total_alignments_in_bins_with_n_variants)
#print(average_alignments_in_bins_with_n_variants)
plt.boxplot(n_alignments_for_n_variants, showfliers=False)
plt.show()
print(total_bins_with_n_variants)
import sys
sys.exit()



average_alignments_in_bins_with_n_variants = total_alignments_in_bins_with_n_variants / (total_bins_with_n_variants + 0.1)

plt.plot(range(0, len(average_alignments_in_bins_with_n_variants)), average_alignments_in_bins_with_n_variants)
plt.show()

plt.figure()
logging.info("Making scatter plot")
plt.scatter(variant_counts, alignment_counts, s=0.5)
plt.show()
logging.info("Done")




