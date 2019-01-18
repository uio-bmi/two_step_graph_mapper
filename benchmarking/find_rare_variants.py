import offsetbasedgraph.vcfmap as vcfmap
import offsetbasedgraph as obg

def get_alle_freq(parts):
    AF = parts[7].split(";")[1].split("=")
    assert AF[0] == "AF", (parts, AF)
    try:
        return max(float(af) for af in AF[1].split(","))
    except:
        print(parts, AF)
        raise


def get_common_variants(vcf_filename):
    lines = (line for line in open(vcf_filename) if not line.startswith("#"))
    parts = (line.split("\t") for line in lines)
    ids = (p for part in parts if get_alle_freq(part)>0.01 for p in part[2].split(";"))
    return set(ids)


def tmp(chromosome, folder, vcf_name):
    base_name = folder + "/" + str(chromosome)
    graph = obg.Graph.from_file(base_name+".nobg")
    seq_graph = obg.SequenceGraph.from_file(base_name + ".nobg.sequences")
    reference = obg.NumpyIndexedInterval.from_file(
        base_name + "_linear_pathv2.interval")


def get_rare_node_ids(all_vcf, individual_vcf, entry_to_edge):
    common_variants = get_common_variants(all_vcf)
    print(len(common_variants))
    vcf_entries = vcfmap.get_vcf_entries(individual_vcf, common_variants)
    variants = (entry_to_edge(entry) for entry in vcf_entries)
    nodes = (node for var in variants if type(var) in (vcfmap.INS, vcfmap.SNP) for node in var.nodes)
    return set(nodes)


def contains_rare_variant(interval, rare_nodes):
    return any(region_path in rare_nodes for region_path in interval.region_paths)


if __name__ == "__main__":
    base_name = "giab_chr6"
    graph = obg.Graph.from_file(base_name + ".nobg")
    seq_graph = obg.SequenceGraph.from_file(base_name + ".nobg.sequences")
    reference = obg.NumpyIndexedInterval.from_file("giab_reference_path_6.intervalcollection.indexed")
    entry_to_edge = vcfmap.entry_to_edge_func(graph, reference, seq_graph)
    nodes = get_rare_node_ids("1000genomes_variants.vcf", "giab_variants.vcf", entry_to_edge)
    print(len(nodes))
