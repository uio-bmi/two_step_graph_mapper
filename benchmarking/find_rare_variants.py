import offsetbasedgraph.vcfmap as vcfmap
import offsetbasedgraph as obg
import logging
logging.basicConfig(level=logging.INFO)

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
    return any(abs(region_path) in rare_nodes for region_path in interval.region_paths)


if __name__ == "__main__":
    from pyvg.conversion import parse_vg_json_alignments
    import sys
    # python3 ../find_rare_variants.py giab_chr6.nobg sim.json giab_reference_path_6.intervalcollection.indexed 1000genomes_variants.vcf giab_variants.vcf

    simulation_graph_file_name = sys.argv[1]
    vg_json_alignments_filename = sys.argv[2]
    simulated_graph_reference_path_file_name = sys.argv[3]
    full_vcf_file_name = sys.argv[4]
    simulation_vcf_file_name = sys.argv[5]


    #base_name = "giab_chr6"
    logging.info("Reading data")
    graph = obg.Graph.from_file(simulation_graph_file_name)
    alignments = parse_vg_json_alignments(vg_json_alignments_filename, graph)
    seq_graph = obg.SequenceGraph.from_file(simulation_graph_file_name + ".sequences")
    reference = obg.NumpyIndexedInterval.from_file(simulated_graph_reference_path_file_name)
    logging.info("Reading vcf")
    entry_to_edge = vcfmap.entry_to_edge_func(graph, reference, seq_graph)
    rare_nodes = get_rare_node_ids(full_vcf_file_name, simulation_vcf_file_name, entry_to_edge)

    i = 0
    for name, interval in alignments:
        if i % 100000 == 0:
            logging.info("%d alignments processed" % i)
        i += 1

        if contains_rare_variant(interval, rare_nodes):
            print(name + "\t" + "1")
        else:
            print(name + "\t" + "0")

    logging.info("%d rare nodes" % len(rare_nodes))
