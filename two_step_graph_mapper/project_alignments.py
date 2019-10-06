from .util import convert_position_on_haplotype_to_position_on_linear_ref
from pysam import AlignmentFile
from tqdm import tqdm
from rough_graph_mapper.util import read_sam_fast, number_of_lines_in_file
import logging
from offsetbasedgraph import NumpyIndexedInterval


def convert_position_on_haplotype_to_position_on_linear_ref(linear_ref_path, haplotype_path, position):
    start_position = position
    haplotype = haplotype_path
    linear_nodes = linear_ref_path.nodes_in_interval()
    try:
        start_node = haplotype.get_node_at_offset(position)
    except IndexError:
        logging.error(
            "Could not find start position for position %d" % position)
        return None

    haplotype_node = haplotype.get_node_at_offset(position)
    while True:
        # Search for a position that matches a node in reference

        if haplotype_node in linear_nodes:
            linear_node = haplotype_node
            break
        haplotype_node += 1

        #position += 5
        # print("On node %d, pos: %d" % (haplotype_node, pos))
        #if abs(position - start_position) > 150:
        #    raise Exception("Cannot find any nodes in linear.")

    assert haplotype.get_offset_at_node(linear_node) < start_position + 150
    linear_offset = linear_ref_path.get_offset_at_node(linear_node)
    if linear_node == start_node:
        add = haplotype.get_node_offset_at_offset(start_position)
        # print("Is start node, add offset: %d" % add)
        linear_offset += add

    return linear_offset


def run_project_alignments(args):
    """ Project mapped sam file"""

    sam = args.sam
    chromosomes = args.chromosomes.split(",")
    graph_dir = args.data_dir

    linear_ref_paths = {}
    haplotype_paths = {}

    #out_sam = AlignmentFile(args.out_sam, "w", template=AlignmentFile(sam))
    out_sam = open(args.out_sam, "w")

    logging.info("Reading linear paths")
    for chromosome in tqdm(chromosomes):
        linear_ref_paths[chromosome] = NumpyIndexedInterval.from_file(graph_dir + chromosome + "_linear_pathv2.interval")
        haplotype_paths[chromosome] = NumpyIndexedInterval.from_file(args.linear_paths_base_name + "_" + chromosome + ".intervalcollection.indexed")

    logging.info("Converting")
    n_unmapped = 0
    for sam_record in tqdm(read_sam_fast(sam), total=number_of_lines_in_file(sam)):
        chromosome = sam_record.chromosome
        if chromosome is None:
            out_sam.write(sam_record.pysam_object)
            n_unmapped += 1
            continue
        #length = len(sam_record.sequence)
        projected_start = convert_position_on_haplotype_to_position_on_linear_ref(linear_ref_paths[chromosome],
                                                                                  haplotype_paths[chromosome],
                                                                                  sam_record.start)
        sam_record.set_start(projected_start)
        if sam_record.text_line is not None:
            out_sam.writelines(["%s\n" % ('\t'.join(sam_record.text_line))])
        else:
            out_sam.write(sam_record.pysam_object)

    logging.info("%d sam records missed chromosome (unmapped)" % n_unmapped)

