from .util import convert_position_on_haplotype_to_position_on_linear_ref
import sys
from pysam import AlignmentFile
from tqdm import tqdm
from rough_graph_mapper.util import read_sam_fast, number_of_lines_in_file
import logging
from offsetbasedgraph import NumpyIndexedInterval
import numpy as np


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

    chromosomes = args.chromosomes.split(",")
    coordinate_maps = {}
    logging.info("Reading coordinate maps")
    for chromosome in tqdm(chromosomes):
        coordinate_maps[chromosome] = np.load("coordinate_map_" + chromosome + ".npy")

    logging.info("Converting")
    n_unmapped = 0

    for line in sys.stdin:
        if line.startswith("@"):
            print(line.strip())
            continue

        l = line.split()
        chromosome = l[2]
        if chromosome not in coordinate_maps:
            n_unmapped += 1
            print(line.strip())
            logging.warning("Could not convert line")
            logging.warning(line)
            continue

        pos = int(l[3])
        try:
            new_pos = int(coordinate_maps[chromosome][pos])
        except IndexError:
            logging.error("Could not get coordinate for position %d on chromosome %s" % (pos, chromosome))
            print(line.strip())
            continue

        l[3] = str(new_pos)
        print('\t'.join(l).strip())

    logging.info("%d sam records missed chromosome (unmapped)" % n_unmapped)

