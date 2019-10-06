import logging


def get_variant_edges(chromosome, graph, sequence_graph, linear_ref_nodes, haplotype_fasta):

    # Traverse
    first_nodes = graph.get_first_blocks()
    assert len(first_nodes) == 1
    logging.info("N nodes in graph: %d" % len(graph.blocks))
    node = first_nodes[0]
    assert node in linear_ref_nodes, "Start node should be in linear ref"
    offset = -1

    return


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

    while True:
        # Search for a position that matches a node in reference
        haplotype_node = haplotype.get_node_at_offset(position)

        if haplotype_node in linear_nodes:
            linear_node = haplotype_node
            break

        position += 5
        # print("On node %d, pos: %d" % (haplotype_node, pos))
        if abs(position - start_position) > 150:
            raise Exception("Cannot find any nodes in linear.")

    assert haplotype.get_offset_at_node(linear_node) < start_position + 150
    linear_offset = linear_ref_path.get_offset_at_node(linear_node)
    if linear_node == start_node:
        add = haplotype.get_node_offset_at_offset(start_position)
        # print("Is start node, add offset: %d" % add)
        linear_offset += add

    return linear_offset
