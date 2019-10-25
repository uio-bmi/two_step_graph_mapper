import logging
from offsetbasedgraph import NumpyIndexedInterval
import numpy as np
from multiprocessing import Process


def make_mapping(chromosome, obg_graph_dir):
    logging.info("Chromosome %s, reading paths" % chromosome)
    path = NumpyIndexedInterval.from_file("predicted_path_" + chromosome + ".intervalcollection.indexed")
    ref = NumpyIndexedInterval.from_file(obg_graph_dir + "/" + chromosome + "_linear_pathv2.interval")

    def fill_zeros_with_last(arr, initial=0):
        ind = np.nonzero(arr)[0]
        cnt = np.cumsum(np.array(arr, dtype=bool))
        return np.where(cnt, arr[ind[cnt-1]], initial)


    logging.info("Making index for chr %s" % chromosome)
    index = ref._node_to_distance[path._distance_to_node - ref.min_node]
    filled_index = fill_zeros_with_last(index)

    # Add node offsets
    node_offsets = np.arange(0, len(path._distance_to_node)) - path._node_to_distance[path._distance_to_node - path.min_node]
    index_with_offsets = filled_index + node_offsets
    file_name = "coordinate_map_" + chromosome

    np.save(file_name, index_with_offsets)
    logging.info("Wrote coordinate map to %s.npy" % file_name)


def make_mappings_parallel(chromosomes, obg_graph_dir):
    processes = []
    for chromosome in chromosomes:
        p = Process(target=make_mapping, args=(chromosome, obg_graph_dir))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    logging.info("Done making all coordinate maps")
