from pyvg.conversion import vg_json_file_to_interval_collection
from offsetbasedgraph import IntervalCollection, Interval
from collections import defaultdict
import pickle
import logging
from scipy.stats import binom_test
import numpy as np


class PathPredicter:
    def __init__(self, alignment_file_name, graph, sequence_graph, chromosome, linear_interval_path, variant_edges,
                 out_file_base_name, linear_ref_bonus=1, max_nodes_to_traverse=None, input_is_edgecounts=False):
        self.input_is_edgecounts = input_is_edgecounts
        self.alignment_file_name = alignment_file_name
        self.graph = graph
        self.sequence_graph = sequence_graph
        self.chromosome = chromosome
        self.linear_interval_path = linear_interval_path
        self.out_file_base_name = out_file_base_name
        self.linear_path_nodes = linear_interval_path.nodes_in_interval()
        self.linear_ref_bonus = linear_ref_bonus
        self.max_nodes_to_traverse = max_nodes_to_traverse
        self.binom_test_cache = {}
        self.alignments = None
        self._read_alignments()
        self.edge_counts = None
        self.variant_edges = variant_edges
        self._get_edge_counts_from_alignments()
        self.predict_path()


    def binom_test(self, x, n, p):
        if (x, n, p) in self.binom_test_cache:
            return self.binom_test_cache[(x, n, p)]
        
        pval = binom_test(x, n, p)
        self.binom_test_cache[(x, n, p)] = pval

        return pval

    def _read_alignments(self):
        if self.input_is_edgecounts:
            self.alignments = pickle.load(open(self.alignment_file_name + "_" + self.chromosome + ".edgecounts", "rb"))
        elif self.alignment_file_name.endswith(".json"):
            logging.info("Reading vg alignments")
            self.alignments = vg_json_file_to_interval_collection(self.alignment_file_name).intervals
            logging.info("Read vg alignments")
        elif self.alignment_file_name.endswith(".graphnodes"):
            self.alignments = (Interval(0, 1, [int(n) for n in line.strip().split()[1].split(",")])
                               for line in open(self.alignment_file_name))
        elif self.alignment_file_name.endswith(".graphalignments"):
            self.alignments = (Interval.from_file_line(line.strip().split("\t")[1]) for line in open(self.alignment_file_name) if line.strip().split("\t")[1] != ".")
        else:
            self.alignments = IntervalCollection.from_file(self.alignment_file_name).intervals

    def _get_edge_counts_from_alignments(self):
        if self.input_is_edgecounts:
            self.edge_counts = defaultdict(int)
            logging.info("Input is edge counts.")
            total_counts = 0
            for edge, counts in self.alignments.items():
                self.edge_counts["%d-%d" % edge] = counts
                total_counts += counts

            logging.info("Average number of reads per edge: %d" % (total_counts / (len(self.alignments.values())+1)))

            return

        logging.info("Getting edge counts")
        edge_counts = defaultdict(int)
        graph_min_node = self.graph.blocks.node_id_offset - 1
        graph_max_node = self.graph.blocks.max_node_id()

        for i, alignment in enumerate(self.alignments):
            if i % 10000000 == 0:
                logging.info("%d alignments processed" % i)

            if abs(alignment.region_paths[0]) < graph_min_node or abs(alignment.region_paths[0]) > graph_max_node:
                continue

            alignment.graph = self.graph
            if alignment.region_paths[0] < 0:
                alignment = alignment.get_reverse()


            prev_node = alignment.region_paths[0]
            for node in alignment.region_paths[1:]:
                edge_id = "%d-%d" % (prev_node, node)
                edge_counts[edge_id] += 1
                prev_node = node

        self.edge_counts = edge_counts

    def predict_path(self):
        logging.info("Using linear bonus %d on chromosome %s" % (self.linear_ref_bonus, self.chromosome))

        logging.info("Using linear out base name %s" % self.out_file_base_name)
        out_file = open("%s_%s.fasta" % (self.out_file_base_name, self.chromosome), "w")

        n_insertions = 0
        n_snps = 0
        n_deletions = 0

        # Traverse
        first_nodes = self.graph.get_first_blocks()
        assert len(first_nodes) == 1

        logging.info("N nodes in graph: %d" % len(self.graph.blocks))

        node = first_nodes[0]
        assert node in self.linear_path_nodes, "Start node should be in linear ref"

        path = []
        n_ambigious = 0
        edges_chosen = set()
        i = 0
        n_special_case = 0
        while True:
            if i % 1000000 == 0:
                logging.info("%d nodes in graph traversed on chrom %s" % (i, self.chromosome))
            i += 1

            if self.max_nodes_to_traverse is not None and i > self.max_nodes_to_traverse:
                logging.warning("Stopped traversing before end because max node to traverse was set")
                break

            path.append(node)

            next_nodes = self.graph.adj_list[node]
            if len(next_nodes) == 0:
                logging.info("Done on node %d" % node)
                break
            elif len(next_nodes) == 1:
                node = next_nodes[0]
            else:
                most_reads = 0

                # Algorithm:
                # if there is any variant edge out with enough reads, follow that
                # if not:
                #  follow edge to linear ref with most reads if we are on a variant node
                #  else: follow linear ref edge


                # Find the variant node with most reads, but only if it has enough reads (if not take next linear ref node)
                next_nodes = self.graph.adj_list[node]
                if len(next_nodes) == 1:
                    most_reads_node = next_nodes[0]
                else:
                    reads_per_node = {n: self.edge_counts["%s-%s" % (node, n)] for n in next_nodes if (node, n) in self.variant_edges}
                    most_reads_nodes = sorted(reads_per_node, key=reads_per_node.get, reverse=True)
                    if len(most_reads_nodes) > 0 and reads_per_node[most_reads_nodes[0]] > self.linear_ref_bonus:
                        most_reads_node = most_reads_nodes[0]
                    else:
                        if node not in self.linear_path_nodes:
                            # On variant node
                            counts = {n: self.edge_counts["%s-%s" % (node, n)] for n in next_nodes if n in self.linear_path_nodes}
                            most_reads_node = sorted(counts, key=counts.get, reverse=True)[0]
                        else:
                            # Not on variant node, choose next node in linear ref path
                            next_linear_ref_node = [n for n in next_nodes if (node, n) not in self.variant_edges]
                            assert len(next_linear_ref_node) == 1
                            most_reads_node = next_linear_ref_node[0]


                assert most_reads_node is not None

                # Decide what kind of variant this is
                if node in self.linear_path_nodes and most_reads_node not in self.linear_path_nodes:
                    if self.graph.adj_list[most_reads_node][0] in self.graph.adj_list[node]:
                        n_insertions += 1
                    else:
                        n_snps += 1

                elif node in self.linear_path_nodes and most_reads_node in self.linear_path_nodes:
                    other_linear_out = [n for n in self.graph.adj_list[node] if n in self.linear_path_nodes and n != most_reads_node and n < most_reads_node]
                    if len(other_linear_out) > 0:
                        n_deletions += 1

                edges_chosen.add("%d-%d" % (node, most_reads_node))
                node = most_reads_node

        # Find statistics of chosen nodes
        nodes_chosen = set(path)
        n_on_linear = len(nodes_chosen.intersection(self.linear_path_nodes))
        n_not_on_linear = len(nodes_chosen) - n_on_linear

        linear_ref_interval = Interval(0, self.graph.blocks[path[-1]].length(), path, self.graph)
        IntervalCollection([linear_ref_interval]).to_file("%s_%s.intervalcollection" % (self.out_file_base_name, self.chromosome),
                                                          text_file=True)

        logging.info("=== STATS FOR CHROMOSOME %s ===" % self.chromosome)
        logging.info("N SNPs found: %d" % n_snps)
        logging.info("N insertions found: %d" % n_insertions)
        logging.info("N deletions found: %d" % n_deletions)
        logging.info("N ambigious choices: %d" % n_ambigious)
        logging.info("Total nodes in linear ref: %d" % len(self.linear_path_nodes))
        logging.info("N nodes chosen that are not in linear ref: %d " % n_not_on_linear)
        logging.info("N nodes chosen that are in linear ref: %d " % n_on_linear)
        logging.info("N special case: %d" % n_special_case)
        logging.info("N nodes in path: %d" % len(path))
        logging.info("Linear path length: %d" % linear_ref_interval.length())

        sequence = self.sequence_graph.get_interval_sequence(linear_ref_interval)

        out_file.writelines([">%s\n" % self.chromosome])
        out_file.writelines([sequence + "\n"])
        out_file.close()

