import logging
logging.basicConfig(level=logging.INFO)
import os
from tqdm import tqdm
import shutil
import subprocess
import argparse
import sys
from multiprocessing import Process
from offsetbasedgraph import Graph, SequenceGraph, NumpyIndexedInterval, IntervalCollection
from .path_predicter import PathPredicter
from rough_graph_mapper.util import run_hybrid_between_bwa_and_minimap, run_bwa_mem
from .project_alignments import run_project_alignments
from pyfaidx import Fasta
from .util import get_variant_edges
import pickle

def main():
    run_argument_parser(sys.argv[1:])


def get_variant_edges_wrapper(args):
    chromosome = args.chromosome
    linear_path = NumpyIndexedInterval.from_file(args.data_dir + "/" + chromosome + "_linear_pathv2.interval")
    linear_ref_nodes = linear_path.nodes_in_interval()
    graph = Graph.from_file(args.data_dir + "/" + chromosome + ".nobg")
    sequence_graph = SequenceGraph.from_file(args.data_dir + "/" + chromosome + ".nobg.sequences")
    haplotype_fasta = Fasta(args.fasta_file)

    edges = get_variant_edges(chromosome, graph, sequence_graph, linear_ref_nodes, haplotype_fasta)

    with open(args.out_file_name, "wb") as f:
        pickle.dump(edges, f)

    logging.info("Wrote edges to file %s" % args.out_file_name)

def run_map_to_path(args):
    n_threads = args.n_threads
    if args.minimap_arguments is None:
        minimap_arguments = "-t %d -k19 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -r50 -p.5 -f90000,180000 -n2 -m20 -s40 -g200 -2K50m --heap-sort=yes -N 7 -a" % n_threads
    else:
        logging.info("Using custom minimap arguments: %s" % args.minimap_arguments)
        minimap_arguments = args.minimap_arguments

    if args.bwa_arguments is None:
        bwa_arguments = "-D 0.05 -h 10000000 -t %d" % args.n_threads
    else:
        logging.info("Using custom bwa arguments: %s" % args.bwa_arguments)
        bwa_arguments = args.bwa_arguments


    if args.skip_mapq_adjustment:
        logging.info("Only running BWA-MEM (not minimap2)")
        run_bwa_mem(args.linear_ref, args.fasta, args.output_file_name, bwa_arguments)
    else:
        run_hybrid_between_bwa_and_minimap(args.linear_ref, args.fasta, args.output_file_name,
                                       bwa_arguments=bwa_arguments,
                                       minimap_arguments=minimap_arguments)


def run_predict_path_single_chromosome(alignment_file_name, chromosome, graph_dir,
                                       linear_ref_bonus, out_file_base_name, max_nodes_to_traverse, input_is_edgecounts):
    sequence_graph = SequenceGraph.from_file(graph_dir + chromosome + ".nobg.sequences")
    graph = Graph.from_file(graph_dir + chromosome + ".nobg")
    linear_path = NumpyIndexedInterval.from_file(graph_dir + "/%s_linear_pathv2.interval" % chromosome)
    PathPredicter(alignment_file_name, graph, sequence_graph, chromosome, linear_path, out_file_base_name,
                  linear_ref_bonus=linear_ref_bonus, max_nodes_to_traverse=max_nodes_to_traverse,
                  input_is_edgecounts=input_is_edgecounts)


def run_bwa_index(fasta_file_name):
    command = "bwa-mem2 index " + fasta_file_name
    subprocess.check_output(command.split())


def run_predict_path(args):
    chromosomes = args.chromosomes.split(",")
    processes = []
    if not os.path.isfile(args.alignments) and not args.input_is_edgecounts:
        logging.error("Input alignments file %s does not exist" % args.alignments)
        sys.exit()

    for chromosome in chromosomes:
        logging.info("Starting process for chromosome %s " % chromosome)
        process = Process(target=run_predict_path_single_chromosome,
                          args=(args.alignments, chromosome, args.data_dir, args.linear_ref_bonus, args.out_file_name, args.max_nodes_to_traverse, args.input_is_edgecounts))
        process.start()
        processes.append(process)

    for process in processes:
        process.join()

    # Merge all fasta files that were produces
    out_fasta = open(args.out_file_name + ".fa", "w")
    logging.info("Merging fasta files")
    for chromosome in tqdm(chromosomes):
        with open(args.out_file_name + "_" + chromosome + ".fasta") as f:
            out_fasta.write(f.read())

    logging.info("Wrote resulting linear reference to %s" % (args.out_file_name + ".fa"))

    # Create indexed intervals for each interval file that was produced
    logging.info("Creating indexed interval for all chromosomes")
    for chromosome in chromosomes:
        file_name = args.out_file_name + "_" + chromosome + ".intervalcollection"
        graph = Graph.from_file(args.data_dir + chromosome + ".nobg")
        intervals = IntervalCollection.from_file(file_name, text_file=True, graph=graph)
        intervals = list(intervals.intervals)
        assert len(intervals) == 1, "Only a single interval in file is supported"
        interval = intervals[0]
        indexed = interval.to_numpy_indexed_interval()
        indexed.to_file(file_name + ".indexed")
        logging.info("Wrote indexed interval to file %s" % file_name + ".indexed")

    from .coordinate_mappings import make_mappings_parallel
    make_mappings_parallel(chromosomes, args.data_dir)

    if not args.skip_bwa_index:
        logging.info("Running bwa index")
        run_bwa_index(args.out_file_name + ".fa")
    else:
        logging.info("Not creating bwa index")

def run_argument_parser(args):

    if shutil.which("bwa-mem2") is None:
        logging.error("BWA MEM (bwa-mem2) cannot be found in path. Make sure BWA MEM 2 is installed.")
        sys.exit()
    if shutil.which("minimap2") is None:
        logging.error("minimap2 cannot be found in path. Make sure minimap2 is installed")
        sys.exit()
    if shutil.which("rough_graph_mapper") is None:
        logging.error("Rough Graph Mapper is not installed. Install first.")

    parser = argparse.ArgumentParser(
        description='Two Step Graph Mapper. Maps to a graph by first predicting a path, and then mapping to that path.',
        prog='two_step_graph_mapper',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers(help='Subcommands')
    subparser_predict = subparsers.add_parser("predict_path", help="Predicting path from graphalignments.")

    subparser_predict.add_argument("-d", "--data-dir", help="", required=True)
    subparser_predict.add_argument("-a", "--alignments", help="A vg json file containing alignments (converted with vg view -aj file.gam > file.json OR a offsetbasedgraph json interval file.", required=True)
    subparser_predict.add_argument("-c", "--chromosomes", help="Comma-separated list of chromosomes", required=True)
    subparser_predict.add_argument("-t", "--n-threads", help="Number of threads to use", type=int, default=8, required=False)
    subparser_predict.add_argument("-l", "--linear-ref-bonus", help="Bonus score for linear reference. Used to favour the linear reference.", default=1, type=int, required=False)
    subparser_predict.add_argument("-o", "--out-file-name", help="Output file name", required=True)
    subparser_predict.add_argument("-s", "--skip-bwa-index", default=False, required=False, help="Set to True to skip creation of bwa index")
    subparser_predict.add_argument("-m", "--max-nodes-to-traverse", type=int, default=None, required=False, help="For debugging/testing. Max number of nodes in graph to traverse.")
    subparser_predict.add_argument("-e", "--input-is-edgecounts", type=bool, default=False, required=False, help="Set to true if input is edgecounts (pickled dict of edgecounts).")
    subparser_predict.set_defaults(func=run_predict_path)

    subparser_map = subparsers.add_parser("map_to_path", help="Map to a predicted path")
    subparser_map.add_argument("-r", "--linear-ref", help="Linear reference fasta file name. Typically the file outputed by predict_path.", required=True)
    subparser_map.add_argument("-f", "--fasta", help="Fasta file to map", required=True)
    subparser_map.add_argument("-o", "--output_file_name", help="Name of output sam file", required=True)
    subparser_map.add_argument("-t", "--n-threads", help="Number of threads to use", type=int, default=8, required=False)
    subparser_map.add_argument("-m", "--minimap-arguments", help="Arguments that will be used with minimap. Only change if you know what you do.", type=str, default=None, required=False)
    subparser_map.add_argument("-b", "--bwa-arguments", help="Arguments that will be used with BWA-MEM. Only change if you know what you do.", type=str, default=None, required=False)
    subparser_map.add_argument("-q", "--skip-mapq-adjustment", help="If set to True, minimap2 will not be run to improve mapq's.", type=bool, default=False, required=False)
    subparser_map.set_defaults(func=run_map_to_path)

    subparser_project = subparsers.add_parser("convert_to_reference_positions", help="Convert positions in a mapped sam file to positions on reference genome")
    subparser_project.add_argument("-s", "--sam", help="Input sam file", required=True)
    subparser_project.add_argument("-d", "--data-dir", help="Data direcoty containing the linear path intervals for the graphs", required=True)
    subparser_project.add_argument("-l", "--linear-paths-base-name", help="Base name for linear paths. Same as --out-file-name used for predict_path.", required=True)
    subparser_project.add_argument("-c", "--chromosomes", required=True)
    subparser_project.add_argument("-o", "--out-sam", required=True, help="Name of sam file to write modified alignments to")
    subparser_project.set_defaults(func=run_project_alignments)

    s = subparsers.add_parser("get_variant_edges", help="Utility functionality: Finds edges supported by a haploid fasta sequence through the graph")
    s.add_argument("-f", "--fasta_file", help="Haplotype fasta file", required=True)
    s.add_argument("-c", "--chromosome", required=True)
    s.add_argument("-d", "--data-dir", required=True)
    s.add_argument("-o", "--out_file_name", required=True)
    s.set_defaults(func=get_variant_edges_wrapper)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)

    if hasattr(args, 'func'):
        args.func(args)

