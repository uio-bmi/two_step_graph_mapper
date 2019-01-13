import logging
import shutil
import subprocess
import argparse
import sys
from multiprocessing import Process
from offsetbasedgraph import Graph, SequenceGraph, NumpyIndexedInterval, IntervalCollection
from .path_predicter import PathPredicter
logging.basicConfig(level=logging.DEBUG)

def main():
    run_argument_parser(sys.argv[1:])


def run_predict_path_single_chromosome(alignment_file_name, chromosome, graph_dir, linear_ref_bonus, out_file_base_name):
    sequence_graph = SequenceGraph.from_file(graph_dir + chromosome + ".nobg.sequences")
    graph = Graph.from_file(graph_dir + chromosome + ".nobg")
    linear_path = NumpyIndexedInterval.from_file(graph_dir + "/%s_linear_pathv2.interval" % chromosome)
    PathPredicter(alignment_file_name, graph, sequence_graph, chromosome, linear_path, out_file_base_name,
                  linear_ref_bonus=linear_ref_bonus)


def run_bwa_index(fasta_file_name):
    command = "bwa index " + fasta_file_name
    subprocess.check_output(command.split())


def run_predict_path(args):
    chromosomes = args.chromosomes.split(",")
    processes = []
    for chromosome in chromosomes:
        logging.info("Starting process for chromosome %s " % chromosome)
        process = Process(target=run_predict_path_single_chromosome,
                          args=(args.alignments, chromosome, args.data_dir, args.linear_ref_bonus, args.out_file_name))
        process.start()
        processes.append(process)

    for process in processes:
        process.join()

    # Merge all fasta files that were produces
    out_fasta = open(args.out_file_name + ".fa", "w")
    for chromosome in chromosomes:
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

    if not args.skip_bwa_index:
        logging.info("Running bwa index")
        run_bwa_index(args.out_file_name + ".fa")
    else:
        logging.info("Not creating bwa index")

def run_argument_parser(args):

    if shutil.which("bwa") is None:
        logging.error("BWA MEM cannot be found in path. Make sure BWA is installed.")
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
    subparser_predict.add_argument("-t", "--n_threads", help="Number of threads to use", type=int, default=8, required=False)
    subparser_predict.add_argument("-l", "--linear-ref-bonus", help="Bonus score for linear reference. Used to favour the linear reference.", default=1, type=int, required=False)
    subparser_predict.add_argument("-o", "--out-file-name", help="Output file name", required=True)
    subparser_predict.add_argument("-s", "--skip-bwa-index", default=False, required=False, help="Set to True to skip creation of bwa index")
    subparser_predict.set_defaults(func=run_predict_path)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)

    if hasattr(args, 'func'):
        args.func(args)
