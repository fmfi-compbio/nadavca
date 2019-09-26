import os
import sys
from nadavca.genome import Genome
from nadavca.kmer_model import KmerModel
from nadavca.read import Read
from nadavca.alignment import ApproximateAligner
from nadavca.estimator import ProbabilityEstimator, Chunk
import nadavca.defaults

import yaml


def estimate_snps(reference_filename,
                  reads,
                  reference=None,
                  config=None,
                  kmer_model=None,
                  bwa_executable=nadavca.defaults.BWA_EXECUTABLE,
                  independent=False,
                  group_name=nadavca.defaults.GROUP_NAME):
    if config is None:
        config = nadavca.defaults.CONFIG_FILE
    if isinstance(config, str):
        try:
            with open(config, 'r') as file:
                config = yaml.load(file)
        except FileNotFoundError:
            sys.stderr.write('failed to load config: {} not found\n'.format(config))
            return None

    if kmer_model is None:
        kmer_model = nadavca.defaults.KMER_MODEL_FILE
    if isinstance(kmer_model, str):
        try:
            kmer_model = KmerModel.load_from_hdf5(kmer_model)
        except FileNotFoundError:
            sys.stderr.write(
                'failed to load k-mer model: {} not found\n'.format(kmer_model))
            return None

    if reference is None:
        try:
            reference = Genome.load_from_fasta(reference_filename)[0].bases
        except FileNotFoundError:
            sys.stderr.write(
                "failed to process: reference {} doesn't exist\n".format(reference_filename))
            return None

    approximate_aligner = ApproximateAligner(bwa_executable, reference, reference_filename)
    estimator = ProbabilityEstimator(kmer_model, approximate_aligner, config)

    if isinstance(reads, str):
        read_basedir = reads
        reads = []
        for file in os.listdir(read_basedir):
            path = os.path.join(read_basedir, file)
            if not os.path.isdir(path) and path.endswith(".fast5"):
                reads.append(path)

    reads = reads[:]
    for i, read in enumerate(reads):
        if isinstance(read, str):
            reads[i] = Read.load_from_fast5(read, group_name)

    Read.normalize_reads(reads)

    if independent:
        result = []
        for read in reads:
            chunks = estimator.estimate_probabilities(reference, [read])
            result.append(chunks[0])
        return result
    else:
        return estimator.estimate_probabilities(reference, reads)

def estimate_snps_command(args):
    try:
        reference = Genome.load_from_fasta(args.reference)[0].bases
    except FileNotFoundError:
        sys.stderr.write("failed to process: reference {} doesn't exist\n".format(args.reference))
        return

    read_filenames = []
    for file in os.listdir(args.read_basedir):
        path = os.path.join(args.read_basedir, file)
        if not os.path.isdir(path) and path.endswith(".fast5"):
            read_filenames.append(path)

    chunks = estimate_snps(reference_filename=args.reference,
                           reads=read_filenames,
                           reference=reference,
                           config=args.configuration,
                           kmer_model=args.kmer_model,
                           bwa_executable=args.bwa_executable,
                           independent=args.independent,
                           group_name=args.group_name)

    if chunks is None:
        return

    if args.independent:
        if args.output:
            if not os.path.exists(args.output):
                os.makedirs(args.output)
            if not os.path.isdir(args.output):
                sys.stderr.write('Failed to create directory {}'
                                 '(maybe a file with that name exists?)\n'.format(args.output))
                return
        for chunk, read_filename in zip(chunks, read_filenames):
            output_file = None
            if args.output:
                basename = os.path.splitext(os.path.basename(read_filename))[0]
                output_file = open(os.path.join(args.output, basename + '.txt'), 'w')
            else:
                output_file = sys.stdout
            Chunk.print_head(output_file)
            chunk.print(output_file, reference)
            if args.output:
                output_file.close()

    else:
        output_file = None
        if args.output:
            output_file = open(args.output, 'w')
        else:
            output_file = sys.stdout
        Chunk.print_head(output_file)
        for chunk in chunks:
            chunk.print(output_file, reference)
        if args.output:
            output_file.close()
