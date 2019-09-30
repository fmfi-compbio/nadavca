import os
import sys
from nadavca.genome import Genome
from nadavca.kmer_model import KmerModel
from nadavca.read import Read
from nadavca.alignment import ApproximateAligner
from nadavca.estimator import ProbabilityEstimator
import nadavca.defaults
import yaml


def align_signal(reference_filename,
                 reads,
                 reference=None,
                 config=nadavca.defaults.CONFIG_FILE,
                 kmer_model=nadavca.defaults.KMER_MODEL_FILE,
                 bwa_executable=nadavca.defaults.BWA_EXECUTABLE,
                 group_name=nadavca.defaults.GROUP_NAME):
    if isinstance(config, str):
        try:
            with open(config, 'r') as file:
                config = yaml.load(file)
        except FileNotFoundError:
            sys.stderr.write('failed to load config: {} not found\n'.format(config))
            return None

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

    return [estimator.get_refined_alignment(reference, read) for read in reads]

def align_signal_command(args):
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

    alignments = align_signal(reference_filename=args.reference,
                              reads=read_filenames,
                              reference=reference,
                              config=args.configuration,
                              kmer_model=args.kmer_model,
                              bwa_executable=args.bwa_executable,
                              group_name=args.group_name)

    if args.output:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        if not os.path.isdir(args.output):
            sys.stderr.write('Failed to create directory {}'
                             '(maybe a file with that name exists?)\n'.format(args.output))
            return

    for alignment, read_filename in zip(alignments, read_filenames):
        output_file = None
        if args.output:
            basename = os.path.splitext(os.path.basename(read_filename))[0]
            output_file = open(os.path.join(args.output, basename + '.txt'), 'w')
        else:
            output_file = sys.stdout
        output_file.write("signal_position\treference_position\n")
        for line in alignment:
            output_file.write("{}\t{}\n".format(line[0], line[1]))
        if args.output:
            output_file.close()
