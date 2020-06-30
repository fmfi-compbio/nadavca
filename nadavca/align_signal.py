import os
import sys
from nadavca.genome import Genome
from nadavca.kmer_model import KmerModel
from nadavca.read import Read
from nadavca.alignment import ApproximateAligner
from nadavca.estimator import ProbabilityEstimator
import nadavca.defaults
import yaml
import numpy as np
from scipy.stats import linregress

def load_model_and_estimator(reference_filename, config=nadavca.defaults.CONFIG_FILE,
                     kmer_model=None,
                     bwa_executable=nadavca.defaults.BWA_EXECUTABLE):
    if kmer_model is None:
        kmer_model=nadavca.defaults.KMER_MODEL_FILE,
    if isinstance(config, str):
        try:
            with open(config, 'r') as file:
                config = yaml.load(file)
        except FileNotFoundError:
            sys.stderr.write('failed to load config: {} not found\n'.format(config))
            return None

    if isinstance(kmer_model, str):
        kmer_model = KmerModel.load_from_hdf5(kmer_model)

    print("kmer model loaded")
    try:
        references = Genome.load_from_fasta(reference_filename)
        references_dict = {r.description[1:]: r.bases for r in references}
    except FileNotFoundError:
        sys.stderr.write(
            "failed to process: reference {} doesn't exist\n".format(reference_filename))
        return None
    print("genome loaded")

    approximate_aligner = ApproximateAligner(bwa_executable, None, reference_filename,
    references_dict)
    return kmer_model, ProbabilityEstimator(kmer_model, approximate_aligner, config)

def align_signal(reference_filename,
                 reads,
                 config=nadavca.defaults.CONFIG_FILE,
                 kmer_model=nadavca.defaults.KMER_MODEL_FILE,
                 bwa_executable=nadavca.defaults.BWA_EXECUTABLE,
                 group_name=nadavca.defaults.GROUP_NAME,
                 renorm_rounds=nadavca.defaults.RENORM_ROUNDS):
    kmer_model, estimator = load_model_and_estimator(reference_filename, config, kmer_model, bwa_executable)

    for read_fn in reads:
        read = Read.load_from_fast5(read_fn, group_name)
        Read.normalize_reads([read])
        apx_alignment, alignment = estimator.get_refined_alignment(read)
        signal_cut = read.normalized_signal[alignment[0][1]:alignment[-1][2]]

        for r in range(renorm_rounds):
            if r % 2 == 0:
                bases_str = apx_alignment.reference_part
                trans = {"A": 0, "C": 1, "G": 2, "T": 3}
                bases_num = [trans[x] for x in bases_str]
                esignal = kmer_model.get_expected_signal(bases_num, [], [])
                expected = np.array(esignal)

                signal_means = []
                expected_expand = []
                al_start = alignment[0][1]
                for i, (_, s, e) in enumerate(alignment):
                    expected_expand += [expected[i]] * (e-s)
                    signal_means.append(np.mean(signal_cut[s-al_start:e-al_start]))

                slope,intercept, _, _, _ = linregress(expected, signal_means)

                read.normalized_signal = (read.normalized_signal - intercept) / slope
                signal_cut = read.normalized_signal[alignment[0][1]:alignment[-1][2]]
            else:
                apx_alignment, alignment = estimator.get_refined_alignment(read)

                signal_cut = read.normalized_signal[alignment[0][1]:alignment[-1][2]]
        yield read, (apx_alignment, alignment)

def write_binary_output(filename, raw_cut, bases, *args):
    np.savez(filename + ".npz", raw_cut, bases, *args)  

def get_output_filename(filename, args):
  base = os.path.splitext(os.path.basename(filename))[0]
  return os.path.join(args.output, base)

def align_signal_command(args):
    read_filenames = []
    for file in os.listdir(args.read_basedir):
        path = os.path.join(args.read_basedir, file)
        if not os.path.isdir(path) and path.endswith(".fast5"):
            read_filenames.append(path)

    alignments = align_signal(reference_filename=args.reference,
                              reads=read_filenames,
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

    for (read, res), read_filename in zip(alignments, read_filenames):
        if res is None:
            continue
        apx_alignment, alignment = res
        raw_cut = read.raw_signal[alignment[0][1]:alignment[-1][1]]
        if len(raw_cut) == 0:
            print("empty raw cut", read_filename)
            continue
#            raise AlignException("empty raw cut")
        alph = "ACGT"
        bases = np.full(raw_cut.shape, 'N')
        last = -47

        for base, start, end in zip(apx_alignment.reference_part, alignment[:-1,1], alignment[1:,1]):
            if start <= last:
                raise AlignException("bad alignment")
            bases[start - alignment[0][1]] = base
            last = start

        write_binary_output(get_output_filename(read_filename, args), raw_cut, bases,
            np.array((str(apx_alignment.reference_range[0]), '-' if apx_alignment.reverse_complement else '+', apx_alignment.contig_name, read.sequence)))
        

        if False:
            output_file = None
            if args.output:
                basename = os.path.splitext(os.path.basename(read_filename))[0]
                output_file = open(os.path.join(args.output, basename + '.txt'), 'w')
            else:
                output_file = sys.stdout
            output_file.write("reference_position\tevent_start\tevent_end\n")
            for line in alignment:
                output_file.write("{}\n".format('\t'.join(line)))
            if args.output:
                output_file.close()
        print("done", read_filename)
