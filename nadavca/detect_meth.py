import os
import sys
from nadavca.genome import Genome
from nadavca.kmer_model import KmerModel
from nadavca.read import Read
from nadavca.alignment import ApproximateAligner
from nadavca.estimator import ProbabilityEstimator, Chunk
from nadavca.align_signal import load_model_and_estimator
import nadavca.defaults
import numpy as np
from scipy.stats import linregress
from scipy import stats
import csv

import yaml



SMALLEST_PVAL = 1e-50

def cdf_scoring(raw, exp):
    z_scores = np.abs(raw.mean() - exp) / 0.35287208
    r_pvals = stats.norm.cdf(-z_scores) * 2.0
    return -np.log(max(SMALLEST_PVAL, r_pvals))

def calculate_meth_scores(signal_cut, alignment, apx_alignment, pattern, kmer_model):
    bases_str = apx_alignment.reference_part
    trans = {"A": 0, "C": 1, "G": 2, "T": 3}
    bases_num = [trans[x] for x in bases_str]
    esignal = kmer_model.get_expected_signal(bases_num, [], [])
    expected = np.array(esignal)
    seq = "".join(bases_str)
    
    features = []
    last_pos = 0
    al_start = alignment[0][1]

    while True:
        pos = seq.find(pattern, last_pos)
        if pos == -1:
            break
        last_pos = pos + 1
        diffs = []
        for i in range(-5, 6):
            if pos + i < 0 or pos + i >= len(alignment):
                continue
            r = signal_cut[alignment[pos+i][1] - al_start : alignment[pos+i][2] - al_start]

            if len(r) == 0:
                break
            diffs.append(cdf_scoring(r, expected[pos+i]))
    
        if len(diffs) == 11:
            features.append((pos, seq[pos-5:pos+6], diffs))
    return features

# TODO: make this configurable
def maxs3(a):
    aa = [a1+a2+a3 for a1,a2,a3 in zip(a, a[1:], a[2:])]
    return max(aa)


def detect_meth(reference_filename,
                reads,
                pattern,
                output,
                config=nadavca.defaults.CONFIG_FILE,
                kmer_model=nadavca.defaults.KMER_MODEL_FILE,
                bwa_executable=nadavca.defaults.BWA_EXECUTABLE,
                group_name=nadavca.defaults.GROUP_NAME,
                renorm_rounds=nadavca.defaults.RENORM_ROUNDS):
    kmer_model, estimator = load_model_and_estimator(reference_filename, config, kmer_model, bwa_executable)

    result_f = open(output, "w") if output is not None else sys.stdout
    result_fc = csv.writer(result_f)
    result_fc.writerow(("Filename", "Position", "Sequence context", "Position scores", "Aggregated score"))

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
        res = calculate_meth_scores(signal_cut, alignment, apx_alignment, pattern, kmer_model)

        for pos, seq, diffs in res:
            result_fc.writerow((read_fn, pos, seq, ",".join(map(str, diffs)), maxs3(diffs)))


def detect_meth_command(args):

    reads = [os.path.join(args.read_basedir, fn) for fn in os.listdir(args.read_basedir) if fn.endswith(".fast5")]

    detect_meth(args.reference, reads, args.pattern, args.output,
                args.configuration, args.kmer_model, args.bwa_executable,
                args.group_name,
                args.renorm_rounds)

    pass
