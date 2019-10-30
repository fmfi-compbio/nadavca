import numpy
from nadavca.alphabet import alphabet, numerical_complement
from nadavca.genome import Genome
import nadavca.dtw


class Chunk:
    def __init__(self, start, end, values, coverage=None):
        self.start = start
        self.end = end
        self.values = values
        self.coverage = coverage
        if coverage is None:
            self.coverage = numpy.ones(end - start, dtype=numpy.int)

    def __lt__(self, other):
        if self.start == other.start:
            return self.end < other.end
        return self.start < other.start

    @staticmethod
    def print_head(file):
        file.write('index\tbase\tcoverage\t{}\n'.format('\t'.join(alphabet)))

    def print(self, file, reference):
        for i in range(self.start, self.end):
            values_string = '\t'.join(map('{:18.16f}'.format, self.values[i - self.start]))
            file.write('{}\t{}\t{}\t{}\n'.format(i,
                                                 reference[i],
                                                 self.coverage[i - self.start],
                                                 values_string))

class ProbabilityEstimator:
    def __init__(self, kmer_model, aligner, config):
        self.kmer_model = kmer_model
        self.aligner = aligner
        self.bandwidth = config['bandwidth']
        self.snp_prior = config['snp_prior_probability']
        self.min_event_length = config['min_alignment_event_length']
        self.distribution_samples_count = config['distribution_samples_count']
        self.noise_sigma = config['noise_sigma']
        self.short_event_acceptability = config['short_event_acceptability']
        self.model_wobbling = config['model_wobbling']
        self.model_transitions = config['model_transitions']
        self.normalization_event_length = config['normalization_event_length']
        self.tweak_signal_normalization = config['tweak_signal_normalization']

    def _normalize_log_likelihoods(self, likelihoods, reference):
        shift = likelihoods[0][reference[0]]
        return (likelihoods - shift) / self.normalization_event_length

    def _get_read_context(self, read, read_sequence_range):
        start, end = read_sequence_range
        k = self.kmer_model.get_k()
        central_position = self.kmer_model.get_central_position()
        context_before = read.sequence[start - central_position : start]
        context_before = Genome.to_numerical(context_before)
        context_after = read.sequence[end : end + k - central_position - 1]
        context_after = Genome.to_numerical(context_after)
        return context_before, context_after

    def _estimate_log_likelihoods(self, reference, read):
        approximate_alignment = self.aligner.get_signal_alignment(read, self.bandwidth)
        if approximate_alignment is None:
            return None

        start_in_reference, end_in_reference = approximate_alignment.reference_range
        reference_part = reference[start_in_reference: end_in_reference]
        if approximate_alignment.reverse_complement:
            reference_part = Genome.reverse_complement(reference_part)
        reference_part = Genome.to_numerical(reference_part)

        start_in_signal, end_in_signal = approximate_alignment.signal_range
        signal = read.normalized_signal[start_in_signal : end_in_signal]
        context_before, context_after = \
            self._get_read_context(read,
                                   approximate_alignment.read_sequence_range)

        if self.tweak_signal_normalization:
            refined_alignment = nadavca.dtw.refine_alignment(
                signal=signal,
                reference=reference_part,
                context_before=context_before,
                context_after=context_after,
                approximate_alignment=approximate_alignment.alignment,
                bandwidth=self.bandwidth,
                short_event_acceptability=self.short_event_acceptability,
                min_event_length=self.min_event_length,
                distribution_samples_count=self.distribution_samples_count,
                noise_sigma=self.noise_sigma,
                kmer_model=self.kmer_model,
                model_transitions=False
            )

            refined_alignment = numpy.array(refined_alignment) + start_in_signal

            expected_signal = self.kmer_model.get_expected_signal(
                reference_part,
                context_before,
                context_after
            )
            read.tweak_signal_normalization(refined_alignment, expected_signal)
            signal = read.tweaked_normalized_signal[start_in_signal : end_in_signal]

        log_likelihoods = nadavca.dtw.estimate_log_likelihoods(
            signal=signal,
            reference=reference_part,
            context_before=context_before,
            context_after=context_after,
            approximate_alignment=approximate_alignment.alignment,
            bandwidth=self.bandwidth,
            short_event_acceptability=self.short_event_acceptability,
            kmer_model=self.kmer_model,
            model_wobbling=self.model_wobbling
        )

        log_likelihoods = numpy.array(log_likelihoods)
        log_likelihoods = self._normalize_log_likelihoods(log_likelihoods, reference_part)

        if approximate_alignment.reverse_complement:
            complement_likelihoods = numpy.zeros(log_likelihoods.shape, dtype=numpy.float)
            for i, line in enumerate(log_likelihoods):
                for j in range(len(alphabet)):
                    complement_likelihoods[i][j] = line[numerical_complement[j]]
            log_likelihoods = numpy.flipud(complement_likelihoods)

        return Chunk(start_in_reference, end_in_reference, log_likelihoods)

    def _corrected_priors(self, context_positions):
        c = len(alphabet) - 1
        p_1 = 1 - self.snp_prior
        p_2 = self.snp_prior / c
        snp_hypothesis_prior = 1 / (p_1 / p_2 + (1 - context_positions) * c)
        nonsnp_hypothesis_prior = 1 - snp_hypothesis_prior * c
        return snp_hypothesis_prior, nonsnp_hypothesis_prior

    def _compute_posterior(self, log_likelihoods, reference):
        probabilities = numpy.zeros(log_likelihoods.shape, dtype=numpy.float)
        k = self.kmer_model.get_k()
        for i, _ in enumerate(probabilities):
            context_start = max(0, i - k + 1)
            context_end = min(i + k, len(probabilities))
            max_likelihood_in_context = numpy.max(log_likelihoods[context_start: context_end])
            snp_hypothesis_prior, nonsnp_hypothesis_prior = self._corrected_priors(
                context_end - context_start - 1)
            for j, base in enumerate(alphabet):
                prior_probability = snp_hypothesis_prior if base != reference[
                    i] else nonsnp_hypothesis_prior
                probabilities[i][j] = numpy.exp(
                    log_likelihoods[i][j] - max_likelihood_in_context) * prior_probability
                if base == reference[i]:
                    for i_2 in range(context_start, context_end):
                        if i_2 == i:
                            continue
                        for j_2, base2 in enumerate(alphabet):
                            if base2 == reference[i_2]:
                                continue
                            probabilities[i][j] += numpy.exp(
                                log_likelihoods[i_2][
                                    j_2] - max_likelihood_in_context) * snp_hypothesis_prior
            probabilities[i] /= sum(probabilities[i])
        return probabilities

    def get_refined_alignment(self, read):
        approximate_alignment = self.aligner.get_signal_alignment(read, self.bandwidth)
        if approximate_alignment is None:
            return None

        start_in_reference, end_in_reference = approximate_alignment.reference_range
        reference_part = Genome.to_numerical(approximate_alignment.reference_part)

        start_in_signal, end_in_signal = approximate_alignment.signal_range
        signal = read.normalized_signal[start_in_signal : end_in_signal]
        context_before, context_after = \
            self._get_read_context(read,
                                   approximate_alignment.read_sequence_range)

        refined_alignment = nadavca.dtw.refine_alignment(
            signal=signal,
            reference=reference_part,
            context_before=context_before,
            context_after=context_after,
            approximate_alignment=approximate_alignment.alignment,
            bandwidth=self.bandwidth,
            short_event_acceptability=self.short_event_acceptability,
            min_event_length=self.min_event_length,
            distribution_samples_count=self.distribution_samples_count,
            noise_sigma = self.noise_sigma,
            kmer_model=self.kmer_model,
            model_transitions=self.model_transitions
        )

        result = numpy.zeros((len(refined_alignment), 3), dtype=int)
        for reference_position, event_range in enumerate(refined_alignment):
            event_start, event_end = event_range
            result[reference_position][1] = event_start + start_in_signal
            result[reference_position][2] = event_end + start_in_signal
            if approximate_alignment.reverse_complement:
                result[reference_position][0] = end_in_reference - reference_position - 1
            else:
                result[reference_position][0] = start_in_reference + reference_position
        return approximate_alignment, result


    def estimate_probabilities(self, reference, reads):
        chunks = []
        for read in reads:
            chunk = self._estimate_log_likelihoods(reference, read)
            if chunk is not None:
                chunks.append(chunk)
        chunks.sort()
        chunk_groups = []
        current_group = []
        current_start = None
        current_end = None
        for i, chunk in enumerate(chunks):
            if current_start is None:
                current_start = chunk.start
                current_end = chunk.end
            current_group.append(chunk)
            current_end = max(current_end, chunk.end)
            if i + 1 >= len(chunks) or chunks[i + 1].start >= current_end:
                chunk_groups.append(((current_start, current_end), current_group))
                current_group = []
                current_start = None
                current_end = None

        result = []
        for group in chunk_groups:
            start, end = group[0]
            group_members = group[1]
            coverage = numpy.zeros(end - start, dtype=numpy.int)
            log_likelihoods = numpy.zeros((end - start, len(alphabet)), dtype=numpy.float)
            for chunk in group_members:
                for i in range(chunk.start, chunk.end):
                    coverage[i - start] += 1
                log_likelihoods[chunk.start - start: chunk.end - start] += chunk.values

            probabilities = self._compute_posterior(log_likelihoods, reference[start:end])

            result.append(Chunk(start, end, probabilities, coverage))
        return result
