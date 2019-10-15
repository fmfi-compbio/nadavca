import statistics
import h5py
import numpy
from scipy import interpolate
from nadavca.genome import Genome


class Read:
    def __init__(self):
        self.raw_signal = None
        self.normalized_signal = None
        self.tweaked_normalized_signal = None
        self.strand = None
        self.fastq = None
        self.sequence_to_signal_mapping = None
        self.sequence = None

    @staticmethod
    def _extract_sequence_to_signal_mapping(events):
        mapping = {}
        sequence_index = 1
        for event in events:
            sequence_index += event['move']
            if event['move'] > 0:
                mapping[sequence_index] = event['start']
        return mapping

    @staticmethod
    def load_from_fast5(filename, basecall_group):
        read = Read()
        with h5py.File(filename, 'r') as file:
            read_group = list(file['Raw/Reads'].values())[0]
            read.raw_signal = numpy.array(read_group['Signal'].value)
            events = file['{}/BaseCalled_template/Events'.format(basecall_group)]
            read.sequence_to_signal_mapping = Read._extract_sequence_to_signal_mapping(events)
            read.fastq = file['{}/BaseCalled_template/Fastq'.format(basecall_group)].value.decode(
                'ascii')
            read.sequence = Genome.create_from_fastq_string(read.fastq)[0].bases
        return read

    @staticmethod
    def normalize_reads(reads):
        total_number_of_values = 0
        for read in reads:
            total_number_of_values += read.raw_signal.shape[0]
        values = numpy.zeros(total_number_of_values, dtype=numpy.float)
        current = 0
        for read in reads:
            length = len(read.raw_signal)
            values[current: current + length] = read.raw_signal
            current += length
        shift = statistics.median(values)
        scale = statistics.median(abs(values - shift))
        for read in reads:
            read.normalized_signal = (read.raw_signal - shift) / scale

    def tweak_signal_normalization(self, alignment, expected_means):
        data = []
        for event, expected_mean in zip(alignment, expected_means):
            mean = numpy.mean(self.normalized_signal[event[0]: event[1]])
            if abs(expected_mean - mean) <= 1:
                data.append((mean, expected_mean))

        data.sort()
        means = [d[0] for d in data]
        expected_means = [d[1] for d in data]
        spline = interpolate.splrep(means, expected_means, s=len(means))
        self.tweaked_normalized_signal = interpolate.splev(self.normalized_signal, spline)
