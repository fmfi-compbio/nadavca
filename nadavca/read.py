import statistics
import h5py
import numpy
from scipy import interpolate
from nadavca.genome import Genome
import numpy as np


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
    def _extract_sequence_to_signal_mapping_from_moves(moves, signal_start, signal_step):
        mapping = {}
        sequence_index = 0
        base_pos = 0
        signal_pos = signal_start

        for row in moves:
            if row == 1:
                mapping[base_pos] = int(signal_pos)
            base_pos += row
            signal_pos += signal_step
        return mapping        

    @staticmethod
    def load_from_fast5(filename, basecall_group, segmentation_group="Analyses/Segmentation_000"):
        read = Read()
        with h5py.File(filename, 'r') as file:
            read_group = list(file['Raw/Reads'].values())[0]
            read.raw_signal = numpy.array(read_group['Signal'][()])
            try:
                events = file['{}/BaseCalled_template/Events'.format(basecall_group)]
                read.sequence_to_signal_mapping = Read._extract_sequence_to_signal_mapping(events)
            except:
                try:
                    moves = file['{}/BaseCalled_template/Move'.format(basecall_group)]
                    signal_start = file["{}/Summary/segmentation".format(segmentation_group)].attrs["first_sample_template"]
                    signal_step = file["{}/Summary/basecall_1d_template".format(basecall_group)].attrs["block_stride"]
                    read.sequence_to_signal_mapping = \
                        Read._extract_sequence_to_signal_mapping_from_moves(moves, signal_start, signal_step)
                except:
                    print("Cannot find Events or Move table from basecaller.")
                    raise
            read.fastq = file['{}/BaseCalled_template/Fastq'.format(basecall_group)][()].decode(
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
            read.normalized_signal = np.clip((read.raw_signal - shift) / scale, -5, 5)

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
