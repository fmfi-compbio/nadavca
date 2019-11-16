from collections import namedtuple
import subprocess
import tempfile
import os
import sys
import numpy
import simplesam
from nadavca.genome import Genome

ApproximateSignalAlignment = namedtuple('ApproximateSignalAlignment',
                                        ['alignment',
                                         'signal_range',
                                         'reference_range',
                                         'read_sequence_range',
                                         'reverse_complement',
                                         'reference_part',
                                         'contig_name'])

class ApproximateAligner:
    # Reference here is only for backwards compatibility with single contig mode
    def __init__(self, bwa_executable, reference, reference_filename, references_dict=None):
        self.bwa_executable = bwa_executable
        assert reference is not None or references_dict is not None, "one of reference or reference_dict must be not None"
        self.reference = reference
        self.references_dict = references_dict
        self.reference_filename = reference_filename

        # A heuristic to check if bwa index has been run already
        if not os.path.isfile(self.reference_filename + '.bwt'):
            subprocess.run([self.bwa_executable, 'index', self.reference_filename],
                           stderr=subprocess.PIPE, check=True)

        self.bwapy_aligner = None
        try:
            from bwapy import BwaAligner
            options = '-x ont2d'
            self.bwapy_aligner = BwaAligner(reference_filename, options=options)
        except ImportError:
            sys.stderr.write("Could't import bwapy, will use bwa executable to align reads\n")

    @staticmethod
    def _parse_cigar(cigar):
        num = 0
        result = []
        for character in cigar:
            if character >= '0' and character <= '9':
                num *= 10
                num += ord(character) - ord('0')
            else:
                result.append((num, character))
                num = 0
        return result

    @staticmethod
    def convert_mapping(base_mapping, read):
        result = []
        for index_in_read, index_in_reference in base_mapping:
            if index_in_read in read.sequence_to_signal_mapping:
                result.append((read.sequence_to_signal_mapping[index_in_read], index_in_reference))
        return numpy.array(result, dtype=numpy.int)

    def _get_reference_contig(self, contig_name):
        if self.references_dict is not None:
            return self.references_dict[contig_name]

        # TODO: remove this afterwards
        return self.reference

    def _get_base_alignment(self, read):
        if self.bwapy_aligner:
            alignments = self.bwapy_aligner.align_seq(''.join(read.sequence))
            if len(alignments) == 0:
                return None
            alignment = alignments[0]
            cigar = alignment.cigar
            is_reverse_complement = alignment.orient == '-'
            mapped_position = alignment.pos
            contig_name = alignment.rname

        else:
            read_fastq_filename = None
            with tempfile.NamedTemporaryFile(mode='w', delete=False, prefix='nadavca_tmp',
                                             suffix='.fastq') as file:
                read_fastq_filename = file.name
                file.write(read.fastq)

            bwa_output_filename = None
            with tempfile.NamedTemporaryFile(delete=True,
                                             prefix='nadavca_tmp',
                                             suffix='.sam') as file:
                bwa_output_filename = file.name

            subprocess.run([self.bwa_executable, 'mem', self.reference_filename,
                            read_fastq_filename, '-o', bwa_output_filename],
                           stderr=subprocess.PIPE,
                           check=True)
            with simplesam.Reader(open(bwa_output_filename, 'r')) as reader:
                sam = reader.next()
                if not sam.mapped:
                    return None
                cigar = sam.cigar
                is_reverse_complement = sam.reverse
                mapped_position = sam.pos - 1
                contig_name = sam.rname

            os.remove(read_fastq_filename)
            os.remove(bwa_output_filename)

        oriented_read = Genome.reverse_complement(
            read.sequence) if is_reverse_complement else read.sequence

        index_in_read = 0
        index_in_reference = mapped_position
        base_mapping = []
        parsed_cigar = self._parse_cigar(cigar)
        reference = self._get_reference_contig(contig_name)

        for num, operation in parsed_cigar:
            if operation == 'S':
                index_in_read += num
            elif operation == 'M':
                for i in range(num):
                    if reference[index_in_reference] == oriented_read[index_in_read]:
                        base_mapping.append((index_in_read, index_in_reference))
                    index_in_read += 1
                    index_in_reference += 1
            elif operation == 'D':
                index_in_reference += num
            elif operation == 'I':
                index_in_read += num
            else:
                raise ValueError('Unknown cigar operation: {}'.format(operation))

        if is_reverse_complement:
            for i, val in enumerate(base_mapping):
                base_mapping[i] = (len(read.sequence) - 1 - val[0],
                                   len(reference) - 1 - val[1])
            base_mapping.reverse()

        return numpy.array(base_mapping, dtype=numpy.int), is_reverse_complement, contig_name

    def get_signal_alignment(self, read, bandwidth):
        base_alignment = self._get_base_alignment(read)
        if base_alignment is None:
            return None

        base_mapping, is_reverse_complement, contig_name = base_alignment
        signal_mapping = self.convert_mapping(base_mapping, read)

        start_in_reference = signal_mapping[0][1]
        end_in_reference = signal_mapping[-1][1] + 1
        signal_mapping[:, 1] -= start_in_reference
        reference = self._get_reference_contig(contig_name)

        if is_reverse_complement:
            start_in_reference, end_in_reference = \
                len(reference) - end_in_reference, \
                len(reference) - start_in_reference
        reference_range = (start_in_reference, end_in_reference)

        start_in_signal = signal_mapping[0][0]
        end_in_signal = signal_mapping[-1][0] + 1
        extended_start_in_signal = max(0, start_in_signal - bandwidth)
        extended_end_in_signal = min(len(read.normalized_signal), end_in_signal + bandwidth)
        signal_range = (extended_start_in_signal, extended_end_in_signal)

        signal_mapping[:, 0] -= extended_start_in_signal

        start_in_read_sequence = base_mapping[0][0]
        end_in_read_sequence = base_mapping[-1][0] + 1
        read_sequence_range = (start_in_read_sequence, end_in_read_sequence)
        reference_part = reference[start_in_reference: end_in_reference]
        if is_reverse_complement:
            reference_part = Genome.reverse_complement(reference_part)

        return ApproximateSignalAlignment(alignment=signal_mapping,
                                          signal_range=signal_range,
                                          reference_range=reference_range,
                                          read_sequence_range=read_sequence_range,
                                          reverse_complement=is_reverse_complement,
                                          reference_part=reference_part,
                                          contig_name=contig_name
                                          )
