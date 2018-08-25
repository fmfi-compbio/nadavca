import subprocess
import tempfile
import os
import simplesam
import numpy
from nadavca.genome import Genome


class ApproximateAligner:
    def __init__(self, bwa_executable):
        self.bwa_executable = bwa_executable

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

    def get_alignment(self, reference, reference_filename, read):
        read_fastq_filename = None
        with tempfile.NamedTemporaryFile(mode='w', delete=False, prefix='nadavca_tmp',
                                         suffix='.fastq') as file:
            read_fastq_filename = file.name
            file.write(read.fastq)

        # A heuristic to check if bwa index has been run already
        if not os.path.isfile(reference_filename + '.bwt'):
            subprocess.run([self.bwa_executable, 'index', reference_filename],
                           stderr=subprocess.PIPE)

        bwa_output_filename = None
        with tempfile.NamedTemporaryFile(delete=True, prefix='nadavca_tmp', suffix='.sam') as file:
            bwa_output_filename = file.name

        subprocess.run([self.bwa_executable, 'mem', reference_filename, read_fastq_filename, '-o',
                        bwa_output_filename],
                       stderr=subprocess.PIPE)
        with simplesam.Reader(open(bwa_output_filename, 'r')) as reader:
            sam = reader.next()
            if not sam.mapped:
                return None
            cigar = sam.cigar
            is_reverse_complement = sam.reverse
            mapped_position = sam.pos - 1

        os.remove(read_fastq_filename)
        os.remove(bwa_output_filename)

        oriented_read = Genome.reverse_complement(
            read.sequence) if is_reverse_complement else read.sequence

        index_in_read = 0
        index_in_reference = mapped_position
        base_mapping = []
        parsed_cigar = self._parse_cigar(cigar)

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
                base_mapping[i] = (len(read.sequence) - 1 - val[0], len(reference) - 1 - val[1])
            base_mapping.reverse()

        return numpy.array(base_mapping, dtype=numpy.int), is_reverse_complement,
