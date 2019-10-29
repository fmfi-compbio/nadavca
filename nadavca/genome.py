from nadavca.alphabet import inv_alphabet, complement
import numpy

class Genome:
    def __init__(self, desc_line=None):
        self.bases = []
        self.description = desc_line

    def _append_line(self, line):
        for base in line:
            self.bases.append(base)

    def _convert_to_numpy(self):
        self.bases = numpy.array(self.bases)

    @staticmethod
    def to_numerical(sequence):
        return numpy.array([inv_alphabet[base] for base in sequence], dtype=numpy.int)

    @staticmethod
    def reverse_complement(sequence):
        res = [complement[x] for x in sequence]
        res.reverse()
        return numpy.array(res)

    @staticmethod
    def load_from_fasta(filename):
        result = []
        with open(filename, 'r') as file:
            current = None
            for line in file:
                if line[0] == '>':
                    current = Genome(line.rstrip().split(' ')[0])
                    result.append(current)
                else:
                    current._append_line(line.rstrip())
        for genome in result:
            genome._convert_to_numpy()
        return result

    @staticmethod
    def create_from_fastq_string(fastq_string):
        result = []
        state = 'balast'
        current = None
        for line in fastq_string.split('\n'):
            if len(line) > 0 and line[0] == '@':
                current = Genome(line)
                result.append(current)
                state = 'sequence'
            elif state == 'sequence':
                current._append_line(line.rstrip())
                state = 'balast'
        for genome in result:
            genome._convert_to_numpy()
        return result

    @staticmethod
    def load_from_fastq(filename):
        with open(filename, 'r') as file:
            return Genome.create_from_fastq_string(file.read())
