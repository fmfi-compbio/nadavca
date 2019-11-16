import h5py
from nadavca.alphabet import alphabet, inv_alphabet
from nadavca.dtw import KmerModel
import numpy as np

def kmer_to_id(kmer):
    result = 0
    for base in kmer:
        result *= len(alphabet)
        result += inv_alphabet[base]
    return result

@staticmethod
def load_kmer_model(filename):
    with h5py.File(filename, "r") as file:
        central_position = file.attrs["central_pos"]
        model_table = file["model"][()]
        mean = np.zeros(len(model_table))
        sigma = np.zeros(len(model_table))
        k = None
        for line in model_table:
            kmer = line[0].decode("ascii")
            kmer_id = kmer_to_id(kmer)
            mean[kmer_id] = line[1]
            sigma[kmer_id] = line[2]
            k = len(kmer)
        return KmerModel(k, central_position, len(alphabet), mean, sigma)

setattr(KmerModel, "load_from_hdf5", load_kmer_model)
