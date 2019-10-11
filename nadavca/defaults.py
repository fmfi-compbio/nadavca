import os
import nadavca

KMER_MODEL_FILE = os.path.join(os.path.dirname(nadavca.__file__), 'default', 'kmer_model.hdf5')
CONFIG_FILE = os.path.join(os.path.dirname(nadavca.__file__), 'default', 'config.yaml')
BWA_EXECUTABLE = 'bwa'
GROUP_NAME = 'Analyses/Basecall_1D_000'
