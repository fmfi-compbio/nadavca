import os
import nadavca

KMER_MODEL_FILE = os.path.join(os.path.dirname(nadavca.__file__), 'default', '10kmer_fact2.h5')
CONFIG_FILE = os.path.join(os.path.dirname(nadavca.__file__), 'default', 'config.yaml')
BWA_EXECUTABLE = 'bwa'
GROUP_NAME = 'Analyses/Basecall_1D_000'
RENORM_ROUNDS = 3
