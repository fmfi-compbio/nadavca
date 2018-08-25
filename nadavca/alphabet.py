alphabet = ['A', 'C', 'G', 'T']
inv_alphabet = {c : i for i, c in enumerate(alphabet)}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
numerical_complement = {i : inv_alphabet[complement[c]] for i, c in enumerate(alphabet)}
