# Nadavca
NAnopore DAta Variant CAller

## Installing

Run:

```
git clone git@github.com:baklazan/nadavca.git
cd nadavca
python3 setup.py build
python3 setup.py install
```

Depending on your configuration, the last step might require root privileges.

After installing, you should be able to run
```
nadavca --help
```

and also import `nadavca` package in python3:
```
python3
...
>>> import nadavca
```

If Nadavca runs, but crashes on unsuccessfully trying to import some package,
please report this as a bug and we will try to fix it. In the meantime, just
manually install that package and try again.

If Nadavca crashes for any other reason, please report that, too.

## Usage

### Calling SNPs

```python
nadavca.estimate_snps(reference_filename,
                      reads,
                      reference=None,
                      config='default/config.yaml',
                      kmer_model='default/kmer_model.hdf5',
                      bwa_executable='bwa',
                      independent=False,
                      group_name='Analyses/Basecall_1D_000')
```

where:

*  `reference_filename` is name of a `.fasta` file with the reference sequence
*  `reads` is either of:
   -  a list of filenames of individual reads
   -  a name of a directory. In that case, all `.fast5` files in the directory will
      be processed
   -  a list of `nadavca.Read` instances
*  If reference sequence has already been loaded into a `nadavca.Genome` instance,
   you can pass it as `reference` argument to prevent Nadavca from unnecessarily
   reloading it. You still need to specify `reference_filename`
*  `config` is either of:
   - a name of a YAML file containing parameters for the SNP-calling algorithm 
     (see [default values](default/config.yaml))
   - a `dict` containing parameters for the SNP-calling algorithm
*  `kmer_model` is either of:
   - a name of a HDF5 file containing expected values of signal for individual 
     k-mers
   - a `nadavca.KmerModel` instance
*  `bwa_executable` is command that runs BWA on your system
*  If `independent` is set to `True`, each read will be treated separately 
   (i.e. Nadavca assumes each read was from a different modification of 
   the reference sequence). Otherwise, information from all reads is combined
   into single consensus score for each position in the sequence.
*  `group_name` is path to the group containing results of basecalling inside
   the `.fast5` files of reads

The return value of `nadavca.estimate_snps()` is a list of 
`nadavca.estimator.Chunk`s. Each `Chunk` contains estimated SNP probabilities
for a contiguous segment of the reference sequence. Bounds of this segment are
indicated by `Chunk`'s `.start` and `.end` attributes.

Estimated probabilities themselves are in `Chunk`'s `.values` attribute,
which is a 2D `numpy` array of dimension `(end - start, 4)`. For each position
in reference sequence between `start` (inclusive) and `end` (exclusive),
it contains 4 numbers: the estimated probability of A, C, G and T, respectively,
on this position.
   
If `independent` was set to `False`, some positions in the reference sequence may
be covered by multiple reads. This information is stored in `.cover` attribute
of `nadavca.estimator.Chunk`.

If `independent` was set to `True`, Nadavca produces a `Chunk` for each read.
If you pass reads to `nadavca.estimate_snps()` as a list, while using `independent=True`,
the order of `Chunk`s corresponds to the order of reads.

There is also a command-line interface for `nadavca.estimate_snps()`:

```
nadavca snp reference reads
```

Run

```
nadavca --help
nadavca snp --help
```

for more details.

### Aligning signal to reference

```
nadavca.align_signal(reference_filename,
                     reads,
                     reference=None,
                     config='default/config.yaml',
                     kmer_model='default/kmer_model.hdf5',
                     bwa_executable='bwa',
                     group_name='Analyses/Basecall_1D_000'):
```

where parameters have the same meaning as in `nadavca.estimate_snps()`.

The return value of `nadavca.align_signal()` is a list of 2D `numpy` arrays.
Each array in the list corresponds to one read, in the same
order as they appeared on the input. Thus, it is advisable to pass reads to
`nadavca.align_signal()` as a list. 

Each array has two columns: the first
column contains positions in the signal, the second column contains 
corresponding positions in the reference sequence. Positions in signal are 
in ascending order. Positions in reference are in ascending or descending
order, for forward strands and reverse strands, respectively.

`nadavca.align_signal()` also has a command-line interface:

```
nadavca align reference reads
```

Run

```
nadavca --help
nadavca align --help
```

for more details.
