# Checkphase
Provided a reference and a query VCF/BCF file, Checkphase compares shared biallelic variants from both files and gives a summary on switch errors (useful for haplotype phasing) and genotype errors (useful for genotype imputation). The tool uses Htslib to load the VCFs and to identify shared sites.

Genotype errors can be presented as *hard* or *soft* errors. Hard errors are based on the hard genotype call while soft errors are based on the imputation dosage. Soft errors are only calculated if allele dosage data is present (i.e. either ADS or HDS tags).

**Note:** The reference and the query files must contain exactly the same samples.

**Note:** A different phase in both input files at a heterozygous site is only considered as a switch error if the phase at the previous shared heterozygous site was equal. Otherwise, the path simply continues and the difference is therefore not another switch error. Futhermore, a phase difference at the first heterozygous position is considered as a switch in the maternal and paternal path (mat/pat switch) and not as a switch error.

The `--stat` switch continuously writes the progress in a file with the provided filename.

The `--dump` switch dumps for each sample the phase error positions to stderr.

### Example output
Comparing two phasing outputs:
```
$ checkphase ref.vcf.gz query.vcf.gz

Input files:
  Reference: ref.vcf.gz
  Query:     query.vcf.gz
  Reference variants (from index): Mref = 8799
  Reference samples:               Nref = 500
  Query variants (from index):   Mquery = 8799
  Query samples:                 Nquery = 500

Passed samples check.
Checking shared variants for phase and genotype errors... done.

Summary:

  Haploid ref samples:      0
  Haploid query samples:    0

  Reference variants:       8799
  Query variants:           8799
  Shared variants:          8799
  Checked variants:         8799

  Excluded shared:          0
    Not biallelic in ref:   0
    Not biallelic in query: 0
    Alleles do not match:   0

  Ref/Alt swaps:            0
  Strand flips:             0
  Ref/Alt + Strand flip:    0

  Missing or unphased sites:
    Total missing in ref:        0
    Total missing in query:      0
    Total unphased het in ref:   0
    Total unphased het in query: 0

  Genotype errors (hard):
    Total genotype errors:     0
    Minimum genotype errors:   0
    Maximum genotype errors:   0
    Average gt err per sample: 0
    Minimum gt error rate:     0
    Maximum gt error rate:     0
    Average gt error rate:     0
    Standard GER deviation:    0
    GER variance:              0

  Switch errors:
    Total switch errors:       27855
    Minimum switch errors:     4
    Maximum switch errors:     379
    Average sw err per sample: 56
    Minimum sw error rate:     0.000454597
    Maximum sw error rate:     0.0430731
    Average sw error rate:     0.0063314
    Standard SER deviation:    0.00318695
    SER variance:              0.000220513
    Switch error free targets: 0
    Mat/Pat switches:          251
```

## Summarize information
The bash script `checkphase_extract` can be used to generate a CSV file from several outputs of Checkphase, e.g. from several chromosome files.

In addition, the bash script `checkphase_summary` summarizes the outputs of several Checkphase runs.

Both scripts use AWK in the background.
