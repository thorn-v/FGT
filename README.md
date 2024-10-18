# FGT

Pairwise Four Gamete Test (FGT) for haploid bi-allelic SNPs.  
(Tested up to 128 individuals/samples and 100,000s of loci. The more loci included, the longer the program will take.)

### Compile

To compile, use the `make` command in the folder with main.c and Makefile.  
To run, use `./bin/4gamete [options] input.txt`

### Input

Expects files with a binary SNP array in the form of loci as rows (with row name column) and samples as columns.   Tab or Space delim.  
Ex.

```text
loci_1 0 1 0 0 1 0 
loci_2 1 1 1 0 1 0 
loci_3 0 0 0 1 1 1 
...
```

This can easily be acheived using bcftools:

`bcftools query -f "%CHROM\_%POS[ %GT]\n" input.vcf.gz > binary_alleles_output.txt`

### Output

This program will write pairs found to stdout by default. To keep the pairs, redirect stdout to a file.  
It will also write statistics to the screen (stderr) including time taken and number of matches.

### Options

-j  number of cores to use (defaults to 1 core with no flag)  
-n  supress written output - do not write out matches to stdout (faster, good for debugging or to get quick statistics)
