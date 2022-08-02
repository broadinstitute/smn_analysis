# SMN_ANALYSIS Tool
smn_analysis is a tool to detect Survival Motor Neuron (SMN) mutations to diagnose Spinal Muscular Atrophy (SMA) by counting the total number of reads and the number of C reads at c.840 position in exon 7.
This tool also works for 

## Dependencies
- Python3
- Python Libraries
- Pandas
- Matplotlib
- Pysam
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- [htslib (bgzip & tabix) and samtools](http://www.htslib.org/download/)
- [gcloud](https://cloud.google.com/sdk/docs/install)
- [gsutil](https://cloud.google.com/storage/docs/gsutil_install#install)
- ...
	
## Installation
Use git to clone the package and install 
```
$ git clone git@github.com:broadinstitute/smn_analysis.git

$ cd smn_analysis/
```
## Usage
```
$ python3 smn_analysis.py -R --reference-fasta -o --output-tsv -v --verbose cram_paths
```

Example: 
```
$ python3 smn_analysis_C.py -R /Users/r/SMA.2022/Homo_sapiens_assembly38.fasta -o c.tsv -v ../MAAC002.cram
```
## Output
Based on the -o/--output-tsv parameter, smn_analysis will output a tsv file with unknown, unaffected and affected samples with SMA.





