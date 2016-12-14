# BLAST\_to\_BED

BLAST\_to\_BED is a simple Python script that parses a BLAST XML file and creates a BED file. One can give it a BLAST XML file to be parsed or give it a FASTA file and path to a nucleotide BLAST database to run BLASTn and parse the results. Basic usage can be found as follows:

```
$ ./BLAST_to_BED.py
usage: BLAST_to_BED.py (-f FASTA FILE | -x XML FILE)
                       [-b BED FILE | -o OUTPUT FILE]
                       [-d REFERENCE BLAST DATABASE] [-e E-VALUE THRESHOLD]
                       [-s MAX HITS] [-m MAX HSPS] [-k]
                       [-X [EXCLUDE CHROMOSOMES [EXCLUDE CHROMOSOMES ...]] |
                       -K [KEEP ONLY CHROMOSOMES [KEEP ONLY CHROMOSOMES ...]]]

Input options:
  Specify the type of input we're working with, '-f | --fasta' requires BLASTn from NCBI to be installed

  -f FASTA FILE, --fasta FASTA FILE
                        Input FASTA file to run BLAST on, incompatible with
                        '-x | --xml'
  -x XML FILE, --xml XML FILE
                        Input BLAST XML file to turn into BED file,
                        incompatible with '-f | --fasta'

Output options:
  Specify how we're writing our output files, create new file or append to preexisting BED file

  -b BED FILE, --bed BED FILE
                        BED file to append results to, incompatible with '-o |
                        --outfile'
  -o OUTPUT FILE, --outfile OUTPUT FILE
                        Name of output file to write to, defaults to
                        '/home/paul/BLAST_to_BED/output.bed', incompatible
                        with '-b | --bed'

BLAST options:
  Options only used when BLASTING, no not need to provide with '-x | --xml'

  -d REFERENCE BLAST DATABASE, --database REFERENCE BLAST DATABASE
                        Reference BLAST database in nucleotide format, used
                        only when running BLAST
  -e E-VALUE THRESHOLD, --evalue E-VALUE THRESHOLD
                        Evalue threshold for BLAST, defaults to '1e-1', used
                        only when running BLAST
  -s MAX HITS, --max-hits MAX HITS
                        Maximum hits per query, defaults to '1', used only
                        when running BLAST
  -m MAX HSPS, --max-hsps MAX HSPS
                        Maximum HSPs per hit, defaults to '1', used only when
                        running BLAST
  -k, --keep-xml        Do we keep the XML results? pas '-k | --keep-xml' to
                        say 'yes', used only when running BLAST

Filtering options:
  Options to filter the resulting BED file, can specify as many or as few chromosome/contig names as you want. You may choose to either exclude or keep, not both. Note: this only affects the immediate output of this script, we will not filter a preexisting BED file

  -X [EXCLUDE CHROMOSOMES [EXCLUDE CHROMOSOMES ...]], --exclude-chrom [EXCLUDE CHROMOSOMES [EXCLUDE CHROMOSOMES ...]]
                        Chromosomes/Contigs to exclude from BED file,
                        everything else is kept, incompatible with '-K |
                        --keep-chrom'
  -K [KEEP ONLY CHROMOSOMES [KEEP ONLY CHROMOSOMES ...]], --keep-chrom [KEEP ONLY CHROMOSOMES [KEEP ONLY CHROMOSOMES ...]]
                        Chromosomes/Contigs to keep from BED file, everything
                        else is excluded, incompatible with '-X | --exclude-
                        chrom'
```

## Inputs

BLAST\_to\_BED requires either a FASTA file or BLAST XML file as input. These two options are mutually exclusive; one cannot give both a FASTA and XML file. Giving a FASTA file will cause BLAST\_to\_BED to run BLASTn against a BLAST database. The resulting XML file will be removed upon completion of the script unles `-k | --keep-xml` is passed to the script

## Outputs

BLAST\_to\_BED creates between one and four output files.

| File name | Contents |
| --------- | -------- |
| *fasta*_*database*.xml | XML results from running BLAST with BLAST\_to\_BED. Generated only if `-k | --keep-xml` is passed to the script.
| *output*.bed | Final four-column BED file (last column is query name), specified by `-o | --outfile` or `-b | --bed` |
| *output*_failed.log | List of queries that failed during the BLAST, generated only if there were failures |
| *output*_filtered.bed | BED file that contains the removed hits from the filtering step, generated only if filtering was specified |

## Filtering

BLAST\_to\_BED can filter the final BED output to exclude or keep only certain chromosomes or contigs, though not both for obvious reasons. To filter, pass `-X | --exclude-chrom` to exclude or `-K | --keep-chrom` to keep and a space-delimited list of chromosome or contigs names exactly as they would be found in the BLAST XML file. If you don't know what the chromosome contig names are, either look in the BLAST XML file (generated with `-k | --keep_xml`) or use grep on the reference FASTA file that was the basis for the BLAST database

```bash
grep '>' reference.fasta
```

## Dependencies

BLAST\_to\_BED depends on the following:
 - [Python 3](https://www.python.org/)
 - [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 - [BioPython](http://biopython.org/wiki/Biopython)
 - [Beautiful Soup 4](https://www.crummy.com/software/BeautifulSoup/)
 - [lxml](http://lxml.de/)

The latter three are all available on [PyPi](https://pypi.python.org/pypi) and can be downloaded using [pip3](https://pip.pypa.io/en/latest/installing/) (included with Python 3.4 or greater)
