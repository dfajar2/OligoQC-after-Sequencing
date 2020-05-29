# OligoQC-after-Sequencing

This script accepts alignment files or PE reads. If the input are alignment 
files, then the script will extract mapped reads not mapping to SNPs locations
and count the most frequent occurrence of the first n bases from each read.
The script will create by default a padded file 100bases +/- the SNP 
location. 

Alternatively, the script can use as input read files. The script will extract
n bases from each read and blast the sequences against the primers DB.

USAGE:
    $ /path/to/script.sh OPTIONS
        Required:

        [ -a Alignments directory ]
        [ -t Target file. bed format ]
        [ -r Reference file Blast DB ]

	or 

	      [ -s Reads directory ]

        [ -p Full path to primers fasta file ]
	      [ -o Output directory  ]

	Optional:
        [ -n increase by n bases the bed regions. Default: 100 bases ]
        [ -b Number of bases to trim. Default: 18 bases. ]
        [ -v verbose ]
