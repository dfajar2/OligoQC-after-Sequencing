#!/bin/bash
set -eo pipefail
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${PATH}
export PATH=/media/software/samtools/1.9/bin:${PATH}
export PATH=/media/software/bedtools/2.27.1/bin:${PATH}
export PATH=/media/software/bedops/2.4.35/bin:${PATH}
export PATH=/media/data/diego/scripts:${PATH}
export PATH=/media/software/fastx_toolkit/0.0.14/bin:${PATH}
export PATH=/usr/games:${PATH}

echo Running pipeline; echo ==================
start=`date +%s`
echo $start

usage="

This script will extract mapped reads not mapping to SNPs locations and 
count the most frequent occurrence of the first n bases from each read.
The script will create by default a padded file 100bases +/- the SNP 
location. 

USAGE:
    $ /path/to/script.sh OPTIONS
        Required:

        [ -a Full path to alignments directory ]
        [ -t Target file. bed format. ]
        [ -p Full path to primers fasta file ]
        [ -r Reference file Blast DB ]
        [ -o Output directory  ]

	Optional:
        [ -n increase by n bases the bed regions. Default: 100 bases ]
        [ -b Number of bases to trim. Default: 18 bases. ]
        [ -v verbose ]
"

n=100 
b=18
m=2000
I=0
v=""
#while getopts "a:t:p:r:o:n:b:m:Iv" options; do
while getopts "a:t:p:r:o:n:b:m:v" options; do
        case "${options}" in
                a)
                        a=${OPTARG} ;;
                t)
                        t=${OPTARG} ;;
                p)
                        p=${OPTARG} ;;
                r)
                        r=${OPTARG} ;;
		o)
			o=${OPTARG} ;;
		n)
			n=${OPTARG} ;;
		b)
			b=${OPTARG} ;;
		m)
			m=${OPTARG} ;;
#                I)
#			I=1 ;;
		v)
			set -xv ;;
		*)
                        echo ${usage}
                        exit 1 ;;
        esac
done

#START LOG
echo "Log file: " | tee ${o}.log


shift $((OPTIND-1))
if [ ${I} -eq 1 ]; then
	echo "IMPORTANT -- ls Running insilico only" ;
	if [ -z "${p}" ] ; then
		echo "ERROR - Missing primers file - fasta format."; echo "$usage"; exit 1
	elif [ -z "${r}" ] ; then
		echo "ERROR - Missing reference blast DB name."; echo "$usage"; exit 1
	elif [ -z "${o}" ] ; then
		echo "ERROR - Missing output directory."; echo "$usage"; exit 1
	fi
elif [ ${I} -eq 0 ]; then
	if [ -z "${a}" ]; then
		echo "ERROR - Missing alignments directory path."; echo "$usage"; exit 1
	elif [ -z "${t}" ]; then
		echo "ERROR - Missing SNP file - bed format."; echo "$usage"; exit 1
	elif [ -z "${p}" ] ; then
		echo "ERROR - Missing primers file - fasta format."; echo "$usage"; exit 1
	elif [ -z "${r}" ] ; then
		echo "ERROR - Missing reference blast DB name."; echo "$usage"; exit 1
	elif [ -z "${o}" ] ; then
		echo "ERROR - Missing output directory."; echo "$usage"; exit 1
	fi
fi

alignments=$a
targets=$t
primers=$p
reference=$r
outdir=$o
number=$n
bases=$b
merge=$m
base_dir=$(pwd)

if [ ${I} -eq 0 ]; then
	echo "Alignments Dir    = "${alignments} | tee -a ${o}.log
	echo "Targets File      = "${targets}    | tee -a ${o}.log
fi

echo "Primers File      = "${primers}    | tee -a ${o}.log
echo "Reference DB      = "${reference}  | tee -a ${o}.log
echo "Output Dir        = "${outdir}     | tee -a ${o}.log
echo "Padded bases      = "${number}     | tee -a ${o}.log
echo "Extracted bases   = "${bases}      | tee -a ${o}.log
echo "Distance insilico = "${merge}      | tee -a ${o}.log
#Check paths

if [ ! -d ${alignments} ]; then
	echo ----- | tee -a ${o}.log ; echo ERROR !!! Alignments Directory do not exist. Check files and paths.  | tee -a ${o}.log ; exit 1
elif [ ! -f ${targets} ] ;then
	echo ----- | tee -a ${o}.log ; echo ERROR !!! SNP Bed File does not exist. Check files and paths.  | tee -a ${o}.log ; exit 1
elif [ ! -f ${primers} ] ;then
	echo ----- | tee -a ${o}.log ; echo ERROR !!! Primers Fasta File does not exist. Check files and paths.  | tee -a ${o}.log ; exit 1
elif [ ! -f ${reference}*nhr ] ;then
	echo ----- | tee -a ${o}.log ; echo ERROR !!! Reference blast DB File does not exist. Check files and paths.  | tee -a ${o}.log ; exit 1
fi

# Check primers file name
if [[ "$primers" == *fasta  ]]
then
	echo Reference is in fasta format  | tee -a ${o}.log 
elif [[ "$primers" == *fa ]]
then
	echo Reference is in fasta format  | tee -a ${o}.log 
elif [[ "$primers" == *fna ]]
then
	echo Reference is in fasta format  | tee -a ${o}.log 
else
	echo What the hell??? Primers file does not have .fasta, .fa nor .fna suffix | tee -a ${o}.log 
	exit 1
fi

echo ============ | tee -a ${o}.log 
if [ ${bases} -lt 0 ]; then
	echo Invalid number of bases to pad bed file. | tee -a ${o}.log
	echo ${bases} bases is an invalid value.  | tee -a ${o}.log
	exit 1;
elif [ ${merge} -lt 0 ]; then
	echo Invalid number of bases to merg blast hits | tee -a ${o}.log
	echo ${merge} bases is an invalid value.  | tee -a ${o}.log
	exit 1;
fi
	
echo ============ | tee -a ${o}.log 
if [ ${I} -eq 0 ] ; then
	ls ${alignments}/*bam | tee -a ${o}.log 
fi

if [ ${number} -gt 0 ];
then
	echo ============ | tee -a ${o}.log ; echo Bed file will be padded by ${number} bases | tee -a ${o}.log 
else
	echo ============ | tee -a ${o}.log ; echo Bed file will not be padded  | tee -a ${o}.log 
fi

echo ============ | tee -a ${o}.log 
echo Output Directory : ${outdir}

#Safe stop
echo --------------- | tee -a ${o}.log 
echo Check files and paths | tee -a ${o}.log 
echo Press Y to continue, any other key to exit. | tee -a ${o}.log 
echo --------------- | tee -a ${o}.log 
read input

if [ "$input" != "Y" ] && [ "$input" != "y" ]; 
then
	echo Exiting... | tee -a ${o}.log   ; 
	echo | tee -a ${o}.log 
	exit 0
fi


if [ -d ${outdir} ]
then
	echo  | tee -a ${o}.log ; echo ---------------; echo WARNING; echo; echo Output folder ${o} already exist!! | tee -a ${o}.log 
	echo Do you want to delete the output folder and run the analysis again? \<Type Y to continue\> | tee -a ${o}.log 
	read input2
	if [ "$input2" != "Y" ] && [ "$input2" != "y" ];
	then
		cd ${outdir}
		echo Current directory: $(pwd)  | tee -a ${base_dir}/${o}.log
	else
		rm -r ${outdir}
		mkdir ${outdir}
		cd ${outdir}
	fi
else
	mkdir ${outdir}
	cd ${outdir}
fi

outdir=$(pwd)


function primersdb(){
	echo =========== | tee -a ${base_dir}/${o}.log
	echo Generate Primers db | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	cp ${primers} . 
	primers_file=$(ls ${primers} | awk -F/ '{print $NF}')
	primers_filename=$(echo ${primers_file} | awk -F".fasta" '{print $1}')
	#ln -s ${primers} .
	primers=${primers_file}
	echo Creating blast DB  | tee -a ${base_dir}/${o}.log
	makeblastdb -in ${primers} -dbtype nucl
	d=${primers}
	
#	echo Check if Reference db exist | tee -a ${base_dir}/${o}.log
#	if [[ "$reference"*nhr == *nhr  ]]
#	then
#		echo Reference has blast index | tee -a ${base_dir}/${o}.log 
#	else
#		echo ERROR - CANNOT CONTINUE!!! - Index Reference file first.  | tee -a ${base_dir}/${o}.log
#		exit 1
#	fi
}

function alignment(){
	echo =========== | tee -a ${base_dir}/${o}.log
	echo Running - function alignment  | tee -a ${base_dir}/${o}.log
	cd ${outdir}
	mkdir alignments ; cd alignments 
	for w in `ls ${alignments}/*.bam | grep -v Unde` ; do ln -s $w . ; done
	cd ..
	
	echo Preparing bed file  | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	if [ ${number} -gt 0 ];
	then
		sort-bed ${targets} | bedops --range ${number} --everything - > regions.bed 
	else
		cp ${targets} regions.bed
	fi
	
	targets="regions.bed"
	
	echo Generate bam of padded regions  | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	for f in `ls alignments/*bam | awk -F"." '{print $1}'`; do echo $f | tee -a ${base_dir}/${o}.log; bedtools intersect -v -a ${f}.bam -b ${targets} -sorted > ${f}_bedtools_v_padded_regions_offtarget.bam ; done
	echo | tee -a ${base_dir}/${o}.log
	echo Extract properly paired from off target regions | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	for f in `ls alignments/*_bedtools_v_padded_regions_offtarget.bam | awk -F"." '{print $1}'`; do echo $f | tee -a ${base_dir}/${o}.log; samtools view -b -f1 -F12 ${f}.bam > ${f}_properly_paired.bam ; done
	echo | tee -a ${base_dir}/${o}.log
	echo Generate R1 and R2 fastq files  | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	for f in `ls alignments/*_properly_paired.bam | awk -F"." '{print $1}'`; do echo $f | tee -a ${base_dir}/${o}.log; bamToFastq -i ${f}.bam -fq ${f}_R1.fastq -fq2 ${f}_R2.fastq ; done
	echo | tee -a ${base_dir}/${o}.log
	echo Generate fasta interleave  | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	for f in `ls alignments/*_properly_paired_R1.fastq | awk -F"_R1.fastq" '{print $1}'`; do echo $f | tee -a ${base_dir}/${o}.log; interleave-fastq ${f}_R1.fastq ${f}_R2.fastq | sed -n '1~4s/^@M/>M/p;2~4p' > ${f}_pairs_full_read.fasta ; done
	echo | tee -a ${base_dir}/${o}.log
	echo Extracting n bases from fasta files | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	for f in `ls alignments/*_pairs_full_read.fasta | awk -F".fasta" '{print $1}'`; do echo $f | tee -a ${base_dir}/${o}.log; fastx_trimmer -f 1 -l ${bases} -m ${bases} -i ${f}.fasta -o ${f}_trimmed_${bases}nt.fasta; done
	cp alignments/*_trimmed_${bases}nt.fasta .
}

function offtarget_to_primersDB(){
	echo =========== | tee -a ${base_dir}/${o}.log
	echo Running - function offtarget_to_primersDB  | tee -a ${base_dir}/${o}.log
	#Blast reads to primers DB
	cd ${outdir} 
	echo | tee -a ${base_dir}/${o}.log
	echo Blast trimmed reads against primers DB | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	for q in `ls *_trimmed_${bases}nt.fasta | awk -F"." '{print $1}'` ; \
		do \
			echo $q  | tee -a ${base_dir}/${o}.log; \
			#blastn -task=blastn -db ${primers} -query ${q}.fasta -outfmt "6 std qlen slen" -num_threads=36 > ${q}_blast_2pdb_temp.out ; \
			blastn -dust no -task blastn-short -db ${primers} -query ${q}.fasta -outfmt "6 std qlen slen" -num_threads=36 > ${q}_blast_2pdb_temp.out ; \
	done
	for blasts in `ls *_blast_2pdb_temp.out | awk -F"_blast_2pdb_temp" '{print $1}'` ; \
		do \
			echo $blasts  | tee -a ${base_dir}/${o}.log; \
			echo filter blast output by pct, length and best hit \(contains 1 and 2 freq\)
			cat ${blasts}_blast_2pdb_temp.out | grep -v ^$ | awk -v b=$bases '$4>=0.8*b && $8>=b-1 {print}'| sed 's/ /\t/g' | sort | awk '!NF || !seen[$1]++' > ${blasts}_reads_to_primersDB_blast.out ; \
	done
	mkdir offtarget_to_primersDB
	mv *offtarget*.* offtarget_to_primersDB
	cd offtarget_to_primersDB
		
	#IMPORTANT - To check!!!
# head offtarget_to_primersDB/merged_orig_bam_sorted_bedtools_v_padded_regions_offtarget_properly_paired_pairs_full_read_trimmed_18nt_reads_to_primersDB_blast.out
# M02484:311:000000000-D7C4P:1:1102:20156:16053/1 cp257_3580_F    100.000 18      0       0       1       18      1       18      2.54e-05        33.7    18      18
# M02484:311:000000000-D7C4P:1:1102:20156:16053/2 cp257_2698_F    100.000 17      0       0       1       18      1       17      8.87e-05        31.9    18      18
# M02484:311:000000000-D7C4P:1:1101:12872:17285/1 cp257_3484_R    94.444  18      1       0       1       18      1       18      0.001   28.3    18      18

	#1st only get best hits individually | tee -a ${base_dir}/${o}.log
	
	for hits in `ls *_reads_to_primersDB_blast.out | awk -F"." '{print $1}'`; do echo ${hits}; awk '{print $2}'  ${hits}.out | sort | uniq -c | sort -rn > ../${hits}_counts.txt ; done
	
		
	# create file with only 2
	echo create file with only 2 | tee -a ${base_dir}/${o}.log
	for l in `ls *_reads_to_primersDB_blast.out | awk -F. '{print $1}'`; do for f in `awk -F/ '{print $1}' ${l}.out | sort | uniq -c | awk '$1==2{print $2}'` ; do grep ${f} ${l}.out ; done > ${l}_only_pairs.out ; done
	
	#Counting Pairs
	echo Counting Pairs	| tee -a ${base_dir}/${o}.log
	for f in `ls *_only_pairs.out | awk -F. '{print $1}'` ; do awk '{if (NR%2==0) print $0,"lololo"; else print $0}' ${f}.out | tr "\n" "\t" | sed 's/lololo/\n/g' | sed 's/ /\t/g'| sed 's/\t\t/\t/g'| sed 's/\t\t/\t/g'| sed 's/\t\t/\t/g'|awk -F" " '{print $2,$16}' | awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }' | sort | uniq -c | sort -rn > ${f}_pairs_counts.txt; done
	
	#Counting individually from pairs
	echo Counting by individually from pairs| tee -a ${base_dir}/${o}.log
	for f in `ls *_only_pairs.out | awk -F. '{print $1}'` ; do awk '{print $2}' ${f}.out | sort | uniq -c | sort -k1rn > ${f}_single_counts.txt ; done
	for f in `ls *_only_pairs.out | awk -F. '{print $1}'` ; do awk '{print $2}' ${f}.out | sort | uniq -c | sort -k1rn > temp_lolo; SUM=$(awk '{sum+=$1}END{print sum}' temp_lolo); awk '{print $2}' ${f}.out | sort | uniq -c | sort -k1rn | awk -v all=$SUM '{print $1/all"\t"$0}' > ${f}_pct.txt; rm temp_lolo; done

#	rename 's/_sorted_bedtools_v_padded_regions_offtarget_properly_paired_pairs_full_read_trimmed_${bases}nt_reads_to_primersDB_blast//g' *txt
	for f in *txt; do mv "$f" "$(echo "$f" | awk -F "_sorted" '{print $1}')_sorted.txt" ; done
	mkdir ../Output_files_to_review
	cp *_sorted.txt ../Output_files_to_review/
	cd ..
}

################
###  BELOW ....
### In silico only primers vs reference NOT USED IN THIS SOP
################


function primers_to_referenceDB(){
	echo =========== | tee -a ${base_dir}/${o}.log
	echo Running - function primers_to_referenceDB  | tee -a ${base_dir}/${o}.log
	#Blast primers to reference genome
	cd ${outdir} 
	echo | tee -a ${base_dir}/${o}.log
	echo Blast primers sequence against reference DB | tee -a ${base_dir}/${o}.log
	echo | tee -a ${base_dir}/${o}.log
	for q in ${primers} ; \
		do \
			echo $q | tee -a ${base_dir}/${o}.log ; \
			#blastn -db ${reference} -query ${primers} -outfmt "6 std qlen slen" -num_threads=36 -task=blastn > ${primers}_blast_2rdb_temp.out ; \
			blastn -dust no -task blastn-short -db ${reference} -query ${primers} -outfmt "6 std qlen slen" -num_threads=36 > ${primers}_blast_2rdb_temp.out ; \
	done
	
	for blasts in `ls *_blast_2rdb_temp.out | awk -F"_blast_2rdb_temp" '{print $1}'` ; \
		do \
			echo $blasts  | tee -a ${base_dir}/${o}.log; \
			# $13-1 account for mismatches allowed at the end
			cat ${primers}_blast_2rdb_temp.out | grep -v ^$ |  awk -v b=$bases '$4>=0.8*b && $8>=$13-1 {print}'| sed 's/ /\t/g'> ${primers_filename}_primers_to_referenceDB_blast.out ; \
	done
}

function count_primers_clusters(){
	echo =========== | tee -a ${base_dir}/${o}.log
	echo Running - function count_primers_clusters  | tee -a ${base_dir}/${o}.log
	
	cd ${outdir} 
	#Blast results to bed format 
	echo Blast results to bed format | tee -a ${base_dir}/${o}.log
	cat ${primers_filename}_primers_to_referenceDB_blast.out | awk '{if ($9<$10) print $2"\t"$9"\t"$10"\t"$1"\t"$2"\t+"; else print $2"\t"$10"\t"$9"\t"$1"\t"$2"\t-"}' | \
	bedtools sort -i - > ${primers_filename}_blast.bed 
	
	
	#New approach
			
	#Merge if primers hits are ${merge}bp or less 
	echo Merge if primers hits are ${merge}bp or less | tee -a ${base_dir}/${o}.log
	# Just 2 by pairs
	bedtools merge -i ${primers_filename}_blast.bed  -d ${merge} -c 4,6,1 -o distinct,collapse,count | awk '$NF==2 && ($3-$2)>150 {print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6}' | grep +,- | awk '{print $5}'| awk -F, '{if (NF==1)print $1","$1; else print $0}' | awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }' | sort | uniq -c | sort -rn > ${primers_filename}_just2_by_pairs_Most_frequent.txt
	# Just 2 split
	bedtools merge -i ${primers_filename}_blast.bed  -d ${merge} -c 4,6,1 -o distinct,collapse,count | awk '$NF==2 && ($3-$2)>150 {print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6}' | grep +,- | awk '{print $5}'| awk -F, '{if (NF==1)print $1","$1; else print $0}' | sed 's/,/\n/g' | awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }'| sort | uniq -c | sort -rn > ${primers_filename}_just2_individually_Most_frequent.txt
	# All regions >2 primers
	bedtools merge -i ${primers_filename}_blast.bed  -d ${merge} -c 4,6,1 -o distinct,collapse,count | awk '$NF>=2 && ($3-$2)>150 {print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6}' | grep +,- | awk '{print $5}'| awk -F, '{if (NF==1)print $1","$1; else print $0}' | sed 's/,/\n/g' | awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }'| sort | uniq -c | sort -rn > ${primers_filename}_all_individually_Most_frequent.txt

	mkdir primers_to_referenceDB
	mv ${primers_filename}* primers_to_referenceDB
	mv primers_to_referenceDB/*Most_frequent.txt .
	}




function iqr(){
echo File "$1"
Q1pos=`cat $1 | awk 'END{u=0.75*(NR+1) ; print u}'`
echo Q1pos "$Q1pos"

if [[ "$Q1pos" =~ ^[0-9]+$ ]]; then
	echo $Q1pos | awk '{print $0 }'
	Q1=`cat $1 | awk -v Q1_a=$Q1pos 'NR==Q1_a{print $1}'`
else
	Q1_1=`echo $Q1pos | awk -F. '{print $1}' `
	Q1_2=`echo $Q1pos | awk -F. '{print ($1)+1}' `
	echo $Q1_1 $Q1_2
	Q11=`cat $1 | awk -v Q1_a=$Q1_1 'NR==Q1_a{print $1}'`
	Q12=`cat $1 | awk -v Q1_b=$Q1_2 'NR==Q1_b{print $1}'`
	echo Q1 $Q12 $Q11
	Q1=`echo $(((Q12+Q11)/2))`
fi

echo $Q1
Q3pos=`cat $1 | awk 'END{u=0.25*(NR+1) ; print u}'`
echo Q3pos "$Q3pos"
if [[ "$Q3pos" =~ ^[0-9]+$ ]]; then
	echo $Q3pos | awk '{print $0 }'
	Q3=`cat $1 | awk -v Q3_a=$Q3pos 'NR==Q3_a{print $1}'`
else
	Q3_1=`echo $Q3pos | awk -F. '{print $1}' `
	Q3_2=`echo $Q3pos | awk -F. '{print ($1)+1}' `
	echo $Q3_1 $Q3_2
	Q31=`cat $1 | awk -v Q3_a=$Q3_1 'NR==Q3_a{print $1}'`
	Q32=`cat $1 | awk -v Q3_b=$Q3_2 'NR==Q3_b{print $1}'`
	echo Q3 $Q32 $Q31
	Q3=`echo $(((Q32+Q31)/2))`
fi

echo $Q3
echo
echo IQR ...
echo IQR = $Q3 - $Q1
IQR=`echo $((Q3-Q1))`
echo
echo "Q3 + ( 1.5 x IQR )"
UP_BOUND=`awk -v a=$Q3 -v b=$IQR 'BEGIN{print a+(b*1.5)}'`
echo $UP_BOUND
cat ${1} | awk -v limit=$UP_BOUND '$1>limit{print $0}' > ${1}_CANDIDATES_IQR.txt
}


if [ ${I} -eq 0 ] ; then
	primersdb
	alignment
	offtarget_to_primersDB
#	primers_to_referenceDB
#	count_primers_clusters
	echo
	echo IMPORTANT - Files to check are located in the Output_files_to_review/ directory.
	echo
	echo ===================
	echo Pipeline Complete!!
elif [ ${I} -eq 1 ]; then
	cp ${primers} . 
	primers_file=$(ls ${primers} | awk -F/ '{print $NF}')
	primers_filename=$(echo ${primers_file} | awk -F".fasta" '{print $1}')
	primers=${primers_file}
	primers_to_referenceDB
	count_primers_clusters
	for f in `ls *_all_individually_Most_frequent.txt` ; do echo Running IQR ${f} ; iqr ${f} ;done
	echo
	echo ===================
	echo Pipeline Complete!!	
fi



