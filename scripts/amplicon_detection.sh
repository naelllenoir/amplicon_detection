#!/usr/bin/bash

set -Eeo pipefail
trap 'echo [ERROR] Error on $0 line $LINENO with the exit code $? >&2' ERR


declare ALI
declare CHR
declare OUT
declare HINT

SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
declare -r SCRIPT_NAME


displayUsage(){
	echo "USAGE :
${SCRIPT_NAME} [options]

OPTIONS :
	-a, --alignment STO	Specify a STOCKHOLM alignment file to use to build profile # should be a choice between this or a single query if no hmmprofile wanted
	-t, --target FASTA	FASTA file containing target chromosomes
	-H, --hints FASTA	FASTA file of consensus protein sequence used by Augustus as a hint # need to make this optional
	-d, --interpro-db	Path to interproscan6 databases
#	-s, --alphafold		Path to alphafold models
	-o, --outdir		Output directory (default to output_files)
	-h, --help			show this help message and exit

DESCRIPTION :
${SCRIPT_NAME} identifies putative duplicated homologs.
/!\ As of 31/08/2026 - Python script necessary for locus reconstruction must be copied in the folder where this script is run # for now - needs further development to use git clone / make.
Extrinsic file for Augustus should also be in the folder where the script is run. /!\
Other options will be added (output filenames, options for the python script (expansion, span), options for augustus, etc...).
"
}

checkCmd(){
	local cmd
	cmd="$1"
	command -v "${cmd}" > /dev/null || { echo "Command ${cmd} not found. Check your installation." >&2 && exit 127 ; }
}

checkEnv(){
	checkCmd realpath
	checkCmd nhmmer
	checkCmd python3
	checkCmd bedtools
	checkCmd augustus
	checkCmd startAlign.pl
	checkCmd align2hints.pl
	checkCmd getAnnoFasta.pl
	checkCmd nextflow
	checkCmd singularity
	checkCmd foldseek
}


testArg(){
	if [[ $2 =~ ^[-] || -z $2 ]]; then
		echo "[ERROR] Missing argument for $1" >&2 && displayUsage && exit 1
	else
		echo "$2"
	fi
}

parseArgs(){
	if(( $# == 0 )); then
		echo -e "[ERROR] No argument provided\n" >&2 && displayUsage && exit 1
	fi
	while(( $# > 0 )); do
  		case $1 in
			-a|--alignment) ALI=$(testArg "$1" "$2"); shift;;
			-t|--target) CHR=$(testArg "$1" "$2"); shift;;
			-H|--hint) HINT=$(testArg "$1" "$2"); shift;;
			-d|--interpro-db) IP=$(testArg "$1" "$2"); shift;;
#			-s|--alphafold) AF=$(testArg "$1" "$2"); shift;;
			-o|--outdir) OUT_DIR=$(testArg "$1" "$2"); shift;;
			-h|--help) displayUsage && exit 0 ;;
			--) shift; break ;;
			-*) echo "[ERROR] ${SCRIPT_NAME} : unrecognized option $1" >&2 && displayUsage && exit 1 ;;
			*) break;;
		esac
			shift
	done
}

checkArgs(){
	if [[ -z ${ALI} ]]; then
		echo "[ERROR] -a|--alignment was not supplied (mandatory file)" >&2 && displayUsage && exit 1
	fi
	if [[ -z ${CHR} ]]; then
		echo "[ERROR] -t|--target was not supplied (mandatory file)" >&2 && displayUsage && exit 1
	fi
	if [[ -z ${HINT} ]]; then
		echo "[ERROR] -H|--hint was not supplied (mandatory file)" >&2 && displayUsage && exit 1
	fi
        if [[ -z ${IP} ]]; then
                echo "[ERROR] -d|--interpro-db was not supplied (mandatory path)" >&2 && displayUsage && exit 1
        fi
 #       if [[ -z ${AF} ]]; then
  #              echo "[ERROR] -s|--alphafold Alphafold models directory was not supplied (mandatory path)" >&2 && displayUsage && exit 1
   #     fi
	if [[ -z "${OUT_DIR}" ]]; then
		mkdir output_files
		OUT_DIR=output_files ;
	elif [[ -e ${OUT_DIR} ]]; then
		echo "${OUT_DIR} already exists. Please delete your existing directory or give another name to the directory you want to create." >&2 && exit
	elif [[ ${OUT_DIR} != output_files ]]; then
		mkdir ${OUT_DIR}
	fi

	ALI=$(realpath -s "${ALI}")
	CHR=$(realpath -s "${CHR}")
        HINT=$(realpath -s "${HINT}")
	IP=$(realpath -s "${IP}")
	OUT_DIR=$(realpath -s "${OUT_DIR}")
#	AF=$(realpath -s "${AF}")
}


hmm_profile(){
	echo -e "Building hmm profile...\n" >> "${LOG}"
	hmmbuild profile.hmm ${ALI} >> "${LOG}"
}

putative_exons(){
	echo -e "\nSearching for putative exons...\n" >> "${LOG}"
	nhmmer --dna -o nhmmer.txt -A msa_nhmmer.txt --tblout table_nhmmer.txt --seed 63 --cpu 10 ${OUT_DIR}/hmm_profile/profile.hmm ${CHR}
        grep -v "#" table_nhmmer.txt | sed -e 's/[[:space:]]\+/,/g' -e 's/,/\t/g' | awk '{ if ($12=="-") { tmp=$7; $7=$8; $8=tmp; tmp2=$9; $9=$10; $10=tmp2 } print }' OFS='\t' | sort -k7,7n | cut -f1-15 >> table_nhmmer_orderedCoordAliStrand_nodesc.txt
}

# Option to add: additional gff3 output file - default output file should be changed to BED
putative_loci(){
	echo -e "Reconstructing putative loci...\n" >> "${LOG}"
	echo -e "$(basename ${CHR%.*})\n---" >> "${NB}"
	python3 ~/amplicon_detection/scripts/loci_reconstruction.py -i ${OUT_DIR}/putative_exons/table_nhmmer_orderedCoordAliStrand_nodesc.txt -o ${OUT_DIR}/putative_loci/putative_loci_nhmmer.gff3 -f "${CHR}" -e 1000 -m 10000 >> "${NB}"
        grep "putative_locus_expanded" ${OUT_DIR}/putative_loci/putative_loci_nhmmer.gff3 | awk -F'\t' 'BEGIN{OFS="\t"}{ match($9,/ID=([^;]+)/,id); print $1,$4-1,$5,id[1],".",$7}' > ${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.bed
}

seq_prediction(){
	echo -e "Starting sequence prediction in putative loci...\n" >> "${LOG}"
	bedtools getfasta -fi ${CHR} -fo  ${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.fasta -bed ${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.bed -s -nameOnly
	startAlign.pl --genome=${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.fasta --prot=${HINT} --prg=exonerate
	align2hints.pl --in=align_exonerate/exonerate.concat.aln --out=prot.hints --prg=exonerate --genome_file=${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.fasta
	augustus --species=human --extrinsicCfgFile=~/amplicon_detection/scripts/extrinsic.MP.cfg --hintsfile=prot.hints --uniqueGeneId=true --alternatives-from-evidence=false --softmasking=off --codingseq=on --protein=on ${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.fasta > augustus_protein_hints_putative_loci_nhmmer_expanded.gff
	getAnnoFasta.pl augustus_protein_hints_putative_loci_nhmmer_expanded.gff
	awk '/^>/ {if (seq) print seq; print; seq=""} /^[^>]/ {seq=seq $0} END {print seq}' augustus_protein_hints_putative_loci_nhmmer_expanded.aa > augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa
	grep -v "#" augustus_protein_hints_putative_loci_nhmmer_expanded.gff | grep -w "gene" > augustus_protein_hints_putative_loci_nhmmer_expanded_gene.gff
	grep '>' augustus_protein_hints_putative_loci_nhmmer_expanded.aa | cut -d'>' -f2 | cut -d'(' -f1 | sort -t _ -k2 -g | xargs -I{} grep -w {} ${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.bed > augustus_protein_hints_putative_loci_nhmmer_expanded.bed
	paste augustus_protein_hints_putative_loci_nhmmer_expanded_gene.gff augustus_protein_hints_putative_loci_nhmmer_expanded.bed | awk 'BEGIN { OFS="\t" } { if ($15=="+") { start=$11+$4; end=$11+$5 } else { start=$12-$5; end=$12-$4 } print $1,$2,$3,start,end,$6,$15,$14,$9 }' > augustus_predicted_sequences.gff3
	awk 'BEGIN {dism=ARGV[2]-ARGV[1]; print "\n\nAUGUSTUS\n----------\nPutative protein-coding sequences predicted by AUGUSTUS:", ARGV[1]"\nPutative non-coding loci dismissed:", dism}' "$(grep '>' augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa | wc -l)" "$(grep '>'  ${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.fasta | wc -l)" >> "${NB}"
}

functional_annot(){
	echo -e "Running Interproscan6 on AUGUSTUS predicted sequences...\n" >> "${LOG}"
	nextflow run ebi-pf-team/interproscan6 -r 6.0.2.2 -profile singularity --input ${OUT_DIR}/augustus_predictions/augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa --datadir ${IP} --interpro 107.0 --applications SUPERFAMILY -ansi-log false >> "${LOG}"
	awk 'BEGIN {dism=ARGV[2]-ARGV[1]; print "\n\nInterproscan6\n----------\nPutative protein-coding sequences predicted by AUGUSTUS with functional domains identified:", ARGV[1]"\nPutative protein-coding sequences predicted by AUGUSTUS without functional domains identified:", dism}' "$(cut -f1 augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa.tsv | sort | uniq | wc -l)" "$(grep '>' ${OUT_DIR}/augustus_predictions/augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa | wc -l)" >> "${NB}"
}

# Need to add a default parameter for min length (0)
annot_filter(){
	echo -e "Filtering Interproscan6 results to only keep sequences with a predicted NAP domain...\n" >> "${LOG}"
	grep -e "IPR002164\|IPR037231" augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa.tsv | cut -f1 | sort | uniq | xargs -I{} grep -A1 {} ${OUT_DIR}/augustus_predictions/augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa > augustus_protein_hints_putative_loci_nhmmer_expanded_filt_annotation.aa
	grep '>' augustus_protein_hints_putative_loci_nhmmer_expanded_filt_annotation.aa | cut -d'>' -f2 | cut -d'(' -f1 | sort -t _ -k2 -g | xargs -I{} grep -w {} ${OUT_DIR}/putative_loci/putative_loci_nhmmer_expanded.bed > nap_loci.bed
	grep -oe 'locus_[0-9]*' nap_loci.bed | xargs -I{} grep -w {} ${OUT_DIR}/augustus_predictions/augustus_predicted_sequences.gff3  | awk '{len=$5-$4+1} {print $0";length="len}' > nap_loci.gff3
	grep "SSF143113" augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa.gff3 | awk -v min=120 '{start=$4; end=$5; size=end-start} { if (size>=min) {print $0,"size="size} }' | cut -d'.' -f1,2 | xargs -I{} grep -w {} ${OUT_DIR}/augustus_predictions/augustus_predicted_sequences.gff3 > "$(basename ${CHR%.*})"_size_filtered_nap_proteins.gff3
	cut -f9 "$(basename ${CHR%.*})"_size_filtered_nap_proteins.gff3 | xargs -I{} grep -A1 -w {} augustus_protein_hints_putative_loci_nhmmer_expanded_filt_annotation.aa > "$(basename ${CHR%.*})"_size_filtered_nap_proteins.aa
	cut -f9 "$(basename ${CHR%.*})"_size_filtered_nap_proteins.gff3 | xargs -I{} grep -w '^{}' augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa.gff3 > "$(basename ${CHR%.*})"_size_filtered_nap_proteins_annotation.gff3
	cat "$(basename ${CHR%.*})"_size_filtered_nap_proteins.aa >> "$(basename ${CHR%.*})"_size_filtered_nap_proteins_annotation.gff3
	awk 'BEGIN {dism=ARGV[2]-ARGV[1]; print "\nPutative protein-coding sequences predicted by AUGUSTUS encoding a NAP domain:", ARGV[1]"\nPutative protein-coding sequences predicted by AUGUSTUS not encoding a NAP domain:", dism}' "$(grep '>' augustus_protein_hints_putative_loci_nhmmer_expanded_filt_annotation.aa | wc -l)" "$(cut -f1 augustus_protein_hints_putative_loci_nhmmer_expanded_oneliner.aa.tsv | sort | uniq | wc -l)" >> "${NB}"
	awk 'BEGIN {dism=ARGV[2]-ARGV[1]; print "\nPutative NAP domain coding sequences passing the NAP domain size threshold:", ARGV[1]"\nPutative NAP domain containing non-coding sequences (not passing the NAP domain size threshold):", dism}' "$(grep '>' $(basename ${CHR%.*})_size_filtered_nap_proteins.aa | wc -l)" "$(grep '>' augustus_protein_hints_putative_loci_nhmmer_expanded_filt_annotation.aa | wc -l)" >> "${NB}"
}

struct_json(){
	for loc in $(grep '>' "$(basename ${CHR%.*})"_size_filtered_nap_proteins.aa); do
		seqname=${loc:1}
		seq=$(awk -v loc=$seqname '$0 == ">"loc {getline;print}' "$(basename ${CHR%.*})"_size_filtered_nap_proteins.aa)
		output=$(echo ${seqname} | grep -oE 'locus_[0-9]+\([+-]\)\.g[0-9]+\.t[0-9]+')
		echo -e '{\n  "name": "'${seqname}'",\n  "sequences": [\n    {\n      "protein": {\n        "id": "A",\n        "sequence": "'${seq}'",\n        "description": "'${seqname}'"\n      }\n    }\n  ],\n  "modelSeeds": [63],\n  "dialect": "alphafold3",\n  "version": 4\n}' > ${OUT_DIR}/structure/json_files/"$(basename ${CHR%.*})"_${output}.json
	done
}

# Will need to add another argument to the script for Alphafold databases and image ad script..
#struct_prediction(){
#parallel "singularity exec --nv --bind ${OUT_DIR}/structure/json_files:/root/af_input --bind ${OUT_DIR}/structure/alphafold_predictions:/root/af_output --bind ${AF}:/root/models --bind /storage/biodatabanks/alphafold3/2026-07-28:/root/public_databases $HOME/alphafold3/alphafold3.sif python $HOME/alphafold3/run_alphafold.py --json_path=/root/af_input/{} --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output" ::: ${OUT_DIR}/structure/json_files/*.json
#}

# Need to add another argument to have one or more reference structure to use for foldseek comparisons - ${SREF} (.cif format)
#struct_comparison(){
#foldseek createdb "${OUT_DIR}"/structure/alphafold_predictions/ "${OUT_DIR}"/structure/folseek_comparisons/folseek_db
#foldseek easy-search ${SREF} "${OUT_DIR}"/structure/folseek_comparisons/folseek_db foldseek_comparison.m8 tmp --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,prob"
#foldseek easy-search ${SREF} "${OUT_DIR}"/structure/folseek_comparisons/folseek_db foldseek_comparison.html tmp --format-mode 3
#}


main(){
	checkEnv
	parseArgs "$@"
	checkArgs
	mkdir -p "${OUT_DIR}"
	cd "${OUT_DIR}"
	NB="${OUT_DIR}"/"$(basename ${CHR%.*})"_counts.txt
	LOG="${OUT_DIR}"/detection_pipeline.log
        mkdir -p "${OUT_DIR}"/hmm_profile
        cd "${OUT_DIR}"/hmm_profile
	hmm_profile
        mkdir -p "${OUT_DIR}"/putative_exons
        cd "${OUT_DIR}"/putative_exons
	putative_exons
        mkdir -p "${OUT_DIR}"/putative_loci
        cd "${OUT_DIR}"
	putative_loci
        mkdir -p "${OUT_DIR}"/augustus_predictions
        cd "${OUT_DIR}"/augustus_predictions
	seq_prediction
        mkdir -p "${OUT_DIR}"/functional_domains
        cd "${OUT_DIR}"/functional_domains
	functional_annot
	annot_filter
	mkdir -p "${OUT_DIR}"/structure/json_files
	struct_json
#	cd "${OUT_DIR}"/structure
#       mkdir -p "${OUT_DIR}"/structure/alphafold_predictions
#	cd "${OUT_DIR}"/structure/alphafold_predictions
#	struct_prediction
#	mkdir -p "${OUT_DIR}"/structure/folseek_comparisons
#	cd "${OUT_DIR}"/structure/folseek_comparisons
#	struct_comparison
}

main "$@"
