MAX_LIBS_NUMBER=5

get_orientation() { 
  gtable_id=$1
  orientation='null'

  for f in $(dx get_details $gtable_id); do
    if [ $orientation == 'now' ]; then        
      orientation=${f//\"/}  # removing double quotes
      break
    fi
    if [ $f == '"pair_orientation":' ]; then
      orientation='now'
    fi
  done

  echo ${orientation,,} # to lowercase
}

process_reads() {
  declare -a reads_array=("${!1}")
  type=$2  # pe, mp or single

  SPADES_READS=
  for index in ${!reads_array[@]}
  do
    if [ $index -eq $MAX_LIBS_NUMBER ]; then          
      break
    fi

    reads_gtable_name=$(dx-jobutil-parse-link "${reads_array[$index]}")
    if [ $type == 'single' ]; then 
      read_fname="unpaired_$index.fastq"

      dx-reads-to-fastq $reads_gtable_name --output $read_fname
      SPADES_READS="$SPADES_READS -s $read_fname"
    else
      lread_fname="left_${type}_${index}.fastq"
      rread_fname="right_${type}_${index}.fastq"

      dx-reads-to-fastq $reads_gtable_name --output $lread_fname --output2 $rread_fname
      SPADES_READS="$SPADES_READS --$type$(($index + 1))-1 $lread_fname --$type$(($index + 1))-2 $rread_fname"
      orientation=$(get_orientation "$reads_gtable_name")
      if ! [ $orientation == 'null' ]; then
          SPADES_READS="$SPADES_READS --$type$(($index + 1))-$orientation"
      fi
    fi
  done

  echo $SPADES_READS
}

main() {
  set -e -x

  # Look for the inputs and build the cmdline

  # Look for reads
  SPADES_READS=
  paired_reads_result=$(process_reads paired_reads[@] 'pe')
  SPADES_READS="$SPADES_READS $paired_reads_result"

  mate_pairs_result=$(process_reads mate_pairs[@] 'mp')
  SPADES_READS="$SPADES_READS $mate_pairs_result"

  unpaired_reads_result=$(process_reads unpaired_reads[@] 'single')
  SPADES_READS="$SPADES_READS $unpaired_reads_result"

  # Check various options
  SPADES_MODE="--disable-gzip-output"
  if $is_single_cell ; then
      SPADES_MODE="$SPADES_MODE --sc"
  fi
  if $is_only_assembler ; then
      SPADES_MODE="$SPADES_MODE --only-assembler"
  fi
  if $is_careful_mode ; then
      SPADES_MODE="$SPADES_MODE --careful"
  fi

  if [ -n "$k" ]; then
      SPADES_MODE="$SPADES_MODE -k $k"
  fi

  if [ -n "$memory" ]; then
      SPADES_MODE="$SPADES_MODE -m $memory"
  fi

  echo "Desired cmdline: $SPADES_READS $SPADES_MODE"
  spades.py -o spades_dnanexus $SPADES_READS $SPADES_MODE

  contigs=$(dx upload spades_dnanexus/contigs.fasta --brief)
  if [ -f spades_dnanexus/scaffolds.fasta ]; then
      scaffolds=$(dx upload spades_dnanexus/scaffolds.fasta --brief)
  else
      echo "Scaffolds not found. Outputting contigs as scaffolds."
      scaffolds=$contigs
  fi  

  dx-jobutil-add-output contigs "$contigs" --class=file
  dx-jobutil-add-output scaffolds "$scaffolds" --class=file
}
