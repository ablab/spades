main() {
  set -e -x

  # Look for the inputs and build the cmdline

  # Look for paired reads
  SPADES_READS=
  for index in ${!left_reads[*]}
  do
      lread_link="${left_reads[$index]}"
      lread_fname="left_$index.fastq"
      rread_link="${right_reads[$index]}"
      rread_fname="right_$index.fastq"

      dx download "$lread_link" -o $lread_fname
      dx download "$rread_link" -o $rread_fname
      SPADES_READS="$SPADES_READS -1 $lread_fname -2 $rread_fname"
  done

  # Look for single reads
  for index in ${!single_reads[*]}
  do
      uread_link="${single_reads[$index]}"
      uread_fname="unpaired_$index.fastq"

      dx download "$uread_link" -o $uread_fname
      SPADES_READS="$SPADES_READS -s $uread_fname"
  done

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
  scaffolds=$(dx upload spades_dnanexus/scaffolds.fasta --brief)

  dx-jobutil-add-output contigs "$contigs" --class=file
  dx-jobutil-add-output scaffolds "$scaffolds" --class=file
}
