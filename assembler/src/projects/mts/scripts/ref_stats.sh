#/bin/bash

if [ "$#" -lt 3 ]; then
    echo "Usage: identify.sh <assemblies_folder> <refs_folder> <out_dir>"
    exit 1
fi

CTG_LENGTH_THR=5000
process_cnt=4
thread_cnt=8
assemblies_folder=$1
refs_folder=$2
#canopy_out=$3
out_dir=$3

folder=$out_dir/metaquast

export LC_ALL=C
mkdir -p $out_dir

~/git/quast/metaquast.py --debug -R $refs_folder -o $out_dir/metaquast $assemblies_folder/*.fasta

#awk ' {print $2,$1} ' $canopy_out | sort > $folder/clusters.txt

rm -rf $out_dir/ref_summary.txt

for ref in $refs_folder/*.fasta ; do
    echo "Processing reference $ref" 
    ref_name=$(basename "$ref")
    ref_name="${ref_name%.*}"

    rm -rf $out_dir/${ref_name}.ctgs

    #for sample in $assemblies_out_dir/sample9.fasta ; do
    for sample in $assemblies_folder/*.fasta ; do 
        sample_name=$(basename "$sample")
        sample_name="${sample_name%.*}"
        aligned=$out_dir/metaquast/quast_corrected_input/${sample_name}_to_${ref_name}.fasta
        ~/git/ngs_scripts/contig_length_filter.py $CTG_LENGTH_THR $aligned $out_dir/long.fasta.tmp
        ~/git/ngs_scripts/contig_info.py $out_dir/long.fasta.tmp $out_dir/ctg.info.tmp
        sed_command="s/ID_/${sample_name}-/g"
        grep -Eo "ID_.*$" $out_dir/ctg.info.tmp | sed -e $sed_command >> $out_dir/${ref_name}.ctgs
        rm $out_dir/long.fasta.tmp
        rm $out_dir/ctg.info.tmp
    done

    sed 's/$/ '"${ref_name}"'/g' $out_dir/${ref_name}.ctgs >> $out_dir/ref_summary.txt

    #sort $out_dir/${ref_name}.ctgs.tmp > $out_dir/${ref_name}.ctgs

    #join $out_dir/${ref_name}.ctgs $out_dir/clusters.txt | awk ' { print $2 } ' | sort | uniq -c | sort -nr | head -10 

    #join $out_dir/${ref_name}.ctgs $out_dir/clusters.txt > $out_dir/join.txt 
    #awk ' { print $2 } ' $out_dir/join.txt | sort | uniq -c | sort -nr | head -10 

    report=$out_dir/metaquast/runs_per_reference/$ref_name/report.txt

    grep "Assembly" $report
    grep "Genome fraction" $report
done

#rm -rf $out_dir
echo "Finished"
