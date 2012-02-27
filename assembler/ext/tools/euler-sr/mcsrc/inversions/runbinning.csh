#!/usr/bin/env csh

pushd /scratch/mchaisso/cftr_ext 
~/fwg_scripts/launch_mysqld.csh   global 
set region = $argv[1]
set dir = $argv[2]
set ext = $argv[3]

pushd /Volumes/external/encode/oriented_july/$dir
~/projects/mcsrc/inversions/db_find_inversions_all.pl $dir $dir.$ext.inversions
popd
pushd /Volumes/external/inversions1_2/$dir
cp /Volumes/external/encode/oriented_july/$dir/$dir.$ext.inversions .
~/projects/mcsrc/inversions/projectionFreeBinning.pl $dir.$ext.inversions sequences.txt fasize $dir.$ext.bins 
~/projects/mcsrc/inversions/GetBinUnionCoordinates.pl $dir.$ext.bins fasize $dir.$ext.coords
mkdir $dir.$ext
~/projects/mcsrc/inversions/ExtractBinSequences.pl $dir.$ext.coords . $dir.$ext/$dir
pushd $dir.$ext
## first mask everything
foreach file (`ls *.fasta`)
  ~/projects/mcsrc/sequtils/powerpc/mskrep $file $file:r.msk.fasta
end

foreach file (`ls *.msk.fasta`)
  ~/projects/mcsrc/inversions/MakeBinDataFile.pl $file cftr_ext $file:r.ortho
  ~/projects/mcsrc/inversions/checkBinRevised.pl cftr_ext $file ../sequences.txt $region ../ ../blastdb $file:r.ortho $file:r.conf.fasta 0.00001
  ~/projects/mcsrc/inversions/FastaBinsToCoordinates.pl $file:r.conf.fasta $file:r.conf.coords
  ~/projects/mcsrc/inversions/ExtractBinSequences.pl $file:r.conf.coords .. -f $file:r:r.unmsk.conf.fasta
  ~/projects/mcsrc/inversions/make_grid2.pl $file:r:r.unmsk.conf.fasta ../sequences.txt $file:r:r.unmsk.chars
#  ~/projects/mcsrc/comparative/i686/makegrid $file:r:r.unmsk.conf.fasta ../sequences.txt $file:r:r.chars $file:r:r.data
end
