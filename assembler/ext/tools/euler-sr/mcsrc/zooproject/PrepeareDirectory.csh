#!/usr/bin/env csh

# $1 = target file

# configure directory structure
mkdir -p lav
mkdir -p net
mkdir -p chain
mkdir -p netlav
mkdir -p fasize
mkdir -p axt 
mkdir -p nib
mkdir -p blastdb
mkdir -p inversions
mkdir -p loci
foreach spec (`cat $1`)
 if (-e $spec.fasta) then
   faSize -detailed=on $spec.fasta > fasize/$spec.size
   faToNib $spec.fasta nib/$spec.nib
   formatdb -p F -i $spec.fasta ; mv $spec.fasta.{nin,nsq,nhr} blastdb; ln $spec.fasta blastdb/
  endif
end 
