#!/usr/bin/env csh
if ($# < 2) then
  echo "usage: ProcessAlignment.csh index1 index2"
  exit
endif
RunBlastz.csh $1.fasta $2.fasta lav/$1.$2.lav "M=2"
RunCmd.csh "lavToAxt lav/$1.$2.lav nib nib" "axt/$1.$2.axt"
RunCmd.csh "axtChain axt/$1.$2.axt  nib nib" "chain/$1.$2.chain"
RunCmd.csh "chainNet chain/$1.$2.chain fasize/$1.size fasize/$2.size net/$1.$2.tnet" "net/$1.$2.qnet"
RunCmd.csh "~/projects/mcsrc/net/${MACHTYPE}/netToLav chain/$1.$2.chain net/$1.$2.tnet $1 fasize/$1.size $2 fasize/$2.size " "netlav/$1.$2.net.lav"
