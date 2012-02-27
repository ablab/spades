#!/usr/bin/env csh

#jif ($argv < 2) then
 # echo usage: CombinePdfs.csh outputfile file1.ps [file2.ps...]
#else
  set outfile = $argv[1]
  shift
  gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=$outfile -dBATCH $argv 
#endif
