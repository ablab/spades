#!/usr/bin/env csh
if ($# < 2) then
 echo "usage: CleanupDirectory.csh targetfile cleaning_commmands"
 exit
endif

set idfile  = $1
set outfile = $2
set nameFile = "~/projects/mcsrc/zooproject/zoonames.txt"
set cwd = `pwd`

foreach id1 (`cat $idfile`)
 foreach id2 (`cat $idfile`)
   @ failed = 0
   if (! -e lav/$id1.$id2.lav ||  -z lav/$id1.$id2.lav || ! -e axt/$id1.$id2.axt || -z axt/$id1.$id2.axt || ! -e chain/$id1.$id2.chain || -z chain/$id1.$id2.chain  ) then
     @ failed = 1
   else
    tail -1 lav/$id1.$id2.lav | grep -q "eof"
      if ($status != 0) then
        @ failed =1
      endif
   endif
   if ($failed == 1) then
     echo "cd $cwd; ~/projects/mcsrc/zooproject/ProcessAlignment.csh $id1 $id2 " >> $outfile
   endif
  end
end

