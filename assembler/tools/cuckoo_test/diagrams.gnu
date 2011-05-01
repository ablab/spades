set term png
FILE_NAME_T='time_insert.tmp'
FILE_NAME_M='memory.tmp'
set output 'time_insert.png'
plot [0:] [0:] FILE_NAME_T using ($2) title 'FILE' with lines
set output 'memory.png'
plot [0:] [0:] FILE_NAME_M using ($2) title 'FILE' with lines
