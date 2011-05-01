set term png
FILE_NAME_T='time_insert.tmp'
FILE_NAME_F='time_find.tmp'
FILE_NAME_M='memory.tmp'
set output 'time_insert.png'
plot [-0.5:4.5] [0:] FILE_NAME_T index 0 using ($2) title 'FILE1' with histograms, FILE_NAME_T index 1 using ($2) title 'FILE2' with histograms, FILE_NAME_T index 2 using ($2) title 'FILE3' with histograms
set output 'time_find.png'
plot [-0.5:4.5] [0:] FILE_NAME_F index 0 using ($2) title 'FILE1' with histograms, FILE_NAME_F index 1 using ($2) title 'FILE2' with histograms, FILE_NAME_F index 2 using ($2) title 'FILE3' with histograms
set output 'memory.png'
plot [-0.5:4.5] [0:] FILE_NAME_M index 0 using ($2) title 'FILE1' with histograms, FILE_NAME_M index 1 using ($2) title 'FILE2' with histograms, FILE_NAME_M index 2 using ($2) title 'FILE3' with histograms

