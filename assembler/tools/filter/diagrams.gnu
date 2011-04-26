set term png
FILE_NAME_T='time.tmp'
FILE_NAME_M='memory.tmp'
set output 'time.png'
plot [-0.5:4.5] [0:] FILE_NAME_T index 0 using ($2)*60+($3) title 'FILE1' with histograms, FILE_NAME_T index 1 using ($2)*60+($3) title 'FILE2' with histograms, FILE_NAME_T index 2 using ($2)*60+($3) title 'FILE3' with histograms
set output 'memory.png'
plot [-0.5:4.5] [0:] FILE_NAME_M index 0 using ($2) title 'FILE1' with histograms, FILE_NAME_M index 1 using ($2) title 'FILE2' with histograms, FILE_NAME_M index 2 using ($2) title 'FILE3' with histograms

