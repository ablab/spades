set term png size 800,600
FILE_NAME_T='time_insert.tmp'
FILE_NAME_F='time_find.tmp'
FILE_NAME_M='memory.tmp'
set style fill solid
set ylabel 'sec'
set output 'time_insert.png'
plot [-0.5:2.5] [0:] FILE_NAME_T index 0 using ($2):xticlabels(3) title MAP1 with histograms, FILE_NAME_T index 1 using ($2):xticlabels(3) title MAP2 with histograms, FILE_NAME_T index 2 using ($2):xticlabels(3) title MAP3 with histograms, FILE_NAME_T index 3 using ($2):xticlabels(3) title MAP4 with histograms, FILE_NAME_T index 4 using ($2):xticlabels(3) title MAP5 with histograms
set output 'time_find.png'
plot [-0.5:2.5] [0:] FILE_NAME_F index 0 using ($2):xticlabels(3) title MAP1 with histograms, FILE_NAME_F index 1 using ($2):xticlabels(3) title MAP2 with histograms, FILE_NAME_F index 2 using ($2):xticlabels(3) title MAP3 with histograms, FILE_NAME_F index 3 using ($2):xticlabels(3) title MAP4 with histograms, FILE_NAME_F index 4 using ($2):xticlabels(3) title MAP5 with histograms
set ylabel 'KB'
set output 'memory.png'
plot [-0.5:2.5] [0:] FILE_NAME_M index 0 using ($2):xticlabels(3) title MAP1 with histograms, FILE_NAME_M index 1 using ($2):xticlabels(3) title MAP2 with histograms, FILE_NAME_M index 2 using ($2):xticlabels(3) title MAP3 with histograms, FILE_NAME_M index 3 using ($2):xticlabels(3) title MAP4 with histograms, FILE_NAME_M index 4 using ($2):xticlabels(3) title MAP5 with histograms

