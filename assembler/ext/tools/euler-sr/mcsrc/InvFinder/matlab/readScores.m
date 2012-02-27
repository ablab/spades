function [len, pos, scores] = readScores(fh)
% len: length of scores
% pos: locations of intervals, 
% scorres: average scores;


len = fscanf(fh, '%d', 1);
pos = fscanf(fh, '%d', 4);
scores = fscanf(fh, '%f', len);

