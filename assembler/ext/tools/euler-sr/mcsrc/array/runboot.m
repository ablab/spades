load clone_array.mat;


expdef = [1 8; 9 16; 17 24];
addpath('/home/mchaisso/projects/mcsrc/array');

[chipdist, chipsup, chipmean, chipstdev, nsamp] = BootstrapDistChip(ca, expdef);


save bootstrap_clone_array.mat chipdist chipsup chipmean chipstdev nsamp;

quit;
