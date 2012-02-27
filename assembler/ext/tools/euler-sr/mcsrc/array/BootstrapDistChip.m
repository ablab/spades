function [chipdist, chipsupp, chipmean, chipstdev, nsamples, chipboot] = BootstrapDistChip(chipdata, expdef)
% BOOTSTRAPDISTCHIP computes the bootstrap distribution several
% experiments on a microarray.
% 
%    [clndist, chpsup, nsamp] = BOOTSTRAPDISTCHIP(chipdata, expdef) 
%       Input:
%       chipdata   The M clones by N experiments log_2 chip values.
%       expdef     Ts the definition of experiments, in the format
%                  [start1 end1; start2 end2; start3 end3]
%       Output:
%       clonedist   The cumulate distribution of mean expression value computed
%                   with the bootstrap.
%       chipsupport The values of the mean expression value.
%       nsamples    The number of samples used to compute the distribution.  
%                   Values that are NaN are not used 
  
chipdist = zeros(size(chipdata,1), size(expdef,1), 100);
chipsupp = zeros(size(chipdata,1), size(expdef,1), 100);
nsamples = zeros(size(chipdata,1), size(expdef,1));
for cloneindex = 1:size(chipdata,1)
  for expindex = 1:size(expdef,1)
    fprintf('ci: %d ei: %d\n', cloneindex, expindex);
    clonevect = SelectExperiment(chipdata, cloneindex, expindex, ...
				 expdef);
    [bootdist, bootsup, mu, sigmasq, muboot] = BootstrapMeanDist(clonevect);
    chipdist(cloneindex,expindex,1:100) = bootdist;
    chipsupp(cloneindex,expindex,1:100) = bootsup;
    nsamples(cloneindex,expindex)       = size(clonevect,2);
    chipmean(cloneindex,expindex) = mu;
    chipstdev(cloneindex, expindex) = sigmasq;
    chipboot(cloneindex, expindex, 1:1000)  = muboot;
  end
end
