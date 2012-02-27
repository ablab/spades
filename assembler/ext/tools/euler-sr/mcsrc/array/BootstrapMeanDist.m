function [cumdist, cumsupport, mu, sigmasq, muBoot] = BoostrapMeanDist(vector)

% configure some parameters

nboot = 1000;
nbins = 100;
sigmasq = NaN;
mu      = NaN;
muBoot  = zeros(1000,1);
% Some error checking.  This entire vector could be bad (or one element)
% if so, don't calculate statistics
if (size(vector,2) <= 1)
  cumsupport = zeros(nbins,1);
  cumdist    = zeros(nbins,1);
  return;
end


% calcuate the mean of the bootstrap samples
muBoot = bootstrp(1000, @mean, vector);

% generate the pdf
[hv, cumsupport] = hist(muBoot, 100);

% integrate the histogram to get the cummulate distribution
cumdist = zeros(nbins,1);
cumdist(1) = hv(1);
for i = 2:nbins
  cumdist(i) = hv(i) + cumdist(i-1);
end

sigmasq = var(muBoot);
mu      = mean(muBoot);

