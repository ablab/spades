function [dist, avg, pcnt, varcnt, expcnt] = CellDistribution(data, step, diagonal)
xmax = max(data(:,1));
xmin = min(data(:,1));
ymax = max(data(:,3));
ymin = min(data(:,3));

xend = floor((xmax )/ step);
yend = floor((ymax - (step + diagonal)) / step);

counts = zeros(xend*yend,1);
i = 1;
for xi = 0:xend-1
  for yi = 0:yend-1
     xs = xi * step ;
     ys = yi * step + step + diagonal;
     np = find(data(:,1) >= xs & ...
 	       data(:,1) <= xs + step & ...
	       data(:,3) >= ys & ...
	       data(:,3) <= ys + step);

     counts(i) = size(np,1);
     i = i + 1;
  end
end
maxc = max(counts);

dist = hist(counts, maxc);

avg  = mean(counts);

npoints = sum(counts);
ncells  = size(counts,1);
maxn = max(counts);
pcnt   = poisspdf([1:maxn], avg);
expcnt = pcnt * ncells;
varcnt = pcnt.*(1-pcnt)*ncells;

