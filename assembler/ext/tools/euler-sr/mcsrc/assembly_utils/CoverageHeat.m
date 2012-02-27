function h = CoverageHeat(cov,gc)
	covn = size(cov,1);
	gcn  = size(gc,1);

	n = min([covn,gcn]);
	h = zeros(101,101);
	maxcov = max(cov);
	maxcov = min([maxcov, 10]);
	for i=1:n
		if (cov(i) < maxcov) 
			x = floor(gc(i)*100) + 1;
			y = floor((cov(i) / maxcov)*100) + 1;
			h(x,y) = h(x,y) + 1;
		end
  end
end
