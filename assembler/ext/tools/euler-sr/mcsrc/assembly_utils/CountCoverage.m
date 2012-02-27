function P = CountCoverage(cov,win) 
	n = size(cov,1);
	mcov = max(cov);
	L = mcov(1);
	% Compute the coverage per position
	P1 = zeros(L+1,1);
	for i=1:n
		P1(cov(i)+1) = P1(cov(i)+1) + 1;
	end

	P = zeros(L+1-win+1,1);
	for i = 1:L+1-win
		P(i) = sum(P1(i:i+win));
		P(i) = P(i) / 50;
	end
end
		
	
	
	
