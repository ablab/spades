function s = CumulativeContigs(lengths, minLength)
	% lengths are the contig lengths
	% ignore all contigs less than minLength
	lengths = flipud(sort(lengths));
  li = find(lengths > minLength);
	nli = size(li,1);
	for i = 1:nli
		s(i) = sum(lengths(1:i));
  end
