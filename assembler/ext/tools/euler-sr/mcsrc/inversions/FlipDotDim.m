function m2 = FlipDotDim(m,l) 
m(:,2) = l - m(:,2);
m(:,4) = l - m(:,4);
m2 = m;
