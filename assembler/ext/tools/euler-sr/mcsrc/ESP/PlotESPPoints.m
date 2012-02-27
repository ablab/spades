function PlotESPPoints(p)
  ffi = find(p(:,2) == 1 & p(:,4) == 1);
  fri = find(p(:,2) == 1 & p(:,4) == -1);
  rfi = find(p(:,2) == -1 & p(:,4) == 1);
  rri = find(p(:,2) == -1 & p(:,4) == -1);

 hold on;
  plot(p(ffi,1), p(ffi,3), 'k+');
  plot(p(fri,1), p(fri,3), 'k>');
  plot(p(rfi,1), p(rfi,3), 'k<');
  plot(p(rri,1), p(rri,3), 'kx');
hold off;
