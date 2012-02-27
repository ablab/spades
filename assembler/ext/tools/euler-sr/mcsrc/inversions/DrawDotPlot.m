function DrawDotPlot(coordsfile)

cfh = fopen(coordsfile, 'r');

while(feof(cfh) == 0)
  dotFile = fscanf(cfh, '%s', [1,1]);
  seqalen = fscanf(cfh, '%d', [1,1]);
  seqblen = fscanf(cfh, '%d', [1,1]);
  imgFile = fscanf(cfh, '%s', [1,1]);
  imgTitle = ReadString(cfh);
  fprintf('got values\n');
  dotFile
  seqalen
  seqblen
  imgFile
  imgTitle
  load(dotFile);
  clf;
  figh = figure(1);
  f1 = gcf;
  if (size(human,1) > 10000) 
    scatter(human(:,1), human(:,3), 1, 'k');
    axis ij;
  else 
    line([human(:,1) human(:,2)]', [human(:,3) human(:,4)]', 'Color', 'k');
    axis([0 seqalen 0 seqblen]);
    axis image;
  end
  title(imgTitle);
  saveas(figh, imgFile, 'jpg');
  clear human;
end
