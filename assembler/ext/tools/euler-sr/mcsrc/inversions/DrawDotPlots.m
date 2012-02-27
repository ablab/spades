function DrawDotPlots(coordsfile)

cfh = fopen(coordsfile, 'r');

while(feof(cfh) == 0)
  hFile = fscanf(cfh, '%s', [1,1]);
  cFile = fscanf(cfh, '%s', [1,1]);
  imgFile = fscanf(cfh, '%s', [1,1]);
  imgTitle = ReadString(cfh);
  fprintf('got input:\n');
  hFile
  cFile
  imgFile
  hsfh = fopen(hFile, 'r');
  hseq = fscanf(hsfh, '%s', [1,1]);
  
  csfh = fopen(cFile, 'r');
  cseq = fscanf(csfh, '%s', [1,1]);
  fprintf('sequenes:\n');
  hseq
  cseq
  
  lhseq = length(hseq)
  lcseq = length(cseq)
  if (lhseq < 20000 && lcseq < 20000) 
    figh = figure(1);
    clf;
    DotPlot(hseq, cseq, figh, imgFile, imgTitle);
  end
  clear hseq;
  clear cseq;
  fclose(hsfh);
  fclose(csfh);
end
