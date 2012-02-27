function DotPlot(seqA, seqB, figh, imgout, imgtitle)

[a, b, coords] = seqdotplotsilent(seqA, seqB, 8,8);

seqBrc = seqrcomplement(seqB);
[a,b,coordsrc] = seqdotplotsilent(seqA, seqBrc, 8,8);

lenb= size(seqB,2);

coordsrcf = FlipDotDim(coordsrc,lenb);
% draw the forward direction
fprintf('coord sizes, coords:\n');
size(coords)
fprintf('coordsrcf:\n');
size(coordsrcf)
DrawDot(coords);
DrawDot(coordsrcf);
axis([0 size(seqA,2) 0 size(seqB,2)])
title(imgtitle);
saveas(figh, imgout, 'jpg');
