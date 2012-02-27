function DrawClusters(coords, clus, w) 
  
for i = 1:size(clus,1)
  j = 1;
  color = random('uniform', 0,1,1,3); 
  ind = find(clus(i,:) ~= -1); 
  nc = size(ind,2);
  if (nc == 1) color = [0,0,0]; end;
  if (nc == 2) color = [0,1,0]; end;
  if (nc == 3) color = [0,0,1]; end;
  if (nc >= 4) color = [1,0,0]; end;
  if (size(ind,2) > 1)
    for ji = 1:size(ind,2)
      cl = clus(i,ind(ji));
      x = coords(cl,1);
      y = coords(cl,3);
      
      line([x-w, x+w], [y-w,y-w],'Color', color);
      line([x+w, x+w], [y-w,y+w],'Color', color);
      line([x-w, x+w], [y+w,y+w],'Color', color);
      line([x-w, x-w], [y+w,y-w],'Color', color);
    end
  end
end
