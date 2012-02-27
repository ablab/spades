function reorder=ReorderArrayData(original, order, map)
  reorder = zeros(size(original));
  maxIndex = max(original(:,1));
  for i = 1:size(original,1)
    reorder(i,:) = original(map(order(i)),:);
  end

  
