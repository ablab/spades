function nameMap=CreateNameMap(clone_array, nameValues)
  nameMap = zeros(1,size(nameValues,1));
  for i=1:size(nameValues,1)
    locs = find(clone_array(:,1) == nameValues(i));
    if (size(locs,1) == 1)
      nameMap(nameValues(i)) = find(clone_array(:,1)==nameValues(i));
    end
  end
