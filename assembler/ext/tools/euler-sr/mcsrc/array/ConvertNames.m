function nameValues=ConvertNames(names)
%  nameValues = zeros(1,size(names,1));
  for i=1:size(names,1)
    val = names{i};
    if (strcmp(val, '') == 0)
      nameValues(i) = str2num(names{i});
    end
  end
