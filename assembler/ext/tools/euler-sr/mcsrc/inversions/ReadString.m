function str = ReadString(fh)
c = 'a';
fprintf('starting, %c', c);
while (strcmp(c,'#') ==0)
  c = fscanf(fh, '%c', [1,1]);
  res = strcmp(c,'#');
end
c = fscanf(fh, '%c', [1,1]);
str = '';
spos = 1
while (feof(fh) == 0 && strcmp(c, '#') == 0)
  str(spos) = c;
  spos = spos + 1;
  c = fscanf(fh, '%c', [1,1]);
end
str
