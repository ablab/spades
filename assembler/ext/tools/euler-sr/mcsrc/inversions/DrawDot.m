function DrawDot(d, color)
if (nargin < 2)
  color = 'k';
end
  
line([d(:,1) d(:,2)]', [d(:,3), d(:,4)]', 'Color', color, 'LineWidth', ...
     1.0);
axis([0 max(d(:,2)) 0  max(d(:,4))]);
