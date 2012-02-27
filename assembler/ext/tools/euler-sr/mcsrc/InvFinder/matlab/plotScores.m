function plotScores(len, pos, scores)
figure(1);
clf;
hold on;
axis([1,len,min(scores)-100, max(scores) + 100]); 
plot(scores(1:pos(1)), 'k', 'MarkerSize', 1);
plot(pos(1):pos(2)-1, scores(pos(1):pos(2)-1), 'b', ...
     'MarkerSize', 1);
plot(pos(2):pos(3)-1, scores(pos(2):pos(3)-1), 'r', ...
     'MarkerSize', 1);
plot(pos(3):pos(4)-1, scores(pos(3):pos(4)-1), 'b', ...
     'MarkerSize', 1);
plot(pos(4):len, scores(pos(4):len), 'k', ...
     'MarkerSize', 1);
hold off;
figure(2);
axis([1,len,min(scores)-100, max(scores) + 100]); 
x1 = 1:len;
xx = 1:0.1:len;
yy = spline(x1,scores, xx);
plot(xx,yy);

%plot(scores(1:pos(1)), 'k.', 'MarkerSize', 1);
%plot(pos(1):pos(2)-1, scores(pos(1):pos(2)-1), 'b.', ...
%     'MarkerSize', 1);
%plot(pos(2):pos(3)-1, scores(pos(2):pos(3)-1), 'r.', ...
%     'MarkerSize', 1);
%plot(pos(3):pos(4)-1, scores(pos(3):pos(4)-1), 'b.', ...
%     'MarkerSize', 1);
%plot(pos(4):len, scores(pos(4):len), 'k.', ...
%     'MarkerSize', 1);
hold off;



