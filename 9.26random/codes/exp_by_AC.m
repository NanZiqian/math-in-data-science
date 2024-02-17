
% x on [0,1], y on [0,5]
x1 = linspace(0,100,1000);
f = 5*exp(-5*x1);
%1000 points
points = rand(1000,2);
% points(:,1) = points(:,1);
points(:,2) = points(:,2)*5;

for ii = 1:1000

    if points(ii,2) > 5*exp(-5*points(ii,1))
        points(ii,1) = -1;
    end
end
%
freq = zeros(1,10);
for ii = 1:1000 % all 1000 points
    if points(ii,1) < 0
        continue
    end
    freq(floor(points(ii,1)/0.1)+1) = freq(floor(points(ii,1)/0.1)+1) + 1;
end

subplot(1,2,1);
plot(0:0.1:0.9,freq);
axis([0,5,0,80]);

subplot(1,2,2);
plot(x1,f);
axis([0,5,0,0.1]);

