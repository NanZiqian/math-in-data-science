
% f = 5e^-5x when x>0
x1 = linspace(0,100,1000);
f = 5*exp(-5*x1);
u = rand(1000,1);
% F^-1 = -1/5*ln(1-y)
x = -0.2*log(1-u);
% [0:0.1:100] 1000 ranges
freq = zeros(1,1000);
for ii = 1:1000
    if x(ii) >= 100
        continue
    end
    freq(round(x(ii)/0.1)+1) = freq(round(x(ii)/0.1)+1) + 1;
end
subplot(1,2,1);

plot(0:0.1:99.9,freq);
axis([0,5,0,200]);
subplot(1,2,2);

plot(x1,f);
axis([0,5,0,0.1]);
