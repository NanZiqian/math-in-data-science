
% |difference of mu, difference of sigma|
Difference = zeros(4,2);

for ii = 1:4
    Num = 10^ii;
    [Difference(ii,1),Difference(ii,2)] = func(Num);
end

x = [1,2,3,4];
figure(1);
scatter(x,Difference(:,1),'r');
xlabel('log_10(number of data)');
ylabel('均值之差');
figure(2);
scatter(x,Difference(:,2),'r');
xlabel('log_10(number of data)');
ylabel('方差之差');

% produce mu and sigma of estimatet data
function [mu,sigma] = func(Num)
    data = sqrt(5)*randn(Num,1)+1;
    [mu,sigma] = normfit(data);
    mu = mu-1;
    sigma = sigma - sqrt(5);
end