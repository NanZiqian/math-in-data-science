
clear
clc
%construct samples
m = 1000;n = 3;k = 2;
distinct = 2;

fea = zeros(m,n);
fea(1:distinct,:) = 10000000*rand(distinct,n);
for i = distinct+1:m
    rand_int = ceil(rand(1,1)*distinct);
    fea(i,:) = fea(rand_int,:);
end

%column normalize
for i = 1:n
    fea(:,i) = fea(:,i) - mean(fea(:,i));
end

% PCA
W = zeros(3,k);
[~,~,V] = svd(fea);
W = V(:,1:k);
Distance = zeros(m,m);
Y = fea*W*W';
% W = {||xi-xj||}
for j = 1:m
    Distance(:,j) = sum((fea-fea(j,:)).^2,2);
end

Distance_Q = zeros(m,m);
for j = 1:m
    Distance_Q(:,j) = sum((Y-Y(j,:)).^2,2);
end

Fnorm_error = norm(Distance_Q-Distance,'fro')^2;
infnorm_error = max(max(abs(Distance-Distance_Q)));
Distance-Distance_Q
