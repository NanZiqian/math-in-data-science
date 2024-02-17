
clear
clc
load mammoth.mat

k = 2;
[m,n] = size(fea);
%column normalize
for i = 1:3
    fea(:,i) = fea(:,i) - mean(fea(:,i));
end

B = fea'*fea;

[V,D] = eig(B);

W = zeros(3,k);
for i = 1:k
    W(:,i) = V(:,end-i+1);
end

Y = fea*W;

%gscatter(Y(:,1),Y(:,2),gnd);

%%
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

error = Distance-Distance_Q;