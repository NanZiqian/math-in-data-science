
% exercise 7.5
clear
clc

D = [0,1.2,3.4,5.1,6.2,7.7;
    1.2,0,2.7,4.8,6.3,8.2;
    3.4,2.7,0,4.4,5.5,3.3;
    5.1,4.8,4.4,0,2.4,2.5;
    6.2,6.3,5.5,2.4,0,0.9;
    7.7,8.2,3.3,2.5,0.9,0
];

[m,~] = size(D);

Cn = eye(m)-ones(m,m)/m;
M = Cn*(D.^2)*Cn/(-2);
[L,V] = eig(M);
eigen_value = diag(V);
[~,id] = sort(eigen_value,"DESCEND");
Q = L(:,id(1:2))*diag(eigen_value(id(1:2)).^(0.5))

scatter(Q(:,1),Q(:,2))
W = zeros(m,m);
% W = {||xi-xj||}
for j = 1:m
    W(:,j) = sum((Q-Q(j,:)).^2,2);
end
error = sqrt(W)-D