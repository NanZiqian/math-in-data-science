
clear
clc

k=3;knn=20;
%% original sample
load CP.mat
figure(1)
subplot(1,2,1)
gscatter(fea(:,1),fea(:,2),gnd)

load SP.mat
subplot(1,2,2)
gscatter(fea(:,1),fea(:,2),gnd)

%% CP.mat
load CP.mat
[m,n] = size(fea);
% sequence distance
D = pdist2(fea,fea);
[knn_dist,knn_index] = mink(D,8);
% Phi
Phi = zeros(m,m);
for i = 1:m
    Phi(i,:) = sum(D<D(:,i),2);
end
Phi = max(Phi,Phi');

sign = Phi<knn;
[~,id] = sort(Phi);% search by column
for i = 1:m
    phi_i_ik(i) = Phi(i,id(k+1,i));
end
A = exp(-Phi.^2./phi_i_ik);
A = A.*sign;
A = min(A,A');

D = sqrt(sum(A,2)).^-1;
% D has been sqrt-ed
D = diag(D);
L = eye(m)-D*A*D;
L = sparse(L);
[Z,~] = eigs(L,2,'smallestabs');
% A = sparse(A);
% L = bsxfun(@ldivide,D,A);
% L = speye(m)-bsxfun(@rdivide,L,D');
% [Z,~] = eigs(L,2,'smallestabs');
Y = D*Z;

figure(2)
gscatter(Y(:,1),Y(:,2),gnd)

%% SP.mat
load CP.mat
[m,n] = size(fea);
% sequence distance
D = pdist2(fea,fea);
[knn_dist,knn_index] = mink(D,8);
% Phi
Phi = zeros(m,m);
for i = 1:m
    Phi(i,:) = sum(D<D(:,i),2);
end
Phi = max(Phi,Phi');

sign = Phi<knn;
[~,id] = sort(Phi);% search by column
for i = 1:m
    phi_i_ik(i) = Phi(i,id(k+1,i));
end
A = exp(-Phi.^2./phi_i_ik);
A = A.*sign;
A = min(A,A');

D = sqrt(sum(A,2)).^-1;
% D has been sqrt-ed
D = diag(D);
L = eye(m)-D*A*D;
L = sparse(L);
[Z,~] = eigs(L,3,'smallestabs');
Y = D*Z;

idx = kmeans(Y,3);
figure(3)
scatter3(Y(:,1),Y(:,2),Y(:,3),80,idx,'filled')

