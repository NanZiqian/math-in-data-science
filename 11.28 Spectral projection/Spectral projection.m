
clear
clc
load mammoth.mat

K = 100;k=3;
N=10000;

%idx(i,j)代表第i个样本第j近的行号,明显idx(i,1)=i
%D(i,j)代表第i个样本和第j近的样本的距离
[idx,D] = knnsearch(fea,fea,"K",K);
A = zeros(N,N);

%xi属于xj的K领域->i in idx(j,1:K)
for j=1:N
    A(:,j) = ismember((1:N)',idx(j,:));
end
for i = 1:N
    B(i,:) = ismember(1:N,idx(i,:));
end
A = A.*B;

%||xi-xj||^2
W = zeros(N,N);
for j = 1:N
    sigmaj = D(j,k);
    sigmai = D(:,k);
    W(:,j)=exp( -sum((fea-fea(j,:)).^2,2)./sigmai/2/sigmaj );
end
W = A.*W;

D = sqrt(sum(W,2));
% D has been sqrt-ed
W = sparse(W);
L = bsxfun(@ldivide,D,W);
L = speye(N)-bsxfun(@rdivide,L,D');
[Z,~] = eigs(L,2,'smallestabs');
Y = diag(D.^(-1))*Z;

gscatter(Y(:,1),Y(:,2),gnd);
