
% p173 exercise 7.1
clear
clc

A = readmatrix("A.csv");
[m,n] = size(A);
[U,Sigma,V] = svd(A);
% third right singular vector
disp("third right singular vector: ")
V(:,3)
% second singular value
disp("second singular value: ")
Sigma(2,2)
% fourth left singular vector
disp("fourth left singular vector: ")
U(:,4)
% rank of A
disp("rank of A: "+rank(Sigma))

% compute Ak
k=2;
Ak = zeros(m,n);
for i = 1:k
    Ak = Ak + U(:,i)*Sigma(i,i)*V(:,i)';
end
% F-norm^2 of A-Ak
F_norm = 0;
for i = 1:n
    F_norm = F_norm + norm(A(:,i)-Ak(:,i),2)^2;
end
F_norm = F_norm
% 2-norm^2 of A-Ak
two_norm = norm(A-Ak,2)^2

% center A, run PCA
A_centered = A - mean(A,2);
[U_centered,Sigma_centered,V_centered] = svd(A_centered);
B = V_centered(:,1:k);
% report
pi_B=A_centered*B*B';
F_norm_PCA = 0;
for i = 1:n
    F_norm_PCA = F_norm_PCA + norm(A_centered(:,i)-pi_B(:,i),2)^2;
end
F_norm_PCA = F_norm_PCA
two_norm_PCA = norm(A_centered-pi_B,2)^2
