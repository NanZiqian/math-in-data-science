
%% generate a linear equation
m = 100;
n = 50;
A = rand(m,n);
[U,S,V] = svd(A,'econ');
temp = diag(S);
temp = temp-min(temp)+1e-16;
S = diag(temp);
A = U*S*V';
x = rand(n,1);
b = A*x;

h = 0.001;

%% solve x by ridge regression
x_lambda = zeros(n,2/h);
i=1;
for lambda = h:h:2
    x_lambda(:,i) = (A'*A+lambda*eye(n))^(-1)*(A'*b);
    i = i + 1;
end
for i = 1:2/h
    resitual(i) = norm(A*x_lambda(:,i)-b,2);
end
[~,number] = min(resitual);
error_lambda = norm(x - x_lambda(:,i),2)
lambda = h*number
x_nolambda = (A'*A)^-1*(A'*b);
error_nolambda = norm(x - x_nolambda,2)
