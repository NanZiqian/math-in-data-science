# 第五次上机课 岭回归


```
error_lambda =

    1.0267


lambda =

   1.0000e-03


error_nolambda =

   13.3053
```

可以发现，有正则化的$x^\lambda$误差明显比没有的要小；但是不知为何，$\lambda$总是取第一个0.001，但是也情有可原，因为稍微改动就可以避免特征值0的出现，而过大的改动会让数值解偏离真值过大。

code:

```

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

```