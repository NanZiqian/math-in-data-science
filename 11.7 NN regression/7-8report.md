# 数据科学的数学基础七、八次上机

南子谦 3210104676

error: 5882.78505194171

我用的sigmoid激活函数，用relu误差更大，到了e4的数量级

code:

```

%% preparation
% 1000 samples
N = 1000;
data_z = rand(N,1);
y = data_z + sin(data_z) + 1e-6*randn(N,1);

w1 = randn(20,1);
w2 = randn(1,20);
b1 = randn(20,1);
b2 = randn(1,1);

yita = 1e-3;

%sigmoid
sigma = @(x) 1./(1+exp(-x));
dsigma = @(x) exp(-x)./(1+exp(-x)).^2;

%Relu
%sigma = @(x) max(x,0);
%dsigma = @(x) double(x>0);

%% backward
for ii = 1:1e3
    dL_dw1 = zeros(20,1);
    dL_dw2 = zeros(1,20);
    dL_db1 = zeros(20,1);
    dL_db2 = 0;
    for i = 1:N
        % sample number i
        % add dL_dw of all samples altogether
        [t1,t2,x1,x2] = forward(data_z(i),w1,w2,b1,b2,sigma);
        dL_db2 = dL_db2 + ((x2-y(i))/N).*dsigma(t2);% delta^2
        dL_dw2 = dL_dw2 + dL_db2*x1';
        dL_db1 = dL_db1 + (w2'*dL_db2).*dsigma(t1);
        dL_dw1 = dL_dw1 + dL_db1*(data_z(i)');
    end
    w1 = w1 - yita*dL_dw1;
    w2 = w2 - yita*dL_dw2;
    b1 = b1 - yita*dL_db1;
    b2 = b2 - yita*dL_db2;
end

error = 0;
for i = 1:N
    error = error + norm(y(i)-forward(data_z(i),w1,w2,b1,b2,sigma),2);
end
error
%% forward
% z is a scalar, x2 is output--a scalar
function [t1,t2,x1,x2] = forward(z,w1,w2,b1,b2,sigma)
    t1 = w1*z + b1;
    x1 = sigma(t1);
    t2 = w2*x1+b2;
    x2 = sigma(t2);
end


```
