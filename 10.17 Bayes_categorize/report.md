# 10.17第四次上机课

### 1. simple bayes

是

p1=1.6875, p2=0.0625, p1 > p2

code:

```

% 大小 颜色 形状 是否好果子
% 大1 红1 圆1 是1
data = [
    0,0,0,0;
    1,1,0,1;
    1,1,1,1;
    1,0,1,1;
    1,0,0,0;
    0,1,1,1;
    0,1,0,0;
    0,0,1,0;
    ];

% m个数据，n个维度
% p_vecx_1(1) = p(x=1)
% p_vecx_1(4) = p(C=1)
[m,n] = size(data);
for i = 1:n-1
    p_vecx_1(i) = sum(data(:,i))/m;
end
p_C_1 = sum(data(:,n))/m;

% here, temp_data only consider Ck=1
temp_data = data;
for i = 1:m
    if temp_data(i,n) == 0
        temp_data(i,1:n-1) = 0
    end
end
temp_m = sum(temp_data(:,n));
for i = 1:n-1
    p_vecx_1_C1(i) = sum(temp_data(:,i))/temp_m;
end

p_vecx_1_C0 = (p_vecx_1-p_C_1*p_vecx_1_C1)/(1-p_C_1);

input = [1,1,1];
p1 = p_vecx_1_C1(1)*p_vecx_1_C1(2)*p_vecx_1_C1(3)*p_C_1/(p_vecx_1(1)*p_vecx_1(2)*p_vecx_1(3));
p2 = p_vecx_1_C0(1)*p_vecx_1_C0(2)*p_vecx_1_C0(3)*(1-p_C_1)/(p_vecx_1(1)*p_vecx_1(2)*p_vecx_1(3));

```

### 2. Gauss bayes

男

ln_p0 = -25.971197781768100
ln_p1 = -13.131434139811583

code:

```

% 身高，体重，脚长，性别
% 男1
data = [
    6, 180, 12, 1;
    5.92,190,11,1;
    5.58,170,12,1;
    5.92,165,10,1;
    5,100,6,0;
    5.5,150,8,0;
    5.42,130,7,0;
    5.75,150,9,0;
];

[m,n] = size(data);
p_C_1 = sum(data(:,n))/m;
n = n-1;

% %Sigma_Omega
% mu = sum(data(:,1:3))/m;
% Sigma_Omega = zeros(n,n);
% for l = 1:m
%     vec_x_mu = data(l,1:3)-mu;
%     Sigma_Omega = Sigma_Omega + vec_x_mu'*vec_x_mu;
% end

%Sigma_1
temp_data = data(1:4,:);
[temp_m,~] = size(temp_data);
mu_1 = sum(temp_data(:,1:3))/temp_m;
Sigma_1 = zeros(n,n);
for l = 1:temp_m
    vec_x_mu = temp_data(l,1:3)-mu_1;
    Sigma_1 = Sigma_1 + vec_x_mu'*vec_x_mu;
end
%Sigma_0
temp_data = data(5:8,:);
[temp_m,~] = size(temp_data);
mu_0 = sum(temp_data(:,1:3))/temp_m;
Sigma_0 = zeros(n,n);
for l = 1:temp_m
    vec_x_mu = temp_data(l,1:3)-mu_0;
    Sigma_0 = Sigma_0 + vec_x_mu'*vec_x_mu;
end
%calculate
Sigma_0 = Sigma_0/3;
Sigma_1 = Sigma_1/3;
x = [6,130,8];
vec_x_mu = x-mu_1;
ln_p1 = -1/2*(log(det(Sigma_1))+vec_x_mu*inv(Sigma_1)*vec_x_mu');
%+log(p_C_1);
%ln_p1 = ln_p1 + 1/2*(log(abs(det(Sigma_Omega)))+vec_x_mu*Sigma_Omega*vec_x_mu'/7);
vec_x_mu = x-mu_0;
ln_p0 = -1/2*(log(det(Sigma_0))+vec_x_mu*Sigma_0^(-1)*vec_x_mu');
%+log(1-p_C_1);
%ln_p0 = ln_p0 + 1/2*(log(abs(det(Sigma_Omega)))+vec_x_mu*Sigma_Omega*vec_x_mu'/7);
```
