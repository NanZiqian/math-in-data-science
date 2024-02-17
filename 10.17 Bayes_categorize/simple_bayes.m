
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
