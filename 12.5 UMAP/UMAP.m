
clear
clc
% 20 items, each has 72 pictures with diff angles
load COIL20.mat

%fea \in Rmn, Y \in Rmk
m = 1440;n = 1024;k = 2;
fea_max_sqaured = norm(fea,'inf')^2;
shrink_rate = 1;
%number of iteration in gredient descent
Iter_Time = 2000;
%K is used to calculate wij
K = 21;
%b is used in calculate vij
b = 1.5;
%epsilon is used in calculate g(Y)'s gredient to yj
epsilon = 0.001;
% UMAP: minimize the sum of similarity change before and after jiangwei

%% similarities of fea
% W = zeros(m,m);
% % W = {||xi-xj||^2 m*m}
% for j = 1:m
%     W(:,j) = sum((fea-fea(j,:)).^2,2);
% end
% 
% %idx(i,j)代表第i个样本第j近的行号,明显idx(i,1)=i
% %D(i,j)代表第i个样本和第j近的样本的距离
% [idx,D] = knnsearch(fea,fea,"K",K);
% 
% sigma = zeros(m,1);
% for i = 1:m
%     fun = @(x)sigma_func(x,i,D(i,:));
%     sigma(i) = fsolve(fun,1);
% end
% sigma = sigma.^(-2);
% W_hat = exp(-W.*sigma);
% W = W_hat + W_hat' - W_hat.*W_hat';

D = pdist2(fea,fea);
[knn_dist,knn_index] = mink(D,21);
knn_dist = knn_dist(2:end,:);
knn_index = knn_index(2:end,:);
n = size(D,1);
log2k = log2(20);
sigma = zeros(n,1);
for ii = 1:n
    f = @(x) sum(exp(-(knn_dist(:,ii)-knn_dist(1,ii))/x))-log2k;
    options = optimset('Display','off');
    sol = fsolve(f,1,options);
    sigma(ii) = sol;
end
sigma(sigma<1) = 1;

W = exp(-bsxfun(@rdivide,D.^2,sigma.^2));
W = W+W'-W.*W'; W(W<eps) = 0;
strong_graph = zeros(n,n);
for ii = 1:n
    strong_graph(ii,knn_index(:,ii)) = 1;
end
strong_graph = min(strong_graph,strong_graph');
strong_graph = strong_graph+diag(ones(n,1));
W = W.*strong_graph;

%% update yi, diff sample influence diffly to yi's gredient
% super coefficient
c1 = 20;c2 = 20;
% initialize
Y = rand(m,k);
figure_num = 1;

for iter = 1:Iter_Time
    alpha = 1 - iter/Iter_Time;
    % each iteration, update y1-ym once
    for i = 1:m
        % calculate f(Y)'s gredient for yi
        gredient = 0;
        % c1 yj in f(Y)'s gredient
        ii = 0;list_num = 10;%temp = 0;
        while ii <= c1
            % select qualified points
            pivot = rand(1,1);
            list = knn_index(1:list_num,i);
            j = get_j(i,m,list);
            
            if pivot*shrink_rate < W(i,j)
                d = norm(Y(i,:)-Y(j,:),2);
                gredient = gredient + 2*b*d^(2*b-2)/(1+d^(2*b))*(Y(i,:)-Y(j,:));
                ii = ii + 1;
            else
                break
%                 if temp > 100
%                     break;
%                 else
%                     temp = temp + 1;
%                     continue;
%                 end
            end
        end
        % calculate g(Y)'s gredient for yi
        ii = 0;%temp = 0;
        while ii <= c2
            pivot = rand(1,1);
            j = get_j(i,m,list);
            if pivot*shrink_rate < 1 - W(i,j)
                d = norm(Y(i,:)-Y(j,:),2);
                gredient = gredient - 2*b/(epsilon+d^2)/(1+d^(2*b))*(Y(i,:)-Y(j,:));
                ii = ii + 1;
            else
                break
%                 if temp > 100
%                     break;
%                 else
%                     temp = temp + 1;
%                     continue;
%                 end
            end
        end
        % apply gredient of f(Y)+g(Y) to yi, total c1 + c2
        % how to avoid gredient explotion?
        % Y(i,:) = Y(i,:) - alpha*gredient/c1/c2;
        if abs(gredient)<fea_max_sqaured
            Y(i,:) = Y(i,:) - alpha*gredient;
        else
            continue
        end
    end
    if rem(iter,100) == 99
        figure(figure_num)
        gscatter(Y(:,1),Y(:,2),gnd)
        pause(1)
        figure_num = figure_num + 1;
    end
end

figure(figure_num)
gscatter(Y(:,1),Y(:,2),gnd)

%% solve sigma
function f = sigma_func(x,i,D)
    f = 1 - log(20);
    for j = 2:20
        f = f + exp(-(D(j+1)-D(2))/x);
    end
end

%% get j,a random int from 1 to m, but != i
% W is sparse, j need to be around i
function j = get_j(i,m,list)
    [list_num,~] = size(list);
    j = list(ceil(rand(1,1)*list_num));%[i-2,i+2]
%     while j == i || j <= 0 || j > m
%         j = i-3+ceil(rand(1,1)*5);
%     end
end
