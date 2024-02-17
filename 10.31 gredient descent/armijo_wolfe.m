
%% main
% result:
% argmin = [1.0000]4*3
% iter_num = [2911	24009	59140
%            2945	20563	50061]
f = @(x) 100*(x(1)^2-x(2))^2+(x(1)-1)^2;
gradf = @(x) [400*x(1)*(x(1)^2-x(2))+2*x(1)-2;-200*(x(1)^2-x(2))];
data_x0 = [[0;0],[1;10],[20;20]];
% armijo [;][;][;];
% wolfe  [;][;][;] 
argmin = zeros(4:3);
% armijo [fisrt,second,third];
% wolfe  [fisrt,second,third]
iter_num = zeros(2:3);
for i = 1:3
    % gredient descent by armijo
    xk = data_x0(:,i);
    iter = 0;
    while norm(gradf(xk),2) >= 1e-5
        alpha = armijo(xk,f,gradf);
        xk = xk-alpha*gradf(xk);
        iter = iter + 1;
    end
    argmin(1:2,i) = xk;
    iter_num(1,i) = iter;

    % wolfe
    xk = data_x0(:,i);
    iter = 0;
    while norm(gradf(xk),2) >= 1e-5
        alpha = wolfe(xk,f,gradf);
        xk = xk-alpha*gradf(xk);
        iter = iter + 1;
    end
    argmin(3:4,i) = xk;
    iter_num(2,i) = iter;
end



% armijo
% xk is a vector
function [alpha] = armijo(xk,f,gradf)
    m=0;
    beta = 0.5;
    sigma = 0.2;
    gk = gradf(xk);
    dk = -gk;
    while 1
        alpha = beta^m;
        if f(xk+alpha*dk) <= f(xk) + sigma*alpha*gk'*dk
            break;
        end
        m = m+1;
    end
end

% wolfe
function [alpha] = wolfe(xk,f,gradf)
    beta1 = 0.5;
    beta2 = 1/3;
    sigma = 0.3;
    rho = 0.2;
    alpha = 1;
    c1 = 0;
    c2 = 0;
    while 1
        gk = gradf(xk);
        dk = -gk;
        %c1 == 是否满足条件1
        c1 = ( f(xk+alpha*dk) <= f(xk)+ rho*alpha*gk'*dk);
        if c1 == 0
            % 不满足1
            alpha = alpha*beta1;
            continue;
        else
            %满足1
            c2 =( gradf(xk+alpha*dk)'*dk >= sigma*gk'*dk );
            if c2 == 0
                %不满足2
                alpha = alpha/beta2;
                continue;
            else
                %满足2
                break;
            end
        end
        %acquire alpha
    end
end


