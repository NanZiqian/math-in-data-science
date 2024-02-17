
Num = 100;
data = randn(2*Num,1);
% mu = 1, sigma = 1
data(1:Num,1) = data(1:Num,1)+1;
% mu = 3, sigma = 1
data(Num:2*Num,1) = data(Num:2*Num,1)+3;

% init, randomly divide Z into two groups
seq = randperm(200);
Xn_1 = zeros(200,1);
Yn_1 = zeros(200,1);
for i = 1:200
    if i <= 100
        Xn_1(i) = data(seq(i));
    else
        Yn_1(i-100) = data(seq(i));
    end
end
[muXn_1,sigmaXn_1] = normfit(Xn_1);
[muYn_1,sigmaYn_1] = normfit(Yn_1);

for j = 1:1000
    [muXn,sigmaXn,muYn,sigmaYn] = iterate(data,muXn_1,muYn_1);
    if abs(muXn-muXn_1)<1e-4 && abs(muYn-muYn_1)<1e-4 && abs(sigmaXn-sigmaXn_1)<1e-4 && abs(sigmaYn-sigmaYn_1)<1e-4
        break;
    end
    muXn_1 = muXn;
    sigmaXn_1 = sigmaXn;
    muYn_1 = muYn;
    sigmaYn_1 = sigmaYn;
end

sprintf('muX: %f\n',muXn)
sprintf('sigmaX: %f\n',sigmaXn)
sprintf('muY: %f\n',muYn)
sprintf('sigmaY: %f\n',sigmaYn)



function [muXn,sigmaXn,muYn,sigmaYn] = iterate(data,muXn_1,muYn_1)
    xj = 1;
    yj = 1;
    for i = 1:200
        d_mu1 = abs(data(i)-muXn_1);% 小则X
        d_mu2 = abs(data(i)-muYn_1);
        if d_mu1 < d_mu2
            Xn(xj) = data(i);
            xj = xj+1;
        elseif d_mu1 > d_mu2
            Yn(yj) = data(i);
            yj = yj+1;
        else
            if rand()>0.5
                Xn(xj) = data(i);
                xj = xj+1;
            else
                Yn(yj) = data(i);
                yj = yj+1;
            end
        end
    end
    [muXn,sigmaXn] = normfit(Xn);
    [muYn,sigmaYn] = normfit(Yn);
end