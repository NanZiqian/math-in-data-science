
A = [1,3;
    -2,1;
    10,6;
    6,-3;
    -1,1
    ];

sigma = [2,4,5];

mu = [0,0;
    3,4;
    -3,2];

result = zeros(5,3);

for i = 1:5
    for j = 1:3
        result(i,j) = g(A(i,:),mu(j,:),sigma(j));
    end
end

W = result./sum(result,2)

function [result] = g(x,mu,sigma)
    result = 1/(2*pi*sigma)*exp(-sigma/2*norm(x-mu,2)^2);
end