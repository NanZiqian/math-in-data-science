
A1 = [4,-1,0;
    -1,4,-1;
    0,-1,4];

A2 = diag(-ones(3,1));

A3 = zeros(3,3);

A = [A1,A2,A3;
    A2,A1,A2;
    A3,A2,A1];

r2 = sqrt(2);

b = [1,r2,1,r2,2,r2,1,r2,1]'/16;

u = A^-1*b;

u_hat = zeros(3,3);

for i = 1:3
    for j = 1:3
        u_hat(i,j) = sin(pi*i/4)*sin(pi*j/4)/pi^2;
    end
end

u_hat2 = zeros(9,1);
for i = 1:3
    u_hat2(3*(i-1)+1:3*i,1) = u_hat(i,:)';
end

max_n = norm(u_hat2-u,inf)
max_2 = norm(u_hat2-u,2)