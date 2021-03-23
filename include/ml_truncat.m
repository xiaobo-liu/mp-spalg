function f = ml_truncat(z,alpha,beta,d)
% A summation algorithm that computes using 200 significant figures the 
% scalar Mittag-Leffler function by trucating the series from definition,
% and then rounds the result to d<=200 significant figures.

d_comput = 200; % digits used in the summation

d_old = mp.Digits(); % 'mp digits' at the start
mp.Digits(d_comput);
alpha = mp(alpha);
beta = mp(beta);
gk = @(k) gamma(mp(alpha*k + beta));
if ~isvector(z)
    z = mp(z);
    j = 0;
    p = mp('1');
    f = mp('0'); 
    tau = p/gamma(beta); % the first term in the series
    while abs(tau)>eps('mp')
        f = f + tau; % f sums each term until 'convergence'
        j = j + 1;
        p = p.*z;
        tau = p/gk(j);
    end
else
    n = length(z);
    f = zeros(n,1,'mp');
    for i=1:n
        z(i) = mp(z(i));
        j = 0;
        p = mp('1');
        tau = p/gamma(beta);
        while abs(tau)>eps('mp')
            f(i) = f(i) + tau;
            j = j + 1;
            p = p.*z(i);
            tau = p/gk(j);
        end
    end
end
mp.Digits(d);
f = mp(f);
mp.Digits(d_old); % return 'mp digits' to the value it had at the start
end