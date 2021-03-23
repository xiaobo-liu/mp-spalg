%TEST   Simple test of FUNM_ND.

format compact
addpath('include');

n = 50;
num = 10;
rng(1)

fprintf('Normwise relative differences between FUNM_ND for the matrix exponential \n')
fprintf('and EXPM, and between FUNM_ND for the matrix square root and SQRTM.\n')
fprintf('The matrices are random matrices of dimension %g\n',n)
fprintf('and the differences should be of order some reasonable multiple of %g.\n',eps/2)
for i = 1:num
    A = rand(n);
    funmA_exp = funm_nd(A,@exp); 
    expmA = expm(A);
    funmA_sqrt = funm_nd(A,@sqrt); 
    sqrtmA = sqrtm(A);
    fprintf('%2.0f: exp%8.1e, ', i, norm(funmA_exp-expmA,'fro')/norm(expmA,'fro'))
    fprintf('sqrt%8.1e\n', norm(funmA_sqrt-sqrtmA,'fro')/norm(sqrtmA,'fro'))
end
