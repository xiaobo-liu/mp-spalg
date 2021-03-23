%test_trimat

% Testing algorithm 4.1 on different matrices using different functions

% The code generates the data reported in Table 4.1 and 4.2

addpath('data');
addpath('external');
addpath('include');
warning('off');

format shorte
format compact

rng(1)

% set the size of the test matrices
m1 = 35;
m2 = 75;

% generate the test matrices
T1_m1 = gallery('kahan',m1);
T1_m2 = gallery('kahan',m2);
T2_m1 = schur(gallery('smoke',m1),'complex');
T2_m2 = schur(gallery('smoke',m2),'complex');
T3_m1 = schur(randn(m1),'complex');
T3_m2 = schur(randn(m2),'complex');
T4_m1 = schur(rand(m1),'complex');
T4_m2 = schur(rand(m2),'complex');
T5_m1 = triu(randn(m1));
T5_m2 = triu(randn(m2));
T6_m1 = triu(rand(m1));
T6_m2 = triu(rand(m2));
T7_m1 = gallery('jordbloc',m1,0.5);
T7_m2 = gallery('jordbloc',m2,0.5);

% initialization

data_exp_sqrt = zeros(14,8);
data_sign_log = zeros(14,8);
data_cos_sin = zeros(14,8);


main_loop = tic; % record the time consumption

fprintf('testing kahan, T1, m=%d...\n',m1)
data_exp_sqrt(1,1:4) = maxerr_trimat(T1_m1,@exp);
data_exp_sqrt(1,5:8) = maxerr_trimat(T1_m1,@sqrt);
data_sign_log(1,1:4) = maxerr_trimat(T1_m1,@sign);
data_sign_log(1,5:8) = maxerr_trimat(T1_m1,@log);
data_cos_sin(1,1:4) = maxerr_trimat(T1_m1,@cos);
data_cos_sin(1,5:8) = maxerr_trimat(T1_m1,@sin);

fprintf('testing kahan, T1, m=%d...\n',m2)
data_exp_sqrt(2,1:4) = maxerr_trimat(T1_m2,@exp);
data_exp_sqrt(2,5:8) = maxerr_trimat(T1_m2,@sqrt);
data_sign_log(2,1:4) = maxerr_trimat(T1_m2,@sign);
data_sign_log(2,5:8) = maxerr_trimat(T1_m2,@log);
data_cos_sin(2,1:4) = maxerr_trimat(T1_m2,@cos);
data_cos_sin(2,5:8) = maxerr_trimat(T1_m2,@sin);

fprintf('testing schur_smoke, T2, m=%d...\n',m1)
data_exp_sqrt(3,1:4) = maxerr_trimat(T2_m1,@exp);
data_exp_sqrt(3,5:8) = maxerr_trimat(T2_m1,@sqrt);
data_sign_log(3,1:4) = maxerr_trimat(T2_m1,@sign);
data_sign_log(3,5:8) = maxerr_trimat(T2_m1,@log);
data_cos_sin(3,1:4) = maxerr_trimat(T2_m1,@cos);
data_cos_sin(3,5:8) = maxerr_trimat(T2_m1,@sin);

fprintf('testing schur_smoke, T2, m=%d...\n',m2)
data_exp_sqrt(4,1:4) = maxerr_trimat(T2_m2,@exp);
data_exp_sqrt(4,5:8) = maxerr_trimat(T2_m2,@sqrt);
data_sign_log(4,1:4) = maxerr_trimat(T2_m2,@sign);
data_sign_log(4,5:8) = maxerr_trimat(T2_m2,@log);
data_cos_sin(4,1:4) = maxerr_trimat(T2_m2,@cos);
data_cos_sin(4,5:8) = maxerr_trimat(T2_m2,@sin);

fprintf('testing schur_randn, T3, m=%d...\n',m1)
data_exp_sqrt(5,1:4) = maxerr_trimat(T3_m1,@exp);
data_exp_sqrt(5,5:8) = maxerr_trimat(T3_m1*(1+1i),@sqrt);
data_sign_log(5,1:4) = maxerr_trimat(T3_m1,@sign);
data_sign_log(5,5:8) = maxerr_trimat(T3_m1*(1+1i),@log);
data_cos_sin(5,1:4) = maxerr_trimat(T3_m1,@cos);
data_cos_sin(5,5:8) = maxerr_trimat(T3_m1,@sin);

fprintf('testing schur_randn, T3, m=%d...\n',m2)
data_exp_sqrt(6,1:4) = maxerr_trimat(T3_m2,@exp);
data_exp_sqrt(6,5:8) = maxerr_trimat(T3_m2*(1+1i),@sqrt);
data_sign_log(6,1:4) = maxerr_trimat(T3_m2,@sign);
data_sign_log(6,5:8) = maxerr_trimat(T3_m2*(1+1i),@log);
data_cos_sin(6,1:4) = maxerr_trimat(T3_m2,@cos);
data_cos_sin(6,5:8) = maxerr_trimat(T3_m2,@sin);

fprintf('testing schur_rand, T4, m=%d...\n',m1)
data_exp_sqrt(7,1:4) = maxerr_trimat(T4_m1,@exp);
data_exp_sqrt(7,5:8) = maxerr_trimat(T4_m1*(1+1i),@sqrt);
data_sign_log(7,1:4) = maxerr_trimat(T4_m1,@sign);
data_sign_log(7,5:8) = maxerr_trimat(T4_m1*(1+1i),@log);
data_cos_sin(7,1:4) = maxerr_trimat(T4_m1,@cos);
data_cos_sin(7,5:8) = maxerr_trimat(T4_m1,@sin);

fprintf('testing schur_rand, T4, m=%d...\n',m2)
data_exp_sqrt(8,1:4) = maxerr_trimat(T4_m2,@exp);
data_exp_sqrt(8,5:8) = maxerr_trimat(T4_m2*(1+1i),@sqrt);
data_sign_log(8,1:4) = maxerr_trimat(T4_m2,@sign);
data_sign_log(8,5:8) = maxerr_trimat(T4_m2*(1+1i),@log);
data_cos_sin(8,1:4) = maxerr_trimat(T4_m2,@cos);
data_cos_sin(8,5:8) = maxerr_trimat(T4_m2,@sin);

fprintf('testing triu_randn, T5, m=%d...\n',m1);
data_exp_sqrt(9,1:4) = maxerr_trimat(T5_m1,@exp);
data_exp_sqrt(9,5:8) = maxerr_trimat(T5_m1*(1+1i),@sqrt);
data_sign_log(9,1:4) = maxerr_trimat(T5_m1,@sign);
data_sign_log(9,5:8) = maxerr_trimat(T5_m1*(1+1i),@log);
data_cos_sin(9,1:4) = maxerr_trimat(T5_m1,@cos);
data_cos_sin(9,5:8) = maxerr_trimat(T5_m1,@sin);

fprintf('testing triu_randn, T5, m=%d...\n',m2)
data_exp_sqrt(10,1:4) = maxerr_trimat(T5_m2,@exp);
data_exp_sqrt(10,5:8) = maxerr_trimat(T5_m2*(1+1i),@sqrt);
data_sign_log(10,1:4) = maxerr_trimat(T5_m2,@sign);
data_sign_log(10,5:8) = maxerr_trimat(T5_m2*(1+1i),@log);
data_cos_sin(10,1:4) = maxerr_trimat(T5_m2,@cos);
data_cos_sin(10,5:8) = maxerr_trimat(T5_m2,@sin);

fprintf('testing triu_rand, T6, m=%d...\n',m1)
data_exp_sqrt(11,1:4) = maxerr_trimat(T6_m1,@exp);
data_exp_sqrt(11,5:8) = maxerr_trimat(T6_m1,@sqrt);
data_sign_log(11,1:4) = maxerr_trimat(T6_m1,@sign);
data_sign_log(11,5:8) = maxerr_trimat(T6_m1,@log);
data_cos_sin(11,1:4) = maxerr_trimat(T6_m1,@cos);
data_cos_sin(11,5:8) = maxerr_trimat(T6_m1,@sin);

fprintf('testing triu_rand, T6, m=%d...\n',m2)
data_exp_sqrt(12,1:4) = maxerr_trimat(T6_m2,@exp);
data_exp_sqrt(12,5:8) = maxerr_trimat(T6_m2,@sqrt);
data_sign_log(12,1:4) = maxerr_trimat(T6_m2,@sign);
data_sign_log(12,5:8) = maxerr_trimat(T6_m2,@log);
data_cos_sin(12,1:4) = maxerr_trimat(T6_m2,@cos);
data_cos_sin(12,5:8) = maxerr_trimat(T6_m2,@sin);

fprintf('testing jordbloc, T7, m=%d...\n',m1)
data_exp_sqrt(13,1:4) = maxerr_trimat(T7_m1,@exp);
data_exp_sqrt(13,5:8) = maxerr_trimat(T7_m1,@sqrt);
data_sign_log(13,1:4) = maxerr_trimat(T7_m1,@sign);
data_sign_log(13,5:8) = maxerr_trimat(T7_m1,@log);
data_cos_sin(13,1:4) = maxerr_trimat(T7_m1,@cos);
data_cos_sin(13,5:8) = maxerr_trimat(T7_m1,@sin);

fprintf('testing jordbloc, T7, m=%d...\n',m2)
data_exp_sqrt(14,1:4) = maxerr_trimat(T7_m2,@exp);
data_exp_sqrt(14,5:8) = maxerr_trimat(T7_m2,@sqrt);
data_sign_log(14,1:4) = maxerr_trimat(T7_m2,@sign);
data_sign_log(14,5:8) = maxerr_trimat(T7_m2,@log);
data_cos_sin(14,1:4) = maxerr_trimat(T7_m2,@cos);
data_cos_sin(14,5:8) = maxerr_trimat(T7_m2,@sin);

fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);

% save the results to a file
close all
myformat = '$T_, m=$ & %8.1e & %8.1e & %8.1e & %8.1e & %8.1e & %8.1e \n';
fid = fopen('data\test_trimat.txt','w');
for k=1:14
    fprintf(fid,myformat,data_exp_sqrt(k,[1 2 3 5 6 7]))
end
for k=1:14
    fprintf(fid,myformat,data_sign_log(k,[1 2 3 5 6 7]))
end
for k=1:14
    fprintf(fid,myformat,data_cos_sin(k,[1 2 3 5 6 7]))
end
fclose(fid);


function data = maxerr_trimat(T,fun)
% This code returns a vector data containing: the maximal relative errors 
% in 10 calls of the algorithms in section 4.1 and 4.2, 
% condition number of the matrix function fun, and the number of 
% digits used in precision uh in the computation.

if isequal(fun,@sign), refm = @signm; condfun = @signm; end
if isequal(fun,@cos), refm = @cosm; condfun = @cosm; end
if isequal(fun,@sin), refm = @sinm; condfun = @sinm; end
if isequal(fun,@cosh), refm = @coshm; condfun = @coshm; end
if isequal(fun,@sinh), refm = @sinhm; condfun = @sinhm; end
if isequal(fun,@exp), refm = @expm; condfun = @expm; end
if isequal(fun,@log), refm = @logm; condfun = @logm; end
if isequal(fun,@sqrt), refm = @sqrtm; condfun = @sqrtm; end

total_samp = 10; % number of samples
err_vec = zeros(total_samp,2);
% num_samp = zeros(1,2);

d_old = mp.Digits();
mp.Digits(200);
refX = double(refm(mp(T)));

u = eps(1/2);
mp.Digits(17); % sinm and cosm can only be called in mp precisions
condfu = double(funm_condest1(mp(T),condfun))*u;
for i=1:total_samp
    [err_vec(i,1),err_vec(i,2),digits] = comput_error(T,fun,refX);
%     if err_vec(i,1)<condfu
%         num_samp(1) = num_samp(1) + 1;
%     end
%     if err_vec(i,2)<condfu
%         num_samp(2) = num_samp(2) + 1;
%     end
end
max_err = max(err_vec);
% ratio = num_samp/total_samp;
data = [max_err,condfu,digits]; % take the maximal error in the table
% data = [max_err,ratio,condfu,digits];
mp.Digits(d_old);
end

function [error_diag,error_full,digits] = comput_error(T,fun,refX)
[F_diag,digits] = trim_diagpertub(T,fun);
F_full = trim_fullpertub(T,fun);
norm_refX = norm(refX,'fro');
error_diag = norm(F_diag-refX,'fro')/norm_refX;
error_full = norm(F_full-refX,'fro')/norm_refX;
end
