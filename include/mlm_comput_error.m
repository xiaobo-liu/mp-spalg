function [error_nd,error_mlm] = mlm_comput_error(A,alpha,beta)
% Copmute the errors produced by funm_nd and mlm for copmuting the 
% matrix Mittag-Leffler functions with two parametes alpha and beta of
% matrix A.

F_nd = mlm_nd(A,alpha,beta);
F_mlm = mlm_mlm(A,alpha,beta);
F = mlm_ref(A,alpha,beta); % reference solution

% compute normwise relative errors in F-norm
F_norm = norm(F,'fro');
error_nd = double(norm(F_nd-F,'fro')/F_norm);
error_mlm = double(norm(F_mlm-F,'fro')/F_norm);
end

function mlm = mlm_ref(A,alpha,beta)
% Function to compute the reference solution for the matrix Mittag-Leffler
% function by using the diagonalization approach at 200 digits precision,
% and then rounding back to double precision.

d_comput = 200; % set the digits to use in diagonalization 
d_old = mp.Digits();
mp.Digits(d_comput);
A = mp(A);

% E.B.Davies' trick: ensures distinct eigenvalues; diagonalisation possible
delA = mp(randn(length(A))); % generate random perturbation
delA = 10^(-d_comput/2)*delA/norm(delA,1);
%   [V, D] = eig(A); if a is not diagonalisable, this may break down
[V, D] = eig(A - delA); % E.B.Davies's trick
% convert back into double presion
mlm = double(V*diag(ml_truncat(diag(D),alpha,beta,d_comput))/V); 
mp.Digits(d_old);
end
