function F = funm_nd_inf(A,fun)
% Algorithm that computes function of the whole Schur factor T by the 
% randomized approximate diagonalization method with a diagonal 
% perturbation 

[m,n] = size(A);
if  ~isfloat(A) || ~ismatrix(A) || m ~= n
   error(message('MATLAB:funm:InputDim'));
end

% First form complex Schur form (if A not already upper triangular).
if isequal(A,triu(A))
   T = A; U = eye(n);
   diagT = diag(T);
else
   [U,T] = schur(A,'complex');
   diagT = diag(T);
end

if isequal(T,diag(diagT)) % Handle special case of diagonal T.
   F = U*diag(fun(diagT))*U';
   return
end

% Calculate F(T)
F = trim_diagpertub(T,fun);

F = U*F*U';

if isreal(A) && norm(imag(F),1) <= 10*n*eps*norm(F,1)
   F = real(F);
end

end