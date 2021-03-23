function F = trim_fullpertub(T,fun)
%TRIM_FULLPERTUB compute functions of triangular matrix using a full 
%perturbation.
%   TRIM_FULLPERTUB(T,fun) computes function fun of a triangular matrix T 
%   using Davies's randomized approximate diagonalization method (with a
%   full perturbation). F approximates fun(T).

[m,n] = size(T);

if  ~isfloat(T) || ~ismatrix(T) || m ~= n
   error(message('MATLAB:funm:InputDim'));
end

if ~matlab.internal.math.isschur(T,'complex')
    error('MATLAB:MatrixNotTriangular',...
        'The input matrix must be triangular.');
end
if isequal(fun,@sign) % Handle special case of f = sign
    fun = @(x) sign(real(x));
end

u = eps(1/2);
% Full perturbation
% E is real if T is real
% if isreal(T) && norm(imag(T),1) <= 10*m*eps*norm(T,1) 
    E = randn(m);
% else
%     E = randn(m) + 1i*randn(m);
% end
max_tij = max(abs(T),[],'all');
E = E/norm(E,'fro')*u*max_tij; % norm of E = u*max_ij(t_ij)


d_old = mp.Digits(); % 'mp digits' at the start
mp.Digits(32); % precision u^2
Ttilde = mp(T)+mp(E);

% Compute V by calling eig
[V,D] = eig(Ttilde); 

F = double(V*diag(fun(diag(D)))/V);
mp.Digits(d_old); % return 'mp digits' to the value it had at the start

if isreal(T) && norm(imag(F),1) <= 10*n*eps*norm(F,1)
   F = real(F);
end
end