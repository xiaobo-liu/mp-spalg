function [F,d_uh] = mlm_trim_diagpertub(T,alpha,beta)
%ML_TRIM_DIAGPERTUB A modified version of TRIM_DIAGPERTUB for the 
%matrix ML function.
%   We need this since the scalar ML function to different digits of 
%   significants may be needed in the computation.
%   ml_trim_tripertub(T,alpha,beta) computes the matrix ML function
%   with two parameters alpha and beta of a triangular matrix T 
%   using the randomized approximate diagonalization method with a
%   diagonal perturbation. F approximates fun(T) and d_uh returns
%   returns digits of the possibly higher than u^2 precision.

d_uh = 32; % the defult output
[m,n] = size(T);
if  ~isfloat(T) || ~ismatrix(T) || m ~= n
   error(message('MATLAB:funm:InputDim'));
end
if ~matlab.internal.math.isschur(T,'complex')
    error('MATLAB:MatrixNotTriangular',...
        'The input matrix must be triangular.');
end

% parameters: thera and delta_1
theta = 0.4;
delta1 = 5e-3;
% % use precision u^2 for computing 1*1 and 2*2 blocks
% fun = @(z) ml_truncat(z,alpha,beta,32);
fun = @(z) ml_truncat(z,alpha,beta,16);

% Handle special case of diagonal T.
diagT = diag(T);
if isequal(T,diag(diagT)) 
   F = double(diag(fun(diagT)));
   return
end

% Handle 1*1 and 2*2 blocks
u = eps(1/2);
if m == 1, F = double(feval(fun,T)); return, end
t11 = T(1,1); t22 = T(2,2);
if m == 2 && abs(t11-t22)>u
    t12 = T(1,2);
    f11 = fun(t11); f22 = fun(t22);
    F = [f11  t12*(f22-f11)/(t22-t11);
         0      f22]; 
     F = double(F);
     return
end

E = diag(randn(m,1));
max_tij = max(abs(T),[],'all');
E = E/norm(E,'fro')*u*max_tij;

d_old = mp.Digits(); % 'mp digits' at the start

% Use precision u^2 to form T + E since E might be 
% too small to be added by T in precision u
d_defult = 32; % precision with unit roundoff u^2
mp.Digits(d_defult);
Ttilde = mp(T)+mp(E);
fun = @(z) ml_truncat(z,alpha,beta,d_defult);

% calculate the largest group size k
ord = blocking(Ttilde,delta1);
k = max(accumarray(ord',1));

if k>1  % evaluate the required precision uh
    diagTtilde = diag(Ttilde);
    D = diag(diagTtilde);
    max_tij_tilde = max(abs(Ttilde-D),[],'all');
    cm = theta*max_tij/sqrt(m);
    uh = cm*u^2/max_tij_tilde*(max_tij_tilde/cm/u+1)^(2-k);
    d_uh = double(ceil(log10(1/uh)));
    d_uh = max(d_uh,d_defult);
    mp.Digits(d_uh);
    fun = @(z) ml_truncat(z,alpha,beta,d_uh);
end

d_uh = mp.Digits;
% increase the precision from u^2 to uh in solving V when necessary
Ttilde = mp(Ttilde); 
[V,Dtilde] = eig(Ttilde); 
F = double(V*diag(fun(diag(Dtilde)))/V);
mp.Digits(d_old); % return 'mp digits' to the value it had at the start

if isreal(T) && norm(imag(F),1) <= 10*n*eps*norm(F,1)
   F = real(F);
end
end

% blocking scheme taken from the MATLAB funm function
function m = blocking(A,delta)
%BLOCKING  Produce blocking pattern for block Parlett recurrence in FUNM.
%   M = BLOCKING(A, DELTA, SHOWPLOT) accepts an upper triangular matrix
%   A and produces a blocking pattern, specified by the vector M,
%   for the block Parlett recurrence.
%   M(i) is the index of the block into which A(i,i) should be placed,
%   for i=1:LENGTH(A).
%   DELTA is a gap parameter (default 0.1) used to determine the blocking.

%   For A coming from a real matrix it should be posible to take
%   advantage of the symmetry about the real axis.  This code does not.

a = diag(A); n = length(a);
m = zeros(1,n); maxM = 0;

if nargin < 2 || isempty(delta), delta = 0.1; end

for i = 1:n

    if m(i) == 0
        m(i) = maxM + 1; % If a(i) hasn`t been assigned to a set
        maxM = maxM + 1; % then make a new set and assign a(i) to it.
    end

    for j = i+1:n
        if m(i) ~= m(j)    % If a(i) and a(j) are not in same set.
            if abs(a(i)-a(j)) <= delta

                if m(j) == 0
                    m(j) = m(i); % If a(j) hasn`t been assigned to a
                                 % set, assign it to the same set as a(i).
                else
                    p = max(m(i),m(j)); q = min(m(i),m(j));
                    m(m==p) = q; % If a(j) has been assigned to a set
                                 % place all the elements in the set
                                 % containing a(j) into the set
                                 % containing a(i) (or vice versa).
                    m(m>p) = m(m>p) -1;
                    maxM = maxM - 1;
                                 % Tidying up. As we have deleted set
                                 % p we reduce the index of the sets
                                 % > p by 1.
                end
            end
        end
    end
end
end
