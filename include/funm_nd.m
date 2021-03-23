function F = funm_nd(A,fun,delta)
%FUNM_ND Evaluate general matrix function. 
%   FUNM_ND(A,FUN) evaluates the function_handle FUN at the square 
%   matrix A.  It is a Schur--Parlett algorithm that computes the
%   nontrivial diagonal blocks in the Schur form using the randomized
%   approximate diagonalization method with a diagonal perturbation.
%   delta is the blocking parameter. If delta is empty the default choice
%   delta = 0.1 is used.  

% Default parameters.
if nargin < 3 || isempty(delta)
    delta = 0.1;
end
ord = [];

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

if isequal(delta,'norm')
    max_tii = norm(diagT,Inf);
    delta = 0.1*max_tii;
end

% Determine reordering of Schur form into block form.
if isempty(ord), ord = blocking(T,delta); end

[ord, ind] = swapping(ord);  % Gives the blocking.
ord = max(ord)-ord+1;        % Since ORDSCHUR puts highest index top left.
[U,T] = ordschur(U,T,ord);

m = length(ind);

% Calculate F(T)
F = zeros(n);
for col=1:m
   j = ind{col};
   F(j,j) = trim_diagpertub(T(j,j),fun);
   for row=col-1:-1:1
      i = ind{row};
      if length(i) == 1 && length(j) == 1
         % Scalar case.
         k = i+1:j-1;
         temp = T(i,j)*(F(i,i) - F(j,j)) + F(i,k)*T(k,j) - T(i,k)*F(k,j);
         F(i,j) = temp/(T(i,i)-T(j,j));
      else
         k = cat(2,ind{row+1:col-1});
         rhs = F(i,i)*T(i,j) - T(i,j)*F(j,j) + F(i,k)*T(k,j) - T(i,k)*F(k,j);
         F(i,j) = sylv_tri(T(i,i),-T(j,j),rhs);
      end
   end
end

F = U*F*U';

if isreal(A) && norm(imag(F),1) <= 10*n*eps*norm(F,1)
   F = real(F);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mm,ind] = swapping(m)
%SWAPPING  Choose confluent permutation ordered by average index.
%   [MM,IND] = SWAPPING(M) takes a vector M containing the integers
%   1:K (some repeated if K < LENGTH(M)), where M(J) is the index of
%   the block into which the element T(J,J) of a Schur form T
%   should be placed.
%   It constructs a vector MM (a permutation of M) such that T(J,J)
%   will be located in the MM(J)'th block counting from the (1,1) position.
%   The algorithm used is to order the blocks by ascending
%   average index in M, which is a heuristic for minimizing the number
%   of swaps required to achieve this confluent permutation.
%   The cell array vector IND defines the resulting block form:
%   IND{i} contains the indices of the i'th block in the permuted form.

mmax = max(m); mm = zeros(size(m));
g = zeros(1,mmax); h = zeros(1,mmax);

for i = 1:mmax
    p = find(m==i);
    h(i) = length(p);
    g(i) = sum(p)/h(i);
end

[~,y] = sort(g);
h = [0 cumsum(h(y))];

ind = cell(mmax,1);
for i = 1:mmax
    mm(m==y(i)) = i;
    ind{i} = h(i)+1:h(i+1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = sylv_tri(T,U,B)
%SYLV_TRI    Solve triangular Sylvester equation.
%   X = SYLV_TRI(T,U,B) solves the Sylvester equation
%   T*X + X*U = B, where T and U are square upper triangular matrices.

m = length(T);
n = length(U);
X = zeros(m,n);
opts.UT = true;

% Forward substitution.
for i = 1:n
    X(:,i) = linsolve(T + U(i,i)*eye(m), B(:,i) - X(:,1:i-1)*U(1:i-1,i), opts);
end
end
