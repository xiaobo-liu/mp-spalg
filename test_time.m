% The code generates a matrix table that contains all the data reported 
% in Table 5.2.
%
% IMFORMATION - This test took about 15 mins on a laptop equipped with an
% Intel i7-6700HQ processor running at 2.60GHz and with 16GB of RAM.

addpath('data');
addpath('external');
addpath('include');
warning('off');

format shorte
format compact

% reproduce the results?
produce_results = false;

% initialization
rng(0)

error_A1_sin = zeros(2, 3);
error_A2_sin = zeros(2, 3);
error_A3_sin = zeros(2, 3);
error_A1_cosh = zeros(2, 3);
error_A2_cosh = zeros(2, 3);
error_A3_cosh = zeros(2, 3);

time_A1_sin = zeros(2, 3);
time_A2_sin = zeros(2, 3);
time_A3_sin = zeros(2, 3);
time_A1_cosh = zeros(2, 3);
time_A2_cosh = zeros(2, 3);
time_A3_cosh = zeros(2, 3);

error_vec_sin = zeros(10,3); % a matrix to store the errors temporarily
error_vec_cosh = zeros(10,3);

size_digits_A1_sin = zeros(2,2);
size_digits_A1_cosh = zeros(2,2);
size_digits_A2_sin = zeros(2,2);
size_digits_A2_cosh = zeros(2,2);
size_digits_A3_sin = zeros(2,2);
size_digits_A3_cosh = zeros(2,2);

table = zeros(12, 8);

if produce_results
    i = 0;
    main_loop = tic; % record the time consumption
    fprintf('n=40\n');
    for n = [40 100]
        i = i+1;
        if i==2, fprintf('n=100\n'); end
        % A = A1
        A = rand(n)/5;
        fprintf('Testing A1...\n');
        for num=1:10
            % f = sin
            [Xfunm, time] = smart_timer(@funm, A, @sin);
            time_A1_sin(i,1) = time_A1_sin(i,1) + time;
            [Xfunm_nd, time] = smart_timer(@funm_nd, A, @sin);
            time_A1_sin(i,2) = time_A1_sin(i,2) + time;
            [Xfunm_nd_inf, time] = smart_timer(@funm_nd_inf, A, @sin);
            time_A1_sin(i,3) = time_A1_sin(i,3) + time;
            old_digits = mp.Digits();
            mp.Digits(100);
            refX_sin = sinm(mp(A));
            refX_cosh = coshm(mp(A));
            mp.Digits(old_digits);
            refX_sin_norm = double(norm(refX_sin,'fro'));
            error_vec_sin(num,1) = norm(Xfunm-refX_sin,'fro')/refX_sin_norm;
            error_vec_sin(num,2) = norm(Xfunm_nd-refX_sin,'fro')/refX_sin_norm;
            error_vec_sin(num,3) = norm(Xfunm_nd_inf-refX_sin,'fro')/refX_sin_norm;
            % f = cosh
            [Xfunm, time] = smart_timer(@funm, A, @cosh);
            time_A1_cosh(i,1) = time_A1_cosh(i,1) + time;
            [Xfunm_nd, time] = smart_timer(@funm_nd, A, @cosh);
            time_A1_cosh(i,2) = time_A1_cosh(i,2) + time;
            [Xfunm_nd_inf, time] = smart_timer(@funm_nd_inf, A, @cosh);
            time_A1_cosh(i,3) = time_A1_cosh(i,3) + time;
            refX_cosh_norm = double(norm(refX_cosh,'fro'));
            error_vec_cosh(num,1) = norm(Xfunm-refX_cosh,'fro')/refX_cosh_norm;
            error_vec_cosh(num,2) = norm(Xfunm_nd-refX_cosh,'fro')/refX_cosh_norm;
            error_vec_cosh(num,3) = norm(Xfunm_nd_inf-refX_cosh,'fro')/refX_cosh_norm;
        end
        time_A1_sin(i,:) = time_A1_sin(i,:)/num; % take the mean time
        time_A1_cosh(i,:) = time_A1_cosh(i,:)/num;
        for j = 1:3
            error_A1_sin(i,j) = max(error_vec_sin(:,j)); % take the maximal error
            error_A1_cosh(i,j) = max(error_vec_cosh(:,j));
        end
        [~,max_digits,max_block] = funm_nd(A,@sin);
        size_digits_A1_sin(i,:) = [max_block, max_digits];
        [~,max_digits,max_block] = funm_nd(A,@cosh);
        size_digits_A1_cosh(i,:) = [max_block, max_digits];
        % A = A2
        A = randn(n)/10;
        fprintf('Testing A2...\n');
        for num=1:10
            % f = sin
            [Xfunm, time] = smart_timer(@funm, A, @sin);
            time_A2_sin(i,1) = time_A2_sin(i,1) + time;
            [Xfunm_nd, time] = smart_timer(@funm_nd, A, @sin);
            time_A2_sin(i,2) = time_A2_sin(i,2) + time;
            [Xfunm_nd_inf, time] = smart_timer(@funm_nd_inf, A, @sin);
            time_A2_sin(i,3) = time_A2_sin(i,3) + time;
            old_digits = mp.Digits();
            mp.Digits(100);
            refX_sin = sinm(mp(A));
            refX_cosh = coshm(mp(A));
            mp.Digits(old_digits);
            refX_sin_norm = double(norm(refX_sin,'fro'));
            error_vec_sin(num,1) = norm(Xfunm-refX_sin,'fro')/refX_sin_norm;
            error_vec_sin(num,2) = norm(Xfunm_nd-refX_sin,'fro')/refX_sin_norm;
            error_vec_sin(num,3) = norm(Xfunm_nd_inf-refX_sin,'fro')/refX_sin_norm;
            % f = cosh
            [Xfunm, time] = smart_timer(@funm, A, @cosh);
            time_A2_cosh(i,1) = time_A2_cosh(i,1) + time;
            [Xfunm_nd, time] = smart_timer(@funm_nd, A, @cosh);
            time_A2_cosh(i,2) = time_A2_cosh(i,2) + time;
            [Xfunm_nd_inf, time] = smart_timer(@funm_nd_inf, A, @cosh);
            time_A2_cosh(i,3) = time_A2_cosh(i,3) + time;
            refX_cosh_norm = double(norm(refX_cosh,'fro'));
            error_vec_cosh(num,1) = norm(Xfunm-refX_cosh,'fro')/refX_cosh_norm;
            error_vec_cosh(num,2) = norm(Xfunm_nd-refX_cosh,'fro')/refX_cosh_norm;
            error_vec_cosh(num,3) = norm(Xfunm_nd_inf-refX_cosh,'fro')/refX_cosh_norm;
        end
        time_A2_sin(i,:) = time_A2_sin(i,:)/num; % take the mean time
        time_A2_cosh(i,:) = time_A2_cosh(i,:)/num;
        for j = 1:3
            error_A2_sin(i,j) = max(error_vec_sin(:,j)); % take the maximal error
            error_A2_cosh(i,j) = max(error_vec_cosh(:,j));
        end
        [~,max_digits,max_block] = funm_nd(A,@sin);
        size_digits_A2_sin(i,:) = [max_block, max_digits];
        [~,max_digits,max_block] = funm_nd(A,@cosh);
        size_digits_A2_cosh(i,:) = [max_block, max_digits];
        % A = A3
        A = gallery('triw',n,-5);
        fprintf('Testing A3...\n');
        for num=1:10
            % f = sin
            [Xfunm, time] = smart_timer(@funm, A, @sin);
            time_A3_sin(i,1) = time_A3_sin(i,1) + time;
            [Xfunm_nd, time] = smart_timer(@funm_nd, A, @sin);
            time_A3_sin(i,2) = time_A3_sin(i,2) + time;
            [Xfunm_nd_inf, time] = smart_timer(@funm_nd_inf, A, @sin);
            time_A3_sin(i,3) = time_A3_sin(i,3) + time;
            old_digits = mp.Digits();
            mp.Digits(100);
            refX_sin = sinm(mp(A));
            refX_cosh = coshm(mp(A));
            mp.Digits(old_digits);
            refX_sin_norm = double(norm(refX_sin,'fro'));
            error_vec_sin(num,1) = norm(Xfunm-refX_sin,'fro')/refX_sin_norm;
            error_vec_sin(num,2) = norm(Xfunm_nd-refX_sin,'fro')/refX_sin_norm;
            error_vec_sin(num,3) = norm(Xfunm_nd_inf-refX_sin,'fro')/refX_sin_norm;
            % f = cosh
            [Xfunm, time] = smart_timer(@funm, A, @cosh);
            time_A3_cosh(i,1) = time_A3_cosh(i,1) + time;
            [Xfunm_nd, time] = smart_timer(@funm_nd, A, @cosh);
            time_A3_cosh(i,2) = time_A3_cosh(i,2) + time;
            [Xfunm_nd_inf, time] = smart_timer(@funm_nd_inf, A, @cosh);
            time_A3_cosh(i,3) = time_A3_cosh(i,3) + time;
            refX_cosh_norm = double(norm(refX_cosh,'fro'));
            error_vec_cosh(num,1) = norm(Xfunm-refX_cosh,'fro')/refX_cosh_norm;
            error_vec_cosh(num,2) = norm(Xfunm_nd-refX_cosh,'fro')/refX_cosh_norm;
            error_vec_cosh(num,3) = norm(Xfunm_nd_inf-refX_cosh,'fro')/refX_cosh_norm;
        end
        time_A3_sin(i,:) = time_A3_sin(i,:)/num; % take the mean time
        time_A3_cosh(i,:) = time_A3_cosh(i,:)/num;
        for j = 1:3
            error_A3_sin(i,j) = max(error_vec_sin(:,j)); % take the maximal error
            error_A3_cosh(i,j) = max(error_vec_cosh(:,j));
        end
        [~,max_digits,max_block] = funm_nd(A,@sin);
        size_digits_A3_sin(i,:) = [max_block, max_digits];
        [~,max_digits,max_block] = funm_nd(A,@cosh);
        size_digits_A3_cosh(i,:) = [max_block, max_digits];
    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
    table = [error_A1_sin(1,:),time_A1_sin(1,:),size_digits_A1_sin(1,:);
             error_A2_sin(1,:),time_A2_sin(1,:),size_digits_A2_sin(1,:);
             error_A3_sin(1,:),time_A3_sin(1,:),size_digits_A3_sin(1,:);
             error_A1_sin(2,:),time_A1_sin(2,:),size_digits_A1_sin(2,:);
             error_A2_sin(2,:),time_A2_sin(2,:),size_digits_A2_sin(2,:);
             error_A3_sin(2,:),time_A3_sin(2,:),size_digits_A3_sin(2,:);
             error_A1_cosh(1,:),time_A1_cosh(1,:),size_digits_A1_cosh(1,:);
             error_A2_cosh(1,:),time_A2_cosh(1,:),size_digits_A2_cosh(1,:);
             error_A3_cosh(1,:),time_A3_cosh(1,:),size_digits_A3_cosh(1,:);
             error_A1_cosh(2,:),time_A1_cosh(2,:),size_digits_A1_cosh(2,:);
             error_A2_cosh(2,:),time_A2_cosh(2,:),size_digits_A2_cosh(2,:);
             error_A3_cosh(2,:),time_A3_cosh(2,:),size_digits_A3_cosh(2,:);];
    % round the times and errors to 2 significant digits
    table(:,1:6) = round(table(:,1:6),2,'significant'); 
    ratios = table([1 2 4 5 7 8 10 11],5) ./ table([1 2 4 5 7 8 10 11],4);
    [ratio_max, index] = max(ratios); 
    save('data/table.mat','table', 'ratio_max', 'index');
else
    load ('data/table.mat');
end

table
ratio_max
index    

function [X, time] = smart_timer(f, varargin)
% This subfunction is due to Massimiliano Fasi.
    init = tic;
    X = f(varargin{:});
    time = toc(init);
    if time < 0.5
      timefun = @()(f(varargin{:}));
      time = timeit(timefun);
    end
end

function [F,max_digits,mb_size] = funm_nd(A,fun,delta)
%FUNM_ND Evaluate general matrix function 
%(a modified version that also reruens the maximal block size and 
% the maximal digits used in computation).
%   FUNM_ND(A,fun,delta), a Schur--Parlett algorithm that computes 
%   the nontrivial diagonal blocks in the Schur form using the 
%   randomized approximate diagonalization method with a triangular 
%   perturbation. 
%
%   F approximates fun(A) and delta is the blocking parameter in the 
%   Schur-Parlett algorithm. If delta is empty the default choice 
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

% Determine reordering of Schur form into block form.
if isempty(ord), ord = blocking(T,delta); end

[ord, ind] = swapping(ord);  % Gives the blocking.
ord = max(ord)-ord+1;        % Since ORDSCHUR puts highest index top left.
[U,T] = ordschur(U,T,ord);

m = length(ind);

% Calculate F(T)
F = zeros(n);
d_uh = zeros(1,m);
for col=1:m
   j = ind{col};
   [F(j,j),d_uh(col)] = trim_diagpertub(T(j,j),fun);
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
max_digits = max(d_uh);
mb_size = max(accumarray(ord',1));
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

