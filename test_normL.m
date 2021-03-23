%test_normL

% test the matrix square root for different matrices.  No extre precision.

% This code generates the data reported in Table 3.1 and 3.2.

format compact
format shorte

rng(1)
warning off MATLAB:nearlySingularMatrix

clear err normL term

u = eps/2;

i = 0;
for lambda = [1 0.5 0.1]
i = i+1;

j = 0;
for n = [10 20 30]
j = j+1;
% change the matrix 
% A = gallery('kahan',n); % T1
% A = schur(gallery('smoke',n),'complex'); % T2 
% A = schur(randn(n),'complex'); % T3
% A = schur(rand(n),'complex'); % T4
% A = triu(randn(n)); % T5 
% A = triu(rand(n)); % T6 
A = gallery('jordbloc',n,lambda); %T7 
% A = gallery('triw',n); % T8
tol = norm(A,'fro')*sqrt(u);

E = randn(n);
E1 = triu(E);
E = tol*E/norm(E,'fro');
E1 = tol*E1/norm(E1,'fro');

mp.Digits(200); 
% [Fx1,alpha,cond_sqrt] = sqrtm(A);
% alpha, cond_sqrt
Fx1 = sqrtm(mp(A));
Fx = double(Fx1);

% Too slow for n = 30:
% normL(i,j) = double( norm(inv(kron(eye(n),Fx1) + kron(Fx1',eye(n))),2) );

% mp.Digits(34);
% mp.Digits(100); %makes no difference.

K = kron(eye(n),Fx) + kron(Fx',eye(n));
normL(i,j) = norm(inv(K),2);

G = randn(n); G = G/norm(G,1); Gt = triu(G);
Lsolves = [norm(K\G(:),1), norm(K\Gt(:),1)]

% term(i,j) = normL(i,j)/norm(Fx,'fro');

Amp = A + E;
[Vx,Dx] = eig(mp(Amp));
condVx = double(cond(Vx))

[V,D] = eig(Amp);
Vd = double(V); Dd = double(D);
condV = double(cond(V))
Fmp = V*diag(sqrt(diag(D)))/V;
F = double(Fmp);
err(i,j) = norm(F-Fx,'fro') / norm(Fx,'fro');

[V,D] = eig(A + E1);  % Triangular pert
condV = double(cond(V))
Fmp = V*diag(sqrt(diag(D)))/V;
F = double(Fmp);
err2(i,j) = norm(F-Fx,'fro') / norm(Fx,'fro');

% Carry out the evaluation at high precision. Doesn't help!
Amp = mp(A) + mp(E);
[V,D] = eig(Amp);
Fmp = V*diag(sqrt(diag(D)))/V;
F = double(Fmp);
err3(i,j) = norm(F-Fx,'fro') / norm(Fx,'fro');

end
end

err
normL
% term
% termXE = term*norm(E,'fro');

% ltx = 1; ex = 1; mth = 1;
% print_matrix(err, [], [], [], ltx, ex, mth)
% fprintf('\n')
% print_matrix(normL, [], [], [], ltx, ex, mth)

warning on MATLAB:nearlySingularMatrix
