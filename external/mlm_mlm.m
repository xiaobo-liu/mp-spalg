function E = mlm_mlm(A,alpha,beta)
%MLM_MLM Algorithm for computing the matrix Mittag-Leffler
%functions by Roberto Garrappa1 and Marina Popolizio.

% Evaluate the Mittag-Leffler (ML) function with two parameters ALPHA and
% BETA of the square matrix argument A.
%
% E = ML(A,ALPHA,BETA) evaluates the ML function with two parameters ALPHA
% and BETA at the square matrix A argument; ALPHA must be any real and 
% positive scalar, BETA any real scalar and A any real or complex square 
% matrix. If the second parameter BETA is missing, it is assumed BETA=1. 
%
% TECHNICAL NOTES:
%
% The ML function on the matrix argument A is evaluated by exploiting the
% Schur-Parlett algorithm and evaluating derivative of the scalar ML
% function by combining, by means of the derivatives balancing technique
% studied in [1], Taylor series, a summation formula based on the Prabhakar 
% function and the numerical inversion of the Laplace transform obtained 
% after generalizing the algorithm described in [2]. For more details users 
% are referred to [1].      
%
% REFERENCES
%
% [1] R. Garrappa and M. Popolizio, Computing the matrix Mittag?Leffler
% function with applications to fractional calculus, submitted
%
% [2] R. Garrappa, Numerical Evaluation of two and three parameter
% Mittag-Leffler functions, SIAM Journal of Numerical Analysis, 2015,
% 53(3), 1350-1369.
%
%  Please, report any problem or comment to : 
%          roberto dot garrappa at uniba dot it
%
%  Copyright (c) 2017
%
%  Authors:
%   Roberto Garrappa (University of Bari, Italy)
%   roberto dot garrappa at uniba dot it
%   Homepage: http://www.dm.uniba.it/Members/garrappa
%
%   Marina Popolizio (University of Salento, Italy)
%   marina dot popolizio at unisalento dot it
%
%
%   Revision: 1.0 - Date: June 5, 2017
% Check inputs
if nargin < 3
    warning('MATLAB:ml_matrix:BetaParameterMissing', ...
        ['The second parameter BETA of the ML function should be given.', ...
        'The value BETA=1 is assumed.']) ;
    beta = 1 ;
    if nargin < 2
        error('MATLAB:ml_matrix:NumberParameters', ...
            'The parameter ALPHA must be specified.');
    end
end
% Check whether the parameters ALPHA, BETA are of scalar type
if length(alpha) > 1 || length(beta) > 1 
    alpha = alpha(1) ; beta = beta(1) ; 
    warning('MATLAB:ml_matrix:ScalarParameters', ...
        ['ALPHA, BETA must be scalar parameters. ', ...
        'Only the first values ALPHA=%f and BETA=%f will be considered. '], ...
        alpha, beta) ;
end
% Check whether the argument is a square matrix
[nA, mA] = size(A) ;
if nA ~= mA
    error('MATLAB:ml_matrix:MatrixNotSquare', ...
        'The matrix argument must be square.');
end
% Check whether the parameters meet the contraints
if real(alpha) <= 0 || ~isreal(alpha) 
    error('MATLAB:ml_matrix:ParametersOutOfRange', ...
        ['Error in one of the parameters of the Mittag-Leffler function. ', ...
        'The first parameter ALPHA must be real and positive. ']) ;
end
if ~isreal(beta) 
    error('MATLAB:ml_matrix:ParametersOutOfRange', ...
        ['Error in one of the parameters of the Mittag-Leffler function. ', ...
        'The second parameter BETA can not be a complex number. ']) ;
end
mldr = @(z,k) mld(z,alpha,beta,k) ;
E = funm(A,mldr) ; 
if isreal(A) 
    E = real(E) ; 
end
end
function Ek = mld(z,alpha,beta,k)
%for j = 1 : length(z)
%fprintf('z(%d)=%f  al=%f  be=%f  k=%d \n',j,z(j),alpha,beta,k)
%end
Ek = zeros(size(z)) ;
tau = 1.0e-14 ; 
max_gamma_arg = 171.624 ;
Jmax = floor((max_gamma_arg - beta)/alpha) ;
z_abs_max = (tau*gamma(alpha*Jmax+beta)/prod(Jmax-(0:k-1)))^(1/(Jmax-k)) ;
i_z_se = abs(z)<z_abs_max ;
z_se = z(i_z_se) ;
E_se = zeros(size(z)) ; Err_Round = zeros(size(z)) ; e = ones(size(z)) ;
[E_se(i_z_se), Err_Round(i_z_se)] = mlds(z_se,alpha,beta,k) ;
i_z_se_accept = (Err_Round./(e+abs(E_se)) < tau) & ~(E_se==0)  ;
Ek(i_z_se_accept) = E_se(i_z_se_accept) ;
i_ze_lt = ~i_z_se_accept ;
if k <= 3, p = 0 ; elseif k <= 7, p = 1; else p = 2 ; end 
c = zeros(k-p+1,k-p+1) ;
c(1,1) = 1 ;
for kk = 1 : k-p
    c(kk+1,1) = (1-alpha*(kk-1)-beta)*c(kk,1) ;
    for j = 1 : kk-1
        c(kk+1,j+1) = c(kk,j) + (j+1-alpha*(kk-1)-beta)*c(kk,j+1) ;
    end
    c(kk+1,kk+1) = 1 ;
end
for j = 0 : k-p
    if abs(c(k-p+1,j+1)) > 1.0e-14
        Ek(i_ze_lt) = Ek(i_ze_lt) + c(k-p+1,j+1)*mldlt(z(i_ze_lt),alpha,(k-p)*alpha+beta-j,p) ;
    end
end
Ek(i_ze_lt) = Ek(i_ze_lt)/alpha^(k-p) ;
end
function E = mldlt(z,alpha,beta,k)
log_epsilon = log(10^(-15)) ; 
E = zeros(size(z)) ;  
for j = 1 : length(z)
    if abs(z(j)) < 1.0e-14
        E(j) = factorial(k)/gamma(alpha*k+beta) ;
    else
        E(j) = LTI(1,z(j),alpha,beta,k,log_epsilon) ;
    end
    if isreal(z(j)), E(j) = real(E(j)) ; end 
end
end
function  E = LTI(t,lambda,alpha,beta,kd,log_epsilon)
theta = angle(lambda) ;
kmin = ceil(-alpha/2 - theta/2/pi) ;
kmax = floor(alpha/2 - theta/2/pi) ;
k_vett = kmin : kmax ;
s_star = abs(lambda)^(1/alpha) * exp(1i*(theta+2*k_vett*pi)/alpha) ;
phi_s_star = (real(s_star)+abs(s_star))/2 ;
[phi_s_star , index_s_star ] = sort(phi_s_star) ;
s_star = s_star(index_s_star) ;
index_save = phi_s_star > 1.0e-15 ;
s_star = s_star(index_save) ;
phi_s_star = phi_s_star(index_save) ;
s_star = [0, s_star] ;
phi_s_star = [0, phi_s_star] ;
J1 = length(s_star) ; J = J1 - 1 ;
if abs(lambda) <= 1.0 && abs(angle(lambda)) > pi*0.9
    p = [ max(0,-2*(-alpha*(kd+1)/2+alpha-beta+1)) , ones(1,J)*(kd+1) ]  ;
else
    p = [ max(0,-2*(alpha-beta+1)) , ones(1,J)*(kd+1) ]  ;
end
q = [ ones(1,J)*(kd+1) , +Inf] ;
phi_s_star = [phi_s_star, +Inf] ;
admissible_regions = find( ...
    (phi_s_star(1:end-1) < (log_epsilon - log(eps))/t) & ...
    (phi_s_star(1:end-1) < phi_s_star(2:end))) ;
JJ1 = admissible_regions(end) ;
mu_vett = ones(1,JJ1)*Inf ;
N_vett = ones(1,JJ1)*Inf ;
h_vett = ones(1,JJ1)*Inf ;
find_region = 0 ;
while ~find_region
    for j1 = admissible_regions
        if j1 < J1
            [muj,hj,Nj] = OptimalParam_RB ...
                (t,phi_s_star(j1),phi_s_star(j1+1),p(j1),q(j1),log_epsilon) ;
        else
            [muj,hj,Nj] = OptimalParam_RU(t,phi_s_star(j1),p(j1),log_epsilon) ;
        end
        mu_vett(j1) = muj ; h_vett(j1) = hj ; N_vett(j1) = Nj ;
    end
    if min(N_vett) > 200
        log_epsilon = log_epsilon + log(10) ;
    else
        find_region = 1 ;
    end
end
[N, iN] = min(N_vett) ; mu = mu_vett(iN) ; h = h_vett(iN) ;
k = -N : N ;
u = h*k ;
z = mu*(1i*u+1).^2 ;
zd = -2*mu*u + 2*mu*1i ;
zexp = exp(z*t) ;
F = z.^(alpha-beta)./(z.^alpha - lambda).^(kd+1).*zd ;
S = zexp.*F ;
Integral = h*sum(S)/2/pi/1i ;
ss_star = s_star(iN+1:end) ;
Residues = mldr(t,ss_star,alpha,beta,kd) ;
E = factorial(kd)*(Integral + sum(Residues)) ;
end
function [muj,hj,Nj] = OptimalParam_RB ...
    (t, phi_s_star_j, phi_s_star_j1, pj, qj, log_epsilon)
log_eps = -36.043653389117154 ; % log(eps)
fac = 1.01 ;
conservative_error_analysis = 0 ;
f_max = exp(log_epsilon - log_eps) ;
sq_phi_star_j = sqrt(phi_s_star_j) ;
threshold = 2*sqrt((log_epsilon - log_eps)/t) ;
sq_phi_star_j1 = min(sqrt(phi_s_star_j1), threshold - sq_phi_star_j) ;
if pj < 1.0e-14 && qj < 1.0e-14
    sq_phibar_star_j = sq_phi_star_j ;
    sq_phibar_star_j1 = sq_phi_star_j1 ;
    adm_region = 1 ;
end
if pj < 1.0e-14 && qj >= 1.0e-14
    sq_phibar_star_j = sq_phi_star_j ;
    if sq_phi_star_j > 0
        f_min = fac*(sq_phi_star_j/(sq_phi_star_j1-sq_phi_star_j))^qj ;
    else
        f_min = fac ;
    end
    if f_min < f_max
        f_bar = f_min + f_min/f_max*(f_max-f_min) ;
        fq = f_bar^(-1/qj) ;
        sq_phibar_star_j1 = (2*sq_phi_star_j1-fq*sq_phi_star_j)/(2+fq) ;
        adm_region = 1 ;
    else
        adm_region = 0 ;
    end
end
if pj >= 1.0e-14 && qj < 1.0e-14
    sq_phibar_star_j1 = sq_phi_star_j1 ;
    f_min = fac*(sq_phi_star_j1/(sq_phi_star_j1-sq_phi_star_j))^pj ;
    if f_min < f_max
        f_bar = f_min + f_min/f_max*(f_max-f_min) ;
        fp = f_bar^(-1/pj) ;
        sq_phibar_star_j = (2*sq_phi_star_j+fp*sq_phi_star_j1)/(2-fp) ;
        adm_region = 1 ;
    else
        adm_region = 0 ;
    end
end
if pj >= 1.0e-14 && qj >= 1.0e-14
    f_min = fac*(sq_phi_star_j+sq_phi_star_j1)/...
        (sq_phi_star_j1-sq_phi_star_j)^max(pj,qj) ;
    if f_min < f_max
        f_min = max(f_min,1.5) ;
        f_bar = f_min + f_min/f_max*(f_max-f_min) ;
        fp = f_bar^(-1/pj) ;
        fq = f_bar^(-1/qj) ;
        if ~conservative_error_analysis
            w = -phi_s_star_j1*t/log_epsilon ;
        else
            w = -2*phi_s_star_j1*t/(log_epsilon-phi_s_star_j1*t) ;
        end
        den = 2+w - (1+w)*fp + fq ;
        sq_phibar_star_j = ((2+w+fq)*sq_phi_star_j + fp*sq_phi_star_j1)/den ;
        sq_phibar_star_j1 = (-(1+w)*fq*sq_phi_star_j ...
            + (2+w-(1+w)*fp)*sq_phi_star_j1)/den ;
        adm_region = 1 ;
    else
        adm_region = 0 ;
    end
end
if adm_region
    log_epsilon = log_epsilon  - log(f_bar) ;
    if ~conservative_error_analysis
        w = -sq_phibar_star_j1^2*t/log_epsilon ;
    else
        w = -2*sq_phibar_star_j1^2*t/(log_epsilon-sq_phibar_star_j1^2*t) ;
    end
    muj = (((1+w)*sq_phibar_star_j + sq_phibar_star_j1)/(2+w))^2 ;
    hj = -2*pi/log_epsilon*(sq_phibar_star_j1-sq_phibar_star_j)...
        /((1+w)*sq_phibar_star_j + sq_phibar_star_j1) ;
    Nj = ceil(sqrt(1-log_epsilon/t/muj)/hj) ;
else
    muj = 0 ; hj = 0 ; Nj = +Inf ;
end
end
function [muj,hj,Nj] = OptimalParam_RU (t, phi_s_star_j, pj, log_epsilon)
sq_phi_s_star_j = sqrt(phi_s_star_j) ;
if phi_s_star_j > 0
    phibar_star_j = phi_s_star_j*1.01 ;
else
    phibar_star_j = 0.01 ;
end
sq_phibar_star_j = sqrt(phibar_star_j) ;
f_min = 1 ; f_max = 10 ; f_tar = 5 ;
stop = 0 ;
while ~stop
    phi_t = phibar_star_j*t ; log_eps_phi_t = log_epsilon/phi_t ;
    Nj = ceil(phi_t/pi*(1 - 3*log_eps_phi_t/2 + sqrt(1-2*log_eps_phi_t))) ;
    A = pi*Nj/phi_t ;
    sq_muj = sq_phibar_star_j*abs(4-A)/abs(7-sqrt(1+12*A)) ;
    fbar = ((sq_phibar_star_j-sq_phi_s_star_j)/sq_muj)^(-pj) ;
    stop = (pj < 1.0e-14) || (f_min < fbar && fbar < f_max) ;
    if ~stop
        sq_phibar_star_j = f_tar^(-1/pj)*sq_muj + sq_phi_s_star_j ;
        phibar_star_j = sq_phibar_star_j^2 ;
    end
end
muj = sq_muj^2 ;
hj = (-3*A - 2 + 2*sqrt(1+12*A))/(4-A)/Nj ;
log_eps = log(eps) ; threshold = (log_epsilon - log_eps)/t ;
if muj > threshold
    if abs(pj) < 1.0e-14 , Q = 0 ; else Q = f_tar^(-1/pj)*sqrt(muj) ; end
    phibar_star_j = (Q + sqrt(phi_s_star_j))^2 ;
    if phibar_star_j < threshold
        w = sqrt(log_eps/(log_eps-log_epsilon)) ;
        u = sqrt(-phibar_star_j*t/log_eps) ;
        muj = threshold ;
        Nj = ceil(w*log_epsilon/2/pi/(u*w-1)) ;
        hj = sqrt(log_eps/(log_eps - log_epsilon))/Nj ;
    else
        Nj = +Inf ; hj = 0 ;
    end
end
end
function R = mldr(t,s,alpha,beta,k)
omega = zeros(1,k+1) ; omega(1) = 1 ;
jj = 0 : k ; p = (alpha-jj) ; pr = 1 ;
for j = 1 : k+1
    pr = pr*p(j) ;
    omega(j+1) = pr/factorial(j) ;
end
Hk = zeros(1, k+1) ; Hk(1) = 1 ;
for j = 1 : k
    ll = 1 : j ;
    Hk(j+1) = -1/alpha*sum(omega(ll+2).*(k*ll/j+1).*Hk(j-ll+1))  ;
end
ck = zeros(1, k+1) ;
for j = 0 : k
    temp = 0 ;
    for l = 0 : k-j
        if l == 0
            p = 1 ;
        else
            ll = 0 : l-1 ;
            p = prod(alpha-beta-ll) ;
        end
        temp = temp + p*Hk(k-j-l+1)/factorial(l) ;
    end
    ck(j+1) = temp/factorial(j) ;
end
R = 1/alpha^(k+1)*exp(t*s).*s.^(1-alpha*k-beta).*polyval(fliplr(ck),s) ;
end
function [E, Err_Round] = mlds(z,al,be,k)
max_gamma_arg = 171.624 ;
Jmax = floor((max_gamma_arg - be)/al) ;
G = gamma(al*(0:Jmax)+be) ;
jj = k : Jmax ;
f = ones(size(jj)) ;
for l = 1 : k
    f = f.*(jj-l+1) ;
end
c = f./G(k+1:Jmax+1) ;
E = zeros(size(z)) ; Err_Round = zeros(size(z)) ;
Err_Round1 = Err_Round ; Err_Round2 = Err_Round ;
for n = 1 : length(z)
    if abs(z(n)) < eps
        E(n) = factorial(k)/G(k+1) ;
    else
        sum_arg = c.*z(n).^(jj-k) ;
        abs_sum_arg = abs(sum_arg) ;
        i_abs_sum_arg = abs_sum_arg > eps/2 ;
        if ~any(i_abs_sum_arg) , i_abs_sum_arg(1) = 1 ; end
        abs_sum_arg = abs_sum_arg(i_abs_sum_arg) ;
        sum_arg = sum_arg(i_abs_sum_arg) ;
        [abs_sum_arg,i_abs_sum_arg] = sort(abs_sum_arg) ;
        sum_arg = sum_arg(i_abs_sum_arg) ;
        if length(sum_arg) == 1
            E(n) = sum_arg ;
        else
            S = cumsum(sum_arg) ; S = S(2:end) ;
            E(n) = S(end) ;
        end
        J = length(sum_arg) - 1 ;
        JJ = [ J , J:-1:1] ;
        Err_Round1(n) = sum(JJ.*abs_sum_arg)*eps ;
        if length(sum_arg) == 1
            Err_Round2(n) = Err_Round1(n) ;
        else
            Err_Round2(n) = sum(abs(S))*eps ;
        end
        Err_Round(n) = exp((log(Err_Round1(n))+log(Err_Round2(n)))/2) ;
    end
    
end
end
