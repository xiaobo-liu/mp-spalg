% The code generates Figure 5.2, where we are comparing the three 
% Schur--Parlett algorithms for different functions f on matrix A from a
% set of 35 matrices.
addpath('data');
addpath('external');
addpath('include');

warning('off');

% reproduce the results?
produce_results = false;

% save the plots?
save_results = true;

% initialization

% marker styles
marker_01   =  '^';
marker_inf  =  'o';
marker_funm  = '+';

% colours
color_01    = [0.8500, 0.3250, 0.0980];
color_inf   = [0, 0.4470, 0.9410];
color_funm  = [0, 0, 0];
% color_cond  = [0.2, 0.8450, 0.9330];
color_cond  = [0, 0, 0];

ls_cond = '-';
% linewidth
lw = 1.3; 
lwcond = 1.5;
lw_ca = 1.3;


msize = 10; % marker size
ca_fsize    = 14;
title_fsize = 18;
lgd_fsize   = 14;

lcn = 'NorthEast';
yaxis_label = 'Normwise relative errors';

rng(1);

n = 32; 
num_mat = testmats_mct(0);

data_sin = zeros(num_mat,4);
data_cos = zeros(num_mat,4);
data_sinh = zeros(num_mat,4);
data_cosh = zeros(num_mat,4);
I = eye(n);

main_loop = tic; % record the time consumption


if produce_results
    for k = 1:num_mat
        % f = sin
        fun = @sin; refm = @sinm;
        A = testmats_mct(k,n);
        old_digits = mp.Digits();
        mp.Digits(200);
        refX = refm(mp(A));
        mp.Digits(old_digits);
        refX_norm = double(norm(refX,'fro'));
        X2 = funm_nd(A,fun);   
        X3 = funm_nd_inf(A,fun);
        X4 = funm(A,fun);
        data_sin(k,2) = norm(X2-refX,'fro')/refX_norm;
        data_sin(k,3) = norm(X3-refX,'fro')/refX_norm;
        data_sin(k,4) = norm(X4-refX,'fro')/refX_norm;
        % compute cond of sin via the chain rule 
        data_sin(k,1) = funm_condest1(A-pi/2*I, @cosm)*eps(1/2);
    
        % f = cos
        fun = @cos; refm = @cosm;
        A = testmats_mct(k,n);
        old_digits = mp.Digits();
        mp.Digits(200);
        refX = refm(mp(A));
        mp.Digits(old_digits);
        refX_norm = double(norm(refX,'fro'));
        X2 = funm_nd(A,fun);   
        X3 = funm_nd_inf(A,fun);
        X4 = funm(A,fun);
        data_cos(k,2) = norm(X2-refX,'fro')/refX_norm;
        data_cos(k,3) = norm(X3-refX,'fro')/refX_norm;
        data_cos(k,4) = norm(X4-refX,'fro')/refX_norm;
        data_cos(k,1) = funm_condest1(A, @cosm)*eps(1/2);
        
        % f = sinh
        fun = @sinh; refm = @sinhm;
        sinhm_condest = @(A) (expm(A)-expm(-A)) / 2; % used in estimating cond
        A = testmats_mct(k,n);
        old_digits = mp.Digits();
        mp.Digits(200);
        refX = refm(mp(A));
        mp.Digits(old_digits);
        refX_norm = double(norm(refX,'fro'));
        X2 = funm_nd(A,fun);   
        X3 = funm_nd_inf(A,fun);
        X4 = funm(A,fun);
        data_sinh(k,2) = norm(X2-refX,'fro')/refX_norm;
        data_sinh(k,3) = norm(X3-refX,'fro')/refX_norm;
        data_sinh(k,4) = norm(X4-refX,'fro')/refX_norm;
        data_sinh(k,1) = funm_condest1(A, sinhm_condest)*eps(1/2);
        
        % f = cosh
        fun = @cosh; refm = @coshm;
        coshm_condest = @(A) (expm(A)+expm(-A)) / 2; % used in estimating cond
        A = testmats_mct(k,n);
        old_digits = mp.Digits();
        mp.Digits(200);
        refX = refm(mp(A));
        mp.Digits(old_digits);
        refX_norm = double(norm(refX,'fro'));
        X2 = funm_nd(A,fun);   
        X3 = funm_nd_inf(A,fun);
        X4 = funm(A,fun);
        data_cosh(k,2) = norm(X2-refX,'fro')/refX_norm;
        data_cosh(k,3) = norm(X3-refX,'fro')/refX_norm;
        data_cosh(k,4) = norm(X4-refX,'fro')/refX_norm;
        data_cosh(k,1) = funm_condest1(A, coshm_condest)*eps(1/2);
        k
    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
    data_sin = sortrows(data_sin,-1);
    data_cos = sortrows(data_cos,-1);
    data_sinh = sortrows(data_sinh,-1);
    data_cosh = sortrows(data_cosh,-1);
    save('data/test_fullmat.mat','data_sin','data_cos','data_sinh','data_cosh');
else
    load('data/test_fullmat.mat');
end

ymin = 1e-18;
ymax = 1e0;

% generate the figures

% #1: f = sin
figure
clf;
hold on
plot(1:num_mat,data_sin(:,4),marker_funm,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_funm);
plot(1:num_mat,data_sin(:,2),marker_01,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_01);
plot(1:num_mat,data_sin(:,3),marker_inf,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_inf);
plot(1:num_mat,data_sin(:,1),ls_cond,'LineWidth',lwcond,'color',color_cond);
hold off
set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('funm', 'funm\_nd', 'funm\_nd\_\infty', 'Location', lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([0,num_mat+1]);
ylim([ymin,ymax]);
xticks([0 5 10 15 20 25 30 35]);
yticks([1e-16,1e-12,1e-8,1e-4,1e0]);
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$f=\sin$','interpreter','latex','FontWeight','normal','fontsize',title_fsize);
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/test_fullmat_sin.pdf');
    export_fig(str)
end

% #2: f = cos
figure
clf;
hold on
plot(1:num_mat,data_cos(:,4),marker_funm,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_funm);
plot(1:num_mat,data_cos(:,2),marker_01,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_01);
plot(1:num_mat,data_cos(:,3),marker_inf,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_inf);
plot(1:num_mat,data_cos(:,1),ls_cond,'LineWidth',lwcond,'color',color_cond);
hold off
set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('funm', 'funm\_nd', 'funm\_nd\_\infty', 'Location', lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([0,num_mat+1]);
ylim([ymin,ymax]);
xticks([0 5 10 15 20 25 30 35]);
yticks([1e-16,1e-12,1e-8,1e-4,1e0]);
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$f=\cos$','interpreter','latex','FontWeight','normal','fontsize',title_fsize);
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/test_fullmat_cos.pdf');
    export_fig(str)
end

% #3: f = sinh
figure
clf;
hold on
plot(1:num_mat,data_sinh(:,4),marker_funm,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_funm);
plot(1:num_mat,data_sinh(:,2),marker_01,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_01);
plot(1:num_mat,data_sinh(:,3),marker_inf,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_inf);
plot(1:num_mat,data_sinh(:,1),ls_cond,'LineWidth',lwcond,'color',color_cond);
hold off
set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('funm', 'funm\_nd', 'funm\_nd\_\infty', 'Location', lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([0,num_mat+1]);
ylim([ymin,ymax]);
xticks([0 5 10 15 20 25 30 35]);
yticks([1e-16,1e-12,1e-8,1e-4,1e0]);
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$f=\sinh$','interpreter','latex','FontWeight','normal','fontsize',title_fsize);
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/test_fullmat_sinh.pdf');
    export_fig(str)
end

% #4: f = cosh
figure
clf;
hold on
plot(1:num_mat,data_cosh(:,4),marker_funm,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_funm);
plot(1:num_mat,data_cosh(:,2),marker_01,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_01);
plot(1:num_mat,data_cosh(:,3),marker_inf,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_inf);
plot(1:num_mat,data_cosh(:,1),ls_cond,'LineWidth',lwcond,'color',color_cond);
hold off
set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('funm', 'funm\_nd', 'funm\_nd\_\infty', 'Location', lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([0,num_mat+1]);
ylim([ymin,ymax]);
xticks([0 5 10 15 20 25 30 35]);
yticks([1e-16,1e-12,1e-8,1e-4,1e0]);
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$f=\cosh$','interpreter','latex','FontWeight','normal','fontsize',title_fsize);
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/test_fullmat_cosh.pdf');
    export_fig(str)
end
