% The code tests the 2nd example in section 6 and generates Figure 6.2.
%
% IMFORMATION - This test took about 14 mins on a laptop equipped with an
% Intel i7-6700HQ processor running at 2.60GHz and with 16GB of RAM.

format compact

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
marker_nd   =  '^:';
marker_mlm  =  'o:';

% colours
color_nd    = [0.8500, 0.3250, 0.0980];
color_mlm   = [0, 0.4470, 0.9410];
% color_cond  = [0.2, 0.8450, 0.9330];
color_cond  = [0, 0, 0];

% line styles

ls_cond = '-';

% linewidth
lw = 1.3; 
lwcond = 1.5;
lw_ca = 1.3;

msize = 10; % marker size
ca_fsize    = 14;
title_fsize = 18;
lgd_fsize   = 14;

ymin = 1e-15;
ymax = 1e-10;

xticks_str = [0:1:10];
yticks_str = [1e-15,1e-14,1e-13,1e-12,1e-11,1e-10];
xmin = 0;
xmax = 10.5;

lcn = 'NorthWest';
yaxis_label = 'Normwise relative errors';


rng(1);

u = eps(1/2);
alpha = 0.8;
beta  = 0.5:0.5:10;
num_beta = length(beta);

% initialization the arrays for storing the errors and conds
error_nd = zeros(2,num_beta);
error_mlm = zeros(2,num_beta);
condu = zeros(2,num_beta);

if produce_results
    main_loop = tic; % record the time consumption
    % testing A_21
    fprintf('Testing A21...\n')
    % generate the test matrix A_{21}
    z = zeros(1,3);
    a = ones(1,6);
    b = 5*a;
    c = -10*ones(1,3);
    D = diag([z a -a b -b c]);
    U = randn(30);
    A21 = U*D/U;
    A = A21;
    i = 1;
    for j = 1:num_beta
        [error_nd(i,j),error_mlm(i,j)] = mlm_comput_error(A,alpha,beta(j));
        condfunc = @(A) mlm_mlm(A,alpha,beta(j));
        condu(i,j) = funm_condest1(A, condfunc)*u;
    end
    % testing A_22
    fprintf('Testing A22...\n')
    % generate the test matrix A_{22}
    a = 0.9 + (1.0-0.9).*rand(1,5);
    b = 1.2 + (1.3-1.2).*rand(1,4);
    c = 0.9 + (1.0-0.9).*rand(1,4) + 1i;
    d = 0.9 + (1.0-0.9).*rand(1,4) - 1i;
    e = 1.4 + (1.5-1.4).*rand(1,3);
    D = diag([a -a b -b c -c d -d e -e]);
    U = randn(40);
    A22 = U*D/U;
    A = A22;
    i = 2;
    for j = 1:num_beta
        [error_nd(i,j),error_mlm(i,j)] = mlm_comput_error(A,alpha,beta(j));
        condfunc = @(A) mlm_mlm(A,alpha,beta(j));
        condu(i,j) = funm_condest1(A, condfunc)*u;
    end    
    fprintf('The test took %.2f minutes.\n', toc(main_loop)/60);
    save('data/exmp2.mat','error_nd','error_mlm','condu');
else
    load('data/exmp2.mat');
end

% generate the figures

% figure for A_{21}
figure
clf;
hold on
plot(beta,error_mlm(1,:),marker_mlm,'MarkerSize',msize,'LineWidth',lw,'color',color_mlm);
plot(beta,error_nd(1,:),marker_nd,'MarkerSize',msize,'LineWidth',lw,'color',color_nd);
plot(beta,condu(1,:),ls_cond,'LineWidth',lwcond,'color',color_cond);
hold off

set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('mlm', 'funm\_nd', 'Location', lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([xmin,xmax]);
ylim([ymin,ymax]);
xticks(xticks_str);
yticks(yticks_str);
xlabel('\beta');
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$A_{21}$','interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/exmp2_A21.pdf');
    export_fig(str)
end

% figure for A_{22}
figure
clf;
hold on
plot(beta,error_mlm(2,:),marker_mlm,'MarkerSize',msize,'LineWidth',lw,'color',color_mlm);
plot(beta,error_nd(2,:),marker_nd,'MarkerSize',msize,'LineWidth',lw,'color',color_nd);
plot(beta,condu(2,:),ls_cond,'LineWidth',lwcond,'color',color_cond);
hold off

set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('mlm', 'funm\_nd', 'Location', lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([xmin,xmax]);
ylim([ymin,ymax]);
xticks(xticks_str);
yticks(yticks_str);
xlabel('\beta');
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$A_{22}$','interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/exmp2_A22.pdf');
    export_fig(str)
end