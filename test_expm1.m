% The code tests the 1st example in section 6 and generates Figure 6.1.
%
% IMFORMATION - This test took about 13 mins on a laptop equipped with an
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
xmin = 0;
xmax = 10.5;

lcn = 'NorthWest';
yaxis_label = 'Normwise relative errors';

rng(1);

u = eps(1/2);

A = -gallery('redheff',20);
alpha = [0.5,0.8];
beta = 0.5:0.5:10;
num_alpha = length(alpha);
num_beta = length(beta);

% initialization the arrays for storing the errors and conds
error_nd = zeros(num_alpha,num_beta);
error_mlm = zeros(num_alpha,num_beta);
condu = zeros(num_alpha,num_beta);

if produce_results
    main_loop = tic; % record the time consumption
    for i =1:num_alpha
        for j = 1:num_beta
            [error_nd(i,j),error_mlm(i,j)] = mlm_comput_error(A,alpha(i),beta(j));
            condfunc = @(A) mlm_mlm(A,alpha(i),beta(j));
            condu(i,j) = funm_condest1(A, condfunc)*u;
        end
    end
    fprintf('The test took %.2f minutes.\n', toc(main_loop)/60);
    save('data/exmp1.mat','error_nd','error_mlm','condu');
else
    load('data/exmp1.mat');
end

% generate the figures

% alpha = 0.5
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
xticks([0:1:10]);
yticks([1e-15,1e-14,1e-13,1e-12,1e-11,1e-10]);
xlabel('\beta');
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$\alpha=0.5$','interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/exmp1_alpha05.pdf');
    export_fig(str)
end

% alpha = 0.8
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
xticks([1:1:10]);
yticks([1e-15,1e-14,1e-13,1e-12,1e-11,1e-10]);
xlabel('\beta');
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$\alpha=0.8$','interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/exmp1_alpha08.pdf');
    export_fig(str)
end