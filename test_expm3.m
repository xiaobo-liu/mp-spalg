% The code tests the 3rd example in section 6 and generates Figure 6.3.
%
% IMFORMATION - This test took about 5 mins on a laptop equipped with an
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
marker_nd   =  '^';
marker_mlm  =  'o';

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


ymin = 1e-18;
ymax = 1e0;
xmin = 0;
xmax = 33;
% xticks_str = [0 4 8 12 16 20 24 28 32];
yticks_str = [1e-17 1e-13 1e-9 1e-5 1e0];


lcn = 'NorthEast';
yaxis_label = 'Normwise relative errors';

% generate the test matrix A_{21}
rng(1);

num_mat = mlm_testmats_mct(0); 
n = 10;

u = eps(1/2);
alpha = 0.8;
beta  = [1.2 8.0];
num_beta = length(beta);

% initialization the arrays for storing the errors and conds
error_nd = zeros(num_mat,num_beta);
error_mlm = zeros(num_mat,num_beta);
condu = zeros(num_mat,num_beta);

if produce_results
    main_loop = tic; % record the time consumption
    for k = 1:num_mat
        for j = 1:num_beta
            A = mlm_testmats_mct(k,n);
            [error_nd(k,j),error_mlm(k,j)] = mlm_comput_error(A,alpha,beta(j));
            condfunc = @(A) mlm_mlm(A,alpha,beta(j));
            condu(k,j) = funm_condest1(A, condfunc)*u;
        end
        fprintf('Testing the %.2dth matrix...\n', k);
    end
    for j=1:num_beta % stores condition numbers
        [condu(:,j),order] = sort(condu(:,j),'descend');
        error_nd(:,j)= error_nd((order),j);
        error_mlm(:,j)= error_mlm((order),j);
    end
    fprintf('The test took %.2f minutes.\n', toc(main_loop)/60);
    save('data/exmp3.mat','num_mat','error_nd','error_mlm','condu');
else
    load('data/exmp3.mat');
end

% map any error less than 1e-17 to 1e-17
miny = 1e-17; 
for j=1:num_beta
    error_nd((error_nd(:,j)<miny),j) = miny;
    error_mlm((error_mlm(:,j)<miny),j) = miny;
end

% figure for beta = 1.2
figure
clf;
hold on
plot(1:num_mat,error_mlm(:,1),marker_mlm,'MarkerSize',msize,'LineWidth',lw,'color',color_mlm);
plot(1:num_mat,error_nd(:,1),marker_nd,'MarkerSize',msize,'LineWidth',lw,'color',color_nd);
plot(1:num_mat,condu(:,1),ls_cond,'LineWidth',lwcond,'color',color_cond);
hold off

set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('mlm', 'funm\_nd', 'Location', lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([xmin,xmax]);
ylim([ymin,ymax]);
% xticks(xticks_str);
yticks(yticks_str);
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$\alpha=0.8$, $\beta=1.2$','interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/exmp3_beta12.pdf');
    export_fig(str)
end

% figure for beta = 8.0
figure
clf;
hold on
plot(1:num_mat,error_mlm(:,2),marker_mlm,'MarkerSize',msize,'LineWidth',lw,'color',color_mlm);
plot(1:num_mat,error_nd(:,2),marker_nd,'MarkerSize',msize,'LineWidth',lw,'color',color_nd);
plot(1:num_mat,condu(:,2),ls_cond,'LineWidth',lwcond,'color',color_cond);
hold off

set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('mlm', 'funm\_nd', 'Location', lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([xmin,xmax]);
ylim([ymin,ymax]);
% xticks(xticks_str);
yticks(yticks_str);
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
title('$\alpha=0.8$, $\beta=8.0$','interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/exmp3_beta80.pdf');
    export_fig(str)
end