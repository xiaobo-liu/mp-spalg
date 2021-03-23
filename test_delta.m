% The code generates Figure 5.1: comparing the errors for algorithms with
% different blocking parameter delta. 

addpath('data');
addpath('external');
addpath('include');
warning('off');

% reproduce the results?
produce_results = true;

% save the plots?
save_results = true;

% initialization

% marker styles
marker_01   =  '^';
marker_02   =  'x';
marker_norm =  'v';
marker_inf  =  'o';

% colours
color_01    = [0.8500, 0.3250, 0.0980];
color_02    = [0, 0.5, 0];
color_norm  = [0.4940, 0.1840, 0.5560]; 
color_inf   = [0, 0.4470, 0.9410];
% color_cond  = [0.2, 0.8450, 0.9330];
color_cond  = [0, 0, 0];

ls_cond = '-';
% linewidth
lw = 2.5; 
lwcond = 2.7;
lw_ca = 2.5;

% marker size
msize = 22;
ca_fsize    = 28;
title_fsize = 37;
lgd_fsize   = 28;

lcn = 'NorthEast';
lcn_small = 'NorthWest'; % location of the legend in the plot for A/1e2
yaxis_label = '';

rng(1);
fun = @sin; if isequal(fun,@sin), refm = @sinm; end

n = 32; 
num_mat = testmats_mct(0);

data = zeros(num_mat,5);
data_small = zeros(num_mat,5);
data_big = zeros(num_mat,5);
I = eye(n);

main_loop = tic; % record the time consumption

if produce_results
    for k = 1:num_mat
        % A from the test set
        A = testmats_mct(k,n);
        old_digits = mp.Digits();
        mp.Digits(200);
        refX = refm(mp(A));
        mp.Digits(old_digits);
        refX_norm = double(norm(refX,'fro'));
        X2 = funm_nd(A,fun);   
        X3 = funm_nd(A,fun,0.2);   
        X4 = funm_nd(A,fun,'norm');   
        X5 = funm_nd_inf(A,fun);
        data(k,2) = norm(X2-refX,'fro')/refX_norm;
        data(k,3) = norm(X3-refX,'fro')/refX_norm;
        data(k,4) = norm(X4-refX,'fro')/refX_norm;
        data(k,5) = norm(X5-refX,'fro')/refX_norm;
        % compute cond of sin via the chain rule 
        data(k,1) = funm_condest1(A-pi/2*I, @cosm)*eps(1/2);
    
        % A*1e-2
        A = testmats_mct(k,n)/1e2;
        mp.Digits(200);
        refX = refm(mp(A));
        mp.Digits(old_digits);
        refX_norm = double(norm(refX,'fro'));
        X2 = funm_nd(A,fun);   
        X3 = funm_nd(A,fun,0.2);   
        X4 = funm_nd(A,fun,'norm');   
        X5 = funm_nd_inf(A,fun);
        data_small(k,2) = norm(X2-refX,'fro')/refX_norm;
        data_small(k,3) = norm(X3-refX,'fro')/refX_norm;
        data_small(k,4) = norm(X4-refX,'fro')/refX_norm;
        data_small(k,5) = norm(X5-refX,'fro')/refX_norm;
        % compute condf
        data_small(k,1) = funm_condest1(A-pi/2*I, @cosm)*eps(1/2);
    
        % figure 3: A*1e2
        A = testmats_mct(k,n)*1e2;
        mp.Digits(200);
        refX = refm(mp(A));
        mp.Digits(old_digits);
        refX_norm = double(norm(refX,'fro'));
        X2 = funm_nd(A,fun);   
        X3 = funm_nd(A,fun,0.2);   
        X4 = funm_nd(A,fun,'norm');   
        X5 = funm_nd_inf(A,fun);
        data_big(k,2) = norm(X2-refX,'fro')/refX_norm;
        data_big(k,3) = norm(X3-refX,'fro')/refX_norm;
        data_big(k,4) = norm(X4-refX,'fro')/refX_norm;
        data_big(k,5) = norm(X5-refX,'fro')/refX_norm;
        % compute condf
        data_big(k,1) = funm_condest1(A-pi/2*I, @cosm)*eps(1/2);
        k
    end
    fprintf('Producing the results took %.2f minutes.\n', toc(main_loop)/60);
    data = sortrows(data,-1);
    data_small = sortrows(data_small,-1);
    data_big = sortrows(data_big,-1);
    save('data/test_dalta.mat','data', 'data_small', 'data_big');
else
    load ('data/test_dalta.mat');
end
condu = data(:,1);
err_01 = data(:,2); 
err_02 = data(:,3);
err_norm = data(:,4);
err_inf = data(:,5);
condu_small = data_small(:,1);
err_01_small = data_small(:,2);
err_02_small = data_small(:,3);
err_norm_small = data_small(:,4);
err_inf_small = data_small(:,5);
condu_big = data_big(:,1);
err_01_big = data_big(:,2);
err_02_big = data_big(:,3);
err_norm_big = data_big(:,4);
err_inf_big = data_big(:,5);
    


% map any error outside [ymin ymax] back to the range
ymin = 1e-18;
ymax = 1e0;

err_01(err_01 > ymax) = ymax;
err_02(err_02 > ymax) = ymax;
err_norm(err_norm > ymax) = ymax;
err_inf(err_inf > ymax) = ymax;
err_01_small(err_01_small > ymax) = ymax;
err_02_small(err_02_small > ymax) = ymax;
err_norm_small(err_norm_small > ymax) = ymax;
err_inf_small(err_inf_small > ymax) = ymax;
err_01_big(err_01_big > ymax) = ymax;
err_02_big(err_02_big > ymax) = ymax;
err_norm_big(err_norm_big > ymax) = ymax;
err_inf_big(err_inf_big > ymax) = ymax;

err_01(err_01 < ymin) = ymin;
err_02(err_02 < ymin) = ymin;
err_norm(err_norm < ymin) = ymin;
err_inf(err_inf < ymin) = ymin;
err_01_small(err_01_small < ymin) = ymin;
err_02_small(err_02_small < ymin) = ymin;
err_norm_small(err_norm_small < ymin) = ymin;
err_inf_small(err_inf_small < ymin) = ymin;
err_01_big(err_01_big < ymin) = ymin;
err_02_big(err_02_big < ymin) = ymin;
err_norm_big(err_norm_big < ymin) = ymin;
err_inf_big(err_inf_big < ymin) = ymin;

% generate the figures
figure
clf;
hold on
plot(1:num_mat,err_01,marker_01,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_01);
plot(1:num_mat,err_02,marker_02,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_02);
plot(1:num_mat,err_norm,marker_norm,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_norm);
plot(1:num_mat,err_inf,marker_inf,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_inf);
plot(1:num_mat,condu,ls_cond,'LineWidth',lwcond,...
    'color',color_cond);
hold off
set(gcf,'outerposition',get(0,'ScreenSize'));
set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('funm\_nd\_0.1', 'funm\_nd\_0.2', 'funm\_nd\_norm',...
    'funm\_nd\_\infty','Location',lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([0,num_mat+1]);
ylim([ymin,ymax]);
yticks([1e-16,1e-12,1e-8,1e-4,1e0]);
title('$A$ from the test set','interpreter','latex','FontWeight','normal','fontsize',title_fsize);
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/test_delta.pdf');
    export_fig(str)
end

figure;
clf;
hold on
plot(1:num_mat,err_01_small,marker_01,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_01);
plot(1:num_mat,err_02_small,marker_02,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_02);
plot(1:num_mat,err_norm_small,marker_norm,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_norm);
plot(1:num_mat,err_inf_small,marker_inf,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_inf);
plot(1:num_mat,condu_small,ls_cond,'LineWidth',lwcond,...
    'color',color_cond);
hold off
set(gcf,'outerposition',get(0,'ScreenSize'));
set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('funm\_nd\_0.1', 'funm\_nd\_0.2', 'funm\_nd\_norm',...
    'funm\_nd\_\infty','Location',lcn_small);
lgd.FontSize = lgd_fsize;
box on
xlim([0,num_mat+1]);
ylim([ymin,ymax]);
yticks([1e-16,1e-12,1e-8,1e-4,1e0]);
title('$A*$1e-2','interpreter','latex','FontWeight','normal','fontsize',title_fsize)
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/test_delta_small.pdf');
    export_fig(str) 
end

figure;
clf;
hold on
plot(1:num_mat,err_01_big,marker_01,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_01);
plot(1:num_mat,err_02_big,marker_02,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_02);
plot(1:num_mat,err_norm_big,marker_norm,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_norm);
plot(1:num_mat,err_inf_big,marker_inf,'MarkerSize',msize,'LineWidth',lw,...
    'color',color_inf);
plot(1:num_mat,condu_big,ls_cond,'LineWidth',lwcond,...
    'color',color_cond);
hold off
set(gcf,'outerposition',get(0,'ScreenSize'));
set(gca,'linewidth',lw_ca)
set(gca,'fontsize',ca_fsize)
set(gca, 'YScale', 'log')
lgd = legend('funm\_nd\_0.1', 'funm\_nd\_0.2', 'funm\_nd\_norm',...
    'funm\_nd\_\infty','Location',lcn);
lgd.FontSize = lgd_fsize;
box on
xlim([0,num_mat+1]);
ylim([ymin,ymax]);
yticks([1e-16,1e-12,1e-8,1e-4,1e0]);
title('$A*$1e2','interpreter','latex','FontWeight','normal','fontsize',title_fsize)
% ylabel(yaxis_label,'interpreter','latex','FontWeight','normal','fontsize',title_fsize)
if save_results
    set(gcf, 'Color', 'w')
    str = ('../figs/test_delta_big.pdf');
    export_fig(str)
end