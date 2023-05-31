clc
clear

%% Parameters and functions
gfuns = functions_given;
pfuns = functions_plot;

%% Parameters
global Record_dt
global area_cal h boundary_O boundary_D boundary_D_2h boundary_H
global MC_std MC_mean dir_data MC_nseq MC_max
gfuns.Para(20221222);

%% Layout settings
% 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
x               =   (area_cal(1,1)+0.5*h):h:(area_cal(1,2)-0.5*h);      nx = length(x);
y               =   (area_cal(2,2)-0.5*h):-h:(area_cal(2,1)+0.5*h);     ny = length(y);
x_ex3           =   (x(1)-3*h):h:(x(end)+3*h);              nx_ex3 = length(x_ex3);
y_ex3           =   (y(1)+3*h):-h:(y(end)-3*h);             ny_ex3 = length(y_ex3);
[bj_ex3ODH]    =   gfuns.Layout(boundary_O,boundary_D,boundary_D_2h,boundary_H,x_ex3,y_ex3);
[bj_calODH]    =   gfuns.Layout(boundary_O,boundary_D,boundary_D_2h,boundary_H,x,y);
cpu_tstart  =   tic;

%% Load data

load([dir_data 'MC' num2str(MC_max) '.mat']);
len_y = size(Record_q1{1,2},1);
len_x = size(Record_q1{1,2},2);
len_t = size(Record_q1{1,2},3);

temp = zeros(len_y,len_x,len_t,MC_max);
for i = 1:MC_max
    temp(:,:,:,i) = Record_q1{i,2};
end
MC_max_mean = mean(temp,4);
MC_max_std = std(temp,0,4);

MC_seq_mean = zeros(len_y,len_x,len_t,length(MC_nseq));
MC_seq_std = zeros(len_y,len_x,len_t,length(MC_nseq));
MC_f_mean = zeros(1,length(MC_nseq));
MC_f_std = zeros(1,length(MC_nseq));

for i = 1:length(MC_nseq)
    load([dir_data 'MC' num2str(MC_nseq(i)) '.mat']);
    temp = zeros(len_y,len_x,len_t,MC_nseq(i));
    for j = 1:MC_nseq(i)
        temp(:,:,:,j) = Record_q1{j,2};
    end
    MC_seq_mean(:,:,:,i) = mean(temp,4);
    MC_f_mean(i) = mean(MC_fin);
    MC_seq_std(:,:,:,i) = std(temp,0,4);
    MC_f_std(i) = std(MC_fin);
    cpu_dt = toc(cpu_tstart);
    fprintf('MC. N = %d. CPUt is %.2f\n',MC_nseq(i),cpu_dt);
end
%% RRMSE
RRMSE_mean = zeros(1,length(MC_nseq));
RRMSE_std = zeros(1,length(MC_nseq));
MC_t = 1:25;
RRMSE_valid = (sum(sum(bj_calODH~=0)))*length(MC_t);

for i = 1:length(MC_nseq)
    RRMSE_mean(i) = sqrt( gfuns.RRMSE( (MC_seq_mean(:,:,MC_t,i) - MC_max_mean(:,:,MC_t)),2,RRMSE_valid) )./gfuns.RRMSE(MC_max_mean(:,:,MC_t),1,RRMSE_valid);
    RRMSE_std(i) = sqrt( gfuns.RRMSE( (MC_seq_std(:,:,MC_t,i) - MC_max_std(:,:,MC_t)),2,RRMSE_valid) )./gfuns.RRMSE(MC_max_std(:,:,MC_t),1,RRMSE_valid);
end

fig = figure();
plot(MC_nseq./100,RRMSE_mean.*100,'.-k'); hold on;
plot(MC_nseq./100,RRMSE_std.*100,'.-r'); hold off;
ylabel('RRMSE (%)');
xlabel('N / 100');
legend('AVE','STD')
xlim([1 20]);
% ylim([0 8]); 
set(fig,'Units','centimeters','Position',[10 10 20 8]);
set(gca,'TickLabelInterpreter','latex');
set(gca, 'LooseInset', [0,0,0.02,0.02]);

fig = figure();
plot(MC_nseq./100,MC_f_mean,'.-k'); hold on;
plot(MC_nseq./100,MC_f_std,'.-r'); hold off;
ylabel('Inflow');
xlabel('N / 100');
legend('AVE','STD')
xlim([1 20]);
% ylim([0 8]); 
set(fig,'Units','centimeters','Position',[10 10 20 8]);
set(gca,'TickLabelInterpreter','latex');
set(gca, 'LooseInset', [0,0,0.02,0.02]);


