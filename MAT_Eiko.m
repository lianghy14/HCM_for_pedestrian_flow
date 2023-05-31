clc
clear

%% Parameters and functions
gfuns   =   functions_given;
pfuns   =   functions_plot;

%% Parameters
global Record_dt
global area_cal h boundary_O boundary_D boundary_D_2h boundary_H
global MC_std MC_mean dir_data MC_nseq

n0 = 1;
gfuns.Para(20230116);

%% Layout settings
% 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
x               =   (area_cal(1,1)+0.5*h):h:(area_cal(1,2)-0.5*h);      nx = length(x);
y               =   (area_cal(2,2)-0.5*h):-h:(area_cal(2,1)+0.5*h);     ny = length(y);
x_ex3           =   (x(1)-3*h):h:(x(end)+3*h);              nx_ex3 = length(x_ex3);
y_ex3           =   (y(1)+3*h):-h:(y(end)-3*h);             ny_ex3 = length(y_ex3);
[bj_ex3ODH]    =   gfuns.Layout(boundary_O,boundary_D,boundary_D_2h,boundary_H,x_ex3,y_ex3);
[bj_calODH]    =   gfuns.Layout(boundary_O,boundary_D,boundary_D_2h,boundary_H,x,y);
cpu_tstart  =   tic;

%% MC simulations
mkdir(dir_data);
h_par = h;

if n0>=2
    load([dir_data 'MC' num2str(MC_nseq(n0-1)) '.mat']);
end

for n = n0:length(MC_nseq)
    MC_n = MC_nseq(n);
    if n == 1
        MC_np = 0;
        MC_fin = randn(MC_n,1).*MC_std+MC_mean;
        Record_q1 = cell(MC_n,2);
    else
        MC_np = MC_nseq(n-1);
        temp = MC_fin;
        MC_fin = zeros(MC_n,1);
        MC_fin(1:MC_np) = temp;
        MC_fin((MC_np+1):end) = randn(MC_n-MC_np,1).*MC_std+MC_mean;
        temp = Record_q1;
        Record_q1 = cell(MC_n,2);
        Record_q1(1:MC_np,:) = temp;
    end
    for i = (MC_np+1):MC_n
        cost = 1./MC_fin(i).*ones(ny_ex3,nx_ex3);
        Record_q1(i,:) = {MC_fin(i),fast_Godunov(cost,h_par,bj_ex3ODH)};
        cpu_dt = toc(cpu_tstart);
        if mod(i,10)==0
            fprintf('MC. N = %d. i = %d. F = %.2f. CPUt is %.2f\n',MC_n,i,MC_fin(i),cpu_dt);
        end
    end
    Record_name = [dir_data 'MC' num2str(MC_n) '.mat'];
    save(Record_name,'Record*','MC_fin');
end