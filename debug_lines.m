%% Parameters and functions
gfuns   =   functions_given;
pfuns   =   functions_plot;
% pfuns.potential(Record_q1{180/Record_dt},bj_calODH,x,y,5,['density. t = ']);
%% Parameters
global Record_dt
global area_cal h boundary_O boundary_D boundary_D_2h boundary_H
global MC_std MC_mean dir_data MC_nseq
% gfuns.Para(20221222);
gfuns.Para(0);

%% Layout settings
% 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
x               =   (area_cal(1,1)+0.5*h):h:(area_cal(1,2)-0.5*h);      nx = length(x);
y               =   (area_cal(2,2)-0.5*h):-h:(area_cal(2,1)+0.5*h);     ny = length(y);
x_ex3           =   (x(1)-3*h):h:(x(end)+3*h);              nx_ex3 = length(x_ex3);
y_ex3           =   (y(1)+3*h):-h:(y(end)-3*h);             ny_ex3 = length(y_ex3);
[bj_ex3ODH]    =   gfuns.Layout(boundary_O,boundary_D,boundary_D_2h,boundary_H,x_ex3,y_ex3);
[bj_calODH]    =   gfuns.Layout(boundary_O,boundary_D,boundary_D_2h,boundary_H,x,y);

mex F90_weno5_r8.f90;

cost = 1.0 ./ max(0,(1.34 .* exp( -0.09*(max(q1_RK1,0.0).^2) ))) + 0.002 .* (max(q1_RK1,0.0)).^2;
potential = fast_WENO3(cost,h,bj_ex3ODH);
potential = fast_Godunov(cost,h,bj_ex3ODH);
[q1_RK2,alpha_LF2,nit,test]  =   F90_weno5_r8(q1_RK1,bj_ex3ODH,F_in,h,tstep);
pfuns.potential(potential,bj_ex3ODH,x_ex3,y_ex3,100,'');


cost = ones(54,104);
potential = fast_WENO3(cost,h,bj_ex3ODH);

[q1_RK1,alpha_LF1,nit,test]  =   F90_weno5(q1,bj_ex2ODH,F_in,h,tstep);

[q1_RK2,alpha_LF2,nit,AAATest]  =   F90_weno5(q1_RK1,bj_ex2ODH,F_in,h,tstep);
pfuns.potential(Record_q1{2}(:,:,13),bj_calODH,x,y,4,'')
contour(x,y,Record_q1{2}(:,:,13)); colorbar;




























