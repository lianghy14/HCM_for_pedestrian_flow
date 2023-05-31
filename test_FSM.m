%% Test fast sweeping method
clc
clear
global iteration
global area_cal x y nx ny
global x_ex2 y_ex2 nx_ex2 ny_ex2
global bj_calODH bj_ex2ODH

global bj_ex2 bj_ex2H bj_ex2D bj_ex2HD bj_ex2HDh bj_ex2O bj_ex2HDhO bj_ex2HDO bj_ex2HO
global bj_ex2x bj_ex2xD bj_ex2y bj_ex2yD 
global bj_midxH bj_midyH bj_midxO bj_midyO
global tStart tEnd method 
global plot_config plot_dt dir_fig
global f1 t_bound h x_max y_max boundary_H boundary_O boundary_D boundary_D_2h boundary_D_outside
global m_ave v_free gama_1 CFL sigma epsilon X_lagrange
global gfuns pfuns
global fast_IC

gfuns = functions_given;
pfuns = functions_plot;
gfuns.Para();
x               =   (area_cal(1,1)+0.5*h):h:(area_cal(1,2)-0.5*h);      nx = length(x);
y               =   (area_cal(2,2)-0.5*h):-h:(area_cal(2,1)+0.5*h);     ny = length(y);
x_ex2           =   (x(1)-2*h):h:(x(end)+2*h);              nx_ex2 = length(x_ex2);
y_ex2           =   (y(1)+2*h):-h:(y(end)-2*h);             ny_ex2 = length(y_ex2);
[bj_ex2ODH]    =   gfuns.Layout(boundary_O,boundary_D,boundary_D_2h,boundary_H,x_ex2,y_ex2);
[bj_calODH]    =   gfuns.Layout(boundary_O,boundary_D,boundary_D_2h,boundary_H,x,y);

% Initialization
fast_IC = 10^12 .* ones(ny_ex2,nx_ex2);
cost = zeros(ny_ex2,nx_ex2);
potential_t = zeros(ny_ex2,nx_ex2);

for i = 1:nx_ex2
    for j = 1:ny_ex2
        cost(j,i) = (pi()./2)* ( (sin( pi()+pi()/2*(x_ex2(i)-1) )).^2 + (sin( pi()+pi./2.*(y_ex2(j)-1) ))^2 )^(1/2);
        potential_t(j,i) = cos( pi()+pi/2*(x_ex2(i)-1) ) + cos( pi()+pi/2*(y_ex2(j)-1) );
    end
end
% 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
potential_t(bj_ex2ODH==0) = 10^12; potential_t(bj_ex2ODH==2) = 10^12;
fast_IC(bj_ex2ODH==1) = potential_t(bj_ex2ODH==1);


% cost = ones(50,50); bj_ex2ODH = zeros(50,50); bj_ex2ODH(3:48,3:48) = 1; h = 1; [potential, p_x, p_y, nit] = F90_fsmgod(cost,bj_ex2ODH,h);
potential = fast_WENO3(cost,h,bj_ex2ODH);
potential(bj_ex2ODH<=2) = potential_t(bj_ex2ODH<=2);
Error_L2 = sqrt(mean(mean((potential - potential_t).^2)));
Error_L2
% potential = potential(3:end-2,3:end-2);pfuns.potential(potential,'potential');
% potential_t = potential_t(3:end-2,3:end-2);pfuns.potential(potential_t,'potential');


