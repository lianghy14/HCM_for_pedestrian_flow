function [ve_x,ve_y] = RHS_FSM(Q_cell)
global area x y nx ny
global tStart tEnd
global plot_config plot_dt dir_fig
global f1 t_bound h x_max y_max boundary_H boundary_O boundary_D boundary_D_2h boundary_D_outside
global m_ave v_free gama_1 CFL sigma epsilon
global gfuns pfuns

% Obtain cost derived by Greenshield's model and discomfort function
% Obtain magnitude of velosity
boundary_judge = ones(ny,nx);
method = 'zero';% 'linear'/'nearest'/'pchip'/'spline'
density_ex = gfuns.Boundary_ex(Q_cell{1},boundary_judge,'x',2,method);
x_ex = (x(1)-2*h):h:(x(end)+2*h); nx_ex = length(x_ex);
velosity_x_ex = gfuns.Velosity(density_ex);
cost_1 = 1./velosity_x_ex; cost_1(isnan(cost_1)) = 10^12;
cost_1 = gfuns.Boundary_value(x_ex,y,cost_1,boundary_H,10^12);
cost_2 = 0;
cost =  cost_1 + cost_2;
% Obtain potential by fast sweeping WENO using Gauss-Seidel iterations
potential = fast_Godunov(cost,boundary_D,x_ex,y);
boundary_judge = ones(ny,nx_ex);
method = 'linear';% 'linear'/'nearest'/'pchip'/'spline'
boundary_judge_H = gfuns.Boundary_value(x_ex,y,boundary_judge,boundary_H,0);
boundary_judge_DH = gfuns.Boundary_value(x_ex,y,boundary_judge_H,boundary_D_outside,0);
potential_x_ex = gfuns.Boundary_ex(potential,boundary_judge_DH,'x',1,method);
potential_y_ex = gfuns.Boundary_ex(potential,boundary_judge_DH,'y',1,method);
potential_x = -(potential_x_ex(:,3:end)-potential_x_ex(:,1:end-2))./2./h;
potential_y = -(potential_y_ex(1:end-2,:)-potential_y_ex(3:end,:))./2./h;
phy_grad = sqrt(potential_x.^2+potential_y.^2);
ve_x = potential_x(:,3:end-2)./phy_grad(:,3:end-2);
ve_y = potential_y(:,3:end-2)./phy_grad(:,3:end-2);
ve_x(isnan(ve_x)) = 0;
ve_y(isnan(ve_y)) = 0;

end