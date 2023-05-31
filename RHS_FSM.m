function [ve_x,ve_y] = RHS_FSM(Q_cell)
global bj_ex2HD
global h ny nx
global gfuns pfuns

% Obtain cost derived by Greenshield's model and discomfort function
% Obtain magnitude of velosity
method = 'nearest';% 'linear'/'nearest'/'pchip'/'spline'
density_ex = gfuns.Boundary_ex(Q_cell{1},ones(ny,nx),'x',2,method);
density_ex = gfuns.Boundary_ex(density_ex,ones(ny,nx+4),'y',2,method);
density_ex = density_ex.*bj_ex2HD;
velosity = gfuns.Velosity(density_ex);
cost_1 = 1./velosity;
cost_2 = 0.02.*density_ex.^2;
cost =  cost_1 + cost_2;
potential = fast_WENO3(cost);

% Obtain potential by fast sweeping WENO using Gauss-Seidel iterations
method = 'linear';% 'linear'/'nearest'/'pchip'/'spline'
potential_x_ex = gfuns.Boundary_ex(potential,bj_calDH,'x',1,method);
potential_y_ex = gfuns.Boundary_ex(potential,bj_calDH,'y',1,method);
potential_x = -(potential_x_ex(:,3:end)-potential_x_ex(:,1:end-2))./2./h;
potential_y = -(potential_y_ex(1:end-2,:)-potential_y_ex(3:end,:))./2./h;
phy_grad = sqrt(potential_x.^2+potential_y.^2);
ve_x = potential_x./phy_grad;
ve_y = potential_y./phy_grad;
ve_x(isnan(ve_x)) = 0;
ve_y(isnan(ve_y)) = 0;

end