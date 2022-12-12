function [Q,alpha] = TVD_RK(Q_cell,t,tstep)
global area x y nx ny
global tStart tEnd
global plot_config plot_dt dir_fig
global f1 t_bound h x_max y_max boundary_H boundary_O boundary_D boundary_D_2h boundary_D_outside
global m_ave v_free gama_1 CFL sigma epsilon
global gfuns pfuns

% Initialization
F_in = gfuns.F_in(t);
Q = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};

% Cost potential: FSM
tstart = tic;
[ve_x,ve_y] = RHS_FSM(Q_cell);
% tend_FSM = toc(tstart);

% Numerical scheme
Q_cell{2} = Q_cell{1}.*gfuns.Velosity(Q_cell{1}).*ve_x;
Q_cell{3} = Q_cell{1}.*gfuns.Velosity(Q_cell{1}).*ve_y;
[alpha,div_F] = Res_LF(Q_cell,F_in);
Q{1} = Q_cell{1}+tstep.*div_F;

alpha = max(alpha);

% tend_LFHE = toc(tstart);
% fprintf('LF, time is %.2f, Computing time (FSM,LF) is %.2f %.2f\n',t,tend_FSM,tend_LFHE);
end

