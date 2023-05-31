function [Q_3,alpha] = TVD_RK(Q_cell,t,tstep)
global h x y nx ny
global gfuns pfuns

%% TVD-RK3
Q_1 = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
Q_2 = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
Q_3 = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};

% STEP 1
F_in = gfuns.F_in(t);
[ve_x,ve_y] = RHS_FSM(Q_cell);
[alpha_1,div_F] = Res_WENO5(Q_cell,F_in,0);
Q_1{1} = Q_cell{1}+tstep.*div_F;
Q_1{2} = Q_1{1}.*gfuns.Velosity(Q_1{1}).*ve_x;
Q_1{3} = Q_1{1}.*gfuns.Velosity(Q_1{1}).*ve_y;

% STEP 2
F_in = gfuns.F_in(t+tstep);
[ve_x,ve_y] = RHS_FSM(Q_1);
[alpha_2,div_F] = Res_WENO5(Q_1,F_in,0);
Q_2{1} = Q_1{1}+tstep.*div_F;
Q_2{2} = Q_1{1}.*gfuns.Velosity(Q_2{1}).*ve_x;
Q_2{3} = Q_1{1}.*gfuns.Velosity(Q_2{1}).*ve_y;

for i = 1:3
    Q_2{i} = 3./4.*Q_cell{i} + 1./4.*Q_2{i};
end

% SREP 3
F_in = gfuns.F_in(t+1/2*tstep);
[ve_x,ve_y] = RHS_FSM(Q_2);
[alpha_3,div_F] = Res_WENO5(Q_2,F_in,0);
Q_3{1} = Q_2{1}+tstep.*div_F;
Q_3{2} = Q_2{1}.*gfuns.Velosity(Q_3{1}).*ve_x;
Q_3{3} = Q_2{1}.*gfuns.Velosity(Q_3{1}).*ve_y;

for i = 1:3
    Q_3{i} = 1./3.*Q_cell{i} + 2./3.*Q_3{i};
end

alpha = max([alpha_1,alpha_2,alpha_3]);

end

