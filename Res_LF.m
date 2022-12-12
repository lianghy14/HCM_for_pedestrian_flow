function [alpha,div] = Res_LF(Q_cell,F_in)
global area x y nx ny
global tStart tEnd
global plot_config plot_dt dir_fig
global f1 t_bound h x_max y_max boundary_H boundary_O boundary_D boundary_D_2h boundary_D_outside
global m_ave v_free gama_1 CFL sigma epsilon
global gfuns pfuns

% Q_vector
n = 1;
boundary_judge = ones(ny,nx);
boundary_judge_H = gfuns.Boundary_value(x,y,boundary_judge,boundary_H,0);
Q_cell_xex = gfuns.Boundary_ex(Q_cell{1},boundary_judge_H,'x',n,'symmetric+');
Q_cell_yex = gfuns.Boundary_ex(Q_cell{1},boundary_judge_H,'y',n,'symmetric+');
F_cell_xex = gfuns.Boundary_ex(Q_cell{2},boundary_judge_H,'x',n,'symmetric-');
G_cell_yex = gfuns.Boundary_ex(Q_cell{3},boundary_judge_H,'y',n,'symmetric-');
Q_cell_xex_D = gfuns.Boundary_ex(Q_cell{1},boundary_judge,'x',n,'nearest');
F_cell_xex_D = gfuns.Boundary_ex(Q_cell{2},boundary_judge,'x',n,'nearest');

x_ex = (x(1)-n*h):h:(x(end)+n*h);
boundary_judge = ones(ny,nx+2);
boundary_judge_D = gfuns.Boundary_value(x_ex,y,boundary_judge,boundary_D,0);
Q_cell_xex = Q_cell_xex.*boundary_judge_D+Q_cell_xex_D.*(1-boundary_judge_D);
F_cell_xex = F_cell_xex.*boundary_judge_D+F_cell_xex_D.*(1-boundary_judge_D);


% Alpha
% +c0.*Q_cell{1}./density_m
velosity_xex = sqrt(2).*F_cell_xex.*Q_cell_xex./sqrt(Q_cell_xex.^4+max(Q_cell_xex,10^-6).^4);
velosity_yex = sqrt(2).*G_cell_yex.*Q_cell_yex./sqrt(Q_cell_yex.^4+max(Q_cell_yex,10^-6).^4);

alpha_F = max(abs(velosity_xex),[],2)*ones(1,nx+2*n);
alpha_G = ones(ny+2*n,1)*max(abs(velosity_yex),[],1);

alpha = max(max(max(alpha_F)),max(max(alpha_G)));
% Obtain deviation of density
% boundary_judge_OD = gfuns.Boundary_value(x,y,boundary_judge_O,boundary_D,0);

Q_cell_xp1 = circshift(Q_cell_xex,[0 -1]);
F_cell_xp1 = circshift(F_cell_xex,[0 -1]);
Q_cell_yp1 = circshift(Q_cell_yex,[1 0]);
G_cell_yp1 = circshift(G_cell_yex,[1 0]);

F_cell_mid = 1./2 .* (F_cell_xex+F_cell_xp1-alpha_F.*(Q_cell_xp1-Q_cell_xex));
G_cell_mid = 1./2 .* (G_cell_yex+G_cell_yp1-alpha_G.*(Q_cell_yp1-Q_cell_yex)); 
F_cell_mid = F_cell_mid(:,1:end-1);
G_cell_mid = G_cell_mid(2:end,:);

F_cell_mid(:,2:end) = gfuns.Boundary_value(x,y,F_cell_mid(:,2:end),boundary_H,0);
F_cell_mid(:,1:end-1) = gfuns.Boundary_value(x,y,F_cell_mid(:,1:end-1),boundary_H,0);
G_cell_mid(1:end-1,:) = gfuns.Boundary_value(x,y,G_cell_mid(1:end-1,:),boundary_H,0);
G_cell_mid(2:end,:) = gfuns.Boundary_value(x,y,G_cell_mid(2:end,:),boundary_H,0);
F_cell_mid  = gfuns.Boundary_value(x_ex,y,F_cell_mid,boundary_O,F_in);

div = - 1 ./ h .* (F_cell_mid(:,2:end) - F_cell_mid(:,1:end-1));
div = div - 1 ./ h .* (G_cell_mid(1:end-1,:) - G_cell_mid(2:end,:));

end