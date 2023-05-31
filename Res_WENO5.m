function [alpha,div] = Res_WENO5(Q_cell,F_in,G_in)
global h x y nx ny
global bj_calH
global bj_ex2xD bj_ex2yD
global bj_midxH bj_midyH bj_midxO bj_midyO
global gfuns pfuns

% Q_vector
n = 2;
Q_cell_xex = gfuns.Boundary_ex(Q_cell{1},bj_calH,'x',n,'symmetric+');
Q_cell_yex = gfuns.Boundary_ex(Q_cell{1},bj_calH,'y',n,'symmetric+');
F_cell_xex = gfuns.Boundary_ex(Q_cell{2},bj_calH,'x',n,'symmetric-');
G_cell_yex = gfuns.Boundary_ex(Q_cell{3},bj_calH,'y',n,'symmetric-');
Q_cell_xex_D = gfuns.Boundary_ex(Q_cell{1},ones(ny,nx),'x',n,'nearest');
F_cell_xex_D = gfuns.Boundary_ex(Q_cell{2},ones(ny,nx),'x',n,'nearest');
G_cell_yex_D = gfuns.Boundary_ex(Q_cell{3},ones(ny,nx),'y',n,'nearest');

Q_cell_xex = Q_cell_xex.*bj_ex2xD+Q_cell_xex_D.*(1-bj_ex2xD);
F_cell_xex = F_cell_xex.*bj_ex2xD+F_cell_xex_D.*(1-bj_ex2xD);
G_cell_yex = G_cell_yex.*bj_ex2yD+G_cell_yex_D.*(1-bj_ex2yD);


% Alpha
% +c0.*Q_cell{1}./density_m
velosity_xex = sqrt(2).*F_cell_xex.*Q_cell_xex./sqrt(Q_cell_xex.^4+max(Q_cell_xex,10^-6).^4);
velosity_yex = sqrt(2).*G_cell_yex.*Q_cell_yex./sqrt(Q_cell_yex.^4+max(Q_cell_yex,10^-6).^4);

alpha_F = max(abs(velosity_xex),[],2)*ones(1,nx+2*n);
alpha_G = ones(ny+2*n,1)*max(abs(velosity_yex),[],1);

alpha = max(max(max(alpha_F)),max(max(alpha_G)));
% Obtain deviation of density

F_cell_right = 1./2 .* (F_cell_xex+alpha_F.*Q_cell_xex);
F_cell_left = 1./2 .* (F_cell_xex-alpha_F.*Q_cell_xex);
G_cell_right = 1./2 .* (G_cell_yex+alpha_G.*Q_cell_yex);
G_cell_left = 1./2 .* (G_cell_yex-alpha_G.*Q_cell_yex);
F_cell_right_mid = Reconstruction_WENO5(F_cell_right,[0 -1]);
F_cell_left_mid = Reconstruction_WENO5(F_cell_left,[0 1]);
G_cell_right_mid = Reconstruction_WENO5(G_cell_right,[1 0]);
G_cell_left_mid = Reconstruction_WENO5(G_cell_left,[-1 0]);
F_cell_mid = (F_cell_right_mid(:,n:(end-n)) + F_cell_left_mid(:,(n+1):(end-n+1))).*bj_midxH;
G_cell_mid = (G_cell_right_mid((n+1):(end-n+1),:) + G_cell_left_mid(n:(end-n),:)).*bj_midyH;

F_cell_mid  = F_cell_mid.*bj_midxO + (1-bj_midxO).*F_in;
G_cell_mid  = G_cell_mid.*bj_midyO + (1-bj_midyO).*G_in;

div = - 1 ./ h .* (F_cell_mid(:,2:end) - F_cell_mid(:,1:end-1));
div = div - 1 ./ h .* (G_cell_mid(1:end-1,:) - G_cell_mid(2:end,:));

end

function flow_p = Reconstruction_WENO5(flow,ind)
global epsilon 
% Compute nvmerical fluxes at cell 'i' interfaces.
% Set Variables
flow_p1 =circshift(flow,ind);
flow_p2 =circshift(flow_p1,ind);
flow_m1 =circshift(flow,-ind);
flow_m2 =circshift(flow_m1,-ind);
% Reconstruction Polynomials
up1 = 1./3.*flow_m2 - 7./6.*flow_m1 + 11./6.*flow;
up2 = -1./6.*flow_m1 + 5./6.*flow + 1./3.*flow_p1;
up3 = 1./3.*flow + 5./6.*flow_p1 - 1./6.*flow_p2;
% Smooth parameters
b1 = 13./12 .* (flow_m2-2.*flow_m1+flow).^2 + 1 ./4 .* (flow_m2-4.*flow_m1+3.*flow).^2; 
b2 = 13./12 .* (flow_m1-2.*flow+flow_p1).^2 + 1 ./4 .* (flow_m1-flow_p1).^2; 
b3 = 13./12 .* (flow-2.*flow_p1+flow_p2).^2 + 1 ./4 .* (3.*flow-4.*flow_p1+flow_p2).^2; 
% Constants
g1 = 1/10; g2 = 3/5; g3 = 3/10;
% weigths
wt1 = g1 ./ (epsilon+b1).^2;
wt2 = g2 ./ (epsilon+b2).^2;
wt3 = g3 ./ (epsilon+b3).^2;
sum_wt = wt1 + wt2 +wt3;
% Non-linear weigths
w1 = wt1 ./ sum_wt;
w2 = wt2 ./ sum_wt;
w3 = wt3 ./ sum_wt;
% WENO polynomial
flow_p = w1.*up1 + w2.*up2 + w3.*up3;
% End of reconstruction.
end