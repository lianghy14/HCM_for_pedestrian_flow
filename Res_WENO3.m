function [alpha,div] = Res_WENO3(Q_cell)
global x y h nx ny gfuns boundary_O boundary_D boundary_H c0 beta density_m d_min;
global D_in
% tstart=tic;
div = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
% Boundary expand
% Q_vector
boundary_judge_x = ones(ny,nx);
boundary_judge_y = ones(ny,nx);
Q_cell_xex = cell(1,3);Q_cell_xex_D = cell(1,3);
Q_cell_yex = cell(1,3);Q_cell_yex_D = cell(1,3);
n = 2;
boundary_judge_H = gfuns.Boundary_value(x,y,boundary_judge_x,boundary_H,0);
Q_cell_xex{1} = gfuns.Boundary_ex(Q_cell{1},boundary_judge_H,'x',n,'symmetric+');
Q_cell_yex{1} = gfuns.Boundary_ex(Q_cell{1},boundary_judge_H,'y',n,'symmetric+');
Q_cell_xex{2} = gfuns.Boundary_ex(Q_cell{2},boundary_judge_H,'x',n,'symmetric-');
Q_cell_yex{2} = gfuns.Boundary_ex(Q_cell{2},boundary_judge_H,'y',n,'symmetric+');
Q_cell_xex{3} = gfuns.Boundary_ex(Q_cell{3},boundary_judge_H,'x',n,'symmetric+');
Q_cell_yex{3} = gfuns.Boundary_ex(Q_cell{3},boundary_judge_H,'y',n,'symmetric-');
for i = 1:3
    Q_cell_xex_D{i} = gfuns.Boundary_ex(Q_cell{i},boundary_judge_x,'x',n,'nearest');
    Q_cell_yex_D{i} = gfuns.Boundary_ex(Q_cell{i},boundary_judge_y,'y',n,'nearest');
end
x_ex = (x(1)-n*h):h:(x(end)+n*h);nx_ex = length(x_ex);
y_ex = (y(1)+n*h):-h:(y(end)-n*h);ny_ex = length(y_ex);
boundary_judge_x = ones(ny,nx_ex);
boundary_judge_y = ones(ny_ex,nx);
boundary_judge_D_x = gfuns.Boundary_value(x_ex,y,boundary_judge_x,boundary_D,0);
boundary_judge_D_y = gfuns.Boundary_value(x,y_ex,boundary_judge_y,boundary_D,0);
% boundary_judge_O = gfuns.Boundary_value(x,y,boundary_judge,boundary_O,0);
for i = 1:3
    Q_cell_xex{i} = Q_cell_xex{i}.*boundary_judge_D_x+Q_cell_xex_D{i}.*(1-boundary_judge_D_x);
    Q_cell_yex{i} = Q_cell_yex{i}.*boundary_judge_D_y+Q_cell_yex_D{i}.*(1-boundary_judge_D_y);
%     Q_cell_xex{i} = Q_cell_xex{i}.*boundary_judge_O+Q_cell_xex_D{i}.*(1-boundary_judge_O);
end
Q_cell_xex{1}  = gfuns.Boundary_value(x_ex,y,Q_cell_xex{1},boundary_O,D_in(1));
Q_cell_xex{2}  = gfuns.Boundary_value(x_ex,y,Q_cell_xex{2},boundary_O,D_in(2));
Q_cell_xex{3}  = gfuns.Boundary_value(x_ex,y,Q_cell_xex{3},boundary_O,D_in(3));
% F_vector
F_cell_ex = cell(1,3);G_cell_ex = cell(1,3);
F_cell_ex{1} = Q_cell_xex{2};
F_cell_ex{2} = Q_cell_xex{2}.^2 ./ Q_cell_xex{1};
F_cell_ex{3} = Q_cell_xex{2}.*Q_cell_xex{3}./ Q_cell_xex{1};
F_cell_ex{2}(isnan(F_cell_ex{2})) = 0;F_cell_ex{3}(isnan(F_cell_ex{3})) = 0;
F_cell_ex{2} = F_cell_ex{2} + gfuns.Force_i(Q_cell_xex{1});
% G_vector
G_cell_ex{1} = Q_cell_yex{3};
G_cell_ex{2} = Q_cell_yex{2}.*Q_cell_yex{3}./ Q_cell_yex{1};
G_cell_ex{3} = Q_cell_yex{3}.^2 ./ Q_cell_yex{1};
G_cell_ex{2}(isnan(G_cell_ex{2})) = 0;G_cell_ex{3}(isnan(G_cell_ex{3})) = 0;
G_cell_ex{3} = G_cell_ex{3} + gfuns.Force_i(Q_cell_yex{1});
% Alpha
% +c0.*(Q_cell_xex{1}./density_m).^beta
% c0.*(cos(pi.*density./density_m)+1)
alpha_F = max(0,max(max(abs(F_cell_ex{1}./Q_cell_xex{1})+c0.*(Q_cell_xex{1}./density_m).^beta )));
alpha_G = max(0,max(max(abs(G_cell_ex{1}./Q_cell_yex{1})+c0.*(Q_cell_yex{1}./density_m).^beta )));
alpha = max(alpha_F,alpha_G);
% Obtain deviation of density
F_cell_right = cell(1,3);F_cell_left = cell(1,3);
G_cell_right = cell(1,3);G_cell_left = cell(1,3);
F_cell_right_mid = cell(1,3);F_cell_left_mid = cell(1,3);
G_cell_right_mid = cell(1,3);G_cell_left_mid = cell(1,3);
F_cell_mid = cell(1,3);G_cell_mid = cell(1,3);
% F_cell_xex_O = cell(1,3);
% for i = 1:3
%     F_cell_xex_O{i} = gfuns.Boundary_ex(F_cell_ex{i},boundary_judge_O,'x',2,'linear');
% end
% for i = 1:3
%     F_cell_ex{i} = F_cell_ex{i}.*boundary_judge_O+F_cell_xex_O{i}(:,3:end-2).*(1-boundary_judge_O);
% end
for i = 1:3
    F_cell_right{i} = 1./2 .* (F_cell_ex{i}+alpha_F.*Q_cell_xex{i});
    F_cell_left{i} = 1./2 .* (F_cell_ex{i}-alpha_F.*Q_cell_xex{i});
    G_cell_right{i} = 1./2 .* (G_cell_ex{i}+alpha_G.*Q_cell_yex{i});
    G_cell_left{i} = 1./2 .* (G_cell_ex{i}-alpha_G.*Q_cell_yex{i});
    F_cell_right_mid{i} = Reconstruction_WENO3(F_cell_right{i},[0 -1]);
    F_cell_left_mid{i} = Reconstruction_WENO3(F_cell_left{i},[0 1]);
    G_cell_right_mid{i} = Reconstruction_WENO3(G_cell_right{i},[1 0]);
    G_cell_left_mid{i} = Reconstruction_WENO3(G_cell_left{i},[-1 0]);
    F_cell_mid{i} = F_cell_right_mid{i}(:,n:(end-n)) + F_cell_left_mid{i}(:,(n+1):(end-n+1));
    G_cell_mid{i} = G_cell_right_mid{i}((n+1):(end-n+1),:) + G_cell_left_mid{i}(n:(end-n),:);
end
F_cell_mid{1}(:,2:end) = gfuns.Boundary_value(x,y,F_cell_mid{1}(:,2:end),boundary_H,0);
F_cell_mid{1}(:,1:end-1) = gfuns.Boundary_value(x,y,F_cell_mid{1}(:,1:end-1),boundary_H,0);
G_cell_mid{1}(1:end-1,:) = gfuns.Boundary_value(x,y,G_cell_mid{1}(1:end-1,:),boundary_H,0);
G_cell_mid{1}(2:end,:) = gfuns.Boundary_value(x,y,G_cell_mid{1}(2:end,:),boundary_H,0);
x_ex = (x(1)-0.5*h):h:(x(end)+0.5*h);
F_cell_mid{1}  = gfuns.Boundary_value(x_ex,y,F_cell_mid{1},boundary_O,D_in(2));
if D_in(1) == 0
    F_cell_mid{2}  = gfuns.Boundary_value(x_ex,y,F_cell_mid{2},boundary_O,0);
    F_cell_mid{3}  = gfuns.Boundary_value(x_ex,y,F_cell_mid{3},boundary_O,0);
else
    F_cell_mid{2} = gfuns.Boundary_value(x_ex,y,F_cell_mid{2},boundary_O,D_in(2).^2 ./ D_in(1)+gfuns.Force_i(D_in(1)));
    F_cell_mid{3} = gfuns.Boundary_value(x_ex,y,F_cell_mid{3},boundary_O,D_in(2).*D_in(3)./D_in(1));
end
for i = 1:3
    div{i} = - 1 ./ h .* (F_cell_mid{i}(:,2:end) - F_cell_mid{i}(:,1:end-1));
    div{i} = div{i} - 1 ./ h .* (G_cell_mid{i}(1:end-1,:) - G_cell_mid{i}(2:end,:));
end
for i = 1:3
    div{i} = gfuns.Boundary_value(x,y,div{i},boundary_H,0);
    div{i} = gfuns.Boundary_value(x,y,div{i},boundary_D,0);
end
% Boundary conditions
% for i = 2:(nx-1)
%     for j = 1:ny
%         if (boundary_judge_H(j,i-1)==0 && div{2}(j,i)<0)
%             div{2}(j,i) = - 1 ./ h .* (F_cell_mid{2}(j,i+1) - gfuns.Force_i(Q_cell{1}(j,i)));
%             div{2}(j,i) = div{2}(j,i) - 1 ./ h .* (G_cell_mid{2}(j,i) - G_cell_mid{2}(j+1,i));
%         else
%             if (boundary_judge_H(j,i+1)==0 && div{2}(j,i)>0)
%                 div{2}(j,i) = - 1 ./ h .* (gfuns.Force_i(Q_cell{1}(j,i)) - F_cell_mid{2}(j,i));
%                 div{2}(j,i) = div{2}(j,i) - 1 ./ h .* (G_cell_mid{2}(j,i) - G_cell_mid{2}(j+1,i));
%             end
%         end
%     end
% end
% for i = 1:nx
%     for j = 2:(ny-1)
%         if (boundary_judge_H(j-1,i)==0 && div{3}(j,i)>0)
%             div{3}(j,i) = - 1 ./ h .* (F_cell_mid{3}(j,i+1) - F_cell_mid{3}(j,i));
%             div{3}(j,i) = div{3}(j,i) - 1 ./ h .* (gfuns.Force_i(Q_cell{1}(j,i)) - G_cell_mid{3}(j+1,i));
%         else
%             if (boundary_judge_H(j+1,i)==0 && div{3}(j,i)<0)
%                 div{3}(j,i) = - 1 ./ h .* (F_cell_mid{3}(j,i+1) - F_cell_mid{3}(j,i));
%                 div{3}(j,i) = div{3}(j,i) - 1 ./ h .* (G_cell_mid{3}(j,i) - gfuns.Force_i(Q_cell{1}(j,i)));
%             end
%         end
%     end
% end
% tend = toc(tstart);
% fprintf('Reconstruction. Computing time is %f\n',tend);
end

function flow_p = Reconstruction_WENO3(flow,ind)
global epsilon 
% Compute nvmerical fluxes at cell 'i' interfaces.
% Set Variables
flow_p1 =circshift(flow,ind);
flow_m1 =circshift(flow,-ind);
% Reconstruction Polynomials
up1 = - 1./2.*flow_m1 + 3./2.*flow;
up2 = 1./2.*flow + 1./2.*flow_p1;
% Smooth parameters
b1 = (flow_m1-flow).^2; 
b2 = (flow-flow_p1).^2; 
% Constants
g1 = 1/3; g2 = 2/3; 
% weigths
wt1 = g1 ./ (epsilon+b1).^2;
wt2 = g2 ./ (epsilon+b2).^2;
sum_wt = wt1 + wt2;
% Non-linear weigths
w1 = wt1 ./ sum_wt;
w2 = wt2 ./ sum_wt;
% WENO polynomial
flow_p = w1.*up1 + w2.*up2;
% End of reconstruction.
end
% function flow_p = Reconstruction_WENO5(flow,ind)
% global epsilon 
% % Compute nvmerical fluxes at cell 'i' interfaces.
% % Set Variables
% flow_p1 =circshift(flow,ind);
% flow_p2 =circshift(flow_p1,ind);
% flow_m1 =circshift(flow,-ind);
% flow_m2 =circshift(flow_m1,-ind);
% % Reconstruction Polynomials
% up1 = 1./3.*flow_m2 - 7./6.*flow_m1 + 11./6.*flow;
% up2 = -1./6.*flow_m1 + 5./6.*flow + 1./3.*flow_p1;
% up3 = 1./3.*flow + 5./6.*flow_p1 - 1./6.*flow_p2;
% % Smooth parameters
% b1 = 13./12 .* (flow_m2-2.*flow_m1+flow).^2 + 1 ./4 .* (flow_m2-4.*flow_m1+3.*flow).^2; 
% b2 = 13./12 .* (flow_m1-2.*flow+flow_p1).^2 + 1 ./4 .* (flow_m1-flow_p1).^2; 
% b3 = 13./12 .* (flow-2.*flow_p1+flow_p2).^2 + 1 ./4 .* (3.*flow-4.*flow_p1+flow_p2).^2; 
% % Constants
% g1 = 1/10; g2 = 3/5; g3 = 3/10;
% % weigths
% wt1 = g1 ./ (epsilon+b1).^2;
% wt2 = g2 ./ (epsilon+b2).^2;
% wt3 = g3 ./ (epsilon+b3).^2;
% sum_wt = wt1 + wt2 +wt3;
% % Non-linear weigths
% w1 = wt1 ./ sum_wt;
% w2 = wt2 ./ sum_wt;
% w3 = wt3 ./ sum_wt;
% % WENO polynomial
% flow_p = w1.*up1 + w2.*up2 + w3.*up3;
% % End of reconstruction.
% end
