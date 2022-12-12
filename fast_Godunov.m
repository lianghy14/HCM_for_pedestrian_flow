function potential_IC = fast_Godunov(cost,boundary,x,y)
% tstart=tic;
global boundary_H
global h 
global sigma gfuns
nx = length(x);
ny = length(y);
boundary_judge = ones(ny,nx);
boundary_judge_H = gfuns.Boundary_value(x,y,boundary_judge,boundary_H,0);
boundary_judge_HF = gfuns.Boundary_value(x,y,boundary_judge_H,boundary,0);

potential_IC = 10^12 .* ones(ny,nx);%coordinate potantial(j,i):x = x(i); y = y(j)
potential_old = zeros(ny,nx);
potential_IC = gfuns.Boundary_value(x,y,potential_IC,boundary_H,10^12);
potential_IC = gfuns.Boundary_value(x,y,potential_IC,boundary,0);
% method = 'linear';
iteration = 0;
GS_it = {[1 nx;1 ny],[nx 1;1 ny],[nx 1;ny 1],[1 nx;ny 1]};
GS_i = [1 1;-1 1;-1 -1;1 -1];
while (norm(potential_IC - potential_old)>=sigma)
    iteration = iteration + 1;
%     potential_xex = gfuns.Boundary_ex(potential_IC,boundary_judge_H,'x',1,method);
%     potential_yex = gfuns.Boundary_ex(potential_IC,boundary_judge_H,'y',1,method);
    potential_old = potential_IC;
    for it = 1:4
        for i = GS_it{it}(2,1):GS_i(it,2):GS_it{it}(2,2)
            for j = GS_it{it}(1,1):GS_i(it,1):GS_it{it}(1,2)
                if boundary_judge_HF(i,j) ~= 0
    %                 potential_x_neighbor = [potential_xex(i,j) potential_xex(i,j+1) potential_xex(i,j+2)];
    %                 potential_y_neighbor = [potential_yex(i,j) potential_yex(i+1,j) potential_yex(i+2,j)];
                    potential_min_x = min([potential_IC(i,min(nx,max(1,j-1))) potential_IC(i,min(nx,max(1,j+1)))]);
                    potential_min_y = min([potential_IC(min(ny,max(1,i+1)),j) potential_IC(min(ny,max(1,i-1)),j)]);
                    if abs(potential_min_x-potential_min_y)-(cost(i,j)*h)>=10^-8
                        potential_IC(i,j) = min([potential_min_x+cost(i,j)*h potential_min_y+cost(i,j)*h]);
                    else
                        potential_IC(i,j)= (potential_min_x+potential_min_y+(2*(cost(i,j)^2)*(h^2)-(potential_min_x-potential_min_y)^2)^0.5)/2;
                    end
                end
            end
        end
    end
    if iteration >= 1000
        debug = 0;
        potential_IC = null;
        return;
    end
end
% tend = toc(tstart);
% fprintf('IC for potential. Computing time is %f. Iteration number is %d\n',tend,iteration);
%% plot
% number = plot_potential(boundary_H,potential_IC,x,y);

end
