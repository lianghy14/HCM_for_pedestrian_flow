function fast_godunov = fast_Godunov(cost,h,bj_ex2ODH)
% tstart=tic;

ny_ex2 = size(cost,1);
nx_ex2 = size(cost,2);
fast_IC = 10^6 .* ones(ny_ex2,nx_ex2);
fast_IC(bj_ex2ODH==1) = 0;
sigma = 10^-12;
potential_old = zeros(ny_ex2,nx_ex2);
iteration = 0;
GS_it = {[1 nx_ex2;1 ny_ex2],[nx_ex2 1;1 ny_ex2],[nx_ex2 1;ny_ex2 1],[1 nx_ex2;ny_ex2 1]};
GS_i = [1 1;-1 1;-1 -1;1 -1];
while (norm(fast_IC - potential_old)>=sigma)
    iteration = iteration + 1;
    potential_old = fast_IC;
    for it = 1:4
        for i = GS_it{it}(2,1):GS_i(it,2):GS_it{it}(2,2)
            for j = GS_it{it}(1,1):GS_i(it,1):GS_it{it}(1,2)
                if bj_ex2ODH(i,j) >2.5
                    potential_min_x = min([fast_IC(i,j-1) fast_IC(i,j+1)]);
                    potential_min_y = min([fast_IC(i+1,j) fast_IC(i-1,j)]);
                    if abs(potential_min_x-potential_min_y)-(cost(i,j)*h)>=0
                        fast_IC_temp = min([potential_min_x potential_min_y])+cost(i,j)*h;
                    else
                        fast_IC_temp= (potential_min_x+potential_min_y+(2*(cost(i,j)^2)*(h^2)-(potential_min_x-potential_min_y)^2)^0.5)/2;
                    end
                    fast_IC(i,j) = min([fast_IC(i,j),fast_IC_temp]);
                end
            end
        end
%         fprintf('Difference is %f. Iteration number is %d\n',norm(fast_IC - potential_old),iteration);
    end
    if iteration >= 1000
        fast_godunov = null;
        return;
    end
end
fast_godunov = fast_IC;
% tend = toc(tstart);
% fprintf('IC for potential. Computing time is %f. Iteration number is %d\n',tend,iteration);
%% plot
% number = plot_potential(boundary_H,potential_IC,x,y);

end
