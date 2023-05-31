function potential = fast_WENO3(cost,h,bj_ex2ODH)

% GS iterations
ny_ex2 = size(cost,1);
nx_ex2 = size(cost,2);
sigma = 10^-11;
GS_it = {[1 nx_ex2;1 ny_ex2],[nx_ex2 1;1 ny_ex2],[nx_ex2 1;ny_ex2 1],[1 nx_ex2;ny_ex2 1]};
GS_i = [1 1;-1 1;-1 -1;1 -1];
% IC for potential using Godunov
potential = fast_Godunov(cost,h,bj_ex2ODH);
potential_old = zeros(ny_ex2,nx_ex2);
% Starting fast sweeping
iteration = 0;
% 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
bj_2 = bj_ex2ODH; bj_2(bj_2==2) = 0;
sigma_old = zeros(50);
weight = zeros(2,ny_ex2,nx_ex2,2);
nit_zhang = 0;

while norm(potential - potential_old)>=sigma
    iteration = iteration + 1;
    fprintf('It No. is %d. Nit No. is %d. Error is %.9f\n',iteration, nit_zhang, norm(potential - potential_old));
    for it = 1:4
        potential_old = potential;
        for i = GS_it{it}(2,1):GS_i(it,2):GS_it{it}(2,2)
            for j = GS_it{it}(1,1):GS_i(it,1):GS_it{it}(1,2)
                if (bj_2(i,j) > 2.5)&&all(bj_2((i-2):(i+2),j) > 2.5)&&all(bj_2(i,(j-2):(j+2)) > 2.5)
                    % Potential xmin and ymin
                    % ((j-2):1:(j+2)) ((i-2):1:(i+2))
                    % ((j-2*GS_i(it,1)):GS_i(it,1):(j+2*GS_i(it,1)))
                    % ((i-2*GS_i(it,2)):GS_i(it,2):(i+2*GS_i(it,2)))
                    potential_x_neighbor = potential( i,((j-2):1:(j+2)) );
                    potential_y_neighbor = potential(  (i-2):1:(i+2) ,j );
                    w1 = reshape(weight(1,i,j,:),[1,2]);
                    w2 = reshape(weight(1,i,j,:),[1,2]);
                    [potential_min_x,w1] = fast_WENO3_reconstruction(nit_zhang,w1,potential_x_neighbor,h);
                    [potential_min_y,w2] = fast_WENO3_reconstruction(nit_zhang,w2,potential_y_neighbor,h);
                    weight(1,i,j,:) = w1;
                    weight(2,i,j,:) = w2;
                    % Solution
                    if (abs(potential_min_x-potential_min_y)-(cost(i,j)*h))>1d-11
                        potential(i,j) = min(potential_min_x,potential_min_y) + cost(i,j)*h;
                    else
                        potential(i,j)= (potential_min_x+potential_min_y+(2.*(cost(i,j).^2).*(h.^2)-(potential_min_x-potential_min_y).^2).^0.5)./2;
                    end
                end
            end
        end
        if (nit_zhang < 50)
            sigma_old(nit_zhang+1) = norm(potential - potential_old);
            if ( abs( log10(sigma_old(nit_zhang+1)) - log10(sum(sigma_old(1:(nit_zhang+1)))/(nit_zhang+1)) ) < 1.0 )
                nit_zhang = nit_zhang+ 1;
            else
                nit_zhang = 0;
            end
        else
            nit_zhang = nit_zhang + 1;
        end

    end

    if iteration == 500
        break
    end
    if iteration >= 1000
        error('FSM error!');
    end
end
% potential = potential(3:end-2,3:end-2);
% tend = toc(tstart);
fprintf('WENO3 for potential. Iteration No. is %d\n',iteration);
% pfuns = functions_plot;
% number = pfuns.potential(boundary_H,potential,x,y);
end


function [potential_min,weight] = fast_WENO3_reconstruction(it_non,weight,potential_neighbor,h)
epsilon = 10^-6;
if (it_non < 50)
    r_m = (epsilon + (potential_neighbor(1)-2.*potential_neighbor(2)+potential_neighbor(3)).^2)./(epsilon + (potential_neighbor(2)-2.*potential_neighbor(3)+potential_neighbor(4)).^2);
    r_p = (epsilon + (potential_neighbor(5)-2.*potential_neighbor(4)+potential_neighbor(3)).^2)./(epsilon + (potential_neighbor(4)-2.*potential_neighbor(3)+potential_neighbor(2)).^2);
    w_m = 1 ./ (1 + 2 .* r_m.^2);
    w_p = 1 ./ (1 + 2 .* r_p.^2);
    weight(1) = (w_m + weight(1)*it_non) / (it_non+1);
    weight(2) = (w_p + weight(2)*it_non) / (it_non+1);
else
    w_m = weight(1);
    w_p = weight(2);
end

potential_upm = (1-w_m).*(potential_neighbor(4)-potential_neighbor(2))./(2.*h) + w_m.*(3.*potential_neighbor(3)-4.*potential_neighbor(2)+potential_neighbor(1))./(2.*h);
potential_upp = (1-w_p).*(potential_neighbor(4)-potential_neighbor(2))./(2.*h) + w_p.*(-3.*potential_neighbor(3)+4.*potential_neighbor(4)-potential_neighbor(5))./(2.*h);

potential_min = min(potential_neighbor(3) - h.*potential_upm, potential_neighbor(3) + h.*potential_upp);

end