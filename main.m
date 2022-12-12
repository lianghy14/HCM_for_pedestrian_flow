clc
clear all
global area x y nx ny
global tStart tEnd
global plot_config plot_dt dir_fig
global f1 t_bound h x_max y_max boundary_H boundary_O boundary_D boundary_D_2h boundary_D_outside
global m_ave v_free gama_1 CFL sigma epsilon
global gfuns pfuns
%% Paramters

gfuns = functions_given; pfuns = functions_plot;

gfuns.Para();
area = [0,x_max;-2,y_max+2]; 
tstep = 0.01;

% Computational domain
x = (area(1,1)+0.5*h):h:(area(1,2)-0.5*h); nx = length(x);
y = (area(2,2)-0.5*h):-h:(area(2,1)+0.5*h);  ny = length(y);
[X,Y] = meshgrid(x,y);

%% IC for density potential
cpu_tstart = tic;
mkdir(dir_fig);
fprintf('Start Computing---------------------------------------------\n');
% Initial Condition
t = tStart;
if tStart == 0
    Q_cell = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
    Q_cell{1} = gfuns.Boundary_value(x,y,Q_cell{1},boundary_H,0);
else
    Q_name = [dir_fig num2str(tStart) '.mat']; 
    load([dir_fig num2str(tStart) '.mat']);
end
alpha = 5;
%% LF scheme
while (t<=tEnd)
%% Plot&Save
    if(mod(t,1)==0)
        if (mod(t,plot_dt)==0) && (plot_config~=0)
            density_f = zeros(ny,nx,8);
            fig = figure(1);
            density = gfuns.Boundary_value(x,y,Q_cell{1},boundary_H,nan);
            subplot(2,1,1); mesh(X,Y,density); colormap Jet; axis([0,x_max,0,y_max,0,6]);
            title(['HCM' ', dx = ',num2str(h),', dy = ',num2str(h),', time: ',num2str(t)]);
            xlabel('x(m)'); ylabel('y(m)'); zlabel('density(ped/m^2)');
            subplot(2,1,2); imagesc(x,y,density,'alphadata',~isnan(density)); colormap Jet; set(gca,'YDir','normal','color',0*[1 1 1]);caxis([0 6]); colorbar; axis([0,x_max,0,y_max]); 
            xlabel('x(m)'); ylabel('y(m)');
            
            set(fig,'unit','centimeters','position',[10 5 15 12]);
            saveas(fig,[dir_fig num2str(t) '.png']);
        end
        tend = toc(cpu_tstart);
        fprintf(['HCM.t = %d. '],t)
        fprintf('Computing time is %f, Alpha is %f\n',tend,alpha);
        Q_name = [dir_fig num2str(t) '.mat'];
        save(Q_name,'Q_cell*');
        tstep = min([CFL*h/alpha, CFL]);
    else
        tstep = min([CFL*h/alpha, ceil(t)-t, CFL]);
    end

    [Q_cell,alpha] = TVD_RK(Q_cell,t,tstep);
    t = t+tstep;
end
