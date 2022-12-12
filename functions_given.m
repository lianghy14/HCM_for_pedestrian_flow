% This function defines all inputs into this HCM simulation

function gfuns = functions_given
gfuns.Para = @Para_def;
gfuns.F_in = @F_in_cal;
gfuns.Velosity = @Velosity_cal;
gfuns.Boundary_ex = @Boundary_ex_cal;
gfuns.Boundary_ex2 = @Boundary_ex2_cal;
gfuns.Boundary_value = @Boundary_value_cal;
end

%% Parameters & Layout
function Para_def()
global tStart tEnd method 
global plot_config plot_dt dir_fig
global f1 t_bound h x_max y_max boundary_H boundary_O boundary_D boundary_D_2h boundary_D_outside
global m_ave v_free gama_1 CFL sigma epsilon

% Simulation time
tStart = 0;
tEnd = 300;
method = 'LF';

% Layout settings
f1 = 1; t_bound = 120;
h = 1; x_max = 100; y_max = 50;
boundary_O = {'Rectangle',1,[-5 0;0 y_max]};
boundary_D = {'Rectangle',1,[x_max 0;x_max+5 y_max]};
boundary_D_2h = {'Rectangle',1,[x_max-h 0;x_max+5 y_max]};
boundary_D_outside = {'Rectangle',1,[x_max+h 0;x_max+5 y_max]};
boundary_H = {'Rectangle',3,[-5 -5;x_max+5 0],[-5 y_max;x_max+5 y_max+5],[50 10;70 30]};

% Parameters
m_ave = 60;
v_free = 1.034;
gama_1 = -0.09;
CFL = 0.2;
sigma = 10^-9;
epsilon = 10^-6;

% plot configration.
plot_config = 1;
plot_dt = 50;
dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study HCM\' 'method' '_h' num2str(h) '_fin' num2str(f1) '/'];

end

%% Fundmental diagram
function Velosity = Velosity_cal(density)
global v_free gama_1
alpha = 2;
Velosity = v_free.*exp(gama_1.*max(0,density).^alpha);
end

%% Flow in
function F_in = F_in_cal(t)
global t_bound f1
switch true
    case t<=60
        F_in = f1*(t / 60);
    case t>60 && t<=t_bound
        F_in = f1;
    case t>t_bound
        F_in = max(0,f1*((t_bound-t+60) / 60));
end
end

%% Boundary expand
function cell_ex = Boundary_ex_cal(cell,cell_judge,dim,n,method)
nx = length(cell(1,:));
ny = length(cell(:,1));
switch dim
    case 'x'
        cell_ex = zeros(ny,nx+2*n);
        cell_judge = [cell_judge,zeros(ny,1)];
        for i = 1:ny
            n_loc = n+1;
            while n_loc<=(nx+n)
                if cell_judge(i,n_loc-n) == 0
                    cell_ex(i,n_loc) = cell(i,n_loc-n);
                    n_loc = n_loc+1;
                else
                    start_loc = n_loc-n;
                    end_loc = find(cell_judge(i,(n_loc-n):(nx+1))==0,1,'first')-1+n_loc-n-1;

                    switch method
                        case 'zero'
                            cell_ex(i,(start_loc):(start_loc+n-1)) = 0;
                            cell_ex(i,(start_loc+n):(end_loc+n)) =  cell(i,(start_loc):(end_loc));
                            cell_ex(i,(end_loc+n+1):(end_loc+2*n)) = 0;
                        case 'symmetric+'
                            cell_ex(i,(start_loc):(start_loc+n-1)) = cell(i,(start_loc+n-1):-1:(start_loc));
                            cell_ex(i,(start_loc+n):(end_loc+n)) =  cell(i,(start_loc):(end_loc));
                            cell_ex(i,(end_loc+n+1):(end_loc+2*n)) = cell(i,(end_loc):-1:(end_loc-n+1));
                        case 'symmetric-'
                            cell_ex(i,(start_loc):(start_loc+n-1)) = - cell(i,(start_loc+n-1):-1:(start_loc));
                            cell_ex(i,(start_loc+n):(end_loc+n)) =  cell(i,(start_loc):(end_loc));
                            cell_ex(i,(end_loc+n+1):(end_loc+2*n)) = - cell(i,(end_loc):-1:(end_loc-n+1));
                        otherwise
                            x = (n+1):(n+end_loc-start_loc+1);
                            x_expand = 1:(n+end_loc-start_loc+1+n);
                            cell_ex(i,(start_loc):(end_loc+2*n)) = interp1(x,cell(i,start_loc:end_loc),x_expand,method,'extrap');
                    end
                    n_loc = end_loc+2*n+1;
                end
            end
        end
    case 'y'
        cell_ex = zeros(ny+2*n,nx);
        cell_judge = [cell_judge;zeros(1,nx)];
        for i = 1:nx
            n_loc = n+1;
            while n_loc<=(ny+n)
                if cell_judge(n_loc-n,i) == 0
                    cell_ex(n_loc,i) = cell(n_loc-n,i);
                    n_loc = n_loc+1;
                else
                    start_loc = n_loc-n;
                    end_loc = find(cell_judge((n_loc-n):(ny+1),i)==0,1,'first')-1+n_loc-n-1;
                    switch method
                        case 'zero'
                            cell_ex((start_loc):(start_loc+n-1),i) = 0;
                            cell_ex((start_loc+n):(end_loc+n),i) =  cell((start_loc):(end_loc),i);
                            cell_ex((end_loc+n+1):(end_loc+2*n),i) = 0;
                        case 'symmetric+'
                            cell_ex((start_loc):(start_loc+n-1),i) = cell((start_loc+n-1):-1:(start_loc),i);
                            cell_ex((start_loc+n):(end_loc+n),i) =  cell((start_loc):(end_loc),i);
                            cell_ex((end_loc+n+1):(end_loc+2*n),i) = cell((end_loc):-1:(end_loc-n+1),i);
                        case 'symmetric-'
                            cell_ex((start_loc):(start_loc+n-1),i) = - cell((start_loc+n-1):-1:(start_loc),i);
                            cell_ex((start_loc+n):(end_loc+n),i) =  cell((start_loc):(end_loc),i);
                            cell_ex((end_loc+n+1):(end_loc+2*n),i) = - cell((end_loc):-1:(end_loc-n+1),i);
                        otherwise
                            y = (n+1):(n+end_loc-start_loc+1);
                            y_expand = 1:(n+end_loc-start_loc+1+n);
                            cell_ex((start_loc):(end_loc+2*n),i) = interp1(y,cell(start_loc:end_loc,i),y_expand,method,'extrap');
                    end
                    n_loc = end_loc+2*n+1;
                end
            end
        end
end

%% Simple zero
% switch dim
%     case 1
%         flow_expand = zeros(ny,nx+2*n);
%         flow_expand(:,(n+1):(n+nx)) = flow;
%         flow_expand(:,[1:n (n+nx+1):(n+nx+n)]) = 0;
%     case 2
%         flow_expand = zeros(ny+2*n,nx);
%         flow_expand((n+1):(n+ny),:) = flow;
%         flow_expand([1:n (n+ny+1):(n+ny+n)],:) = 0;
% end;

end

function cell_ex2 = Boundary_ex2_cal(cell,cell_judge,n,method)
nx = length(cell(1,:));
ny = length(cell(:,1));
cell_judge(cell_judge~=0) = 1;

cell_vector = zeros(sum(sum(cell_judge)),3);
count = 1;
for i = 1:ny
    for j = 1:nx
        if cell_judge(i,j) ~= 0
            cell_vector(count,:) = [j,i,cell(i,j)];
            count = count+1;
        end
    end
end
[X,Y] = meshgrid((1-n):(nx+n),(1-n):(ny+n));
cell_ex2 = griddata(cell_vector(:,1),cell_vector(:,2),cell_vector(:,3),X,Y,method);
end

%% Boundary value
function cell = Boundary_value_cal(x,y,cell,boundary,value)
n = 1;
while(n<=length(boundary))
    switch boundary{n}
        case 'Rectangle'
            for k = (n+2):(n+1+boundary{n+1})
                for i = length(cell(1,:)):-1:1
                    for j = 1:length(cell(:,1))
                        if (x(i)>=boundary{k}(1,1))&&(x(i)<=boundary{k}(2,1))&&(y(j)>=boundary{k}(1,2))&&(y(j)<=boundary{k}(2,2))
                                cell(j,i) = value;
                        end
                    end
                end
            end
            n = n + boundary{n+1}+2;
        case 'Circle'
            for k = (n+2):(n+1+boundary{n+1})
                for i = length(cell(1,:)):-1:1
                    for j = 1:length(cell(:,1))
                        if ((x(i)-boundary{k}(1))^2+(y(j)-boundary{k}(2))^2) <= boundary{k}(3)^2
                                cell(j,i) = value;
                        end
                    end
                end
            end
            n = n + boundary{n+1}+2;
    end
end
end
