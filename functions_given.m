% This function defines all inputs into this HCM simulation

function gfuns = functions_given
gfuns.Para            =   @Para_cal;
gfuns.Layout            =   @Layout_cal;
gfuns.F_in              =   @F_in_cal;
gfuns.Boundary_value    =   @Boundary_value_cal;
gfuns.RRMSE = @RRMSE_cal;
end
%% Parameters
function Para_cal(day)
global Record_dt
global area_cal h boundary_O boundary_D boundary_D_2h boundary_H
global MC_std MC_mean dir_data MC_nseq MC_max

switch day
    case 20230116
        Record_dt = 30;
        MC_std = 0.20;
        MC_mean = 1;
        MC_nseq = 400:400:6000;
        MC_max = 10000;
%         MC_nseq = MC_max;
        h = 1; 
        x_max = 100;  y_max = 50; area_cal = [0,x_max;0,y_max]; 
        sext = 0.1*h;  ext  =  5*h+sext;
        boundary_O      =   {'Rectangle',1, [-ext -sext;sext y_max+sext]};
        boundary_D      =   {'Rectangle',1, [x_max-sext 10+sext;x_max+ext y_max-10-sext]};
        boundary_D_2h   =   {'Rectangle',1, [x_max-sext-h 10+sext;x_max+ext y_max-10-sext]};
        boundary_H      =   {'Rectangle',5, [-ext -ext;x_max+ext -sext],[-ext y_max+sext;x_max+ext y_max+ext],[50-sext 10-sext;70+sext 30+sext],...
                            [x_max-sext y_max-10-sext;x_max+ext y_max+sext],[x_max-sext -sext;x_max+ext 10+sext]};
        dir_data = ['C:\Users\HOWIE-PC\Desktop\Case study\MC simulations Eiko\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
%         dir_data = ['C:\Users\Howie\Desktop\Case study\MC simulations Eiko\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
%         dir_data = ['C:\Users\HLiang\Desktop\Case study\MC simulations Eiko\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
    case 20230102
        Record_dt = 30;
        MC_std = 0.10;
        MC_mean = 1;
        MC_nseq = [20:20:80,100:100:1000];
        MC_max = 1000;
%         MC_nseq = MC_max;
        h = 1; 
        x_max = 100;  y_max = 50; area_cal = [0,x_max;0,y_max]; 
        sext = 0.1*h;  ext  =  5*h+sext;
        boundary_O      =   {'Rectangle',1, [-ext -sext;sext y_max+sext]};
        boundary_D      =   {'Rectangle',1, [x_max-sext 10+sext;x_max+ext y_max-10-sext]};
        boundary_D_2h   =   {'Rectangle',1, [x_max-sext-h 10+sext;x_max+ext y_max-10-sext]};
        boundary_H      =   {'Rectangle',5, [-ext -ext;x_max+ext -sext],[-ext y_max+sext;x_max+ext y_max+ext],[50-sext 10-sext;70+sext 30+sext],...
                            [x_max-sext y_max-10-sext;x_max+ext y_max+sext],[x_max-sext -sext;x_max+ext 10+sext]};
%         dir_data = ['C:\Users\HOWIE-PC\Desktop\Case study\MC simulations\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
        dir_data = ['C:\Users\Howie\Desktop\Case study\MC simulations\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
%         dir_data = ['C:\Users\HLiang\Desktop\Case study\MC simulations\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
    case 20221222
        Record_dt = 30;
        MC_std = 0.1;
        MC_mean = 1;
        MC_nseq = [20:20:80,100:100:2000];
        MC_max = 2000;
%         MC_nseq = MC_max;
        h = 1; 
        x_max = 100;  y_max = 50; area_cal = [0,x_max;0,y_max]; 
        sext = 0.1*h;  ext  =  5*h+sext;
        boundary_O      =   {'Rectangle',1, [-ext -sext;sext y_max+sext]};
        boundary_D      =   {'Rectangle',1, [x_max-sext 10+sext;x_max+ext y_max-10-sext]};
        boundary_D_2h   =   {'Rectangle',1, [x_max-sext-h 10+sext;x_max+ext y_max-10-sext]};
        boundary_H      =   {'Rectangle',5, [-ext -ext;x_max+ext -sext],[-ext y_max+sext;x_max+ext y_max+ext],[50-sext 10-sext;70+sext 30+sext],...
                            [x_max-sext y_max-10-sext;x_max+ext y_max+sext],[x_max-sext -sext;x_max+ext 10+sext]};
%         dir_data = ['C:\Users\HOWIE-PC\Desktop\Case study\MC simulations\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
        dir_data = ['C:\Users\HLiang\Desktop\Case study\MC simulations\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
    case 0
        Record_dt = 30;
        MC_std = 0.10;
        MC_mean = 1;
        MC_nseq = [20:20:80,100:100:1000];
        MC_max = 1000;
        MC_nseq = MC_max;
        h = 1; 
        x_max = 100;  y_max = 50; area_cal = [0,x_max;0,y_max]; 
        sext = 0.1*h;  ext  =  5*h+sext;
        boundary_O      =   {'Rectangle',1, [-ext -sext;sext y_max+sext]};
        boundary_D      =   {'Rectangle',1, [x_max-sext 10+sext;x_max+ext y_max-10-sext]};
        boundary_D_2h   =   {'Rectangle',1, [x_max-sext-h 10+sext;x_max+ext y_max-10-sext]};
        boundary_H      =   {'Rectangle',5, [-ext -ext;x_max+ext -sext],[-ext y_max+sext;x_max+ext y_max+ext],[50-sext 10-sext;70+sext 30+sext],...
                            [x_max-sext y_max-10-sext;x_max+ext y_max+sext],[x_max-sext -sext;x_max+ext 10+sext]};
%         dir_data = ['C:\Users\HOWIE-PC\Desktop\Case study\MC simulations\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
        dir_data = ['C:\Users\HLiang\Desktop\Case study\MC simulations cc\' 'H' num2str(h) '_AVE' num2str(MC_mean) '_STD' num2str(MC_std) '\'];
    case -1
        h = 0.2; 
        x_max = 2;  y_max = 2; area_cal = [0,2;0,2]; 
        sext = 0.1*h;  ext  =  5*h+sext;
        
        
end
end
%% Layout
function [bj_cal] = Layout_cal(boundary_O,boundary_D,boundary_D_2h,boundary_H,x,y)
    bj_cal = ones(length(y),length(x)).*9;
    bj_cal = Boundary_value_cal(x,y,bj_cal,boundary_D_2h,3);
    bj_cal = Boundary_value_cal(x,y,bj_cal,boundary_O,2);
    bj_cal = Boundary_value_cal(x,y,bj_cal,boundary_D,1);
    bj_cal = Boundary_value_cal(x,y,bj_cal,boundary_H,0);
end
  
%% Flow in
function F_in = F_in_cal(t,f1)
syms x
F_in = zeros(1,3);
switch true
    case t<=60
        F_in(2) = f1*(t / 60);
        F_in(3) = 0;
        F_in(1) = vpa(solve(1.34*x*exp(-0.09*(x^2))==F_in(2),x));
    case t>60
        F_in(2) = max(0,f1*((120-t) / 60));
        F_in(3) = 0;
        F_in(1) = vpa(solve(1.34*x*exp(-0.09*(x^2))==F_in(2),x));

end
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

%% RRMSE
function y = RRMSE_cal(x,order,n)
    y = sum(sum(sum( x.^order )))./n;
end
