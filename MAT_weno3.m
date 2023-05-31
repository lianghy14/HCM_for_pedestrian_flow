function Record_q1 = MAT_weno3(f_in,bj_ex3ODH,h,Record_dt)

%% Simulation period
tStart  =   0;
tEnd    =   300;
gfuns   =   functions_given;
pfuns   =   functions_plot;

%% Record data
Record_q1   =   zeros(size(bj_ex3ODH,1)-6,size(bj_ex3ODH,2)-6,ceil(tEnd/Record_dt)+1);

%% Initialization
tstep       =   0.01;
CFL         =   0.2;
cpu_tstart  =   tic;
t           =   tStart;
q1          =   zeros(size(bj_ex3ODH));
alpha_LF    =   5;
rn = 0;

%% TVD-RK scheme
while (t<=tEnd)

    if(mod(t,1)==0)
        if (mod(t,Record_dt)==0)
            cpu_dt              =   toc(cpu_tstart);
            rn = rn+1;
            Record_q1(:,:,rn)   =   q1(4:end-3,4:end-3);
        end        
        tstep           =   min([CFL*h/alpha_LF, CFL]);
    else
        tstep           =   min([CFL*h/alpha_LF, ceil(t)-t, CFL]);
    end
    %% TVD-RK3
    
    % STEP 1
    F_in = gfuns.F_in(t,f_in);
    [q1_RK1,alpha_LF1,nit,test]  =  F90_weno5_cc(q1,bj_ex3ODH,F_in,h,tstep);
    if any(any(isnan(q1_RK1))) || (nit>1999)
        error('FSM error!');
    end
    % STEP 2
    F_in = gfuns.F_in(t+tstep,f_in);
    [q1_RK2,alpha_LF2,nit,test]  =   F90_weno5_cc(q1_RK1,bj_ex3ODH,F_in,h,tstep);
    if any(any(isnan(q1_RK2))) || (nit>1999)
        error('FSM error!');
    end
    q1_RK2  =   3./4.*q1 + 1./4.*q1_RK2;

    % SREP 3
    F_in = gfuns.F_in(t+1/2*tstep,f_in);
    [q1_RK3,alpha_LF3,nit,test]  =   F90_weno5_cc(q1_RK2,bj_ex3ODH,F_in,h,tstep);
    if any(any(isnan(q1_RK3))) || (nit>1999)
        error('FSM error!');
    end
    q1      =   1./3.*q1 + 2./3.*q1_RK3;

    alpha_LF = max([alpha_LF1,alpha_LF2,alpha_LF3]);
    
    if (mod(t,Record_dt)==0)
        fprintf('HCM. T = %d. CPUt is %.2f. FSM iterations is %d.\n',t,cpu_dt,nit);
    end   
    t = t+tstep;

end

%% Error

if sum(sum(Record_q1(:,:,rn)))*h*h > 0.1
    error('tEnd is not enough!');
end
end