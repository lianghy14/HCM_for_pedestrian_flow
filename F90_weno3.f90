#include "fintrf.h" !必须有的头文件,里面有mxGetPr, mxGetM, mxGetN,mxCreateDoubleMatrix等函数的申明

Subroutine mexFunction(OutSum,OutVar,InSum,InVar)!函数接口名称必须为mexFunction,
    !OutSum:输出参数个数, OutVar:输出参数数组指针, InSum:输入参数个数, InVar:输入参数数组指针

    Integer InSum,OutSum
    mwPointer InVar(*),OutVar(*)                           !mwPointer专门用于表示指针变量,这个不能随意用Integer代替
    mwPointer mxGetPr, mxGetM, mxGetN, mxCreateDoubleMatrix !这个对返回指针函数的再次申明,
    integer, parameter              ::  acc = 8
    Real(kind = acc),Allocatable    ::  q1(:,:),q2(:,:),boundary(:,:)
    real(kind = acc)                ::  h,tstep,f_in(2),alpha_LF
    Integer                         ::  ny,nx

    If(InSum/=5)Then
        call mexErrMsgIdAndTxt('MATLAB:InputTooBig','Less than 5 inputs.')
        Return
    EndIf

    ny = mxGetM(InVar(1))!获取第1个输入参数的行数
    nx = mxGetN(InVar(1))!获取第1个输入参数的列数
    Allocate(q1(ny,nx),q2(ny,nx),boundary(ny,nx))

    Call mxCopyPtrToReal8(mxGetPr(InVar(1)),q1,ny*nx)
    Call mxCopyPtrToReal8(mxGetPr(InVar(2)),boundary,ny*nx)
    Call mxCopyPtrToReal8(mxGetPr(InVar(3)),f_in,2)
    Call mxCopyPtrToReal8(mxGetPr(InVar(4)),h,1)
    Call mxCopyPtrToReal8(mxGetPr(InVar(5)),tstep,1)

    Call WENO3_euler(q1,q2,boundary,f_in,h,tstep,ny,nx,alpha_LF)

    OutVar(1) = mxCreateDoubleMatrix(ny,nx,0)!给返回参数分配内存
    OutVar(2) = mxCreateDoubleMatrix(1,1,0)!给返回参数分配内存

    Call mxCopyReal8ToPtr(q2,mxGetPr(OutVar(1)),ny*nx)
    Call mxCopyReal8ToPtr(alpha_LF,mxGetPr(OutVar(2)),1)

    DeAllocate(q1,q2,boundary)!释放临时分配的内存
    Return

End SubRoutine


SubRoutine WENO3_euler(q1,q2,bj_real,f_in,h,tstep,ny,nx,alpha_LF)

    implicit none
    integer, parameter          ::  acc = 8
    real(kind=acc),parameter    ::  vel_f = 1.34, gama = -0.09
    integer, Intent(In)         ::  ny, nx
    real(kind=acc), Intent(In)  ::  q1(ny,nx), f_in(2), h, tstep, bj_real(ny,nx)
    real(kind=acc), Intent(Out) ::  q2(ny,nx), alpha_LF
    real(kind=acc)              ::  cost(ny,nx),potential(ny,nx),p_x(ny,nx),p_y(ny,nx),vel_mag(ny,nx),&
                                    &vel_x(ny,nx),vel_y(ny,nx),flow_x(ny,nx),flow_y(ny,nx),q_xex(ny,nx),&
                                    &q_yex(ny,nx),flow_xp(ny,nx),flow_yp(ny,nx),flow_xm(ny,nx),flow_ym(ny,nx),&
                                    &flow_xcom(ny,nx),flow_ycom(ny,nx)
    integer                     ::  i, j, n, boundary(ny,nx), nit

    boundary    =   NINT(bj_real)

    ! Obtain velosity & cost
    vel_mag =   vel_f * exp( gama*(max(q1,0.0)**2) )
    cost    =   1.0 / vel_mag + 0.02 * (max(q1,0.0))**2

    ! Obtain potential gradients
    call F90_fsmgod(cost,boundary,potential,p_x,p_y,h,ny,nx,nit)

    ! Obtain velosity & flow on grids
    vel_x   =   - p_x * vel_mag
    vel_y   =   - p_y * vel_mag
    flow_x  =   vel_x * q1
    flow_y  =   vel_y * q1
    q_xex   =   q1
    q_yex   =   q1

    ! ghost points
    ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
    ! 2
    do i = 1,ny,1
        do j = 1,nx,1
            if (boundary(i,j)>2) then
                if (boundary(i,j+1)==2) then
                    flow_x(i,j+1:j+2:1)     = f_in(1)
                    q_xex(i,j+1:j+2:1)      = q_xex(i,j)
                end if
                if (boundary(i,j-1)==2) then
                    flow_x(i,j-2:j-1:1)     = f_in(1)
                    q_xex(i,j-2:j-1:1)      = q_xex(i,j)
                end if
                if (boundary(i+1,j)==2) then
                    flow_y(i+1:i+2:1,j)     = f_in(2)
                    q_yex(i+1:i+2:1,j)      = q_yex(i,j)
                end if
                if (boundary(i-1,j)==2) then
                    flow_y(i-2:i-1:1,j)     = f_in(2)
                    q_yex(i-2:i-1:1,j)      = q_yex(i,j)
                end if
            end if
        end do
    end do
    ! 1
    do i = 1,ny,1
        do j = 1,nx,1
            if (boundary(i,j)>2) then
                if (boundary(i,j+1)==1) then
                    flow_x(i,j+1:j+2:1)     = flow_x(i,j)
                    q_xex(i,j+1:j+2:1)      = q_xex(i,j)
                end if
                if (boundary(i,j-1)==1) then
                    flow_x(i,j-2:j-1:1)     = flow_x(i,j)
                    q_xex(i,j-2:j-1:1)      = q_xex(i,j)
                end if
                if (boundary(i+1,j)==1) then
                    flow_y(i+1:i+2:1,j)     = flow_y(i,j)
                    q_yex(i+1:i+2:1,j)      = q_yex(i,j)
                end if
                if (boundary(i-1,j)==1) then
                    flow_y(i-2:i-1:1,j)     = flow_y(i,j)
                    q_yex(i-2:i-1:1,j)      = q_yex(i,j)
                end if
            end if
        end do
    end do
    ! 0
    do i = 1,ny,1
        do j = 1,nx,1
            if (boundary(i,j)>2) then
                if (boundary(i,j+1)==0) then
                    flow_x(i,j+1)           =   - flow_x(i,j)
                    flow_x(i,j+2)           =   - flow_x(i,j-1)
                    q_xex(i,j+1:j+2:1)      =   q_xex(i,j:j-1:-1)
                end if
                if (boundary(i,j-1)==0) then
                    flow_x(i,j-1)           =   - flow_x(i,j)
                    flow_x(i,j-2)           =   - flow_x(i,j+1)
                    q_xex(i,j-1:j-2:-1)     =   q_xex(i,j:j+1:1)
                end if
                if (boundary(i+1,j)==0) then
                    flow_y(i+1,j)           =   - flow_y(i,j)
                    flow_y(i+2,j)           =   - flow_y(i-1,j)
                    q_yex(i+1:i+2:1,j)      =   q_yex(i:i-1:-1,j)
                end if
                if (boundary(i-1,j)==0) then
                    flow_y(i-1,j)           =   - flow_y(i,j)
                    flow_y(i-2,j)           =   - flow_y(i+1,j)
                    q_yex(i-1:i-2:-1,j)     =   q_yex(i:i+1:1,j)
                end if
            end if
        end do
    end do

    ! Flow x/y upwind/downwind
    alpha_LF    =   max(maxval(maxval(abs(vel_x),2),1),maxval(maxval(abs(vel_y),1),1))
    flow_xp     =   0.5*flow_x + 0.5*spread(maxval(abs(vel_x),2),2,nx)*q_xex
    flow_xm     =   0.5*flow_x - 0.5*spread(maxval(abs(vel_x),2),2,nx)*q_xex
    flow_yp     =   0.5*flow_y + 0.5*spread(maxval(abs(vel_y),1),1,ny)*q_yex
    flow_ym     =   0.5*flow_y - 0.5*spread(maxval(abs(vel_y),1),1,ny)*q_yex
    
    call WENO_Res_3(flow_xp, 2, 1,ny,nx)
    call WENO_Res_3(flow_xm, 2,-1,ny,nx)
    call WENO_Res_3(flow_yp, 1,-1,ny,nx)
    call WENO_Res_3(flow_ym, 1, 1,ny,nx)
    
    ! Combination
    flow_xcom = flow_xp + eoshift(flow_xm, 1,0.0,2)
    flow_ycom = flow_yp + eoshift(flow_ym,-1,0.0,1)
    
    ! boundary value
    ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
    do i = 1,ny,1
        do j = 1,nx-1,1
            if ((boundary(i,j)==2).OR.(boundary(i,j+1)==2)) then
                flow_xcom(i,j) = f_in(1)
            end if
            if ((boundary(i,j)==0).OR.(boundary(i,j+1)==0)) then
                flow_xcom(i,j) = 0
            end if
        end do
    end do
    do i = 2,ny,1
        do j = 1,nx,1
            if ((boundary(i-1,j)==2).OR.(boundary(i,j)==2)) then
                flow_ycom(i,j) = f_in(2)
            end if
            if ((boundary(i-1,j)==0).OR.(boundary(i,j)==0)) then
                flow_ycom(i,j) = 0
            end if
        end do
    end do

    ! Divergence
    q2  =   q1  -1.0/h *tstep * (flow_xcom - eoshift(flow_xcom,-1,0.0,2) +&
                                &flow_ycom - eoshift(flow_ycom,1,0.0,1))
    q2((/1:2,ny-1:ny/),:) = 0.0
    q2(:,(/1:2,nx-1:nx/)) = 0.0

    return
end subroutine WENO3_euler



subroutine WENO_Res_3(flow, dim, dir, ny, nx)
    implicit none
    integer, parameter          ::  acc = 8
    real, parameter             ::  epsilon = 1.0d-6
    integer, Intent(In)         ::  ny, nx, dim, dir
    real(kind=acc)              ::  flow(ny,nx)
    real(kind=acc)              ::  flow_shift(3,ny,nx), up(2, ny,nx), b(2, ny,nx), g(2) = [1./3, 2./3],&
                                    &w(2,ny,nx), sum_w(ny,nx)
    ! Set Variables
    flow_shift(1,:,:) = eoshift(flow,dir,0.0,dim)
    flow_shift(2,:,:) = flow
    flow_shift(3,:,:) = eoshift(flow,-dir,0.0,dim)
    ! Reconstruction Polynomials
    up(1,:,:) = -1./2*flow_shift(1,:,:) + 3./2*flow_shift(2,:,:)
    up(2,:,:) =  1./2*flow_shift(2,:,:) + 1./2*flow_shift(3,:,:)
    ! Smooth parameters
    b(1,:,:) = (flow_shift(1,:,:)-flow_shift(2,:,:))**2
    b(2,:,:) = (flow_shift(2,:,:)-flow_shift(3,:,:))**2
    ! weigths
    w(1,:,:) = g(1) / (epsilon+b(1,:,:)**2)
    w(2,:,:) = g(2) / (epsilon+b(2,:,:)**2)
    sum_w = w(1,:,:) + w(2,:,:)
    ! Non-linear weigths
    w(1,:,:) = w(1,:,:) / sum_w
    w(2,:,:) = w(2,:,:) / sum_w
    ! WENO polynomial
    flow = w(1,:,:)*up(1,:,:) + w(2,:,:)*up(2,:,:)
    return
end subroutine WENO_Res_3


subroutine F90_fsmgod(cost,boundary,potential,p_x,p_y,h,ny,nx,nit)
    implicit none
    integer, parameter          ::  acc = 8
    real, parameter             ::  sigma = 1.0d-9,  p_max = 1.0d12
    integer, Intent(In)         ::  ny, nx,boundary(ny,nx)
    real(kind=acc), Intent(In)  ::  cost(ny,nx),h
    real(kind=acc), Intent(Out) ::  potential(ny,nx),p_x(ny,nx),p_y(ny,nx)
    integer                     ::  i, j, n, nit_c = 1, Gsit(4,2,2), GSi(4,2),nit
    real(kind=acc)              ::  p_old(ny,nx), pot_x, pot_y, p_t, p_temp

    Gsit    = reshape(  (/1 , nx, nx, 1 ,&
                            ny, ny, 1 , 1 ,&
                            nx, 1 , 1 , nx,&
                            1 , 1 , ny, ny/),(/4,2,2/))
    GSi     = reshape(  (/ 1, -1, -1,  1,&
                        -1, -1,  1,  1/),(/4,2/))
    p_x = 0.0
    p_y = 0.0
    potential = p_max
    potential = merge(0.0,potential,boundary==1)
    p_old = 0.0
    nit = 0
    ! Obtain potential by Fast Godunov
    do nit = 1,1000,1
        if ( ( sum(abs(potential - p_old)) ) >= sigma ) then
            p_old = potential
            do n = 1,4,1
                do j = GSit(n,2,1), GSit(n,2,2), GSi(n,2)
                    do i = GSit(n,1,1), GSit(n,1,2), GSi(n,1)

                        ! 0 is bound_H, 1 is bound_D, 2 is bound_O
                        if (boundary(j,i)>=3) then
                            
                            pot_x = min(potential(j,max(1,i-1)), potential(j,min(nx,i+1)))
                            pot_y = min(potential(max(1,j-1),i), potential(min(ny,j+1),i))
                            if (abs(pot_x-pot_y)>=(cost(j,i))) then
                                p_temp = (min(pot_x, pot_y) + cost(j,i)*h)
                            else 
                                p_temp = (pot_x+pot_y + (2*(cost(j,i)*h)**2-(pot_x-pot_y)**2)**0.5) / 2.0
                            end if
                            potential(j,i) = min(potential(j,i),p_temp)

                            
                        end if 

                    end do
                end do
            end do
        else
            do j = 1,ny,1
                do i = 1,nx,1
                    ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D

                    ! original
                    if (boundary(j,i)>2) then
                            
                        pot_x = min(potential(j,max(1,i-1)), potential(j,min(nx,i+1)))
                        pot_y = min(potential(max(1,j-1),i), potential(min(ny,j+1),i))
                        if (abs(pot_x-pot_y)>=(cost(j,i))) then
                            p_temp = (min(pot_x, pot_y) + cost(j,i)*h)
                        else 
                            p_temp = (pot_x+pot_y + (2*(cost(j,i)*h)**2-(pot_x-pot_y)**2)**0.5) / 2.0
                        end if
                        
                        if (potential(j,max(1,i-1)) >= potential(j,min(nx,i+1))) then
                            p_x(j,i) = (potential(j,min(nx,i+1)) - potential(j,i)) / h
                        else
                            p_x(j,i) = (potential(j,i) - potential(j,max(1,i-1))) / h
                        end if
                        if (potential(min(ny,j+1),i) >= potential(max(1,j-1),i)) then
                            p_y(j,i) = (potential(max(1,j-1),i) - potential(j,i)) / h
                        else
                            p_y(j,i) = (potential(j,i) - potential(min(ny,j+1),i)) / h
                        end if
                        potential(j,i) = min(potential(j,i),p_temp)

                    end if 

                    ! central
                    ! if (boundary(j,i)>2) then

                    !     if (((boundary(j,i-1)/=2).AND.(boundary(j,i-1)/=0)) .AND. ((boundary(j,i+1)/=2).AND.(boundary(j,i+1)/=0))) then
                    !         p_x(j,i) = (potential(j,i+1) - potential(j,i-1)) / h / 2
                    !     else
                    !         if (((boundary(j,i-1)==2).OR.(boundary(j,i-1)==0)) .AND. ((boundary(j,i+1)/=2).AND.(boundary(j,i+1)/=0))) then
                    !             p_x(j,i) = (potential(j,i+1) - potential(j,i)) / h
                    !         end if
                    !         if (((boundary(j,i+1)==2).OR.(boundary(j,i+1)==0)) .AND. ((boundary(j,i-1)/=2).AND.(boundary(j,i-1)/=0))) then
                    !             p_x(j,i) = (potential(j,i) - potential(j,i-1)) / h
                    !         end if
                    !     end if
                        
                    !     if (((boundary(j-1,i)/=2).AND.(boundary(j-1,i)/=0)) .AND. ((boundary(j+1,i)/=2).AND.(boundary(j+1,i)/=0))) then
                    !         p_y(j,i) = (potential(j-1,i) - potential(j+1,i)) / h / 2
                    !     else
                    !         if (((boundary(j-1,i)==2).OR.(boundary(j-1,i)==0)) .AND. ((boundary(j+1,i)/=2).AND.(boundary(j+1,i)/=0))) then
                    !             p_y(j,i) = (potential(j,i) - potential(j+1,i)) / h
                    !         end if
                    !         if (((boundary(j+1,i)==2).OR.(boundary(j+1,i)==0)) .AND. ((boundary(j-1,i)/=2).AND.(boundary(j-1,i)/=0))) then
                    !             p_y(j,i) = (potential(j-1,i) - potential(j,i)) / h
                    !         end if
                    !     end if
                    ! else
                    !     if (boundary(j,i)==1) then
                    !         if (boundary(j,i-1)>2) then
                    !             p_x(j,i) = (potential(j,i) - potential(j,i-1)) / h
                    !         end if
                    !         if (boundary(j,i+1)>2) then
                    !             p_x(j,i) = (potential(j,i+1) - potential(j,i)) / h
                    !         end if
                    !         if (boundary(j+1,i)>2) then
                    !             p_y(j,i) = (potential(j,i) - potential(j+1,i)) / h
                    !         end if
                    !         if (boundary(j-1,i)>2) then
                    !             p_y(j,i) = (potential(j-1,i) - potential(j,i)) / h
                    !         end if
                    !     end if
                    ! end if

                end do
            end do
            exit
        end if
        
    end	do

    return
end subroutine F90_fsmgod