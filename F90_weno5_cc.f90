#include "fintrf.h" !必须有的头文件,里面有mxGetPr, mxGetM, mxGetN,mxCreateDoubleMatrix等函数的申明

Subroutine mexFunction(OutSum,OutVar,InSum,InVar)!函数接口名称必须为mexFunction,
    !OutSum:输出参数个数, OutVar:输出参数数组指针, InSum:输入参数个数, InVar:输入参数数组指针

    Integer InSum,OutSum
    mwPointer InVar(*),OutVar(*)                           !mwPointer专门用于表示指针变量,这个不能随意用Integer代替
    mwPointer mxGetPr, mxGetM, mxGetN, mxCreateDoubleMatrix !这个对返回指针函数的再次申明,
    integer, parameter              ::  acc = 8
    Real(kind = acc),Allocatable    ::  q1(:,:),q2(:,:),boundary(:,:),potential(:,:)
    real(kind = acc)                ::  h,tstep,f_in(3),alpha_LF,nit_r
    Integer                         ::  ny,nx,nit

    If(InSum/=5)Then
        call mexErrMsgIdAndTxt('MATLAB:InputTooBig','Less than 5 inputs.')
        Return
    EndIf

    ny = mxGetM(InVar(1))!获取第1个输入参数的行数
    nx = mxGetN(InVar(1))!获取第1个输入参数的列数
    Allocate(q1(ny,nx),q2(ny,nx),boundary(ny,nx),potential(ny,nx))

    Call mxCopyPtrToReal8(mxGetPr(InVar(1)),q1,ny*nx)
    Call mxCopyPtrToReal8(mxGetPr(InVar(2)),boundary,ny*nx)
    Call mxCopyPtrToReal8(mxGetPr(InVar(3)),f_in,3)
    Call mxCopyPtrToReal8(mxGetPr(InVar(4)),h,1)
    Call mxCopyPtrToReal8(mxGetPr(InVar(5)),tstep,1)

    Call WENO5_euler(q1,q2,boundary,f_in,h,tstep,ny,nx,alpha_LF,nit,potential)
    nit_r = real(nit)

    OutVar(1) = mxCreateDoubleMatrix(ny,nx,0)!给返回参数分配内存
    OutVar(2) = mxCreateDoubleMatrix(1,1,0)!给返回参数分配内存
    OutVar(3) = mxCreateDoubleMatrix(1,1,0)!给返回参数分配内存
    OutVar(4) = mxCreateDoubleMatrix(ny,nx,0)!给返回参数分配内存

    Call mxCopyReal8ToPtr(q2,mxGetPr(OutVar(1)),ny*nx)
    Call mxCopyReal8ToPtr(alpha_LF,mxGetPr(OutVar(2)),1)
    Call mxCopyReal8ToPtr(nit_r,mxGetPr(OutVar(3)),1)
    Call mxCopyReal8ToPtr(potential,mxGetPr(OutVar(4)),ny*nx)

    DeAllocate(q1,q2,boundary,potential)!释放临时分配的内存

    Return

End SubRoutine


SubRoutine WENO5_euler(q1,q2,bj_real,f_in,h,tstep,ny,nx,alpha_LF,nit,test)

    implicit none
    integer, parameter          ::  acc = 8
    real(kind=acc),parameter    ::  vel_f = 1.34d0, gama = -0.09d0
    integer, Intent(In)         ::  ny, nx
    real(kind=acc), Intent(In)  ::  q1(ny,nx), f_in(3), h, tstep, bj_real(ny,nx)
    real(kind=acc), Intent(Out) ::  q2(ny,nx), alpha_LF, test(ny,nx)
    real(kind=acc)              ::  cost(ny,nx),potential(ny,nx),p_x(ny,nx),p_y(ny,nx),vel_mag(ny,nx),&
                                    &vel_x(ny,nx),vel_y(ny,nx),flow_x(ny,nx),flow_y(ny,nx),q_xex(ny,nx),&
                                    &q_yex(ny,nx),flow_xp(ny,nx),flow_yp(ny,nx),flow_xm(ny,nx),flow_ym(ny,nx),&
                                    &flow_xcom(ny,nx),flow_ycom(ny,nx),q_ex(ny,nx)
    integer                     ::  i, j, n, boundary(ny,nx), nit

    boundary    =   NINT(bj_real)

    ! Obtain velosity & cost
    ! Obtain potential gradients
    q_ex   =   q1
    do i = 1,ny,1
        do j = 1,nx,1
            if (boundary(i,j)>2) then
                if (boundary(i,j+1)==2) then
                    q_ex(i,j+1:j+3:1)      = f_in(1)
                end if
                if (boundary(i,j-1)==2) then
                    q_ex(i,j-3:j-1:1)      = f_in(1)
                end if
                if (boundary(i+1,j)==2) then
                    q_ex(i+1:i+3:1,j)      = f_in(1)
                end if
                if (boundary(i-1,j)==2) then
                    q_ex(i-3:i-1:1,j)      = f_in(1)
                end if
            end if
        end do
    end do
    vel_mag =   vel_f * exp( gama*(max(q_ex,0d0)**2) )
    cost    =   1d0
    vel_mag =   merge(0d0,vel_mag, (boundary == 0).OR.(boundary == 1))
    ! call F90_fsmweno3(cost,boundary,potential,p_x,p_y,h,ny,nx,nit,test)
    call F90_fsmgod(cost,boundary,potential,p_x,p_y,h,ny,nx,nit,test)
    ! test = potential

    ! p_x = p_x / max(1.0d-9,cost)
    ! p_y = p_y / max(1.0d-9,cost)
    p_x = p_x / max(1d-9,(p_x**2 + p_y**2)**0.5)
    p_y = p_y / max(1d-9,(p_x**2 + p_y**2)**0.5)

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
        ! do i = 1,ny,1
        !     do j = 1,nx,1
        !         if (boundary(i,j)>2) then
        !             if (boundary(i,j+1)==2) then
        !                 flow_x(i,j+1:j+3:1)     = f_in(2)
        !                 q_xex(i,j+1:j+3:1)      = f_in(1)
        !             end if
        !             if (boundary(i,j-1)==2) then
        !                 flow_x(i,j-3:j-1:1)     = f_in(2)
        !                 q_xex(i,j-3:j-1:1)      = f_in(1)
        !             end if
        !             if (boundary(i+1,j)==2) then
        !                 flow_y(i+1:i+3:1,j)     = f_in(3)
        !                 q_yex(i+1:i+3:1,j)      = f_in(1)
        !             end if
        !             if (boundary(i-1,j)==2) then
        !                 flow_y(i-3:i-1:1,j)     = f_in(3)
        !                 q_yex(i-3:i-1:1,j)      = f_in(1)
        !             end if
        !         end if
        !     end do
        ! end do
        ! do i = 1,ny,1
        !     do j = 1,nx,1
        !         if (boundary(i,j)>2) then
        !             if (boundary(i,j+1)==2) then
        !                 flow_x(i,j+1)     = 3*flow_x(i,j) - 3*flow_x(i,j-1) + flow_x(i,j-2)
        !                 flow_x(i,j+2)     = 3*flow_x(i,j+1) - 3*flow_x(i,j) + flow_x(i,j-1)
        !                 flow_x(i,j+3)     = 3*flow_x(i,j+2) - 3*flow_x(i,j+1) + flow_x(i,j)
        !                 q_xex(i,j+1:j+3:1)      = q_xex(i,j)
        !             end if
        !             if (boundary(i,j-1)==2) then
        !                 flow_x(i,j-1)     = 3*flow_x(i,j) - 3*flow_x(i,j+1) + flow_x(i,j+2)
        !                 flow_x(i,j-2)     = 3*flow_x(i,j-1) - 3*flow_x(i,j) + flow_x(i,j+1)
        !                 flow_x(i,j-3)     = 3*flow_x(i,j-2) - 3*flow_x(i,j-1) + flow_x(i,j)
        !                 q_xex(i,j-3:j-1:1)      = q_xex(i,j)
        !             end if
        !             if (boundary(i+1,j)==2) then
        !                 flow_y(i+1,j)     = 3*flow_y(i,j) - 3*flow_y(i-1,j) + flow_y(i-2,j)
        !                 flow_y(i+2,j)     = 3*flow_y(i+1,j) - 3*flow_y(i,j) + flow_y(i-1,j)
        !                 flow_y(i+3,j)     = 3*flow_y(i+2,j) - 3*flow_y(i+1,j) + flow_y(i,j)
        !                 q_yex(i+1:i+3:1,j)      = q_yex(i,j)
        !             end if
        !             if (boundary(i-1,j)==2) then
        !                 flow_y(i-1,j)     = 3*flow_y(i,j) - 3*flow_y(i+1,j) + flow_y(i+2,j)
        !                 flow_y(i-2,j)     = 3*flow_y(i-1,j) - 3*flow_y(i,j) + flow_y(i+1,j)
        !                 flow_y(i-3,j)     = 3*flow_y(i-2,j) - 3*flow_y(i-1,j) + flow_y(i,j)
        !                 q_yex(i-3:i-1:1,j)      = q_yex(i,j)
        !             end if
        !         end if
        !     end do
        ! end do
        ! do i = 1,ny,1
        !     do j = 1,nx,1
        !         if (boundary(i,j)>2) then
        !             if (boundary(i,j+1)==2) then
        !                 flow_x(i,j+1:j+3:1)     = flow_x(i,j)
        !                 q_xex(i,j+1:j+3:1)      = q_xex(i,j)
        !             end if
        !             if (boundary(i,j-1)==2) then
        !                 flow_x(i,j-3:j-1:1)     = flow_x(i,j)
        !                 q_xex(i,j-3:j-1:1)      = q_xex(i,j)
        !             end if
        !             if (boundary(i+1,j)==2) then
        !                 flow_y(i+1:i+3:1,j)     = flow_y(i,j)
        !                 q_yex(i+1:i+3:1,j)      = q_yex(i,j)
        !             end if
        !             if (boundary(i-1,j)==2) then
        !                 flow_y(i-3:i-1:1,j)     = flow_y(i,j)
        !                 q_yex(i-3:i-1:1,j)      = q_yex(i,j)
        !             end if
        !         end if
        !     end do
        ! end do

        ! 1
        ! do i = 1,ny,1
        !     do j = 1,nx,1
        !         if (boundary(i,j)>2) then
        !             if (boundary(i,j+1)==1) then
        !                 flow_x(i,j+1)       = 3*flow_x(i,j) - 3*flow_x(i,j-1) + flow_x(i,j-2)
        !                 flow_x(i,j+2)       = 3*flow_x(i,j+1) - 3*flow_x(i,j) + flow_x(i,j-1)
        !                 flow_x(i,j+3)       = 3*flow_x(i,j+2) - 3*flow_x(i,j+1) + flow_x(i,j)
        !                 q_xex(i,j+1:j+3:1)  = q_xex(i,j)
        !             end if
        !             if (boundary(i,j-1)==1) then
        !                 flow_x(i,j-1)       = 3*flow_x(i,j) - 3*flow_x(i,j+1) + flow_x(i,j+2)
        !                 flow_x(i,j-2)       = 3*flow_x(i,j-1) - 3*flow_x(i,j) + flow_x(i,j+1)
        !                 flow_x(i,j-3)       = 3*flow_x(i,j-2) - 3*flow_x(i,j-1) + flow_x(i,j)
        !                 q_xex(i,j-3:j-1:1)      = q_xex(i,j)
        !             end if
        !             if (boundary(i+1,j)==1) then
        !                 flow_y(i+1,j)       = 3*flow_y(i,j) - 3*flow_y(i-1,j) + flow_y(i-2,j)
        !                 flow_y(i+2,j)       = 3*flow_y(i+1,j) - 3*flow_y(i,j) + flow_y(i-1,j)
        !                 flow_y(i+3,j)       = 3*flow_y(i+2,j) - 3*flow_y(i+1,j) + flow_y(i,j)
        !                 q_yex(i+1:i+3:1,j)      = q_yex(i,j)
        !             end if
        !             if (boundary(i-1,j)==1) then
        !                 flow_y(i-1,j)       = 3*flow_y(i,j) - 3*flow_y(i+1,j) + flow_y(i+2,j)
        !                 flow_y(i-2,j)       = 3*flow_y(i-1,j) - 3*flow_y(i,j) + flow_y(i+1,j)
        !                 flow_y(i-3,j)       = 3*flow_y(i-2,j) - 3*flow_y(i-1,j) + flow_y(i,j)
        !                 q_yex(i-3:i-1:1,j)      = q_yex(i,j)
        !             end if
        !         end if
        !     end do
        ! end do
        ! do i = 1,ny,1
        !     do j = 1,nx,1
        !         if (boundary(i,j)>2) then
        !             if (boundary(i,j+1)==1) then
        !                 q_xex(i,j+1:j+3:1)      = q_xex(i,j)
        !                 flow_x(i,j+1:j+3:1)     = vel_f*q_xex(i,j)
        !             end if
        !             if (boundary(i,j-1)==1) then
        !                 q_xex(i,j-3:j-1:1)      = q_xex(i,j)
        !                 flow_x(i,j-3:j-1:1)     = - vel_f*q_xex(i,j)
        !             end if
        !             if (boundary(i+1,j)==1) then
        !                 q_yex(i+1:i+3:1,j)      = q_yex(i,j)
        !                 flow_y(i+1:i+3:1,j)     = vel_f*q_yex(i,j)
        !             end if
        !             if (boundary(i-1,j)==1) then
        !                 q_yex(i-3:i-1:1,j)      = q_yex(i,j)
        !                 flow_y(i-3:i-1:1,j)     = - vel_f*q_yex(i,j)
        !             end if
        !         end if
        !     end do
        ! end do
        do i = 1,ny,1
            do j = 1,nx,1
                if (boundary(i,j)>2) then
                    if (boundary(i,j+1)==1) then
                        q_xex(i,j+1:j+3:1)      = q_xex(i,j)
                        flow_x(i,j+1:j+3:1)     = flow_x(i,j)
                    end if
                    if (boundary(i,j-1)==1) then
                        q_xex(i,j-3:j-1:1)      = q_xex(i,j)
                        flow_x(i,j-3:j-1:1)     = flow_x(i,j)
                    end if
                    if (boundary(i+1,j)==1) then
                        q_yex(i+1:i+3:1,j)      = q_yex(i,j)
                        flow_y(i+1:i+3:1,j)     = flow_y(i,j)
                    end if
                    if (boundary(i-1,j)==1) then
                        q_yex(i-3:i-1:1,j)      = q_yex(i,j)
                        flow_y(i-3:i-1:1,j)     = flow_y(i,j)
                    end if
                end if
            end do
        end do

    ! Flow x/y upwind/downwind
    alpha_LF    =   max(maxval(maxval(abs(vel_mag),2),1),maxval(maxval(abs(vel_mag),1),1))
    flow_xp     =   5d-1*flow_x + 5d-1*spread(maxval(abs(vel_mag),2),2,nx)*q_xex
    flow_xm     =   5d-1*flow_x - 5d-1*spread(maxval(abs(vel_mag),2),2,nx)*q_xex
    flow_yp     =   5d-1*flow_y + 5d-1*spread(maxval(abs(vel_mag),1),1,ny)*q_yex
    flow_ym     =   5d-1*flow_y - 5d-1*spread(maxval(abs(vel_mag),1),1,ny)*q_yex

    call WENO_ghost(flow_xp,flow_yp,boundary,ny,nx)
    call WENO_ghost(flow_xm,flow_ym,boundary,ny,nx)

    call WENO_Res_5(flow_xp, 2, 1,ny,nx)
    call WENO_Res_5(flow_xm, 2,-1,ny,nx)
    call WENO_Res_5(flow_yp, 1,-1,ny,nx)
    call WENO_Res_5(flow_ym, 1, 1,ny,nx)
    
    ! Combination
    flow_xcom = flow_xp + eoshift(flow_xm, 1,0d0,2)
    flow_ycom = flow_yp + eoshift(flow_ym,-1,0d0,1)

    ! boundary value
    ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
    do i = 1,ny,1
        do j = 1,nx-1,1
            if ((boundary(i,j)==2).OR.(boundary(i,j+1)==2)) then
                flow_xcom(i,j) = f_in(2)
            end if
            if ((boundary(i,j)==0).OR.(boundary(i,j+1)==0)) then
                flow_xcom(i,j) = 0d0
            end if
        end do
    end do
    do i = 2,ny,1
        do j = 1,nx,1
            if ((boundary(i-1,j)==2).OR.(boundary(i,j)==2)) then
                flow_ycom(i,j) = f_in(3)
            end if
            if ((boundary(i-1,j)==0).OR.(boundary(i,j)==0)) then
                flow_ycom(i,j) = 0d0
            end if
        end do
    end do

    ! Divergence
    q2  =   q1  -1d0/h *tstep * (flow_xcom - eoshift(flow_xcom,-1,0d0,2) +&
                                &flow_ycom - eoshift(flow_ycom, 1,0d0,1))
    q2((/1:3,ny-2:ny/),:) = 0d0
    q2(:,(/1:3,nx-2:nx/)) = 0d0

    return
end subroutine WENO5_euler

subroutine WENO_ghost(flow_x,flow_y,boundary,ny,nx)
    implicit none
    integer, parameter          ::  acc = 8
    integer, Intent(In)         ::  ny, nx, boundary(ny,nx)
    real(kind=acc)              ::  flow_x(ny,nx),flow_y(ny,nx)
    integer                     ::  i,j
    ! ghost points
    ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
    ! 2
    ! do i = 1,ny,1
    !     do j = 1,nx,1
    !         if (boundary(i,j)>2) then
    !             if (boundary(i,j+1)==2) then
    !                 flow_x(i,j+1)       = 2d0*flow_x(i,j) - flow_x(i,j-1)
    !                 flow_x(i,j+2)       = 2d0*flow_x(i,j+1) - flow_x(i,j)
    !                 flow_x(i,j+3)       = 2d0*flow_x(i,j+2) - flow_x(i,j+1)
    !             end if
    !             if (boundary(i,j-1)==2) then
    !                 flow_x(i,j-1)       = 2d0*flow_x(i,j) - flow_x(i,j+1)
    !                 flow_x(i,j-2)       = 2d0*flow_x(i,j-1) - flow_x(i,j)
    !                 flow_x(i,j-3)       = 2d0*flow_x(i,j-2) - flow_x(i,j-1)
    !             end if
    !             if (boundary(i+1,j)==2) then
    !                 flow_y(i+1,j)       = 2d0*flow_y(i,j) - flow_y(i-1,j)
    !                 flow_y(i+2,j)       = 2d0*flow_y(i+1,j) - flow_y(i,j)
    !                 flow_y(i+3,j)       = 2d0*flow_y(i+2,j) - flow_y(i+1,j)
    !             end if
    !             if (boundary(i-1,j)==2) then
    !                 flow_y(i-1,j)       = 2d0*flow_y(i,j) - flow_y(i+1,j)
    !                 flow_y(i-2,j)       = 2d0*flow_y(i-1,j) - flow_y(i,j)
    !                 flow_y(i-3,j)       = 2d0*flow_y(i-2,j) - flow_y(i-1,j)
    !             end if
    !         end if
    !     end do
    ! end do
    ! do i = 1,ny,1
    !     do j = 1,nx,1
    !         if (boundary(i,j)>2) then
    !             if (boundary(i,j+1)==2) then
    !                 flow_x(i,j+1)       = 3*flow_x(i,j) - 3*flow_x(i,j-1) + flow_x(i,j-2)
    !                 flow_x(i,j+2)       = 3*flow_x(i,j+1) - 3*flow_x(i,j) + flow_x(i,j-1)
    !                 flow_x(i,j+3)       = 3*flow_x(i,j+2) - 3*flow_x(i,j+1) + flow_x(i,j)
    !             end if
    !             if (boundary(i,j-1)==2) then
    !                 flow_x(i,j-1)       = 3*flow_x(i,j) - 3*flow_x(i,j+1) + flow_x(i,j+2)
    !                 flow_x(i,j-2)       = 3*flow_x(i,j-1) - 3*flow_x(i,j) + flow_x(i,j+1)
    !                 flow_x(i,j-3)       = 3*flow_x(i,j-2) - 3*flow_x(i,j-1) + flow_x(i,j)
    !             end if
    !             if (boundary(i+1,j)==2) then
    !                 flow_y(i+1,j)       = 3*flow_y(i,j) - 3*flow_y(i-1,j) + flow_y(i-2,j)
    !                 flow_y(i+2,j)       = 3*flow_y(i+1,j) - 3*flow_y(i,j) + flow_y(i-1,j)
    !                 flow_y(i+3,j)       = 3*flow_y(i+2,j) - 3*flow_y(i+1,j) + flow_y(i,j)
    !             end if
    !             if (boundary(i-1,j)==2) then
    !                 flow_y(i-1,j)       = 3*flow_y(i,j) - 3*flow_y(i+1,j) + flow_y(i+2,j)
    !                 flow_y(i-2,j)       = 3*flow_y(i-1,j) - 3*flow_y(i,j) + flow_y(i+1,j)
    !                 flow_y(i-3,j)       = 3*flow_y(i-2,j) - 3*flow_y(i-1,j) + flow_y(i,j)
    !             end if
    !         end if
    !     end do
    ! end do
    ! do i = 1,ny,1
    !     do j = 1,nx,1
    !         if (boundary(i,j)>2) then
    !             if (boundary(i,j+1)==2) then
    !                 flow_x(i,j+1:j+3:1)     = flow_x(i,j)
    !             end if
    !             if (boundary(i,j-1)==2) then
    !                 flow_x(i,j-3:j-1:1)     = flow_x(i,j)
    !             end if
    !             if (boundary(i+1,j)==2) then
    !                 flow_y(i+1:i+3:1,j)     = flow_y(i,j)
    !             end if
    !             if (boundary(i-1,j)==2) then
    !                 flow_y(i-3:i-1:1,j)     = flow_y(i,j)
    !             end if
    !         end if
    !     end do
    ! end do
    ! 1
        ! do i = 1,ny,1
        !     do j = 1,nx,1
        !         if (boundary(i,j)>2) then
        !             if (boundary(i,j+1)==1) then
        !                 flow_x(i,j+1)       = 3*flow_x(i,j) - 3*flow_x(i,j-1) + flow_x(i,j-2)
        !                 flow_x(i,j+2)       = 3*flow_x(i,j+1) - 3*flow_x(i,j) + flow_x(i,j-1)
        !                 flow_x(i,j+3)       = 3*flow_x(i,j+2) - 3*flow_x(i,j+1) + flow_x(i,j)
        !             end if
        !             if (boundary(i,j-1)==1) then
        !                 flow_x(i,j-1)       = 3*flow_x(i,j) - 3*flow_x(i,j+1) + flow_x(i,j+2)
        !                 flow_x(i,j-2)       = 3*flow_x(i,j-1) - 3*flow_x(i,j) + flow_x(i,j+1)
        !                 flow_x(i,j-3)       = 3*flow_x(i,j-2) - 3*flow_x(i,j-1) + flow_x(i,j)
        !             end if
        !             if (boundary(i+1,j)==1) then
        !                 flow_y(i+1,j)       = 3*flow_y(i,j) - 3*flow_y(i-1,j) + flow_y(i-2,j)
        !                 flow_y(i+2,j)       = 3*flow_y(i+1,j) - 3*flow_y(i,j) + flow_y(i-1,j)
        !                 flow_y(i+3,j)       = 3*flow_y(i+2,j) - 3*flow_y(i+1,j) + flow_y(i,j)
        !             end if
        !             if (boundary(i-1,j)==1) then
        !                 flow_y(i-1,j)       = 3*flow_y(i,j) - 3*flow_y(i+1,j) + flow_y(i+2,j)
        !                 flow_y(i-2,j)       = 3*flow_y(i-1,j) - 3*flow_y(i,j) + flow_y(i+1,j)
        !                 flow_y(i-3,j)       = 3*flow_y(i-2,j) - 3*flow_y(i-1,j) + flow_y(i,j)
        !             end if
        !         end if
        !     end do
        ! end do

    ! 0
        do i = 1,ny,1
            do j = 1,nx,1
                if (boundary(i,j)>2) then
                    if (boundary(i,j+1)==0) then
                        flow_x(i,j+1)           =   - flow_x(i,j)
                        flow_x(i,j+2)           =   - flow_x(i,j-1)
                        flow_x(i,j+3)           =   - flow_x(i,j-2)
                    end if
                    if (boundary(i,j-1)==0) then
                        flow_x(i,j-1)           =   - flow_x(i,j)
                        flow_x(i,j-2)           =   - flow_x(i,j+1)
                        flow_x(i,j-3)           =   - flow_x(i,j+2)
                    end if
                    if (boundary(i+1,j)==0) then
                        flow_y(i+1,j)           =   - flow_y(i,j)
                        flow_y(i+2,j)           =   - flow_y(i-1,j)
                        flow_y(i+3,j)           =   - flow_y(i-2,j)
                    end if
                    if (boundary(i-1,j)==0) then
                        flow_y(i-1,j)           =   - flow_y(i,j)
                        flow_y(i-2,j)           =   - flow_y(i+1,j)
                        flow_y(i-3,j)           =   - flow_y(i+2,j)
                    end if
                end if
            end do
        end do
    
    return
end subroutine WENO_ghost

subroutine WENO_Res_5(flow, dim, dir, ny, nx)
    implicit none
    integer, parameter          ::  acc = 8
    real(kind=acc), parameter   ::  epsilon = 1d-6
    integer, Intent(In)         ::  ny, nx, dim, dir
    real(kind=acc)              ::  flow(ny,nx)
    real(kind=acc)              ::  flow_shift(5,ny,nx), up(3, ny,nx), b(3, ny,nx),g(3),&
                                    &w(3,ny,nx), sum_w(ny,nx)
    ! Set Variables
    g(:) = (/1d0/10, 3d0/5, 3d0/10/)

    if ((dim == 2) .AND. (dir == 1)) then
        flow_shift(1,:,:) = eoshift(flow,2,flow(:,1),2)
        flow_shift(2,:,:) = eoshift(flow,1,flow(:,1),2)
        flow_shift(3,:,:) = flow
        flow_shift(4,:,:) = eoshift(flow,-1,flow(:,nx),2)
        flow_shift(5,:,:) = eoshift(flow,-2,flow(:,nx),2)
    else if((dim == 2) .AND. (dir == -1)) then
        flow_shift(1,:,:) = eoshift(flow,-2,flow(:,nx),2)
        flow_shift(2,:,:) = eoshift(flow,-1,flow(:,nx),2)
        flow_shift(3,:,:) = flow
        flow_shift(4,:,:) = eoshift(flow,1,flow(:,1),2)
        flow_shift(5,:,:) = eoshift(flow,2,flow(:,1),2)
    else if((dim == 1) .AND. (dir == -1)) then
        flow_shift(1,:,:) = eoshift(flow,-2,flow(ny,:),1)
        flow_shift(2,:,:) = eoshift(flow,-1,flow(ny,:),1)
        flow_shift(3,:,:) = flow
        flow_shift(4,:,:) = eoshift(flow,1,flow(1,:),1)
        flow_shift(5,:,:) = eoshift(flow,2,flow(1,:),1)
    else if((dim == 1) .AND. (dir == 1)) then
        flow_shift(1,:,:) = eoshift(flow,2,flow(1,:),1)
        flow_shift(2,:,:) = eoshift(flow,1,flow(1,:),1)
        flow_shift(3,:,:) = flow
        flow_shift(4,:,:) = eoshift(flow,-1,flow(ny,:),1)
        flow_shift(5,:,:) = eoshift(flow,-2,flow(ny,:),1)
    end if
                                    
    ! Reconstruction Polynomials
    up(1,:,:) =  1d0/3*flow_shift(1,:,:) - 7d0/6*flow_shift(2,:,:) + 11d0/6*flow_shift(3,:,:)
    up(2,:,:) = -1d0/6*flow_shift(2,:,:) + 5d0/6*flow_shift(3,:,:) + 1d0/3*flow_shift(4,:,:)
    up(3,:,:) =  1d0/3*flow_shift(3,:,:) + 5d0/6*flow_shift(4,:,:) - 1d0/6*flow_shift(5,:,:)
    ! Smooth parameters
    b(1,:,:) = 13d0/12 * (flow_shift(1,:,:)-2d0*flow_shift(2,:,:)+flow_shift(3,:,:))**2 +&
                1d0/4 * (flow_shift(1,:,:)-4d0*flow_shift(2,:,:)+3d0*flow_shift(3,:,:))**2
    b(2,:,:) = 13d0/12 * (flow_shift(2,:,:)-2d0*flow_shift(3,:,:)+flow_shift(4,:,:))**2 +&
                1d0/4* (flow_shift(2,:,:)-flow_shift(4,:,:))**2
    b(3,:,:) = 13d0/12 * (flow_shift(3,:,:)-2d0*flow_shift(4,:,:)+flow_shift(5,:,:))**2 +&
                1d0/4 * (3*flow_shift(3,:,:)-4d0*flow_shift(4,:,:)+flow_shift(5,:,:))**2
    ! weigths
    w(1,:,:) = g(1) / ((epsilon + b(1,:,:))**2)
    w(2,:,:) = g(2) / ((epsilon + b(2,:,:))**2)
    w(3,:,:) = g(3) / ((epsilon + b(3,:,:))**2)
    sum_w = w(1,:,:) + w(2,:,:) + w(3,:,:)
    ! Non-linear weigths
    w(1,:,:) = w(1,:,:) / sum_w
    w(2,:,:) = w(2,:,:) / sum_w
    w(3,:,:) = w(3,:,:) / sum_w
    ! WENO polynomial
    flow = w(1,:,:)*up(1,:,:) + w(2,:,:)*up(2,:,:) + w(3,:,:)*up(3,:,:)
    return
end subroutine WENO_Res_5

subroutine F90_fsmweno3(cost,boundary,potential,p_x,p_y,h,ny,nx,nit,test)
    implicit none
    integer, parameter          ::  acc = 8
    real(kind=acc), parameter   ::  sigma = 1d-11,  p_max = 1d6
    integer, Intent(In)         ::  ny, nx,boundary(ny,nx)
    real(kind=acc), Intent(In)  ::  cost(ny,nx),h
    real(kind=acc), Intent(Out) ::  potential(ny,nx),p_x(ny,nx),p_y(ny,nx),test(ny,nx)
    integer                     ::  i, j, n, nit_c = 1, Gsit(4,2,2), GSi(4,2),nit,bj_fsm(ny,nx),nit_zhang
    real(kind=acc)                ::  p_old(ny,nx), pot_x, pot_y, p_t, p_temp, pot_xs(2), pot_ys(2),&
                                    &sigma_old(200),weight(2,ny,nx,2), p_ste(5)

    ! Gsit    = reshape(  (/1 , nx, nx, 1 ,&
    !                         1, 1, ny , ny ,&
    !                         nx, 1 , 1 , nx,&
    !                         ny , ny , 1, 1/),(/4,2,2/))
    Gsit    = reshape(  (/1 , nx, nx, 1 ,&
                            1, 1, ny , ny ,&
                            nx, 1 , 1 , nx,&
                            ny, ny , 1, 1/),(/4,2,2/))
    GSi     = reshape(  (/ 1, -1, -1,  1,&
                        1, 1,  -1,  -1/),(/4,2/))
    bj_fsm = boundary
    bj_fsm = merge(0,boundary,((boundary==1).OR.(boundary==0)))
    call F90_fsmgod(cost,boundary,potential,p_x,p_y,h,ny,nx,nit)
    p_old = 0d0
    nit = 0
    nit_zhang = 0
    weight = 0d0
    sigma_old = 0d0
    test = 0d0
    do nit = 1,2000,1
        if ( ( norm2(potential - p_old) ) >= sigma ) then
            do n = 1,4,1
                p_old = potential
                do j = GSit(n,2,1), GSit(n,2,2), GSi(n,2)
                    do i = GSit(n,1,1), GSit(n,1,2), GSi(n,1)

                        ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is Dh
                        if ((bj_fsm(j,i)>0).AND.(i>2).AND.(i<(nx-2)).AND.(j>2).AND.(j<(ny-2))) then
                            if ((all(bj_fsm(j,(i-2):(i+2)) > 0)) .AND. (all(bj_fsm((j+2):(j-2):-1,i) >0))) then

                                pot_x = minval(fast_WENO_Res(nit_zhang,weight(1,j,i,:),potential(j,(i-2):(i+2)),h))
                                pot_y = minval(fast_WENO_Res(nit_zhang,weight(2,j,i,:),potential((j+2):(j-2):-1,i),h))
                            
                                if (abs(pot_x-pot_y) > (cost(j,i)*h)) then
                                    potential(j,i) = (min(pot_x, pot_y) + cost(j,i)*h)
                                else 
                                    potential(j,i) = (pot_x+pot_y + sqrt(2d0*(cost(j,i)*h)**2-(pot_x-pot_y)**2)) / 2
                                end if

                            ! else if ((any(boundary(j,(i-2):(i+2)) == 2)) .OR. (any(boundary((j+2):(j-2):-1,i) == 2))) then

                            !     if (any(boundary(j,(i-2):(i+2)) == 2)) then
                            !         pot_x = minval(fast_WENO_Int(nit_zhang,weight(1,j,i,:),potential(j,(i-2):(i+2)),boundary(j,(i-2):(i+2)),h))
                            !     else
                            !         pot_x = minval(fast_WENO_Res(nit_zhang,weight(1,j,i,:),potential(j,(i-2):(i+2)),h))
                            !     end if
                            !     if (any(boundary((j+2):(j-2):-1,i) == 2)) then
                            !         pot_y = minval(fast_WENO_Int(nit_zhang,weight(2,j,i,:),potential((j+2):(j-2):-1,i),boundary((j+2):(j-2):-1,i),h))
                            !     else
                            !         pot_y = minval(fast_WENO_Res(nit_zhang,weight(2,j,i,:),potential((j+2):(j-2):-1,i),h))
                            !     end if
                            !     if (abs(pot_x-pot_y) > (cost(j,i)*h)) then
                            !         potential(j,i) = (min(pot_x, pot_y) + cost(j,i)*h)
                            !     else 
                            !         potential(j,i) = (pot_x+pot_y + sqrt(2d0*(cost(j,i)*h)**2-(pot_x-pot_y)**2)) / 2
                            !     end if

                            end if
                                                     
                        end if

                    end do
                end do

                if ((nit_zhang < 200).AND.(nit>300))  then
                    sigma_old(nit_zhang+1) = norm2(potential - p_old)
                    if ( abs( log10(sigma_old(nit_zhang+1)) - log10(sum(sigma_old(1:(nit_zhang+1)))/(nit_zhang+1)) ) < 1.0 ) then
                        nit_zhang = nit_zhang+ 1
                    else
                        nit_zhang = 0
                    end if
                else
                    nit_zhang = nit_zhang + 1
                end if          

            end do

        else
                do j = 1, ny, 1
                    do i = 1, nx, 1
                        ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is Dh
                            if ((bj_fsm(j,i)>0).AND.(i>2).AND.(i<(nx-2)).AND.(j>2).AND.(j<(ny-2))) then

                                if ((all(bj_fsm(j,(i-2):(i+2)) > 0)) .AND. (all(bj_fsm((j+2):(j-2):-1,i) >0))) then


                                    pot_xs = fast_WENO_Res(nit_zhang,weight(1,j,i,:),potential(j,(i-2):(i+2)),h)
                                    pot_ys = fast_WENO_Res(nit_zhang,weight(2,j,i,:),potential((j+2):(j-2):-1,i),h)
                                
                                    pot_x = minval(pot_xs)
                                    pot_y = minval(pot_ys)
                                    if (abs(pot_x-pot_y) > (cost(j,i)*h)) then
                                        p_temp = (min(pot_x, pot_y) + cost(j,i)*h)
                                    else 
                                        p_temp = (pot_x+pot_y + sqrt(2d0 * (cost(j,i)*h)**2-(pot_x-pot_y)**2)) / 2
                                    end if
                                    
                                    if (pot_xs(1) >= pot_xs(2)) then
                                        p_x(j,i) = (pot_xs(2) - potential(j,i)) / h
                                    else
                                        p_x(j,i) = (potential(j,i) - pot_xs(1)) / h
                                    end if
                                    if (pot_ys(1) >= pot_ys(2)) then
                                        p_y(j,i) = (pot_ys(2) - potential(j,i)) / h
                                    else
                                        p_y(j,i) = (potential(j,i) - pot_ys(1)) / h
                                    end if
                                    potential(j,i) = p_temp

                                ! else if ((any(boundary(j,(i-2):(i+2)) == 2)) .OR. (any(boundary((j+2):(j-2):-1,i) == 2))) then

                                !     if (any(boundary(j,(i-2):(i+2)) == 2)) then
                                !         pot_xs = (fast_WENO_Int(nit_zhang,weight(1,j,i,:),potential(j,(i-2):(i+2)),boundary(j,(i-2):(i+2)),h))
                                !     else
                                !         pot_xs = (fast_WENO_Res(nit_zhang,weight(1,j,i,:),potential(j,(i-2):(i+2)),h))
                                !     end if
                                !     if (any(boundary((j+2):(j-2):-1,i) == 2)) then
                                !         pot_ys = (fast_WENO_Int(nit_zhang,weight(2,j,i,:),potential((j+2):(j-2):-1,i),boundary((j+2):(j-2):-1,i),h))
                                !     else
                                !         pot_ys = (fast_WENO_Res(nit_zhang,weight(2,j,i,:),potential((j+2):(j-2):-1,i),h))
                                !     end if

                                !     pot_x = minval(pot_xs)
                                !     pot_y = minval(pot_ys)
                                !     if (abs(pot_x-pot_y) > (cost(j,i)*h)) then
                                !         p_temp = (min(pot_x, pot_y) + cost(j,i)*h)
                                !     else 
                                !         p_temp = (pot_x+pot_y + sqrt(2d0 * (cost(j,i)*h)**2-(pot_x-pot_y)**2)) / 2
                                !     end if
                                    
                                !     if (pot_xs(1) >= pot_xs(2)) then
                                !         p_x(j,i) = (pot_xs(2) - potential(j,i)) / h
                                !     else
                                !         p_x(j,i) = (potential(j,i) - pot_xs(1)) / h
                                !     end if
                                !     if (pot_ys(1) >= pot_ys(2)) then
                                !         p_y(j,i) = (pot_ys(2) - potential(j,i)) / h
                                !     else
                                !         p_y(j,i) = (potential(j,i) - pot_ys(1)) / h
                                !     end if
                                !     potential(j,i) = p_temp
    
                                end if
                            end if

                    end do
                end do
            
            exit
        end if
    end	do
    test = potential-p_old

    return
contains

! Fast WENO reconstruction
    function fast_WENO_Res(it_non,weight,x_pot,h)
        implicit none
        integer, parameter          ::  acc = 8
        real(kind=acc), parameter   ::  epsilon = 1d-6
        real(kind=acc)              ::  fast_WENO_Res(2), weight(2), x_pot(5), r_m, r_p, w_m, w_p, pot_m, pot_p,h
        integer                     ::  it_non
        if (it_non < 200) then
            r_m = (epsilon + (x_pot(1)-2*x_pot(2)+x_pot(3))**2) / (epsilon + (x_pot(2)-2*x_pot(3)+x_pot(4))**2)
            r_p = (epsilon + (x_pot(5)-2*x_pot(4)+x_pot(3))**2) / (epsilon + (x_pot(4)-2*x_pot(3)+x_pot(2))**2)
            w_m = 1d0 / (1d0 + 2d0 * (r_m**2))
            w_p = 1d0 / (1d0 + 2d0 * (r_p**2))
            weight(1) = (w_m + weight(1)*it_non) / (it_non+1)
            weight(2) = (w_p + weight(2)*it_non) / (it_non+1)
        else
            w_m = weight(1)
            w_p = weight(2)
        end if
        pot_m = (1d0-w_m)*(x_pot(4)-x_pot(2))/(2d0*h) + w_m*( 3d0*x_pot(3)-4d0*x_pot(2)+x_pot(1))/(2d0*h)
        pot_p = (1d0-w_p)*(x_pot(4)-x_pot(2))/(2d0*h) + w_p*(-3d0*x_pot(3)+4d0*x_pot(4)-x_pot(5))/(2d0*h)
        fast_WENO_Res(1) = x_pot(3)-h*pot_m
        fast_WENO_Res(2) = x_pot(3)+h*pot_p
        return
    end function fast_WENO_Res

! Fast WENO reconstruction with interpolation
    function fast_WENO_Int(it_non,weight,x_old,bj_pot,h)
        implicit none
        integer, parameter          ::  acc = 8
        real(kind=acc) , parameter  ::  epsilon = 1d-6
        real(kind=acc)              ::  fast_WENO_Int(2),weight(2),x_pot(5),x_old(5), r_m, r_p, w_m, w_p, pot_m, pot_p,h
        integer                     ::  bj_pot(5),it_non
        x_pot = x_old
        ! interp
            ! if (bj_pot(2)==2) then
            !     x_pot(1:2) = x_pot(3)
            ! else if (bj_pot(1)==2) then
            !     x_pot(1) = x_pot(2)
            ! else if (bj_pot(4)==2) then
            !     x_pot(4:5) = x_pot(3)
            ! else
            !     x_pot(5) = x_pot(4)
            ! end if
            ! if (bj_pot(2)==2) then
            !     x_pot(2) = 2d0*x_pot(3) - x_pot(4)
            !     x_pot(1) = 2d0*x_pot(2) - x_pot(3)
            ! else if (bj_pot(1)==2) then
            !     x_pot(1) = 2d0*x_pot(2) - x_pot(3)
            ! else if (bj_pot(4)==2) then
            !     x_pot(4) = 2d0*x_pot(3) - x_pot(2)
            !     x_pot(5) = 2d0*x_pot(4) - x_pot(3)
            ! else
            !     x_pot(5) = 2d0*x_pot(4) - x_pot(3)
            ! end if
            if (bj_pot(2)==2) then
                x_pot(2) = 3*x_pot(3) - 3*x_pot(4) + x_pot(5)
                x_pot(1) = 3*x_pot(2) - 3*x_pot(3) + x_pot(4)
            else if (bj_pot(1)==2) then
                x_pot(1) = 3*x_pot(2) - 3*x_pot(3) + x_pot(4)
            else if (bj_pot(4)==2) then
                x_pot(4) = 3*x_pot(3) - 3*x_pot(2) + x_pot(1)
                x_pot(5) = 3*x_pot(4) - 3*x_pot(3) + x_pot(2)
            else
                x_pot(5) = 3*x_pot(4) - 3*x_pot(3) + x_pot(2)
            end if
        
        if (it_non < 200) then
            r_m = (epsilon + (x_pot(1)-2d0*x_pot(2)+x_pot(3))**2)/(epsilon + (x_pot(2)-2d0*x_pot(3)+x_pot(4))**2)
            r_p = (epsilon + (x_pot(5)-2d0*x_pot(4)+x_pot(3))**2)/(epsilon + (x_pot(4)-2d0*x_pot(3)+x_pot(2))**2)
            w_m = 1d0 / (1d0 + 2d0 * (r_m**2))
            w_p = 1d0 / (1d0 + 2d0 * (r_p**2))
            weight(1) = (w_m + weight(1)*it_non) / (it_non+1)
            weight(2) = (w_p + weight(2)*it_non) / (it_non+1)
        else
            w_m = weight(1)
            w_p = weight(2)
        end if

        pot_m = (1d0-w_m)*(x_pot(4)-x_pot(2))/(2d0*h) + w_m*( 3d0*x_pot(3)-4d0*x_pot(2)+x_pot(1))/(2d0*h)
        pot_p = (1d0-w_p)*(x_pot(4)-x_pot(2))/(2d0*h) + w_p*(-3d0*x_pot(3)+4d0*x_pot(4)-x_pot(5))/(2d0*h)
        fast_WENO_Int(1) = x_pot(3)-h*pot_m
        fast_WENO_Int(2) = x_pot(3)+h*pot_p
        return
    end function fast_WENO_Int

end subroutine F90_fsmweno3

subroutine F90_fsmgod(cost,boundary,potential,p_x,p_y,h,ny,nx,nit)
    implicit none
    integer, parameter          ::  acc = 8
    real(kind=acc), parameter   ::  sigma = 1d-12,  p_max = 1d6
    integer, Intent(In)         ::  ny, nx,boundary(ny,nx)
    real(kind=acc), Intent(In)  ::  cost(ny,nx),h
    real(kind=acc), Intent(Out) ::  potential(ny,nx),p_x(ny,nx),p_y(ny,nx)
    integer                     ::  i, j, n, nit_c = 1, Gsit(4,2,2), GSi(4,2),nit, bj_fsm(ny,nx)
    real(kind=acc)              ::  p_old(ny,nx), pot_x, pot_y, p_t, p_temp

    Gsit    = reshape(  (/1 , nx, nx, 1 ,&
                            ny, ny, 1 , 1 ,&
                            nx, 1 , 1 , nx,&
                            1 , 1 , ny, ny/),(/4,2,2/))
    GSi     = reshape(  (/ 1, -1, -1,  1,&
                        -1, -1,  1,  1/),(/4,2/))
    p_x = 0d0
    p_y = 0d0
    bj_fsm = boundary
    potential = p_max
    potential = merge(0d0,potential,boundary==1)
    bj_fsm = merge(0,boundary,((boundary==0).OR.(boundary==1)))
    p_old = 0d0
    nit = 0
    ! Obtain potential by Fast Godunov
    do nit = 1,10000,1
        if ( ( sum(abs(potential - p_old)) ) >= sigma ) then
            p_old = potential
            do n = 1,4,1
                do j = GSit(n,2,1), GSit(n,2,2), GSi(n,2)
                    do i = GSit(n,1,1), GSit(n,1,2), GSi(n,1)

                        ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
                        if (bj_fsm(j,i) > 0) then
                                
                            pot_x = min(potential(j,max(1,i-1)), potential(j,min(nx,i+1)))
                            pot_y = min(potential(max(1,j-1),i), potential(min(ny,j+1),i))
                            if (abs(pot_x-pot_y) > (cost(j,i)*h)) then
                                p_temp = (min(pot_x, pot_y) + cost(j,i)*h)
                            else 
                                p_temp = (pot_x+pot_y + sqrt(2d0*(cost(j,i)*h)**2-(pot_x-pot_y)**2)) / 2d0
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
                    end do                  
                end do
            end do
        else
            exit
        end if
        
    end	do

    return
end subroutine F90_fsmgod