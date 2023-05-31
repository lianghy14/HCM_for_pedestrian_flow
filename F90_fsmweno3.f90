#include "fintrf.h" !必须有的头文件,里面有mxGetPr, mxGetM, mxGetN,mxCreateDoubleMatrix等函数的申明

Subroutine mexFunction(OutSum,OutVar,InSum,InVar)!函数接口名称必须为mexFunction,
    !OutSum:输出参数个数, OutVar:输出参数数组指针, InSum:输入参数个数, InVar:输入参数数组指针

    Integer InSum,OutSum
    mwPointer InVar(*),OutVar(*)                           !mwPointer专门用于表示指针变量,这个不能随意用Integer代替
    mwPointer mxGetPr, mxGetM, mxGetN, mxCreateDoubleMatrix !这个对返回指针函数的再次申明,
    integer, parameter              ::  acc = 8
    integer, Allocatable            ::  boundary(:,:)
    Real(kind = acc),Allocatable    ::  cost(:,:),potential(:,:),p_x(:,:),p_y(:,:)
    integer                         ::  nit,ny,nx
    real(kind = acc)                ::  h

    If(InSum/=3)Then
        call mexErrMsgIdAndTxt('MATLAB:InputTooBig','Less than 3 inputs.')
        Return
    EndIf

    ny = mxGetM(InVar(1))!获取第1个输入参数的行数
    nx = mxGetN(InVar(1))!获取第1个输入参数的列数
    Allocate(cost(ny,nx),boundary(ny,nx),potential(ny,nx),p_x(ny,nx),p_y(ny,nx))

    Call mxCopyPtrToReal8(mxGetPr(InVar(1)),cost,nx*ny)
    Call mxCopyPtrToReal8(mxGetPr(InVar(2)),boundary,nx*ny)
    Call mxCopyPtrToReal8(mxGetPr(InVar(3)),h,1)

    Call F90_fsmweno3(cost,boundary,potential,p_x,p_y,h,ny,nx,nit)!调用内部函数

    OutVar(1) = mxCreateDoubleMatrix(ny,nx,0)!给返回参数分配内存
    OutVar(2) = mxCreateDoubleMatrix(ny,nx,0)!给返回参数分配内存
    OutVar(3) = mxCreateDoubleMatrix(ny,nx,0)!给返回参数分配内存
    OutVar(4) = mxCreateDoubleMatrix(1,1,0)!给返回参数分配内存

    Call mxCopyReal8ToPtr(potential,mxGetPr(OutVar(1)),nx*ny)
    Call mxCopyReal8ToPtr(p_x,mxGetPr(OutVar(2)),nx*ny)
    Call mxCopyReal8ToPtr(p_y,mxGetPr(OutVar(3)),nx*ny)
    Call mxCopyReal8ToPtr(nit,mxGetPr(OutVar(4)),1)

    DeAllocate(cost,boundary,potential,p_x,p_y)!释放临时分配的内存
    Return

End SubRoutine

subroutine F90_fsmweno3(cost,boundary,potential,p_x,p_y,h,ny,nx,nit)
    implicit none
    integer, parameter          ::  acc = 8
    real, parameter             ::  sigma = 1.0d-9,  p_max = 1.0d12
    integer, Intent(In)         ::  ny, nx, boundary(ny,nx)
    integer, Intent(Out)        ::  nit
    real(kind=acc), Intent(In)  ::  cost(ny,nx),h
    real(kind=acc), Intent(Out) ::  potential(ny,nx),p_x(ny,nx),p_y(ny,nx)
    integer                     ::  i, j, n, nit_c = 1, Gsit(4,2,2), GSi(4,2)
    real(kind=acc)              ::  p_GOD(ny,nx), p_old(ny,nx), pot_x, pot_y, weight = 0.0, sigma_old(50)

    Gsit    = reshape(  (/1 , nx, nx, 1 ,&
                        ny, ny, 1 , 1 ,&
                        nx, 1 , 1 , nx,&
                        1 , 1 , ny, ny/),(/4,2,2/))
    GSi     = reshape(  (/ 1, -1, -1,  1,&
                        -1, -1,  1,  1/),(/4,2/))
    Call F90_fsmgod(cost,boundary,potential,h,ny,nx)
    p_old = 0.0
    p_x = 0.0
    p_y = 0.0
    nit = 0
    do while ( ( sum(abs(potential - p_old)) ) >= sigma )
        nit = nit + 1
        p_old = potential

        do n = 1,4,1
            do j = GSit(n,2,1), GSit(n,2,2), GSi(n,2)
                do i = GSit(n,1,1), GSit(n,1,2), GSi(n,1)

                    ! 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
                    if (boundary(j,i)>3.5) then

                        pot_x=minval(WENO3_Res(h,potential_tem(j,i-2:(i+2):1)))
                        pot_y=minval(WENO3_Res(h,potential_tem((j+2):(j-2):-1,i)))
                        if (abs(pot_x-pot_y)>=(cost(j,i)*h)) then
                            potential(j,i) = (min(pot_x, pot_y) + cost(j,i)*h)
                        else 
                            potential(j,i) = (pot_x+pot_y + (2*(cost(j,i)**2)*(h**2)-(pot_x-pot_y)**2)**0.5) / 2.0
                        end if

                    end if

                end do
            end do
        end do

        ! Convergence method introduced by Zhang
        if (nit_c < 50) then
            sigma_old(nit_c) = log10(sum(abs(potential - p_old)))
            if ( abs(log10(sum(abs(potential - p_old)))-sum(sigma_old(1:nit_c))/nit_c) < 1.0 ) then
                nit_c = nit_c + 1
            else
                nit_c = 1
            end if
        else
            nit_c = nit_c + 1
        end if

        if (nit >= 1000) then
            potential = abs(potential - p_old)
            exit
        end if

    end do

    return
end subroutine F90_fsmweno3

! Fast WENO3 reconstruction
function WENO3_Res(h,x_pot)
    real(kind=acc)  ::  h, WENO3_Res(2), x_pot(:), r_m, r_p, w_m, w_p, pot_m, pot_p
    integer         ::  n, j, i

    !if (it_non <= 50) then
        r_m = (1.0d-6 + (x_pot(1)-2*x_pot(2)+x_pot(3))**2)/(1.0d-6 + (x_pot(2)-2*x_pot(3)+x_pot(4))**2)
        r_p = (1.0d-6 + (x_pot(5)-2*x_pot(4)+x_pot(3))**2)/(1.0d-6 + (x_pot(4)-2*x_pot(3)+x_pot(2))**2)
        w_m = 1./(1 + 2 * r_m**2)
        w_p = 1./(1 + 2 * r_p**2)
        !weight(it_non,n,j,i,1) = w_m
        !weight(it_non,n,j,i,2) = w_p
    !else
        !w_m = sum(weight(:,n,j,i,1)) / 50
        !w_p = sum(weight(:,n,j,i,2)) / 50
    !end if
    pot_m = (1.-w_m)*(x_pot(4)-x_pot(2))/(2.*h) + w_m*( 3*x_pot(3)-4*x_pot(2)+x_pot(1))/(2.*h)
    pot_p = (1.-w_p)*(x_pot(4)-x_pot(2))/(2.*h) + w_p*(-3*x_pot(3)+4*x_pot(4)-x_pot(5))/(2.*h)
    WENO3_Res(1) = x_pot(3)-h*pot_m
    WENO3_Res(2) = x_pot(3)+h*pot_p

end function WENO3_Res

! First-Order Godunov
subroutine F90_fsmgod(cost,boundary,p_GOD,h,ny,nx)
    implicit none
    integer, parameter          ::  acc = 8
    real, parameter             ::  sigma = 1.0d-9,  p_max = 1.0d12
    integer, Intent(In)         ::  ny, nx, boundary(ny,nx)
    real(kind=acc), Intent(In)  ::  cost(ny,nx),h
    real(kind=acc), Intent(Out) ::  p_GOD(ny,nx)
    integer                     ::  i, j, n, nit = 0, nit_c = 1, Gsit(4,2,2), GSi(4,2)
    real(kind=acc)              ::  p_old(ny,nx), pot_x, pot_y

    Gsit    = reshape(  (/1 , nx, nx, 1 ,&
                            ny, ny, 1 , 1 ,&
                            nx, 1 , 1 , nx,&
                            1 , 1 , ny, ny/),(/4,2,2/))
    GSi     = reshape(  (/ 1, -1, -1,  1,&
                        -1, -1,  1,  1/),(/4,2/))

    p_GOD = p_max * boundary
    p_old = 0.0

    ! Obtain p_GOD by Fast Godunov
    do while ( ( sum(abs(p_GOD - p_old)) ) >= sigma )
        nit = nit + 1
        p_old = p_GOD
        
        do n = 1,4,1
            do j = GSit(n,2,1), GSit(n,2,2), GSi(n,2)
                do i = GSit(n,1,1), GSit(n,1,2), GSi(n,1)

                    ! 0 is bound_H, 1 is bound_D, 2 is bound_O
                    if (boundary(j,i)>2.5) then
                        
                        pot_x = min(p_GOD(j,max(1,i-1)), p_GOD(j,i), p_GOD(j,min(nx,i+1)))
                        pot_y = min(p_GOD(max(1,j-1),i), p_GOD(j,i), p_GOD(min(ny,j+1),i))
                        if (abs(pot_x-pot_y)>=(cost(j,i))) then
                            p_GOD(j,i) = (min(pot_x, pot_y) + cost(j,i)*h)
                        else 
                            p_GOD(j,i) = (pot_x+pot_y + (2*(cost(j,i)*h)**2-(pot_x-pot_y)**2)**0.5) / 2.0
                        end if
                        
                    end if 

                end do
            end do
        end do
        
        if (nit >= 1000) then
            p_GOD = nit
            exit
        end if
    end	do
end subroutine F90_fsmgod