module global
    implicit none
    integer, parameter :: nx = 100, ny = 50, nB = 3, acc = 8, fast_WENO_switch = 1, nH = 1, &
                        method_de = 1, method_pe = 3, nx_e = nx+2*nB, ny_e = ny+2*nB
    integer :: it_non = 1, iteration = 0
    real(kind=acc), parameter :: zero = 0.0, one = 1.0, xmax = 100d0, ymax = 50d0,&
            h = 1d0, CFL = 2d-1, epsilon = 1.0d-6, sigma = 1.0d-9, potmax = 1.0d12
    real(kind=acc) :: t,cpu_t0,cpu_t1,x(nx_e),y(ny_e),flow_in,potential(ny_e,nx_e),weight(50,8,ny_e,nx_e,2)
    real(kind=acc) :: boundary_O(1,2,2) = reshape([-nB*h, 0.5*h, nH*h, ymax-nH*h], [1,2,2]), &
            boundary_D_2h(1,2,2) = reshape([xmax-2*h, xmax+nB*h, -nB*h, ymax+nB*h], [1,2,2]), &
            boundary_D(2,2,2) = reshape([xmax-h,   xmax-h, &
                                         xmax,     xmax, &
                                         5d0,      30d0, &
                                         20d0,     50d0], [2,2,2]), &
            boundary_H(6,2,2) = reshape([-nB*h,     -nB*h,     xmax-nH*h, xmax-nH*h, xmax-nH*h, 40d0, &
                                         xmax+nB*h, xmax+nB*h, xmax+nB*h, xmax+nB*h, xmax+nB*h, 60d0, &
                                         -nB*h,     ymax-nH*h, zero,      20d0,      45d0,      10d0, &
                                         nH*h,      ymax+nB*h, 5d0,       30d0,      50d0,      30d0],[6,2,2])
            !boundary_D(1,2,2) = reshape([xmax-h, xmax+nB*h, nH*h, ymax-nH*h], [1,2,2]), &
            !boundary_H(2,2,2) = reshape([-nB*h, -nB*h, xmax+nB*h, xmax+nB*h, -nB*h, ymax-nH*h, nH*h, ymax+nB*h], [2,2,2])
end module global

program main
    use global
    implicit none
    integer :: i
    real(kind=acc) :: dt, t_end, div_d(ny,nx), density(ny,nx), density_e(ny_e,nx_e)
    real(kind=acc), allocatable :: density_1(:,:), density_2(:,:), density_3(:,:)
    character*30 :: savefile
    dt = CFL * h / 2.0
    t_end = 200d0
    ! Compute vector arrays
    do i = 1,nx_e
        x(i) = one/2*h + h*(i-1-nB)
    end do
    do i = 1,ny_e
        y(i) = ymax - one/2*h - h*(i-1-nB)
    end do
    ! IC for density and potential
    density = zero;
    potential = 1000d0;
    call boundary_value(x,y,potential,boundary_D,zero);
    call boundary_value(x,y,potential,boundary_H,potmax);
    ! WENO Scheme
    print *, "Start Computing ... ..."
    call cpu_time(cpu_t0)
    do i = 0, int(t_end/dt), 1
        ! Third order TVD
        allocate(density_1(ny,nx), density_2(ny,nx), density_3(ny,nx))
        density_1 = density

        call boundary_expand(density_1,1,nB,method_de)
        call boundary_expand(density_1,2,nB,method_de)
        t = (i) * dt
        call div_density(density_1, div_d)
        density_2 = density + dt * div_d

        call boundary_expand(density_2,1,nB,method_de)
        call boundary_expand(density_2,2,nB,method_de)
        t = (i + 1) * dt
        call div_density(density_2, div_d)
        density_3 = 0.75 * density + 0.25 * ( density_2(1+nB:ny+nB,1+nB:nx+nB) + dt * div_d)

        call boundary_expand(density_3,1,nB,method_de)
        call boundary_expand(density_3,2,nB,method_de)
        !call boundary_value(x,y,density_3,boundary_D,zero);
        t = (i + 0.5) * dt
        call div_density(density_3, div_d)
        density = 1.0/3 * density + 2.0/3 * ( density_3(1+nB:ny+nB,1+nB:nx+nB) + dt * div_d)
        deallocate(density_1,density_2,density_3)
        ! Save density
        if (mod(i*dt,1.0) == 0) then
            if (fast_WENO_switch==0) then
                savefile = "result_g/density"
            else
                savefile = "result_f/density"
            end if
            write(savefile(15:18),"(i4.4)") int(i*dt)
            savefile(19:22) = ".txt"
            call output(density, savefile)
        end if
    end do

contains
    ! Fifth-Order WENO reconstruction
    subroutine div_density(density, div)
        use global
        implicit none
        integer i, j, n
        real(kind=acc) :: div(:,:), density(:,:), cost(ny_e,nx_e), velosity(ny_e,nx_e),flow_mag(ny_e,nx_e),&
                flow_x(ny_e,nx_e), flow_y(ny_e,nx_e), xnew(nx+1), ynew(ny+1), &
                flow_xcom(ny_e,nx+1), flow_ycom(ny+1,nx_e), boundary_judge(ny_e,nx_e) = zero, pot_x(2), pot_y(2),&
                flow_xm(ny_e,nx_e), flow_xp(ny_e,nx_e), flow_ym(ny_e,nx_e), flow_yp(ny_e,nx_e)
        real(kind=acc), allocatable :: pot_e(:,:)
        
        ! Flow in: flow_in
        if (t<=60) then
            flow_in = t / 12.0
        else 
            if (t<=120) then
                flow_in = 10 - t / 12.0
            else 
                flow_in = 0.0
            end if
        end if
        ! Flow location: x_new, y_new
        do i = 1,nx+1
            xnew(i) = h*(i-1)
        end do
        do i = 1,ny+1
            ynew(i) = ymax - h*(i-1)
        end do
        ! Obtain Velosity & Cost & Flow_mag & Potential
        do i = 1,size(density,2),1
            do j = 1,size(density,1),1
                velosity(j,i) = min(2d0,2*(1.0 - max(density(j,i),zero)/10.0))
                cost(j,i) = 1.0 / velosity(j,i) + 0.002 * (max(zero,density(j,i)) )**2
                flow_mag(j,i) = velosity(j,i) * density(j,i)
            end do
        end do
        call fast_WENO(cost)
        ! Split flow_mag
        allocate(pot_e(ny_e,nx_e))
        pot_e = potential
        call boundary_expand(pot_e,1,2,method_pe)
        call boundary_expand(pot_e,2,2,method_pe)
        call boundary_value(x,y,boundary_judge,boundary_H,one)
        if (t==t_end) then
            call output(pot_e,"debug/debug00.txt")
            call output(boundary_judge,"debug/debug0.txt")
        end if
        it_non = 1
        do i = 1,nx_e
            do j = 1,ny_e
                if (boundary_judge(j,i)==zero) then
                    pot_x = fast_WENO_Res(1,j,i,pot_e(j+2,i:i+4))
                    pot_y = fast_WENO_Res(1,j,i,pot_e(j+4:j:-1,i+2))
                    pot_x(1) = ( pot_e(j+2,i+2) - pot_x(1)) / h
                    pot_x(2) = ( pot_x(2) - pot_e(j+2,i+2)) / h
                    pot_y(1) = ( pot_e(j+2,i+2) - pot_y(1)) / h
                    pot_y(2) = ( pot_y(2) - pot_e(j+2,i+2)) / h
                    if (abs(pot_x(1) - pot_x(2)) > cost(j,i)) then
                        if (abs(pot_x(1)) >= abs(pot_x(2))) then
                            flow_x(j,i) = - pot_x(2) * flow_mag(j,i) / cost(j,i)
                        else
                            flow_x(j,i) = - pot_x(1) * flow_mag(j,i) / cost(j,i)
                        end if
                    else
                        flow_x(j,i) = - sum(pot_x) / 2 * flow_mag(j,i) / cost(j,i)
                    end if
                    if (abs(pot_y(1) - pot_y(2)) > cost(j,i)) then
                        if (abs(pot_y(1)) >= abs(pot_y(2))) then
                            flow_y(j,i) = - pot_y(2) * flow_mag(j,i) / cost(j,i)
                        else
                            flow_y(j,i) = - pot_y(1) * flow_mag(j,i) / cost(j,i)
                        end if
                    else
                        flow_y(j,i) = - sum(pot_y) / 2  * flow_mag(j,i) / cost(j,i)
                    end if
                else 
                    flow_x(j,i) = 0
                    flow_y(j,i) = 0
                end if
            end do
        end do
        call boundary_value(x, y, flow_x, boundary_O, flow_in)
        call boundary_value(x, y, flow_y, boundary_O, zero)
        call boundary_value(x, y, flow_x, boundary_H, zero)
        call boundary_value(x, y, flow_y, boundary_H, zero)
        if (t==t_end) then
            call output(flow_x,"debug/debug01.txt")
            call output(flow_y,"debug/debug02.txt")
        end if
        ! Flow x/y upwind/downwind
        flow_xp = 0.5*flow_x + 0.5*spread(maxval(abs(velosity),2),2,nx_e)*density
        flow_xm = 0.5*flow_x - 0.5*spread(maxval(abs(velosity),2),2,nx_e)*density
        flow_yp = 0.5*flow_y + 0.5*spread(maxval(abs(velosity),1),1,ny_e)*density
        flow_ym = 0.5*flow_y - 0.5*spread(maxval(abs(velosity),1),1,ny_e)*density
        if (t==t_end) then
            call output(flow_xp,"debug/debug21.txt")
            call output(flow_xm,"debug/debug22.txt")
            call output(flow_yp,"debug/debug23.txt")
            call output(flow_ym,"debug/debug24.txt")
        end if
        call WENO_Res(flow_xp, 2, 1)
        call WENO_Res(flow_xm, 2,-1)
        call WENO_Res(flow_yp, 1, 1)
        call WENO_Res(flow_ym, 1,-1)
        if (t==t_end) then
            call output(flow_xp,"debug/debug31.txt")
            call output(flow_xm,"debug/debug32.txt")
            call output(flow_yp,"debug/debug33.txt")
            call output(flow_ym,"debug/debug34.txt")
        end if
        ! Combination
        flow_xcom = flow_xp(:,nB:nx+nB) + flow_xm(:,(nB+1):nx+(nB+1))
        flow_ycom = flow_yp((nB+1):ny+(nB+1),:) + flow_ym(nB:ny+nB,:)
        ! call boundary_value(xnew, y, flow_xcom, boundary_O, flow_in)
        flow_xcom(:,nx+1) = flow_xcom(:,nx)
        call boundary_value(x, ynew, flow_ycom, boundary_D, zero)
        call boundary_value(xnew, y, flow_xcom, boundary_H, zero)
        call boundary_value(x, ynew, flow_ycom, boundary_H, zero)
        if (t==t_end) then
            call output(flow_xcom,"debug/debug41.txt")
            call output(flow_ycom,"debug/debug42.txt")
        end if
        ! Divergence
        div = -1./h * (flow_xcom((nB+1):ny+nB,2:nx+1) - flow_xcom((nB+1):ny+nB,1:nx) +&
                        flow_ycom(1:ny,(nB+1):nx+nB) - flow_ycom(2:ny+1,(nB+1):nx+nB))
        !div = max(div, -density/dt)        
        if (t==t_end) then
            call output(div,"debug/debug5.txt")
            call output(div*dt+density,"debug/debug6.txt")
        end if
        deallocate(pot_e)
        return
    end subroutine div_density

    ! WENO reconstruction: 1-y, 2-x, --up wind, +-down wind
    subroutine WENO_Res(flow, dim, dir)
        use global
        integer :: dim, dir
        real(kind=acc) :: flow(:,:)
        real(kind=acc) :: flow_shift(5,size(flow,1),size(flow,2)), up(3, size(flow,1),size(flow,2)),&
            b(3, size(flow,1),size(flow,2)), g(3) = [1./10, 3./5, 3./10],&
            w(3,size(flow,1),size(flow,2)), sum_w(size(flow,1),size(flow,2))
        ! Set Variables
        flow_shift(1,:,:) = eoshift(flow,2*dir,zero,dim)
        flow_shift(2,:,:) = eoshift(flow,dir,zero,dim)
        flow_shift(3,:,:) = flow
        flow_shift(4,:,:) = eoshift(flow,-dir,zero,dim)
        flow_shift(5,:,:) = eoshift(flow,-2*dir,zero,dim)
        ! Reconstruction Polynomials
        up(1,:,:) =  1./3*flow_shift(1,:,:) - 7./6*flow_shift(2,:,:) + 11./6*flow_shift(3,:,:)
        up(2,:,:) = -1./6*flow_shift(2,:,:) + 5./6*flow_shift(3,:,:) + 1./3*flow_shift(4,:,:)
        up(3,:,:) =  1./3*flow_shift(3,:,:) + 5./6*flow_shift(4,:,:) - 1./6*flow_shift(5,:,:)
        ! Smooth parameters
        b(1,:,:) = 13./12 * (flow_shift(1,:,:)-2*flow_shift(2,:,:)+flow_shift(3,:,:))**2 +&
             1./4 * (flow_shift(1,:,:)-4*flow_shift(2,:,:)+3*flow_shift(3,:,:))**2
        b(2,:,:) = 13./12 * (flow_shift(2,:,:)-2*flow_shift(3,:,:)+flow_shift(4,:,:))**2 +&
             1./4 * (flow_shift(2,:,:)-flow_shift(4,:,:))**2
        b(3,:,:) = 13./12 * (flow_shift(3,:,:)-2*flow_shift(4,:,:)+flow_shift(5,:,:))**2 +&
             1./4 * (3*flow_shift(3,:,:)-4*flow_shift(4,:,:)+flow_shift(5,:,:))**2
        ! weigths
        w(1,:,:) = g(1) / ((epsilon+b(1,:,:))**2)
        w(2,:,:) = g(2) / ((epsilon+b(2,:,:))**2)
        w(3,:,:) = g(3) / ((epsilon+b(3,:,:))**2)
        sum_w = w(1,:,:) + w(2,:,:) + w(3,:,:)
        ! Non-linear weigths
        w(1,:,:) = w(1,:,:) / sum_w
        w(2,:,:) = w(2,:,:) / sum_w
        w(3,:,:) = w(3,:,:) / sum_w
        ! WENO polynomial
        flow = w(1,:,:)*up(1,:,:) + w(2,:,:)*up(2,:,:) + w(3,:,:)*up(3,:,:)
        return
    end subroutine WENO_Res

    subroutine WENO_Res_3(flow, dim, dir)
        use global
        integer :: dim, dir
        real(kind=acc) :: flow(:,:)
        real(kind=acc) :: flow_shift(3,size(flow,1),size(flow,2)), up(2, size(flow,1),size(flow,2)),&
            b(2, size(flow,1),size(flow,2)), g(2) = [1./3, 2./3],&
            w(2,size(flow,1),size(flow,2)), sum_w(size(flow,1),size(flow,2))
        ! Set Variables
        flow_shift(1,:,:) = eoshift(flow,dir,zero,dim)
        flow_shift(2,:,:) = flow
        flow_shift(3,:,:) = eoshift(flow,-dir,zero,dim)
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

    ! Third-Order Fast Sweeping WENO
    subroutine fast_WENO(cost)
        use global
        implicit none
        real(kind=acc) :: sigma_old(50), cost(:,:), potential_old(ny_e,nx_e), pot_x, pot_y,&
                        boundary_judge(ny_e,nx_e) = zero, pot_tem
        real(kind=acc), allocatable :: potential_tem(:,:)
        integer :: GSit(4,2,2), GSi(4,2), i, j, n
        data Gsit /1 , nx_e, nx_e, 1 ,&
                   ny_e, ny_e, 1 , 1 ,&
                   nx_e, 1 , 1 , nx_e,&
                   1 , 1 , ny_e, ny_e/,&
             GSi  / 1, -1, -1,  1,&
                   -1, -1,  1,  1/
        call fast_Godunov(cost,1,potential)
        ! Obtain potential by Fast Sweeping WENO
        if (fast_WENO_switch==1) then
        call boundary_value(x,y,boundary_judge,boundary_D_2h,one)
        call boundary_value(x,y,boundary_judge,boundary_H,one)
        potential_old = zero
        iteration = 0
        it_non = 1
        weight = zero
        do while ( ( sum(abs(potential - potential_old)) ) >= sigma )
            iteration = iteration + 1
            potential_old = potential
            do n = 1,4,1
                allocate(potential_tem(ny_e,nx_e))
                potential_tem = potential
                call boundary_expand(potential_tem,1,2,method_pe)
                call boundary_expand(potential_tem,2,2,method_pe)
                do j = GSit(n,2,1), GSit(n,2,2), GSi(n,2)
                    do i = GSit(n,1,1), GSit(n,1,2), GSi(n,1)
                        if (boundary_judge(j,i)==zero) then
                            pot_x=minval(fast_WENO_Res(2*n-1,j,i,potential_tem(j+2,(i):(i+4))))
                            pot_y=minval(fast_WENO_Res(2*n  ,j,i,potential_tem((j+4):(j):-1,i+2)))
                            if (abs(pot_x-pot_y)>=(cost(j,i)*h)) then
                                pot_tem = (min(pot_x, pot_y) + cost(j,i)*h)
                            else 
                                pot_tem = (pot_x+pot_y + (2*(cost(j,i)**2)*(h**2)-(pot_x-pot_y)**2)**0.5) / 2.0
                            end if
                            ! potential_tem(j+2,i+2) = pot_tem
                            potential_tem(j+2,i+2) = min(pot_tem, potential_tem(j+2,i+2))
                        end if
                    end do
                end do
                potential = potential_tem(3:ny_e+2,3:nx_e+2)
                deallocate(potential_tem)

                if ((iteration == 10000).and.(n == 1)) then
                    call output(potential, "debug/debug_potential_10000_1.txt")
                end if
                if ((iteration == 10000).and.(n == 2)) then
                    call output(potential, "debug/debug_potential_10000_2.txt")
                end if
                if ((iteration == 10000).and.(n == 3)) then
                    call output(potential, "debug/debug_potential_10000_3.txt")
                end if
                if ((iteration == 10000).and.(n == 4)) then
                    call output(potential, "debug/debug_potential_10000_4.txt")
                end if

            end do

            ! Convergence method introduced by Zhang
            if (it_non < 50) then
                sigma_old(it_non) = log10(sum(abs(potential - potential_old)))
                if ( abs( log10(sum(abs(potential - potential_old))) -sum(sigma_old(1:it_non))/it_non ) < one ) then
                    it_non = it_non + 1
                else
                    it_non = 1
                end if
            else
                it_non = it_non + 1
            end if
            if (iteration>=1000) then
                if (mod(iteration,1000)==0) then
                    print *, "Debug Iterations:", iteration, it_non, sum(abs(potential - potential_old))
                end if
                if (iteration == 10000) then
                    call output(potential, "debug/debug_potential_10000.txt")
                    call output(cost, "debug/debug_cost.txt")
                end if
                if (iteration == 10001) then
                    call output(potential, "debug/debug_potential_10001.txt")
                end if
                if (iteration == 10002) then
                    call output(potential, "debug/debug_potential_10002.txt")
                    print *, "Fast Sweeping WENO didn't converge! ", sum(abs(potential - potential_old))
                    exit
                end if
            end if

        end do

        end if
        call cpu_time(cpu_t1)
        print *,"Simulation time:", t, "Iterations:", iteration, it_non, "Computing time:", cpu_t1 - cpu_t0
        return
    end subroutine fast_WENO

    ! Fast WENO reconstruction
    function fast_WENO_Res(n,j,i,x_pot)
        use global
        real(kind=acc) :: fast_WENO_Res(2), x_pot(:), r_m, r_p, w_m, w_p, pot_m, pot_p
        integer :: n, j, i
        !if (it_non <= 50) then
            r_m = (epsilon + (x_pot(1)-2*x_pot(2)+x_pot(3))**2)/(epsilon + (x_pot(2)-2*x_pot(3)+x_pot(4))**2)
            r_p = (epsilon + (x_pot(5)-2*x_pot(4)+x_pot(3))**2)/(epsilon + (x_pot(4)-2*x_pot(3)+x_pot(2))**2)
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
        fast_WENO_Res(1) = x_pot(3)-h*pot_m
        fast_WENO_Res(2) = x_pot(3)+h*pot_p
    end function fast_WENO_Res

    ! First-Order Godunov
    subroutine fast_Godunov(cost,h_G,potential_G)
        use global
        implicit none
        integer :: h_G, GSit(4,2,2), GSit_G(4,2,2), GSi(4,2), i, j, n
        real(kind=acc), allocatable :: potential_tem(:,:)
        real(kind=acc) :: potential_G(:,:), cost(:,:), pot_x, pot_y, x_G(h_G*nx_e), y_G(h_G*ny_e), &
                        boundary_judge_G(h_G*ny_e,h_G*nx_e), potential_old_G(h_G*ny_e,h_G*nx_e)
        data Gsit /1 , nx_e, nx_e, 1 ,&
                   ny_e, ny_e, 1 , 1 ,&
                   nx_e, 1 , 1 , nx_e,&
                   1 , 1 , ny_e, ny_e/,&
             GSi  / 1, -1, -1,  1,&
                   -1, -1,  1,  1/
        ! Obtain potential IC by Fast Godunov
        do i = 1,h_G*nx_e
            x_G(i) = 1.0/2*h/h_G + h/h_G*(i-1-nB)
        end do
        do i = 1,h_G*ny_e
            y_G(i) = 50 - 1.0/2*h/h_G - h/h_G*(i-1-nB)
        end do
        GSit_G = Gsit
        do n = 1,4,1
            do i = 1,2,1
                do j = 1,2,1
                    if (GSit_G(n,i,j)/=1) then
                        GSit_G(n,i,j) = GSit_G(n,i,j)*h_G
                    end if
                end do
            end do
        end do 

        boundary_judge_G = zero
        call boundary_value(x_G,y_G,boundary_judge_G,boundary_D,one)
        call boundary_value(x_G,y_G,boundary_judge_G,boundary_H,one)
        potential_G = 1000d0
        call boundary_value(x_G,y_G,potential_G,boundary_D,zero)
        call boundary_value(x_G,y_G,potential_G,boundary_H,potmax)
        potential_old_G = zero
        iteration = 0
        do while ( ( sum(abs(potential_G - potential_old_G)) ) >= sigma )
            iteration = iteration + 1
            allocate(potential_tem(ny_e, nx_e))
            potential_tem = potential_G
            potential_old_G = potential_G
            call boundary_expand(potential_tem,1,1,5)
            call boundary_expand(potential_tem,2,1,5)
            do n = 1,4,1
                do j = GSit_G(n,2,1), GSit_G(n,2,2), GSi(n,2)
                    do i = GSit_G(n,1,1), GSit_G(n,1,2), GSi(n,1)
                        if (boundary_judge_G(j,i)/=one) then
                            pot_x = min(potential_tem(j+1,i), potential_tem(j+1,i+1), potential_tem(j+1,i+2))
                            pot_y = min(potential_tem(j,i+1), potential_tem(j+1,i+1), potential_tem(j+2,i+1))
                            if (abs(pot_x-pot_y)>=(cost(ceiling(1.0*j/h_G),ceiling(1.0*i/h_G))*h/h_G)) then
                                potential_tem(j+1,i+1) = (min(pot_x, pot_y) + cost(ceiling(1.0*j/h_G),ceiling(1.0*i/h_G))*h/h_G)
                            else 
                                potential_tem(j+1,i+1) = (pot_x+pot_y + (2*(cost(ceiling(1.0*j/h_G),ceiling(1.0*i/h_G))*h/h_G)&
                                **2-(pot_x-pot_y)**2)**0.5) / 2.0
                            end if
                        end if 
                    end do
                end do
            end do
            potential_G = potential_tem(2:(h_G*ny_e+1), 2:(h_G*nx_e+1))
            deallocate(potential_tem)
            if (iteration >= 1000) then
                print *, "Fast Godunov didn't converge! "
                exit
            end if
        end do
        return
    end subroutine fast_Godunov

    ! Set boundary values
    subroutine boundary_value(x,y,m,boundary,value)
        implicit none
        real(kind=acc) :: x(:), y(:), m(:,:), boundary(:,:,:), value
        integer :: n, i, j
        do n = 1,size(boundary,1),1
            do i = 1,size(x),1
                do j = size(y),1,-1
                    if ( (x(i)>=boundary(n,1,1)).and.(x(i)<=boundary(n,2,1)).and.&
                    (y(j)>=boundary(n,1,2)).and.(y(j)<=boundary(n,2,2)) ) then
                        m(j,i) = value
                    end if
                end do
            end do
        end do
        return
    end subroutine boundary_value

    ! Expand boundary
    subroutine boundary_expand(m,dim,n,method)
        implicit none
        real(kind=acc),allocatable :: m(:,:)
        real(kind=acc) ::tem(size(m,1), size(m,2))
        integer :: dim, n, method, i, lenx, leny
        lenx = size(m,2)
        leny = size(m,1)
        tem = m
        deallocate(m)
        select case (dim)
            case (1)
                allocate( m(leny, lenx+2*n) )
                m(:,(n+1):(lenx+n)) = tem
                do i = 1,n,1
                    select case (method)
                        case (0)
                            m(:,n+1-i) = 0d0
                            m(:,lenx+n+i) = 0d0
                        case (1)
                            m(:,n+1-i) = m(:,n+2-i)
                            m(:,lenx+n+i) = m(:,lenx+n+i-1)
                        case (2)
                            m(:,n+1-i) = max(0d0, 2*m(:,n+2-i) - m(:,n+3-i))
                            m(:,lenx+n+i) = max(0d0, 2*m(:,lenx+n+i-1) - m(:,lenx+n+i-2))
                        case (3)
                            m(:,n+1-i) = max(0d0, 3*m(:,n+2-i) - 3*m(:,n+3-i) + m(:,n+4-i))
                            m(:,lenx+n+i) = max(0d0, 3*m(:,lenx+n+i-1) - 3*m(:,lenx+n+i-2) + m(:,lenx+n+i-3))
                        case (4)
                            m(:,n+1-i) = max(0d0, 4*m(:,n+2-i) - 6*m(:,n+3-i) + 4*m(:,n+4-i) - m(:,n+5-i))
                            m(:,lenx+n+i) = max(0d0, 4*m(:,lenx+n+i-1) - 6*m(:,lenx+n+i-2) + 4*m(:,lenx+n+i-3) - m(:,lenx+n+i-4))
                        case (5)
                            m(:,n+1-i) = potmax
                            m(:,lenx+n+i) = potmax
                    end select
                end do
            case (2)
                allocate( m(leny+2*n, lenx) )
                m((n+1):(leny+n),:) = tem
                do i = 1,n,1
                    select case (method)
                        case (0)
                            m(n+1-i,:) = 0d0
                            m(leny+n+i,:) = 0d0
                        case (1)
                            m(n+1-i,:) = m(n+2-i,:)
                            m(leny+n+i,:) = m(leny+n+i-1,:)
                        case (2)
                            m(n+1-i,:) = max(0d0, 2*m(n+2-i,:) - m(n+3-i,:))
                            m(leny+n+i,:) = max(0d0, 2*m(leny+n+i-1,:) - m(leny+n+i-2,:))
                        case (3)
                            m(n+1-i,:) = max(0d0, 3*m(n+2-i,:) - 3*m(n+3-i,:) + m(n+4-i,:))
                            m(leny+n+i,:) = max(0d0, 3*m(leny+n+i-1,:) - 3*m(leny+n+i-2,:) + m(leny+n+i-3,:))
                        case (4)
                            m(n+1-i,:) = max(0d0, 4*m(n+2-i,:) - 6*m(n+3-i,:) + 4*m(n+4-i,:) - m(n+5-i,:))
                            m(leny+n+i,:) = max(0d0, 4*m(leny+n+i-1,:) - 6*m(leny+n+i-2,:) + 4*m(leny+n+i-3,:) - m(leny+n+i-4,:))
                        case (5)
                            m(n+1-i,:) = potmax
                            m(leny+n+i,:) = potmax
                    end select
                end do
        end select
        return
    end subroutine boundary_expand

    ! Output
    subroutine output(m,str)
        implicit none
        real(kind=acc) :: m(:,:)
        character(len = *) :: str
        open(1, file = str, status = 'replace')
        write (1, *), m
        close(1)
        return
    end subroutine output
    
end program main


 ! if ((boundary(j,i)==3)) then
                            !     if ((i>2).AND.(i<(nx-1))) then
                            !         p_ste = potential(j,(i-2):(i+2))
                            !     else if (i==1) then
                            !         p_ste(3:5)  = potential(j,(i):(i+2))
                            !         p_ste(2)    = p_max !3*p_ste(3)-3*p_ste(4)+p_ste(5)
                            !         p_ste(1)    = p_max !3*p_ste(2)-3*p_ste(3)+p_ste(4)
                            !     else if (i==2) then
                            !         p_ste(2:5)  = potential(j,(i-1):(i+2))
                            !         p_ste(1)    = p_max !3*p_ste(2)-3*p_ste(3)+p_ste(4)
                            !     else if (i==(nx-1)) then
                            !         p_ste(1:4)  = potential(j,(i-2):(i+1))
                            !         p_ste(5)    = p_max !3*p_ste(4)-3*p_ste(3)+p_ste(2)
                            !     else if (i==nx) then
                            !         p_ste(1:3)  = potential(j,(i-2):(i))
                            !         p_ste(4)    = p_max !3*p_ste(3)-3*p_ste(2)+p_ste(1)
                            !         p_ste(5)    = p_max !3*p_ste(4)-3*p_ste(3)+p_ste(2)
                            !     end if
                            !     pot_xs = fast_WENO_Res(0,weight(1,j,i,:),p_ste,h)

                            !     if ((j>2).AND.(j<(ny-1))) then
                            !         p_ste = potential((j+2):(j-2):-1,i)
                            !     else if (j==ny) then
                            !         p_ste(3:5)  = potential((j):(j-2):-1,i)
                            !         p_ste(2)    = p_max !3*p_ste(3)-3*p_ste(4)+p_ste(5)
                            !         p_ste(1)    = p_max !3*p_ste(2)-3*p_ste(3)+p_ste(4)
                            !     else if (j==(ny-1)) then
                            !         p_ste(2:5)  = potential((j+1):(j-2):-1,i)
                            !         p_ste(1)    = p_max !3*p_ste(2)-3*p_ste(3)+p_ste(4)
                            !     else if (j==2) then
                            !         p_ste(1:4)  = potential((j+2):(j-1):-1,i)
                            !         p_ste(5)    = p_max !3*p_ste(4)-3*p_ste(3)+p_ste(2)
                            !     else if (j==1) then
                            !         p_ste(1:3)  = potential((j+2):(j):-1,i)
                            !         p_ste(4)    = p_max !3*p_ste(3)-3*p_ste(2)+p_ste(1)
                            !         p_ste(5)    = p_max !3*p_ste(4)-3*p_ste(3)+p_ste(2)
                            !     end if
                            !     pot_ys = fast_WENO_Res(0,weight(2,j,i,:),p_ste,h)

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
                            ! end if 


                            ! if (any(boundary(j,(i-2):(i+2)) == 2)) then
                            !     pot_x = minval(fast_WENO_Int(nit_zhang,weight(1,j,i,:),potential(j,(i-2):(i+2)),boundary(j,(i-2):(i+2)),h))
                            ! else
                            !     pot_x = minval(fast_WENO_Res(nit_zhang,weight(1,j,i,:),potential(j,(i-2):(i+2)),h))
                            ! end if
                            ! if (any(boundary((j+2):(j-2):-1,i) == 2)) then
                            !     pot_y = minval(fast_WENO_Int(nit_zhang,weight(2,j,i,:),potential((j+2):(j-2):-1,i),boundary((j+2):(j-2):-1,i),h))
                            ! else
                            !     pot_y = minval(fast_WENO_Res(nit_zhang,weight(2,j,i,:),potential((j+2):(j-2):-1,i),h))
                            ! end if

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