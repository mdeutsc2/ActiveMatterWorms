program lbm
    !importing concrete type library
    use, intrinsic :: iso_fortran_env

    implicit none

    logical,parameter :: save_vtk = .true.
    integer(int64),parameter :: save_interval = 500
    
    ! PARAMS
    integer(int64),parameter :: nx = 801 !domain size
    integer(int64),parameter :: ny = 201
    integer(int64),parameter :: nsteps = 20000
    real(real64),parameter :: rho0 = 1.0d0
    real(real64),parameter :: niu = 0.01!viscosity
    integer(int64),parameter :: cy = 1! whether or not to have a cylindrical obstacle
    real(real64),parameter,dimension(3) :: cy_para = (/160.d0,100.d0,20.d0/)

    ! ARRAYS
    real(real64),dimension(0:nx,0:ny) :: rho
    real(real64) :: vel(0:nx,0:ny,2)
    integer(int64),dimension(0:nx,0:ny) :: mask
    real(real64),dimension(9,0:nx,0:ny) :: f_old
    real(real64),dimension(9,0:nx,0:ny) :: f_new
    real(real64),dimension(9) :: w
    real(real64),dimension(2,9) :: e
    integer(int64),dimension(0:3) :: bc_type
    real(real64),dimension(2,0:3) :: bc_value

    ! VARIABLES
    real(real64) ::tau,inv_tau,elapsed,ftime,stime
    integer(int64) :: istep,i
    integer(int32) :: vtkfileunit


    
    w = (/ 4.d0/9.d0, 1.d0/9.d0,1.d0/9.d0,1.d0/9.d0,1.d0/9.d0,1.d0/36.d0,1.d0/36.d0,1.d0/36.d0,1.d0/36.d0/)
    e(1,:) = (/0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0/) ! x
    e(2,:) = (/0.d0, 0.d0, 1.d0, 0.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0/) ! y
    ! ! LID-driven cavity
    ! ![left,top,right,bottom] =0, Dirichlet, =1, Neumann
    ! bc_type = (/0,0,0,0/)
    ! !if bc_type = 0 we need to specify a velocity in bc_value
    ! bc_value(1,:) = (/0.d0, 0.1d0, 0.d0, 0.d0/)
    ! bc_value(2,:) = (/0.d0, 0.d0, 0.d0, 0.d0/)

    ! Karman Vortex Street
    ![left,top,right,bottom] =0, Dirichlet, =1, Neumann
    bc_type = (/0,0,1,0/)
    !if bc_type = 0 we need to specify a velocity in bc_value
    bc_value(1,:) = (/0.1d0, 0.0d0, 0.d0, 0.d0/)
    bc_value(2,:) = (/0.d0, 0.d0, 0.d0, 0.d0/)

    ! INITIALIZATION
    tau = 3.0*niu + 0.5
    inv_tau = 1.0/tau
    print*,0.1*cy_para(3)/niu
    call lbm_init(nx,ny,vel,rho,f_new,f_old,e,w,mask,cy,cy_para,rho0)
    istep = 0
    call write_out(vel,vtkfileunit,istep,nx,ny)
    print *, "istep, max(rho), max(vx), max(vy)"
    do istep = 1,nsteps
        call cpu_time(stime)
        call collide_and_stream(nx,ny,f_new,f_old,e,inv_tau)

        call update_macro_vars(nx,ny,f_new,f_old,vel,rho,e)

        call apply_bc(nx,ny,vel,bc_type,bc_value,rho,e,w,mask,cy,cy_para,f_old)
        call cpu_time(ftime)
        elapsed = elapsed + (ftime-stime)
        !print *, istep,maxval(rho),maxval(vel(:,:,1)),maxval(vel(:,:,2))
        if (save_vtk .eqv. .true.) then
            if (mod(istep,save_interval/5) .eq. 0) then
                print *, istep,maxval(rho),maxval(vel(:,:,1)),maxval(vel(:,:,2)),elapsed/istep
            endif
            if (mod(istep,save_interval) .eq. 0) then
                call write_out(vel,vtkfileunit,istep,nx,ny)
            endif
        endif
    enddo
    
    contains
        subroutine write_out(vel,vtkfileunit,istep,nx,ny)
            real(real64), intent(in out) :: vel(:,:,:)
            integer(int32), intent(in) :: vtkfileunit
            integer(int64), intent(in) :: istep,nx,ny
            integer(int64) :: i,j
            character(len=12):: filename
            
            write(filename, '(A3,I5.5,A4)')'vel',istep,'.vtk'
            print *, filename
            open(unit=vtkfileunit,file=filename,status='unknown')
            ! writing the header
            write(vtkfileunit,'(A)')trim('# vtk DataFile Version 2.0')
            write(vtkfileunit,'(A)')trim(filename)
            write(vtkfileunit,'(A)')'ASCII'
            write(vtkfileunit,'(A)')'DATASET STRUCTURED_POINTS'
            write(vtkfileunit,'(A,1X,I3,1X,I3,1X,I1)')'DIMENSIONS',nx,ny,1
            write(vtkfileunit,'(A)')'ORIGIN 0 0 0'
            write(vtkfileunit,'(A)')'SPACING 1 1 1'
            write(vtkfileunit,'(A,1X,I6)')'POINT_DATA ',nx*ny
            write(vtkfileunit,'(A)')'VECTORS velocity_field float'
            do j = 1,ny
                do i = 1,nx
                    write(vtkfileunit,'(*(es13.6,1x))')vel(i,j,1),vel(i,j,2),0.0d0
                enddo
            enddo
        end subroutine write_out
        
        subroutine lbm_init(nx,ny,vel,rho,f_new,f_old,e,w,mask,cy,cy_para,rho0)
            integer(int64),intent(in) :: nx,ny
            real(real64),intent(in out) :: vel(0:nx,0:ny,2),rho(0:nx,0:ny),f_new(9,0:nx,0:ny),f_old(9,0:nx,0:ny),e(:,:),w(:)
            integer(int64),intent(in out) :: mask(0:nx,0:ny)
            integer(int64),intent(in) :: cy
            real(real64),intent(in) :: cy_para(:),rho0
            integer(int64) :: i,j,k
            real(real64) :: r
            print*, nx,ny
            do i = 0,nx
                do j = 0,ny
                    !print*, i,j
                    vel(i,j,1) = 0.d0
                    vel(i,j,2) = 0.d0
                    rho(i,j) = rho0
                    mask(i,j) = 0
                    do k = 1,9
                        f_new(k,i,j) = f_eq(nx,ny,i,j,k,e,vel,w,rho)
                        f_old(k,i,j) = f_new(k,i,j)
                    enddo
                    r =  (real(i,real64) - cy_para(1))**2.0 + (real(j,real64) - cy_para(2))**2.0 
                    if ((cy .eq. 1) .and. (r .le. cy_para(3))) then
                        mask(i,j) = 1
                    endif
                enddo
            enddo
            return
        end subroutine lbm_init
        
        pure subroutine collide_and_stream(nx,ny,f_new,f_old,e,inv_tau)
            integer(int64),intent(in) :: nx,ny
            real(real64),intent(in out) :: f_new(9,0:nx,0:ny),f_old(9,0:nx,0:ny),e(:,:)
            real(real64),intent(in) :: inv_tau
            integer(int64)::i,j,k,ip,jp
            
            do k = 1,9
                do j = 1,ny-1
                    do i = 1,nx-1
                        ip = i - int(e(1,k),int64)
                        jp = j - int(e(2,k),int64)
                        f_new(k,i,j) = (1.0-inv_tau)*f_old(k,ip,jp) + f_eq(nx,ny,ip,jp,k,e,vel,w,rho)*inv_tau
                    enddo
                enddo
            enddo
            return
        end subroutine collide_and_stream

        pure subroutine update_macro_vars(nx,ny,f_new,f_old,vel,rho,e)
            integer(int64),intent(in) :: nx,ny
            real(real64),intent(in out) :: f_new(9,0:nx,0:ny),f_old(9,0:nx,0:ny),vel(0:nx,0:ny,2),rho(0:nx,0:ny)
            real(real64),intent(in) :: e(:,:)
            integer(int64)::i,j,k

            do j = 1,ny-1
                do i = 1,nx-1
                    rho(i,j) = 0.0d0
                    vel(i,j,1) = 0.0d0
                    vel(i,j,2) = 0.0d0
                    do k = 1,9
                        f_old(k,i,j) = f_new(k,i,j)
                        rho(i,j) = rho(i,j) + f_new(k,i,j)
                        vel(i,j,1) = vel(i,j,1) + e(1,k)*f_new(k,i,j)
                        vel(i,j,2) = vel(i,j,2) + e(2,k)*f_new(k,i,j)
                    enddo
                    vel(i,j,1) = vel(i,j,1) / rho(i,j)
                    vel(i,j,2) = vel(i,j,2) / rho(i,j)
                enddo
            enddo
            return
        end subroutine update_macro_vars

        pure subroutine apply_bc_core(nx,ny,outer,dr,ibc,jbc,inb,jnb,vel,bc_type,bc_value,rho,e,w,f_old)
            integer(int64),intent(in) :: nx,ny
            integer(int64),intent(in) :: outer,dr,ibc,jbc,inb,jnb,bc_type(0:3)
            real(real64),intent(in out) :: vel(0:nx,0:ny,2),bc_value(2,0:3),rho(0:nx,0:ny),e(:,:),w(:),f_old(9,0:nx,0:ny)
            integer(int64) :: k
            if (outer .eq. 1) then !handle outer boundary
                ! bc_type = 0 (Dirichlet) velocity is set as constant value given in bc_value
                ! e.g. if bc_value at boundary is zero, then there is a no-slip condition.
                if (bc_type(dr) .eq. 0) then
                    vel(ibc,jbc,1) = bc_value(1,dr)
                    vel(ibc,jbc,2) = bc_value(2,dr)
                !bc_type = 1 (Neumann), the derivative of velocity in boundary normal direction is set to zero
                ! this is known as bounce-back
                else if (bc_type(dr) .eq. 1) then
                    vel(ibc,jbc,1) = vel(inb,jnb,1)
                    vel(ibc,jbc,2) = vel(inb,jnb,2)
                endif
                rho(ibc,jbc) = rho(inb,jnb)
                do k = 1,9
                    f_old(k,ibc,jbc) = f_eq(nx,ny,ibc,jbc,k,e,vel,w,rho) - f_eq(nx,ny,inb,jnb,k,e,vel,w,rho) + f_old(k,inb,jnb)
                enddo
            endif
            return
        end subroutine apply_bc_core

        pure subroutine apply_bc(nx,ny,vel,bc_type,bc_value,rho,e,w,mask,cy,cy_para,f_old)
            integer(int64),intent(in) :: nx,ny
            real(real64),intent(in out) :: vel(0:nx,0:ny,2),rho(0:nx,0:ny),w(:),bc_value(2,0:3),e(:,:),f_old(9,0:nx,0:ny)
            integer(int64),intent(in out) :: bc_type(0:3),mask(0:nx,0:ny)
            integer(int64),intent(in) :: cy
            real(real64),intent(in) :: cy_para(:)
            integer(int64) :: outer,dr,ibc,jbc,inb,jnb
            integer(int64) :: i,j
            !dr= [left,top,right,bottom] bc_type=0, Dirichlet, bc_type=1, Neumann
            !left and right bc
            do j = 1,ny-1
                ! left: dr = 0, ibc=0, jbc = j, inb = 1, jnb = j
                outer = 1
                dr = 0
                ibc = 0
                jbc = j
                inb = 1
                jnb = j
                call apply_bc_core(nx,ny,outer,dr,ibc,jbc,inb,jnb,vel,bc_type,bc_value,rho,e,w,f_old)
                            
                !right: dr = 2, ibc = nx-1, jbc = j, inb = nx-2, jnb = j
                outer = 1
                dr = 2
                ibc = nx-1
                jbc = j
                inb = nx-2
                jnb = j
                call apply_bc_core(nx,ny,outer,dr,ibc,jbc,inb,jnb,vel,bc_type,bc_value,rho,e,w,f_old)
            enddo
            !top and bottom bc
            do i = 0,nx
                !top: dr = 1, ibc=i, jbc=ny-1, inb = i, jnb = ny-2
                outer = 1
                dr = 1
                ibc = i
                jbc = ny-1
                inb = i
                jnb = ny-2
                call apply_bc_core(nx,ny,outer,dr,ibc,jbc,inb,jnb,vel,bc_type,bc_value,rho,e,w,f_old)
                            
                !bottom: dr = 3, ibc = 1, jbc = 0, inb = i, jnb = 1
                outer = 1
                dr = 3
                ibc = 1
                jbc = 0
                inb = i
                jnb = 1
                call apply_bc_core(nx,ny,outer,dr,ibc,jbc,inb,jnb,vel,bc_type,bc_value,rho,e,w,f_old)
            enddo
                            
            !cylindrical obstacle
            do j = 0,ny
                do i = 0,nx
                    if ((cy .eq. 1) .and. (mask(i,j) .eq. 1)) then
                        vel(i,j,1) = 0.d0 !velocity is zero at solid boundary
                        vel(i,j,2) = 0.d0
                        inb = 0
                        jnb = 0
                        if (real(i,real64) .gt. cy_para(1))then
                            inb = i+1
                        else
                            inb = i-1
                        endif
                        if (real(j,real64) .gt. cy_para(2)) then
                            jnb = j+1
                        else
                            jnb = j-1
                        endif
                        outer = 0
                        dr = 0
                        call apply_bc_core(nx,ny,outer,dr,i,j,inb,jnb,vel,bc_type,bc_value,rho,e,w,f_old)
                    endif
                enddo
            enddo
            return
        end subroutine apply_bc

        pure function f_eq(nx,ny,i,j,k,e,vel,w,rho) result(f)
            integer(int64),intent(in) :: nx,ny
            real(real64),intent(in) :: e(:,:),vel(0:nx,0:ny,2),w(:),rho(0:nx,0:ny)
            integer(int64),intent(in) :: i,j,k
            real(real64) :: f,eu,uv
            eu = e(1,k)*vel(i,j,1) + e(2,k)*vel(i,j,2)
            uv = vel(i,j,1)*vel(i,j,1) + vel(i,j,2)*vel(i,j,2)
            f = w(k)*rho(i,j)*(1.0 + 3.0*eu + 4.5*eu*eu - 1.5*uv)
        end function f_eq
end program lbm
