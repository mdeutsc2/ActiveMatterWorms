module fluid_mod
   !importing concrete type library
   use, intrinsic :: iso_fortran_env
   !use active_mod, only: nxcell, nycell
   implicit none

   !DATA DECLARATIONS
   real(real64), dimension(:, :), allocatable, save :: rho
   real(real64), dimension(:, :, :), allocatable, save :: f_old, f_new, vel,force, S
   integer(int64), dimension(:, :), allocatable, save :: mask
   !interpolated velocities of particles from active matter sim
   !real(real64), dimension(:, :), allocatable, save :: intvx, intvy

   real(real64), dimension(9), save :: w
   real(real64), dimension(2, 9), save :: e
   integer(int64), dimension(0:3), save :: bc_type
   real(real64), dimension(2, 0:3), save :: bc_value
   integer(int64), save :: nx, ny, centerx, centery,res
   real(real64), save :: tau, inv_tau
   integer(int32), save :: vtkfileunit
   logical, parameter :: save_vtk = .true.

contains

   subroutine init_fluid(nxin, nyin, rho0, niu,res)
      integer(int64), intent(in) :: nxin, nyin, res
      real(real64), intent(in) :: rho0, niu
      integer(int64) :: stat, memsize, i, j, k, istep, dr
      real(real64) :: r
      character(100) :: errmsg

      w = (/4.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0/)
      e(1, :) = (/0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0/) ! x
      e(2, :) = (/0.d0, 0.d0, 1.d0, 0.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0/) ! y
      ![left,top,right,bottom] =0, Dirichlet, =1, Neumann
      bc_type = (/0, 0, 0, 0/)
      !if bc_type = 0 we need to specify a velocity in bc_value
      bc_value(1, :) = (/0.d0, 0.0d0, 0.d0, 0.d0/)
      bc_value(2, :) = (/0.d0, 0.d0, 0.d0, 0.d0/)

      ! if the cylinder radius is r, have the domain be 2r and set the mask
      !nx = 2*dr+1
      !ny = 2*dr+1
      !centerx = nx/2
      nx = nxin
      ny = nyin
      centerx = (nx - 1*res)/2
      !centery = ny/2
      centery = (ny - 1*res)/2
      dr = centerx - 1*res
      ! now that we have the domain size, allocate the arrays
      allocate (rho(0:nx, 0:ny), stat=stat, errmsg=errmsg)
      allocate (vel(0:nx, 0:ny, 2), stat=stat, errmsg=errmsg)
      allocate (force(0:nx,0:ny,2), stat=stat, errmsg=errmsg) ! fluid-bead force
      allocate (f_new(9, 0:nx, 0:ny), stat=stat, errmsg=errmsg)
      allocate (f_old(9, 0:nx, 0:ny), stat=stat, errmsg=errmsg)
      allocate (S(9, 0:nx, 0:ny), stat=stat, errmsg=errmsg)
      allocate (mask(0:nx, 0:ny), stat=stat, errmsg=errmsg)
      if (stat .gt. 0) then
         error stop errmsg
      end if
      memsize = sizeof(rho) + sizeof(vel) + sizeof(f_new) + sizeof(f_old) + &
                sizeof(f_new) + sizeof(mask) + sizeof(w) + sizeof(e) + &
                sizeof(force) + sizeof(S)
      print *, "fluid allocations done,", float(memsize)/1e6, " mb"

      ! test if centerx or centery is off-by-1 that there is enough space
      !if (((centerx .le. dr) .or. (centery .le. dr)) .or. (((nx-centerx) .le. dr) .or. ((ny-centery) .le. dr))) then
      tau = 3.0*niu + 0.5
      inv_tau = 1.0/tau
      print *, "Total fluid size:", nx, ny
      print *, "boundary radius:", dr, "centered", centerx, centery
      !initializing values for distributions, density and velocity
      do i = 0, nx
         do j = 0, ny
            !print*, i,j
            vel(i, j, 1) = 0.d0
            vel(i, j, 2) = 0.d0
            rho(i, j) = rho0
            !mask(i,j) = 0 ! remove for circular boundary
            do k = 1, 9
               f_new(k, i, j) = f_eq(nx, ny, i, j, k, e, vel, w, rho)
               f_old(k, i, j) = f_new(k, i, j)
               S(k,i,j) = 0.0d0
            end do
            r = (real(i, real64) - real(centerx, real64))**2 + (real(j, real64) - real(centery, real64))**2
            if (r .ge. (dr**2)) then
               mask(i, j) = 1
            end if
         end do
      end do
      !vel(centerx:centerx+40,centery:centery+60,2) = -0.1
      !vel(centerx:centerx+40,centery:centery-60,2) = 0.1
      istep = 0
      call write_out(istep, nx, ny)!, vel)
   end subroutine init_fluid

   subroutine fluid_step(istep, save_interval)
      integer(int64), intent(in) :: istep, save_interval

      !print*,istep,sum(f_new),sum(vel(:,:,1)),sum(vel(:,:,2))
      call collide_and_stream(nx, ny, f_new, f_old, e, inv_tau)

      call update_macro_vars(nx, ny, f_new, f_old, vel, rho, e)

      call apply_bc(nx, ny, vel, bc_type, bc_value, rho, e, w, mask, f_old)
      
      ! if (save_vtk) then
      !    if (mod(istep, save_interval) .eq. 0) then
      !       call write_out(istep, nx, ny, vel)
      !    end if
      ! end if
   end subroutine fluid_step

   subroutine write_out(istep, nx, ny)!, vel)
      integer(int64), intent(in) :: istep, nx, ny
      !real(real64), intent(in) :: vel(0:nx, 0:ny, 2)
      integer(int64) :: i, j,ic
      character(len=14) :: filename

      write (filename, '(A3,I6.6,A4)') 'fvel', istep, '.vtk'
      !print *, filename,nx,ny,sum(vel(:, :, 1)), sum(vel(:, :, 2))
      open (unit=vtkfileunit, file=filename, status='unknown')
      ! writing the header
      write (vtkfileunit, '(A)') trim('# vtk DataFile Version 2.0')
      write (vtkfileunit, '(A)') trim(filename)
      write (vtkfileunit, '(A)') 'ASCII'
      write (vtkfileunit, '(A)') 'DATASET STRUCTURED_POINTS'
      write (vtkfileunit, '(A,1X,I3,1X,I3,1X,I1)') 'DIMENSIONS', nx, ny, 1
      write (vtkfileunit, '(A)') 'ORIGIN 0 0 0'
      write (vtkfileunit, '(A)') 'SPACING 1 1 1'
      write (vtkfileunit, '(A,1X,I6)') 'POINT_DATA ', nx*ny
      write (vtkfileunit, '(A)') 'VECTORS velocity_field float'
      ic = 9
      do j = 0,ny
         do i = 0,nx
            write (vtkfileunit, '(*(es13.6,1x))') vel(i, j, 1), vel(i, j, 2), 0.0d0
            ic = ic + 1
         end do
      end do
      print *, filename," lines written",ic
   end subroutine write_out

   subroutine collide_and_stream(nx, ny, f_new, f_old, e, inv_tau)
      integer(int64), intent(in) :: nx, ny
      real(real64), intent(in out) :: f_new(9, 0:nx, 0:ny), f_old(9, 0:nx, 0:ny), e(:, :)
      real(real64), intent(in) :: inv_tau
      integer(int64)::i, j, k, ip, jp
      !print*,"collide_and_stream start",sum(f_new),sum(vel)
#ifdef _OPENMP
!$omp parallel do
#endif
      do k = 1, 9
         do j = 1, ny - 1
            do i = 1, nx - 1
               ip = i - int(e(1, k), int64)
               jp = j - int(e(2, k), int64)
               f_new(k, i, j) = (1.0 - inv_tau)*f_old(k, ip, jp) + f_eq(nx, ny, ip, jp, k, e, vel, w, rho)*inv_tau + S(k,i,j)
            end do
         end do
      end do
#ifdef _OPENMP
!$omp end parallel do
#endif
      !print*,"collide_and_stream end",sum(f_new),sum(vel)
      return
   end subroutine collide_and_stream

   subroutine update_macro_vars(nx, ny, f_new, f_old, vel, rho, e)
      integer(int64), intent(in) :: nx, ny
      real(real64), intent(in out) :: f_new(9, 0:nx, 0:ny), f_old(9, 0:nx, 0:ny), vel(0:nx, 0:ny, 2), rho(0:nx, 0:ny)
      real(real64), intent(in) :: e(:, :)
      integer(int64)::i, j, k
      !print*,"macro_vars start",sum(f_new),sum(vel)
#ifdef _OPENMP
!$omp parallel do
#endif
      do j = 1, ny - 1
         do i = 1, nx - 1
            rho(i, j) = 0.0d0
            vel(i, j, 1) = 0.0d0
            vel(i, j, 2) = 0.0d0
            do k = 1, 9
               f_old(k, i, j) = f_new(k, i, j)
               rho(i, j) = rho(i, j) + f_new(k, i, j)
               vel(i, j, 1) = vel(i, j, 1) + e(1, k)*f_new(k, i, j)
               vel(i, j, 2) = vel(i, j, 2) + e(2, k)*f_new(k, i, j)
            end do
            vel(i, j, 1) = vel(i, j, 1)/rho(i, j) + force(i,j,1)/(2*rho(i,j))
            vel(i, j, 2) = vel(i, j, 2)/rho(i, j) + force(i,j,2)/(2*rho(i,j))
         end do
      end do
#ifdef _OPENMP
!$omp end parallel do
#endif
      !print*,"macro_vars end",sum(f_new),sum(vel)
      return
   end subroutine update_macro_vars

   pure subroutine apply_bc_core(nx, ny, outer, dr, ibc, jbc, inb, jnb, vel, &
                                 bc_type, bc_value, rho, e, w, f_old)
      integer(int64), intent(in) :: nx, ny
      integer(int64), intent(in) :: outer, dr, ibc, jbc, inb, jnb, bc_type(0:3)
      real(real64), intent(in out) :: vel(0:nx, 0:ny, 2), bc_value(2, 0:3), rho(0:nx, 0:ny)
      real(real64), intent(in out) ::  e(:, :), w(:), f_old(9, 0:nx, 0:ny)
      integer(int64) :: k
      if (outer .eq. 1) then !handle outer boundary
         if (bc_type(dr) .eq. 0) then
            vel(ibc, jbc, 1) = bc_value(1, dr)
            vel(ibc, jbc, 2) = bc_value(2, dr)
         else if (bc_type(dr) .eq. 1) then
            vel(ibc, jbc, 1) = vel(inb, jnb, 1)
            vel(ibc, jbc, 2) = vel(inb, jnb, 2)
         end if
         rho(ibc, jbc) = rho(inb, jnb)
         do k = 1, 9
            f_old(k, ibc, jbc) = f_eq(nx, ny, ibc, jbc, k, e, vel, w, rho) - &
                                 f_eq(nx, ny, inb, jnb, k, e, vel, w, rho) + f_old(k, inb, jnb)
         end do
      end if
      return
   end subroutine apply_bc_core

   subroutine apply_bc(nx, ny, vel, bc_type, bc_value, rho, e, w, mask, f_old)
      integer(int64), intent(in) :: nx, ny
      real(real64), intent(in out) :: vel(0:nx, 0:ny, 2), rho(0:nx, 0:ny), w(:)
      real(real64), intent(in out) :: bc_value(2, 0:3), e(:, :), f_old(9, 0:nx, 0:ny)
      integer(int64), intent(in out) :: bc_type(0:3), mask(0:nx, 0:ny)
      integer(int64) :: outer, dr, ibc, jbc, inb, jnb
      integer(int64) :: i, j, centerx, centery
      !dr= [left,top,right,bottom] bc_type=0, Dirichlet, bc_type=1, Neumann
      ! for the circular active matter boundary, bc_type=0 always
      ! and we modify the cylindrical boundary routine
      centerx = nx/2
      centery = ny/2
      !print*,"bc start",sum(f_new),sum(vel)
      do j = 1, ny - 1
         do i = 1, nx - 1
            if (mask(i, j) .eq. 1) then
               vel(i, j, 1) = 0.d0 !velocity is zero at solid boundary
               vel(i, j, 2) = 0.d0
               inb = 0
               jnb = 0
               if (real(i, real64) .gt. dr + centerx) then
                  inb = i - 1
               else
                  inb = i + 1
               end if
               if (real(j, real64) .gt. dr + centery) then
                  jnb = j - 1
               else
                  jnb = j + 1
               end if
               outer = 0
               dr = 0
               call apply_bc_core(nx, ny, outer, dr, i, j, inb, jnb, vel, bc_type, bc_value, rho, e, w, f_old)
            end if
         end do
      end do
      !print*,"bc start",sum(f_new),sum(vel)
      return
   end subroutine apply_bc

   pure function f_eq(nx, ny, i, j, k, e, vel, w, rho) result(f)
      integer(int64), intent(in) :: nx, ny
      real(real64), intent(in) :: e(:, :), vel(0:nx, 0:ny, 2), w(:), rho(0:nx, 0:ny)
      integer(int64), intent(in) :: i, j, k
      real(real64) :: f, eu, uv
      eu = e(1, k)*vel(i, j, 1) + e(2, k)*vel(i, j, 2)
      uv = vel(i, j, 1)*vel(i, j, 1) + vel(i, j, 2)*vel(i, j, 2)
      f = w(k)*rho(i, j)*(1.0 + 3.0*eu + 4.5*eu*eu - 1.5*uv)
   end function f_eq

end module fluid_mod
