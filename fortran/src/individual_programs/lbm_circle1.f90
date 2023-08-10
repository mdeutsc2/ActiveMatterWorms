program lbm
   !importing concrete type library
   use, intrinsic :: iso_fortran_env

   implicit none

   logical, parameter :: save_vtk = .true.
   integer(int64), parameter :: save_interval = 2000

   ! PARAMS
   integer(int64), parameter :: dr = 200
   integer(int64), parameter :: resolution = 1
   integer(int64), parameter :: nsteps = 20000
   real(real64), parameter :: rho0 = 1.0d0
   real(real64), parameter :: niu = 0.01!viscosity
   integer(int64), parameter :: cy = 0! whether or not to have a cylindrical obstacle
   real(real64), parameter, dimension(3) :: cy_para = (/160.d0, 100.d0, 20.d0/)

   ! ARRAYS
   real(real64), dimension(:, :), allocatable :: rho
   real(real64), dimension(:, :, :), allocatable :: f_old, f_new, vel
   integer(int64), dimension(:, :), allocatable :: mask
   real(real64), dimension(9) :: w
   real(real64), dimension(2, 9) :: e
   integer(int64), dimension(0:3) :: bc_type
   real(real64), dimension(2, 0:3) :: bc_value

   ! VARIABLES
   real(real64) ::tau, inv_tau, elapsed, ftime, stime, r
   integer(int64) :: istep, i, j, k
   integer(int32) :: vtkfileunit
   integer(int64) :: nx, ny, centerx, centery

   w = (/4.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/9.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0/)
   e(1, :) = (/0.d0, 1.d0, 0.d0, -1.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0/) ! x
   e(2, :) = (/0.d0, 0.d0, 1.d0, 0.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0/) ! y
   ! ! LID-driven cavity
   ! ![left,top,right,bottom] =0, Dirichlet, =1, Neumann
   ! bc_type = (/0,0,0,0/)
   ! !if bc_type = 0 we need to specify a velocity in bc_value
   ! bc_value(1,:) = (/0.d0, 0.1d0, 0.d0, 0.d0/)
   ! bc_value(2,:) = (/0.d0, 0.d0, 0.d0, 0.d0/)

   ! Karman Vortex Street
   ![left,top,right,bottom] =0, Dirichlet, =1, Neumann
   bc_type = (/0, 0, 0, 0/)
   !if bc_type = 0 we need to specify a velocity in bc_value
   bc_value(1, :) = (/0.d0, 0.0d0, 0.d0, 0.d0/)
   bc_value(2, :) = (/0.d0, 0.d0, 0.d0, 0.d0/)

   ! INITIALIZATION
   ! if the cylinder radius is r, have the domain be 2r and set the mask
   nx = 2*dr + 1
   ny = 2*dr + 1
   centerx = nx/2
   centery = ny/2
   ! now that we have the domain size, allocate the arrays
   allocate (rho(1:nx + 1, 1:ny + 1))
   allocate (vel(1:nx + 1, 1:ny + 1, 2))
   allocate (f_new(9, 1:nx + 1, 1:ny + 1))
   allocate (f_old(9, 1:nx + 1, 1:ny + 1))
   allocate (mask(1:nx + 1, 1:ny + 1))

   ! test if centerx or centery is off-by-1 that there is enough space
   !if (((centerx .le. dr) .or. (centery .le. dr)) .or. (((nx-centerx) .le. dr) .or. ((ny-centery) .le. dr))) then
   print *, nx, ny, centerx, centery, dr

   tau = 3.0*niu + 0.5
   inv_tau = 1.0/tau
   call lbm_init(nx, ny, vel, rho, f_new, f_old, e, w, mask, cy, cy_para, rho0, dr)
   !setting a chunk of vx to 0.1 to test the circular bc
   vel(centerx:centerx + 40, centery:centery + 60, 2) = -0.1
   vel(centerx:centerx + 40, centery:centery - 60, 2) = 0.1
   istep = 0
   call write_out(vel, vtkfileunit, istep, nx, ny)
   !print *, "istep, max(rho), max(vx), max(vy)"
   print *, "istep, sum(rho), sum(vx), sum(vy), iter/s"
   do istep = 1, nsteps
      call cpu_time(stime)
      call collide_and_stream(nx, ny, f_new, f_old, e, inv_tau)

      call update_macro_vars(nx, ny, f_new, f_old, vel, rho, e)

      call apply_bc(nx, ny, vel, bc_type, bc_value, rho, e, w, mask, cy, cy_para, f_old)
      call cpu_time(ftime)
      elapsed = elapsed + (ftime - stime)
      !print *, istep,sum(rho),sum(vel(:,:,1)),sum(vel(:,:,2)),elapsed/istep
      ! if (mod(istep,save_interval/5) .eq. 0) then
      !     print *, istep,maxval(rho),maxval(vel(:,:,1)),maxval(vel(:,:,2)),elapsed/istep
      ! endif
      if (save_vtk .eqv. .true.) then
         if (mod(istep, save_interval) .eq. 0) then
            call write_out(vel, vtkfileunit, istep, nx, ny)
         end if
      end if
   end do

contains
   subroutine write_out(vel, vtkfileunit, istep, nx, ny)
      real(real64), intent(in out) :: vel(:, :, :)
      integer(int32), intent(in) :: vtkfileunit
      integer(int64), intent(in) :: istep, nx, ny
      integer(int64) :: i, j
      character(len=12):: filename

      write (filename, '(A3,I5.5,A4)') 'vel', istep, '.vtk'
      print *, filename, sum(vel(:, :, 1)), sum(vel(:, :, 2))
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
      do j = 1, ny + 1
         do i = 1, nx + 1
            write (vtkfileunit, '(*(es13.6,1x))') vel(i, j, 1), vel(i, j, 2), 0.0d0
         end do
      end do
   end subroutine write_out

   subroutine lbm_init(nx, ny, vel, rho, f_new, f_old, e, w, mask, cy, cy_para, rho0, dr)
      integer(int64), intent(in) :: nx, ny, dr
      real(real64), intent(in out) :: vel(1:nx + 1, 1:ny + 1, 2), rho(1:nx + 1, 1:ny + 1)
      real(real64), intent(in out) :: f_new(9, 1:nx + 1, 1:ny + 1), f_old(9, 1:nx + 1, 1:ny + 1), e(:, :), w(:)
      integer(int64), intent(in out) :: mask(0:nx, 0:ny)
      integer(int64), intent(in) :: cy
      real(real64), intent(in) :: cy_para(:), rho0
      integer(int64) :: i, j, k
      real(real64) :: r
      print *, nx, ny
      do i = 1, nx + 1
         do j = 1, ny + 1
            !print*, i,j
            vel(i, j, 1) = 0.d0
            vel(i, j, 2) = 0.d0
            rho(i, j) = rho0
            mask(i, j) = 0
            do k = 1, 9
               f_new(k, i, j) = f_eq(nx, ny, i, j, k, e, vel, w, rho)
               f_old(k, i, j) = f_new(k, i, j)
            end do
            r = (real(i, real64) - real(centerx, real64))**2 + (real(j, real64) - real(centery, real64))**2
            if (r .ge. dr**2) then
               mask(i, j) = 1
            end if
            r = (real(i, real64) - cy_para(1))**2.0 + (real(j, real64) - cy_para(2))**2.0
            if ((cy .eq. 1) .and. (r .le. cy_para(3))) then
               mask(i, j) = 1
            end if
         end do
      end do
      return
   end subroutine lbm_init

   pure subroutine collide_and_stream(nx, ny, f_new, f_old, e, inv_tau)
      integer(int64), intent(in) :: nx, ny
      real(real64), intent(in out) :: f_new(9, 1:nx + 1, 1:ny + 1), f_old(9, 1:nx + 1, 1:ny + 1), e(:, :)
      real(real64), intent(in) :: inv_tau
      integer(int64)::i, j, k, ip, jp

      do k = 1, 9
         do j = 2, ny
            do i = 2, nx
               ip = i - int(e(1, k), int64)
               jp = j - int(e(2, k), int64)
               f_new(k, i, j) = (1.0 - inv_tau)*f_old(k, ip, jp) + f_eq(nx, ny, ip, jp, k, e, vel, w, rho)*inv_tau
            end do
         end do
      end do
      return
   end subroutine collide_and_stream

   pure subroutine update_macro_vars(nx, ny, f_new, f_old, vel, rho, e)
      integer(int64), intent(in) :: nx, ny
      real(real64), intent(in out) :: vel(1:nx + 1, 1:ny + 1, 2), rho(1:nx + 1, 1:ny + 1)
      real(real64), intent(in out) :: f_new(9, 1:nx + 1, 1:ny + 1), f_old(9, 1:nx + 1, 1:ny + 1)
      real(real64), intent(in) :: e(:, :)
      integer(int64)::i, j, k

      do j = 2, ny
         do i = 2, nx
            rho(i, j) = 0.0d0
            vel(i, j, 1) = 0.0d0
            vel(i, j, 2) = 0.0d0
            do k = 1, 9
               f_old(k, i, j) = f_new(k, i, j)
               rho(i, j) = rho(i, j) + f_new(k, i, j)
               vel(i, j, 1) = vel(i, j, 1) + e(1, k)*f_new(k, i, j)
               vel(i, j, 2) = vel(i, j, 2) + e(2, k)*f_new(k, i, j)
            end do
            vel(i, j, 1) = vel(i, j, 1)/rho(i, j)
            vel(i, j, 2) = vel(i, j, 2)/rho(i, j)
         end do
      end do
      return
   end subroutine update_macro_vars

   pure subroutine apply_bc_core(nx, ny, outer, dr, ibc, jbc, inb, jnb, vel, bc_type, bc_value, rho, e, w, f_old)
      integer(int64), intent(in) :: nx, ny
      integer(int64), intent(in) :: outer, dr, ibc, jbc, inb, jnb, bc_type(0:3)
      real(real64), intent(in out) :: vel(1:nx + 1, 1:ny + 1, 2), bc_value(2, 0:3)
      real(real64), intent(in out) :: rho(1:nx + 1, 1:ny + 1), e(:, :), w(:), f_old(9, 1:nx + 1, 1:ny + 1)
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
     f_old(k, ibc, jbc) = f_eq(nx, ny, ibc, jbc, k, e, vel, w, rho) - f_eq(nx, ny, inb, jnb, k, e, vel, w, rho) + f_old(k, inb, jnb)
         end do
      end if
      return
   end subroutine apply_bc_core

   pure subroutine apply_bc(nx, ny, vel, bc_type, bc_value, rho, e, w, mask, cy, cy_para, f_old)
      integer(int64), intent(in) :: nx, ny
      real(real64), intent(in out) :: vel(1:nx + 1, 1:ny + 1, 2), rho(1:nx + 1, 1:ny + 1), w(:)
      real(real64), intent(in out) :: bc_value(2, 0:3), e(:, :), f_old(9, 1:nx + 1, 1:ny + 1)
      integer(int64), intent(in out) :: bc_type(0:3), mask(0:nx, 0:ny)
      integer(int64), intent(in) :: cy
      real(real64), intent(in) :: cy_para(:)
      integer(int64) :: outer, dr, ibc, jbc, inb, jnb
      integer(int64) :: i, j, centerx, centery
      !dr= [left,top,right,bottom] bc_type=0, Dirichlet, bc_type=1, Neumann
      ! for the circular active matter boundary, bc_type=0 always
      ! and we modify the cylindrical boundary routine
      centerx = nx/2
      centery = ny/2
      do j = 2, ny
         do i = 2, nx
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

      return
   end subroutine apply_bc

   pure function f_eq(nx, ny, i, j, k, e, vel, w, rho) result(f)
      integer(int64), intent(in) :: nx, ny
      real(real64), intent(in) :: e(:, :), vel(1:nx + 1, 1:ny + 1, 2), w(:), rho(1:nx + 1, 1:ny + 1)
      integer(int64), intent(in) :: i, j, k
      real(real64) :: f, eu, uv
      eu = e(1, k)*vel(i, j, 1) + e(2, k)*vel(i, j, 2)
      uv = vel(i, j, 1)*vel(i, j, 1) + vel(i, j, 2)*vel(i, j, 2)
      f = w(k)*rho(i, j)*(1.0 + 3.0*eu + 4.5*eu*eu - 1.5*uv)
   end function f_eq
end program lbm
