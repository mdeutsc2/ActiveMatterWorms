module active_mod
   !importing concrete type library
   use, intrinsic :: iso_fortran_env
   use fluid_mod
   implicit none

   !DATA DECLARATIONS
   real(real64), dimension(:, :), allocatable, save :: x, y, vx, vy, vxave, vyave, fx, fy, fxold, fyold
   real(real64), dimension(:), allocatable, save :: savex, savey
   real(real64), dimension(:,:), allocatable,save :: intvx,intvy
   integer(int64), dimension(:), allocatable, save :: ireverse, ddx, ddy, hhead, ipointto
   integer(int64), dimension(:, :), allocatable, save :: nnab

   real(real64), save :: r2cut, r2cutsmall, rwall, r2inside, density, hx, hy, hxo2, hyo2
   integer(int64), save :: nxcell, nycell, ncells, iseed, nstepso500, ireadsaved, iwalldrive
   real(real64), save :: dcell, gnoise, dto2, dt2o2, length2, lengthmax, thetanow, a, rmin
   integer(int64), save :: scell, scnab
   real(real64), save :: twopi, pi, pio4
   ! IO variables and format statements
   integer(int32), save :: xyzfileunit, activefileunit
   !character(len=:), allocatable, save :: xyzfmt, datfmt

contains

   subroutine init_worms(nsteps, dt, nworms, np, rcut, rcutsmall, length0, fluid_cpl,res)
      integer(int64), intent(in) :: nworms, np, nsteps,res
      real(real64), intent(in) :: dt, rcut, rcutsmall, length0
      logical, intent(in) :: fluid_cpl
      integer(int64) :: stat, memsize, iw, ip, i
      real(real64) :: rand1, dth, xangle, r
      character(100) :: errmsg

      allocate (x(nworms, np), stat=stat, errmsg=errmsg)
      allocate (y(nworms, np), stat=stat, errmsg=errmsg)
      allocate (vx(nworms, np), stat=stat, errmsg=errmsg)
      allocate (vy(nworms, np), stat=stat, errmsg=errmsg)
      allocate (vxave(nworms, np), stat=stat, errmsg=errmsg)
      allocate (vyave(nworms, np), stat=stat, errmsg=errmsg)
      allocate (fx(nworms, np), stat=stat, errmsg=errmsg)
      allocate (fy(nworms, np), stat=stat, errmsg=errmsg)
      allocate (fxold(nworms, np), stat=stat, errmsg=errmsg)
      allocate (fyold(nworms, np), stat=stat, errmsg=errmsg)

      allocate (savex(np), stat=stat, errmsg=errmsg)
      allocate (savey(np), stat=stat, errmsg=errmsg)
      allocate (ireverse(nworms), stat=stat, errmsg=errmsg)
      allocate (ddx(9), stat=stat, errmsg=errmsg)
      allocate (ddy(9), stat=stat, errmsg=errmsg)
      allocate (hhead(504*504), stat=stat, errmsg=errmsg)
      allocate (ipointto(nworms*np), stat=stat, errmsg=errmsg)
      allocate (nnab(nworms, np), stat=stat, errmsg=errmsg)

      if (stat .gt. 0) then
         error stop errmsg
      end if

      !xyzfmt = '(A2,1x,f8.4,1x,f8.4,1x,f8.4)'
      r2cut = rcut*rcut
      r2cutsmall = rcutsmall**2
      rwall = 125.d0*rcutsmall*dsqrt(2.d0)
      twopi = 8.d0*datan(1.d0)
      pi = 0.5d0*twopi
      pio4 = pi*0.25d0
      density = nworms*np/(pi*rwall**2)
      hx = 2.d0*rwall + 1.d0
      hy = hx
      hyo2 = 0.5d0*hy
      hxo2 = 0.5d0*hx
      nxcell = (hx/rcutsmall) - 1
      nycell = (hy/rcutsmall) - 1
      dcell = hx/float(nxcell)
      ncells = nxcell*nycell

      if(fluid_cpl) then ! allocating interpolated velocities
         allocate (intvx(nxcell*res+1,nycell*res+1),stat=stat, errmsg=errmsg)
         allocate (intvy(nxcell*res+1,nycell*res+1),stat=stat, errmsg=errmsg)
         if (stat .gt. 0) then
            error stop errmsg
         end if
         intvx(:,:) = 0.0d0
         intvy(:,:) = 0.0d0
      end if

      memsize = sizeof(x) + sizeof(y) + sizeof(vx) + sizeof(vy) + &
               sizeof(vxave) + sizeof(vyave) + sizeof(fx) + sizeof(fy) + &
               sizeof(fxold) + sizeof(fyold) + sizeof(savex) + &
               sizeof(savey) + sizeof(ireverse) + sizeof(ddx) + &
               sizeof(ddy) + sizeof(hhead) + sizeof(ipointto) + &
               sizeof(nnab) + sizeof(intvx) + sizeof(intvy)

      print *, "active matter allocations done,", memsize/1e6, " mb"

      iseed = 6487735
      call random_init(repeatable=.true., image_distinct=.false.)
      !ss=rand(iseed)
      !ss=rand()
      nstepso500 = nsteps/500
      gnoise = .8d0/dsqrt(10.d0)*0.8d0
      dt2o2 = dt*dt*0.5d0
      dto2 = dt*0.5d0

      length2 = 2.d0*length0
      lengthmax = (length0)*float(np - 1)

      !array for locating neighbor cells
      ddx(1) = 1
      ddy(1) = 0
      ddx(2) = 1
      ddy(2) = 1
      ddx(3) = 0
      ddy(3) = 1
      ddx(4) = -1
      ddy(4) = 1
      ddx(5) = -1
      ddy(5) = 0
      ddx(6) = -1
      ddy(6) = -1
      ddx(7) = 0
      ddy(7) = -1
      ddx(8) = 1
      ddy(8) = -1
      ddx(9) = 0
      ddy(9) = 0

      r2inside = (rwall - rcutsmall)**2
      thetanow = 5.d0*pi
      a = 0.24
      rmin = a*thetanow
      iwalldrive = 1
      print *, "nworms", nworms
      print *, "np", np
      print *, "rwall", rwall
      print *, "nxcells", nxcell, "nycells", nycell
      print *, "density", density
      !setting up the worms
      do iw = 1, nworms
         ireverse(iw) = 0
         call random_number(rand1)
         if (rand1 .le. 0.5d0) then
            ireverse(iw) = 1
         end if
         do i = 1, np
            r = a*thetanow
            dth = length0/r
            thetanow = thetanow + dth
            x(iw, i) = hxo2 + r*dcos(thetanow)
            y(iw, i) = hyo2 + r*dsin(thetanow)
            xangle = atan2(y(iw, i) - hyo2, x(iw, i) - hxo2)
            !give them an initial velocity going around the circle
            !vx(iw,i)=-dsin(xangle)*0.2d0*float(ireverse(iw))
            !vy(iw,i)=dcos(xangle)*0.2d0*float(ireverse(iw))
            !vx(iw,i)=-dsin(xangle)*0.2d0
            !vy(iw,i)=dcos(xangle)*0.2d0
            vx(iw, i) = 0.d0
            vy(iw, i) = 0.d0
            fx(iw, i) = 0.d0
            fy(iw, i) = 0.d0
         end do
         thetanow = thetanow + 4.d0*dth
      end do

      !reverse some of the worms and give them crazy colors
      do iw = 1, nworms
         call random_number(rand1)
         if (ireverse(iw) .eq. 1) then
            do i = 1, np
               savex(i) = x(iw, i)
               savey(i) = y(iw, i)
            end do
            do ip = 1, np
               x(iw, ip) = savex(np + 1 - ip)
               y(iw, ip) = savey(np + 1 - ip)
            end do
         end if
      end do
      print *, minval(x), maxval(x)
      print *, minval(y), maxval(y)
      print *, "Init done"
   end subroutine init_worms

   subroutine write_xyz(istep,nworms, np)
      integer(int64), intent(in) :: istep,nworms, np
      !integer(int32), intent(in) :: xyzfileunit
      integer(int64) :: iw, i,ic
      real(real64) :: dx, dy, xang, rx, ry, dot
      character(len=18) :: filename
      write (filename, '(A7,I7.7,A4)') 'amatter',istep,".xyz"
      open (unit=xyzfileunit, file=filename, status='unknown')

      write (xyzfileunit, *) nworms*np + 4
      write (xyzfileunit, *) "# 0"
      ic = 2
      do iw = 1, nworms
         dx = x(iw, 1) - hxo2
         dy = y(iw, 1) - hyo2
         xang = atan2(dy, dx)
         rx = -dsin(xang)
         ry = dcos(xang)
         dot = (x(iw, 1) - x(iw, np))*rx + (y(iw, 1) - y(iw, np))*ry
         if (dot .ge. 0.d0) then
            do i = 1, np
               write (xyzfileunit, '(A2,1x,f8.4,1x,f8.4,1x,f8.4)') 'A ', x(iw, i), y(iw, i), 0.0d0
               ic = ic + 1
            end do
         else
            do i = 1, np
               write (xyzfileunit, '(A2,1x,f8.4,1x,f8.4,1x,f8.4)') 'B ', x(iw, i), y(iw, i), 0.0d0
               ic = ic + 1
            end do
         end if
      end do
      write (xyzfileunit, '(A2,1x,f8.4,1x,f8.4,1x,f8.4)') 'E ', hxo2 - rwall, hyo2 - rwall, 0.0d0
      write (xyzfileunit, '(A2,1x,f8.4,1x,f8.4,1x,f8.4)') 'E ', hxo2 - rwall, hyo2 + rwall, 0.0d0
      write (xyzfileunit, '(A2,1x,f8.4,1x,f8.4,1x,f8.4)') 'E ', hxo2 + rwall, hyo2 - rwall, 0.0d0
      write (xyzfileunit, '(A2,1x,f8.4,1x,f8.4,1x,f8.4)') 'E ', hxo2 + rwall, hyo2 + rwall, 0.0d0
      print *, filename," lines written",ic
   end subroutine write_xyz

   subroutine write_active_vel(istep, nxin, nyin)
      integer(int64), intent(in) :: istep, nxin, nyin
      !real(real64), intent(in) :: intvx(nx,ny), intvy(nx,ny)
      integer(int64) :: i, j,ic
      character(len=15) :: filename
      ic = 0
      write (filename, '(A3,I7.7,A4)') 'avel', istep, '.vtk'
      open (unit=activefileunit, file=filename, status='unknown')
      ! writing the header
      write (activefileunit, '(A)') trim('# vtk DataFile Version 2.0')
      write (activefileunit, '(A)') trim(filename)
      write (activefileunit, '(A)') 'ASCII'
      write (activefileunit, '(A)') 'DATASET STRUCTURED_POINTS'
      write (activefileunit, '(A,1X,I3,1X,I3,1X,I1)') 'DIMENSIONS', nxin, nyin, 1
      write (activefileunit, '(A)') 'ORIGIN 0 0 0'
      write (activefileunit, '(A)') 'SPACING 1 1 1'
      write (activefileunit, '(A,1X,I6)') 'POINT_DATA ', nxin*nyin
      write (activefileunit, '(A)') 'VECTORS velocity_field float'
      ic = ic + 9
      do j = 1, nyin
         do i = 1, nxin
            write (activefileunit, '(*(es15.8,1x))') intvx(i, j), intvy(i, j), 0.0d0
            !print*,i*j,intvx(i,j),intvy(i,j)
            ic = ic + 1
         end do
      end do
      print *, filename," lines written",ic
   end subroutine write_active_vel

   subroutine update_pos(nworms, np, dt)
      integer(int64), intent(in) :: nworms, np
      real(real64), intent(in) :: dt
      integer(int64) :: iw, i
      do iw = 1, nworms
         do i = 1, np
            x(iw, i) = x(iw, i) + vx(iw, i)*dt + fx(iw, i)*dt2o2
            y(iw, i) = y(iw, i) + vy(iw, i)*dt + fy(iw, i)*dt2o2
            fxold(iw, i) = fx(iw, i)
            fyold(iw, i) = fy(iw, i)
         end do
      end do
   end subroutine update_pos

   subroutine update_vel(nworms, np)
      integer(int64), intent(in) :: nworms, np
      integer(int64) :: iw, i
      do iw = 1, nworms
         do i = 1, np
            vx(iw, i) = vx(iw, i) + dto2*(fx(iw, i) + fxold(iw, i))
            vy(iw, i) = vy(iw, i) + dto2*(fy(iw, i) + fyold(iw, i))
            vxave(iw, i) = vxave(iw, i)/nnab(iw, i)
            vyave(iw, i) = vyave(iw, i)/nnab(iw, i)
         end do
      end do
   end subroutine update_vel

   subroutine calc_forces(nworms, np, kspring, kbend, length0)
      integer(int64), intent(in) :: nworms, np
      real(real64), intent(in) :: kspring, kbend, length0
      real(real64):: rsq, rand1, rand2, v1, v2, fac, g1, th
      real(real64):: dx, dy, r, ff, ffx, ffy, dot, f2x, f2y, f3x, f3y, f4x, f4y
      real(real64)::x2, y2, x3, y3, x4, y4, y23, y34, x23, x34, r23, r34, cosvalue, sinvalue
      integer(int64) :: iw, i, ip1, i2, i3, i4
      !zero out the force arrays and add Gaussian noise
      rsq = 0.0d0
      do iw = 1, nworms
         do i = 1, np
            do while (rsq .ge. 0.999 .or. rsq .le. 0.001)
               call random_number(rand1)
               call random_number(rand2)
               v1 = 2.0*rand1 - 1.d0
               v2 = 2.0*rand2 - 1.d0
               rsq = v1*v1 + v2*v2
            end do
            call random_number(rand1)
            fac = dsqrt(-2.0*dlog(rsq)/rsq)
            g1 = v1*fac*gnoise
            th = rand1*twopi
            fx(iw, i) = g1*dcos(th)
            fy(iw, i) = g1*dsin(th)
         end do
      end do
      !first set of springs nearest neighbor springs
      do iw = 1, nworms
         do i = 1, np - 1
            ip1 = i + 1

            dx = x(iw, ip1) - x(iw, i)
            dy = y(iw, ip1) - y(iw, i)
            r = dsqrt(dx*dx + dy*dy)

            ff = -kspring*(r - length0)/r
            ffx = ff*dx
            ffy = ff*dy
            fx(iw, ip1) = fx(iw, ip1) + ffx
            fx(iw, i) = fx(iw, i) - ffx
            fy(iw, ip1) = fy(iw, ip1) + ffy
            fy(iw, i) = fy(iw, i) - ffy
         end do
      end do
      !bond bending terms
      do iw = 1, nworms
         do i2 = 1, np - 2
            i3 = i2 + 1
            i4 = i2 + 2
            !print*, i2,i3,i4
            x2 = x(iw, i2)
            y2 = y(iw, i2)
            x3 = x(iw, i3)
            y3 = y(iw, i3)
            x4 = x(iw, i4)
            y4 = y(iw, i4)
            y23 = y3 - y2
            y34 = y4 - y3
            x23 = x3 - x2
            x34 = x4 - x3
            r23 = dsqrt(x23*x23 + y23*y23)
            r34 = dsqrt(x34*x34 + y34*y34)

            cosvalue = abs(x23*x34 + y23*y34)/(r23*r34) 
            !cosvalue = (x23*x34 + y23*y34)/(r23*r34)
            !print*, x23,x34,y23,y34,r23,r34
            if (cosvalue .lt. 0.d0) then
               cosvalue = 0.0d0
               !print *, x23, x34, y23, y34, r23, r34
               !print *, "cosvalue .lt. 0.0d0"
               !stop
            end if
            if (cosvalue .gt. 1.d0) then
               cosvalue = 1.d0
            end if
            sinvalue = dsqrt(1.d0 - cosvalue*cosvalue)

            ff = -kbend*sinvalue/(r23*r34)

            dot = x23*x34 + y23*y34
            fac = dot/(r23*r23)

            f2x = ff*(x34 - fac*x23)
            f2y = ff*(y34 - fac*y23)

            fac = dot/(r34*r34)
            f4x = ff*(fac*x34 - x23)
            f4y = ff*(fac*y34 - y23)
            f3x = -f2x - f4x
            f3y = -f2y - f4y

            fx(iw, i2) = fx(iw, i2) + f2x
            fy(iw, i2) = fy(iw, i2) + f2y

            fx(iw, i3) = fx(iw, i3) + f3x
            fy(iw, i3) = fy(iw, i3) + f3y

            fx(iw, i4) = fx(iw, i4) + f4x
            fy(iw, i4) = fy(iw, i4) + f4y
         end do
      end do
   end subroutine calc_forces

   subroutine worm_wall(nworms, np, fdogic, fdogicwall, diss)
      integer(int64), intent(in) :: nworms, np
      real(real64), intent(in) :: fdogic, fdogicwall, diss
      real(real64) :: dx, dy, r, r2, th, xwall, ywall, rr2, ffor, dxi, dyi, ri, dxj, dyj, ffx, ffy
      integer(int64) :: iw, i, ip1
      !put worm-wall interactions here, and dissipation force proportional to velocity
      do iw = 1, nworms
         do i = 1, np
            !dissipation proportional to v relative to local average
            ! TODO swap vxave with intvx
            fx(iw, i) = fx(iw, i) - diss*(vx(iw, i) - vxave(iw, i))
            fy(iw, i) = fy(iw, i) - diss*(vy(iw, i) - vyave(iw, i))
            !now that we have used them, zero out vxave and vyave, recalculate below
            vxave(iw, i) = vx(iw, i)
            vyave(iw, i) = vy(iw, i)
            nnab(iw, i) = 1
            !calculate distance to the center
            dx = x(iw, i) - hxo2
            dy = y(iw, i) - hyo2
            r2 = (dx*dx + dy*dy)
            !if close enough to the wall, calculate wall forces
            !use the short cut-off
            if (r2 .ge. r2inside) then
               !find the nearest spot on the wall
               th = atan2(dy, dx)
               xwall = hxo2 + rwall*dcos(th)
               ywall = hyo2 + rwall*dsin(th)
               dx = xwall - x(iw, i)
               dy = ywall - y(iw, i)
               rr2 = dx*dx + dy*dy
               r = dsqrt(rr2)
               !ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
               ffor = -48.d0*rr2**(-7) + 24.d0*rr2**(-4)
               fx(iw, i) = fx(iw, i) + ffor*dx
               fy(iw, i) = fy(iw, i) + ffor*dy
               
               iwalldrive = 1
               !Turning on dogic drive with the wall!!!
               if (iwalldrive .eq. 1) then
                  !first calculate unit vector along the worm
                  ip1 = i + 1
                  if (ip1 .le. np) then
                     dxi = x(iw, ip1) - x(iw, i)
                     dyi = y(iw, ip1) - y(iw, i)
                  else
                     dxi = x(iw, i) - x(iw, i - 1)
                     dyi = y(iw, i) - y(iw, i - 1)
                  end if
                  !make it a unit vector
                  ri = dsqrt(dxi*dxi + dyi*dyi)
                  dxi = dxi/ri
                  dyi = dyi/ri
                  !calculate the unit vector along the wall
                  dxj = -dsin(th)
                  dyj = dcos(th)

                  !if the vectors are not antiparallel, reverse the vector along the wall
                  if (dxi*dxj + dyi*dyj .gt. 0.d0) then
                     dxj = -dxj
                     dyj = -dyj
                  end if

                  !if the two vectors have any component pointing in opposite directions
                  if (dxi*dxj + dyi*dyj .lt. 0.d0) then
                     !Find the direction for the force...
                     dx = (dxi - dxj)/2.d0
                     dy = (dyi - dyj)/2.d0
                     !normalize the direction vector
                     ri = dsqrt(dx*dx + dy*dy)
                     dx = dx/ri
                     dy = dy/ri

                     !turn on extra-strong driving force
                     ffx = fdogicwall*dx
                     ffy = fdogicwall*dy
                     fx(iw, i) = fx(iw, i) + ffx
                     fy(iw, i) = fy(iw, i) + ffy
                  end if
               end if
            end if
         end do
      end do
   end subroutine worm_wall

   subroutine cell_sort(nworms, np, fdogic, fdep)
      integer(int64), intent(in) :: nworms, np
      real(real64), intent(in) :: fdogic, fdep
      real(real64) :: dddx, dddy, r2, riijj, ffor, ffx, ffy, dxi, dyi, dxj, dyj, ri, rj, r, dx, dy
      integer(int64) :: iworm, jworm, ip, jp, ii, jj, kk, icell, jcell, i, ip1, jp1, scell, idir, scnab, inogo, icnab, jcnab

      do iworm = 1, nworms
         do ip = 1, np
            ii = (iworm - 1)*np + ip
            icell = 1 + floor(x(iworm, ip)/dcell)
            jcell = 1 + floor(y(iworm, ip)/dcell)
            !pid(icell-1,jcell-1)%row =
            if (icell .gt. nxcell .or. icell .lt. 1) then
               print *, 'nxcell=', nxcell
               print *, 'icell out of bounds', iworm, ip, x(iworm, ip)
               print *, 'icell=', icell, 'dcell=', dcell
               do i = 1, np
                  print *, 'check positions:', i, x(iworm, i)
               end do
               stop
            end if
            if (jcell .gt. nycell .or. jcell .lt. 1) then
               print *, 'nycell=', nycell
               print *, 'jcell out of bounds', iworm, ip, y(iworm, ip)
               print *, 'jcell=', jcell, 'dcell=', dcell
               do i = 1, np
                  print *, 'check positions:', i, x(iworm, i)
               end do
               stop
            end if
            scell = (icell) + (jcell - 1)*nxcell !1d-indexing for 2d cells
            if (scell .gt. ncells .or. scell .lt. 1) then
               print *, 'scell out of bounds, i=', i
            end if

            ipointto(ii) = hhead(scell)
            hhead(scell) = ii
         end do
      end do
      do icell = 1, nxcell
         do jcell = 1, nycell
            scell = icell + (jcell - 1)*nxcell
            if (hhead(scell) .ne. -1) then
               !there are particles in the cell called scell so
               !lets check all the neighbor cells

               do idir = 1, 9
                  icnab = icell + ddx(idir)
                  !if(icnab.gt.nxcell)goto 1303
                  if (icnab .gt. nxcell) cycle
                  !if(icnab.eq.0)goto 1303
                  if (icnab .eq. 0) cycle
                  jcnab = jcell + ddy(idir)
                  !if(jcnab.gt.nycell)goto 1303
                  if (jcnab .gt. nycell) cycle
                  !if(jcnab.eq.0)goto 1303
                  if (jcnab .eq. 0) cycle

                  scnab = icnab + (jcnab - 1)*nxcell !1d neighbor
                  if (hhead(scnab) .ne. -1) then
                     !there are particles in the cell called scnab
                     ii = hhead(scell) ! ii is the # of the head particle

                     do while (ii .gt. 0)
                        iworm = 1 + int((ii - 1)/np) ! find which worm ii is in
                        !print*, iworm
                        ip = ii - np*(iworm - 1) ! which particle in the worm is ii?
                        jj = hhead(scnab) ! head particle of neighboring cell

                        do while (jj .gt. 0)
                           jworm = 1 + int((jj - 1)/np)
                           !print*, jworm
                           jp = jj - np*(jworm - 1)
                           inogo = 0
                           if (iworm .eq. jworm .and. abs(ip - jp) .le. 2) then
                              inogo = 1 !on the same worm and close means no interaction calculted here
                           end if
                           if (ii .lt. jj .and. inogo .eq. 0) then

                              dddx = x(jworm, jp) - x(iworm, ip)
                              dddy = y(jworm, jp) - y(iworm, ip)
                              r2 = dddx**2 + dddy**2
                              riijj = dsqrt(r2)
                              ! add attractive force fdep between all pairs
                              if (r2 .le. r2cutsmall) then
                                 ffor = -48.d0*r2**(-7) + 24.d0*r2**(-4) + fdep/riijj
                                 ffx = ffor*dddx
                                 ffy = ffor*dddy
                                 fx(iworm, ip) = fx(iworm, ip) + ffx
                                 fx(jworm, jp) = fx(jworm, jp) - ffx
                                 fy(iworm, ip) = fy(iworm, ip) + ffy
                                 fy(jworm, jp) = fy(jworm, jp) - ffy

                                 !take these neighbors into account in calculating vxave and vyave

                                 vxave(iworm, ip) = vxave(iworm, ip) + vx(jworm, jp)
                                 vyave(iworm, ip) = vyave(iworm, ip) + vy(jworm, jp)
                                 nnab(iworm, ip) = nnab(iworm, ip) + 1
                                 vxave(jworm, jp) = vxave(jworm, jp) + vx(iworm, ip)
                                 vyave(jworm, jp) = vyave(jworm, jp) + vy(iworm, ip)
                                 nnab(jworm, jp) = nnab(jworm, jp) + 1

                                 ! add 'dogic drive' to interacting pairs
                                 !first calculate unit vectors along each worm
                                 ip1 = ip + 1
                                 if (ip1 .le. np) then
                                    dxi = x(iworm, ip1) - x(iworm, ip)
                                    dyi = y(iworm, ip1) - y(iworm, ip)
                                 else
                                    dxi = x(iworm, ip) - x(iworm, ip - 1)
                                    dyi = y(iworm, ip) - y(iworm, ip - 1)
                                 end if

                                 jp1 = jp + 1
                                 if (jp1 .le. np) then
                                    dxj = x(jworm, jp1) - x(jworm, jp)
                                    dyj = y(jworm, jp1) - y(jworm, jp)
                                 else
                                    dxj = x(jworm, jp) - x(jworm, jp - 1)
                                    dyj = y(jworm, jp) - y(jworm, jp - 1)
                                 end if
                                 !if the two vectors have any component pointing in opposite directions
                                 if (dxi*dxj + dyi*dyj .le. 0.d0) then
                                    !normalize those vectors to make them unit vectors
                                    ri = dsqrt(dxi*dxi + dyi*dyi)
                                    dxi = dxi/ri
                                    dyi = dyi/ri

                                    rj = dsqrt(dxj*dxj + dyj*dyj)
                                    dxj = dxj/rj
                                    dyj = dyj/rj
                                    !now they are both unit vectors. Find the direction for the force...

                                    dx = (dxi - dxj)/2.d0
                                    dy = (dyi - dyj)/2.d0

                                    !normalize

                                    r = dsqrt(dx*dx + dy*dy)
                                    dx = dx/r
                                    dy = dy/r

                                    !add an extra attractive component where kinesin drive is present

                                    ffx = fdogic*(dx) + 0.7d0*dddx/riijj
                                    ffy = fdogic*(dy) + 0.7d0*dddy/riijj

                                    !ffx=fdogic*(dx)
                                    !ffy=fdogic*(dy)

                                    fx(iworm, ip) = fx(iworm, ip) + ffx
                                    fx(jworm, jp) = fx(jworm, jp) - ffx
                                    fy(iworm, ip) = fy(iworm, ip) + ffy
                                    fy(jworm, jp) = fy(jworm, jp) - ffy

                                 end if
                              end if
                           end if
                           jj = ipointto(jj)
                        end do
                        ii = ipointto(ii)
                     end do
                  end if
               end do
            end if
         end do
      end do
   end subroutine cell_sort

   subroutine bead_fluid_interaction(nworms,np,res,eps)!,intvx,intvy)
      !use fluid_mod, only : vel,force
      integer(int64), intent(in) :: nworms,np,res
      real(real64), intent(in) :: eps
      real(real64) :: dcell2,cs2,inv_cs2
      integer(int64) :: iworm,ip,ii,icell,jcell
      integer(int64) :: i,j,k
      ! fluid resolution must be smaller than the
      cs2 = 1.0d0/3.0d0 
      inv_cs2 = 1/cs2
      force(:,:,1) = 0.0d0
      force(:,:,2) = 0.0d0
      intvx(:,:) = 0.0d0
      intvy(:,:) = 0.0d0
      dcell2 = (hx/float(nxcell*res))
      do iworm = 1, nworms
         do ip = 1, np
            ii = (iworm - 1)*np + ip
            icell = 1 + floor(x(iworm, ip)/dcell2)
            jcell = 1 + floor(y(iworm, ip)/dcell2)
            !calculating velocity of the fluid 
            intvx(icell,jcell) = intvx(icell,jcell) + vx(iworm,ip)
            intvy(icell,jcell) = intvx(icell,jcell) + vy(iworm,ip)
            force(icell,jcell,1) = force(icell,jcell,1) + eps*(vx(iworm,ip) - vel(icell,jcell,1))
            force(icell,jcell,2) = force(icell,jcell,2) + eps*(vy(iworm,ip) - vel(icell,jcell,2))
            !fx(iworm,ip) = fx(iworm,ip) + force(icell,jcell,1)
            !fy(iworm,ip) = fy(iworm,ip) + force(icell,jcell,2)
         enddo
      enddo
      ! calculate the source term to add into the lattice boltzmann eq
      do i = 0,nxcell
         do j = 0, nycell
            do k = 1,9
               S(k,i,j) = w(k)*( (force(i,j,1)*e(1,k) + force(i,j,2)*e(2,k))*inv_cs2 + &
               2*inv_cs2*force(i,j,1)*intvy(i,j)*(e(1,k)*e(1,k) - cs2) + 2*inv_cs2*force(i,j,2)*intvx(i,j)*(e(2,k)*e(2,k) - cs2))
            end do   
         end do
      end do
      !print*, sum(force(:,:,1)),sum(force(:,:,2))
   end subroutine bead_fluid_interaction

   subroutine finalize(fluid_cpl, calc_time, io_time)
      real(real64), intent(in) :: calc_time, io_time
      logical, intent(in) :: fluid_cpl
      deallocate (x)
      deallocate (y)
      deallocate (vx)
      deallocate (vy)
      deallocate (vxave)
      deallocate (vyave)
      deallocate (fx)
      deallocate (fy)
      deallocate (fxold)
      deallocate (fyold)
      if (fluid_cpl) then
         deallocate (intvx)
         deallocate (intvy)
      end if

      deallocate (savex)
      deallocate (savey)
      deallocate (ireverse)
      deallocate (ddx)
      deallocate (ddy)
      deallocate (hhead)
      deallocate (ipointto)
      deallocate (nnab)
      !dclose (unit=xyzfileunit)
      print *, "DONE!"
      print *, "Total time", calc_time + io_time
      print *, "Calculations", calc_time, calc_time/(calc_time + io_time)
      print *, "Input/Output", io_time, io_time/(calc_time + io_time)
   end subroutine finalize
end module active_mod
