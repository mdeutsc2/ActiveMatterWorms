program amatter
   !import Fortran 2008 kinds
   use, intrinsic :: iso_fortran_env
   implicit none
   !implicit real*8(a-h,o-z)

   ! PARAMETERS
   integer(int64),parameter :: np=400
   integer(int64),parameter :: nworms=250
   integer(int64),parameter :: nsteps=2000000

   real(real64),parameter :: fdogic=0.06d0
   real(real64),parameter :: fdogicwall=0.0000d0
   real(real64),parameter :: fdep=1.0d0
   real(real64),parameter :: fdepwall=0.0d0
   real(real64),parameter :: diss=0.08d0
   real(real64),parameter :: dt=0.02d0
   real(real64),parameter :: kspring=57.146436d0
   real(real64),parameter :: kbend=40.0d0
   real(real64),parameter :: length0=0.8d0

   ! ARRAY ALLOCATION
   real(real64),dimension(nworms,np) ::  x,y,vx,vy,vxave,vyave,fx,fy,fxold,fyold
   real(real64),dimension(np) :: savex,savey
   integer(int64),dimension(nworms) :: ireverse
   integer(int64),dimension(9) :: ddx,ddy
   integer(int64),dimension(504*504) :: hhead
   integer(int64),dimension(nworms*np) :: ipointo
   integer(int64),dimension(nworms,np) :: nnab

   ! VARIABLES
   real(real64) :: rcut,r2cut,rcutsmall,r2cutsmall,rwall,r2inside,density,hx,hy,hxo2,hyo2
   integer(int64) :: nxcell,nycell,ncells,iseed,nstepso500,ireadsaved,iwalldrive
   real(real64) :: dcell,gnoise,dto2,dt2o2,length2,lengthmax,thetanow,a,rmin,t
   integer(int64) scell,scnab
   real(real64) twopi,pi,pio4
   !variables that are inside loops
   real(real64) :: r,dth,xangle,dx,dy,xang,rx,ry,dot,rsq,v1,v2,fac,g1,th,r2
   real(real64) :: ff,ffx,ffy,f4x,f4y,f3x,f3y,f2x,f2y,ffor
   real(real64) :: x2,x3,x4,y2,y3,y4,x23,x34,y23,y34,r23,r34,cosvalue,sinvalue
   real(real64) :: xwall,ywall,rr2,dxi,dyi,ri,rj,dxj,dyj,dddx,dddy,riijj
   integer(int64) :: ip1,jp1,i2,i3,i4,icell,jcell,icnab,jcnab,inogo,idir
   real(real64) :: rand1,rand2

   !iterator variables
   integer(int64) :: iw,iworm,jworm,ip,jp,i,itime,ii,jj

   ! IO variables and format statements
   integer(int32) :: xyzfileunit
   character(len=:), allocatable :: xyzfmt,datfmt
   xyzfmt = '(A2,1x,f8.4,1x,f8.4,1x,f8.4)'

   rcut=2.5d0
   r2cut=rcut*rcut
   rcutsmall=2.0d0**(1.0/6.0)
	r2cutsmall=rcutsmall**2
  rwall=125.d0*rcutsmall*dsqrt(2.d0)

	print *, 'rwall=',rwall
	twopi=8.d0*datan(1.d0)
  pi=0.5d0*twopi
  pio4=pi*0.25d0
  print *, 'nworms= ',nworms
	density= nworms*np/(pi*rwall**2)
	print *, 'density=',density


  hx=2.d0*rwall+1.d0
	hy=hx
  hyo2=0.5d0*hy
  hxo2=0.5d0*hx

  nxcell=(hx/rcutsmall)-1
  nycell=(hy/rcutsmall)-1
	dcell=hx/float(nxcell)
  !		write(6,*)'nxcell=',nxcell
  !		write(6,*)'dcell=',dcell
  !		write(6,*)'rcutsmall=',rcutsmall

  ncells=nxcell*nycell

  iseed=6487735
  call random_init(repeatable=.true.,image_distinct=.false.)
	!ss=rand(iseed)
	!ss=rand()
  nstepso500=nsteps/500


	open(unit=xyzfileunit,file='amatter.xyz',status='unknown')


  gnoise=.8d0/dsqrt(10.d0)*0.8d0
	dt2o2=dt*dt*0.5d0
	dto2=dt*0.5d0


  length2=2.d0*length0
	lengthmax=(length0)*float(np-1)

  !array for locating neighbor cells
  ddx(1)=1
  ddy(1)=0
  ddx(2)=1
  ddy(2)=1
  ddx(3)=0
  ddy(3)=1
  ddx(4)=-1
  ddy(4)=1
  ddx(5)=-1
  ddy(5)=0
  ddx(6)=-1
  ddy(6)=-1
  ddx(7)=0
  ddy(7)=-1
  ddx(8)=1
  ddy(8)=-1
  ddx(9)=0
  ddy(9)=0

  r2inside = (rwall-rcutsmall)**2
  thetanow=5.d0*pi
  a=0.24
  rmin=a*thetanow

  !setting up the worms
	do iw=1,nworms
    ireverse(iw)=0
    call random_number(rand1)
		if(rand1.le.0.5d0) then
     ireverse(iw)=1
    endif
    do i=1,np
       r=a*thetanow
       dth=length0/r
       thetanow=thetanow+dth
       x(iw,i)=hxo2+r*dcos(thetanow)
       y(iw,i)=hyo2+r*dsin(thetanow)
       xangle=atan2(y(iw,i)-hyo2,x(iw,i)-hxo2)
       !give them an initial velocity going around the circle
       !vx(iw,i)=-dsin(xangle)*0.2d0*float(ireverse(iw))
       !vy(iw,i)=dcos(xangle)*0.2d0*float(ireverse(iw))
       !vx(iw,i)=-dsin(xangle)*0.2d0
       !vy(iw,i)=dcos(xangle)*0.2d0
       vx(iw,i)=0.d0
       vy(iw,i)=0.d0
       fx(iw,i)=0.d0
       fy(iw,i)=0.d0
    enddo
    thetanow=thetanow+4.d0*dth
  enddo

  !reverse some of the worms and give them crazy colors
	do iw=1,nworms
    call random_number(rand1)
	  if(ireverse(iw).eq.1)then
      do i=1,np
         savex(i)=x(iw,i)
         savey(i)=y(iw,i)
      enddo
      do ip=1,np
         x(iw,ip)=savex(np+1-ip)
         y(iw,ip)=savey(np+1-ip)
      enddo
     endif
  enddo

  !ireadsaved=0
	!if(ireadsaved.eq.1)then
  !  open(unit=21,file='corah04.dat',status='old')
	!  read(21,*)nworms,np
	!	read(21,*)diss,gnoise,fdogic,rwall,fdep
	!  do iw=1,nworms
  !    do i=1,np
  !       read(21,datfmt)x(iw,i),y(iw,i),vx(iw,i),vy(iw,i)
  !    enddo
  !  enddo
  !  close(unit=21)
	!	r2inside=(rwall-rcutsmall)**2
  !endif

  write(xyzfileunit,*)nworms*np+4
	write(xyzfileunit,*)"# 0"
 	do iw=1,nworms
     dx=x(iw,1)-hxo2
     dy=y(iw,1)-hyo2
     xang=atan2(dy,dx)
     rx=-dsin(xang)
     ry=dcos(xang)
     dot=(x(iw,1)-x(iw,np))*rx+(y(iw,1)-y(iw,np))*ry
     if(dot.ge.0.d0)then
        do i=1,np
           write(xyzfileunit,xyzfmt)'A ',x(iw,i),y(iw,i),0.0d0
        enddo
     else
        do i=1,np
           write(xyzfileunit,xyzfmt)'B ',x(iw,i),y(iw,i),0.0d0
        enddo
     endif
  enddo

  write(xyzfileunit,xyzfmt)'E ',hxo2-rwall,hyo2-rwall,0.0d0
  write(xyzfileunit,xyzfmt)'E ',hxo2-rwall,hyo2+rwall,0.0d0
  write(xyzfileunit,xyzfmt)'E ',hxo2+rwall,hyo2-rwall,0.0d0
  write(xyzfileunit,xyzfmt)'E ',hxo2+rwall,hyo2+rwall,0.0d0

  iwalldrive=1


  !MAIN LOOP
	do itime=1,nsteps
    t=float(itime)*dt
    print*,itime
    !first update positions and store old forces
    call update_pos(nworms,np,x,y,vx,vy,fx,fy,fxold,fyold,dt,dt2o2)

    !calculate forces
    !zero out the force arrays and add Gaussian noise
    rsq = 0.0d0
    do iw=1,nworms
       do i=1,np
          do while(rsq.ge.0.999.or.rsq.le.0.001)
             call random_number(rand1)
             call random_number(rand2)
             v1=2.0*rand1-1.d0
             v2=2.0*rand2-1.d0
             rsq=v1*v1+v2*v2
          enddo
          call random_number(rand1)
          fac=dsqrt(-2.0*dlog(rsq)/rsq)
          g1=v1*fac*gnoise
          th=rand1*twopi
          fx(iw,i)=g1*dcos(th)
          fy(iw,i)=g1*dsin(th)
       enddo
    enddo

    !first set of springs nearest neighbor springs
    do iw=1,nworms
       do i=1,np-1
          ip1=i+1

          dx=x(iw,ip1)-x(iw,i)
          dy=y(iw,ip1)-y(iw,i)
          r=dsqrt(dx*dx+dy*dy)

          ff=-kspring*(r-length0)/r
          ffx=ff*dx
          ffy=ff*dy
          fx(iw,ip1)=fx(iw,ip1)+ffx
          fx(iw,i)=fx(iw,i)-ffx
          fy(iw,ip1)=fy(iw,ip1)+ffy
          fy(iw,i)=fy(iw,i)-ffy
       enddo
    enddo



    !bond bending terms
    do iw=1,nworms
       do i2=1,np-2
          i3=i2+1
          i4=i2+2
          !print*, i2,i3,i4
          x2=x(iw,i2)
          y2=y(iw,i2)
          x3=x(iw,i3)
          y3=y(iw,i3)
          x4=x(iw,i4)
          y4=y(iw,i4)
          y23=y3-y2
          y34=y4-y3
          x23=x3-x2
          x34=x4-x3
          r23=dsqrt(x23*x23+y23*y23)
          r34=dsqrt(x34*x34+y34*y34)

          cosvalue=(x23*x34+y23*y34)/(r23*r34)
          print*, x23,x34,y23,y34,r23,r34
          if (cosvalue .lt. 0.d0) then
            print*, x23,x34,y23,y34,r23,r34
            print*, "cosvalue .lt. 0.0d0"
            stop
          endif
          if(cosvalue.gt.1.d0) then
             cosvalue=1.d0
          endif
          sinvalue=dsqrt(1.d0-cosvalue*cosvalue)

          ff=-kbend*sinvalue/(r23*r34)

          dot=x23*x34+y23*y34
          fac=dot/(r23*r23)

          f2x=ff*(x34-fac*x23)
          f2y=ff*(y34-fac*y23)

          fac=dot/(r34*r34)
          f4x=ff*(fac*x34-x23)
          f4y=ff*(fac*y34-y23)
          f3x=-f2x-f4x
          f3y=-f2y-f4y


          fx(iw,i2)=fx(iw,i2)+f2x
          fy(iw,i2)=fy(iw,i2)+f2y

          fx(iw,i3)=fx(iw,i3)+f3x
          fy(iw,i3)=fy(iw,i3)+f3y

          fx(iw,i4)=fx(iw,i4)+f4x
          fy(iw,i4)=fy(iw,i4)+f4y
       enddo
    enddo



    !intra-worm interactions, 3rd nbors and beyond, short range repulsive interactions
    !adding fdep as well

!        do 290 iw=1,nworms
!         do 291 i1=1,np-3
!           do 292 i2=i1+3,np
!	      dx=x(iw,i1)-x(iw,i2)
!	      dy=y(iw,i1)-y(iw,i2)
!	      rr=dx*dx+dy*dy
!		  r=dsqrt(rr)

!	      if(rr.le.r2cutsmall)then
!   	        ff = (48.0/rr**7)-(24.0/rr**4)+fdep/r
!	        fx(iw,i1)=fx(iw,i1)+ff*dx
!		    fx(iw,i2)=fx(iw,i2)-ff*dx
!	        fy(iw,i1)=fy(iw,i1)+ff*dy
!		    fy(iw,i2)=fy(iw,i2)-ff*dy
!          endif
!      continue
!      continue
!      continue



    !put worm-wall interactions here, and dissipation force proportional to velocity
    do iw=1,nworms
       do i=1,np
          !dissipation proportional to v relative to local average
          fx(iw,i)=fx(iw,i)-diss*(vx(iw,i)-vxave(iw,i))
          fy(iw,i)=fy(iw,i)-diss*(vy(iw,i)-vyave(iw,i))
          !now that we have used them, zero out vxave and vyave, recalculate below
          vxave(iw,i)=vx(iw,i)
          vyave(iw,i)=vy(iw,i)
          nnab(iw,i)=1
          !calculate distance to the center
          dx=x(iw,i)-hxo2
          dy=y(iw,i)-hyo2
          r2=(dx*dx+dy*dy)
          !if close enough to the wall, calculate wall forces
          !use the short cut-off
          if(r2.ge.r2inside)then
             !find the nearest spot on the wall
             th=atan2(dy,dx)
             xwall=hxo2+rwall*dcos(th)
             ywall=hyo2+rwall*dsin(th)
             dx=xwall-x(iw,i)
             dy=ywall-y(iw,i)
             rr2=dx*dx+dy*dy
             r=dsqrt(rr2)
             !ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
             ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)
             fx(iw,i)=fx(iw,i)+ffor*dx
             fy(iw,i)=fy(iw,i)+ffor*dy

             !Turning on dogic drive with the wall!!!
             if(iwalldrive.eq.1)then
                !first calculate unit vector along the worm
                ip1=i+1
                if(ip1.le.np)then
                   dxi=x(iw,ip1)-x(iw,i)
                   dyi=y(iw,ip1)-y(iw,i)
                else
                   dxi=x(iw,i)-x(iw,i-1)
                   dyi=y(iw,i)-y(iw,i-1)
                endif
                !make it a unit vector
                ri=dsqrt(dxi*dxi+dyi*dyi)
                dxi=dxi/ri
                dyi=dyi/ri
                !calculate the unit vector along the wall
                dxj=-dsin(th)
                dyj= dcos(th)

                !if the vectors are not antiparallel, reverse the vector along the wall
                if(dxi*dxj+dyi*dyj.gt.0.d0)then
                   dxj=-dxj
                   dyj=-dyj
                endif

                !if the two vectors have any component pointing in opposite directions
                if(dxi*dxj+dyi*dyj.lt.0.d0)then
                   !Find the direction for the force...
                   dx=(dxi-dxj)/2.d0
                   dy=(dyi-dyj)/2.d0
                   !normalize the direction vector
                   ri=dsqrt(dx*dx+dy*dy)
                   dx=dx/ri
                   dy=dy/ri

                   !turn on extra-strong driving force
                   ffx=fdogicwall*dx
                   ffy=fdogicwall*dy
                   fx(iw,i)=fx(iw,i)+ffx
                   fy(iw,i)=fy(iw,i)+ffy
                endif
             endif
          endif
       enddo
    enddo

    !add calculation of vxave and vyave
    do ii=1,ncells
       hhead(ii)=-1
    enddo

    !assign particles to their cells
    do iworm=1,nworms
       do ip=1,np
          ii=(iworm-1)*np+ip
          icell=1+floor(x(iworm,ip)/dcell)
          jcell=1+floor(y(iworm,ip)/dcell)
          if(icell.gt.nxcell.or.icell.lt.1)then
             print *, 'nxcell=',nxcell
             print *, 'icell out of bounds',iworm,ip,x(iworm,ip)
             print *, 'icell=',icell, 'dcell=',dcell
             do i=1,np
                print *, 'check positions:',i,x(iworm,i)
             enddo
             stop
          endif
          if(jcell.gt.nycell.or.jcell.lt.1)then
             print *, 'nycell=',nycell
             print *, 'jcell out of bounds',iworm,ip,y(iworm,ip)
             print *, 'jcell=',jcell, 'dcell=',dcell
             do i = 1,np
                print *, 'check positions:',i,x(iworm,i)
             enddo
             stop
          endif
          scell=(icell)+(jcell-1)*nxcell !1d-indexing for 2d cells
          if(scell.gt.ncells.or.scell.lt.1)then
             print *, 'scell out of bounds, i=',i
          endif

          ipointo(ii)=hhead(scell)
          hhead(scell)=ii
       enddo
    enddo

    do icell=1,nxcell
       do jcell=1,nycell
          scell=icell+(jcell-1)*nxcell
          if(hhead(scell).ne.-1)then
             !there are particles in the cell called scell so
             !lets check all the neighbor cells

             do idir=1,9
                icnab=icell+ddx(idir)
                !if(icnab.gt.nxcell)goto 1303
                if(icnab.gt.nxcell) cycle
                !if(icnab.eq.0)goto 1303
                if(icnab.eq.0) cycle
                jcnab=jcell+ddy(idir)
                !if(jcnab.gt.nycell)goto 1303
                if(jcnab.gt.nycell) cycle
                !if(jcnab.eq.0)goto 1303
                if(jcnab.eq.0) cycle

                scnab=icnab+(jcnab-1)*nxcell !1d neighbor
                if(hhead(scnab).ne.-1)then
                   !there are particles in the cell called scnab
                   ii=hhead(scell) ! ii is the # of the head particle

                   do while(ii.gt.0)
                      iworm=1+int((ii-1)/np) ! find which worm ii is in
                      !print*, iworm
                      ip=ii-np*(iworm-1) ! which particle in the worm is ii?
                      jj=hhead(scnab) ! head particle of neighboring cell

                      do while(jj.gt.0)
                         jworm=1+int((jj-1)/np)
                         !print*, jworm
                         jp=jj-np*(jworm-1)
                         inogo=0
                         if(iworm.eq.jworm.and.abs(ip-jp).le.2) then
                            inogo=1 !on the same worm and close means no interaction calculted here
                         endif
                         if (ii.lt.jj.and.inogo.eq.0) then

                            dddx=x(jworm,jp)-x(iworm,ip)
                            dddy=y(jworm,jp)-y(iworm,ip)
                            r2=dddx**2+dddy**2
                            riijj=dsqrt(r2)
                            ! add attractive force fdep between all pairs
                            if(r2.le.r2cutsmall)then
                               ffor=-48.d0*r2**(-7)+24.d0*r2**(-4)+fdep/riijj
                               ffx=ffor*dddx
                               ffy=ffor*dddy
                               fx(iworm,ip)=fx(iworm,ip)+ffx
                               fx(jworm,jp)=fx(jworm,jp)-ffx
                               fy(iworm,ip)=fy(iworm,ip)+ffy
                               fy(jworm,jp)=fy(jworm,jp)-ffy

                               !take these neighbors into account in calculating vxave and vyave

                               vxave(iworm,ip)=vxave(iworm,ip)+vx(jworm,jp)
                               vyave(iworm,ip)=vyave(iworm,ip)+vy(jworm,jp)
                               nnab(iworm,ip)=nnab(iworm,ip)+1
                               vxave(jworm,jp)=vxave(jworm,jp)+vx(iworm,ip)
                               vyave(jworm,jp)=vyave(jworm,jp)+vy(iworm,ip)
                               nnab(jworm,jp)=nnab(jworm,jp)+1

                               ! add 'dogic drive' to interacting pairs
                               !first calculate unit vectors along each worm
                               ip1=ip+1
                               if(ip1.le.np)then
                                  dxi=x(iworm,ip1)-x(iworm,ip)
                                  dyi=y(iworm,ip1)-y(iworm,ip)
                               else
                                  dxi=x(iworm,ip)-x(iworm,ip-1)
                                  dyi=y(iworm,ip)-y(iworm,ip-1)
                               endif

                               jp1=jp+1
                               if(jp1.le.np)then
                                  dxj=x(jworm,jp1)-x(jworm,jp)
                                  dyj=y(jworm,jp1)-y(jworm,jp)
                               else
                                  dxj=x(jworm,jp)-x(jworm,jp-1)
                                  dyj=y(jworm,jp)-y(jworm,jp-1)
                               endif
                               !if the two vectors have any component pointing in opposite directions
                               if(dxi*dxj+dyi*dyj.le.0.d0)then
                                  !normalize those vectors to make them unit vectors
                                  ri=dsqrt(dxi*dxi+dyi*dyi)
                                  dxi=dxi/ri
                                  dyi=dyi/ri

                                  rj=dsqrt(dxj*dxj+dyj*dyj)
                                  dxj=dxj/rj
                                  dyj=dyj/rj
                                  !now they are both unit vectors. Find the direction for the force...

                                  dx=(dxi-dxj)/2.d0
                                  dy=(dyi-dyj)/2.d0

                                  !normalize

                                  r=dsqrt(dx*dx+dy*dy)
                                  dx=dx/r
                                  dy=dy/r

                                  !add an extra attractive component where kinesin drive is present

                                  ffx=fdogic*(dx)+0.7d0*dddx/riijj
                                  ffy=fdogic*(dy)+0.7d0*dddy/riijj

                                  !ffx=fdogic*(dx)
                                  !ffy=fdogic*(dy)


                                  fx(iworm,ip)=fx(iworm,ip)+ffx
                                  fx(jworm,jp)=fx(jworm,jp)-ffx
                                  fy(iworm,ip)=fy(iworm,ip)+ffy
                                  fy(jworm,jp)=fy(jworm,jp)-ffy

                               endif
                            endif
                         endif
                         jj=ipointo(jj)
                      enddo
                      ii=ipointo(ii)
                   enddo
                endif
             enddo
          endif
       enddo
    enddo

    !update velocities and normalize average velocities (they are one time step behind)
    call update_vel(nworms,np,vx,vy,vxave,vyave,fx,fy,fxold,fyold,nnab,dt2o2)

   if(mod(itime,2000).eq.0)then
       ! writing out the positions of the worms
       write(xyzfileunit,*)nworms*np+4
       write(xyzfileunit,*)"# 0"
       do iw=1,nworms
          dx=x(iw,1)-hxo2
          dy=y(iw,1)-hyo2
          xang=atan2(dy,dx)
          rx=-dsin(xang)
          ry=dcos(xang)
          dot=(x(iw,1)-x(iw,np))*rx+(y(iw,1)-y(iw,np))*ry
          if(dot.ge.0.d0)then
             do i=1,np
                write(xyzfileunit,xyzfmt)'A ',x(iw,i),y(iw,i),0.0d0
             enddo
          else
             do i=1,np
                write(xyzfileunit,xyzfmt)'B ',x(iw,i),y(iw,i),0.0d0
             enddo
          endif
       enddo
       write(xyzfileunit,xyzfmt)'E ',hxo2-rwall,hyo2-rwall,0.0d0
       write(xyzfileunit,xyzfmt)'E ',hxo2-rwall,hyo2+rwall,0.0d0
       write(xyzfileunit,xyzfmt)'E ',hxo2+rwall,hyo2-rwall,0.0d0
       write(xyzfileunit,xyzfmt)'E ',hxo2+rwall,hyo2+rwall,0.0d0
   endif


  enddo ! end of main loop
  close(unit=xyzfileunit)
  print *, "DONE!"

  contains
    subroutine update_pos(nworms,np,x,y,vx,vy,fx,fy,fxold,fyold,dt,dt2o2)
      integer(int64),intent(in) :: nworms,np
      real(real64),intent(in) :: vx(:,:),vy(:,:),fx(:,:),fy(:,:)
      real(real64),intent(in out) :: x(:,:),y(:,:),fxold(:,:),fyold(:,:)
      real(real64),intent(in) :: dt,dt2o2
      integer(int64) :: iw,i
      do iw=1,nworms
         do i=1,np
            x(iw,i)=x(iw,i)+vx(iw,i)*dt+fx(iw,i)*dt2o2
            y(iw,i)=y(iw,i)+vy(iw,i)*dt+fy(iw,i)*dt2o2
            fxold(iw,i)=fx(iw,i)
            fyold(iw,i)=fy(iw,i)
         enddo
      enddo
    end subroutine update_pos

    subroutine update_vel(nworms,np,vx,vy,vxave,vyave,fx,fy,fxold,fyold,nnab,dt2o2)
      integer(int64),intent(in) :: nworms,np
      real(real64), intent(in) :: dt2o2
      real(real64), intent(in out) :: vx(:,:),vy(:,:),vxave(:,:),vyave(:,:)
      real(real64), intent(in) :: fx(:,:),fy(:,:),fxold(:,:),fyold(:,:)
      integer(int64),intent(in) :: nnab(:,:)
      integer(int64) :: iw,i
      do iw=1,nworms
         do i=1,np
            vx(iw,i)=vx(iw,i)+dto2*(fx(iw,i)+fxold(iw,i))
            vy(iw,i)=vy(iw,i)+dto2*(fy(iw,i)+fyold(iw,i))
            vxave(iw,i)=vxave(iw,i)/nnab(iw,i)
            vyave(iw,i)=vyave(iw,i)/nnab(iw,i)
         enddo
      enddo
    end subroutine update_vel

end program amatter
