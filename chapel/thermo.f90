 program thermo
      implicit real*8(a-h,o-z)
      real*8 x1,y1,x2,y2,length0,fx1,fy1,fx2,fy2
      real*8 fxold1,fyold1,fxold2,fyold2,zzero
      real*8 vx1,vy1,vx2,vy2,r2,r,kspring,dt,w,mass,period,pi
      integer np
!     choose parameters
      length0=4.d0
        kspring=7.d0
        mass=1.d0
        pi=4.d0*atan(1.d0)
        w=dsqrt(kspring/mass)
        period=2.d0*pi/w
        dt=period/100.d0
        nsteps=1000
        np=2
        zzero=0.d0  
!     set initial condition     

        x1=-0.5*length0-0.2
        x2=0.5*length0+0.2
        y1=0.d0
        y2=0.d0
        vx1=0.d0
        vy1=0.1d0
        vx2=0.d0
        vy2=-0.1d0  
!     force calculation before main loop
 
        fx1=0.d0
        fy1=0.d0
        fx2=0.d0
        fy2=0.d0
        dx=(x2-x1)
        dy=(y2-y1)
        r2=dx*dx+dy*dy
        r=dsqrt(r2)
        ff=-kspring*(r-length0)
        ffx=ff*dx/r
        ffy=ff*dy/r
        fx1=fx1-ffx
        fx2=fx2+ffx
        fy1=fy1-ffy
       fy2=fy2+ffy

!     output file setup and first data
      open(unit=10,file='twomass01.xyz',status='unknown')
        write(10,*)np  
        write(10,*)
        write(10,*)'A ',x1,y1,zzero
        write(10,*)'A ',x2,y2,zzero
      
!     main loop on steps

      do 999 istep=1,nsteps

!     update positions

      x1=x1+vx1*dt+0.5*(fx1/mass)*dt*dt
      y1=y1+vy1*dt+0.5*(fy1/mass)*dt*dt
      x2=x2+vx2*dt+0.5*(fx2/mass)*dt*dt
      y2=y2+vy2*dt+0.5*(fy2/mass)*dt*dt

!     store old forces

      fx1old=fx1
      fy1old=fy1
      fx2old=fx2
      fy2old=fy2

!     calculate forces

        fx1=0.d0
        fy1=0.d0
        fx2=0.d0
        fy2=0.d0
        dx=(x2-x1)
        dy=(y2-y1)
        r2=dx*dx+dy*dy
        r=dsqrt(r2)
        ff=-kspring*(r-length0)
        ffx=ff*dx/r
        ffy=ff*dy/r
        fx1=fx1-ffx
        fx2=fx2+ffx
        fy1=fy1-ffy
        fy2=fy2+ffy

!     update velocities

      vx1=vx1+0.5*(fx1+fx1old)*dt/mass
        vy1=vy1+0.5*(fy1+fy1old)*dt/mass
      vx2=vx2+0.5*(fx2+fx2old)*dt/mass
        vy2=vy2+0.5*(fy2+fy2old)*dt/mass 
            
        xmom=mass*vx1+mass*vx2
        ymom=mass*vy1+mass*vy2
        epot=0.5*kspring*(r-length0)*(r-length0)
        ekin=0.5*mass*(vx1*vx1+vy1*vy1+vx2*vx2+vy2*vy2)
        etot=epot+ekin
        write(6,1000)istep*dt,xmom,ymom,epot,ekin,etot
1000  format(f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4)

      if(mod(istep,10).eq.0)then
         write(10,*)np
         write(10,*)
         write(10,*)'A ',x1,y1,zzero
         write(10,*)'A ',x2,y2,zzero
      endif
  
999   continue
      close(unit=10)
      stop
end program thermo
