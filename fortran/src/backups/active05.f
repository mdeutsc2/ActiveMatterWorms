program active05 
	  implicit real*8(a-h,o-z)
	  parameter(np=400,nworms=250,nsteps=2000000)

	  dimension x(nworms,np),y(nworms,np)
      dimension vx(nworms,np),vy(nworms,np)
      dimension vxave(nworms,np),vyave(nworms,np)
      dimension fx(nworms,np),fy(nworms,np)
	  dimension fxold(nworms,np),fyold(nworms,np)

	  real*8 kspring,kbend
	  real*8 length0,length2,lengthmax,dt,t,hy,hx
	  real*8 savex(np),savey(np)
	  real*8 twopi,pi
	integer ireverse(nworms)
	integer ddx(9),ddy(9)

	integer hhead(504*504)
	integer ipointo(nworms*np)
	integer scell,scnab
	integer nnab(nworms,np)
	character*1 col(nworms)
	character*1 collist(10)




c       setting up dissipation relative to local average velocity summed over neighbors



c       SWITCHING ON the Dogic drive at the wall
c       SWITCHING OFF depletion force at the wall, just short range attraction
c       using higher value of kinesin drive fdogic
c       increasing the density
c       adding fdep to intraworm

c
c       switching to velocity verlet
	  zero=0.d0
c	  rcut=2.d0**(1./6.)
        rcut=2.5d0
        r2cut=rcut*rcut
	    rcutsmall=2.d0**(1./6.)
	    r2cutsmall=rcutsmall**2

            rwall=125.d0*rcutsmall*dsqrt(2.d0)

		  write(6,*)'rwall=',rwall

	    twopi=8.d0*datan(1.d0)
        pi=0.5d0*twopi
		pio4=pi*0.25d0



		 write(6,*)'nworms= ',nworms


		  density= nworms*np/(pi*rwall**2)
		  write(6,*)'density=',density


          hx=2.d0*rwall+1.d0
		  hy=hx



		hyo2=0.5d0*hy
		hxo2=0.5d0*hx




        nxcell=(hx/rcutsmall)-1
        nycell=(hy/rcutsmall)-1
		dcell=hx/float(nxcell)
c		write(6,*)'nxcell=',nxcell
c		write(6,*)'dcell=',dcell
c		write(6,*)'rcutsmall=',rcutsmall

        ncells=nxcell*nycell




        collist(1)='A'
        collist(2)='B'
        collist(3)='C'
        collist(4)='D'
        collist(5)='E'
        collist(6)='F'
        collist(7)='G'
        collist(8)='H'
	    collist(9)='I'
        collist(10)='J'


c         18a
               iseed=6487735

	   ss=rand(iseed)
	   ss=rand()

        nstepso500=nsteps/500


	   open(unit=10,file='active05.xyz',status='unknown')
	   open(unit=16,file='active05-debug.dat',status='unknown')
	    npm1=np-1

        npm2=np-2




          fdogic=0.06d0
c		fdogicwall=0.06d0
		fdogicwall=0.04d0

		 fdep=1.d0
c         fdep=1.d0
		 fdepwall=0.d0

c          diss=0.004
           diss=0.08
		gnoise=.80/dsqrt(10.d0)
                gnoise=gnoise*0.8d0




          dt=0.02
c           dt=0.02
	dt2o2=dt*dt*0.5d0
	dto2=dt*0.5d0

	    kspring=57.146436d0
        kbend=40.0d0


	length0=0.8d0
        length2=2.d0*length0
	lengthmax=(length0)*float(np-1)





c       array for locating neighbor cells
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




		r2inside=(rwall-rcutsmall)**2


	  thetanow=5.d0*pi

c	  a=0.24
      a=0.16
        rmin=a*thetanow




		do 210 iw=1,nworms
		ireverse(iw)=0
		if(rand().le.0.5d0)ireverse(iw)=1

       do 211 i=1,np

 	   r=a*thetanow
       dth=length0/r
	   thetanow=thetanow+dth
	   x(iw,i)=hxo2+r*dcos(thetanow)
	   y(iw,i)=hyo2+r*dsin(thetanow)
	   xangle=atan2(y(iw,i)-hyo2,x(iw,i)-hxo2)
c       give them an initial velocity going around the circle
c	   vx(iw,i)=-dsin(xangle)*0.2d0*float(ireverse(iw))
c	   vy(iw,i)=dcos(xangle)*0.2d0*float(ireverse(iw))
c	   vx(iw,i)=-dsin(xangle)*0.2d0
c	   vy(iw,i)=dcos(xangle)*0.2d0
       vx(iw,i)=0.d0
	   vy(iw,i)=0.d0
	   fx(iw,i)=0.d0
	   fy(iw,i)=0.d0
211   continue
      thetanow=thetanow+4.d0*dth


210   continue

c       reverse some of the worms and give them crazy colors
	do 799 iw=1,nworms
			col(iw)=collist(1+int(rand()*5))
c      rr=dsqrt((x(iw,1)-hxo2)**2+(y(iw,1)-hyo2)**2)
c	  if(rr.gt.0.95*rwall)then
c	     ireverse(iw)=0
c	  else
c	     ireverse(iw)=1
c	  endif
	  if(ireverse(iw).eq.1)then
	  do 801 i=1,np
	   savex(i)=x(iw,i)
	   savey(i)=y(iw,i)
801	  continue

	  do 802 ip=1,np
	   x(iw,ip)=savex(np+1-ip)
	   y(iw,ip)=savey(np+1-ip)
802	  continue
	  endif

799	  continue

        ireadsaved=0
		nworms1=0
		if(ireadsaved.eq.1)then
		open(unit=21,file='active05.dat',status='old')
	    read(21,*)nworms1,nnp
		read(21,*)diss,gnoise,fdogic,rwall,fdep
	    do 1701 iw=1,nworms
	    do 1702 i=1,np
	    read(21,7000)x(iw,i),y(iw,i),vx(iw,i),vy(iw,i)

1702     continue
1701     continue
        close(unit=21)
		r2inside=(rwall-rcutsmall)**2

c		 write(6,1000)'A ',x(1,1),y(1,1),zero
1000     format(A2,1x,f8.4,1x,f8.4,1x,f8.4)

        endif


	write(10,*)nworms*np+4
	write(10,*)"# 0"
c	do 5 iw=1,nworms
c	  do 6 ip=1,np
c	  write(10,*)col(iw),x(iw,ip),y(iw,ip),zero
c6     continue
c5     continue


 	  do 5 iw=1,nworms
      dx=x(iw,1)-hxo2
	  dy=y(iw,1)-hyo2
      xang=atan2(dy,dx)
      rx=-dsin(xang)
      ry=dcos(xang)
      dot=(x(iw,1)-x(iw,np))*rx+(y(iw,1)-y(iw,np))*ry
      if(dot.ge.0.d0)then
 	  do 6 i=1,np
	  write(10,1000)'A ',x(iw,i),y(iw,i),zero
6     continue
	  else
	   do 7 i=1,np
	   write(10,1000)'B ',x(iw,i),y(iw,i),zero
7	   continue
	  endif
5	continue





       write(10,1000)'E ',hxo2-rwall,hyo2-rwall,zero
       write(10,1000)'E ',hxo2-rwall,hyo2+rwall,zero
	   write(10,1000)'E ',hxo2+rwall,hyo2-rwall,zero
	   write(10,1000)'E ',hxo2+rwall,hyo2+rwall,zero


	  iwalldrive=1


c   MAIN LOOP



	    do 999 itime=1,nsteps

	   if(mod(itime,200).eq.0)then
             write(6,*)itime
           endif
	t=float(itime)*dt

c   first update positions and store old forces



	do 29 iw=1,nworms

 	do 30 i=1,np

 	  x(iw,i)=x(iw,i)+vx(iw,i)*dt+fx(iw,i)*dt2o2

c	  if(x(iw,i).lt.0.d0.or.x(iw,i).gt.hx)then
c	      write(6,*)'x out of bounds'
c		  write(6,*)'iw = ',iw
c		  write(6,*)'i=',i
c		  close(unit=10)
c          stop
c      endif
   	  y(iw,i)=y(iw,i)+vy(iw,i)*dt+fy(iw,i)*dt2o2
c	  	  if(y(iw,i).lt.0.d0.or.y(iw,i).gt.hy)then
c	      write(6,*)'y out of bounds'
c		  write(6,*)'iw = ',iw
c		  write(6,*)'i=',i
c		  close(unit=10)
c          stop
c      endif
      fxold(iw,i)=fx(iw,i)
	  fyold(iw,i)=fy(iw,i)
 30     continue

 29	continue



c       calculate forces

c       zero out the force arrays and add Gaussian noise
	do 24 iw=1,nworms
	do 25 i=1,np
888      continue
	   v1=2.0*rand()-1.d0
	   v2=2.0*rand()-1.d0
	   rsq=v1*v1+v2*v2
	   if(rsq.ge.0.999.or.rsq.le.0.001)goto 888
	   fac=dsqrt(-2.0*dlog(rsq)/rsq)
	   g1=v1*fac*gnoise
c	   g2=v2*fac

       th=rand()*twopi
	   fx(iw,i)=g1*dcos(th)
	   fy(iw,i)=g1*dsin(th)
25       continue
24	continue


c       first set of springs
	do 99 iw=1,nworms
	do 110 i=1,np-1
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
110      continue

 99	continue



c         bond bending terms

         do 120 iw=1,nworms
		  do 121 i2=1,np-2
                        i3=i2+1
			i4=i2+2
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
                         if(cosvalue.gt.1.d0)cosvalue=1.d0
                        sinvalue=dsqrt(1.d0-cosvalue*cosvalue)
c			thetanow=dacos(cosvalue)
c                         sinvalue=(x23*y34-y23*x34)/(r23*r34)

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



121      continue
120    continue



c       intra-worm interactions, 3rd nbors and beyond, short range repulsive interactions
c        adding fdep as well

c        do 290 iw=1,nworms
c         do 291 i1=1,np-3
c           do 292 i2=i1+3,np
c	      dx=x(iw,i1)-x(iw,i2)
c	      dy=y(iw,i1)-y(iw,i2)
c	      rr=dx*dx+dy*dy
c		  r=dsqrt(rr)

c	      if(rr.le.r2cutsmall)then
c   	        ff = (48.0/rr**7)-(24.0/rr**4)+fdep/r
c	        fx(iw,i1)=fx(iw,i1)+ff*dx
c		    fx(iw,i2)=fx(iw,i2)-ff*dx
c	        fy(iw,i1)=fy(iw,i1)+ff*dy
c		    fy(iw,i2)=fy(iw,i2)-ff*dy
c          endif
c292      continue
c291      continue
c290      continue



c       put worm-wall interactions here, and dissipation force proportional to velocity


	do 1300 iw=1,nworms
	 do 1301 i=1,np
c     dissipation proportional to v relative to local average
      fx(iw,i)=fx(iw,i)-diss*(vx(iw,i)-vxave(iw,i))
	  fy(iw,i)=fy(iw,i)-diss*(vy(iw,i)-vyave(iw,i))
c      now that we have used them, zero out vxave and vyave, recalculate below
	  vxave(iw,i)=vx(iw,i)
	  vyave(iw,i)=vy(iw,i)
	  nnab(iw,i)=1
c     calculate distance to the center
	  dx=x(iw,i)-hxo2
      dy=y(iw,i)-hyo2

	  r2=(dx*dx+dy*dy)
c     if close enough to the wall, calculate wall forces
c     use the short cut-off
      if(r2.ge.r2inside)then
c     find the nearest spot on the wall
		  	th=atan2(dy,dx)
			xwall=hxo2+rwall*dcos(th)
			ywall=hyo2+rwall*dsin(th)
			dx=xwall-x(iw,i)
			dy=ywall-y(iw,i)
			rr2=dx*dx+dy*dy
			r=dsqrt(rr2)
c			ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)+fdepwall/r
			ffor=-48.d0*rr2**(-7)+24.d0*rr2**(-4)
            fx(iw,i)=fx(iw,i)+ffor*dx
			fy(iw,i)=fy(iw,i)+ffor*dy

c     Turning on dogic drive with the wall!!!
           if(iwalldrive.eq.1)then

cc               first calculate unit vector along the worm
				ip1=i+1
				if(ip1.le.np)then
	            dxi=x(iw,ip1)-x(iw,i)
 	            dyi=y(iw,ip1)-y(iw,i)
				else
				dxi=x(iw,i)-x(iw,i-1)
				dyi=y(iw,i)-y(iw,i-1)
				endif

cc               make it a unit vector
				ri=dsqrt(dxi*dxi+dyi*dyi)
				dxi=dxi/ri
				dyi=dyi/ri
cc               calculate the unit vector along the wall
                dxj=-dsin(th)
 				dyj= dcos(th)

cc               if the vectors are not antiparallel, reverse the vector along the wall
				if(dxi*dxj+dyi*dyj.gt.0.d0)then
                   dxj=-dxj
                   dyj=-dyj
                 endif

cc               if the two vectors have any component pointing in opposite directions

                 if(dxi*dxj+dyi*dyj.lt.0.d0)then

cc               Find the direction for the force...

                 dx=(dxi-dxj)/2.d0
                 dy=(dyi-dyj)/2.d0
c                normalize the direction vector
                 ri=dsqrt(dx*dx+dy*dy)
				 dx=dx/ri
				 dy=dy/ri

c      turn on extra-strong driving force
                 ffx=fdogicwall*dx
                 ffy=fdogicwall*dy
                 fx(iw,i)=fx(iw,i)+ffx
				 fy(iw,i)=fy(iw,i)+ffy
				 endif
			endif
         endif
1301	 continue
1300	continue

c        add calculation of vxave and vyave

c         write(16,*)'forces flag 4: ',fx(313,18),fy(313,18)
c       every cell is initially empty; hhead is -1
         do 381 ii=1,ncells
           hhead(ii)=-1
381      continue



c       assign particles to their cells
        do 382 iworm=1,nworms
          do 383 ip=1,np
             ii=(iworm-1)*np+ip
c           iwormcalc=1+int((ii-1)/np)
c           ipcalc=ii-np*(iworm-1)
c            write(6,*)'compare iworm: ',iworm,iwormcalc
c            write(6,*)'compare ip: ',ip,ipcalc


           icell=1+floor(x(iworm,ip)/dcell)
           jcell=1+floor(y(iworm,ip)/dcell)
           if(icell.gt.nxcell.or.icell.lt.1)then
            write(6,*)'nxcell=',nxcell
            write(6,*)'icell out of bounds',iworm,ip,x(iworm,ip)
            write(6,*)'icell=',icell, 'dcell=',dcell

			do 1501 i=1,np
			  write(6,*)'check positions:',i,x(iworm,i)
1501        continue
            stop
           endif
           if(jcell.gt.nycell.or.jcell.lt.1)then
            write(6,*)'jcell out of bounds',iworm,ip,y(iworm,ip)
            stop
           endif
           scell=(icell)+(jcell-1)*nxcell

           if(scell.gt.ncells.or.scell.lt.1)then
               write(6,*)'scell out of bounds, i=',i
           endif

           ipointo(ii)=hhead(scell)
           hhead(scell)=ii
383      continue
382      continue


       do 301 icell=1,nxcell
       do 302 jcell=1,nycell
         scell=icell+(jcell-1)*nxcell
         if(hhead(scell).ne.-1)then
c          there are particles in the cell called scell so
c          let's check all the neighbor cells

           do 303 idir=1,9
            icnab=icell+ddx(idir)
            if(icnab.gt.nxcell)goto 1303
            if(icnab.eq.0)goto 1303
            jcnab=jcell+ddy(idir)
            if(jcnab.gt.nycell)goto 1303
            if(jcnab.eq.0)goto 1303

            scnab=icnab+(jcnab-1)*nxcell
            if(hhead(scnab).ne.-1)then
c           there are particles in the cell called scnab

             ii=hhead(scell)

             do 310 while(ii>0)
             iworm=1+int((ii-1)/np)
             ip=ii-np*(iworm-1)
              jj=hhead(scnab)

              do 311 while(jj>0)
              jworm=1+int((jj-1)/np)
              jp=jj-np*(jworm-1)
c              write(6,*)'interact: ',ii,jj
               inogo=0
			   if(iworm.eq.jworm.and.abs(ip-jp).le.2)inogo=1
               if(ii.lt.jj.and.inogo.eq.0)then

                dddx=x(jworm,jp)-x(iworm,ip)
                dddy=y(jworm,jp)-y(iworm,ip)
                r2=dddx**2+dddy**2
				riijj=dsqrt(r2)



c               add attractive force fdep between all pairs
                if(r2.le.r2cutsmall)then

				 ffor=-48.d0*r2**(-7)+24.d0*r2**(-4)+fdep/riijj
                 ffx=ffor*dddx
                 ffy=ffor*dddy



                 fx(iworm,ip)=fx(iworm,ip)+ffx
                 fx(jworm,jp)=fx(jworm,jp)-ffx
                 fy(iworm,ip)=fy(iworm,ip)+ffy
                 fy(jworm,jp)=fy(jworm,jp)-ffy

c				 take these neighbors into account in calculating vxave and vyave

		vxave(iworm,ip)=vxave(iworm,ip)+vx(jworm,jp)
		vyave(iworm,ip)=vyave(iworm,ip)+vy(jworm,jp)
				 nnab(iworm,ip)=nnab(iworm,ip)+1
				 vxave(jworm,jp)=vxave(jworm,jp)+vx(iworm,ip)
				 vyave(jworm,jp)=vyave(jworm,jp)+vy(iworm,ip)
				 nnab(jworm,jp)=nnab(jworm,jp)+1


c                add 'dogic drive' to interacting pairs

c               first calculate unit vectors along each worm
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
c               if the two vectors have any component pointing in opposite directions
				if(dxi*dxj+dyi*dyj.le.0.d0)then


c               normalize those vectors to make them unit vectors
				ri=dsqrt(dxi*dxi+dyi*dyi)
				dxi=dxi/ri
				dyi=dyi/ri

				rj=dsqrt(dxj*dxj+dyj*dyj)
				dxj=dxj/rj
				dyj=dyj/rj
c               now they are both unit vectors. Find the direction for the force...

                dx=(dxi-dxj)/2.d0
                dy=(dyi-dyj)/2.d0

c               normalize

                r=dsqrt(dx*dx+dy*dy)
				dx=dx/r
				dy=dy/r

c                add an extra attractive component where kinesin drive is present

                ffx=fdogic*(dx)+0.7d0*dddx/riijj
                ffy=fdogic*(dy)+0.7d0*dddy/riijj

c                ffx=fdogic*(dx)
c                ffy=fdogic*(dy)


                 fx(iworm,ip)=fx(iworm,ip)+ffx
                 fx(jworm,jp)=fx(jworm,jp)-ffx
                 fy(iworm,ip)=fy(iworm,ip)+ffy
                 fy(jworm,jp)=fy(jworm,jp)-ffy

				 endif


                endif
c               "if(r2.le.r2cut)"
               endif
c              "if(ii.lt.jj)"
               jj=ipointo(jj)
311           enddo
             ii=ipointo(ii)
310          enddo
             endif
c            "if(hhead(scnab).ne.-1)"
1303         continue
303         continue
            endif
c           "if(hhead(scell).ne.-1)"

302      enddo
301       enddo


c        write(16,*)'forces flag end: ',fx(313,18),fy(313,18)


c   update velocities and normalize average velocities (they are one time step behind)




	do 19 iw=1,nworms

	do 20 i=1,np

	vx(iw,i)=vx(iw,i)+dto2*(fx(iw,i)+fxold(iw,i))
	vy(iw,i)=vy(iw,i)+dto2*(fy(iw,i)+fyold(iw,i))
	vxave(iw,i)=vxave(iw,i)/nnab(iw,i)
	vyave(iw,i)=vyave(iw,i)/nnab(iw,i)

20	continue

19	continue





	if(mod(itime,2000).eq.0)then

c	write(10,*)nworms*np
c	write(10,*)"# 0"
c	do 35 iw=1,nworms
c	  do 36 ip=1,np
c	  write(10,*)col(iw),x(iw,ip),y(iw,ip),zero
c36     continue
c35     continue

	write(10,*)nworms*np+4
	write(10,*)"# 0"
	do 205 iw=1,nworms

      dx=x(iw,1)-hxo2
	  dy=y(iw,1)-hyo2
      xang=atan2(dy,dx)
      rx=-dsin(xang)
	  ry=dcos(xang)
       dot=(x(iw,1)-x(iw,np))*rx+(y(iw,1)-y(iw,np))*ry
       if(dot.ge.0.d0)then
 	  do 206 i=1,np
	  write(10,1000)'A ',x(iw,i),y(iw,i),zero
206     continue
	  else
	   do 207 i=1,np
	   write(10,1000)'B ',x(iw,i),y(iw,i),zero
207	   continue
	  endif
205	  continue
       write(10,1000)'E ',hxo2-rwall,hyo2-rwall,zero
       write(10,1000)'E ',hxo2-rwall,hyo2+rwall,zero
	   write(10,1000)'E ',hxo2+rwall,hyo2-rwall,zero
	   write(10,1000)'E ',hxo2+rwall,hyo2+rwall,zero

      endif

	  if(mod(itime,5000).eq.0)then
	   open(unit=21,file='active04.dat',status='unknown')
	    write(21,*)nworms,np
		write(21,*)diss,gnoise,fdogic,rwall,fdep
	    do 901 iw=1,nworms
	    do 902 i=1,np
	    write(21,7000)x(iw,i),y(iw,i),vx(iw,i),vy(iw,i)

902     continue
901     continue
        close(unit=21)
	  endif

999     continue




	close(unit=10)


	    open(unit=21,file='active04.dat',status='unknown')
	    write(21,*)nworms,np
		write(21,*)diss,gnoise,fdogic,rwall,fdep
	    do 701 iw=1,nworms
	    do 702 i=1,np
	    write(21,7000)x(iw,i),y(iw,i),vx(iw,i),vy(iw,i)
7000    format(1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6)
702     continue
701     continue
        close(unit=21)



	stop
	end
