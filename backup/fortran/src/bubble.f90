	implicit real*8(a-h,o-z)
	parameter(np=18,nworms=60,nbubble=250,nsteps=4000000)
        parameter(iwait=1000000)
	dimension x(nworms,np),y(nworms,np),
     &   vx(nworms,np),vy(nworms,np),
     &   fx(nworms,np),fy(nworms,np)
       
	real*8 randomfmax(nworms)
	real*8 kspring,k2spring
	real*8 length0,length2,lengthmax,
     &  lengthnpm1,dt,t,mu,hy,hx,savex(np),savey(np)
	real*8 x0(nworms),y0(nworms),theta(nworms),w(nworms),twopi
        real*8 bubblex(nbubble),bubbley(nbubble) // particles in boundary
        real*8 bubblevx(nbubble),bubblevy(nbubble)
        real*8 bubblefx(nbubble),bubblefy(nbubble)
        real*8 lbubble,kbubble
        integer ddx(9),ddy(9)
	integer col(nworms)
        integer hhead(100*100)
        integer ipointo(nworms*np)
        integer scell,scnab
        
c       changing the driving force mechanism
c       these chains are driven along their length in the direction
c       of bond orientation in a polar fashion
c       extra random force on the head renewed at every delta t
c       there is also a little disorder in the initial config
c       size 100 by 100 with 560 chains


	twopi=8.d0*datan(1.d0)
        pi=0.5d0*twopi

        iseed=4381321
	ss=rand(iseed)
	ss=rand()

        nstepso1000=nsteps/1000

	open(unit=10,file='bubble3.xyz',status='unknown')
	npm1=np-1
        nbubblem1=nbubble-1
        npm2=np-2
	dt=0.015d0
	mu=0.15d0
	kspring=57.146436d0
        kbubble=kspring*0.5
	k2spring=kspring/2.d0
	knpm1spring=kspring

	length0=0.8d0
        length2=2.d0*length0
	lengthmax=(length0)*float(np-1)
        lbubble=1.d0
	zero=0.d0
	rcut=2.d0**(1./6.)
	r2cut=rcut*rcut

	hy=100.d0*rcut
	hyo2=0.5d0*hy

	hx=100.d0*rcut
	hxo2=0.5d0*hx
        nxcell=100
        nycell=100
        ncells=nxcell*nycell

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

        dcell=rcut

	iw=0
	do 200 ii=2,4
           shifter=rand()
	   do 201 jj=21,40
	      iw=iw+1
	      randomfmax(iw)=0.2d0
c	      randomfmax(iw)=0.2d0+0.1*(rand()-0.5)
	      x0(iw)=0.001+ii*hx/7.
	      y0(iw)=0.001+jj*hy/80.+shifter
	theta(iw)=0.d0
	col(iw)=1+int(rand()*5)
 201	continue
 200	continue
	do 2 iw=1,nworms
	do 1 i=1,np
	   x(iw,i)=x0(iw)+float(i-1)*length0*dcos(theta(iw))
	   if(x(iw,i).gt.hx)x(iw,i)=x(iw,i)-hx
	   if(x(iw,i).lt.0.d0)x(iw,i)=x(iw,i)+hx

	   y(iw,i)=y0(iw)+float(i-1)*length0*dsin(theta(iw))
	   if(y(iw,i).gt.hy)y(iw,i)=y(iw,i)-hy
	   if(y(iw,i).lt.0.d0)y(iw,i)=y(iw,i)+hy

	   vx(iw,i)=0.d0
	   vy(iw,i)=0.d0
	   fx(iw,i)=0.d0
	   fy(iw,i)=0.d0
1       continue
2	continue


c       reverse half of the worms
	do 799 iw=1,nworms
	   if(rand().le.0.5)then
c           if(rand().lt.0.d0)then
	do 801 i=1,np
	   savex(i)=x(iw,i)
	   savey(i)=y(iw,i)
801	continue

	do 802 ip=1,np
	   x(iw,ip)=savex(np+1-ip)
	   y(iw,ip)=savey(np+1-ip)
802	continue
	endif

799	continue


        rbubble=nbubble*lbubble/twopi
        write(6,*)'rbubble=',rbubble
        
        do 810 i=1,nbubble // puts the bubbles in space
        tth=i*twopi/float(nbubble)
        bubblex(i)=hxo2+rbubble*dcos(tth)
        bubbley(i)=hyo2+rbubble*dsin(tth)
c        write(6,*)bubblex(i),bubbley(i)
810     continue
        
        


	write(10,*)nworms*np+nbubble+4
	write(10,*)"# 0"
	do 5 iw=1,nworms

	   if(col(iw).eq.1)then
	do 6 i=1,np
	write(10,*)'A ',x(iw,i),y(iw,i),zero
6       continue

	elseif(col(iw).eq.2)then
	do 7 i=1,np
	write(10,*)'B ',x(iw,i),y(iw,i),zero
 7	continue	  

	elseif(col(iw).eq.3)then
	do 8 i=1,np
	write(10,*)'C ',x(iw,i),y(iw,i),zero
8	continue	  

	elseif(col(iw).eq.4)then
	do 9 i=1,np
	write(10,*)'D ',x(iw,i),y(iw,i),zero
9	continue	  

	elseif(col(iw).eq.5)then
	do 10 i=1,np
	write(10,*)'E ',x(iw,i),y(iw,i),zero
 10	continue	  
	  
	endif

5	continue

        do 11 ib=1,nbubble
          write(10,*)'F ',bubblex(ib),bubbley(ib),zero
11      continue



        write(10,*)'E ',zero,hy,zero
        write(10,*)'E ',hx,hy,zero
	write(10,*)'E ',zero,zero,zero
	write(10,*)'E ',hx,zero,zero

        
               

	do 999 itime=1,nsteps
	   if(mod(itime,100).eq.0)then
             write(6,*)itime
           endif
	t=float(itime)*dt
         if(itime.gt.iwait.and.mod(itime,2500).eq.0)then
          lbubble=lbubble*0.99
c          kbubble=kbubble/0.995
         endif

c       calculate forces

c       zero out the force arrays
	do 24 iw=1,nworms
	do 25 i=1,np
	fx(iw,i)=0.d0
	fy(iw,i)=0.d0
25       continue
24	continue

        do 26 ib=1,nbubble
        bubblefx(ib)=0.d0
        bubblefy(ib)=0.d0
26      continue       

c       first set of springs and driving forces
	do 99 iw=1,nworms
	do 110 i=1,np-1
        ip1=i+1
	dx=x(iw,ip1)-x(iw,i)
	if(dx.gt.hxo2)dx=dx-hx
	if(dx.lt.-hxo2)dx=dx+hx

	dy=y(iw,ip1)-y(iw,i)

	if(dy.gt.hyo2)dy=dy-hy
	if(dy.lt.-hyo2)dy=dy+hy

	r=dsqrt(dx*dx+dy*dy)
        fx(iw,i)=fx(iw,i)+randomfmax(iw)*dx/r
        fy(iw,i)=fy(iw,i)+randomfmax(iw)*dy/r
	ff=-kspring*(r-length0)/r
	ffx=ff*dx
	ffy=ff*dy
	fx(iw,ip1)=fx(iw,ip1)+ffx
	fx(iw,i)=fx(iw,i)-ffx
	fy(iw,ip1)=fy(iw,ip1)+ffy
	fy(iw,i)=fy(iw,i)-ffy
110      continue


c       spring from head to toe
c	dx=x(iw,np)-x(iw,1)
c	if(dx.gt.hxo2)dx=dx-hx
c	if(dx.lt.-hxo2)dx=dx+hx

c	dy=y(iw,np)-y(iw,1)

c	if(dy.gt.hyo2)dy=dy-hy
c	if(dy.lt.-hyo2)dy=dy+hy

c	r=dsqrt(dx*dx+dy*dy)
c	ff=-kextra*(r-lengthmax)/r
c	ffx=ff*dx
c	ffy=ff*dy
c	fx(iw,np)=fx(iw,np)+ffx
c	fx(iw,1)=fx(iw,1)-ffx
c	fy(iw,np)=fy(iw,np)+ffy
c	fy(iw,1)=fy(iw,1)-ffy


c       IMPORTANT
c       add extra random force on the head

        ffextra=rand()*5.d0
        tth=rand()*twopi
        fx(iw,np)=fx(iw,np)+ffextra*dcos(tth)
        fy(iw,np)=fy(iw,np)+ffextra*dsin(tth)



c       second set of springs

	do 12 i=1,npm2
        ip2=i+2
	dx=x(iw,ip2)-x(iw,i)
	if(dx.gt.hxo2)dx=dx-hx
	if(dx.lt.-hxo2)dx=dx+hx

	dy=y(iw,ip2)-y(iw,i)

	if(dy.gt.hyo2)dy=dy-hy
	if(dy.lt.-hyo2)dy=dy+hy

	r=dsqrt(dx*dx+dy*dy)
	ff=-k2spring*(r-length2)/r
	ffx=ff*dx
	ffy=ff*dy
	fx(iw,ip2)=fx(iw,ip2)+ffx
	fx(iw,i)=fx(iw,i)-ffx
	fy(iw,ip2)=fy(iw,ip2)+ffy
	fy(iw,i)=fy(iw,i)-ffy
12      continue


 99	continue

c       intra-worm interactions, 3rd nbors and beyond
        do 290 iw=1,nworms
         do 291 i1=1,np-3
           do 292 i2=i1+3,np
	      dx=x(iw,i1)-x(iw,i2)
	if(dx.gt.hxo2)dx=dx-hx
	if(dx.lt.-hxo2)dx=dx+hx

	      dy=y(iw,i1)-y(iw,i2)

	if(dy.gt.hyo2)dy=dy-hy
	if(dy.lt.-hyo2)dy=dy+hy

	      rr=dx*dx+dy*dy
	      if(rr.le.r2cut)then
   	        ff = (48.0/rr**7)-(24.0/rr**4)	      
	        fx(iw,i1)=fx(iw,i1)+ff*dx
		fx(iw,i2)=fx(iw,i2)-ff*dx
	        fy(iw,i1)=fy(iw,i1)+ff*dy
		fy(iw,i2)=fy(iw,i2)-ff*dy
            endif
292      continue
291      continue
290      continue

c       bubble-bubble interactions
        do 400 i=1,nbubble 
        ip1=i+1
        if(i.eq.nbubble)ip1=1
	dx=bubblex(ip1)-bubblex(i)
	if(dx.gt.hxo2)dx=dx-hx
	if(dx.lt.-hxo2)dx=dx+hx

	dy=bubbley(ip1)-bubbley(i)

	if(dy.gt.hyo2)dy=dy-hy
	if(dy.lt.-hyo2)dy=dy+hy

	r=dsqrt(dx*dx+dy*dy)
	ff=-kbubble*(r-lbubble)/r
	ffx=ff*dx
	ffy=ff*dy
	bubblefx(ip1)=bubblefx(ip1)+ffx
	bubblefx(i)=bubblefx(i)-ffx
	bubblefy(ip1)=bubblefy(ip1)+ffy
	bubblefy(i)=bubblefy(i)-ffy
400     continue

c       put worm-bubble interactions here


	do 1300 iw=1,nworms
	 do 1301 i=1,np
	  do 1310 ib=1,nbubble
              
	  dx=x(iw,i)-bubblex(ib)
  	  if(dx.gt.hxo2)dx=dx-hx !periodic boundaries
          if(dx.lt.-hxo2)dx=dx+hx
          dy=y(iw,i)-bubbley(ib)

	if(dy.gt.hyo2)dy=dy-hy
	if(dy.lt.-hyo2)dy=dy+hy

	      rr=dx*dx+dy*dy
	      if(rr.le.r2cut)then
    	        ff = (48.0/rr**7)-(24.0/rr**4)	      
	        fx(iw,i)=fx(iw,i)+ff*dx
		bubblefx(ib)=bubblefx(ib)-ff*dx
	        fy(iw,i)=fy(iw,i)+ff*dy
		bubblefy(ib)=bubblefy(ib)-ff*dy
             endif
1310	  continue
1301	 continue
1300	continue
        
        
      
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
            write(6,*)'icell out of bounds',iworm,ip,x(iworm,ip)
            write(6,*)'icell=',icell, 'dcell=',dcell
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
            if(icnab.gt.nxcell)icnab=1
            if(icnab.eq.0)icnab=nxcell            
            jcnab=jcell+ddy(idir)
            if(jcnab.gt.nycell)jcnab=1
            if(jcnab.eq.0)jcnab=nycell

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
               if(ii.lt.jj.and.iworm.ne.jworm)then
                dx=x(jworm,jp)-x(iworm,ip)
                if(dx.gt.hxo2)dx=dx-hx
                if(dx.lt.-hxo2)dx=dx+hx
                dy=y(jworm,jp)-y(iworm,ip)
                if(dy.gt.hyo2)dy=dy-hy
                if(dy.lt.-hyo2)dy=dy+hy
                r2=dx**2+dy**2
                if(r2.le.r2cut)then
                 ffor=-48.d0*r2**(-7)+24.d0*r2**(-4)
                 ffx=ffor*dx
                 ffy=ffor*dy
                 fx(iworm,ip)=fx(iworm,ip)+ffx
                 fx(jworm,jp)=fx(jworm,jp)-ffx
                 fy(iworm,ip)=fy(iworm,ip)+ffy
                 fy(jworm,jp)=fy(jworm,jp)-ffy
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
303         continue
            endif
c           "if(hhead(scell).ne.-1)"         

302      enddo
301       enddo


        


	do 19 iw=1,nworms
	do 20 i=1,np
	vx(iw,i)=mu*fx(iw,i)
	vy(iw,i)=mu*fy(iw,i)
20	continue
19	continue

        do 21 ib=1,nbubble
        bubblevx(ib)=mu*bubblefx(ib)
        bubblevy(ib)=mu*bubblefy(ib)
21      continue

	do 29 iw=1,nworms
 	do 30 i=1,np
 	x(iw,i)=x(iw,i)+vx(iw,i)*dt
	if(x(iw,i).gt.hx)x(iw,i)=x(iw,i)-hx
	if(x(iw,i).lt.0.d0)x(iw,i)=x(iw,i)+hx
 	y(iw,i)=y(iw,i)+vy(iw,i)*dt
	if(y(iw,i).gt.hy)y(iw,i)=y(iw,i)-hy
	if(y(iw,i).lt.0.d0)y(iw,i)=y(iw,i)+hy
 30     continue
 29	continue

        if(itime.gt.iwait)then
        do 31 ib=1,nbubble
        bubblex(ib)=bubblex(ib)+bubblevx(ib)*dt
	if(bubblex(ib).gt.hx)bubblex(ib)=bubblex(ib)-hx
	if(bubblex(ib).lt.0.d0)bubblex(ib)=bubblex(ib)+hx
        bubbley(ib)=bubbley(ib)+bubblevy(ib)*dt
	if(bubbley(ib).gt.hy)bubbley(ib)=bubbley(ib)-hy
	if(bubbley(ib).lt.0.d0)bubbley(ib)=bubbley(ib)+hy
31      continue
        endif

	if(mod(itime,nstepso1000).eq.0)then

	write(10,*)np*nworms+nbubble
	write(10,*)"#",itime

	do 39 iw=1,nworms

        if(col(iw).eq.1)then
	do 506 i=1,np
	write(10,*)'A ',x(iw,i),y(iw,i),zero
506     continue

	elseif(col(iw).eq.2)then
	do 507 i=1,np
	write(10,*)'B ',x(iw,i),y(iw,i),zero
 507	continue	  

	elseif(col(iw).eq.3)then
	do 508 i=1,np
	write(10,*)'C ',x(iw,i),y(iw,i),zero
508	continue	  

	elseif(col(iw).eq.4)then
	do 509 i=1,np
	write(10,*)'D ',x(iw,i),y(iw,i),zero
509	continue	  

	elseif(col(iw).eq.5)then
	do 510 i=1,np
	write(10,*)'E ',x(iw,i),y(iw,i),zero
510	continue	  

	endif
39	continue

        do 520 ib=1,nbubble
        write(10,*)'F ',bubblex(ib),bubbley(ib),zero
520     continue
       




        endif
	   

999     continue



	close(unit=10)

	stop
	end
	
