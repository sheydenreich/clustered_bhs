        implicit real*8 (A-H,O-Z)
        dimension a(-100:4050,-100:4050), b(-100:4050,-100:4050)
	DIMENSION xa1(100000), xa2(100000),eps(100000)
	dimension rowax1(20000),rowax2(20000),roway1(20000),roway2(20000)
	dimension colax1(20000),colax2(20000),colay1(20000),colay2(20000)
  dimension xh1(20,1000), xh2(20,1000)


	character(len=100)::inname
	character(len=100)::outname
	character(len=100)::paramname
  character(len=20)::thalonumber
  character(len=20)::tconcentration

	call get_command_argument(1,inname)
	call get_command_argument(2,outname)
  call get_command_argument(3,paramname)

  call get_command_argument(4,thalonumber)
  call get_command_argument(5,tconcentration)


  read(thalonumber,*) halonumber
  read(tconcentration,*) concentration

  thres = 2*2*4*4

	counter1 = 0
	counter2 = 0

	open(1,file=inname) ! Input file with stars coordinates and masses (x1,x2,mass)
	open(2,file=paramname) ! Input file with the pattern parameters
	open(3,file=outname) ! Output pattern

  open(4,file='testhalo.dat')
  do j=1,halonumber
        do i=1,concentration
          read(4,*) rx1,rx2
          xh1(j,i) = rx1
          xh2(j,i) = rx2
        end do
  end do

	! Input parameters
	read(2,*) yl,ny ! half size of the pattern side in Einstein radius units; size in pixels of the pattern
	read(2,*) xl,nx ! half size of the shooting region in Einstein radius units; number of rays to be shoot in one dimension
	read(2,*) rkstars,rks,rgs ! kappa and gamma of the smooth distribution of matter
	read(2,*) threshold ! threshold to check the connectivity of the tetragons

  

	allow=1.0
	! rkstars=0.019
	rktot=rkstars+rks
	xl1=allow*1.5*yl/abs(1-rktot-rgs)
	xl2=1.5*yl/abs(1-rktot+rgs)

	! Computation of the number of rays per unlensed pixel
	ns=ny-1
	xs=2*xl/(nx-1)
	ys=2*yl/(ny-1)
	raypix=float(nx-1)*float(nx-1)/xl/xl/float(ny-1)/float(ny-1)*yl*yl


	! Input of stars coordinates and masses (x1,x2,mass)
	allow=100
	nstars=0
	do i=1,100000
	   read(1,*,end=1) rx1,rx2,re
	   if ((abs(rx1).le.xl1*allow).and.(abs(rx2).le.(xl2*allow))) then
	   nstars=nstars+1
	   xa1(nstars)=rx1
	   xa2(nstars)=rx2
	   eps(nstars)=re
	   end if
	end do
1	continue
	n=i-1


	! Computation of the mean mass and kappa of the stars
	rmasa=0.
	do k=1,nstars
	   rmasa=rmasa+eps(k)
	end do
	rmeanm=rmasa/nstars
	rk=rmasa*3.14159/xl/xl/4.
	rkcircle=2.*rk/3.14159
	n=nstars


	f1=1.-rks-rgs
	f2=1.-rks+rgs
	kconex=0

	! Lattice generation
	! First row
	krow=0
	do x1=-xl1,xl1,xs
	      krow=krow+1
	      y1=x1
	      x2=-xl2
	      y2=-xl2
	      y1=y1*f1
	      y2=y2*f2
	      do k=1,n
	         d2=(x1-xa1(k))*(x1-xa1(k))+(x2-xa2(k))*(x2-xa2(k))
           if (d2.gt.thres) then
             y1=y1-eps(k)*((x1-xa1(k))/(d2))
  	         y2=y2-eps(k)*((x2-xa2(k))/(d2))
          else
            y1t = 0
            y2t = 0
            call halodeflection(concentration,halonumber,x1,x2,k,xa1(k),xa2(k),xh1,xh2,y1t,y2t)
	    y1 = y1+y1t
	    y2 = y2+y2t
          end if
	      end do
	      rowax1(krow)=x1
	      rowax2(krow)=-xl2
	      roway1(krow)=(y1+yl)/ys+1
	      roway2(krow)=(y2+yl)/ys+1
	end do



	!First column
	kcol=0
	do x2=-xl2,xl2,xs
	      kcol=kcol+1
	      x1=-xl1
	      y1=-xl1
	      y2=x2
	      y1=y1*f1
	      y2=y2*f2
	      do k=1,n
	         d2=(x1-xa1(k))*(x1-xa1(k))+(x2-xa2(k))*(x2-xa2(k))
           if (d2.gt.thres) then
  	         y1=y1-eps(k)*((x1-xa1(k))/(d2))
  	         y2=y2-eps(k)*((x2-xa2(k))/(d2))
          else
            y1t = 0
            y2t = 0
            call halodeflection(concentration,halonumber,x1,x2,k,xa1(k),xa2(k),xh1,xh2,y1t,y2t)
	    y1 = y1+y1t
	    y2 = y2+y2t
          end if
	      end do
	      colax1(kcol)=-xl1
	      colax2(kcol)=x2
	      colay1(kcol)=(y1+yl)/ys+1
	      colay2(kcol)=(y2+yl)/ys+1
	end do

	ksup=0
	kiter=0
	kintro=0
	krow=1
	kshear=0
	ipercent=0
	niter=4.*xl1*xl2/xs/xs
	do x1=-xl1+xs,xl1,xs
	   krow=krow+1
	   vbx1=rowax1(krow)!-1)
	   vbx2=rowax2(krow)!-1)
	   vby1=roway1(krow)!-1)
	   vby2=roway2(krow)!-1)
	   kcol=1
	   do x2=-xl2+xs,xl2,xs
	      kiter=kiter+1
	      kcol=kcol+1
	      vax1=colax1(kcol-1)
	      vax2=colax2(kcol-1)
	      vay1=colay1(kcol-1)
	      vay2=colay2(kcol-1)
	      vdx1=colax1(kcol)
	      vdx2=colax2(kcol)
	      vdy1=colay1(kcol)
	      vdy2=colay2(kcol)
	      y1=x1
	      y2=x2
	      y1=y1*f1
	      y2=y2*f2
	      do k=1,n
	         d2=(x1-xa1(k))*(x1-xa1(k))+(x2-xa2(k))*(x2-xa2(k))
           if (d2.gt.thres) then
		counter1 = counter1 + 1
  	         y1=y1-eps(k)*((x1-xa1(k))/(d2))
  	         y2=y2-eps(k)*((x2-xa2(k))/(d2))
          else
		counter2 = counter2 + 1
            y1t = 0
            y2t = 0
            call halodeflection(concentration,halonumber,x1,x2,k,xa1(k),xa2(k),xh1,xh2,y1t,y2t)
	    y1 = y1+y1t
	    y2 = y2+y2t
          end if
	      end do
	      vcx1=x1
	      vcx2=x2
	      vcy1=(y1+yl)/ys+1
	      vcy2=(y2+yl)/ys+1
       	imin=min(vay1,vby1,vcy1,vdy1)+0.5
	jmin=min(vay2,vby2,vcy2,vdy2)+0.5
	imax=max(vay1,vby1,vcy1,vdy1)+0.5
	jmax=max(vay2,vby2,vcy2,vdy2)+0.5
	rnorm=0
	      if ((imin.gt.0).and.(jmin.gt.0).and.(imax.le.ns).and.(jmax.le.ns)) then
	do in=imin,imax
	  do jn=jmin,jmax
	     call sum02(in,jn,vay1,vay2,vby1,vby2,sab)
	     call sum02(in,jn,vby1,vby2,vcy1,vcy2,sbc)
	     call sum02(in,jn,vcy1,vcy2,vdy1,vdy2,scd)
	     call sum02(in,jn,vdy1,vdy2,vay1,vay2,sda)
	     b(in,jn)=sab+sbc+scd+sda
	     rnorm=rnorm+abs(b(in,jn))
	  end do
	end do

	rnormant=rnorm
	if (rnorm.le.0.0000001) rnorm=1.
	if (abs(rnorm).le.threshold) then
	   ic1=-1
      call corte(vay1,vay2,vby1,vby2,vcy1,vcy2,vdy1,vdy2,ic1,aby1,aby2)
           ic2=-1
      call corte(vay1,vay2,vdy1,vdy2,vcy1,vcy2,vby1,vby2,ic2,ady1,ady2)

	   if (ic1.eq.1) then
	   kconex=kconex+1
	   rnorm=0
	   do in=imin,imax
	     do jn=jmin,jmax
	       call sum02(in,jn,vay1,vay2,aby1,aby2,saab)
	       call sum02(in,jn,aby1,aby2,vdy1,vdy2,sabd)
	       call sum02(in,jn,vdy1,vdy2,vay1,vay2,sda)
	       call sum02(in,jn,aby1,aby2,vby1,vby2,sabb)
	       call sum02(in,jn,vby1,vby2,vcy1,vcy2,sbc)
	       call sum02(in,jn,vcy1,vcy2,aby1,aby2,scab)
	       b(in,jn)=abs(saab+sabd+sda)+abs(sabb+sbc+scab)
	       rnorm=rnorm+b(in,jn)
	     end do
	   end do
	   do in=imin,imax
	     do jn=jmin,jmax
	        a(in,jn)=a(in,jn)+abs(b(in,jn)/rnorm)
	     end do
	   end do
	   else if (ic2.eq.1) then
	   kconex=kconex+1
	   rnorm=0
	   do in=imin,imax
	     do jn=jmin,jmax
	       call sum02(in,jn,vay1,vay2,ady1,ady2,saad)
	       call sum02(in,jn,ady1,ady2,vby1,vby2,sadb)
	       call sum02(in,jn,vby1,vby2,vay1,vay2,sba)
	       call sum02(in,jn,ady1,ady2,vdy1,vdy2,sadd)
	       call sum02(in,jn,vdy1,vdy2,vcy1,vcy2,sdc)
	       call sum02(in,jn,vcy1,vcy2,ady1,ady2,scad)
	       b(in,jn)=abs(saad+sadb+sba)+abs(sadd+sdc+scad)
	       rnorm=rnorm+abs(b(in,jn))
	     end do
	   end do
	   do in=imin,imax
	     do jn=jmin,jmax
	        a(in,jn)=a(in,jn)+abs(b(in,jn)/rnorm)
	     end do
	   end do
	   else
	   do in=imin,imax
	     do jn=jmin,jmax
	        a(in,jn)=a(in,jn)+abs(b(in,jn)/rnorm)
	     end do
	   end do
	   end if
	else
	do in=imin,imax
	  do jn=jmin,jmax
	    a(in,jn)=a(in,jn)+abs(b(in,jn)/rnorm)
	  end do
	end do
	end if
	      end if
	      colax1(kcol-1)=vbx1
	      colax2(kcol-1)=vbx2
	      colay1(kcol-1)=vby1
	      colay2(kcol-1)=vby2
	      vbx1=vcx1
	      vbx2=vcx2
	      vby1=vcy1
	      vby2=vcy2
	  end do
	      colax1(kcol)=vcx1
	      colax2(kcol)=vcx2
	      colay1(kcol)=vcy1
	      colay2(kcol)=vcy2
	end do
	rmagni=0
	do i=1,ny
	   do j=1,ny
	      rmagni=rmagni+a(i,j)/raypix
	      write (3,*) a(i,j)/raypix
	   end do
	end do
	rmagni=rmagni/ny/ny
	rmagnitheor=1./((1-rktot-rgs)*(1-rktot+rgs))
	rmagnimass=1./((1-(rk+rks)-rgs)*(1-(rk+rks)+rgs))
	ratio = counter1/(counter1+counter2)
	print*, ratio
	end

	subroutine sum02(i,j,vay1,vay2,vby1,vby2,sab)
        implicit real*8 (A-H,O-Z)
	sgab=(vby1-vay1)/abs(vby1-vay1)
	if (sgab.ge.0) then
	   xi=max(vay1,i-0.5)
	   xs=min(vby1,i+0.5)
	else
	   xi=min(vay1,i+0.5)
	   xs=max(vby1,i-0.5)
	end if
	if (sgab*(xs-xi).le.0) then
	   sab=0
	   return
	end if
	a=(vby2-vay2)/(vby1-vay1)
	b=vay2-a*vay1
	fxi=b+a*xi
	fxs=b+a*xs
	xmin=(j-0.5-b)/a
	xplu=(j+0.5-b)/a
	if (fxi.ge.(j+0.5)) then
	   if (fxs.ge.(j+0.5)) then
	      sab=xs-xi
	      return
	   else if (fxs.ge.(j-0.5)) then
	      sab=xplu-xi+0.5*a*(xs*xs-xplu*xplu)
	      sab=sab+(b-j+0.5)*(xs-xplu)
	      return
	   else
	      sab=xplu-xi+0.5*a*(xmin*xmin-xplu*xplu)
	      sab=sab+(b-j+0.5)*(xmin-xplu)
	      return
	   end if
	else if (fxi.ge.(j-0.5)) then
	   if (fxs.ge.(j+0.5)) then
	      sab=xs-xplu+0.5*a*(xplu*xplu-xi*xi)
	      sab=sab+(b-j+0.5)*(xplu-xi)
	      return
	   else if (fxs.ge.(j-0.5)) then
	      sab=0.5*a*(xs*xs-xi*xi)
	      sab=sab+(b-j+0.5)*(xs-xi)
	      return
	   else
	      sab=0.5*a*(xmin*xmin-xi*xi)
	      sab=sab+(b-j+0.5)*(xmin-xi)
	      return
	   end if
	else
	   if (fxs.ge.(j+0.5)) then
	      sab=xs-xplu+0.5*a*(xplu*xplu-xmin*xmin)
	      sab=sab+(b-j+0.5)*(xplu-xmin)
	      return
	   else if (fxs.ge.(j-0.5)) then
	      sab=0.5*a*(xs*xs-xmin*xmin)
	      sab=sab+(b-j+0.5)*(xs-xmin)
	      return
	   else
	      sab=0
	      return
	   end if
	end if
	print*,"ERROR---------------------------"
	return
	end

	subroutine sum03(i,j,va1,va2,vb1,vb2,sab)
        implicit real*8 (A-H,O-Z)
	if (vb1.lt.va1) then
	   vay1=vb1
	   vay2=vb2
	   vby1=va1
	   vby2=va2
	else
	   vay1=va1
	   vay2=va2
	   vby1=vb1
	   vby2=vb2
	end if
	a=(vby2-vay2)/(vby1-vay1)
	b=vay2-a*vay1
	if (a.gt.0) then
	   xmin=(j-0.5-b)/a
	   xplu=(j+0.5-b)/a
	   xi=min(i+0.5,max(i-0.5,xmin,vay1))
	   xf=min(i+0.5,max(i-0.5,xmin,vby1))
	   xl=max(xi,min(xf,xplu))
	   sab=0.5*a*(xl*xl-xi*xi)+(b-j+0.5)*(xl-xi)
	   sab=sab+xf-xl
	   if (vb1.lt.va1) sab=-sab
	   return
 	else if (a.lt.0) then
	   xmin=(j-0.5-b)/a
	   xplu=(j+0.5-b)/a
	   xi=max(i-0.5,min(i+0.5,xmin,vay1))
	   xf=max(i-0.5,min(i+0.5,xmin,vby1))
	   xl=max(xi,min(xf,xplu))
	   sab=0.5*a*(xf*xf-xl*xl)+(b-j+0.5)*(xf-xl)
	   sab=sab+xl-xi
	   if (vb1.lt.va1) sab=-sab
	   return
	else
	   xi=min(i+0.5,max(i-0.5,vay1))
	   xf=min(i+0.5,max(i-0.5,vby1))
	   fun=min(b,j+0.5)
	   if (b.lt.(j-0.5)) then
	      sab=0
	   else
	      sab=(fun-j+0.5)*(xf-xi)
	   end if
	   if (vb1.lt.va1) sab=-sab
	   return
	end if
	end

	subroutine norm02(vay1,vay2,vby1,vby2,sab)
        implicit real*8 (A-H,O-Z)
	a=(vby2-vay2)/(vby1-vay1)
	b=vay2-a*vay1
	sab=0.5*a*(vby1*vby1-vay1*vay1)+b*(vby1-vay1)
	return
	end

	subroutine corte(vay1,vay2,vby1,vby2,vcy1,vcy2,vdy1,vdy2,ic,xc,yc)
        implicit real*8 (A-H,O-Z)
	a=(vby2-vay2)/(vby1-vay1)
	b=vay2-a*vay1
	c=(vdy2-vcy2)/(vdy1-vcy1)
	d=vcy2-c*vcy1
	xc=(b-d)/(c-a)
	yc=a*xc+b
	rl=(xc-vay1)/(vby1-vay1)
	rm=(xc-vcy1)/(vdy1-vcy1)
	if ((rl.gt.0).and.(rl.lt.1).and.(rm.gt.0).and.(rm.lt.1)) ic=1
	return
	end

  subroutine halodeflection(concentration,halonumber,x1,x2,k,xt1,xt2,xh1,xh2,y1r,y2r)
        implicit real*8 (A-H,O-Z)
    real*8 xh1(int(halonumber),int(concentration))
    real*8 xh2(int(halonumber),int(concentration))
    integer::k2
    k2 = mod(k,int(halonumber))+1

    bhmass = 1/concentration
    y1r=0
    y2r=0
    conc2 = int(concentration)
    do j2=1,conc2
      dt = (x1-xt1-xh1(k2,j2))*(x1-xt1-xh1(k2,j2))+(x2-xt2-xh2(k2,j2))*(x2-xt2-xh2(k2,j2))
      y1r=y1r-bhmass*((x1-xt1-xh1(k2,j2))/(dt))
      y2r=y2r-bhmass*((x2-xt2-xh2(k2,j2))/(dt))
    end do
    return
    end
