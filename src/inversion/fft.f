ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine trans_3dto2d(seis,nt,dt,fmin,fmax)
	parameter(maxsample=4096)
	parameter(pi=3.141592654)
	dimension seis(*)
	complex sf(4*maxsample),xtmp
	nnn=1
	do i=1,20
	   if(nnn.lt.2*nt)then
	      nnn=nnn*2
	   endif
	enddo
	do i=1,nt
	   sf(i)=cmplx(seis(i),0.0)
	enddo
	do i=nt+1,nnn
	   sf(i)=cmplx(0.0,0.0)
	enddo
	call fft(sf,nnn,1)
	xtmp=cmplx(0.5,0.5)
	ii1=nint(fmin*nnn*dt)+1
	ii2=nint(fmax*nnn*dt)+1
	ii2=min(ii2,nnn/2+1)
	do i=1,ii1-1
	   sf(i)=cmplx(0.0,0.0)
	enddo
	do i=ii1,ii2
	   w1=sqrt((i-1)/real(nnn*dt)*2.0*pi)
	   sf(i)=sf(i)*xtmp/w1
	enddo
	do i=ii2+1,nnn/2+1
	   sf(i)=cmplx(0.0,0.0)
	enddo
	ii3=nnn-ii2+2
	ii4=nnn-ii1+2
	do i=nnn/2+2,ii3-1
	   sf(i)=cmplx(0.0,0.0)
        enddo
	do i=ii3,ii4
	   w1=sqrt((nnn+1-i)/real(nnn*dt)*2.0*pi)
           sf(i)=sf(i)*xtmp/w1
        enddo
	do i=ii4+1,nnn
	   sf(i)=cmplx(0.0,0.0)
        enddo
	call fft(sf,nnn,-1)
	do i=1,nt
	   tmp=float(i-1)*dt
	   tmp=sqrt(tmp) 
	   seis(i)=real(sf(i))*tmp
	enddo
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine lowpass(trace,nt,dt,fmax)
	parameter(maxsample=16384,ntaper_max=1000)
	real trace(*)
	real taper1(ntaper_max),taper2(ntaper_max)
	complex x(4*maxsample)
c
	pi=acos(-1.0)
c
	nfft=1
	do i=1,20
	   if(nfft.lt.2*nt)then
	      nfft=nfft*2
	   endif
	enddo
	dfreq=1.0/float(nfft)/dt
        f2=fmax
	if2=1+nint(f2/dfreq)
	if2=min(if2,1+nfft/2)
	if3=nfft-if2+2
	if4=nfft
	ilength=if2
	ratio=0.05
	call get_taper(taper1,taper2,ntaper,nt,ratio)
c
	do it=1,ntaper
	   x(it)=cmplx(trace(it)*taper1(it),0.0)
	enddo
	do it=ntaper+1,nt-ntaper
	   x(it)=cmplx(trace(it),0.0)
	enddo
	do it=nt-ntaper+1,nt
	   x(it)=cmplx(trace(it)*taper2(it-nt+ntaper),0.0)
	enddo
	do it=nt+1,nfft
	   x(it)=cmplx(0.0,0.0)
	enddo
	call fft_new(x,nfft,1)
	ratio=0.05
	call get_taper(taper1,taper2,ntaper,nt,ratio)
c	write(*,*)ntaper*dfreq
c	do it=if2-ntaper,if3-1
c	   x(it)=cmplx(0.0,0.0)
c        enddo
c	   do it=if2-ntaper+1,if2
c	      x(it)=x(it)*taper2(it-if2+1)
c	   enddo
c	   do it=if3,if3+ntaper-1
c	      x(it)=x(it)*taper1(it-if3+ntaper)
c	   enddo
	do it=if2+ntaper,if3-ntaper
	   x(it)=cmplx(0.0,0.0)
        enddo
	   do it=if2,if2+ntaper-1
	      x(it)=x(it)*taper2(it-if2+1)
	   enddo
	   do it=if3-ntaper+1,if3
	      x(it)=x(it)*taper1(it-if3+ntaper)
	   enddo
	call fft_new(x,nfft,-1)
	do it=1,nt
	   trace(it)=real(x(it))
	enddo
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	subroutine bandpass(trace,nt,dt,f1,f2)
	parameter(maxsample=16384,ntaper_max=1000)
	real trace(*)
	real taper1(ntaper_max),taper2(ntaper_max)
	complex x(4*maxsample)
c
	pi=acos(-1.0)
c
	nfft=1
	do i=1,20
	   if(nfft.lt.2*nt)then
	      nfft=nfft*2
	   endif
	enddo
	dfreq=1.0/float(nfft)/dt
	if(f1.gt.f2)then
	   tmp=f2
	   f2=f1
	   f1=tmp
	endif
	if1=1+nint(f1/dfreq)
	if1=max(if1,1)
	if2=1+nint(f2/dfreq)
	if2=min(if2,1+nfft/2)
	if3=nfft-if2+2
	if4=nfft-if1+2
	ilength=if2-if1+1
	ratio=0.05
	call get_taper(taper1,taper2,ntaper,nt,ratio)
c
	do it=1,ntaper
	   x(it)=cmplx(trace(it)*taper1(it),0.0)
	enddo
	do it=ntaper+1,nt-ntaper
	   x(it)=cmplx(trace(it),0.0)
	enddo
	do it=nt-ntaper+1,nt
	   x(it)=cmplx(trace(it)*taper2(it-nt+ntaper),0.0)
	enddo
	do it=nt+1,nfft
	   x(it)=cmplx(0.0,0.0)
	enddo
	call fft_new(x,nfft,1)
c	call get_taper(taper1,taper2,ntaper,ilength,ratio)
	do it=1,if1-ntaper
	   x(it)=cmplx(0.0,0.0)
	enddo
	do it=if2+ntaper,if3-ntaper
	   x(it)=cmplx(0.0,0.0)
        enddo
	do it=if4+ntaper,nfft
           x(it)=cmplx(0.0,0.0)
        enddo
	   do it=if1-ntaper+1,if1
	      x(it)=x(it)*taper1(it-if1+ntaper)
	   enddo
	   do it=if2,if2+ntaper-1
	      x(it)=x(it)*taper2(it-if2+1)
	   enddo
	   do it=if3-ntaper+1,if3
	      x(it)=x(it)*taper1(it-if3+ntaper)
	   enddo
	   do it=if4,if4+ntaper-1
	      x(it)=x(it)*taper2(it-if4+1)
	   enddo
	call fft_new(x,nfft,-1)
	do it=1,nt
	   trace(it)=real(x(it))
	enddo
	return
	end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine get_taper(taper1,taper2,ntaper,nt,ratio)
        real taper1(*),taper2(*)
        pi=acos(-1.0)
        ntaper=int(float(nt)*ratio)
        do i=1,ntaper
           taper1(i)=1.0-cos(0.5*pi*float(i-1)/float(ntaper))**2
           taper2(ntaper-i+1)=taper1(i)
        enddo
        return
        end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine get_taper2(taper,nt,ratio)
        real taper(*)
        pi=acos(-1.0)
        ntaper=int(float(nt)*ratio)
        do i=1,ntaper
           taper(i)=1.0-cos(0.5*pi*float(i-1)/float(ntaper))**2
           taper(nt-i+1)=taper(i)
        enddo
	do i=ntaper+1,nt-ntaper
	   taper(i)=1.0
	enddo
        return
        end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine bandpass2(trace,nt,dt,f1,f2)
	parameter(maxsample=4096,ntaper=3)
	real trace(*)
	real taper(20)
	complex x(4*maxsample)
c
	pi=acos(-1.0)
	do i=1,ntaper
	   ss=float(i-1)*pi/2.0/float(ntaper)
           taper(i)=1.0-cos(ss)*cos(ss)
	enddo
c
	nfft=1
	do i=1,20
	   if(nfft.lt.2*nt)then
	      nfft=nfft*2
	   endif
	enddo
	dfreq=1.0/float(nfft)/dt
	if(f1.gt.f2)then
	   tmp=f2
	   f2=f1
	   f1=tmp
	endif
	if1=1+nint(f1/dfreq)
	if1=max(if1,1)
	if2=1+nint(f2/dfreq)
	if2=min(if2,1+nfft/2)
	if3=nfft-if2+1
	if4=nfft-if1+1
c
c	do it=1,ntaper
c	   x(it)=cmplx(trace(it)*taper(it),0.0)
c	enddo
c	do it=ntaper+1,nt-ntaper
c	   x(it)=cmplx(trace(it),0.0)
c	enddo
c	do it=nt-ntaper+1,nt
c	   x(it)=cmplx(trace(it)*taper(nt-it+1),0.0)
c	enddo
c	do it=nt+1,nfft
c	   x(it)=cmplx(0.0,0.0)
c	enddo
	do it=1,nt
	   x(it)=cmplx(trace(it),0.0)
	enddo
	do it=nt+1,nfft
	   x(it)=cmplx(0.0,0.0)
       	enddo
	call fft(x,nfft,1)
c	do it=1,nfft
c	   write(*,*)it,abs(x(it))
c	enddo
c	pause
	do it=1,if1-ntaper-1
	   x(it)=cmplx(0.0,0.0)
	enddo
	do it=max(1,if1-ntaper),if1-1
	   x(it)=x(it)*taper(it-if1+ntaper+1)
	enddo
	do it=if2+1,min(if2+ntaper,nfft/2)
	   x(it)=x(it)*taper(if2+ntaper+1-it)
	enddo
	do it=if2+ntaper+1,nfft/2
	   x(it)=cmplx(0.0,0.0)
	enddo
	do it=nfft/2+1,if3-ntaper-1
	   x(it)=cmplx(0.0,0.0)
	enddo
	do it=max(if3-ntaper,nfft/2+1),if3-1
	   x(it)=x(it)*taper(it-if3+ntaper+1)
	enddo
	do it=if4+1,min(if4+ntaper,nfft)
	   x(it)=x(it)*taper(if4+ntaper+1-it)
	enddo
	do it=if4+ntaper+1,nfft
	   x(it)=cmplx(0.0,0.0)
	enddo
c	call fft_time(x,nfft,-1)
	call fft(x,nfft,-1)
	do it=1,nt
	   trace(it)=real(x(it))
	enddo
	return
	end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_bandpass(filter,xf1,xf2,xf3,xf4,nnn_kkk,dt)
	real filter(*)
        pi=acos(-1.0)
        df=1.0/float(nnn_kkk)/dt
        i1=NINT(xf1/df)+1
        i2=NINT(xf2/df)+1
        i3=NINT(xf3/df)+1
        i4=NINT(xf4/df)+1
        do i=1,i1
           filter(i)=0.0
        enddo
        do i=nnn_kkk,nnn_kkk-i1+2,-1
           filter(i)=0.0
        enddo
        do i=i4,nnn_kkk-i4+2
           filter(i)=0.0
        enddo
        do i=i2,i3
	   filter(i)=1.0
	enddo
	do i=nnn_kkk-i2+2,nnn_kkk-i3+2,-1
	   filter(i)=1.0
	enddo
 	c=pi*0.5/float(i2-i1)
        do i=i1+1,i2-1
	   s=sin(c*float(i-i1))
           filter(i)=s*s
	enddo
	do i=nnn_kkk-i1+1,nnn_kkk-i2+3,-1
           filter(i)=filter(nnn_kkk-i+2)
	enddo
	c=pi*0.5/float(i4-i3)
	do i=i3+1,i4-1
	   s=sin(c*float(i4-i))
	   filter(i)=s*s
	enddo
	do i=nnn_kkk-i3+1,nnn_kkk-i4+3,-1
	   filter(i)=filter(nnn_kkk-i+2)
	enddo
	return
	end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fft_new(x,kk,key)
c
        complex x(*),t
        real*8 pi,s,a01,a02,a03,b01,b02,b03
        pi=acos(-1.0)
c
        j=1
        s=-pi*key
        kk2=kk/2
        do 14 i=1,kk
        if (i.gt.j) goto 9
        t=x(j)
        x(j)=x(i)
        x(i)=t
9       m=kk2
10      if (j.le.m) goto 14
        j=j-m
        m=m/2
        if(m.ge.1) goto 10
14      j=j+m
        l=1
16      l1=l*2
        a03=s/l
        a01=1.0D0
        a02=0.0D0
        b01=dcos(a03)
        b02=dsin(a03)
        do 23 m=1,l
        do 22 i=m,kk,l1
        ip=i+l
        t=cmplx(a01,a02)*x(ip)
        x(ip)=x(i)-t
22      x(i)=x(i)+t
        b03=a01
        a01=a01*b01+a02*b02
23      a02=-b03*b02+b01*a02
        l=l1
        if(l.lt.kk)go to 16
        if(key.eq.-1)then
          a01=1.D0/real(kk)
          do i=1,kk
            x(i)=x(i)*a01
          enddo
        endif
c
        return
        end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   ******************
c                   * Subroutine FFT *
c                   ******************
c
        subroutine fft(x,kk,key)
c
        complex x(*),a,t,w
	double precision pi
	pi=dacos(-1.0d0)
c
        j=1
        s=sqrt(real(kk)**(key-1))
        do 14 i=1,kk
        if (i.gt.j) goto 9
        t=s*x(j)
        x(j)=s*x(i)
        x(i)=t
9       m=kk/2
10      if (j.le.m) goto 14
        j=j-m
        m=m/2
        if(m.ge.1) goto 10
14      j=j+m
        l=1
16      d=2*l
        do 23 m=1,l
        a=-(0.0,1.0)*pi*key*(m-1)/l
        w=cexp(a)
        do 23 i=m,kk,d
        t=w*x(i+l)
        x(i+l)=x(i)-t
23      x(i)=x(i)+t
        l=d
        if(l.lt.kk) goto 16
c
        return
        end
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c                   ******************
c                   * Subroutine FFT *
c                   ******************
c
        subroutine fft_time(x,kk,key)
c
        complex x(*),a,t,w
	double precision pikey,tmp
c
	do i=kk/2+2,kk
	   x(i)=x(kk+2-i)
	enddo
c
	pikey=-dacos(-1.0d0)*key
c
        j=1
        s=sqrt(real(kk)**(key-1))
        do 14 i=1,kk
        if (i.gt.j) goto 9
        t=s*x(j)
        x(j)=s*x(i)
        x(i)=t
9       m=kk/2
10      if (j.le.m) goto 14
        j=j-m
        m=m/2
        if(m.ge.1) goto 10
14      j=j+m
        l=1
16      d=2*l
	tmp=pikey/l
        do 23 m=1,l
        a=cmplx(0.0,tmp*(m-1))
        w=cexp(a)
        do 23 i=m,kk,d
        t=w*x(i+l)
        x(i+l)=x(i)-t
23      x(i)=x(i)+t
        l=d
        if(l.lt.kk) goto 16
c
	do i=2,kk
	   x(i)=x(i)*2.0
	enddo
        return
        end
