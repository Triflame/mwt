ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_input_parameter

      implicit none
      include 'wtw_model.para'
      character*200 line(100),line_tmp
      integer i, k, ind, is_suc, is_sou, is_den
      real x(4)
c     --------------------------------------------------------
      open(11,file='WTW_INPUT.DAT',status='old')
      do i=1,1000
        read(11,'(a)',end=111)line_tmp
        ind=index(line_tmp,'#')
        if(ind.ne.0)then
          do k=ind,200
            line_tmp(k:k)=' ' 
          enddo
        endif
        line(i)=line_tmp
      enddo
111   continue
      close(11)
      i=i-1

c  get input file names
      call get_str_new(coord_temp,'COORD_TEMP',line,i,is_suc)
      if(is_suc.eq.0)then
        write(*,*)'Please add input item COORD_TEMP.'
        stop
      endif
      call get_str_new(vel_in,'VEL_IN',line,i,is_suc)
      if(is_suc.eq.0)then
        write(*,*)'Please add input item VEL_IN.'
        stop
      endif
      call get_num_new(x,1,'IS_SU_VELOCITY',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(28)=0
      else
        data_int(28)=int(x(1))
      endif
      call get_str_new(csg_pre,'CSG_PRE_IN',line,i,is_suc)
      if(is_suc.eq.0)then
        write(*,*)'Please add input item CSG_PRE_IN.'
        stop
      endif
      call get_str_new(csg_suffix,'CSG_SUFFIX_IN',line,i,is_suc)
      call get_str_new(source_in,'SOURCE_IN',line,i,is_sou)
      call get_str_new(density_in,'DENSITY_IN',line,i,is_den)
      call get_num_new(x,1,'IS_SU_DENSITY',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(29)=0
      else
        data_int(29)=int(x(1))
      endif
c get output file names
      call get_str_new(wtw_ou,'WTW_OU',line,i,is_suc)
      if(is_suc.eq.0)wtw_ou='WTW_OU'
      call get_str_new(res_wave,'RES_WAVEFORM',line,i,is_suc)
      if(is_suc.eq.0)res_wave='RES_WAVE'
      call get_str_new(syn_ou_pre,'SYN_OU_PRE',line,i,is_suc)
      if(is_suc.eq.0)syn_ou_pre='SYN_OU'
      call get_str_new(syn_ou_suffix,'SYN_OU_SUFFIX',line,i,is_suc)
      if(is_suc.eq.0)syn_ou_suffix=''
      call get_str_new(vel_ou_pre,'VEL_OU_PRE',line,i,is_suc)
      if(is_suc.eq.0)vel_ou_pre='VEL_OU'
      call get_str_new(vel_ou_suffix,'VEL_OU_SUFFIX',line,i,is_suc)
      if(is_suc.eq.0)vel_ou_suffix=''
      call get_str_new(grad_ou,'GRAD_OU',line,i,is_suc)
      if(is_suc.eq.0)grad_ou='GRAD_OU'

c model parameters
      call get_num_new(x,3,'Nx_IN,Nz_IN,Dx_IN',line,i,is_suc)
      if(is_suc.eq.0)then
        write(*,*)'Please input Nx_IN,Nz_IN,Dx_IN.'
        stop
      endif
      data_int(1)=int(x(1))
      data_int(2)=int(x(2))
      data_real(1)=x(3)
      if(data_int(1).gt.nx_mod0)then
        write(*,*)'nx_in(',data_int(1),')>nx_mod0(',nx_mod0,')'
        stop
      endif
      if(data_int(2).gt.nz_mod0)then
        write(*,*)'nz_in(',data_int(2),')>nz_mod0(',nz_mod0,')'
        stop
      endif
      call get_num_new(x,2,'XMIN_IN,ZMIN_IN',line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(25)=0.0
        data_real(26)=0.0
      else
         data_real(25)=x(1)
         data_real(26)=x(2)
      endif
      call get_num_new(x,1,'XMAX',line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(19)=data_real(25)+float(data_int(1)-1)*data_real(1)
      else
        data_real(19)=x(1)
      endif
      call get_num_new(x,1,'ZMAX',line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(21)=data_real(26)+float(data_int(2)-1)*data_real(1)
      else
        data_real(21)=x(1)
      endif

c data parameters
      call get_num_new(x,3,'Nt_IN,Dt_IN,IS_SU_CSG',line,i,is_suc)
      if(is_suc.eq.0)then
        write(*,*)'Please input Nt_IN,Dt_IN,IS_SU_CSG.'
        stop
      endif
      data_int(3)=int(x(1))
      data_real(2)=x(2)
      data_int(4)=int(x(3))
      if(data_int(3).gt.nt0)then
        write(*,*)'nt_in(',data_int(3),')>nt0(',nt0,')'
        stop
      endif
      call get_num_new(x,2,'Dx_OUT,Dt_OUT',line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(3)=data_real(1)
        data_real(4)=data_real(2)
      else
        data_real(3)=x(1)
        data_real(4)=x(2)
      endif
      call get_num_new(x,1,'ISMARINE',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(5)=0
      else
        data_int(5)=int(x(1))
      endif
      call get_num_new(x,2,'IS_RESET_SURFACE,IS_RESET_SG',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(6)=0
        data_int(7)=0
      else
        data_int(6)=int(x(1))
        data_int(7)=int(x(2))
      endif
      call get_num_new(x,2,'DEPSOU,DEPREC',line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(5)=0.0
        data_real(6)=0.0
      else
        data_real(5)=x(1)
        data_real(6)=x(2)
      endif

c source parameters
      call get_num_new(x,1,'NSHOT_INT',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(27)=1
      else
        data_int(27)=int(x(1))
      endif
      call get_num_new(x,1,'NTPAD',line,i,is_suc)
      if(is_suc.eq.0)then
           data_int(26)=0
        else
           data_int(26)=int(x(1))
        endif
      call get_num_new(x,4,'NOPSOU,FREQ,NTPSOU,DTSOU',line,i,is_suc)
      if(is_suc.eq.0)then
        write(*,*)'Please input NOPSOU,FREQ,NTPSOU,DTSOU.'
        stop
      endif
      data_int(8)=int(x(1))
      data_real(7)=x(2)
      data_int(9)=int(x(3))
      data_real(24)=x(4)
      if(int(x(1)).ne.0.and.is_sou.eq.0)then
        write(*,*)'Please input SOURCE_IN or change NOPSOU to 1 or 2'
        stop
      endif

c density parameters
      call get_num_new(x,4,'NDENS,DENCO,CP00,DENP',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(10)=1
c        data_real(8)=1.0
c        data_real(9)=1.0
c        data_real(10)=1.0
        data_real(8)=310.0
        data_real(9)=1.0
        data_real(10)=0.25
      else
        data_int(10)=int(x(1))
        data_real(8)=x(2)
        data_real(9)=x(3)
        data_real(10)=x(4)
      endif
      if(data_int(10).eq.0.and.is_den.eq.0)then
        write(*,*)'Please input DENSITY_IN or change NDENS to 1 or 2'
        stop
      endif

c free surface
      call get_num_new(x,1,'NFREE',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(11)=0
      else
        data_int(11)=int(x(1))
      endif

c control parameters
      call get_num_new(x,2,'ITERATION,ITER00',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(12)=200
        data_int(13)=1
      else
        data_int(12)=int(x(1))
        data_int(13)=int(x(2))
      endif
      call get_num_new(x,1,'NITE',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(14)=20
      else
        data_int(14)=int(x(1))
      endif
      call get_num_new(x,1,'IS_USE_CON1',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(15)=0
      else
        data_int(15)=int(x(1))
      endif
      call get_num_new(x,3,'DEPTH_CON1,VMIN_CON1,VMAX_CON1',
     1           line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(11)=0.0
        data_real(12)=1500.0
        data_real(13)=1500.0
      else
        data_real(11)=x(1)
        data_real(12)=x(2)
        data_real(13)=x(3)
      endif
      call get_num_new(x,1,'IS_USE_WATER_BOTTOM',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(25)=0
      else
        data_int(25)=int(x(1))
      endif
      if( data_int(25).eq.1)then
        call get_str_new(waterBottomFile,'WATER_BOTTOM_FILE',line,
     1            i,is_suc)
        if(is_suc.eq.0)then
          write(*,*)'Please input DENSITY_IN or change NDENS to 1 or 2'
          stop
        endif
      endif
    
      call get_num_new(x,4,'DVMIN,DVMAX,VMIN,VMAX',line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(14)=6.0
        data_real(15)=30.0
        data_real(16)=700.0
        data_real(17)=7000.0
      else
        data_real(14)=x(1)
        data_real(15)=x(2)
        data_real(16)=x(3)
        data_real(17)=x(4)
      endif
      call get_num_new(x,2,'XMIN,ZMIN',line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(18)=0.0
        data_real(20)=0.0
      else
        data_real(18)=x(1)
        data_real(20)=x(2)
      endif
      call get_num_new(x,2,'OFFMIN,OFFMAX',line,i,is_suc)
      if(is_suc.eq.0)then
        data_real(22)=0.0
        data_real(23)=-1.0
      else
        data_real(22)=x(1)
        data_real(23)=x(2)
      endif
      call get_num_new(x,3,'ISPRE,ISNORM,ISUSECONJUGATE',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(16)=1
        data_int(17)=21
        data_int(18)=1
      else
        data_int(16)=int(x(1))
        data_int(17)=int(x(2))
        data_int(18)=int(x(3))
      endif
      call get_num_new(x,2,'ISPRECONDITION,ISEPOW',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(19)=1
        data_int(24)=1
      else
        data_int(19)=int(x(1))
        data_int(24)=int(x(2))
      endif
      call get_num_new(x,1,'INDEX_SHOT_OUTPUT',line,i,is_suc)
      if(is_suc.eq.0)then
        data_int(20)=1
      else
         data_int(20)=int(x(1))
      endif

c Nx_OUT, Nz_OUT
      data_int(21)=nint((data_real(19)-data_real(18))/data_real(3))+1
      data_int(22)=nint((data_real(21)-data_real(20))/data_real(3))+1
c      data_int(21)=nint(float(data_int(1)-1)*data_real(1)
c     1          /data_real(3))+1
c      data_int(22)=nint(float(data_int(2)-1)*data_real(1)
c     1          /data_real(3))+1
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_geometry

      implicit none
      include 'wtw_model.para'
      real xwater(nx_mod),zwater(nx_mod)
      real xwater2(nx_mod),zwater2(nx_mod)
      real x000, z000, x1, z1, x2, z2, tt, tmp1, tmp2, tmp3, xx01,
     1  xx02, xx, zz, zerodis, xcor, zcor
      integer ind(nx_mod), ig, i, ii, jj, isin, j, ix, iz, kk, nwp

      real sqrt, min, max
      intrinsic sqrt, min, max
c     --------------------------------------------------------
      do is=1,ns0
         ngis(is)=0
         nxss(is)=0
         nzss(is)=0
         xss(is)=0.0
         zss(is)=0.0
         do ig=1,ng0
            nxg(ig,is)=0
            nzg(ig,is)=0
            xgg(ig,is)=0.0
            zgg(ig,is)=0.0
            tr(ig,is)=-1.0
            offset(ig,is)=-1.0
         enddo
      enddo
      ns=0
      xmin=data_real(18)
      xmax=data_real(19)
      zmin=data_real(20)
      zmax=data_real(21)
      x000=xmin
      z000=zmin
      offmin=data_real(22)
      offmax=data_real(23)
      dx=data_real(3)

      write(*,*) xmin,xmax,zmin,zmax,offmin,offmax

      open(unit=31,file=coord_temp)
      rewind(31)

      do i=1,ns0*ng0
         read(31,*,end=101)ii,jj,x1,z1,x2,z2,tt
         tmp1=x1-x2
         tmp2=z1-z2
         tmp3=sqrt(tmp1*tmp1+tmp2*tmp2)
c         if((x1.ge.xmin.and.x1.le.xmax).and.(x2.ge.xmin.and.x2.le.xmax)
c     1        .and.((offmax.lt.offmin).or.
c     2        (tmp3.ge.offmin.and.tmp3.le.offmax)))then
         if((x1.ge.xmin.and.x1.le.xmax).and.(x2.ge.xmin.and.x2.le.xmax)
     >      )then
            isin=0
            do j=1,ns
               if(ii.eq.indshot(j))then
                  isin=j
                  go to 102
               endif
            enddo
102         continue
            if(isin.ne.0)then
               ngis(isin)=max(ngis(isin),jj)
               xgg(jj,isin)=x2-xmin
               zgg(jj,isin)=z2-zmin
               tr(jj,isin)=tt
               offset(jj,isin)=tmp3
            else
               ns=ns+1
               indshot(ns)=ii
               xss(ns)=x1-xmin
               zss(ns)=z1-zmin
               ngis(ns)=max(ngis(ns),jj)
               xgg(jj,ns)=x2-xmin
               zgg(jj,ns)=z2-zmin
               tr(jj,ns)=tt
               offset(jj,ns)=tmp3
            endif
         endif
      enddo
101   continue
      close(31)
      data_int(23)=ns

      do is=1,ns
         xx01=xss(is)
         xx02=xss(is)
         do ig=1,ngis(is)
            if(offset(ig,is).gt.-0.5.and.tr(ig,is).ge.0.0)then
               xx01=min(xx01,xgg(ig,is))
               xx02=max(xx02,xgg(ig,is))
            endif
         enddo
         ix1_shot(is)=nint(xx01/dx)+1-npad_shot
         ix2_shot(is)=nint(xx02/dx)+1+npad_shot
c         if(is.eq.241)then
c           write(*,*)is,ix1_shot(is),ix2_shot(is),nx_shot
c         endif
      enddo
      if(offmax.lt.offmin)then
         offmin=1.0e+10
         offmax=0.0
         do is=1,ns
            do ig=1,ngis(is)
               if(offset(ig,is).gt.-0.5.and.tr(ig,is).ge.0.0)then
                  offmin=min(offmin,offset(ig,is))
                  offmax=max(offmax,offset(ig,is))
               endif
            enddo
         enddo
         data_real(22)=offmin
         data_real(23)=offmax
      endif
      xmax=xmax-xmin
      zmax=zmax-zmin
      xmin=0.0
      zmin=0.0
      data_real(18)=xmin
      data_real(19)=xmax
      data_real(20)=zmin
      data_real(21)=zmax
      call get_surf

c      do ix=1,nx_tot
c         nsurf_tot(ix)=1
c      enddo

      is_marine=data_int(5)
      is_reset_surface=data_int(6)
      is_reset_sg=data_int(7)
      depsou=data_real(5)
      deprec=data_real(6)
c      if(is_marine.eq.1)then
c	 if(is_reset_surface.eq.1)then
c	    do ix=1,nx_tot
c	       nsurf_tot(ix)=1
c	    enddo
c	 endif
c	 if(is_reset_sg.eq.1)then
c	    do is=1,ns
c	       zss(is)=zss(is)+depsou
c	       do ig=1,ngis(is)
c	          zgg(ig,is)=zgg(ig,is)+deprec
c	       enddo
c	    enddo
c	 endif
c      endif
c
      do is=1,ns
	 ix=nint(xss(is)/dx)+1
	 iz=nint(zss(is)/dx)+1
	 if(iz.le.nsurf_tot(ix))then
	    iz=nsurf_tot(ix)+1
	 endif
         nxss(is)=ix-ix1_shot(is)+na_new+1
	 nzss(is)=iz+na_new
	 do ig=1,ngis(is)
            ix=nint(xgg(ig,is)/dx)+1
	    iz=nint(zgg(ig,is)/dx)+1
	    if(iz.le.nsurf_tot(ix))then
               iz=nsurf_tot(ix)+1
            endif
	    nxg(ig,is)=ix-ix1_shot(is)+na_new+1
	    nzg(ig,is)=iz+na_new
	 enddo
      enddo
c
      if( data_int(25).eq.1 )then
	 open(40,file=waterBottomFile)
	 rewind(40)
	 do i=1,10000
	    read(40,*,err=601,end=601)xx,zz
	    xwater(i)=xx-x000
	    zwater(i)=zz-z000
	    nwp=i
         enddo
601	 continue
	 close(40)
	 call sort(nwp,xwater,ind)
	 kk=1
	 xwater2(kk)=xwater(ind(1))
	 zwater2(kk)=zwater(ind(1))
         zerodis=0.01*dx
	 do i=2,nwp
	    if(abs(xwater(ind(i))-xwater2(kk)).gt.zerodis)then
	       kk=kk+1
	       xwater2(kk)=xwater(ind(i))
	       zwater2(kk)=zwater(ind(i))
	    endif
	 enddo
         nx_tot=data_int(21)
	 do i=1,nx_tot
	    xcor=float(i-1)*dx
	    call get1D(zcor,xcor,zwater2,xwater2,kk)
	    nwater(i)=nint(zcor/dx)+1
	 enddo
      endif

      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_surf

      implicit none
      include 'wtw_model.para'
      real eps_tmp
      parameter(eps_tmp=0.01)
      integer ind(ns0*(ng0+1))
      real xsurf(ns0*(ng0+1)),zsurf(ns0*(ng0+1)),xsurf2(ns0*(ng0+1)),
     1  zsurf2(ns0*(ng0+1))
      integer i, j, k, kk
      real zerodis, xcor, zcor

      zerodis=eps_tmp*dx
      k=0
      do i=1,ns
        k=k+1
        xsurf2(k)=xss(i)
        zsurf2(k)=zss(i)
        do j=1,ngis(i)
          if(offset(j,i).ge.-0.5.and.tr(j,i).ge.0.0)then
            k=k+1
            xsurf2(k)=xgg(j,i)
            zsurf2(k)=zgg(j,i)
          endif
        enddo
      enddo
      call sort(k,xsurf2,ind)
      kk=1
      xsurf(kk)=xsurf2(ind(1))
      zsurf(kk)=zsurf2(ind(1))
      do i=2,k
        if(abs(xsurf2(ind(i))-xsurf(kk)).gt.zerodis)then
          kk=kk+1
          xsurf(kk)=xsurf2(ind(i))
          zsurf(kk)=zsurf2(ind(i))
        endif
      enddo
      nx_tot=data_int(21)
      do i=1,nx_tot
        xcor=float(i-1)*dx
        call get1D(zcor,xcor,zsurf,xsurf,kk)
        nsurf_tot(i)=nint(zcor/dx)+1
      enddo

      open(20,file='surf.txt',form='formatted')
      do i=1,nx_tot
        write(20,*) nsurf_tot(i)
      enddo
      close(20)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_velocity

      implicit none
      include 'wtw_model.para'
c     ---------------------------------------------------------
      integer nx1, nz1, iii, isd, nx2, nz2, ix, iz
      real zz, xx, value
      real vel0(nx_mod0,nz_mod0)
      real den0(nx_mod0,nz_mod0)

      nx1=data_int(1)
      nz1=data_int(2)
      dx1=data_real(1)
      if(data_int(28).eq.0)then
        write(*,*) nx1, nx_mod0, data_int(21)
        call bird(vel0,vel_in,nx_mod0,nz_mod0,iii)
!        call bird3(vel0,vel_in,1,nx1,1,nz1,nx_mod0,nz_mod0,iii)
      else
        call surd11(vel0,vel_in,nx_mod0,nz_mod0,nx1,nz1,iii)
      endif
      if(iii.lt.nx1)then
        write(*,*)'Something wrong when reading file:',vel_in
        write(*,*)' only ',iii,' out of ',nx1,' were read.'
        stop
      endif
      ndens=data_int(10)
      if(ndens.eq.0)then
        if(data_int(29).eq.0)then
          call bird3(den0,density_in,1,nx1,1,nz1,nx_mod0,nz_mod0,iii)
        else
          call surd11(den0,density_in,nx1,nz1,nx_mod0,nz_mod0,iii)
        endif
        if(iii.lt.nx1)then
          write(*,*)'Something wrong when reading file:',density_in
          write(*,*)' only ',iii,' out of ',nx1,' were read.'
          stop
        endif
      endif
      if(is_reset_sg.eq.1)then
        isd=nint(depsou/dx)
      else
        isd=0
      endif
      nx2=data_int(21)
      nz2=data_int(22)+isd
      data_int(22)=nz2
      dx=data_real(3)
c	write(*,*)nx1,nz1,dx1,nx2,nz2,data_real(3)
      do iz=isd+1,nz2
        zz=float(iz-isd-1)*dx+data_real(26)-data_real(20)
        do ix=1,nx2
          xx=float(ix-1)*dx+data_real(25)-data_real(18)
          call get2Dvalue(value,xx,zz,vel0,nx_mod0,nz_mod0,dx1,dx1)
            velocity(ix,iz)=value
          if(ndens.eq.0)then
            call get2Dvalue(value,xx,zz,den0,nx_mod0,nz_mod0,dx1,dx1)
            density(ix,iz)=value
          endif
        enddo
      enddo
      do iz=1,isd
        do ix=1,nx2
          velocity(ix,iz)=velocity(ix,isd+1)
        enddo
      enddo
      if(ndens.eq.0)then
        do iz=1,isd
          do ix=1,nx2
            density(ix,iz)=density(ix,isd+1)
          enddo
        enddo
      endif

c Chaiwoot: set velocity
c      do ix=1,nx2
c        do iz=1,nsurf_tot(ix)
c          velocity(ix,iz)=vel0(ix+160,iz)
c        enddo
c      enddo

      write(*,*) nx2, nz2, nx_mod, nz_mod
      call biwt3(velocity,'VEL_INT',1,nx2,1,nz2, nx_mod,nz_mod)

      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source2
c
      include 'wtw_model.para'
      character*200 file_tmp
      real s(nt_max2),tdump(nt_max2),tdump2(nt_max2)
c       ---------------------------------------------------------
c Making a source wavelet 
c       ---------------------------------------------------------
      nt6=nts*6
      if(nt6.gt.nt_work) nt6=nt_work
      do iijj=is_first,is_end
        do it=1,nt_work
          sou(iijj-is_first+1,it)=0.
        enddo
      enddo
c
      do it=1,nt_work
        tdump2(it)=0.0
      enddo
c
      if(nopsou.eq.0) then
        b=(3.1415926*vm)**2
        a=2.*b
        do i=1,nts+1
          t2=((i-1)*dt)**2
          tdump2(nts+i)=(1.-a*t2)*exp(-b*t2)
          tdump2(nts+2-i)=tdump2(nts+i)
        enddo
        do iijj=is_first,is_end
          sou(iijj-is_first+1,1+ntpad)=tdump2(1)
          do it=2,2*nts
            sou(iijj-is_first+1,it+ntpad)=tdump2(it)
          enddo
        enddo
        nwave=max(2*nts,ntpsou)
      else if( nopsou.eq.2)then
        do it=1,nt_max2
          tdump(it)=float(it-1)*dt_sou
        enddo
        do iijj=is_first,is_end
          do i=1,nt_work
            sou(iijj-is_first+1,i)=0.
          enddo
        enddo
        do iijj=is_first,is_end
c          is=1+(iijj-1)*nshot_int
          is=ns-(iijj-1)*nshot_int
          call filename(file_tmp,source_in,indshot(is),'')
          open(66,file=file_tmp)
          do it=1,nt
            read(66,*,err=502,end=502)xxx,sss
            s(it)=sss
            nntsou=nt
          enddo
502       continue
          close(66)
          if(abs(dt_sou-dt).lt.0.0001)then
            do it=1,nntsou
              sou(iijj-is_first+1,it)=s(it)
            enddo
            do it=1+nntsou,nt_work
              sou(iijj-is_first+1,it)=0.0
            enddo
          else
            nntsou_new=nint((nntsou-1)*dt_sou/dt)+1
            nntsou_new=min(nt_work,nntsou_new)
            do it=1,nntsou_new
              tt=float(it-1)*dt
              call get1D(sss,tt,s,tdump,nntsou)
              sou(iijj-is_first+1,it)=sss
            enddo
            do it=1+nntsou_new,nt_work
              sou(iijj-is_first+1,it)=0.0
            enddo
          endif
        enddo
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_message

      implicit none
      include 'wtw_model.para'
      real dt_sou

      open(unit=30,file=wtw_ou)
      close(30,status='delete')
      open(unit=30,file=wtw_ou)
      open(unit=31,file=res_wave)
      write(30,*)'COORD_TEMP=',coord_temp
      write(30,*)'VEL_IN=',vel_in
      write(30,*)'IS_SU_VELOCITY=',data_int(28)
      write(30,*)'Nx_IN,Nz_IN,Dx_IN=',nx_in,nz_in,dx1
      write(30,*)'CSG_PRE_IN=',csg_pre
      write(30,*)'CSG_SUFFIX_IN=',csg_suffix
      write(30,*)'Nt_IN,Dt_IN,IS_SU_CSG=',nt,dt1,is_su_csg
      write(30,*)'Nx_OUT,Nz_OUT,Dx_OUT=',nx_tot,nz_tot,dx
      write(30,*)'Nt_OUT,Dt_OUT=',nt_work,dt
      write(30,*)'SOURCE_IN=',source_in
      write(30,*)'DENSITY_IN=',density_in
      write(30,*)'IS_SU_DENSITY=',data_int(29)
      write(30,*)'WTW_OU=',wtw_ou
      write(30,*)'RES_WAVEFORM=',res_wave
      write(30,*)'SYN_OU_PRE=',syn_ou_pre
      write(30,*)'SYN_OU_SUFFIX=',syn_ou_suffix
      write(30,*)'VEL_OU_PRE=',vel_ou_pre
      write(30,*)'VEL_OU_SUFFIX=',vel_ou_suffix
      write(30,*)'GRAD_OU=',grad_ou
      write(30,*)'ISMARINE,IS_RESET_SURFACE,IS_RESET_SG=',
     1		is_marine,is_reset_surface,is_reset_sg
      write(30,*)'DEPSOU,DEPREC=',depsou,deprec
      write(30,*)'NOPSOU,FREQ,NTPSOU,DTSOU=',nopsou,vm,ntpsou,dt_sou
      write(30,*)'NTPAD=',ntpad
      write(30,*)'NDENS,DENCO,CP00,DENP=',ndens,denco,cp00,denp
      write(30,*)'NFREE=',nfree
      write(30,*)'ITERATION,ITER00=',nit,iter00
      write(30,*)'NITE=',nite
      write(30,*)'IS_USE_CON1,DEPTH_CON1,VMIN_CON1,VMAX_CON1=',
     1	is_use_constraint1,depth_constraint1,vmin_constraint1,
     2  vmax_constraint1
      write(30,*)'DVMIN,DVMAX,VMIN,VMAX=',dvmin,dvmax,vmin,vmax
      write(30,*)'XMIN,XMAX,ZMIN,ZMAX=',xmin,xmax,zmin,zmax
      write(30,*)'OFFMIN,OFFMAX=',offmin,offmax
      write(30,*)'ISPRE,ISNORM,ISUSECONJUGATE,ISPRECONDITION=',
     1	ispre,isusenorm,isuseconjugate,isuseprecondition
      write(30,*)'NSHOT_INT=',nshot_int
      write(30,*)'INDEX_SHOT_OUTPUT=',ishot_out
      write(30,*)'Ns=',ns

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_assigned(is_first,is_end,is1,is2,ip,ime)

      implicit none
      integer ip, is_first, is_end, is1, is2, ime, nnn, ii1, ii2
c
c subroutine to assign job from is1 to is2 to Node No. ime in ip nodesc
c
      if(ip.eq.1)then
        is_first=is1
        is_end=is2
        return
      endif
      if(ip.le.0.or.ime.lt.0)then
        write(*,*)'Something wrong! ip should be positive and
     1          ime should be nonnegative.'
        return
      endif
      nnn=is2-is1+1
      ii1=int(nnn/ip)
      if(ii1.eq.0)then
        is_first=is1+ime
        is_end=is_first
        return
      endif
      ii2=nnn-ii1*ip
      if(ime.lt.(ip-ii2))then
        is_first=is1+ime*ii1
        is_end=is1+(ime+1)*ii1-1
      else 
        is_first=is1+ii2+ime*ii1+ime-ip
        is_end=is1+ii2+(ime+1)*ii1+ime-ip
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function mod2(i,j)

      implicit none
      integer i, j, mod2
      if(j.eq.1) then
        mod2=0
        return
      endif
      if(i.ge.0)then
        mod2=mod(i,j)
      else
        mod2=j+mod(i,j)
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sort(n,arr,ind)

      implicit none
      real aln2i, tiny
      parameter (aln2i=1./0.69314718,tiny=1.e-5)
      integer arr(1),ind(1)
      integer i, lognb2, n, m, nn, j, k, l, n1

      lognb2=int(alog(float(n))*aln2i+tiny)
      do i=1,n
        ind(i)=i
      enddo
      m=n
      do 12 nn=1,lognb2
        m=m/2
        k=n-m
        do 11 j=1,k
          i=j
3         continue
          l=i+m
          if(arr(ind(l)).lt.arr(ind(i)))then
            n1=ind(i)
            ind(i)=ind(l)
            ind(l)=n1
            i=i-m
            if(i.ge.1)goto 3
          endif
11      continue
12    continue
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 1D linear interpolation

      subroutine get1D(y,x,ycor,xcor,n)

      implicit none
      integer n, i
      real x, y, ycor(1), xcor(1), a, b

      if(xcor(1).le.xcor(n))then
        if(x.le.xcor(1))then
          y=ycor(1)
          return
        else if(x.ge.xcor(n))then
          y=ycor(n)
          return
        else
          do i=2,n
            if(x.le.xcor(i))then
              a=(xcor(i)-x)/(xcor(i)-xcor(i-1))
              b=1.0-a
              y=ycor(i)*b+ycor(i-1)*a
              return
            endif
          enddo
        endif
      else
        if(x.ge.xcor(1))then
          y=ycor(1)
          return
        else if(x.le.xcor(n))then
          y=ycor(n)
          return
        else
          do i=2,n
            if(x.ge.xcor(i))then
              a=(xcor(i)-x)/(xcor(i)-xcor(i-1))
              b=1.0-a
              y=ycor(i)*b+ycor(i-1)*a
              return
            endif
          enddo
        endif
      endif
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get2Dvalue(value,x,y,data,nx,ny,dx,dy)

      implicit none
      integer nx, ny, i11, i12, i21, i22, i31, i32, i41, i42
      real x, y, value, data(nx,ny), dx, dy, a, b, c, d

      i11=int(x/dx)+1
      i12=int(y/dy)+1
      if(i11.le.1)i11=1
      if(i11.ge.nx)i11=nx
      if(i12.le.1)i12=1
      if(i12.ge.ny)i12=ny
      i21=i11+1
      if(i21.ge.nx)i21=nx
      i22=i12
      i31=i21
      i32=i12+1
      if(i32.ge.ny)i32=ny
      i41=i11
      i42=i32
      a=x/dx-int(x/dx)
      b=1.0-a
      c=y/dy-int(y/dy)
      d=1.0-c
      value=d*(b*data(i11,i12)+a*data(i21,i22))
     1       +c*(b*data(i41,i42)+a*data(i31,i32))
      end
