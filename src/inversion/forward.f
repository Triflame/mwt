      subroutine forw_back(c41,c42,iii_nite,resd1,resd2)

      implicit none
      include 'wtw_model.para'
      character*200 file_tmp
      integer i, j, iss
      real seis(2000,50), c41, c42, iii_nite, resd1, resd2, resid

      resd1=0.0
      resd2=0.0

      do iss=is_first,is_end
c        is=1+(iss-1)*nshot_int
        is=ns-(iss-1)*nshot_int
        is00=iss-is_first+1
        nxs=nxss(is)
        nzs=nzss(is)
        ng=ngis(is)

        call c1c2cl_shot(is,c41,c42)
        call forw(c41,c42,1)
c        call mute(wc,nt_max,ngis(is))
        call norm12
        if(indshot(is).eq.ishot_out)then
          call filename(file_tmp,
     1          syn_ou_pre,indshot(is),syn_ou_suffix)
          call suwt22(wc,file_tmp,1.0,-1.0*ntpad*dt,1.0,dt,ngis(is),
     1         nt_work,ng0,nt_max)
        endif

        call residual(resid,iii_nite,resd1,resd2)
        call back(c41,c42)
888     continue
      enddo

999   continue
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fun_step(alp,c41,c42,iii_nite)

      include 'wtw_model.para'
      real fun_step
      character*200 file_tmp

      s1=1.0/vmin
      s2=1.0/vmax
      do ix=1,nx_tot
        do iz=nsurf_tot(ix),nz_tot
          tmp=vdump(ix,iz)-d2(ix,iz)*alp
          tmp=min(tmp,s1)
          tmp=max(tmp,s2)
          velocity(ix,iz)=1.0/tmp
        enddo
      enddo
      call pad_model

c      call c1c2cl(c41,c42)
      fun_step=0.0

      do iss=iter-1-iter00,iter-1-iter00+iii_nite
     1          *(nite-1),iii_nite
        iss2=is_first+mod2(iss,(is_end-is_first+1))
c        is=1+(iss2-1)*nshot_int
        is=ns-(iss2-1)*nshot_int
        is00=iss2-is_first+1
        nxs=nxss(is)
        nzs=nzss(is)
        ng=ngis(is)
        call c1c2cl_shot(is,c41,c42)
        call forw(c41,c42,0)
        call norm12
        call residual22(rrr)
        fun_step=fun_step+rrr
c        write(*,*)'step:is,rrr,fun_step=',is,rrr,fun_step
      enddo
c      write(*,*)'is,alp,fun_step=',is,alp,fun_step

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine forw(c41,c42,is_store_boundary)

      implicit none
      include 'wtw_model.para'
      character*200 file_tmp
      integer ix, iz, it, is_store_boundary, ixx, izz, k, jj, j, iz1,
     1  iz2, iz3
      real c41, c42
c     ----------------------------------------------------------------
c     Forward modeling to get synthetic seismograms
c     -----------------------------------------------------------------
c
      do iz=1,nz
        do ix=1,nx
          p2(ix,iz)=0.
          u1(ix,iz)=0.
          w1(ix,iz)=0.
        enddo
      enddo
      do ix=1,na_pml
        do iz=1,nz
          p2xL(ix,iz)=0.0
          p2xR(ix,iz)=0.0
          p2yL(ix,iz)=0.0
          p2yR(ix,iz)=0.0
        enddo
      enddo
      do ix=1,nx
        do iz=1,na_pml
          p2xT(ix,iz)=0.0
          p2yT(ix,iz)=0.0
          p2xB(ix,iz)=0.0
          p2yB(ix,iz)=0.0
        enddo
      enddo

c Main Loop
c      write(*,*) 'nt_work = ', nt_work
c      do 10 it=1,nt_work+timeshift
      do 10 it=1,nt_work
        do iz=nzm10,nzm20
          do ix=nxm10,nxm20
            p2(ix,iz)=p2(ix,iz)+c2(ix,iz)*(
     1        c41*( u1(ix  ,iz)-u1(ix-1,iz)+w1(ix,iz  )-w1(ix,iz-1))
     2       +c42*( u1(ix+1,iz)-u1(ix-2,iz)+w1(ix,iz+1)-w1(ix,iz-2)))
          enddo
        enddo
c
c adding source
c
c        p2(nxs,nzs)=p2(nxs,nzs)+sou(is00,it)*c2(nxs,nzs)*0.001*dx
        p2(nxs,nzs)=p2(nxs,nzs)+sou(1,it)*c2(nxs,nzs)*0.001*dx
c
c p2x
c
        do iz=na_pml,2,-1
          do ix=nxm10,nxm20
            p2xT(ix,iz)=p2xT(ix,iz)+c2(ix,iz)*(
     1       c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
          enddo
        enddo
        do iz=2,na_pml
          do ix=na_pml,3,-1
            p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+c2(ix,iz)*(
     1       c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz))) 
          enddo
          ix=2
          p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+
     1          c2(ix,iz)*(u1(ix,iz)-u1(ix-1,iz))
        enddo
        do ix=nxm20+1,nx-1
          ixx=nx-ix+1
          do iz=2,na_pml
            p2xR(ixx,iz)=p2xR(ixx,iz)*exR(ixx,iz)+c2(ix,iz)*(
     1       c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
          enddo
        enddo
        do iz=nzm20+1,nz-1
           izz=nz-iz+1
	   do ix=nxm10,nxm20
	      p2xB(ix,izz)=p2xB(ix,izz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
	do ix=na_pml,3,-1
	   do iz=na_pml+1,nz-1
	      p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
	ix=2
	do iz=na_pml+1,nz-1
	   p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+
     1		c2(ix,iz)*(u1(ix,iz)-u1(ix-1,iz))
	enddo
	do ix=nxm20+1,nx-1
	   ixx=nx-ix+1
	   do iz=na_pml+1,nz-1
	      p2xR(ixx,iz)=p2xR(ixx,iz)*exR(ixx,iz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
c
c p2y
c
	   do iz=na_pml,3,-1
	   do ix=2,nx-1
	      p2yT(ix,iz)=p2yT(ix,iz)*eyT(ix,iz)+c2(ix,iz)*(
     1	   c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
	   enddo
	   enddo
	   iz=2
	   do ix=2,nx-1
	   p2yT(ix,iz)=p2yT(ix,iz)*eyT(ix,iz)+
     1		c2(ix,iz)*(w1(ix,iz)-w1(ix,iz-1))
	   enddo
	do iz=nzm20+1,nz-1
	   izz=nz-iz+1
           do ix=2,nx-1
              p2yB(ix,izz)=p2yB(ix,izz)*eyB(ix,izz)+c2(ix,iz)*(
     1     c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
           enddo
        enddo
	do ix=na_pml,2,-1
	   do iz=nzm10,nzm20
	      p2yL(ix,iz)=p2yL(ix,iz)+c2(ix,iz)*(
     1	   c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
	   enddo
	enddo
        do ix=nxm20+1,nx-1
	   ixx=nx-ix+1
           do iz=nzm10,nzm20
              p2yR(ixx,iz)=p2yR(ixx,iz)+c2(ix,iz)*(
     1     c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
           enddo
        enddo
c
	do iz=2,na_pml
	   do ix=2,na_pml
	      p2(ix,iz)=p2xL(ix,iz)+p2yT(ix,iz)
	   enddo
	   do ix=nxm10,nxm20
	      p2(ix,iz)=p2xT(ix,iz)+p2yT(ix,iz)
	   enddo
	   do ix=nxm20+1,nx
	      ixx=nx-ix+1
	      p2(ix,iz)=p2xR(ixx,iz)+p2yT(ix,iz)
	   enddo
	enddo
	do iz=nzm20+1,nz
	   izz=nz-iz+1
           do ix=2,na_pml
              p2(ix,iz)=p2xL(ix,iz)+p2yB(ix,izz)
           enddo
           do ix=nxm10,nxm20
              p2(ix,iz)=p2xB(ix,izz)+p2yB(ix,izz)
           enddo
           do ix=nxm20+1,nx
              ixx=nx-ix+1
              p2(ix,iz)=p2xR(ixx,iz)+p2yB(ix,izz)
           enddo
        enddo   
	do iz=nzm10,nzm20
	   do ix=2,na_pml
	      p2(ix,iz)=p2xL(ix,iz)+p2yL(ix,iz)
	   enddo
	   do ix=nxm20+1,nx
              ixx=nx-ix+1
	      p2(ix,iz)=p2xR(ixx,iz)+p2yR(ixx,iz)
	   enddo
	enddo
c
        if(nfree.eq.1)then
           do ix=1,nx
              p2(ix,nsurf(ix))=0.0
c              p2(ix,nsurf(ix)-1)=p2(ix,nsurf(ix)+1)

c Changed by Chaiwoot (06/25/07)
              p2(ix,nsurf(ix)-1) = -p2(ix,nsurf(ix)+1)
           enddo
        endif

        if(ispre.eq.1)then
c
c  output the seismogram (pressure)
c
c          if (it .gt. timeshift) then
          do k=1,ng
            if(itmmm(k,is00).gt.0)then
            wc(it,k)= p2(nxg(k,is),nzg(k,is))
c            wc(it-timeshift,k)= p2(nxg(k,is),nzg(k,is))
            endif
          enddo
c          endif
        endif

      do 4222 iz=nzm10,nzm20
      do 4222 ix=nxm10,nxm20
      u1(ix,iz)=u1(ix,iz)+  
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
4222  continue
      do 4221 iz=nzm10,nzm20
      do 4221 ix=nxm10,nxm20
      w1(ix,iz)=w1(ix,iz)+ 
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )  
4221  continue

c
c u1
c
	do iz=na_pml,2,-1
	   do ix=nxm10,nxm20
	      u1(ix,iz)=u1(ix,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
	   enddo
	enddo
	do iz=nzm20+1,nz-1
	   do ix=nxm10,nxm20
	      u1(ix,iz)=u1(ix,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
	   enddo
	enddo
	do ix=na_pml,2,-1
	   do iz=2,na_pml
	      u1(ix,iz)=u1(ix,iz)*exL(ix,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
	   enddo
	   do iz=na_pml+1,nz-1
	      u1(ix,iz)=u1(ix,iz)*exL(ix,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )  
	   enddo
	enddo
	do ix=nxm20+1,nx-2
	   ixx=nx-ix+1
           do iz=2,na_pml
              u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
           do iz=na_pml+1,nz-1
              u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
        enddo
	ix=nx-1
	ixx=nx-ix+1
	do iz=2,na_pml
	   u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+
     1            cl(ix,iz)*(p2(ix+1,iz)-p2(ix,iz))
	enddo
	do iz=na_pml+1,nz-1
	   u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+
     1            cl(ix,iz)*(p2(ix+1,iz)-p2(ix,iz))
	enddo
c
c w1
c
	do iz=na_pml,2,-1
	   do ix=2,nx-1
              w1(ix,iz)=w1(ix,iz)*eyT(ix,iz)+ 
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )  
	   enddo
	enddo
	do iz=nzm20+1,nz-2
	   izz=nz-iz+1
	   do ix=2,nx-1
	      w1(ix,iz)=w1(ix,iz)*eyB(ix,izz)+
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
        enddo
	iz=nz-1
	izz=nz-iz+1
	do ix=2,nx-1
	   w1(ix,iz)=w1(ix,iz)*eyB(ix,izz)+
     1 cl(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
	enddo
	do ix=na_pml,2,-1
	   do iz=nzm10,nzm20
	      w1(ix,iz)=w1(ix,iz)+
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
        enddo
	do ix=nxm20+1,nx
	   do iz=nzm10,nzm20
	      w1(ix,iz)=w1(ix,iz)+
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
	enddo
c
        if(nfree.eq.1)then
           do ix=1,nx
              iz=nsurf(ix)-1
              w1(ix,iz) = w1(ix,nsurf(ix)+1)
c              w1(ix,iz)=w1(ix,iz)+
c     1           cl41(ix,iz)*(-p2(ix,iz))
c     2          +cl42(ix,iz)*(p2(ix,iz+2)+p2(ix,iz+3))
           enddo
        end if
c
	if(ispre.ne.1)then
c
c  output the seismogram (particle velocity)
c
          do k=1,ng
	    if(itmmm(k,is00).gt.0)then
            wc(it,k)= w1(nxg(k,is),nzg(k,is))
	    endif
	  enddo
	endif

      if(is_store_boundary.eq.1) then
c
c     store records and boundary value
c
        jj=0
        do j=na+1,nx-na
          bord4(jj+1,it)=p2(j,nz  -na)  
          bord4(jj+2,it)=p2(j,nz-1-na)  
          bord4(jj+3,it)=w1(j,nz-1-na)
          jj=jj+3
        enddo
        do j=k0+1,nz-na
          bord4(jj+1,it)=p2(na+1,j)  
          bord4(jj+2,it)=p2(na+2,j)
          bord4(jj+3,it)=u1(na+1,j)  
          bord4(jj+4,it)=p2(nx  -na,j)  
          bord4(jj+5,it)=p2(nx-1-na,j)  
          bord4(jj+6,it)=u1(nx-1-na,j)
          jj=jj+6
        enddo
        do ix=na+1,nx-na
          iz1=nsurf(ix)
          iz2=iz1+1
          iz3=iz1+2
          bord4(jj+1,it)=p2(ix,iz1)
          bord4(jj+2,it)=p2(ix,iz2)
          bord4(jj+3,it)=p2(ix,iz3)
          bord4(jj+4,it)=u1(ix,iz1)
          bord4(jj+5,it)=u1(ix,iz2)
          bord4(jj+6,it)=w1(ix,iz1)
          jj=jj+6
        enddo

        endif
c End of Main Loop
c
10    continue

      do it=nt_work-timeshift,nt_work
        do k=1,ng
          wc(it,k)=0.0
        enddo
      enddo

999   continue
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine back(c41,c42)

      implicit none
      include 'wtw_model.para'
c      -----------------------------------------------------------------
c      Resortion of the forward field, back projection of the pseudo-traveltime
c      residuals, calculation of the gradient g1
c      ------------------------------------------------------------------
c
      real pp1(nx0,nz0),pp2(nx0,nz0),pp0(nx0,nz0),
     1     uu1(nx0,nz0),uu2(nx0,nz0),uu0(nx0,nz0),
     2     ww1(nx0,nz0),ww2(nx0,nz0),ww0(nx0,nz0)
      character*200 file_tmp
      integer nt_work0, ks, ke, js, je, ix, iz, k, it, jj, j, iz1, iz2,
     1  iz3, ixx, izz
      real c41, c42, rt, aa, dd, bb, cc
      integer min
      intrinsic min

      if(isusenorm.eq.1)then
        nt_work0=2*nt_work-100
      else
        nt_work0=nt_work
      endif

      rt=1./2./dt
      ks=nzm1
      ke=nzm2
      js= nxm1
      je=nxm2
      do iz=1,nz
        do ix=1,nx
          pp1(ix,iz)=p2(ix,iz)
          uu1(ix,iz)=u1(ix,iz)
          ww1(ix,iz)=w1(ix,iz)
          p2(ix,iz)=0.
          u1(ix,iz)=0.
          w1(ix,iz)=0.
          pp0(ix,iz)=0.0
          uu0(ix,iz)=0.0
          ww0(ix,iz)=0.0
        enddo
      enddo

      do ix=1,na_pml
        do iz=1,nz
          p2xL(ix,iz)=0.0
          p2xR(ix,iz)=0.0
          p2yL(ix,iz)=0.0
          p2yR(ix,iz)=0.0
        enddo
      enddo
      do ix=1,nx
        do iz=1,na_pml
          p2xT(ix,iz)=0.0
          p2yT(ix,iz)=0.0
          p2xB(ix,iz)=0.0
          p2yB(ix,iz)=0.0
        enddo
      enddo

c ispre = 1 for pressure field
      if(ispre.eq.1)then
        do k=1,ng
          if(itmmm(k,is00).gt.0)then
c dw is wavefield residual
            p2(nxg(k,is),nzg(k,is))=dw(nt_work0,k)
          endif
        enddo
      else
        do k=1,ng
          if(itmmm(k,is00).gt.0)then
            w1(nxg(k,is),nzg(k,is))=dw(nt_work0,k)
          endif
        enddo
      endif

c Main loop for back propagation of residual
      do 10 it=nt_work0,2,-1
        if(it.le.nt_work) then
          do ix=na+2,nx-na-2
            do iz=nsurf(ix)+2,nz-na-2
              uu2(ix,iz)=uu1(ix,iz)-  
     1          cl41(ix,iz)*( pp1(ix+1,iz)-pp1(ix  ,iz) )
     2         -cl42(ix,iz)*( pp1(ix+2,iz)-pp1(ix-1,iz) )
            enddo
          enddo

          do ix=na+3,nx-na-2
            do iz=nsurf(ix)+1,nz-na-2
              ww2(ix,iz)=ww1(ix,iz)-  
     1          cl41(ix,iz)*( pp1(ix,iz+1)-pp1(ix,iz  ) )
     2         -cl42(ix,iz)*( pp1(ix,iz+2)-pp1(ix,iz-1) )
            enddo
          enddo

          jj=0
          do j=na+1,nx-na
            ww2(j,nz-na-1)=bord4(jj+3,it-1)
            jj=jj+3
          enddo

          do j=k0+1,nz-na
            uu2(na+1,j)=bord4(jj+3,it-1)
            uu2(nx-na-1,j)=bord4(jj+6,it-1)
            jj=jj+6
          enddo

          do ix=na+1,nx-na
            iz1=nsurf(ix)
            iz2=iz1+1
            iz3=iz1+2
            uu2(ix,iz1)=bord4(jj+4,it-1)
            uu2(ix,iz2)=bord4(jj+5,it-1)
            ww2(ix,iz1)=bord4(jj+6,it-1)
            jj=jj+6
          enddo

c          pp1(nxs,nzs)=pp1(nxs,nzs)-sou(is00,it)
c     1      *c2(nxs,nzs)*0.001*dx
          pp1(nxs,nzs)=pp1(nxs,nzs)-sou(1,it)
     1      *c2(nxs,nzs)*0.001*dx

          do ix=na+3,nx-2-na
            do iz=nsurf(ix)+3,nz-2-na
              pp2(ix,iz)=pp1(ix,iz)-c2(ix,iz)*(
     1  c41*( uu2(ix  ,iz)-uu2(ix-1,iz)+ww2(ix,iz  )-ww2(ix,iz-1) )
     2 +c42*( uu2(ix+1,iz)-uu2(ix-2,iz)+ww2(ix,iz+1)-ww2(ix,iz-2) ) )
            enddo
          enddo

          jj=0
          do j=na+1,nx-na
            pp2(j,nz-na  )=bord4(jj+1,it-1)
            pp2(j,nz-na-1)=bord4(jj+2,it-1)
            jj=jj+3
          enddo

          do j=k0+1,nz-na
            pp2(na+1,j)=bord4(jj+1,it-1)
            pp2(na+2,j)=bord4(jj+2,it-1)
            pp2(nx-na  ,j)=bord4(jj+4,it-1)
            pp2(nx-na-1,j)=bord4(jj+5,it-1)
            jj=jj+6
          enddo

          do ix=na+1,nx-na
            iz1=nsurf(ix)
            iz2=iz1+1
            iz3=iz1+2
            pp2(ix,iz1)=bord4(jj+1,it-1)
            pp2(ix,iz2)=bord4(jj+2,it-1)
            pp2(ix,iz3)=bord4(jj+3,it-1)
            jj=jj+6
          enddo

c          pp1(nxs,nzs)=pp1(nxs,nzs)+sou(is00,it)
c     1        *c2(nxs,nzs)*0.001*dx
          pp1(nxs,nzs)=pp1(nxs,nzs)+sou(1,it)
     1        *c2(nxs,nzs)*0.001*dx

c For particle velocity
          if(it.lt.nt_work-1.and.it.gt.2.and.ispre.ne.1)then
            do ix=js+npad_shot,je-npad_shot
              ixx=ix-nxm1+ix1_shot(is)
              ixx=max(1,ixx)
              ixx=min(ixx,nx_tot)
              do iz=nsurf(ix),ke
                aa=p2(ix,iz)*(pp0(ix,iz)-pp2(ix,iz))/2.0
                dd=den(ix,iz)*c1(ix,iz)
                bb=dd*u1(ix,iz)*(dd*uu2(ix,iz)-dd*uu0(ix,iz))/2.0
                cc=dd*w1(ix,iz)*(dd*ww2(ix,iz)-dd*ww0(ix,iz))/2.0
                g1(ixx,iz-nzm1+1)=g1(ixx,iz-nzm1+1)+aa*9.0+(bb+cc)
              enddo
            enddo
          endif

c For pressure
          if(it.lt.nt_work-1.and.it.gt.2.and.ispre.eq.1)then
            do ix=js+npad_shot,je-npad_shot
              ixx=ix-nxm1+ix1_shot(is)
              ixx=max(1,ixx)
              ixx=min(ixx,nx_tot)
              do iz=nsurf(ix),ke
                aa=p2(ix,iz)*(pp0(ix,iz)-pp2(ix,iz))/2.0
                dd=den(ix,iz)*c1(ix,iz)
                bb=dd*u1(ix,iz)*(dd*uu2(ix,iz)-dd*uu0(ix,iz))/2.0
                cc=dd*w1(ix,iz)*(dd*ww2(ix,iz)-dd*ww0(ix,iz))/2.0
c	         g1(ix,iz)=d2(ix,iz)+aa/100000.0
                g1(ixx,iz-nzm1+1)=g1(ixx,iz-nzm1+1)+aa*9.0+(bb+cc)
              enddo
            enddo
          endif

          do iz=k0+1,nz-na
            do ix=na+1,nx-na
              pp0(ix,iz)=pp1(ix,iz)
              pp1(ix,iz)=pp2(ix,iz)
              uu0(ix,iz)=uu1(ix,iz)
              uu1(ix,iz)=uu2(ix,iz)
              ww0(ix,iz)=ww1(ix,iz)
              ww1(ix,iz)=ww2(ix,iz)
            enddo
          enddo
        endif
c
c Back propagation of the residuals calculated in sub. residual
c
      do iz=nzm10,nzm20
        do ix=nxm10,nxm20
          u1(ix,iz)=u1(ix,iz)+
     1      cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2     +cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
        enddo
      enddo
      do iz=nzm10,nzm20
        do ix=nxm10,nxm20
          w1(ix,iz)=w1(ix,iz)+
     1      cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2     +cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
        enddo
      enddo

      if(ispre.ne.1)then
c
c  particle velocity residual
        do k=1,ng
          if(itmmm(k,is00).gt.0)then
            w1(nxg(k,is),nzg(k,is))=w1(nxg(k,is),nzg(k,is))
     1       +dw(it,k)*wweight(k,is00)*c2(nxg(k,is),nzg(k,is))*dx
          endif
        enddo
      endif
c
c u1
c
        do iz=na_pml,2,-1
           do ix=nxm10,nxm20
              u1(ix,iz)=u1(ix,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
        enddo
        do iz=nzm20+1,nz-1
           do ix=nxm10,nxm20
              u1(ix,iz)=u1(ix,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
        enddo
        do ix=na_pml,2,-1
           do iz=2,na_pml
              u1(ix,iz)=u1(ix,iz)*exL(ix,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
           do iz=na_pml+1,nz-1
              u1(ix,iz)=u1(ix,iz)*exL(ix,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
        enddo
        do ix=nxm20+1,nx-2
           ixx=nx-ix+1
           do iz=2,na_pml
              u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
           do iz=na_pml+1,nz-1
              u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+
     1 cl41(ix,iz)*( p2(ix+1,iz)-p2(ix  ,iz) )
     2+cl42(ix,iz)*( p2(ix+2,iz)-p2(ix-1,iz) )
           enddo
        enddo
        ix=nx-1
        ixx=nx-ix+1
        do iz=2,na_pml
           u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+
     1            cl(ix,iz)*(p2(ix+1,iz)-p2(ix,iz))
        enddo
        do iz=na_pml+1,nz-1
           u1(ix,iz)=u1(ix,iz)*exR(ixx,iz)+
     1            cl(ix,iz)*(p2(ix+1,iz)-p2(ix,iz))
        enddo
c
c w1
c
        do iz=na_pml,2,-1
           do ix=2,nx-1
              w1(ix,iz)=w1(ix,iz)*eyT(ix,iz)+
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
        enddo
        do iz=nzm20+1,nz-2
           izz=nz-iz+1
           do ix=2,nx-1
              w1(ix,iz)=w1(ix,iz)*eyB(ix,izz)+
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
        enddo
        iz=nz-1
        izz=nz-iz+1
        do ix=2,nx-1
           w1(ix,iz)=w1(ix,iz)*eyB(ix,izz)+
     1 cl(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
        enddo
        do ix=na_pml,2,-1
           do iz=nzm10,nzm20
              w1(ix,iz)=w1(ix,iz)+
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
        enddo
        do ix=nxm20+1,nx-1
           do iz=nzm10,nzm20
              w1(ix,iz)=w1(ix,iz)+
     1 cl41(ix,iz)*( p2(ix,iz+1)-p2(ix,iz  ) )
     2+cl42(ix,iz)*( p2(ix,iz+2)-p2(ix,iz-1) )
           enddo
        enddo
c
        if(nfree.eq.1)then
           do ix=1,nx
              iz=nsurf(ix)-1
              w1(ix,iz) = w1(ix,nsurf(ix)+1)
c              w1(ix,iz)=w1(ix,iz)+
c     1           cl41(ix,iz)*(-p2(ix,iz))
c     2          +cl42(ix,iz)*(p2(ix,iz+2)+p2(ix,iz+3))
           enddo
        end if
c
      	do iz=nzm10,nzm20
      	do ix=nxm10,nxm20
      p2(ix,iz)=p2(ix,iz)+c2(ix,iz)*(
     1  c41*( u1(ix  ,iz)-u1(ix-1,iz)+w1(ix,iz  )-w1(ix,iz-1) )
     2 +c42*( u1(ix+1,iz)-u1(ix-2,iz)+w1(ix,iz+1)-w1(ix,iz-2) ) )
	enddo
	enddo
c
	if(ispre.eq.1)then
c
c  pressure residual
          do k=1,ng
	    if(itmmm(k,is00).gt.0)then
               p2(nxg(k,is),nzg(k,is))=p2(nxg(k,is),nzg(k,is))
     1      +dw(it,k)*c2(nxg(k,is),nzg(k,is))*dx*wweight(k,is00)
	    endif
          enddo
	endif
c
c p2x
c
	do iz=na_pml,2,-1
	   do ix=nxm10,nxm20
	      p2xT(ix,iz)=p2xT(ix,iz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
        do iz=nzm20+1,nz-1
           izz=nz-iz+1
	   do ix=nxm10,nxm20
	      p2xB(ix,izz)=p2xB(ix,izz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
	do ix=na_pml,3,-1
	   do iz=2,na_pml
	      p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	   do iz=na_pml+1,nz-1
	      p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
	ix=2
	do iz=2,na_pml
	   p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+
     1		c2(ix,iz)*(u1(ix,iz)-u1(ix-1,iz))
	enddo
	do iz=na_pml+1,nz-1
	   p2xL(ix,iz)=p2xL(ix,iz)*exL(ix,iz)+
     1		c2(ix,iz)*(u1(ix,iz)-u1(ix-1,iz))
	enddo
	do ix=nxm20+1,nx-1
	   ixx=nx-ix+1
	   do iz=2,na_pml
	      p2xR(ixx,iz)=p2xR(ixx,iz)*exR(ixx,iz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	   do iz=na_pml+1,nz-1
	      p2xR(ixx,iz)=p2xR(ixx,iz)*exR(ixx,iz)+c2(ix,iz)*(
     1	   c41*(u1(ix,iz)-u1(ix-1,iz))+c42*(u1(ix+1,iz)-u1(ix-2,iz)))
	   enddo
	enddo
c
c p2y
c
	do iz=na_pml,3,-1
	   do ix=2,nx-1
	      p2yT(ix,iz)=p2yT(ix,iz)*eyT(ix,iz)+c2(ix,iz)*(
     1	   c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
	   enddo
	enddo
	iz=2
	do ix=2,nx-1
	   p2yT(ix,iz)=p2yT(ix,iz)*eyT(ix,iz)+
     1		c2(ix,iz)*(w1(ix,iz)-w1(ix,iz-1))
	enddo
	do iz=nzm20+1,nz-1
	   izz=nz-iz+1
           do ix=2,nx-1
              p2yB(ix,izz)=p2yB(ix,izz)*eyB(ix,izz)+c2(ix,iz)*(
     1     c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
           enddo
        enddo
	do ix=na_pml,2,-1
	   do iz=nzm10,nzm20
	      p2yL(ix,iz)=p2yL(ix,iz)+c2(ix,iz)*(
     1	   c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
	   enddo
	enddo
        do ix=nxm20+1,nx-1
	   ixx=nx-ix+1
           do iz=nzm10,nzm20
              p2yR(ixx,iz)=p2yR(ixx,iz)+c2(ix,iz)*(
     1     c41*(w1(ix,iz)-w1(ix,iz-1))+c42*(w1(ix,iz+1)-w1(ix,iz-2)))
           enddo
        enddo

        do iz=2,na_pml
           do ix=2,na_pml
              p2(ix,iz)=p2xL(ix,iz)+p2yT(ix,iz)
           enddo
           do ix=nxm10,nxm20
              p2(ix,iz)=p2xT(ix,iz)+p2yT(ix,iz)
           enddo
           do ix=nxm20+1,nx
              ixx=nx-ix+1
              p2(ix,iz)=p2xR(ixx,iz)+p2yT(ix,iz)
           enddo
        enddo
        do iz=nzm20+1,nz
           izz=nz-iz+1
           do ix=2,na_pml
              p2(ix,iz)=p2xL(ix,iz)+p2yB(ix,izz)
           enddo
           do ix=nxm10,nxm20
              p2(ix,iz)=p2xB(ix,izz)+p2yB(ix,izz)
           enddo
           do ix=nxm20+1,nx
              ixx=nx-ix+1
              p2(ix,iz)=p2xR(ixx,iz)+p2yB(ix,izz)
           enddo
        enddo
        do iz=nzm10,nzm20
           do ix=2,na_pml
              p2(ix,iz)=p2xL(ix,iz)+p2yL(ix,iz)
           enddo
           do ix=nxm20+1,nx
              ixx=nx-ix+1
              p2(ix,iz)=p2xR(ixx,iz)+p2yR(ixx,iz)
           enddo
        enddo

        if(nfree.eq.1)then
           do ix=1,nx
              p2(ix,nsurf(ix))=0.0
c              p2(ix,nsurf(ix)-1)=p2(ix,nsurf(ix)+1)

c Changed by Chaiwoot
              p2(ix,nsurf(ix)-1) = -p2(ix,nsurf(ix)+1)
           enddo
        endif
10    continue

      end

