      subroutine forw_back(c41,c42,iii_nite,resd1,resd2)

      implicit none
      include 'wtw_model.para'
      character*200 file_tmp
      integer i, j, iss
      real seis(2000,50), c41, c42, iii_nite, resd1, resd2, resid

      resd1=0.0
      resd2=0.0

      do iss=is_first,is_end
        is=ns-(iss-1)*nshot_int
        is00=iss-is_first+1
        nxs=nxss(is)
        nzs=nzss(is)
        ng=ngis(is)

        call c1c2cl_shot(is,c41,c42)
        call forw(c41,c42,1)
c        call mute(wc,nt_work,ngis(is))
        call filename(file_tmp,
     1          syn_ou_pre,indshot(is),syn_ou_suffix)
        call suwt22(wc,file_tmp,1.0,-1.0*ntpad*dt,1.0,dt,ngis(is),
     1         nt_work,ng0,nt_max)
!        open(10,file=file_tmp,access='direct',recl=4*nt_work)
!        do i=1,ngis(is)
!          write(10,rec=i) (wc(j,i),j=1,nt_work)
!        enddo
!        close(10)
      enddo

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

