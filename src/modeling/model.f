      subroutine pad_model

      implicit none
      include 'wtw_model.para'
      integer ix, iz
      real v1

      do ix=1,nx_tot
        do iz=1,nsurf_tot(ix)-1
          velocity(ix,iz)=airvel
        enddo
        do iz=nsurf_tot(ix),nz_tot
          v1 = max(velocity(ix,iz),vmin)
          velocity(ix,iz) = min(v1,vmax)
        enddo
      enddo

c  constraint III: shallow layer (like water layer)
c
c      if(is_use_constraint1.eq.1)then
c        n_constraint1=nint(depth_constraint1/dx)
c        do i=1,nx_tot
c          do j=nsurf_tot(i),nsurf_tot(i)+n_constraint1
c            velocity(i,j)=max(velocity(i,j),vmin_constraint1)
c            velocity(i,j)=min(velocity(i,j),vmax_constraint1)
c          enddo
c        enddo
c      endif

c      if(is_use_water_bottom.eq.1)then
c        do i=1,nx_tot
c          do j=nsurf_tot(i),nsurf_tot(i)+nwater(i)-1
c            velocity(i,j)=max(velocity(i,j),1490.0)
c            velocity(i,j)=min(velocity(i,j),1500.0)
c          enddo
c        enddo
c      endif

      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine c1c2cl_shot(ishot,c41,c42)

      implicit none
      include 'wtw_model.para'
      integer nxm, ishot, ix, iz, ixx, i, j, izz, n_constraint1
      real c41, c42, tmp0, alpha, tmp, alpha_max0, const, factor_att

      nxm=ix2_shot(ishot)-ix1_shot(ishot)+1
      nxm1=na_new+1
      nxm2=nxm1+nxm-1
      nxm10=nxm1-na_pad
      nxm20=nxm2+na_pad
      nx=nxm+2*na_new

      do ix=nxm1,nxm2
        ixx=ix-nxm1+ix1_shot(ishot)
        if(ixx.lt.1)then
          ixx=1
        else if(ixx.gt.nx_tot)then
          ixx=nx_tot
        endif
        do iz=nzm1,nzm2
          c1(ix,iz)=velocity(ixx,iz-nzm1+1)
        enddo
      enddo

      do ix=1,nx
        ixx=ix-nxm1+ix1_shot(ishot)
        if(ixx.lt.1)then
          ixx=1
        else if(ixx.gt.nx_tot)then
          ixx=nx_tot
        endif
        nsurf(ix)=nsurf_tot(ixx)+na_new
      enddo
c
c  constraint I: vmin,vmax
c
      do ix=nxm1,nxm2
        do iz=nsurf(ix),nz
          c1(ix,iz)=max(c1(ix,iz),vmin)
          c1(ix,iz)=min(c1(ix,iz),vmax)
        enddo
      enddo
c
c  constraint III: shallow layer (like water layer)
c
      if(is_use_constraint1.eq.1)then
        n_constraint1=nint(depth_constraint1/dx)
        do i=nxm1,nxm2
          do j=nsurf(i),nsurf(i)+n_constraint1
            c1(i,j)=max(c1(i,j),vmin_constraint1)
            c1(i,j)=min(c1(i,j),vmax_constraint1)
          enddo
        enddo
      endif
      if(is_use_water_bottom.eq.1)then
        do i=nxm1,nxm2
          ixx=i-nxm1+ix1_shot(ishot)
          if(ixx.lt.1)then
            ixx=1
          else if(ixx.gt.nx_tot)then
            ixx=nx_tot
          endif
          do j=nsurf(i),nsurf(i)+nwater(ixx)-1
            c1(i,j)=max(c1(i,j),1490.0)
            c1(i,j)=min(c1(i,j),1500.0)
          enddo
        enddo
      endif
c
c  pad the boundaries
c
      do ix=nxm1,nxm2
        do iz=nzm2+1,nz
          c1(ix,iz)=c1(ix,nzm2)
        enddo
        do iz=nzm1-1,1,-1
          c1(ix,iz)=c1(ix,nzm1)
        enddo
      enddo
      do iz=1,nz
        do ix=nxm1-1,1,-1
          c1(ix,iz)=c1(nxm1,iz)
        enddo
        do ix=nxm2+1,nx
          c1(ix,iz)=c1(nxm2,iz)
        enddo
      enddo
c
c  constraint II: above the surface
c
      if(nfree.ne.0)then
        do ix=1,nx
          do iz=1,nsurf(ix)
            c1(ix,iz)=0.0
          enddo
        enddo
      else
        do ix=1,nx
          do iz=1,nsurf(ix)-1
            c1(ix,iz)=c1(ix,nsurf(ix))
          enddo
        enddo
      endif
c
c  density
c
      if(ndens.eq.1) then
        do ix=1,nx
          do iz=nsurf(ix)+1,nz
            if(is_marine.eq.0.or.(is_marine.eq.1.and.c1(ix,iz)
     1                 .gt.1500.0))then
              den(ix,iz)=denco*(c1(ix,iz)/cp00)**denp
            else
              den(ix,iz)=1027.0
            endif

c Chaiwoot: For testing
c            den(ix,iz) = 1.0

            cl(ix,iz)=dtx/den(ix,iz)
            cl41(ix,iz)=cl(ix,iz)*c41
            cl42(ix,iz)=cl(ix,iz)*c42
          enddo
          do iz=1,nsurf(ix)
            den(ix,iz)=den(ix,nsurf(ix)+1)
            cl(ix,iz)=cl(ix,nsurf(ix)+1)
            cl41(ix,iz)=cl(ix,iz)*c41
            cl42(ix,iz)=cl(ix,iz)*c42
          enddo
        enddo
      else
        do ix=nxm1,nxm2
          ixx=ix-nxm1+ix1_shot(ishot)
          if(ixx.lt.1)then
            ixx=1
          else if(ixx.gt.nx_tot)then
            ixx=nx_tot
          endif
          do iz=nzm1,nzm2
            den(ix,iz)=density(ixx,iz-nzm1+1)
            write(*,*) den(ix,iz)
          enddo
        enddo
        do ix=nxm1,nxm2
          do iz=nzm2+1,nz
            den(ix,iz)=den(ix,nzm2)
          enddo
          do iz=nzm1-1,1,-1
            den(ix,iz)=den(ix,nzm1)
          enddo
        enddo
        do iz=1,nz
          do ix=nxm1-1,1,-1
            den(ix,iz)=den(nxm1,iz)
          enddo
          do ix=nxm2+1,nx
            den(ix,iz)=den(nxm2,iz)
          enddo
        enddo
        do ix=1,nx
          do iz=1,nsurf(ix)
            den(ix,iz)=den(ix,nsurf(ix)+1)
          enddo
        enddo
        do ix=1,nx
          do iz=1,nz
            cl(ix,iz)=dtx/den(ix,iz)
            cl41(ix,iz)=cl(ix,iz)*c41
            cl42(ix,iz)=cl(ix,iz)*c42
          enddo
        enddo
      endif
c
      do iz=1,nz
        do ix=1,nx
          c2(ix,iz)= (c1(ix,iz)**2)*den(ix,iz)*dtx
        enddo
      enddo
c
c  exL
c
      factor_att=0.1
      const=2.302585*dtx
      alpha_max0=sqrt(3.0)*vmax*(8.0/15.0-0.03*na_pml
     1        +na_pml*na_pml/1500.0)*dtx
      do ix=1,na_pml
        tmp=(na_pml+1-ix)/float(na_pml+1)
        tmp=tmp*tmp
        do iz=1,nz
          alpha=alpha_max0*tmp
          exL(ix,iz)=1.0-alpha
        enddo
      enddo
c
c  exR
c
      do ix=nxm20+1,nx
        ixx=nx-ix+1
        tmp=(na_pml+1-ixx)/float(na_pml+1)
        tmp=tmp*tmp
        do iz=1,nz
          alpha=alpha_max0*tmp
          exR(ixx,iz)=1.0-alpha
        enddo
      enddo
c
c  eyT
c
      do iz=1,na_pml
        tmp0=(na_pml+1.0-iz)/float(na_pml+1)
        tmp0=tmp0*tmp0
        do ix=1,nx
          alpha=alpha_max0*tmp0
          eyT(ix,iz)=1.0-alpha
        enddo
      enddo
c
c  eyB
c
      do iz=nzm20+1,nz
        izz=nz-iz+1
        tmp0=(na_pml+1.0-izz)/float(na_pml+1)
        tmp0=tmp0*tmp0
        do ix=1,nx
          alpha=alpha_max0*tmp0
          eyB(ix,izz)=1.0-alpha
        enddo
      enddo

999   continue
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine d2pre

      implicit none
      include 'wtw_model.para'
      integer ix, iz, n_constraint1, i, j
      real xx1

c  precondition
c
      if(isuseprecondition.eq.1)then
        do ix=1,nx_tot
          do iz=nsurf_tot(ix),nz_tot
            d2(ix,iz)=d2(ix,iz)*wafic(ix,iz)
          enddo
        enddo
      endif
      if(nfree.ne.0)then
        do ix=1,nx_tot
          d2(ix,nsurf_tot(ix))=d2(ix,nsurf(ix)+1)
        enddo
      endif

c  de/dv change to de/ds
c
      do ix=1,nx_tot
        do iz=nsurf_tot(ix),nz_tot
          xx1=1000.0/velocity(ix,iz)
          if(ispre.ne.1)then
            d2(ix,iz)=-d2(ix,iz)*xx1
          else
            d2(ix,iz)=d2(ix,iz)*xx1
          endif
        enddo
      enddo
c
c free surface boundary constraints
c
      do ix=1,nx_tot
        do iz=1,nsurf_tot(ix)-1
          d2(ix,iz)=0.0
        enddo
      enddo
c
c constraint: shallow layer
c
      if(is_use_constraint1.eq.1)then
        n_constraint1=nint(depth_constraint1/dx)
        do i=1,nx_tot
          do j=1,nsurf_tot(i)+n_constraint1
            if(abs(vmin_constraint1-velocity(i,j)).lt.0.1.and.
     1           d2(i,j).lt.0.0)then
              d2(i,j)=0.0
            endif
            if(abs(vmax_constraint1-velocity(i,j)).lt.0.1.and.
     1           d2(i,j).gt.0.0)then
              d2(i,j)=0.0
            endif
          enddo
        enddo
      endif
      if(is_use_water_bottom.eq.1)then
        do i=1,nx_tot
          do j=1,nwater(i)
            if(velocity(i,j).le.1490.0.and.
     1           d2(i,j).lt.0.0)then
              d2(i,j)=0.0
            endif
            if(velocity(i,j).ge.1500.0.and.
     1           d2(i,j).gt.0.0)then
              d2(i,j)=0.0
            endif
          enddo
        enddo
      endif

      end

