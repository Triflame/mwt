      program ewt_cdp

      implicit none
      include './mpif.h'
      include 'global.h'

      integer nump,my_rank,ierr,status(MPI_STATUS_SIZE)
      include 'wtw_model.para'
      real gold, glimit, tiny, cgold, zero
      parameter(gold=1.618034,glimit=100.0,tiny=1.0e-20)
      parameter(cgold=0.3819660,zero=1.0e-10)
      character*200 file_tmp
      integer dump3(ns0*ng0), it, iss, iijj, ntsou_old, ntsou,
     1  ntsou_new, iii_nite, ix, iz, ig, kk0, nstop, j, k, kk,
     2  iii, ii, nn, myrank, nair
      real ccc(ns0), tt, rxx, sdis, szz, rdis, aamax, tt0, xamp1,
     1  xamp2, c41, c42, resddd, sss, ss2, resd1, resd2, smax, smin,
     2  gmax, gmin, gmean, ssg, ssmin, ssmax, ax, bx, fa, aaa, fb0,
     3  step_ini, fb, fc0, ux, ff, alp, s1, s2, tmp, gg, dgg, gam,
     4  fu0, per, fu, fc, cx, rzz, sxx, xx, ss, time_air, 
     5  slowmax, slowmin, dt_sou
      real al(2),aa(2)
      real dump1(nx_mod*nz_mod),dump11(nx_mod*nz_mod),dump2(ns0*ng0)
      real dump4(ng0,nt0),dump5(nt_max2),tdump(nt_max2)

! For Residual
      integer nresd_repeat, nrepeat_max
      real resd0, resd_perc, resd_diff

      real fun_step
      external fun_step

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nump,ierr)
      myrank = my_rank

! Residual is considered repeated if it decreases less than 5%
      nrepeat_max = 10
      nresd_repeat = 0
      resd_perc = 1.0

      if(my_rank.eq.0) then
        call get_input_parameter
        call get_geometry
        call get_velocity
      endif

      call MPI_BCAST(data_int,27,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(data_real,24,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(source_in,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(syn_ou_pre,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(syn_ou_suffix,200,MPI_CHARACTER,0,MPI_COMM_WORLD
     1          ,ierr)
      call MPI_BCAST(csg_pre,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(csg_suffix,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      nx_in=data_int(1)
      nz_in=data_int(2)
      dx1=data_real(1)
      nt=data_int(3)
      dt1=data_real(2)
      is_su_csg=data_int(4)
      dx=data_real(3)
      dt=data_real(4)
      is_marine=data_int(5)
      is_reset_surface=data_int(6)
      is_reset_sg=data_int(7)
      depsou=data_real(5)
      deprec=data_real(6)
      nopsou=data_int(8)
      vm=data_real(7)
      ntpsou=data_int(9)
      ndens=data_int(10)
      denco=data_real(8)
      cp00=data_real(9)
      denp=data_real(10)
      nfree=data_int(11)
      nit=data_int(12)
      iter00=data_int(13)
      nite=data_int(14)
      is_use_constraint1=data_int(15)
      depth_constraint1=data_real(11)
      vmin_constraint1=data_real(12)
      vmax_constraint1=data_real(13)
      dvmin=data_real(14)
      dvmax=data_real(15)
      vmin=data_real(16)
      vmax=data_real(17)
      xmin=data_real(18)
      xmax=data_real(19)
      zmin=data_real(20)
      zmax=data_real(21)
      offmin=data_real(22)
      offmax=data_real(23)
      dt_sou=data_real(24)
      ispre=data_int(16)
      isusenorm=data_int(17)
      isuseconjugate=data_int(18)
      isuseprecondition=data_int(19)
      ishot_out=data_int(20)
      nx_tot=data_int(21)
      nz_tot=data_int(22)
      ns=data_int(23)
      isepow=data_int(24)
      is_use_water_bottom=data_int(25)
      ntpad=data_int(26)
      nshot_int=data_int(27)

      if (my_rank.eq.0) then
        write(*,*) 'xmin = ', xmin, ', xmax = ', xmax
        write(*,*) 'offmin = ', offmin, ', offmax = ', offmax
      endif

      dtx=dt/dx
      dxm=dx
      nn=1
      nts=nint(1./vm/dt)
      na=na_store
      k0=na_store

      if(nx_tot.gt.nx_mod) then
        write(*,*)'nx_mod(',nx_mod,')<nx_tot(',nx_tot,')'
        stop
      endif

      if(nz_tot.gt.nz_mod) then
        write(*,*)'nz_mod(',nz_mod,')<nz_tot(',nz_tot,')'
        stop
      endif

      nzm1=na_new+1
      nzm2=nzm1+nz_tot-1
      nzm10=nzm1-na_pad
      nzm20=nzm2+na_pad
      nz=nz_tot+2*na_new

      ns_work=int((ns-1)/nshot_int)+1
      call get_assigned(is_first,is_end,1,ns_work,nump,my_rank)

      call MPI_BCAST( indshot,ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( ix2_shot,ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( ix1_shot,ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( xss,ns,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( zss,ns,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( nxss,ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( nzss,ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( ngis,ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( nsurf_tot,nx_tot,MPI_INTEGER,0
     1             ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( nwater,nx_tot,MPI_INTEGER,0
     1             ,MPI_COMM_WORLD,ierr)

      call MPI_BCAST( xgg,ns0*ng0,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( zgg,ns0*ng0,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( tr,ns0*ng0,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( offset,ns0*ng0,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( nxg,ns0*ng0,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( nzg,ns0*ng0,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( velocity,nx_mod*nz_mod,MPI_REAL,0,
     1              MPI_COMM_WORLD,ierr)

      slowmin=1.0/vmax
      slowmax=1.0/vmin
      do ix=1,nx_tot
        do iz=1,nz_tot
          velocity(ix,iz)=min(vmax,velocity(ix,iz))
          velocity(ix,iz)=max(vmin,velocity(ix,iz))
        enddo
      enddo
      call pad_model
      if(ndens.eq.0) then
        call MPI_BCAST( density,nx_mod*nz_mod,MPI_REAL,0,
     1              MPI_COMM_WORLD,ierr)
      endif

      nt_work=0

c is_first = the starting shot number for a computing node
c is_end   = the end shot number for a computing node
      do iss=is_first,is_end
        is=ns-(iss-1)*nshot_int
        is00=iss-is_first+1
        do k=1,ngis(is)

          itmmm(k,is00)=0

c Chaiwoot: this condition is quite important!!!
!          if(offset(k,is).gt.-0.5.and.tr(k,is).ge.0.0) then
          if(offset(k,is).lt.offmax.and.offset(k,is).gt.offmin
     >         .and.tr(k,is).ge.0.0) then
            itmmm(k,is00)=min(nt_max,nint(tr(k,is)/dt)+1+
     >        ntpsou+ntpad)

c            if(is_marine.eq.0) then
              time_air=offset(k,is)/airvel
              nair=nint(time_air/dt)+1+ntpad
              itmmm(k,is00)=min(nair,itmmm(k,is00))
c            endif

c Chaiwoot: Added only for testing: use this will disable end-trace muting
c            itmmm(k,is00)=nt_max
c            itmmm(k,is00)=1500

          endif
          if (itmmm(k,is00).gt.0.0) then
c            itmmm(k,is00) = 3000
            itmmm(k,is00) = nt_max
          endif

c          nt_work=max(nt_work,itmmm(k,is00)+4*ntpsou)
          nt_work=max(nt_work,itmmm(k,is00))
        enddo
      enddo
      if(abs(dt1-dt).lt.1.0e-6) then
c        nt_work=min(nt,nt_work+100)
        nt_work=min(nt,nt_work)
      else
c        nt_work=min(nt_max,nt_work+100)
        nt_work=min(nt_max,nt_work)
      endif

      do iss=is_first,is_end
        is=ns-(iss-1)*nshot_int
        is00=iss-is_first+1
        ss=-1.0
        ii=0
        do ig=1,ngis(is)
          if(itmmm(ig,is00).gt.0.and.offset(ig,is).gt.0.0) then
            if(ii.eq.0) then
              ss=offset(ig,is)
              ii=ig
            else
              if(offset(ig,is).lt.ss) then
                ss=offset(ig,is)
                ii=ig
              endif
            endif
          endif
        enddo
        indgeop(is00)=ii
      enddo

c Read the source wavelet file if existed
      if(nopsou.eq.1.and.my_rank.eq.0) then
        do it=1,nt_work
          dump5(it)=0.0
        enddo

! Chaiwoot: Read binary source file
        open(66,file=source_in,access='direct',recl=i4*1)
        do it=1,nt_work
          read(66,rec=it,err=5020) sss
          dump5(it)=sss
          data_int(1)=it
        enddo

c        open(66,file=source_in,access='direct',recl=355)
c        read(66,rec=1) (dump5(it),it=1,355)
c        data_int(1)=355
c        do it=1,nt_work
c          read(66,*,err=5020,end=5020)xx,sss
c          dump5(it)=sss
c          data_int(1)=it
c        enddo
5020    continue
        close(66)
      endif

c Broadcast the source wavelet
      if(nopsou.eq.1) then
        call MPI_BCAST(data_int,1,MPI_INTEGER,0,
     1          MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dump5,data_int(1),MPI_REAL,0,
     1          MPI_COMM_WORLD,ierr)
        if(abs(dt_sou-dt).lt.0.0001) then
          do it=1,data_int(1)
            do iijj=is_first,is_end
              sou(iijj-is_first+1,it)=dump5(it)
            enddo
          enddo
          do it=1+data_int(1),nt_work
            do iijj=is_first,is_end
              sou(iijj-is_first+1,it)=0.0
            enddo
          enddo
        else
          ntsou_old=data_int(1)
          ntsou_new=nint((ntsou_old-1)*dt_sou/dt)+1
          ntsou_new=min(ntsou_new,nt_work)
          do it=1,ntsou_old
            tdump(it)=float(it-1)*dt_sou
          enddo
          do it=1,ntsou_new
            tt=float(it-1)*dt
            call get1D(sss,tt,dump5,tdump,ntsou_old)
            do iijj=is_first,is_end
              sou(iijj-is_first+1,it)=sss
            enddo
          enddo
          do it=ntsou_new+1,nt_work
            do iijj=is_first,is_end
              sou(iijj-is_first+1,it)=0.0
            enddo
          enddo
        endif
      else

c Generate a source wavelet using Ricker wavelet
c        call source2
        call source
c        call source_new
      endif

c Write some input parameters into a file
      if(my_rank.eq.0) then
        call get_message
        call flush(30)

        write(*,*) 'nt_work = ', nt_work
        open(10,file='src.dat',access='direct',form='unformatted',
     >    recl=i4*nt_work)
        write(10,rec=1) (sou(1,it),it=1,nt_work)
        close(10)
      endif

      nite=min(nite,is_end-is_first+1)
      iii_nite=(is_end-is_first+1)/nite

      do iz=1,nz_tot
        zzcor(iz)=float(iz-1)*dx
      enddo
      do ix=1,nx_tot
        xxcor(ix)=float(ix-1)*dx
      enddo

      if (isuseprecondition.eq.1) then
        do iz=1,nz_mod
          do ix=1,nx_mod
            g1(ix,iz)=0.0
            wafic(ix,iz)=0.0
          enddo
        enddo
        do iss=is_first,is_end
          is=ns-(iss-1)*nshot_int
          do ig=1,ngis(is)
            if(offset(ig,is).ge.-0.5.and.tr(ig,is).ge.0.0) then
              do iz=1,nz_tot
                szz=zss(is)-zzcor(iz)
                rzz=zgg(ig,is)-zzcor(iz)
                do ix=ix1_shot(is)+npad_shot,ix2_shot(is)-npad_shot
                  sxx=xss(is)-xxcor(ix)
                  rxx=xgg(ig,is)-xxcor(ix)
                  sdis=sqrt(sxx*sxx+szz*szz)
                  rdis=sqrt(rxx*rxx+rzz*rzz)
                  g1(ix,iz)=g1(ix,iz)+1.0/sqrt((sdis+dx)*(rdis+dx))
                enddo
              enddo
            endif
          enddo
        enddo
        call MPI_ALLREDUCE(g1,wafic,nx_mod*nz_mod,MPI_REAL,
     1                  MPI_SUM,MPI_COMM_WORLD,ierr)
        aamax=0.0
        do iz=1,nz_tot
          do ix=1,nx_tot
            aamax=max(aamax,wafic(ix,iz))
          enddo
        enddo
        do iz=1,nz_tot
          do ix=1,nx_tot
            if(wafic(ix,iz).gt.0.0) then
              wafic(ix,iz)=aamax/wafic(ix,iz)
            endif
          enddo
        enddo
      endif

      if(my_rank.eq.0.and.isuseprecondition.eq.1) then
        call suwt11(wafic,'WAFIC.DAT',0.0,0.0,dx,dx,nx_tot,
     1          nz_tot,nx_mod,nz_mod)
      endif

      do it=1,nt
        tdump(it)=float(it-1)*dt1
      enddo

c Read CSG data files
      do iss=is_first,is_end
        is=ns-(iss-1)*nshot_int
        call filename(file_tmp,csg_pre,indshot(is),csg_suffix)

        if(is_su_csg.eq.1) then
c         Read CSG file in SU format
          call surd11(dump4,file_tmp,ng0,nt0,ngis(is),nt,iii)
        else
c         Read CSG file in binary format
          call bird3(dump4,file_tmp,1,ngis(is),1,nt,ng0,nt0,iii)
        endif

        if(iii.lt.ngis(is))then
          write(*,*)'Something wrong when reading file:',file_tmp
          write(*,*)'Only ',iii,' out of ',ngis(is),' were read.'
          stop
        endif

        is00=iss-is_first+1
        if(abs(dt1-dt).lt.1.0e-6)then
          do ig=1,ngis(is)
            do it=1,nt_work
              wr(it,ig,is00)=dump4(ig,it)
            enddo

c Chaiwoot: take into account the wavelet shift
c            do it=1,timeshift
c              wr(it,ig,is00) = 0.0
c            enddo
c            do it=1,nt_work-timeshift
c              wr(it+timeshift,ig,is00)=dump4(ig,it)
c            enddo
          enddo
        else
          do ig=1,ngis(is)
            do it=1,nt
              dump5(it)=dump4(ig,it)
            enddo
            do it=1,nt_work
              tt=float(it-1)*dt
              call get1D(sss,tt,dump5,tdump,nt)
              wr(it,ig,is00)=sss
            enddo

c Chaiwoot: take into account the wavelet shift
c            do it=1,timeshift
c              wr(it,ig,is00) = 0.0
c            enddo
c            do it=1,nt_work-timeshift
c              tt=float(it-1)*dt
c              call get1D(sss,tt,dump5,tdump,nt)
c              wr(it+timeshift,ig,is00)=sss
c            enddo
          enddo
        endif
      enddo

      if(isusenorm.eq.0.or.isusenorm.eq.1.or.isusenorm.eq.2) then
        do iss=is_first,is_end
          is=ns-(iss-1)*nshot_int
          is00=iss-is_first+1

          if (1.eq.0) then
          do ig=1,ngis(is)
            do it=1,nt_work
              wc(it,ig) = wr(it,ig,is00)
            enddo
          enddo
          call mute(wc, nt_max, ngis(is))
          do ig=1,ngis(is)
            do it=1,nt_work
              wr(it,ig,is00) = wc(it,ig)
            enddo
          enddo
          endif
          sss=0.0
          do ig=1,ngis(is)
            if(itmmm(ig,is00).gt.0)then
c Chaiwoot: testing
c              kk=itmmm(ig,is00)
              kk=nt_work

              tt0=min(tr(ig,is)*0.9,tr(ig,is)-2.0/vm)
              kk0=int(tt0/dt)+1100

c Chaiwoot: Added for testing: use this when there is no first-break pick.
              kk0 = 0
c              kk0 = timeshift

c Mute traces before the direct wave arrivals
              do it=1,kk0
                wr(it,ig,is00)=0.0
              enddo
c Mute traces after early arrivals
              do it=kk+1,nt_work
                wr(it,ig,is00)=0.0
              enddo

c Find the maximum value of the trace for trace normalization
              do it=1,kk
                sss=max(sss,abs(wr(it,ig,is00)))
              enddo
            else
c Mute entire traces
              do it=1,nt_work
                wr(it,ig,is00)=0.0
              enddo
            endif
          enddo
          sss=1.0/sss

c Normalize traces
          do ig=1,ngis(is)
            if(itmmm(ig,is00).gt.0)then
c              do it=1,itmmm(ig,is00)
              do it=1,nt_work
                wr(it,ig,is00)=wr(it,ig,is00)*sss
              enddo
            endif
          enddo

          do ig=1,ngis(is)
c           Find zero traces
            if(itmmm(ig,is00).gt.0)then
              xamp1=wr(1,ig,is00)
              xamp2=wr(1,ig,is00)
              do it=2,itmmm(ig,is00)
                xamp1=min(xamp1,wr(it,ig,is00))
                xamp2=max(xamp2,wr(it,ig,is00))
              enddo
              if(xamp2-xamp1.lt.1.0e-10) itmmm(ig,is00)=0
            endif

            if(itmmm(ig,is00).gt.0)then
              wweight(ig,is00)=1.0
              if(isusenorm.eq.1)then
                ss2=max(abs(xamp1),abs(xamp2))
                ss2=1.0/ss2
                do it=1,itmmm(ig,is00)
                  wr(it,ig,is00)=wr(it,ig,is00)*ss2
                enddo
              endif
            endif
          enddo

          if (1.eq.0) then
          do ig=1,ngis(is)
            do it=1,nt_work
              wc(it,ig) = wr(it,ig,is00)
            enddo
          enddo
          call mute(wc, nt_max, ngis(is))
          do ig=1,ngis(is)
            do it=1,nt_work
              wr(it,ig,is00) = wc(it,ig)
            enddo
          enddo
          endif

          if (indshot(is).eq.ishot_out) then
            open(10,file='shot1.bin',access='direct',recl=nt_work*i4,
     >        status='replace')
            do ig=1,ngis(1)
              write(10,rec=ig) (wr(it,ig,is00),it=1,nt_work)
            enddo
            close(10)
!            call suwt22(wr(1,1,is00),'csg1.su',1.0,-1.0*ntpad*dt,1.0,
!     >         dt,ngis(is), nt_work,ng0,nt_max)
          endif
        enddo
      endif

c 4-order FD coefficients
      c41=9./8.
      c42=-1./24.

      do iz=1,nz_mod
        do ix=1,nx_mod
          g1(ix,iz)=0.
          d2(ix,iz)=0.
        enddo
      enddo
      pur=0.

      iter=iter00-1
      call forw_back(c41,c42,iii_nite,resd1,resd2)
      resddd=0.0
      resd0 = 1.0e10
      aaa=0.0

c Sum all the gradient from all CPU

      call MPI_ALLREDUCE(g1,d2,nx_mod*nz_mod,MPI_REAL,
     1                  MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(resd1,resddd,1,MPI_REAL,
     1                  MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(resd2,aaa,1,MPI_REAL,
     1                  MPI_SUM,MPI_COMM_WORLD,ierr)
      call d2pre

c Chaiwoot: for testing
c        do iz=1,4
c          do ix=1,nx_mod
c            d2(ix,iz)=0.0
c          enddo
c        enddo

      if(my_rank.eq.0)then
        write(30,*)'iter=',iter00-1,' resddd=',resddd

        resd_diff = abs(resddd-resd0)*100.0/resd0
        if (resd_diff .lt. resd_perc) then
          nresd_repeat = nresd_repeat + 1
        endif
        resd0 = resddd
        write(*,*)'iter=',iter00-1,' resddd=',resddd,
     >    ' diff=',resd_diff, nresd_repeat

c       Output the gradient
        call suwt11(d2,grad_ou,0.0,0.0,dx,dx,nx_tot,
     1          nz_tot,nx_mod,nz_mod)
      endif

c Chaiwoot: mute the gradient above the free surface
      do ix=1,nx_tot
        do iz=1,nsurf_tot(ix)-1
            d2(ix,iz) = 0.0
        enddo
      enddo

      do iz=1,nz_tot
        do ix=1,nx_tot
          g1tmp(ix,iz)=d2(ix,iz)
          g1old(ix,iz)=g1tmp(ix,iz)
        enddo
      enddo

ccccccccccccccccccccccc
c
c MAIN INVERSION LOOP
c
ccccccccccccccccccccccc
      step_ini=5.0
      do 800 iter=iter00,nit

c  search for a step length
c
        nstop=0
        do ix=1,nx_tot
          do iz=nsurf_tot(ix),nz_tot
            ss=max(velocity(ix,iz),vmin)
            ss=min(velocity(ix,iz),vmax)
            vdump(ix,iz)=1.0/ss
          enddo
        enddo

        if(is_use_water_bottom.eq.1) then
          gmax=0.0
          ssmax=0.0
          ssmin=1.0
          gmean=0.0
          ssg=0.0
          do j=1,nx_tot
            do k=nwater(j),nz_tot
              gmax=max(gmax,abs(d2(j,k)))
              ssmax=max(ssmax,vdump(j,k))
              ssmin=min(ssmin,vdump(j,k))
              gmean=gmean+abs(d2(j,k))
              ssg=ssg+1.0
            enddo
          enddo
        else
          gmax=0.0
          gmean=0.0
          ssg=0.0
          ssmax=0.0
          ssmin=1.0
          do j=1,nx_tot
            do k=nsurf_tot(j),nz_tot
              gmax=max(gmax,abs(d2(j,k)))
              ssmax=max(ssmax,vdump(j,k))
              ssmin=min(ssmin,vdump(j,k))
              gmean=gmean+abs(d2(j,k))
              ssg=ssg+1.0
            enddo
          enddo
        endif
        gmean=gmean/ssg*2.0
c        pur=ssmax/100.0/gmax
        pur=ssmax/100.0/gmean
        nstop=0
  
        ax=0.0
        bx=5.0
        fa=aaa
        fb0=fun_step(pur*bx,c41,c42,iii_nite)
        fb=0.0
        call MPI_ALLREDUCE(fb0,fb,1,MPI_REAL,
     1                  MPI_SUM,MPI_COMM_WORLD,ierr)
        if(my_rank.eq.0)then
          write(30,*)'ax,bx,fa,fb=',ax,bx,fa,fb
        endif

        if(fb.gt.fa)then
699       fc=fb
          cx=bx
          bx=bx/2.0
          fb0=fun_step(pur*bx,c41,c42,iii_nite)
          fb=0.0
          call MPI_ALLREDUCE(fb0,fb,1,MPI_REAL,
     1                  MPI_SUM,MPI_COMM_WORLD,ierr)
c          if(bx.lt.1.0e-3) stop
          if(fb.gt.fa.and.bx.gt.0.1)then
            go to 699
          endif
        else
          cx=bx*2.0
          fc0=fun_step(pur*cx,c41,c42,iii_nite)
          fc=0.0
          call MPI_ALLREDUCE(fc0,fc,1,MPI_REAL,
     1                  MPI_SUM,MPI_COMM_WORLD,ierr)
        endif

        if(my_rank.eq.0)then
          write(30,*)'ax,bx,cx,fa,fb,fc=',ax,bx,cx,fa,fb,fc
        endif
c Quadratic interpolation model for line search
        ux=(fa*bx*bx-fc*bx*bx-fa*cx*cx+fb*cx*cx)/
     1                  (fa*bx-fc*bx-fa*cx+fb*cx)*0.5

        if(ux.lt.0.0)then
          ff=fb
          per=bx
          if(fb.gt.fc)then
            ff=fc
            per=cx
          endif
          go to 605
        endif

        if(ux.gt.20.0) ux=20.0
        fu0=fun_step(pur*ux,c41,c42,iii_nite)
        fu=0.0
        call MPI_ALLREDUCE(fu0,fu,1,MPI_REAL,
     1                  MPI_SUM,MPI_COMM_WORLD,ierr)
        if(my_rank.eq.0)then
          write(30,*)'uu,fun=',ux,fu
        endif
        ff=fb
        per=bx
        if(fb.gt.fc)then
          ff=fc
          per=cx
        endif
        if(fu.lt.ff)then
          ff=fu
          per=ux
        endif
605     continue
c
        if(my_rank.eq.0)then
          write(30,*)'per,ff=',per,ff
          call flush(30)
        endif
c Update velocity model
        alp=pur*per
        s1=1.0/vmin
        s2=1.0/vmax

c We should constraint the gradient at the top part of the model
c due to the high sensitivity
        do ix=1,nx_tot
          do iz=nsurf_tot(ix),nz_tot
            tmp=vdump(ix,iz)-d2(ix,iz)*alp
            tmp=min(tmp,s1)
            tmp=max(tmp,s2)
            velocity(ix,iz)=1.0/tmp
          enddo
        enddo
        call pad_model

        if(my_rank.eq.0)then
          call filename(file_tmp,vel_ou_pre,iter,vel_ou_suffix)
          call suwt11(velocity,file_tmp,0.0,0.0,dx,dx,nx_tot,
     1          nz_tot,nx_mod,nz_mod)
        endif

        do iz=1,nz_mod
          do ix=1,nx_mod
            g1(ix,iz)=0.0
            d2(ix,iz)=0.0
          enddo
        enddo

        call forw_back(c41,c42,
     1          iii_nite,resd1,resd2)
        call MPI_ALLREDUCE(g1,d2,nx_mod*nz_mod,MPI_REAL,
     1          MPI_SUM,MPI_COMM_WORLD,ierr)
        resddd=0.0

        call MPI_ALLREDUCE(resd1,resddd,1,MPI_REAL,
     1          MPI_SUM,MPI_COMM_WORLD,ierr)
        aaa=0.0
        call MPI_ALLREDUCE(resd2,aaa,1,MPI_REAL,
     1          MPI_SUM,MPI_COMM_WORLD,ierr)

c Chaiwoot: mute the gradient above the free surface
        do ix=1,nx_tot
          do iz=1,nsurf_tot(ix)-1
              d2(ix,iz) = 0.0
          enddo
        enddo

c Change the derivative of gradient to be repected to slowness
c instead of velocity
        call d2pre

        if(my_rank.eq.0)then
          write(30,*)'iter=',iter,' resddd=',resddd
          write(31,*) iter, resddd
c          write(*,*)'iter=',iter,' resddd=',resddd

          resd_diff = (resd0-resddd)*100.0/resd0
          if (abs(resd_diff) .lt. resd_perc) then
            nresd_repeat = nresd_repeat + 1
          else if (nresd_repeat .gt. 0) then
            nresd_repeat = 0
          endif
          resd0 = resddd
          write(*,*)'iter=',iter,' resddd=',resddd,
     >      ' diff=',resd_diff, nresd_repeat, ', window=',nt_work
          call flush(30)
        endif

        if (nresd_repeat .eq. nrepeat_max) goto 999

        if (nresd_repeat .eq. nrepeat_max .and. .false.) then
          nt_work = min(nt_max, nt_work+200)
          do iss=is_first,is_end
            is=ns-(iss-1)*nshot_int
            is00=iss-is_first+1
 
c Chaiwoot: Increasing window size when the residual is not changing significantly 
            do k=1,ngis(is)
              itmmm(k,is00)=min(nt_work,itmmm(k,is00)+200)
            enddo
          enddo
          nresd_repeat = 0
        endif

        if(isuseconjugate.eq.1.or.(isuseconjugate.eq.2.and.
     1         (iter-iter/3*3).ne.0))then
          dgg=0.0
          gg=0.0
          do ix=1,nx_tot
            do iz=nsurf_tot(ix),nz_tot
              gg=gg+g1tmp(ix,iz)*g1tmp(ix,iz)
c              dgg=dgg+d2(ix,iz)*d2(ix,iz)
              dgg=dgg+(d2(ix,iz)-g1tmp(ix,iz))*d2(ix,iz)
            enddo
          enddo
          dgg=max(dgg,0.0)
          if(gg.le.1.0e-16) stop
          gam=dgg/gg
          do ix=1,nx_tot
            do iz=nsurf_tot(ix),nz_tot
              g1tmp(ix,iz)=d2(ix,iz)
              g1old(ix,iz)=g1tmp(ix,iz)+gam*g1old(ix,iz)
              d2(ix,iz)=g1old(ix,iz)
            enddo
          enddo
        else if(isuseconjugate.eq.2.and.(iter-iter/3*3).eq.0)then
          do ix=1,nx_tot
            do iz=nsurf_tot(ix),nz_tot
              g1tmp(ix,iz)=d2(ix,iz)
              g1old(ix,iz)=g1tmp(ix,iz)
            enddo
          enddo
        endif

c Chaiwoot: for testing
c        do iz=1,4
c          do ix=1,nx_mod
c            d2(ix,iz)=0.0
c          enddo
c        enddo

c Output the gradient d2
        do ix=1,nx_tot
          do iz=1,nsurf_tot(ix)-1
              d2(ix,iz) = 0.0
          enddo
        enddo
        if(my_rank.eq.0)then
c          open(20,file='surf.txt',form='formatted')
c          do ix=1,nx_tot
c            write(20,*) nsurf_tot(ix)
c          enddo
c          close(20)
          call suwt11(d2,grad_ou,0.0,0.0,dx,dx,nx_tot,
     1          nz_tot,nx_mod,nz_mod)
        endif

c        goto 999
800   continue

      if(my_rank.eq.0)then
        close(30)
        close(31)
      endif

999   continue
      call MPI_FINALIZE(ierr)

      end

