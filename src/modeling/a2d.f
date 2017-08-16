      program ewt_cdp

      implicit none
      include 'mpif.h'
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
      nrepeat_max = 5
      nresd_repeat = 0
      resd_perc = 5.0

      if(my_rank.eq.0) then
        call get_input_parameter
        call get_geometry
        call get_velocity
      endif

      call MPI_BCAST(data_int,27,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(data_real,24,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(source_in,200,MPI_CHARACTER,0,MPI_COMM_WORLD
     1          ,ierr)
      call MPI_BCAST(syn_ou_pre,200,MPI_CHARACTER,0,MPI_COMM_WORLD
     1          ,ierr)
      call MPI_BCAST(syn_ou_suffix,200,MPI_CHARACTER,0,MPI_COMM_WORLD
     1          ,ierr)
      call MPI_BCAST(csg_pre,200,MPI_CHARACTER,0,MPI_COMM_WORLD
     1          ,ierr)
      call MPI_BCAST(csg_suffix,200,MPI_CHARACTER,0,MPI_COMM_WORLD
     1          ,ierr)

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

            if(is_marine.eq.0) then
              time_air=offset(k,is)/airvel
              nair=nint(time_air/dt)+1+ntpad
              itmmm(k,is00)=min(nair,itmmm(k,is00))
            endif

c Chaiwoot: Added only for testing: use this will disable end-trace muting
c            itmmm(k,is00)=nt_max
c            itmmm(k,is00)=1500

          endif
          if (itmmm(k,is00).gt.0.0) then
c            itmmm(k,is00) = 5000
            itmmm(k,is00) = nt_max
          endif

c          nt_work=max(nt_work,itmmm(k,is00)+4*ntpsou)
          nt_work=max(nt_work,itmmm(k,is00))
        enddo
      enddo
      if(abs(dt1-dt).lt.1.0e-6) then
        nt_work=min(nt,nt_work+100)
      else
        nt_work=min(nt_max,nt_work+100)
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

999   continue
      call MPI_FINALIZE(ierr)

      end

