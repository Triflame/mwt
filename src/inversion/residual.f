      subroutine norm12

      implicit none
      include 'wtw_model.para'
      integer ig, it, kk, kk0
      real ss, ss0, tt0
c     --------------------------------------------------------------
c     normalize seismograms
c     ---------------------------------------------------------------
c
      ss=0.0
      do ig=1,ngis(is)
        if(itmmm(ig,is00).gt.0)then

c         For muting data before first-arrival
          tt0 = min(tr(ig,is)*0.9,tr(ig,is)-2.0/vm)
          kk0 = int(tt0/dt)+1+1100
          kk0 = 0
          do it=1,kk0
            wc(it,ig)=0.0
          enddo

c         For muting calculated seismograms after early arrivals
          wweight(ig,is00)=1.0
          kk=itmmm(ig,is00)

          do it=kk+1,nt_work
            wc(it,ig)=0.0
          enddo

          if(isusenorm.eq.1)then
            ss=0.0
            do it=1,kk
              ss=max(ss,abs(wc(it,ig)))
            enddo
            ss0=ss*ss
            if(ss.lt.1.0e-16)then
              ss=0.0
            else
              ss=1.0/ss
            endif
            do it=1,kk
              wc(it,ig)=wc(it,ig)*ss
            enddo
          else
            do it=1,kk
              ss=max(ss,abs(wc(it,ig)))
            enddo
          endif
        else
          do it=1,nt_work
            wc(it,ig) = 0.0
          enddo
        endif
      enddo

      if(isusenorm.ne.1)then
        ss=1.0/ss
        do ig=1,ngis(is)
          if(itmmm(ig,is00).gt.0)then
            do it=1,itmmm(ig,is00)
              wc(it,ig)=wc(it,ig)*ss
            enddo
          endif
        enddo
      endif

c Write muted computed CSG
c      call suwt22(wc,'mute.su',1.0,-1.0*ntpad*dt,1.0,dt,ngis(is),
c     1         nt_work,ng0,nt_max)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine residual(resid,iii_nite,resd1,resd2)

      implicit none
      include 'wtw_model.para'
      integer iii_nite, k, it, nt_work0, iss, ig, ii2, ii
      real resid, resd1, resd2, tmp, tt, s1, s2, s3, s4
      integer mod2
      external mod2
      real sqrt
      intrinsic sqrt

c     ----------------------------------------------------------
c     Calculating pseudo-traveltime residuals
c     -----------------------------------------------------------
c

c      write(*,*) 'rank ', myrank, ', here', isusenorm
      if(isusenorm.eq.1)then
        nt_work0=2*nt_work-100
      else
        nt_work0=nt_work
      endif

      do k=1,ng0
        do it=1,nt_work0
          dw(it,k)=0.0
        enddo
      enddo
      resid=0.0
      tt=1.0

c     L-2 Norm Residual
      if(isusenorm.eq.2)then
        do ig=1,ngis(is)
          if(itmmm(ig,is00).gt.0)then
            s1=0.0
            s2=0.0
            s3=0.0
            do it=1,itmmm(ig,is00)
              if(isepow.eq.1)then
                tt=float(it-1)*dt
                tt=exp(-(tt/(tr(ig,is)+dt))**2.)
              endif
              s1=s1+wc(it,ig)*wc(it,ig)*tt*tt
              s2=s2+wr(it,ig,is00)*wc(it,ig)*tt*tt
              s3=s3+wr(it,ig,is00)*wr(it,ig,is00)*tt*tt
            enddo
            s1=sqrt(s1)
            s3=sqrt(s3)
            if(s1.gt.0.000001)then
              s1=1.0/(s1)
              s3=1.0/(s3)
              s4=s1*s1*s1
              do it=1,itmmm(ig,is00)
                if(isepow.eq.1)then
                  tt=float(it-1)*dt
                  tt=exp(-(tt/(tr(ig,is)+dt))**2.)
                endif
                dw(it,ig)=(s1*wr(it,ig,is00)*tt*tt-s2*s4
     1                    *wc(it,ig)*tt)*s3
              enddo
              tmp=s2*s1*s3
              resid=resid-tmp
            endif
          endif
        enddo
      else if(isusenorm.eq.0.or.isusenorm.eq.1)then
        do ig=1,ngis(is)
          if(itmmm(ig,is00).gt.0)then
            do it=1,itmmm(ig,is00)
              if(isepow.eq.1)then
                tt=float(it-1)*dt
                tt=exp(-(tt/(tr(ig,is)+dt))**2.)
              endif
              tmp=wr(it,ig,is00)-wc(it,ig)
              dw(it,ig)=tmp*tt
              resid=resid+tmp*tmp*tt
            enddo
          endif
        enddo
      endif

      resd1=resd1+resid
      do iss=iter-iter00,iter-iter00+iii_nite*(nite-1),iii_nite
        ii2=is_first+mod2(iss,(is_end-is_first+1))
        ii=ns-(ii2-1)*nshot_int
        if( ii.eq.is ) then
          resd2=resd2+resid
          go to 111
        endif
      enddo
111   continue
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine residual22(resid)

      include 'wtw_model.para'

      resid=0.0
      tt=1.0
      if(isusenorm.eq.2)then
        do ig=1,ngis(is)
          if(itmmm(ig,is00).gt.0)then
            s1=0.0
            s2=0.0
            s3=0.0
            do it=1,itmmm(ig,is00)
              if(isepow.eq.1)then
                tt=float(it-1)*dt
                tt=exp(-(tt/(tr(ig,is)+dt))**2.)
              endif
              s1=s1+wc(it,ig)*wc(it,ig)*tt*tt
              s2=s2+wr(it,ig,is00)*wc(it,ig)*tt*tt
              s3=s3+wr(it,ig,is00)*wr(it,ig,is00)*tt*tt
            enddo
            if(sqrt(s1).gt.0.000001)then
              s4=sqrt(s1)*sqrt(s3)
              resid=resid-s2/s4
            endif
          endif
        enddo
      else if(isusenorm.eq.0.or.isusenorm.eq.1)then
        do ig=1,ngis(is)
          if(itmmm(ig,is00).gt.0)then
            do it=1,itmmm(ig,is00)
              if(isepow.eq.1)then
                tt=float(it-1)*dt
                tt=exp(-(tt/(tr(ig,is)+dt))**2.)
              endif
              tmp=wr(it,ig,is00)-wc(it,ig)
              resid=resid+tmp*tmp*tt
            enddo
          endif
        enddo
      endif
      
      end

