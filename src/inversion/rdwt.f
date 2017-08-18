      subroutine bird(trace,inname,nx,nz,iii)

      implicit none
      include 'global.h'
      integer nx,nz,iii, j, itrace
      real trace(nx,nz)
      character*(*) inname

      open(41,file=inname,access='direct',recl=nz*i4)
      iii=0
      do itrace = 1, nx
        read(41,rec=itrace,err=111) (trace(itrace,j),j=1,nz)
        iii=iii+1
      enddo
111   continue
      close(41)
      end

!----------------------------------------------------------------------
      subroutine bird3(trace,inname,nxm1,nxm2,nzm1,nzm2,mx,mz,iii)

      implicit none
      include 'global.h'
      integer nxm1,nxm2,nzm1,nzm2,mx,mz,iii, j, itrace, nlen
      real trace(mx,mz)
      character*(*) inname

      nlen = nzm2-nzm1+1
      open(41,file=inname,access='direct',recl=nlen*i4)
      iii=0
      do itrace = nxm1, nxm2
        read(41,rec=itrace,err=111)
     1   (trace(itrace+nxm1-1,j),j=nzm1,nzm2)
        iii=iii+1
      enddo
c      do itrace = 1, nxm2-nxm1+1
c        read(41,rec=itrace,err=111)
c     1   (trace(itrace+nxm1-1,j),j=nzm1,nzm2)
c        iii=iii+1
c      enddo
111   continue
      close(41)
      end

!----------------------------------------------------------------------
      subroutine biwt3(trace,inname,nxm1,nxm2,nzm1,nzm2,mx,mz)

      implicit none
      include 'global.h'
      integer nxm1,nxm2,nzm1,nzm2,mx,mz, nlen, itrace, j
      real trace(mx,mz)
      character*(*) inname

      nlen = nzm2-nzm1+1
      open(40,file=inname,access='direct',recl=nlen*i4)
      close(40,status='delete')
      open(40,file=inname,access='direct',recl=nlen*i4)
      do itrace=1,nxm2-nxm1+1
        write(40,rec=itrace)
     1   (trace(itrace+nxm1-1,j),j=nzm1,nzm2)
      enddo
      close(40)
      end

!----------------------------------------------------------------------
      subroutine surd11(trace,inname,mx,mz,mx1,mz1,iii)

      implicit none
      include 'global.h'
      include "./segy.h"
      integer mx1,mz1,mx,mz, nhead, nlen, i, j, k, itrace, iii
      real trace(mx,mz)
      character*(*) inname
c
c Include the SEGY-header format variables
c
      nhead = 60
      nlen = mz1 + nhead
      open(40,file=inname,access='direct',recl=nlen*i4)

      iii=0
      do 10 itrace = 1, mx1
        read(40,rec=itrace,err=111)
     1   tracl,tracr,fldr,tracf,ep,cdp,cdpt,
     1   trid,nvs,nhs,duse,offset,gelev,selev,sdepth,
     1   gdel,sdel,swdep,gwdep,scalel,scalco,
     1   sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat,
     1   tstat,laga,lagb,delrt,muts,mute,ns,dt,
     1   gain,igc,igi,corr,sfs,sfe,slen,styp,stas,
     1   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,
     1   lcs,hcs,year,day,hour,minute,sec,timbas,
     1   trwf,grnors,grnofr,grnlof,gaps,otrav,
     1   d1,f1,d2,f2,ungpow,unscale,mark,(unass(k),k=1,17),
     1   (trace(itrace,j),j=1,mz1)
        iii=iii+1
10    continue
111   continue
      close(40)
      end

!----------------------------------------------------------------------
      subroutine suwt11(trace,inname,xmin,zmin,dx,dz,mx1,mz1,mx,mz)

      implicit none
      include 'global.h'
      include "./segy.h"
      integer mx1,mz1,mx,mz, nhead, nlen, i, j, k, itrace
      real trace(mx,mz), xmin,zmin,dx,dz
      character*(*) inname
c
c Include the SEGY-header format variables
c
      nhead = 60
      nlen = mz1 + nhead
      open(40,file=inname,access='direct',recl=nlen*i4)
      close(40,status='delete')
      call init_segy(
     1   tracl,tracr,fldr,tracf,ep,cdp,cdpt,
     1   trid,nvs,nhs,duse,offset,gelev,selev,sdepth,
     1   gdel,sdel,swdep,gwdep,scalel,scalco,
     1   sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat,
     1   tstat,laga,lagb,delrt,muts,mute,ns,dt,
     1   gain,igc,igi,corr,sfs,sfe,slen,styp,stas,
     1   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,
     1   lcs,hcs,year,day,hour,minute,sec,timbas,
     1   trwf,grnors,grnofr,grnlof,gaps,otrav,
     1   d1,f1,d2,f2,ungpow,unscale,mark,unass)

      ns = mz1
      ep = 1
      fldr = 1
c********************************************************
c  sampling interval "dt"  should be manually adjusted
c  to get the right value
c
      dt = int(dz*1000*1000)
      d1 = dz
      d2 = dx
      f1 = (zmin)
      f2 = (xmin)
      tracl = 0
      tracr = 0

      open(40,file=inname,access='direct',recl=nlen*i4) ! bug is here
      do itrace = 1,mx1
        tracl = tracl + 1
        tracr = tracr + 1
         cdp = itrace
        cdpt = itrace

        write(40,rec=itrace)
     1   tracl,tracr,fldr,tracf,ep,cdp,cdpt,
     1   trid,nvs,nhs,duse,offset,gelev,selev,sdepth,
     1   gdel,sdel,swdep,gwdep,scalel,scalco,
     1   sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat,
     1   tstat,laga,lagb,delrt,muts,mute,ns,dt,
     1   gain,igc,igi,corr,sfs,sfe,slen,styp,stas,
     1   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,
     1   lcs,hcs,year,day,hour,minute,sec,timbas,
     1   trwf,grnors,grnofr,grnlof,gaps,otrav,
     1   d1,f1,d2,f2,ungpow,unscale,mark,(unass(k),k=1,17),
     1   (trace(itrace,j),j=1,ns)
      enddo
      close(40)
      end

!----------------------------------------------------------------------
      subroutine suwt22(trace,inname,xmin,zmin,dx,dz,mx1,mz1,mx,mz)

      include 'global.h'
      include "./segy.h"
      real trace(mz,mx)
      character*(*) inname

c Include the SEGY-header format variables
c
      nhead = 60
      nlen = mz1 + nhead
      open(40,file=inname,access='direct',recl=nlen*i4)
      close(40,status='delete')
      call init_segy(
     1   tracl,tracr,fldr,tracf,ep,cdp,cdpt,
     1   trid,nvs,nhs,duse,offset,gelev,selev,sdepth,
     1   gdel,sdel,swdep,gwdep,scalel,scalco,
     1   sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat,
     1   tstat,laga,lagb,delrt,muts,mute,ns,dt,
     1   gain,igc,igi,corr,sfs,sfe,slen,styp,stas,
     1   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,
     1   lcs,hcs,year,day,hour,minute,sec,timbas,
     1   trwf,grnors,grnofr,grnlof,gaps,otrav,
     1   d1,f1,d2,f2,ungpow,unscale,mark,unass)

c in linux cluster nlen becomes zero after here, so read data gives
c errors, so re-number nlen as above
c add by Tieyuan March 2017
      nhead = 60
      nlen = mz1 + nhead

      ns = mz1
      ep = 1
      fldr = 1
c********************************************************
c  sampling interval "dt"  should be manually adjusted
c  to get the right value
c
      dt = int(dz*1000*1000)
      d1 = dz
      d2 = dx
      f1 = (zmin)
      f2 = (xmin)
      tracl = 0
      tracr = 0
      
      open(40,file=inname,access='direct',recl=nlen*i4)
      do itrace = 1,mx1
        tracl = tracl + 1
        tracr = tracr + 1
        cdp = itrace
        cdpt = itrace

        write(40,rec=itrace)
     1   tracl,tracr,fldr,tracf,ep,cdp,cdpt,
     1   trid,nvs,nhs,duse,offset,gelev,selev,sdepth,
     1   gdel,sdel,swdep,gwdep,scalel,scalco,
     1   sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat,
     1   tstat,laga,lagb,delrt,muts,mute,ns,dt,
     1   gain,igc,igi,corr,sfs,sfe,slen,styp,stas,
     1   stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,
     1   lcs,hcs,year,day,hour,minute,sec,timbas,
     1   trwf,grnors,grnofr,grnlof,gaps,otrav,
     1   d1,f1,d2,f2,ungpow,unscale,mark,(unass(k),k=1,17),
     1   (trace(j,itrace),j=1,ns)
      enddo
      close(40)
      end

