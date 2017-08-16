      subroutine init_segy(
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

      integer*2 unass(17)

      tracl=0
      tracr=0
      fldr=0
      tracf=0
      ep=0
      cdp=0
      cdpt=0
      trid=0
      nvs=0
      nhs=0
      duse=0
      offset=0
      gelev=0
      selev=0
      sdepth=0
      gdel=0
      sdel=0
      swdep=0
      gwdep=0
      scalel=0
      scalco=0
      sx=0
      sy=0
      gx=0
      gy=0
      counit=0
      wevel=0
      swevel=0
      sut=0
      gut=0
      sstat=0
      gstat=0
      tstat=0
      laga=0
      lagb=0
      delrt=0
      muts=0
      mute=0
      ns=0
      dt=0
      gain=0
      igc=0
      igi=0
      corr=0
      sfs=0
      sfe=0
      slen=0
      styp=0
      stas=0
      stae=0
      tatyp=0
      afilf=0
      afils=0
      nofilf=0
      nofils=0
      lcf=0
      hcf=0
      lcs=0
      hcs=0
      year=0
      day=0
      hour=0
      minute=0
      sec=0
      timbas=0
      trwf=0
      grnors=0
      grnofr=0
      grnlof=0
      gaps=0
      otrav=0
      d1=0
      f1=0
      d2=0
      f2=0
      ungpow=0
      unscale=0
      mark=0
      do i=1,17
        unass(i) = 0
      enddo

      end

