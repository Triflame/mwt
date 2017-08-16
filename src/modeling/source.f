      subroutine source
!     Creaing a Ricker source wavelet 

      implicit none
      include 'wtw_model.para'
      integer it, is1
      character*200 file_tmp
      real tshift,pi2,smax,b,const,tim1,tim2,u,amp,
     >     tdump(nt_max2),tdump2(nt_max2)

      tshift = sqrt(2.0)/vm
c      tshift = 5.0/vm
c      tshift = 3.0/vm
c      timeshift = nint(tshift/dt)
      timeshift = 0
      pi2 = sqrt(pi)/2.0
      b = sqrt(6.0)/(pi*vm)
      const = 2.0*sqrt(6.0)/b

      smax = 0.0
      do it=1,nt_work
        tim1 = real(it-1)*dt
        tim2 = tim1 - tshift
        u = const*tim2
        amp = ((u*u)/4.0-0.5)*pi2*(exp(-u*u/4.0))
        sou(1,it) = -amp
        if (smax .lt. abs(amp)) then
          smax = abs(amp)
        endif
      enddo

      smax = smax*2
      do it=1,nt_work
        sou(1,it) = sou(1,it)/smax
      enddo

c      do is1=is_first+1,is_end
c        do it=1,nt_work
c          sou(is1-is_first,it) = sou(1,it)
c        enddo
c      enddo
c      write(*,*) 'nt_work = ', nt_work

        end
!----------------------------------------------------------------------
      subroutine source_new
! Ricker is the 2nd derivative of a gaussian
!  f:   float array to contain the wavelet
!  n :  number of data points returned
!  fpeak: peak frequency (khz). --> sqrt(2)/tau for Ricker
!  dt:  time increment (ms)
! sign: initial polarity
! significant energy to twice fpeak
! peak frequency = sqrt(2)/tau
! the total length of f[i] should be 2*nts+1

      implicit none
      include 'wtw_model.para'
      integer it, is1, is_shift, sign1, i
      character*200 file_tmp
      real tshift,pi2,smax,a,b,const,tim1,tim2,u,amp,t,t2
      real sqrt
      intrinsic sqrt

!      timeshift = 44
      timeshift = 0
      nts = int((sqrt(2.0)/vm)/dt+0.5)
      b = (pi*vm*1000)
      b = b*b
      a = 2.*b
      sign1 = 1

      do i=1,nts+1
       t = (i-1)*dt*0.001
       t2 = t*t
       sou(1,nts+i) = sign1*(1.-a*t2)*exp(-b*t2)
       sou(1,nts+2-i) = sou(1,nts+i)
      enddo
!      n = 2*nts+1;

c      write(*,*) 'nt_work = ', nt_work
c      open(10,file='src.dat',access='direct',form='unformatted',
c     >  recl=nt_work)
c      write(10,rec=1) (sou(1,it),it=1,nt_work)
c      close(10)

      end

