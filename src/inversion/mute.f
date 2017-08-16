      subroutine mute(wc, nt, ng)

      implicit none
      integer ntaper
      parameter(ntaper=50)
      integer nt, ng, it, ig, iw, nsample, window, tmute
      real wc(nt,ng), trace(nt), dtrace(nt), trace_max, trace_mean(nt),
     >  previous, current, next, threshold, taper(ntaper), ss, pi

      pi = 3.1415926
      do it=1,ntaper
        ss = float(it-1)*pi/2.0/float(ntaper)
        taper(it)=1.0-cos(ss)*cos(ss)
      enddo

      window = 200
      threshold = 0.000001
      do ig=1,ng
        trace_max = 0.0
        do it=1,nt
          trace(it) = wc(it,ig)
          if (trace_max .lt. abs(trace(it))) trace_max = abs(trace(it))
        enddo
        if (trace_max .lt. 1.0e-10) goto 222

        do it=1,nt
          trace(it) = trace(it)/trace_max
        enddo

c       Calculate window mean
        trace_max = 0.0
        do it=1,nt
          nsample = 0
          trace_mean(it) = 0.0
          do iw=max(1,it-window/2),min(nt,it+window/2)
            nsample = nsample + 1
            trace_mean(it) = trace_mean(it) + trace(iw)
          enddo
          trace_mean(it) = trace_mean(it)/real(nsample)
          if (trace_max .lt. abs(trace_mean(it)))
     >      trace_max = abs(trace_mean(it))
        enddo
        do it=1,nt
          trace_mean(it) = trace_mean(it)/trace_max
        enddo

        tmute = 0
        do it=2,nt-1
          previous = trace_mean(it-1)
          current = trace_mean(it)
          next = trace_mean(it+1)
          if (current .gt. previous .and. current .gt. next .and. 
     >        current .gt. threshold) then
            tmute = it
            goto 111
          endif
        enddo
111     continue
c        tmute = tmute + 150
        tmute = tmute + 320
c       Mute trace
        tmute = 800
        do it=1,tmute
          wc(it,ig) = 0.0
        enddo
        do it=tmute+1,tmute+ntaper
          wc(it,ig) = wc(it,ig)*taper(it-tmute)
        enddo
222     continue
      enddo

      end

