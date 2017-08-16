c Read EWT output file and convert to residual file
c
c Chaiwoot Boonyasiriwat (July 17, 2007)

      program read_residual

      implicit none

      character(len=200) input, output, line
      integer iIter, iRes, iter
      real residual

      if (iarg() .lt. 2) then
        write(*,*) 'Usage: read_residual.x input output'
        goto 999
      endif
      call getarg(1, input)
      call getarg(2, output)

c      write(*,*) input,output
      open(10,file=input)
      open(20,file=output)
10    continue
        read(10,'(a)',end=20) line
c        write(*,*) line
        iIter = index(line,'iter=')
        iRes = index(line,'resddd=')
c        write(*,*) iIter, iRes
        if (iIter.eq.0 .or. iRes.eq.0) goto 10
c        write(*,*) line(iIter+5:iRes-1)
c        write(*,*) line(iRes+7:200)
        read(line(iIter+5:iRes-1),*) iter
        read(line(iRes+7:200),*) residual
        write(20,*) iter, residual
      goto 10
20    continue
      close(10)
      close(20)
999   continue
      end

