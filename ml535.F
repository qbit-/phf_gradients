

#ifdef GAUSS_LINDA
      subroutine Real_Main
#else
      program ML535
#endif

      implicit none

      real*8      work(1)

      integer     in, iout, ipunch
      common      /io/ in, iout, ipunch

#ifdef _ALIGN_CORE_
      common      /gxwork/ work
#endif

#ifdef P64
      integer*8   ioff
#else
      integer     ioff
#endif

      integer     istat, mdv
      integer     gsetjmp

      in = 5
      iout = 6
      ipunch = 7
      istat = gsetjmp (0)

      if ( istat .eq. 0 ) then
        call initscm (1, 0, 0, work, ioff, mdv)
        call phfgrad (work(ioff+1), mdv)
      else
        call prtstat (istat, 'ML535')
      endif


      end


