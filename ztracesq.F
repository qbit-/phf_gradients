

      subroutine ztracesq (a, n, re, im, add)

      implicit none

C +------------------------------------------------------------+
C |  ztracesq   --   RSS,  06.2013                             |
C |                                                            |
C |                                                            |
C |  Computes traces of square complex matrices.               |
C |                                                            |
C +------------------------------------------------------------+

C     input / output variables

C       a      - imput matrix
C       n      - size of the matrix
C       re     - real part of trace [out]
C       im     - imaginary part of trace [out]
C      add     - whether to accumulate the values in Re and Im

      integer  n, i
      real*8   a(2,n,n), re, im
      logical  add
      
      if ( .not. add ) then
         re = 0.0D0
         im = 0.0D0
      endif
      
      do 10 i = 1, n
         re = re + a(1,i,i)
         im = im + a(2,i,i)
 10   continue
      
      return
      end
