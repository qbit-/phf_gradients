

      subroutine evalpg (iout, iprint, v, lenv, itype, ispin,
     $     nbasis, na, nb, ngrdt, igrdt, indarr, grda, grdb, grdg,  
     $     mo, mot, ng, n1, n2, pg)

      implicit none

C +------------------------------------------------------------+
C |  evalpg   --   RSS,  06.2013                               |
C |                                                            |
C |                                                            |
C |  Build the transition density matrices, P(g).              |
C |                                                            |
C |  ** NOTE: The order of complex conjugation operators is    |
C |           inverted when forming P(g) comparing to N(g).    |
C |           For example, if we built                         |
C |           N(g) = <0| K! R I |0>^(-1) =                     |
C |           ( C^T . R . C )^(-1), then we use it to build    |
C |           P(g) = R I |0> <0| K! / <0| K! R I |0> =         |
C |           R . C . N . C^T                                  |    
C |                                                            |
C | Memory requirement:                                        |
C |             szmo = 2*nel*2*nbasis or less                  |
C |             szp  = 2*2*2*nbasis or less /buffer for P/     |
C |             32*nbsq if spin projection is on               |
C +------------------------------------------------------------+

C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       ncik   - number of basis states for complex conj projection
C       nbasis - size of the basis
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       ngrdt  - number of the grid points
C       igrdt  - index of the grid point number of grid points
C       indarr - mapping array igrdt => igrd?
C       grd?   - integration grid for each variable
C       mo     - occupied part of the MO matrix in AO basis (C)
C       mot    - transpose of the occupied part of the MO matrix in AO basis (C!)
C       s      - atomic overlap matrix
C       ng     - overlap matrix N(g)
C       n1     - right complex conjugation operator (1 or 2)
C       n2     - left  complex conjugation operator (1 or 2)
C       pg     - transition density matrix P(g) matrices [ out ]

      real*8      grda(*), grdb(*), grdg(*)
      integer     iout, iprint, lenv, itype, ncik, nbasis, na, nb, ngrdt
      integer     igrdt
      logical     ispin

      real*8      v(*), mo(*), mot(*), ng(*), pg(*)
      integer     indarr(ngrdt,*)

C     other variables

      integer     nel, nbsq, szn, szp, n1, n2
      integer     jxmo, jxmoa, jxmob, jpbuf, jmltr, jend

      integer     igrda, igrdb, igrdg
      real*8      anga, angb, angg
      real*8      junk

C     functions called

      integer     intowp

C     constants

      complex*16  zero, one

      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)

C     Define useful quantities:
C       - size of matrices N(g) in complex numbers (sznc)
C       - size of matrices N(g) in real numbers (szn)

      nel  = na + nb
      nbsq = nbasis*nbasis

      if ( itype .eq. 1 ) then
         szn = 2*na*na
         szp = 2*nbsq

      elseif ( itype .eq. 2 ) then 
         szn = 2*(na*na + nb*nb)
         szp = 2*2*nbsq

      elseif ( itype .eq. 3 ) then
         szn = 2*nel*nel
         szp = 2*2*2*nbsq
      endif

C     Allocate memory for:
C       - scratch space for multiplication of MO and N (jxmo)
C       - scratch space for mltdg function (jmltr) (only if ispin = .true.)
C       - scratch space to  keep R(g) (jxrg) (only if ispin = .true.)
C       - buffer for P(g) (jpbuf)
      
      if ( itype .eq. 1 ) then
         jxmo  = 1
         jpbuf = jxmo  + 2*na*nbasis
         jend  = jpbuf + szp

      elseif ( itype .eq. 2 ) then 
         jxmoa = 1
         jxmob = jxmoa + 2*na*nbasis
         jpbuf = jxmob + 2*nb*nbasis
         jend  = jpbuf + szp

      elseif ( itype .eq. 3 ) then
         jxmo  = 1
         jpbuf = jxmo  + 2*nel*2*nbasis
         jmltr = jpbuf + szp
         jend  = jmltr + 32*nbsq

      endif
      
      call tstcor (jend-1, lenv, 'formpg')

C     Get grid point coordinates

      igrda = indarr(igrdt,1)
      igrdb = indarr(igrdt,2)
      igrdg = indarr(igrdt,3)
      
      anga = grda(igrda)
      angb = grdb(igrdb)
      angg = grdg(igrdg)

C     Compute matrices P(g)


C     Form I |0> N = T1
      if ( n2 .eq. 1 ) then
         if ( itype .eq. 1 ) then
C           ..C . N = T1 aa block 
            call zgemm ('n', 'n', nbasis, na, na, one,
     $           mo, nbasis, ng, na, zero, v(jxmo),
     $           nbasis)
            
         elseif ( itype .eq. 2 ) then
C           ..C . N = T1  aa block
            call zgemm ('n', 'n', nbasis, na, na, one,
     $           mo, nbasis, ng, na, zero, v(jxmoa),
     $           nbasis)
C           ..C . N = T1  bb block
            call zgemm ('n', 'n', nbasis, nb, nb, one,
     $           mo(2*na*nbasis+1), nbasis, ng(2*na*na + 1),
     $           nb, zero, v(jxmob), nbasis)
            
         elseif ( itype .eq. 3 ) then
C           ..C . N = T1 all blocks
            call zgemm ('n', 'n', 2*nbasis, nel, nel, one,
     $           mo, 2*nbasis, ng, nel, zero, v(jxmo),
     $           2*nbasis)
         endif
         
C     Form K |0> N = T1
      else
         if ( itype .eq. 1 ) then
C           ..C^T . N = T1 aa block 
            call zgemm ('t', 'n', nbasis, na, na, one,
     $           mot, na, ng, na, zero, v(jxmo),
     $           nbasis)
            
         elseif ( itype .eq. 2 ) then
C           ..C^T . N = T1  aa block
            call zgemm ('t', 'n', nbasis, na, na, one,
     $           mot, na, ng, na, zero, v(jxmoa),
     $           nbasis)
C           ..C^T . N = T1  bb block
            call zgemm ('t', 'n', nbasis, nb, nb, one,
     $           mot(2*na*nbasis+1), nb, ng(2*na*na + 1), nb, zero, 
     $           v(jxmob), nbasis)
            
         elseif ( itype .eq. 3 ) then
C           ..C^T . N = T1 all blocks
            call zgemm ('t', 'n', 2*nbasis, nel, nel, one,
     $           mot, nel, ng, nel, zero, v(jxmo),
     $           2*nbasis)
         endif
      endif
            
C     Form P = T1 . C! = X |0> N <0| I! 
      if ( n1 .eq. 1 ) then
         if ( itype .eq. 1 ) then
C           ..T1 . C! = P aa block
            call zgemm ('n', 'n', nbasis, nbasis, na, one,
     $           v(jxmo), nbasis, mot, na, zero, v(jpbuf),
     $           nbasis)
            
         elseif ( itype .eq. 2 ) then
C           ..T1 . C! = T2 aa block
            call zgemm ('n', 'n', nbasis, nbasis, na, one,
     $           v(jxmoa), nbasis, mot, na, zero, v(jpbuf),
     $           nbasis)
C           ..T1 . C! = T2  bb block
            call zgemm ('n', 'n', nbasis, nbasis, nb, one,
     $           v(jxmob), nbasis, mot(2*na*nbasis+1), nb, 
     $           zero, v(jpbuf + 2*nbsq), nbasis)
            
         elseif ( itype .eq. 3 ) then
C           ..T1 . C! = T2 all blocks
            call zgemm ('n', 'n', 2*nbasis, 2*nbasis, nel, one,
     $           v(jxmo), 2*nbasis, mot, nel, zero, v(jpbuf),
     $           2*nbasis)
         endif
               
C     Form P = T1 . C^T = X |0> R(g) <0| K!
      else
         if ( itype .eq. 1 ) then
C           ..T1 . C^T = T2 aa block
            call zgemm ('n', 't', nbasis, nbasis, na, one,
     $           v(jxmo), nbasis, mo, nbasis, zero, v(jpbuf),
     $           nbasis)
            
         elseif ( itype .eq. 2 ) then
C           ..T1 . C^T = T2 aa block
            call zgemm ('n', 't', nbasis, nbasis, na, one,
     $           v(jxmoa), nbasis, mo, nbasis, zero, v(jpbuf),
     $           nbasis)
C           ..T1 . C^T = T2  bb block
            call zgemm ('n', 't', nbasis, nbasis, nb, one,
     $           v(jxmob), nbasis, mo(2*na*nbasis+1), 
     $           nbasis, zero, v(jpbuf + 2*nbsq), nbasis)
                  
         elseif ( itype .eq. 3 ) then
C           ..T1 . C^T = T2 all blocks
            call zgemm ('n', 't', 2*nbasis, 2*nbasis, nel, one,
     $           v(jxmo), 2*nbasis, mo, 2*nbasis, zero, v(jpbuf),
     $           2*nbasis)
         endif
      endif
      
C     ...If spin projection is done, apply rotation
      if ( ispin ) then
         call mltdg (iout, iprint, v(jmltr), 32*nbsq, nbasis, 
     $        v(jpbuf), junk, anga, angb, angg, 1)
      endif

C     Save the results in the output array

      call amove (szp, v(jpbuf), pg)
            
      end
