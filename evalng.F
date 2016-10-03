

      subroutine evalng (iout, iprint, v, lenv, itype, ispin,     
     $     nbasis, na, nb, ngrdt, igrdt, indarr, grda, grdb, grdg,  
     $     mo, mot, s, n1, n2, ng, detng)

      implicit none


C +------------------------------------------------------------+
C |  evalng   --   RSS,  05.2013                               |
C |                                                            |
C |                                                            |
C |  Build the overlap matrices, N(g), necessary to            |
C |  evaluate the transition density matrices, and their       |
C |  determinants, needed for integration over the grid.       |
C |                                                            |
C |  Given a rotation matrix R(g), related to the spin         |
C |  projection operator with an angle g, and given a          |
C |  MO coefficients rectangular matrices D^*, the matrix      |
C |  N(g) can take the form:                                   |
C |                                                            |
C |  if nu = I, mu = I  ( <0| I! R(g) I |0> )                  |
C |                                                            |
C |    [N(g)]^{-1}  =  C! .  R(g) . C                          |
C |                                                            |
C |  if nu = K, mu = I  ( <0| I! R(Omega) K |0> )              |
C |                                                            |
C |    [N(g)]^{-1}  =  C! . R(g) . C*                          |
C |                                                            |
C |  if nu = I, mu = K  ( <0| K! R(Omega) I |0> )              |
C |                                                            |
C |    [N(g)]^{-1}  =  C^T  . R(g) . C                         |
C |                                                            |
C |  if nu = K, mu = K  ( <0| K! R(Omega) K |0> )              |
C |                                                            |
C |    [N(g)]^{-1}  =  C^T  . R(g) . C*                        |
C |                                                            |
C |  where nu, mu = I, K are operators related to complex      |
C |  conjugation restoration.                                  |
C |                                                            |
C |  Here, T denotes the transpose of a matrix, * denotes      |
C |  complex conjugation. The MO coefficients matrices  are    |
C |  rectangular matrices of size nbasis x norbs. The R(g)     |
C |  matrices are always nbasis x nbasis.                      |
C |  This routine computes one or a tuple of four matrices     |
C |  N(g) (if complex conjugation projection is used)  and     |
C |  places them to the array ng                               |
C |                                                            |
C |  Memory requirement:                                       |      
C |              szmo = 2*nel*nbasis or less                   |
C |            + szn  = 2*nel*nel or less /buffer for N/       |
C |            + szn                          /LAPACK/         |
C |            + 2*(na + nb)                  /pivot/          |
C |            + 8*nbsq + 32*nbsq if spin projection is on     |  
C +------------------------------------------------------------+

C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
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
C       n1     - left complex conjugation operator (1 or 2)
C       n2     - right  complex conjugation operator (1 or 2)
C       ng     - overlap matrix N(g) [ out ]
C       detng  - the determinant of matrix N(g) (complex) [ out ]

      real*8      grda(*), grdb(*), grdg(*)
      integer     iout, iprint, lenv, itype, ncik, nbasis, na, nb, ngrdt
      integer     igrdt
      logical     ispin

      real*8      v(*),  mo(*), mot(*), s(*), ng(*)
      integer     indarr(ngrdt,*)
      complex*16  detng

C     other variables

      real*8      adet, pdet, adetb, pdetb
      integer     nel, nbsq, szn, sznc, n1, n2
      integer     rgmode, info
      integer     jxmo, jxmoa, jxmob, jxrg, jnbuf, jpiv, jmltr, jlscr
      integer     jend

      integer     igrda, igrdb, igrdg
      real*8      anga, angb, angg

C     functions called

      integer     intowp

C     debug 
      
      integer     jj
C     constants

      complex*16  zero, one, im

      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)

C     Define useful quantities:
C       - size of matrices N(g) in complex numbers (sznc)
C       - size of matrices N(g) in real numbers (szn)

      nel  = na + nb
      nbsq = nbasis*nbasis

      if ( itype .eq. 1 ) then
         sznc  = na*na
         szn   = 2*sznc
 
      elseif ( itype .eq. 2 ) then 
         sznc  = (na*na + nb*nb)
         szn   = 2*sznc
 
      elseif ( itype .eq. 3 ) then
         sznc  = nel*nel
         szn   = 2*sznc
      endif

C     Allocate memory for:
C       - scratch space for multiplication of S and MO (jxmo, jxmoa, jxmob)
C       - scratch space for mltdg function (jmltr) (only if ispin = .true.)
C       - scratch space to  keep S*R(g) product (jxrg) (only if ispin = .true.)
C       - sctarch space for LAPACK operations (jlscr)
C       - buffer for N(g) (jnbuf)
C       - pivoting vector to compute matrix inverse (jpiv)
      
      if ( itype .eq. 1 ) then
         jxmo  = 1
         jlscr = jxmo  + 2*na*nbasis
         jnbuf = jlscr + szn
         jpiv  = jnbuf + szn
         jend  = jpiv  + intowp(na)

      elseif ( itype .eq. 2 ) then 
         jxmoa = 1
         jxmob = jxmoa + 2*na*nbasis 
         jlscr = jxmob  + 2*nb*nbasis
         jnbuf = jlscr + szn
         jpiv  = jnbuf + szn
         jend  = jpiv  + intowp(nel)

      elseif ( itype .eq. 3 ) then
         jxmo  = 1
         jxrg  = jxmo  + 2*nel*2*nbasis
         jmltr = jxrg  + 8*nbsq
         jlscr = jmltr + 32*nbsq
         jnbuf = jlscr + szn
         jpiv  = jnbuf + szn
         jend  = jpiv  + intowp(nel)

      endif
      
      call tstcor (jend-1, lenv, 'formng')

C     Get grid point coordinates

      igrda = indarr(igrdt,1)
      igrdb = indarr(igrdt,2)
      igrdg = indarr(igrdt,3)
      
      anga = grda(igrda)
      angb = grdb(igrdb)
      angg = grdg(igrdg)

C     form a full spin-orbital overlap/rotation matrix if
C     all spin blocks are present

      if ( itype .eq. 3 ) then
         if ( ispin ) then
C           ..Will use S*R(g) rotation matrix
            rgmode = 10
         else
C     .     .Will use only overlap S*I
            rgmode = 11
         endif
               
         call mltdg (iout, iprint, v(jmltr), 32*nbsq, nbasis, 
     $        v(jxrg), s, anga, angb, angg, rgmode)

      endif

C     Compute matrices N(g)

C     Form <0| I! R(g) = T1
      if ( n1 .eq. 1 ) then
         if ( itype .eq. 1 ) then
C           ..C! . S = T1 aa block 
            call zgemm ('n', 'n', na, nbasis, nbasis, one,
     $           mot, na, s, nbasis, zero, v(jxmo),
     $           na)
            
         elseif ( itype .eq. 2 ) then
C           ..C! . S = T1  aa block
            call zgemm ('n', 'n', na, nbasis, nbasis, one,
     $           mot, na, s, nbasis, zero, v(jxmoa),
     $           na)

C           ..C! . S = T1  bb block
            call zgemm ('n', 'n', nb, nbasis, nbasis, one,
     $           mot(2*na*nbasis+1), nb, s, nbasis, zero, 
     $           v(jxmob), nb)
            
         elseif ( itype .eq. 3 ) then
C           ..C! . S = T1 all blocks
            call zgemm ('n', 'n', nel, 2*nbasis, 2*nbasis, one,
     $           mot, nel, v(jxrg), 2*nbasis, zero, v(jxmo),
     $           nel)
         endif
               
C     Form <0| K! R(g) = T1
      else
         if ( itype .eq. 1 ) then
C           ..C^T . S = T1 aa block 
            call zgemm ('t', 'n', na, nbasis, nbasis, one,
     $           mo, nbasis, s, nbasis, zero, v(jxmo),
     $           na)
                  
         elseif ( itype .eq. 2 ) then
C           ..C^T . S = T1  aa block
            call zgemm ('t', 'n', na, nbasis, nbasis, one,
     $           mo, nbasis, s, nbasis, zero, v(jxmoa),
     $           na)
C           ..C^T . S = T1  bb block
            call zgemm ('t', 'n', nb, nbasis, nbasis, one,
     $           mo(2*na*nbasis+1), nbasis, s, nbasis, zero, 
     $           v(jxmob), nb)
            
         elseif ( itype .eq. 3 ) then
C           ..C^T . S = T1 all blocks
            call zgemm ('t', 'n', nel, 2*nbasis, 2*nbasis, one,
     $           mo, 2*nbasis, v(jxrg), 2*nbasis, zero, v(jxmo),
     $           nel)
         endif
      endif
         
C     Form T1 . C = <0| I! R(g) I |0>
      if ( n2 .eq. 1 ) then
         if ( itype .eq. 1 ) then
C           ..T1 . C = T2 aa block
            call zgemm ('n', 'n', na, na, nbasis, one,
     $           v(jxmo), na, mo, nbasis, zero, v(jnbuf),
     $           na)
            
         elseif ( itype .eq. 2 ) then
C           ..T1 . C = T2 aa block
            write (iout, *) 'Orbitals'
            write (iout, '(F22.18)') (mo(jj), jj = 1,
     $           2*nel*nbasis)

            call zgemm ('n', 'n', na, na, nbasis, one,
     $           v(jxmoa), na, mo, nbasis, zero, v(jnbuf),
     $           na)
C           ..T1 . C = T2  bb block
            call zgemm ('n', 'n', nb, nb, nbasis, one,
     $           v(jxmob), nb, mo(2*na*nbasis+1), nbasis, 
     $           zero, v(jnbuf + 2*na*na),nb)
            
         elseif ( itype .eq. 3 ) then
C           ..T1 . C = T2 all blocks
            call zgemm ('n', 'n', nel, nel, 2*nbasis, one,
     $           v(jxmo), nel, mo, 2*nbasis, zero, v(jnbuf),
     $           nel)
         endif

C           Form T1 . C* = <0| I! R(g) K |0>
      else
         if ( itype .eq. 1 ) then
C           ..T1 . C* = T2 aa block
            call zgemm ('n', 't', na, na, nbasis, one,
     $           v(jxmo), na, mot, na, zero, v(jnbuf),
     $           na)
            
         elseif ( itype .eq. 2 ) then
C           ..T1 . C* = T2 aa block
            call zgemm ('n', 't', na, na, nbasis, one,
     $           v(jxmoa), na, mot, na, zero, v(jnbuf),
     $           na)
C           ..T1 . C* = T2  bb block
            call zgemm ('n', 't', nb, nb, nbasis, one,
     $           v(jxmob), nb, mot(2*na*nbasis+1), nb, 
     $           zero, v(jnbuf + 2*na*na), nb)
                  
         elseif ( itype .eq. 3 ) then
C           ..T1 . C* = T2 all blocks
            call zgemm ('n', 't', nel, nel, 2*nbasis, one,
     $           v(jxmo), nel, mot, nel, zero, v(jnbuf),
     $           nel)
         endif
      endif

C     Evaluate the matrix N(g):
C     N(g) = ( T2 )^{-1}.
      if ( itype .eq. 1 ) then
         call zgetrf (na, na, v(jnbuf), na, v(jpiv), info)
         if ( info .ne. 0 ) call gauerr ('info .ne. 0 in formng')
         
         call detmat (na, v(jnbuf), v(jpiv), adet, pdet,
     $        .false.)
         
         adet = adet*adet
         pdet = 2*pdet

         call zgetri (na, v(jnbuf), na, v(jpiv), v(jlscr), 
     $        na*na, info)
         if ( info .ne. 0 ) call gauerr ('info .ne. 0 in formng')
         
      elseif ( itype .eq. 2 ) then
         call zgetrf (na, na, v(jnbuf), na, v(jpiv), info)
         if ( info .ne. 0 ) call gauerr ('info .ne. 0 in formng')
         
         call detmat (na, v(jnbuf), v(jpiv), adet, pdet)
         
         call zgetri (na, v(jnbuf), na, v(jpiv), v(jlscr), 
     $        na*na, info)
         if ( info .ne. 0 ) call gauerr ('info .ne. 0 in formng')
         
         if ( nb .ge. 1 ) then
            call zgetrf (nb, nb, v(jnbuf + 2*na*na), nb, 
     $           v(jpiv + na), info)
            if ( info .ne. 0 ) call gauerr ('info .ne. 0 in formng')
            
            call detmat (nb, v(jnbuf + 2*na*na), v(jpiv + na), 
     $           adetb, pdetb)
            adet = adet * adetb
            pdet = pdet + pdetb

            call zgetri (nb, v(jnbuf + 2*na*na), nb, 
     $           v(jpiv + na), v(jlscr), nb*nb, info)
            if ( info .ne. 0 ) call gauerr ('info .ne. 0 in formng')
         endif
         
      elseif ( itype .eq. 3 ) then
         call zgetrf (nel, nel, v(jnbuf), nel, v(jpiv), info)
         if ( info .ne. 0 ) call gauerr ('info .ne. 0 in formng')
         
         call detmat (nel, v(jnbuf), v(jpiv), adet, pdet)
         
         call zgetri (nel, v(jnbuf), nel, v(jpiv), v(jlscr), 
     $        nel*nel, info)
         if ( info .ne. 0 ) call gauerr ('info .ne. 0 in formng')
      endif
      
C     Save the results in the output array
      
      call amove (szn, v(jnbuf), ng)
      detng = (cmplx (adet) * exp (im * cmplx (pdet)))

      end
