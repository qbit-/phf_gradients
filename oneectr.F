


      subroutine oneectr (iout, iprint, v, lenv, itype,
     $     ipflag, allowp, jsym2e, nopuse, nop1, neqbas, neqshl, neqsh2,
     $     rotall, neqall, accdes, pint, lint, irwmol, irwf)

      implicit none


C +------------------------------------------------------+
C |  oneectr   --   RSS,  07.2013                        |
C |                                                      |
C |                                                      |
C |  Contract integrated L-term and integrated           |
C |  transition density matrix with one particle         |
C |  hamiltonian derivatives (interface to PRISM).       |
C |                                                      |
C |  Note: the one-particle hamiltonian and overlap are  |
C |  real symmetric matrces. The only nonzero            |
C |  contributions to forces come from real symmetric    |
C |  and imaginary symmetric parts of P and L.           |
C |  As we don't form P using full projection, it is a   |
C |  general matrix (non-hermitian), same applies to L.  |
C |  Here we simply drop imaginary contributions.        |
C +------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       itype  - decides wich spin blocks need to be kept
C       pint   - integrated transition density matrix
C       lint   - integrated L term
C       irwmol - RWF file with /mol/ common block
C       irwf   - RWF file containing forces [ out, incremented ]

C     include Gaussian common blocks

#include "commonmol.inc"

C     input variables

      real*8      v(*), pint(*), lint(*)
      integer     iout, iprint, lenv, itype
      integer     irwmol, irwf

C       .. prism variables .. [ input variables ]

      real*8      rotall(*), accdes
      integer     ipflag, nopuse, nop1, jsym2e
      integer     neqbas(nbasis,*), neqshl(*), neqsh2(*), neqall(*)
      logical     allowp(50)

C       .. pbc variables ..
      
      integer     npdir
      integer     jpbcsta, jpbcend
      real*8      celvec(3,3)

C     other variables

      integer    nbas6d, ntt, nbsq, ntt6d, nbsq6d, nat3, nat3pe
      integer    iprtf
      integer    jpaa, jpbb, jlaa, jlbb
      integer    jaars, jbbrs, jaais, jbbis
      integer    jfxyz, jfxyzo, jfxyzi, jxrs, jxra, jxis, jxia, jscr
      integer    jend, mdv
      integer    ijunk
      real*8     junk

C     functions

      integer     intowp, itqry

      jend = 1

C     retrieve PBC information

      call getpbc(npdir, celvec)
      jpbcsta = jend
      jpbcend = jpbcsta
      call rdipbc(1, npdir, 0, jpbcsta, jpbcend, v, lenv)
      jend = jpbcend

C     restore /mol/ common block

      call rwmol(2, irwmol)

C     retrieve length of cartesian basis

      call getnb6 (nbas6d)

C     printing level for PRISM
      
      iprtf = iprint - 2
      if ( iprtf .lt. 0 ) iprtf = 0

C     Construct useful quantities. 

      ntt   = nbasis * (nbasis + 1) / 2  
      nbsq  = nbasis * nbasis 

      ntt6d   = nbas6d * (nbas6d + 1) / 2  
      nbsq6d  = nbas6d * nbas6d 

      nat3   = natoms * 3
      nat3pe = nat3 + 3

C     %%%%%%%%%%%%%%%%%%%%%%%%
C     %  Memory allocatioin  %
C     %%%%%%%%%%%%%%%%%%%%%%%%

C      Memory allocation
C     ===================

C      Allocate memory for:
C       - jxaa  - alpha part of matrix
C       - jxbb  - beta  part of matrix
C       - jaars - real symmetric alpha part of matrix
C       - jbbrs - real symmetric beta  part of matrix
C       - jfxyz - forces on atoms
C       - jfxyzo - previous value of forces
C       - jscr   - scratch space
      
      jpaa  = jend
      jpbb  = jpaa + 2*nbsq
      jlaa  = jpbb + 2*nbsq
      jlbb  = jlaa + 2*nbsq
      jend  = jlbb + 2*nbsq

      jaars = jend
      jbbrs = jaars + 2*ntt6d
      jaais = jbbrs + 2*ntt6d
      jbbis = jaais + 2*ntt6d
      jxrs  = jbbis + 2*ntt6d
      jxra  = jxrs  + ntt6d
      jxis  = jxra  + ntt6d
      jxia  = jxis  + ntt6d
      jend  = jxia + ntt6d

      jfxyz  = jend
      jfxyzo = jfxyz  + nat3pe
      jfxyzi = jfxyzo + nat3
      jscr   = jfxyzi + nat3
      jend   = jscr   + 2*nbsq

      call tstcor (jend-1, lenv, 'oneectr')

      mdv = lenv - jend + 1

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %   Prepare calculation    %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     Load previous part of gradient

C      call fileio (2, -irwf, nat3, v(jfxyz), 0)
      call aclear (nat3pe, v(jfxyz))
      call aclear (nat3, v(jfxyzi))

      if ( iprint .gt. 3 ) then
         call dmpfrc(iout, 'input to oneectr', 0, natoms, v(jfxyz))
      endif
      
      call amove (nat3, v(jfxyz+3), v(jfxyzo))

C     Extract spin block from P and L matrices

      if ( itype .eq. 1 ) then
         call amove (2*nbsq, pint, v(jpaa))
         call amove (2*nbsq, pint, v(jpbb))
         call amove (2*nbsq, lint, v(jlaa))
         call amove (2*nbsq, lint, v(jlbb))

      elseif ( itype .eq. 2 ) then
         call amove (2*nbsq, pint, v(jpaa))
         call amove (2*nbsq, pint(2*nbsq + 1), v(jpbb))
         call amove (2*nbsq, lint, v(jlaa))
         call amove (2*nbsq, lint(2*nbsq + 1), v(jlbb))

      elseif ( itype .eq. 3 ) then
         call dmblck (iout, iprint, v(jend), mdv, nbasis, pint,
     $        v(jpaa), v(jscr), v(jscr), v(jpbb), 2, 1)
         call dmblck (iout, iprint, v(jend), mdv, nbasis, lint,
     $        v(jlaa), v(jscr), v(jscr), v(jlbb), 2, 1)         
      endif

 
C     %%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Perform contractions  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%

C     Contract L with dS/.dx and P with dT/.dx..
C     ..alpha part..
      call sptblk (iout, iprint, v(jend), mdv, nbasis, v(jlaa), v(jxrs),
     $     v(jxra), v(jxis), v(jxia), 1)
      call amove (ntt, v(jxrs), v(jaars))
      call amove (ntt, v(jxis), v(jaais))

      call sptblk (iout, iprint, v(jend), mdv, nbasis, v(jpaa), v(jxrs),
     $     v(jxra), v(jxis), v(jxia), 1)
      call amove (ntt, v(jxrs), v(jaars+ntt))
      call amove (ntt, v(jxis), v(jaais+ntt))

C     ..beta part.. 
      call sptblk (iout, iprint, v(jend), mdv, nbasis, v(jlbb), v(jxrs),
     $     v(jxra), v(jxis), v(jxia), 1)
      call amove (ntt, v(jxrs), v(jbbrs))
      call amove (ntt, v(jxis), v(jbbis))

      call sptblk (iout, iprint, v(jend), mdv, nbasis, v(jpbb), v(jxrs),
     $     v(jxra), v(jxis), v(jxia), 1) 
      call amove (ntt, v(jxrs), v(jbbrs+ntt))
      call amove (ntt, v(jxis), v(jbbis+ntt))

C     ..add alpha and beta parts as 1e hamiltonian is spin 
C     independent and oneeli seem not to look at beta part anyway

      call aadd (2*ntt, v(jaars), v(jbbrs), v(jaars))
      call aadd (2*ntt, v(jaais), v(jbbis), v(jaais))
      
      call oneeli(iout, iprtf, 1, 1, 101, 1, 1,
     $  0, 0, 0, natoms, 0, v(jaars), v(jbbrs), nbasis, 2, 0,
     $  ian, c, .false., iattyp, atmchg, natoms, ijunk, junk, 0, 
     $  v(jfxyz), junk, accdes, 0, 0, ijunk, junk, ipflag, allowp, 
     $  .false., v(jpbcsta),celvec,
     $  v(jend),v(jend),mdv)

      if ( iprint .gt. 3 ) then
         call dmpfrc(iout, 'trace(dS/.dx*L + dT/.dx*P)', 0, natoms,
     $        v(jfxyz))
      endif
      call aadd (nat3, v(jfxyz+3), v(jfxyzo), v(jfxyzo))

C      call oneeli(iout, iprtf, 1, 1, 101, 1, 1,
C     $  0, 0, 0, natoms, 0, v(jaais), v(jbbis), nbasis, 2, 0,
C     $  ian, c, .false., iattyp, atmchg, natoms, ijunk, junk, 0, 
C     $  v(jfxyz), junk, accdes, 0, 0, ijunk, junk, ipflag, allowp, 
C     $  .false., v(jpbcsta),celvec,
C     $  v(jend),v(jend),mdv)

C      if ( iprint .gt. 3 ) then
C         call dmpfrc(iout, 'trace(dS/.dx*L + dT/.dx*P) imag', 0, natoms,
C     $        v(jfxyz))
C      endif
C      call aadd (nat3, v(jfxyz+3), v(jfxyzi), v(jfxyzi))

C     Contract P with V..

      call sptblk (iout, iprint, v(jend), mdv, nbasis, v(jpaa), v(jxrs),
     $     v(jxra), v(jxis), v(jxia), 1)
      call amove (ntt, v(jxrs), v(jaars))
      call amove (ntt, v(jxis), v(jaais))

      call sptblk (iout, iprint, v(jend), mdv, nbasis, v(jpbb), v(jxrs),
     $     v(jxra), v(jxis), v(jxia), 1)
      call amove (ntt, v(jxrs), v(jbbrs))
      call amove (ntt, v(jxis), v(jbbis))

C     ..add alpha and beta parts as 1e hamiltonian is spin 
C     independent and oneeli seem not to look at beta part anyway

      call aadd (2*ntt, v(jaars), v(jbbrs), v(jaars))
      call aadd (2*ntt, v(jaais), v(jbbis), v(jaais))
      
      call oneeli(iout, iprtf, 1, 3, 101, 1, 1,
     $  0, 0, 0, natoms, 0, v(jaars), v(jbbrs), nbasis, 1, 0,
     $  ian, c, .false., iattyp, atmchg, natoms, ijunk, junk, 0, 
     $  v(jfxyz), junk, accdes, 0, 0, ijunk, junk, ipflag, allowp, 
     $  .false., v(jpbcsta),celvec,
     $  v(jend),v(jend),mdv)

      if ( iprint .gt. 3 ) then
         call dmpfrc(iout, 'trace(dV/.dx * P)', 0, natoms, v(jfxyz))
      endif
      call aadd (nat3, v(jfxyz+3), v(jfxyzo), v(jfxyzo))

      if ( iprint .gt. 2 ) then
         call dmpfrc(iout, '1e contribution w/o Vnn', 1, natoms,
     $        v(jfxyzo))
      endif

C      call oneeli(iout, iprtf, 1, 3, 101, 1, 1,
C     $  0, 0, 0, natoms, 0, v(jaais), v(jbbis), nbasis, 1, 0,
C     $  ian, c, .false., iattyp, atmchg, natoms, ijunk, junk, 0, 
C     $  v(jfxyz), junk, accdes, 0, 0, ijunk, junk, ipflag, allowp, 
C     $  .false., v(jpbcsta),celvec,
C     $  v(jend),v(jend),mdv)

C      if ( iprint .gt. 3 ) then
C         call dmpfrc(iout, 'trace(dV/.dx * P) imag',0, natoms, v(jfxyz))
C      endif
C      call aadd (nat3, v(jfxyz+3), v(jfxyzi), v(jfxyzi))

      if ( iprint .gt. 3 ) then
         call dmpfrc(iout, '1e contribution w/o Vnn', 1, natoms, 
     $        v(jfxyzo))
         call dmpfrc(iout, '1e contribution w/o Vnn, imag', 1, natoms, 
     $        v(jfxyzi))
      endif
         
C     add nuclear repulsion energy 

      call repnn (iout, iprint, .true., 0, 0, .true., natoms,
     $  iattyp, atmchg, c, junk, 0, v(jfxyz), v(jend), lenv - jend + 1)

      if ( iprint .gt. 3 ) then
         call dmpfrc(iout, 'Vnr', 0, natoms, v(jfxyz))
      endif
      call aadd (nat3, v(jfxyz+3), v(jfxyzo), v(jfxyzo))

C     Save forces

      call dmpfrc(iout, 'PHF grad: One-particle final forces:', 1,
     $     natoms, v(jfxyzo))

      call conddf (nat3, irwf)
      call fileio (1, -irwf, nat3, v(jfxyzo), 0)

      return
      end
