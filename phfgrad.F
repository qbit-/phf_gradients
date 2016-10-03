
      subroutine phfgrad (v, ngot)

      implicit none


C +------------------------------------------------------------+
C |  PHF-GRAD   ---   the Scuseria group, 07.2011              |
C |                                                            |
C |                                                   v0.2 -   |
C |                                                            |
C |    ( link 535, part of the Gaussian system of programs )   |
C |                                                            |
C |                                                            |
C |  This link computes projected (spin projection, complex    |
C |  conjugation restoration) Hartree-Fock (PHF) analytic      |
C |  gradients.                                                |
C |                                                            |
C |                                                            |
C |  Please contact Prof. Gustavo Scuseria for any further     |
C |  questions about the method and / or the use of the code.  |
C |  Contact information is available at:                      |
C |                                                            |
C |    http://scuseria.rice.edu/                               |
C |                                                            |
C |                                                            |
C |  IMPORTANT NOTE:                                           |
C |                                                            |
C |    We do not guarantee that this software is bug-free.     |
C |    The use of this code is sole responsibility of the      |
C |    user.                                                   |
C |                                                            |
C +------------------------------------------------------------+

C +------------------------------------------------------------+
C |                                                            |
C |  Description of IOps                                       |
C |  ===================                                       |
C |                                                            |
C |                                                            |
C |  iop(6) - printing level                                   |
C |     = 0, [default] print results at convergence            |
C |     = 1, print energy at every iteration                   |
C |     = 2, print CI matrices at every iteration              |
C |     = 3, print initial and converged density matrices      |
C |     = 4, print density matrix at every iteration           |
C |     = 5, print (gauge angle) matrices at every iteration   |
C |     ...                                                    |
C |     = 10, print two-electron integrals                     |
C |                                                            |
C |  iop(8) - number of grid points in integral over alpha     |
C |           for spin projection (GHF)                        |
C |     = 0  [defaults to 9]                                   |
C |     = N, use N grid points (trapezoidal rule)              |
C |                                                            |
C |  iop(7) - iphf type of PHF calculation                     |
C |     = 1, closed-shell HF + complex conjugation projection  |
C |     = 2, open-shell HF   + complex conjugation projection  |
C |     = 3, general HF      + complex conjugation projection  |
C |     = 4, open-shell HF   + spin projection                 |
C |     = 5, open-shell HF   + spin + complex projection       |
C |     = 6, general HF      + spin projection                 |
C |     = 7, general HF      + spin + complex projection       |
C |                                                            |
C |  iop(9) - number of grid points in integral over beta      |
C |           for spin projection (UHF, GHF)                   |
C |     = 0  [defaults to 8]                                   |
C |     = N, use N grid points (Gauss-Legendre quadrature)     |
C |                                                            |
C |  iop(10) - number of grid points in integral over gamma    |
C |            for spin projection (GHF)                       |
C |     = 0  [defaults to 9]                                   |
C |     = N, use N grid points (trapezoidal rule)              |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |       TEMPORARY IOPS. WILL PASS THIS INFO VIA RWF          | 
C |                                                            |
C |  iop(11) - S value to use for spin projection              |
C |     =  0  [defaults to S = Sz = (Na - Nb)/2]               |
C |     =  N  project into S = N                               |
C |     = -N  project into S = N/2                             |
C |                                                            |
C |     ** Note that S >= Sz = (Na - Nb) / 2.                  |
C |        For half-integer projection, only odd values        |
C |        should be used.                                     |
C |                                                            |
C |  iop(12) - PHF Energy (eref)                               |
C +------------------------------------------------------------+


C +------------------------------------------------------------+
C |  phfgrad  ---   RSS3, 05.2013                              |
C |                                                            |
C |  Driver routine for PHF derivatives calculation.           |
C |                                                            |
C +------------------------------------------------------------+



C     input / output variables

      real*8      v(*)
      integer     ngot


C     include some Gaussian files

#include "commonb.inc"
#include "commonmol.inc"
#include "repall.inc"


C     iop common block

      integer     iop, maxiop
      parameter   ( maxiop = DEFMAXIOP )

      common      /iop/ iop(maxiop)

C     other common blocks

      integer     in, iout, ipunch
      integer     ibf, isym2e

      common      /io/ in, iout, ipunch
      common      /ibf/ ibf(30)

      equivalence (isym2e, ibf(30))

C     general array

      integer     lengen
      parameter   ( lengen = 55 )

      real*8      dgen(lengen)


C     prism variables

      integer     ifmm, ipflag, fmflag, fmflg1, nfxflg, ihmeth, iseall
      integer     lseall, momega, nomega
      logical     allowp(50), fmm
      parameter   ( momega = DEFMOMEGA, nomega = DEFNOMEGA )

      real*8      omega(momega,nomega)


C     general options variables

      real*8      cnvdef, acurcy, accdes
      integer     iprint, iphf, itype
      integer     ngrda, ngrdadf, ngrdamx
      integer     ngrdb, ngrdbdf, ngrdbmx
      integer     ngrdg, ngrdgdf, ngrdgmx
      integer     sval, sval1, ncik
      logical     icmplx, ispuhf, ispghf, isint

      parameter   ( ngrdadf = 9, ngrdamx = 99 )
      parameter   ( ngrdbdf = 8, ngrdbmx = 99 )
      parameter   ( ngrdgdf = 9, ngrdgmx = 99 )

C     memory allocation variables

      integer     jstrt, jend, lenv, mdv
      integer     ngrdt, szmat, nat3
      integer     jgrda, jgrdb, jgrdg, jwgta, jwgtb, jwgtg
      integer     indarr

C     PHF gradients RWF files
C       - irwpa1 - density matrix stored as spin blocks (lower triangular)
C                  in case of PRHF it this is a single block,
C                  in case of PUHF - two blocks, PGHF - whole matrix
C       - irwca1 - spin orbitals rearranged as spin block
C       - irwt   - transition density matrices at each grid point
C       - irwo   - overlap between the initial and transformed wavefunction
C                  at each grid point
C       - irwd   - RWF file with 2e forces on grid

      integer     irwpa1, irwca1
      integer     irwt, irwo, irwg, irwl, irwd

C      parameter   ( irwt = 820, irwo = 821, irwg = 822, irwl = 823, 
C     $     irwd = 824)

C     PHF RWF files
C       - irwciv - CI coefficients of the CI-like problem associated with PHF

      integer     irwciv
      parameter   ( irwciv = 806 )

C     Gaussian RWF files
C       - irwgen - general array (see Link 1 for details)
C       - irwibf - file with IBF array
C       - irweig - MO orbital energies
C       - irwc?  - MO (x-spin) coefficients
C       - irwp?  - density matrix (x-spin)
C       - irws   - overlap matrix
C       - irwx   - transformation matrix (=S^(-1/2))
C       - irwh   - core Hamiltonian matrix
C       - irwf   - cartesian forces on atoms, FXYZ
C       - irwmol - /mol/ common block
      
      integer     irwgen, irwibf, irweig
      integer     irwca, irwcb, irwpa, irwpb, irws, irwx, irwh, irwf
      integer     irwmol

      parameter   ( irwgen = 501, irwibf = 508, irweig = 522 )
      parameter   ( irwca = 524, irwcb = 526, irwpa = 528, irwpb = 530 )
      parameter   ( irws = 514, irwx = 685, irwh = 515, irwf = 584 )
      parameter   ( irwmol = 997 )

C       - neq, neqshl, neqsh2 - files useful for symmetry

      integer     neq, neqshl, neqsh2

      parameter   ( neq = 580, neqshl = 565, neqsh2 = 726 )


C     symmetry related variables

      integer     jsym2e
      integer     nop1, nop2, nopuse
      integer     lenneq, lennes, lenne2
      integer     jneq, jneqsh, jneqs2


C     other variables

      real*8      eref, enuc, etot
      integer     nel, iopcl, junk
      integer     irws1, irwx1, irwh1

C     functions called

      real*8      scfacc, decacc
      integer     intowp, itqry
      integer     igetrw
      integer     nproc

C     format statements

 1001 format (' PHF gradients were invoked with the following options:')
 1010 format ('   method selection:')
 1011 format (' iprint = ', I2)
 1012 format (' iphf   =  1, [ RHF + complex projection ]')
 1013 format (' iphf   =  2, [ UHF + complex projection ]')
 1014 format (' iphf   =  3, [ GHF + complex projection ]')
 1015 format (' iphf   =  4, [ UHF + spin projection ]')
 1016 format (' iphf   =  5, [ UHF + spin + complex projection ]')
 1017 format (' iphf   =  6, [ GHF + spin projection ]')
 1018 format (' iphf   =  7, [ GHF + spin + complex projection ]')
 1020 format ('   integration grid  options:')
 1021 format (' ngrda  = ', I6, ', ngrdb  = ', I6, ', ngrdg  = ', I6)
 1030 format (' Energies received:')
 1031 format (5X, '  E(total) = ', F20.12, 3X, ' E(elec) = ',
     $     F20.12, 3X,)
 1035 format (' requested integral accuracy   = ', 1P, E14.4)
 1040 format ('   use of symmetry:')
 1041 format (' two-electron integral symmetry not used')
 1042 format (' two-electron integrals replicated using symmetry')

      call drum (v, ngot)

      write (iout, *) ' '
      write (iout, *) ' +-------------------------------------------+'
      write (iout, *) ' |  LINK 535 - Projected Hartree-Fock (PHF)  |'
      write (iout, *) ' |     Analytic gradients                    |'
      write (iout, *) ' +-------------------------------------------+'
      write (iout, *) ' '

C     Set up prism and FMM control flags.

C       The calculation will crash if a semi-empirical Hamiltonian was
C       requested. FMM is set to default, but we will also ignore it.

      ifmm = 0
      iseall = 0

      call setpfl (iout, iprint, ifmm, ipflag, allowp, fmm, fmflag,
     $     fmflg1, nfxflg, ihmeth, omega, iseall, lseall, jstrt, v,
     $     ngot)

      if ( jstrt .lt. 1 ) jstrt = 1

      jend = jstrt
      lenv = ngot - jstrt + 1


C     Retrieve nuclear repulsion and total PHF energy and
C     from general array.
      call fileio (2, -irwgen, lengen, dgen, 0)

      enuc = dgen(41)
      etot = dgen(32)
      eref = etot - enuc

      write (iout, 1030)
      write (iout, 1031) etot, eref

C     Read iopcl from ILSW file.

C       iopcl = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)
C             Other values are incorrect after the PHF code run

      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in phfgrad.')
      elseif ( iopcl .lt. 2 ) then
        call gauerr ('Incorrect iopcl in phfgrad.')         
      endif


C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Read options and initialize internal variables  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      write (iout, *)
      write (iout, 1001)
      write (iout, *)


C     Printing level.
C     ===============
 
      iprint = iop(6)

      if ( iprint .lt. 0 .or. iprint .gt. 10 ) iprint = 0


C     Type of PHF calculation
C     =======================

C       iphf - type of PHF calculation

C       iphf = 0,  [defaults to 1, 2, or 3, depending on Gaussian guess]
C            = 1,  closed-shell HF (RHF) + complex conjugation restoration
C            = 2,  open-shell HF   (UHF) + complex conjugation restoration
C            = 3,  general HF      (GHF) + complex conjugation restoration
C            = 4,  open-shell HF   (UHF) + spin projection
C            = 5,  open-shell HF   (UHF) + spin + complex projection
C            = 6,  general HF      (GHF) + spin projection
C            = 7,  general HF      (GHF) + spin + complex projection

      iphf = iop(7)

      if ( iphf .le. 0 .or. iphf .gt. 7 ) then
        call gauerr ('Incorrect iphf in phfgrad.')
      endif


C     Decipher iphf.

C       icmplx - whether complex conjugation projection is on
C       ispuhf - whether UHF-based spin projection is used
C       ispghf - whether GHF-based spin projection is used
C       itype  - decides which spin blocks of transition density
C                matrices are stored

C       itype  = 1,  closed shell [alpha-alpha block]
C              = 2,  open shell [alpha-alpha, beta-beta blocks]
C              = 3,  general [all spin blocks are active]

      icmplx = iphf .le. 3 .or. iphf .eq. 5 .or. iphf .eq. 7
      ispuhf = iphf .eq. 4 .or. iphf .eq. 5
      ispghf = iphf .eq. 6 .or. iphf .eq. 7

      itype = 1

      if ( iphf .eq. 2 ) itype = 2
      if ( iphf .ge. 3 ) itype = 3


      write (iout, 1010)
      write (iout, 1011) iprint

      if ( iphf .eq. 1 ) write (iout, 1012)
      if ( iphf .eq. 2 ) write (iout, 1013)
      if ( iphf .eq. 3 ) write (iout, 1014)
      if ( iphf .eq. 4 ) write (iout, 1015)
      if ( iphf .eq. 5 ) write (iout, 1016)
      if ( iphf .eq. 6 ) write (iout, 1017)
      if ( iphf .eq. 7 ) write (iout, 1018)

      write (iout, *)


C     Integration grid
C     ================

C       ngrda - number of grid points to use for alpha integration
C       ngrdb - number of grid points to use for beta integration
C       ngrdg - number of grid points to use for gamma integration

      ngrda = iop(8)
      ngrdb = iop(9)
      ngrdg = iop(10)

      if ( .not. ( ispghf ) ) then
        ngrda = 1
        ngrdg = 1
      endif

      if ( .not. ( ispuhf .or. ispghf ) ) ngrdb = 1

      if ( ngrda .le. 0 ) ngrda = ngrdadf
      if ( ngrdb .le. 0 ) ngrdb = ngrdbdf
      if ( ngrdg .le. 0 ) ngrdg = ngrdgdf

      ngrda = min (ngrda, ngrdamx)
      ngrdb = min (ngrdb, ngrdbmx)
      ngrdg = min (ngrdg, ngrdgmx)

C     Spin projection
C     ===============

C       sval  - quantum number s on which to project
C       isint - .true. for integer s, .false. for half-integer s

C       sval is determined by iop(11), which selects the value of s
C       desired. Positive numbers should be used for integer s, negative
C       numbers for half-integer s. Note that s >= Sz = (Na - Nb) / 2,
C       where Na, Nb are the number of alpha and beta electrons,
C       respectively.

C       If iop (11) = 0, then sval = Sz = (Na - Nb) / 2.


      if ( abs (iop(11)) .gt. 0 ) sval1 = iop(11)

      if ( iop(11) .eq. 0 .or. .not. ( ispuhf .or. ispghf ) ) then
        if ( mod (nae-nbe, 2) .eq. 0 ) then
          sval = (nae-nbe)/2
          isint = .true.
        else
          sval = nae-nbe
          isint = .false.
        endif

      elseif ( mod (nae-nbe, 2) .eq. 0 ) then
        if ( sval1 .ge. 0 .and. sval1 .ge. (nae-nbe)/2 ) then
          sval = sval1
          isint = .true.
        else
          call gauerr ('Inconsistent multiplicity with iop (11)')
        endif

      else
        if ( sval1 .lt. 0 .and. -sval1 .ge. nae-nbe ) then
          sval = -sval1
          isint = .false.
        else
          call gauerr ('Inconsistent multiplicity with iop (11)')
        endif
      endif

C     Accuracy requested from PRISM
C     =============================

      cnvdef = scfacc (0)

      if ( iop(14) .ge. 0 ) then
        acurcy = scfacc (iop(14))
      else
        acurcy = scfacc (0)
      endif

      accdes = decacc (0, 0.0d0, .false.)

      if ( acurcy .lt. cnvdef ) then
        accdes = accdes * acurcy / cnvdef
      endif

      write (iout, 1035) accdes

C     Use of symmetry
C     ===============

C     We currently set jsym2e = 2 in the code (users cannot change that).

C       jsym2e = 2 uses integral symmetry by replicating integrals, but
C       does not force the density or Fock matrices to be symmetric.
C
C     FIXME: symmetry seems to be prohibited entirely for derivatives. (jsym2e=0)
C            All symmetry related code should be cleaned from this routine
C
      jsym2e = 0

C     Load IBF array.

      call fileio (2, -irwibf, intowp (30), ibf, 0)

      if ( isym2e .eq. 0 ) jsym2e = 0


      write (iout, 1040)

      if ( jsym2e .eq. 0 ) write (iout, 1041)
      if ( jsym2e .eq. 2 ) write (iout, 1042)

      write (iout, *)


C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %    Prepare calculation. Allocate memory and set up    %
C     %    needed quantities.                                 %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     Retrieve symmetry information  
C     =============================

C     Get number of symmetry operations.
C       - nop1,    order of the concise abelian subgroup
C       - nop2,    order of the full abelian subgroup
C       - nopuse,  actual oder to be used

      call getnop (nop1, nop2)
      call pccki (0, junk, neqall, mxatso, nopall, natoms, nopall)

      nopuse = nopall

C     Memory allocation for symmetry related quantities.

      lenneq = 0
      lennes = 0
      lenne2 = 0

      if ( isym2e .eq. 1 ) then
        lenneq = intowp (nop1*nbasis)
        lennes = max (itqry (neqshl), 0)
        lenne2 = max (itqry (neqsh2), 0)
      endif

      jneq   = jend
      jneqsh = jneq   + lenneq
      jneqs2 = jneqsh + lennes
      jend   = jneqs2 + lenne2

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'phfgrad')


C     Fill arrays with symmetry related quantities.

      if ( lenneq .ne. 0 ) then
        call fileio (2, -neq, lenneq, v(jneq), 0)
      endif

      if ( lennes .ne. 0 ) then
        call fileio (2, -neqshl, lennes, v(jneqsh), 0)
      endif

      if ( lenne2 .ne. 0 ) then
        call fileio (2, -neqsh2, lenne2, v(jneqs2), 0)
      endif


C     Define  number of grid points (ngrdt)

      ngrdt = ngrda * ngrdb * ngrdg
      
      nat3  = natoms*3

      ncik = 1
      if ( icmplx ) ncik = 2

      if ( itype .eq. 1 ) then
        szmat   = 2*nbasis*nbasis

      elseif ( itype .eq. 2 ) then
        szmat   = 2*2*nbasis*nbasis

      elseif ( itype .eq. 3 ) then
        szmat   = 2*2*2*nbasis*nbasis

      endif

C     Allocate RWF files

C       ..overlap and transition density matrices..

      irwt =  igetrw (ngrdt*szmat*ncik*ncik, 0)
      irwo =  igetrw (ngrdt*2*ncik*ncik, 0)

C       ..transition Fock matrices and L-matrices..

      irwg = igetrw (ngrdt*szmat*ncik*ncik, 0)
      irwl = igetrw (ngrdt*szmat*ncik*ncik, 0)

C       ..2e forces on grid..
      irwd =  igetrw (ngrdt*nat3*2*ncik*ncik, 0)

C     Allocate memory for the grid 
C     ===============================
        
C     Allocate space for the integration grid (see setgrd)
C       - array with indicies for each spin angle (indarr) 
C       - grid for each angle (jgrd?)
C       - integration weight for each angle (jwgt?)

      indarr = jend
      jgrda  = indarr + intowp (3*ngrdt)
      jwgta  = jgrda  + ngrda
      jgrdb  = jwgta  + ngrda
      jwgtb  = jgrdb  + ngrdb
      jgrdg  = jwgtb  + ngrdb
      jwgtg  = jgrdg  + ngrdg
      jend   = jwgtg  + ngrdg

      call tstcor (jend-1, lenv, 'phfgrad')
      
C     Make sure that the overlap matrix, transformation matrix and core
C     Hamiltonian matrix are in the form we expect them to be.
C     =================================================================

      call fixmat (iout, iprint, v(jend), mdv, iopcl, nbasis, 
     $     nbsuse, irws, irwx, irwh, irws1, irwx1, irwh1)

C     Rearrange MOs and density matrix to spin blocks
C     if GHF-type SCF is used
C     ================================================

      call fixgaufmt (iout, iprint, v(jend), mdv, iopcl, nbasis,
     $     nbsuse, irwpa, irwca, irwpa1, irwca1)


C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Perform a gradient calculation      %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     Fill the arrays defining the integration grid
C     =============================================

      call setgrd (iout, iprint, v(jend), mdv, ngrdt,
     $     ngrda, ngrdb, ngrdg, v(indarr), v(jgrda), 
     $     v(jwgta), v(jgrdb), v(jwgtb), v(jgrdg), 
     $     v(jwgtg))


C     Form the transition density matrices, 
C     overlap matrices N(g) and their traces n(g)
C     ===========================================

      call formden (iout, iprint, v(jend), mdv, iphf, nbasis, 
     $     nbsuse, nae, nbe, ngrdt, v(indarr), v(jgrda), v(jgrdb),
     $     v(jgrdg), irwca1, irwcb, irws1, irwt, irwo)

C     Form the 1e part of the gradient
C     ================================

      call frc1e (iout, iprint, v(jend), mdv, iphf, sval, nae, nbe, 
     $     isint, nbasis,eref,ngrda, ngrdb, ngrdg, v(indarr), v(jgrda),
     $     v(jgrdb), v(jgrdg), v(jwgta), v(jwgtb), v(jwgtg), allowp,
     $     ipflag, jsym2e, nopuse, nop1, v(jneq), v(jneqsh),
     $     v(jneqs2), rotall, neqall, accdes, irwmol, irwh1, irws1, 
     $     irwciv, irwo, irwt, irwg, irwl, irwf)

C     Form the 2e part of the gradient
C     ================================

      call frc2e (iout, iprint, v(jend), mdv, iphf, sval, nae, nbe, 
     $     isint, natoms, nbasis, ngrda, ngrdb, ngrdg, v(indarr),
     $     v(jgrda), v(jgrdb), v(jgrdg), v(jwgta), v(jwgtb), v(jwgtg),
     $     ifmm, iseall, jsym2e, accdes, neq, neqshl, neqsh2,
     $     irwmol, irwh1, irws1, irwciv, irwo, irwt, irwd, irwf)


C     Free scratch RWF files.

      call freerw (irwt)
      call freerw (irwo)
      call freerw (irwg)
      call freerw (irwl)
      call freerw (irwd)

      if ( irws1 .ne. irws ) call freerw (irws1)
      if ( irwx1 .ne. irwx ) call freerw (irwx1)
      if ( irwh1 .ne. irwh ) call freerw (irwh1)
      if ( irwpa1 .ne. irwpa ) call freerw (irwpa1)
      if ( irwca1 .ne. irwca ) call freerw (irwca1)

      call chainX (0)


      return
      end


