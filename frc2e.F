


      subroutine frc2e (iout, iprint, v, lenv, iphf, sval, na, nb, 
     $     isint, natoms, nbasis, ngrda, ngrdb, ngrdg, indarr, grda,
     $     grdb, grdg, wgta, wgtb, wgtg, ifmm, iseall, jsym2e, accdes,
     $     irwneq, irwneqshl, irwneqsh2, irwmol, irwh, irws, irwciv,
     $     irwo, irwt, irwd, irwf)

      implicit none


C +------------------------------------------------------+
C |  frc2e   --   RSS,  07.2013                          |
C |                                                      |
C |                                                      |
C |  Compute the two-particle part of the gradient       |
C |                                                      |
C +------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       itype  - decides wich spn blocks need to be kept
C       ncik   - number of basis states for complex conj projection
C       natoms - number of atoms 
C       nbasis - size of the basis (pure spherical harmonics)
C       ngrdX  - number of grid points for variable X
C       ndarr  - index array for grid (see setgrd)
C       ifmm   - flag to set up ab-initio calculation. See phfgrad
C       iseall - see phfgrad
C       jsym2e - controls the use of symmetry
C       irwneq    - RWF file with shell equvalence data
C       irwneqshl - RWF file with shell equvalence data
C       irwneqsh2 - RWF file with shell equvalence data
C       irwmol - Gaussian RWF containing /mol/ common block
C       irwh   - RWF file with the one particle Hamiltonian
C       irwt   - RWF file with the transition density matrices P(g)
C       irwd   - RWF file with 2e forces on grid [ out ]
C       irwf   - RWF file containing the gradient [ out ]
C       accdes - accuracy desired in PRISM

C     input variables

      integer     iout, iprint, lenv, iphf, natoms, nbasis
      integer     sval, na, nb 
      integer     ngrda, ngrdb, ngrdg
      integer     ifmm, iseall, jsym2e
      integer     irwmol, irwh, irws, irwt, irwo, irwf, irwciv
      integer     irwneq, irwneqshl, irwneqsh2
      integer     indarr(*)
      logical     isint

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      real*8      v(*), accdes

C     other varables 

      integer     jend, mdv
      integer     jciv, jfxyz, jfxyzi, jfxyzo, jscr, jdmt, jdind
      integer     ncis, ncik, nci, ngrdt, itype, nbas6d
      integer     nbsq, ntt, nat3, nmatsv
      integer     irwd
      logical     ispuhf, ispghf, icmplx, ispin

C     functions
      
      integer     igetrw
      complex*16  trcabc

C     %%%%%%%%%%%%%%%%%
C     %   Initialize  %
C     %%%%%%%%%%%%%%%%%

C     Decipher iphf.

C       icmplx - whether complex conjugation projection is on
C       ispuhf - whether UHF-based spin projection is used
C       ispghf - whether GHF-based spin projection is used
C       ispin  - whether spin projection is on
C       itype  - decides which spin blocks of transition density
C                matrices are stored

C       itype  = 1,  closed shell [alpha-alpha block]
C              = 2,  open shell [alpha-alpha, beta-beta blocks]
C              = 3,  general [all spin blocks are active]

      icmplx = iphf .le. 3 .or. iphf .eq. 5 .or. iphf .eq. 7
      ispuhf = iphf .eq. 4 .or. iphf .eq. 5
      ispghf = iphf .eq. 6 .or. iphf .eq. 7

      ispin  = iphf .ge. 4

      itype = 1

      if ( iphf .eq. 2 ) itype = 2
      if ( iphf .ge. 3 ) itype = 3

C     Error checking.

      if ( iphf .lt. 1 .or. iphf .gt. 7 ) then
        call gauerr ('Unrecognized iphf in frc1e.')
      endif

C     %%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Prepare calculation  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%

C     Construct useful quantities
C     ===========================

      ntt  = nbasis * (nbasis + 1) / 2
      nbsq = nbasis * nbasis
      nat3 = natoms * 3

C     Determine size of CI problem associated with projection
C       .. number of CI vectors for spin projection (>1 only if GHF)

      if ( ispghf ) then
        if ( isint ) then
          ncis = 2*sval + 1
        else
          ncis = sval + 1
        endif
      else
        ncis = 1
      endif

      nci = ncis

C       .. double the number of CI states if complex conjugation
C          restoration was performed

      ncik = 1
      if ( icmplx ) ncik = 2

      nci  = nci*ncik
      
C     Total number of grid points

      ngrdt = ngrda * ngrdb * ngrdg

C     Total number of matrices to store if
C        complex conjugation projection was performed

      nmatsv = ngrdt * ncik * ncik

C     %%%%%%%%%%%%%%%%%%%%%%%%
C     %  Memory allocatioin  %
C     %%%%%%%%%%%%%%%%%%%%%%%%

C      Memory allocation
C     ===================
C        - jciv, vector of the CI coefficients
C        - jfxyz, integrated forces d(g)
C        - jfxyzo, previous forces value
C        - jscr, scratch for forces, possibly complex

      jciv   = 1
      jfxyz  = jciv   + 2*nci
      jfxyzo = jfxyz  + nat3
      jfxyzi = jfxyzo + nat3
      jscr   = jfxyzi + nat3
      jend   = jscr + nat3*2

C      - jdmt - small Wigner d matrices (if spin projection)
C      - jdind - indicies for small Wigner d matrix (if spin projection)

      if ( ispin ) then
         jdmt  = jend
         jdind = jdmt  + ngrdb * ncis * ncis
         jend  = jdind + ncis
      endif

      mdv   = lenv - jend + 1

      call tstcor (jend - 1, lenv, 'frc1e')

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Compute 2-electron part of gradient %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     Form trace(dG(g)/.dx . P(g)) = d(g)
C     ===================================

      call evaldg (iout, iprint, v(jend), mdv, itype, ncik, 
     $     ngrdt, ifmm, iseall, jsym2e, accdes, irwneq, irwneqshl, 
     $     irwneqsh2, irwmol, irwt, irwd)

C     Load the vector of CI coefficients
      
      call fileio (2, -irwciv, 2*nci, v(jciv), 0)

C     Prepare Wigner's rotation matrices (being integration 
C     weights in spin projection)

      if ( ispin ) then
         call wigarr (iout, iprint, sval, isint, na, nb, ncis, ngrdb,
     $        grdb, v(jdmt), v(jdind))
      endif

C     Integrate forces d(g) over the grid
C     =====================================

      call intdrv (iout, iprint, v(jend), mdv, iphf, ncis, ncik,  
     $     ngrda, ngrdb, ngrdg, ngrdt, indarr, grda, wgta, grdb, wgtb,  
     $     grdg, wgtg, v(jdmt), v(jdind), v(jciv), irwo, irwd,
     $     2*nat3, v(jscr))

C     Extract real part of force

      call vc2ri (nat3, v(jscr), v(jfxyz), v(jfxyzi))

C      call dmpfrc(iout, '2e contribution,real:', 1, natoms,
C     $        v(jfxyz))

C      call dmpfrc(iout, '2e contribution, imag:', 1, natoms,
C     $        v(jfxyzi))

C     Recover previous value and update forces

      call fileio (2, -irwf, nat3, v(jfxyzo), 0)

C      call dmpfrc(iout, 'input to frc2e', 1, natoms, v(jfxyzo))
      
      call aadd (nat3, v(jfxyz), v(jfxyzo), v(jfxyz))

      write (iout, *) '' 
      call dmpfrc(iout, 'Final forces:', 1, natoms,
     $        v(jfxyz))

C     Save      

      call fileio (1, -irwf, nat3, v(jfxyz), 0)

      return 
      end
