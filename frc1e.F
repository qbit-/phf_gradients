


      subroutine frc1e (iout, iprint, v, lenv, iphf, sval, na, nb, 
     $     isint, nbasis, eref, ngrda, ngrdb, ngrdg, indarr, grda, grdb,
     $     grdg, wgta, wgtb, wgtg, allowp, ipflag, jsym2e, nopuse,
     $     nop1, neqbas, neqshl, neqsh2, rotall, neqall, accdes,
     $     irwmol, irwh, irws, irwciv, irwo, irwt, irwg, irwl, irwf)

      implicit none


C +------------------------------------------------------+
C |  frc1e   --   RSS,  06.2013                          |
C |                                                      |
C |                                                      |
C |  Compute the one-particle part of the gradient       |
C |                                                      |
C +------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       itype  - decides wich spn blocks need to be kept
C       ncik   - number of basis states for complex conj projection
C       nbasis - size of the basis (pure spherical harmonics)
C       eref   - reference electronic energy
C       ngrdt  - total number of grid points
C       irwmol - Gaussian RWF containing /mol/ common block
C       irwh   - RWF file with the one particle Hamiltonian
C       irwt   - RWF file with the transition density matrices P(g)
C       irwg   - RWF file with transition Fock matrices  [ out ]
C       irwl   - RWF file with L-terms for 1e derivatives  [ out ]
C       irwf   - RWF file containing the gradient [ out ]
C       accdes - accuracy desired in PRISM

      integer     iout, iprint, lenv, iphf, nbasis
      integer     sval, na, nb 
      integer     ngrda, ngrdb, ngrdg
      integer     irwmol, irwh, irws, irwt, irwo, irwf, irwciv
      integer     indarr(*)
      logical     isint

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      real*8      v(*), eref, accdes

C       .. prism variables .. [ input variables ]

      real*8      rotall(*)
      integer     ipflag, nopuse, nop1, jsym2e
      integer     neqbas(nbasis,*), neqshl(*), neqsh2(*), neqall(*)
      logical     allowp(50)

C     other varables 

      integer     jend, mdv
      integer     jciv, jpint, jlint, jdmt, jdind
      integer     jh, jhf, jscr, js, jsf
      integer     ncis, ncik, nci, ngrdt, itype, nbas6d
      integer     nbsq, ntt, nmatsv, szmat
      integer     irwg, irwl
      logical     ispuhf, ispghf, icmplx, ispin

      complex*16  enr, ove

C     constants
      
      complex*16  two

C     debug 
      integer jtst

C     functions
      
      integer     igetrw
      complex*16  trcabc

      two = cmplx (2.0d0, 0.0d0)

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
      nbsq = nbasis*nbasis

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

C     size of matrices to work with
      
      if ( itype .eq. 1 ) then
         szmat = 2 * nbsq
      elseif ( itype .eq. 2 ) then
         szmat = 2 * 2 * nbsq
      elseif ( itype .eq. 3 ) then
         szmat = 2 * 4 * nbsq
      endif

C     %%%%%%%%%%%%%%%%%%%%%%%%
C     %  Memory allocatioin  %
C     %%%%%%%%%%%%%%%%%%%%%%%%

C      Memory allocation
C     ===================
C        - jciv, vector of the CI coefficients
C        - jpint, integrated transition density matrix P(g)
C        - jlint, integrated matrix L(g)

      jciv  = 1
      jpint = jciv  + 2*nci
      jlint = jpint + szmat
      jend  = jlint + szmat

C      - jciv - vector of CI coefficients
C      - jdmt - small Wigner d matrices (if spin projection)
C      - jdind - indicies for small Wigner d matrix (if spin projection)

      if ( ispin ) then
         jdmt  = jend
         jdind = jdmt  + ngrdb * ncis * ncis
         jend  = jdind + ncis
      endif

C     debug...
C         - alpha part of 1-particle Hamiltonian (jh)
C         - full hamiltonian
C         - scratch space

      jh   = jend 
      jhf  = jh  + 2*ntt
      js   = jhf + szmat
      jsf  = js  + 2*ntt
      jend = jsf + szmat

      if ( itype .eq. 3 ) then 
         jscr = jend 
         jend = jscr + 2*nbsq
      endif

C     end debug...

      mdv   = lenv - jend + 1

      call tstcor (jend - 1, lenv, 'frc1e')

C     debug 

      jtst = jend 
      jend = jtst + 2*ncik*ncik

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Compute 1-electron part of gradient %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     Form G(x) - the two-electron part of the Hamiltonian
C     ====================================================

      call evalgg (iout, iprint, v(jend), mdv, itype, ncik, 
     $     nbasis, ngrdt, irwt, irwg, allowp,
     $     ipflag, jsym2e, nopuse, nop1, neqbas, neqshl, neqsh2, rotall,
     $     neqall, accdes)

C     Form L(x) - the term to contract wth overlap derivatives
C     ========================================================

      call evallg (iout, iprint, v(jend), mdv, itype, ncik, nbasis,
     $     eref, ngrdt, irwh, irwg, irwt, irwl)

C     Load the vector of CI coefficients
      
      call fileio (2, -irwciv, 2*nci, v(jciv), 0)

C     Prepare Wigner's rotation matrices (being integration 
C     weights in spin projection)

      if ( ispin ) then
         call wigarr (iout, iprint, sval, isint, na, nb, ncis, ngrdb,
     $        grdb, v(jdmt), v(jdind))
      endif

      call intdrv (iout, iprint, v(jend), mdv, iphf, ncis, ncik,  
     $     ngrda, ngrdb, ngrdg, ngrdt, indarr, grda, wgta, grdb, wgtb,  
     $     grdg, wgtg, v(jdmt), v(jdind), v(jciv), irwo, 824,
     $     2, v(jtst))

C     Integrate P(g) and L(g) over the grid
C     =====================================

      call intdrv (iout, iprint, v(jend), mdv, iphf, ncis, ncik,  
     $     ngrda, ngrdb, ngrdg, ngrdt, indarr, grda, wgta, grdb, wgtb,  
     $     grdg, wgtg, v(jdmt), v(jdind), v(jciv), irwo, irwt,
     $     szmat, v(jpint))

      call intdrv (iout, iprint, v(jend), mdv, iphf, ncis, ncik,  
     $     ngrda, ngrdb, ngrdg, ngrdt, indarr, grda, wgta, grdb, wgtb,  
     $     grdg, wgtg, v(jdmt), v(jdind), v(jciv), irwo, irwl,
     $     szmat, v(jlint))


C     Load core Hamiltonian

      if ( itype .eq. 1 ) then
         call fileio (2, -irwh, 2*ntt, v(jh), 0)
         call zsquare (v(jh), v(jhf), nbasis, nbasis, 0)
         call fileio (2, -irws, 2*ntt, v(js), 0)
         call zsquare (v(js), v(jsf), nbasis, nbasis, 0)

      elseif ( itype .eq. 2 ) then
         call fileio (2, -irwh, 2*ntt, v(jh), 0)
         call zsquare (v(jh), v(jhf), nbasis, nbasis, 0)
         call fileio (2, -irwh, 2*ntt, v(jh), 0)
         call zsquare (v(jh), v(jhf+2*nbsq), nbasis, nbasis, 0)

         call fileio (2, -irws, 2*ntt, v(js), 0)
         call zsquare (v(js), v(jsf), nbasis, nbasis, 0)
         call fileio (2, -irws, 2*ntt, v(js), 0)
         call zsquare (v(js), v(jsf+2*nbsq), nbasis, nbasis, 0)

      elseif ( itype .eq. 3 ) then
         call fileio (2, -irwh, 2*ntt, v(jh), 0)
         call aclear (2*nbsq, v(jscr))
         call dmblck (iout, iprint, v(jend), mdv, nbasis, v(jhf), 
     $        v(jh), v(jscr), v(jscr), v(jh), 1, -2)

         call fileio (2, -irws, 2*ntt, v(js), 0)
         call aclear (2*nbsq, v(jscr))
         call dmblck (iout, iprint, v(jend), mdv, nbasis, v(jsf), 
     $        v(js), v(jscr), v(jscr), v(js), 1, -2)

      endif

      call vconjg (szmat, v(jhf))
      call vconjg (szmat, v(jsf))
      
      if ( itype .eq. 1 ) then
         enr = trcabc(1, nbasis, v(jpint), v(jhf))
         ove = trcabc(1, nbasis, v(jpint), v(jsf))
         enr = enr * two

      elseif ( itype .eq. 2 ) then
         enr = trcabc(2, nbasis, v(jpint), v(jhf))
         ove = trcabc(2, nbasis, v(jpint), v(jsf))

      elseif ( itype .eq. 3 ) then
         enr = trcabc(1, 2*nbasis, v(jpint), v(jhf))
         ove = trcabc(1, 2*nbasis, v(jpint), v(jsf))

      endif
      
      call vconjg (szmat, v(jhf))
      call vconjg (szmat, v(jsf))

      write (iout, *) ''
      write (iout, *) 'Correctness test:'
      write (iout, *) '   Hsp = ', enr
      write (iout, *) '   Hph = ', cmplx(v(jtst), v(jtst+1))
      write (iout, *) '   dE = ', enr + cmplx(v(jtst), v(jtst+1)) - eref
      write (iout, *) '   Overlap = ', ove
      write (iout, *) ''

C     Contract P and L with overlap and hamiltonian derivatives
C     =========================================================

      call oneectr (iout, iprint, v(jend), mdv, itype,
     $     ipflag, allowp, jsym2e, nopuse, nop1, neqbas, neqshl, neqsh2,
     $     rotall, neqall, accdes, v(jpint), v(jlint), irwmol, irwf)


      return 
      end
