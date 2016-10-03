


      subroutine intgrd (iout, iprint, v, lenv, iphf, sval, na, nb,  
     $     isint, ngrda, ngrdb, ngrdg, indarr, grda, grdb, grdg,
     $     wgta, wgtb, wgtg, irwciv, irwo, irwf, nat3)

      implicit none

C +------------------------------------------------------------+
C |  intgrd   --   RSS, 06.2013                                |
C |                                                            |
C |                                                            |
C |                                                            |
C |  Evaluate the forces as a sum of projected values with     |
C |  respective CI coefficients. The CI coefficients are an    |
C |  eigenvector of small CI-like problem, associated with     |
C |  PHF method (see Carlos' paper). They come  along with     |
C |  the solution of the PHF procedure                         |
C |  (see link 530).                                           |
C |                                                            |
C |  The whole theory is too cumbersome to explain in this     |
C |  note and can be found in:                                 | 
C |                                                            |
C |                                                            |
C |  In short, here we form projected values of force          | 
C |  (as a sum for complex conjugation or as an integral,      |
C |  which is discretised, for spin projection).               |
C |  The summation or integration has to be                    |  
C |  done with a specific weight, wich depends on the symmetry |
C |  underlaying the projector. In case of spin conjugation    |
C |  the weight is a Wigner D matrix, in case of complex       |
C |  conjugation the weights are not known analytically and    |
C |  are included in the associated CI problem.                | 
C |  After the projected values were formed we sum each        |
C |  projected value of                                        |
C |  forces (or any other array of scalars) with respective    |
C |  CI coefficients to get the answer.                        |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       sval   - S^2 value to project to
C       isint  - whether s is integer
C       ngrda  - number of grid points in alpha integration
C       ngrdb  - number of grid points in beta integration
C       ngrdg  - number of grid points in gamma integration
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       wgtX   - array with quadrature weights for each variable
C       irwciv - RWF file with vector of CI coefficients associated with proojectors 
C       irwo  - RWF file with overlaps x(g)
C       irwf   - RWF file with forces at each grid point
C       nat3   - number of the elements of gradient

C     input variables 

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      integer     iout, iprint, lenv, iphf, sval, na, nb, nat3
      logical     isint
      integer     ngrda, ngrdb, ngrdg
      integer     indarr(*)
      integer     irwciv, irwf, irwo

      real*8      v(*)

C     other variables
      
      integer    ncis, ncik, nci, ngrdt
      logical    icmplx, ispuhf, ispghf, ispin
      integer    jciv, jdmt, jdind, jfxyz, jend, mdv

C     Decipher iphf.

C       icmplx - whether complex conjugation projection is on
C       ispuhf - whether UHF-based spin projection is used
C       ispghf - whether GHF-based spin projection is used
C       ispin  - whether spin projection is on

      icmplx = iphf .le. 3 .or. iphf .eq. 5 .or. iphf .eq. 7
      ispuhf = iphf .eq. 4 .or. iphf .eq. 5
      ispghf = iphf .eq. 6 .or. iphf .eq. 7

      ispin = iphf .ge. 4
      
C     Error checking.

      if ( iphf .lt. 1 .or. iphf .gt. 7 ) then
        call gauerr ('Unrecognized iphf in intgrd.')
      endif

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

C     Allocate memory for:
C     ====================      
C        - jciv, vector of the CI coefficients
C        - jfxyz, forces
      jciv  = 1
      jfxyz = jciv + 2*nci
      jend  = jfxyz + nat3 * 2

C      - jciv - vector of CI coefficients
C      - jdmt - small Wigner d matrices (for spin projection)
C      - jdind - indicies for small Wigner d

      if ( ispin ) then
         jdmt  = jend
         jdind = jdmt  + ngrdb * ncis * ncis
         jend  = jdind + ncis
      endif

      mdv   = lenv - jend + 1

      call tstcor (jend-1, lenv, 'intgrd')

C     Load the vector of CI coefficients
      
      call fileio (2, -irwciv, 2*nci, v(jciv), 0)

C     Prepare Wigner's rotation matrices (being integration 
C     weights in spin projection)

      if ( ispin ) then
         call wigarr (iout, iprint, sval, isint, na, nb, ncis, ngrdb,
     $        grdb, v(jdmt), v(jdind))
      endif

C     Perform the integration / summation

      call intdrv (iout, iprint, v(jend), mdv, iphf, ncis, ncik,  
     $     ngrda, ngrdb, ngrdg, ngrdt, indarr, grda, wgta, grdb, wgtb,  
     $     grdg, wgtg, v(jdmt), v(jdind), v(jciv), irwo, irwf,
     $     nat3, v(jfxyz))

C     debug
      write (iout, *) 'HSP'
      write (iout, *) v(jfxyz), v(jfxyz+1)
      return
      end












