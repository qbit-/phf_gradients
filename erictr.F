

      subroutine erictr (iout, iprint, v, lenv, junk, iopcl, nmts,
     $     nmtt, nmta, xa, xb, ya, yb, ipflag, allowp, jsym2e, nsymop,
     $     nop1, neqbas, neqshl, neqsh2, rotop, neqatm, accdes)

      implicit none


C +------------------------------------------------------------+
C |  erictr   --   CAJH, 07.2011                               |
C |                                                            |
C |                                   (based on PHFB version)  |
C |                                                            |
C |                                                            |
C |  This is the routine that calls PRISM (FoFDir), and        |
C |  should be called by formgg. It contracts P(g) with        |
C |  electron repulsion integrals.                             |
C |                                                            |
C |  The xa and xb arrays contain, on input, the alpha and     |
C |  beta density matrices to be contracted (they are          |
C |  DESTROYED). The ya and yb arrays are returned on output   |
C |  with the two-electron matrices. ALL MATRICES ARE STORED   |
C |  IN LOWER TRIANGULAR FORM.                                 |
C |                                                            |
C |  The integers nmts, nmtt, and nmta control the number of   |
C |  singlet, triplet, and antisymmetrices matrices,           |
C |  respectively. (Singlet matrices are contracted with       |
C |  Coulomb and exchange integrals, whereas triplet           |
C |  matrices are contracted only with exchange integrals).    |
C |                                                            |
C |  Note, however, that all matrices should be packed into    |
C |  the arrays xa, xb, with the following structure:          |
C |                                                            |
C |    . 1 ... nmts matrices should be singlet                 |
C |    . nmts+1 ... nmts+nmtt should be triplet                |
C |    . nmts+nmtt+1 ... nmts+nmtt+nmta should be antisymm     |
C |                                                            |
C |  All the logic of converting the output arrays ya and      |
C |  yb into G(g) should be done elsewhere.                    |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  A final cautionary note. The xa, xb, ya, and yb matrices  |
C |  are packed as a series of (nmts+nmtt+nmta) matrices,      |
C |  each of size NTT = nbasis*(nbasis+1)/2. However, the      |
C |  arrays xa, xb, ya, and yb should be allocated with size   |
C |  (nmts+nmtt+nmta)*NTT6D, where                             |
C |                                                            |
C |    NTT6D = nbas6d*(nbas6d+1)/2.                            |
C |                                                            |
C |  This is a PRISM requirement...                            |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       nbasis - number of basis functions (pure) !!!!!!junk
C       iopcl  - IOpCl (see FoFDir for details)
C       nmts   - number of singlet matrices to contract
C       nmts   - number of triplet matrices to contract
C       nmta   - number of antisymmetric matrices to contract
C       xa     - set of alpha density matrices
C       xb     - set of beta density matrices
C       ya     - set of alpha Hamiltonian matrices [ out ]
C       yb     - set of beta Hamiltonian matrices [ out ]

C     Include Gaussian common blocks 

#include "commonmol.inc"

C     Input variables

      real*8      v(*), accdes
      real*8      xa(*), xb(*), ya(*), yb(*)
      integer     iout, iprint, lenv
      integer     iopcl, nmts, nmtt, nmta

C       .. prism variables .. [ input variables ]

      real*8      rotop(*)
      integer     ipflag, nsymop, nop1, jsym2e
      integer     lseall, momega, nomega
      parameter   ( momega = DEFMOMEGA, nomega = DEFNOMEGA )
      real*8      omega(momega,nomega)
      integer     neqbas(nbasis,*), neqshl(*), neqsh2(*), neqatm(*)
      logical     allowp(50), fmm
      integer     jpflsta
      integer     ifmm, iseall
      integer     fmflag, fmflg1, nfxflg, ihmeth

C       .. PBC variables ..

      integer     npdir, jpbcsta, jpbcend
      real*8      celvec(3,3)

C     other variables

      real*8      junk
      integer     jend
      integer     irwmol, icntrl
      integer     ijunk
      integer     nmat, nmats, nmatt, nmatd
      logical     initf, dopurf, dopurp, keepp
      complex*16  one, cjunk

      ijunk = 0
      junk = 0.0d0
      jend = 1
      irwmol = 997
      one = cmplx(1.0d0,0.0d0)
      cjunk = cmplx(0.0d0,0.0d0)

C     set up PBC information
      
      call getpbc(npdir, celvec)
      jpbcsta = jend
      jpbcend = jpbcsta
      call rdipbc(1, npdir, 0, jpbcsta, jpbcend, v, lenv)
      jend = jpbcend

C     restore /mol/ common block

      call rwmol(2, irwmol)

C     set up PRISM control flags

      ifmm = 0
      iseall = 0

      jpflsta = jend

      call setpfl (iout, iprint, ifmm, ipflag, allowp, fmm, fmflag,
     $     fmflg1, nfxflg, ihmeth, omega, iseall,lseall,jpflsta,v(jend),
     $     lenv-jend+1)

C     Set up some options before calling FoFDir.

      icntrl = 0
      initf  = .true.
      dopurp = .true.
      dopurf = .true.
      keepp  = .false.

      nmat  = nmts + nmtt + nmta
      nmats = nmts
      nmatt = nmtt
      nmatd = nmat


C     Call FoFCou.

      call fofcou(iout,iprint,1,iopcl,icntrl,1,nomega,
     $  omega,fmm,fmflag,fmflg1,nfxflg,ipflag,allowp,lseall,lseall,
     $  1,natoms,keepp,initf,dopurp,dopurf,.false.,0,accdes,
     $  one,0,nmat,nmats,nmatt,nmatd,0,0,nbasis,
     $  jsym2e,nsymop,nop1,neqbas,neqshl,neqsh2,rotop,neqatm,0,0,xa,
     $  xb,junk,junk,junk,junk,junk,junk,junk,0,ya,yb,0,
     $  junk,junk,natoms,ian,atmchg,c,iattyp,iatfrg,junk,junk,
     $  v(jpbcsta),celvec,cjunk,cjunk,junk,0,0,junk,junk,ijunk,
     $  ijunk,junk,v(jend),v(jend),lenv-jend+1)

      return
      end


