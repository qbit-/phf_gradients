

      subroutine derictr (iout, iprint, v, lenv, iopcl, nmts,
     $     nmtt, nmta, pa, pb, fa, fb, ipflag, fmm, fmflag, fmflg1, 
     $     nfxflg, lseall, nomega, omega, allowp, jsym2e,
     $     nopuse, nop1, neq, neqsh, neqs2, rotall, neqall, pbc, 
     $     celvec, accdes)

      implicit none

C +------------------------------------------------------------+
C |  derctr   --   RSS, 07.2013                                |
C |                                                            |
C |                                   (based on PHFB version)  |
C |                                                            |
C |                                                            |
C |  This is the routine that calls PRISM (FofCou), and        |
C |  should be called by evaldg. It contracts P(g) with        |
C |  derivatives of electron repulsion integrals to form       |
C |  forces.                                                   |
C |                                                            |
C |  The pa and pb arrays contain, on input, the alpha and     |
C |  beta density matrices to be contracted (they are          |
C |  DESTROYED). The fa and fb arrays are returned on output   |
C |  with the 3*natoms two-electron matrices. ALL MATRICES     |
C |  ARE STORED IN LOWER TRIANGULAR FORM.                      |
C |                                                            |
C |  The integers nmts, nmtt, and nmta control the number of   |
C |  singlet, triplet, and antisymmetrices matrices,           |
C |  respectively. (Singlet matrices are contracted with       |
C |  Coulomb and exchange integrals, whereas triplet           |
C |  matrices are contracted only with exchange integrals).    |
C |                                                            |
C |  Note, however, that all matrices should be packed into    |
C |  the arrays pa, pb, with the following structure:          |
C |                                                            |
C |    . 1 ... nmts matrices should be singlet                 |
C |    . nmts+1 ... nmts+nmtt should be triplet                |
C |    . nmts+nmtt+1 ... nmts+nmtt+nmta should be antisymm     |
C |                                                            |
C |  All the logic of converting the output arrays fa and      |
C |  fb into dG(g)/da should be done elsewhere.                |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  A final cautionary note. The pa, pb, fa, and fb matrices  |
C |  are packed as a series of (nmts+nmtt+nmta) matrices,      |
C |  each of size NTT = nbasis*(nbasis+1)/2. However, the      |
C |  arrays pa, pb, fa, and fb should be allocated with size   |
C |  (nmts+nmtt+nmta)*NTT6D, where                             |
C |                                                            |
C |    NTT6D = nbas6d*(nbas6d+1)/2.                            |
C |                                                            |
C |  This is a PRISM requirement...                            |
C |                                                            |
C +------------------------------------------------------------+

C     Include Gaussian common blocks 

#include "commonmol.inc"


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iopcl  - IOpCl (see FoFDir for details)
C       nmts   - number of singlet matrices to contract
C       nmts   - number of triplet matrices to contract
C       nmta   - number of antisymmetric matrices to contract
C       pa     - set of alpha density matrices
C       pb     - set of beta density matrices
C       fa     - set of alpha Hamiltonian matrices [ out ]
C       fb     - set of beta Hamiltonian matrices [ out ]

C     Input variables

      real*8      v(*), accdes
      real*8      pa(*), pb(*), fa(*), fb(*)
      integer     iout, iprint, lenv
      integer     iopcl, jsym2e, nmts, nmtt, nmta
      logical     fmm

C       .. prism variables [ input variables ]..

      integer     ipflag, fmflag, fmflg1, nfxflg
      integer     lseall, nomega
      logical     allowp(*)
      real*8      omega(*)

C       .. symmetry related variables [ input variables ]..

      integer     neq(*), neqsh(*), neqs2(*), neqall(*)
      real*8      rotall(*)
      integer     nopuse, nop1

C      .. PBC related variables ..
      
      real*8      pbc(*), celvec(*)

C     Other variables

      integer     icntrl
      integer     ijunk
      integer     nmats, nmatt, nmat, nmatd
      logical     keepp, dopurp, dopurf, initf

C     Constants

      real*8     junk
      complex*16 zero, one, cjunk

      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)

      ijunk = 0
      junk  = 0.0d0
      cjunk = zero

C     Call FofCou..

      icntrl =  3721
      keepp  = .false.
      dopurp = .true.
      dopurf = .true.
      initf  = .true.
      nmats  = nmts
      nmatt  = nmtt
      nmat   = nmts + nmtt + nmta
      nmatd  = nmat
          
      call fofcou(iout,iprint,1,iopcl,icntrl,1,nomega,
     $  omega,fmm,fmflag,fmflg1,nfxflg,ipflag,allowp,lseall,lseall,
     $  1,natoms,keepp,initf,dopurp,dopurf,.false.,0,accdes,
     $  one,0,nmat,nmats,nmatt,nmatd,0,0,nbasis,
     $  jsym2e,nopuse,nop1,neq,neqsh,neqs2,rotall,neqall,0,0,pa,
     $  pb,junk,junk,junk,junk,junk,junk,junk,0,junk,junk,0,
     $  fa,fb,natoms,ian,atmchg,c,iattyp,iatfrg,junk,junk,pbc,
     $  celvec,cjunk,cjunk,junk,0,0,junk,junk,ijunk,
     $  ijunk,junk,v(1),v(1),lenv)

      return
      end


