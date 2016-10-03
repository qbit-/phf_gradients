      

      subroutine evaldg (iout, iprint, v, lenv, itype, ncik, 
     $     ngrdt, ifmm, iseall, jsym2e, accdes, irwneq, irwneqshl, 
     $     irwneqsh2, irwmol, irwt, irwd)

      implicit none

C +------------------------------------------------------------+
C |  evaldg   --   RSS, 07.2013                                |
C |                                                            |
C |                                   (based on CAJH's formgg) |
C |                                                            |
C |                                                            |
C |  Form the matrices dG(g)/da by contracting P(g) with       |
C |  derivatives of electron-repulsion integrals               |
C |  (Coulomb and exchange), then contract them with P(g)      |
C |  to form a two-body contribution to forces d(g)            |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  On input, P(g) are stored in AO basis                     |
C |  Before sending them to PRISM,                             |
C |  we need to decompose them into simmetric/antisymmetric/   |
C |  triplet parts in PRISM terminology                        |
C |                                                            |
C |  The matrices P(g) are prepared in batches, then we call   |
C |  PRISM. The number of density matrices                     |
C |  sent per batch is decided here based on the amount of     |
C |  memory available (see below).                             |
C |                                                            |
C |                                                            |    
C +------------------------------------------------------------+
C     include Gaussian common blocks

#include "commonmol.inc"
#include "repall.inc"

C     input / output variables

C       v       - scratch array
C       lenv    - size of scratch array
C       itype   - decides which spin blocks of transition density
C                matrices are stored
C       ncik    - number of basis states in complex conjugation CI problem
C       ngrdt   - total number of grid points
C       ifmm    - PRISM control flag
C       iseall  - PRISM control flag
C       jsym2e  - controls the use of symmetry in atomic integral evaluation in PRISM
C       accdes  - desired accuracy of atomic integrals in PRISM
C       irwneq  - RWF file with shell equvalence data
C       irwneqshl - RWF file with shell equvalence data
C       irwneqsh2 - RWF file with shell equvalence data
C       irwmol  - Gaussian /mol/ common block
C       irwt    - RWF file with transition density matrices P(g)
C       irwd    - RWF file with forces d(g) [ out ]

      real*8      v(*), accdes
      integer     iout, iprint, lenv, itype, ncik
      integer     ifmm, iseall, jsym2e
      integer     ngrdt
      integer     irwneq, irwneqshl, irwneqsh2, irwmol 
      integer     irwt, irwd

C       .. prism variables ..

      integer     ipflag, fmflag, fmflg1, nfxflg, ihmeth
      integer     lseall, momega, nomega
      logical     allowp(50), fmm
      parameter   ( momega = DEFMOMEGA, nomega = DEFNOMEGA )
      real*8      omega(momega,nomega)
      integer     jpflsta

C       .. PBC variables ..

      integer     npdir, jpbcsta, jpbcend
      real*8      celvec(3,3)

C       .. symmetry related variables ..

      integer     nop1, nop2, nopuse
      integer     lenneq, lennes, lenne2
      integer     jneq, jneqsh, jneqs2

C     other variables

      real*8      junk 
      integer     jend, mdv
      integer     iprtf, iopcl, nbas6d  
      integer     ntt, ntt6d, nbsq, szmat, nat3
      integer     lenp, lenio, npmem, npmax, npuse, ipar
      integer     mblklen, blklen
      integer     szbfmat, szbfgg, szbfntt, szpgin, szdgout
      integer     szpas, szpat, szpaa, szpbs, szpbt, szpba
      integer     szfas, szfat, szfaa, szfbs, szfbt, szfba
      integer     szextrap, szextraf, szextra
      integer     lscr1, lscr2, lscr3, lscr 
      integer     jpgin, jdgout, jbfmat, jbfgg, jbfntt, jscr 
      integer     jpas, jpat, jpaa, jfas, jfat, jfaa
      integer     jpbs, jpbt, jpba, jfbs, jfbt, jfba
      integer     ipgin, idgout, ibfmat, ibfgg, ibfntt, iscr
      integer     ipas, ipat, ipaa, ifas, ifat, ifaa
      integer     ipbs, ipbt, ipba, ifbs, ifbt, ifba
      integer     iscrrs, iscris, iscrra, iscria
      integer     indrs, indra, indis, india
      integer     izaa, izbb, izba, izab
      integer     nmaxs, nmaxt, nmaxa, szg
      integer     ngrdtt, iblock
      integer     ip, iat, igrdt, tid, ii, jj
      integer     inddg
      logical     usemt, dumpmt, dumpfk

      complex*16  dg

C     functions 
      
      integer     intowp, itqry, nproc, lappar
      integer     omp_get_thread_num  
      complex*16  trcabc

C     constants

      complex*16  zero, two, half
      
      zero = cmplx (0.0d0, 0.0d0)
      half = cmplx (0.5d0, 0.0d0)

C     %%%%%%%%%%%%%%%%%%%%%%
C     %   Initialization   %
C     %%%%%%%%%%%%%%%%%%%%%%

C     Initialize data structures to call PRISM...

C     Specify iopcl.

C       iopcl = 0,  when only alpha matrices are sent to PRISM
C             = 1,  when both alpha and beta matrices are sent to PRISM

      if ( itype .eq. 1 ) then
        iopcl = 0
      elseif ( itype .eq. 2 .or. itype .eq. 3 ) then
        iopcl = 1
      endif
   
C     restore /mol/ common block

      call rwmol(2, irwmol)
      
C     retrieve PBC information

      jend = 1

      call getpbc(npdir, celvec)
      jpbcsta = jend
      jpbcend = jpbcsta
      call rdipbc(1, npdir, 0, jpbcsta, jpbcend, v, lenv)
      jend = jpbcend

C     set PRISM control flags

      jpflsta = jend
      call setpfl (iout, iprint, ifmm, ipflag, allowp, fmm, fmflag,
     $     fmflg1, nfxflg, ihmeth, omega, iseall,lseall,jpflsta,v(jend),
     $     lenv-jend+1)

      jend = jpflsta + jend

C     Retrieve symmetry information 

C     Get number of symmetry operations.
C       - nop1,    order of the concise abelian subgroup
C       - nop2,    order of the full abelian subgroup
C       - nopuse,  actual oder to be used

      call getnop (nop1, nop2)
      call pccki (0, junk, neqall, mxatso, nopall, natoms, nopall)

      nopuse = nopall

C     Set up arrays with symmetry information

      lenneq = 0
      lennes = 0
      lenne2 = 0

      if ( jsym2e .eq. 1 ) then
        lenneq = intowp (nop1*nbasis)
        lennes = max (itqry (irwneqshl), 0)
        lenne2 = max (itqry (irwneqsh2), 0)
      endif

      jneq   = jend
      jneqsh = jneq   + lenneq
      jneqs2 = jneqsh + lennes
      jend   = jneqs2 + lenne2

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'evaldg')

C     Fill arrays with symmetry related quantities.

      if ( lenneq .ne. 0 ) then
        call fileio (2, -irwneq, lenneq, v(jneq), 0)
      endif

      if ( lennes .ne. 0 ) then
        call fileio (2, -irwneqshl, lennes, v(jneqsh), 0)
      endif

      if ( lenne2 .ne. 0 ) then
        call fileio (2, -irwneqsh2, lenne2, v(jneqs2), 0)
      endif
      
C     size of basis in cartesian functions

      call getnb6 (nbas6d)

C     Construct useful variables.

      ntt   = nbasis * (nbasis + 1) / 2
      ntt6d = nbas6d * (nbas6d + 1) / 2

      nbsq = nbasis * nbasis

      if ( itype .eq. 1 ) then
         szmat = 2*nbsq
      elseif ( itype .eq. 2 ) then
         szmat = 2*2*nbsq
      elseif ( itype .eq. 3 ) then
         szmat = 2*2*2*nbsq
      endif

      nat3 = natoms*3
 
C     Printing level for PRISM...

      iprtf = iprint - 2
      if ( iprtf .lt. 0 ) iprtf = 0

C     %%%%%%%%%%%%%%%%%%%%%%%
C     %  Memory allocation. %
C     %%%%%%%%%%%%%%%%%%%%%%%

C       Any matrix can be decomposed to two symmetric and two antisymmetric
C       matrices (real and imaginary).
C       Here we specify the maximum number of singlet, triplet, and antisymmetric
C       blocks (per grid point) for each itype case.

C       The spin blocks of P(g) can be decomposed into:

C              spin block     # singlet   # triplet   # antisymm
C              -------------------------------------------------
C                [ aa ]           2                        2
C                [ ab ]                        2           2
C                [ ba ]                        2           2
C                [ bb ]           2                        2

C         ** The [aa] and [bb] blocks are 'singlets' because they need to
C            be contracted with Coulomb and exchange integrals, while the
C            [ab] and [ba] blocks only need to be contracted with exchange
C            integrals. All antisymmetric matrices are considered as
C            antisymmetric.

C       Furthermore, for itype = 2, 3, we will send both alpha and beta
C       matrices, while for itype = 1 we will only send alpha matrices:

C         [ aa ] block is sent as alpha
C         [ ab ] block is sent as alpha
C         [ bb ] block is sent as beta
C         [ ba ] block is sent as beta


      if ( itype .eq. 1 .or. itype .eq. 2 ) then
        nmaxs = 2
        nmaxt = 0
        nmaxa = 2

      elseif ( itype .eq. 3 ) then
        nmaxs = 2
        nmaxt = 2
        nmaxa = 4
      endif

C     Prepare memory for parallel execution
C     =====================================

C       - lenp,  total length of scratch space per processor
C       - lenio, length of the total output per one grid point
C

C     Scratch space per process (lenp):
C     ==================================
C       - szbfpaa, length of buffers for a matrix P(g) decomposed to spin blocks
C       - four scratch lower triangular arrays to decompose transition 
C                  density matrix before feeding it to PRISM
      
      szbfmat = szmat
      szbfgg  = szmat
      szbfntt = 4 * ntt

C       - scratch array to call dmblck, and sptblk

C       - lscr1, length of scratch space required for sptblk, mode = 1 (split)
C       - lscr2, length of scratch space required for sptblk, mode = 2 (combine)
C       - lscr3, length of scratch space required for dmblck (only if all spin
C                blocks are active) 
C       - lscr is the total length of scratch space

      lscr1 = 2*nbsq + 4*ntt
      lscr2 = 3*nbsq

      lscr  = max (lscr1, lscr2)

      if ( itype .eq. 3 ) then
        lscr3 = 16*nbsq
        lscr = max (lscr, lscr3)
      endif

C     Total scratch:

      lenp = szbfmat + szbfgg + szbfntt + lscr

C     Space needed for IO to/from parallel region per grid point (lenio):
C     ===================================================================

C       - szpgin  - buffer length for transition density matrix
C       - szpXX, szfXX  - buffer length for decomposed density matrices
C                         and Fock matrices
C           Naming convention:
C             szpas - density, alpha - singlet
C             szpat - density, alpha - triplet
C             szpaa - density, alpha - antisymmetric
C             ...

C       ** Note that PRISM requires that all matrices be allocated ntt6d,
C          as opposed to ntt. However, all matrices (singlet, triplet,
C          antisymmetric) are sent as a single block. That's why
C          allocation may look tricky.

C       - szdgout - buffer for forces d(g) at current grid point (complex)
C       - szextra - an extra amount we need to provide for PRISM to do
C          sperical/cartesian conversion of matrices. This is a PRISM
C          requirement. Each lower triangular matrix is converted. 
      szpgin  = szmat
      szpas   = ntt * nmaxs
      szpat   = ntt * nmaxt
      szpaa   = ntt * nmaxa
      
      if ( itype .eq. 2 .or. itype .eq. 3 ) then
         szpbs   = szpas
         szpbt   = szpat
         szpba   = szpaa
      else
         szpbs = 0
         szpbt = 0
         szpba = 0
      endif

      szfas   = ntt * nmaxs * nat3
      szfat   = ntt * nmaxt * nat3
      szfaa   = ntt * nmaxa * nat3

      if ( itype .eq. 2 .or. itype .eq. 3 ) then
         szfbs   = szfas
         szfbt   = szfat
         szfba   = szfaa
      else
         szfbs = 0
         szfbt = 0
         szfba = 0
      endif

      szextrap = (ntt6d - ntt)*(nmaxs+nmaxt+nmaxa)*(nat3 + 1)
      szextraf = (ntt6d - ntt)*(nmaxs+nmaxt+nmaxa)*nat3
      szextra  = szextrap + szextraf

      if ( itype .eq. 2 .or. itype .eq. 3 ) then
         szextra = 2*szextra
      endif

      szdgout = 2 * nat3

C     Total IO:
      
      lenio = szpgin + szpas + szpat + szpaa + szpbs + szpbt + 
     $     + szpba + szfas + szfat + szfaa + szfbs + szfbt + szfba +
     $     szdgout + szextra

C     ...define how many threads with 1 block of input we can run
C
C       According to our experience with PRISM, it 
C       may internally make a copy of the input array of matrices.
C       We need to provide to PRISM at least twice the amount 
C       of memory we used for sending density matrices. 

C       This is the upper bound of the number of threads

      npmem = (lenv - jend + 1) / (3*(lenp + lenio))

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in evaldg for a single CPU.')
      endif

C     Total number of points
      
      ngrdtt = ngrdt * ncik * ncik
      
C     ...define how many threads are allowed
C     ...check if LAPACK is allowed to run in parallel

      ipar  = lappar (0)
      if ( ipar .eq. 1 ) then
         npmax = nproc (0)
      else
         npmax = 1
      endif

C     In case of parallel execution find the number of blocks per thread
C     we can afford.
C     We substract the length of scratch from the total amount.
C     Then we divide by the total number of threads as each thread 
C     replicates the input array of matrices. 

      if ( npmem .gt. npmax ) then
         npuse   = npmax
         mblklen = min (ngrdtt,(lenv - jend + 1 - lenp*npuse) / 
     $        ((npuse+1)*lenio))
      else
         npuse   = npmem
         mblklen = npmem
      endif

C     Allocate memory for the scratch per thread:

C     SCRATCH

      jbfmat  = jend
      jbfgg   = jbfmat + npuse*szbfmat 
      jbfntt  = jbfgg  + npuse*szbfgg
      jscr    = jbfntt + npuse*szbfntt
      jend    = jscr   + npuse*lscr

C     Allocate memory for the I/O scratch per block:
C
C       - jpgin   - transitioin density matrices
C       - jpXX, jfXX  - arrays storing all the matrices to send to PRISM:
C             p - density matrices (sent to PRISM)
C             f - Fock matrices (returned from PRISM)

C IO
      jdgout = jend
      jend   = jdgout + mblklen*szdgout

      jpgin  = jend 
      jpas   = jpgin  + mblklen*szpgin
      jpat   = jpas   + mblklen*szpas
      jpaa   = jpat   + mblklen*szpat
      jend   = jpaa   + mblklen*szpaa + szextrap

      call aclear (mblklen*(szpas+szpat+szpaa+szextrap), v(jpas))

      if ( itype .eq. 2 .or. itype .eq. 3 ) then
         jpbs = jend
         jpbt = jpbs + mblklen*szpbs
         jpba = jpbt + mblklen*szpbt
         jend = jpba + mblklen*szpba + szextrap
         call aclear (mblklen*(szpbs+szpbt+szpba+szextrap), v(jpbs))

      endif

      jfas   = jend
      jfat   = jfas  + mblklen*szfas
      jfaa   = jfat  + mblklen*szfat
      jend   = jfaa  + mblklen*szfaa + szextraf

      if ( itype .eq. 2 .or. itype .eq. 3 ) then
         jfbs = jend
         jfbt = jfbs + mblklen*szfbs
         jfba = jfbt + mblklen*szfbt
         jend = jfba + mblklen*szfba + szextraf
      endif
C     debug
      jend = jend

      call tstcor (jend - 1, lenv, 'evaldg')

      mdv = lenv - jend + 1

C     %%%%%%%%%%%%%%%
C     %  Form d(g)  %
C     %%%%%%%%%%%%%%%

C     Contract P(g) with two-electron integrals...

C     Loop over batches.

      usemt = .false.
      dumpmt = .false.
      do 101 iblock = 1, ngrdtt, mblklen

        blklen = min (mblklen, ngrdtt - iblock + 1)

C     Load  P(g) in blocks

        if (.not. usemt) then
           call fileio (2, -irwt, blklen*szpgin, v(jpgin),
     $          (iblock-1)*szpgin)
        else
           open(111,file='/home/shiva/git/ela/r/denin.txt')
           do 115 ii = 1,szpgin*blklen
              read(111,*) v(jpgin + ii-1)
 115       continue
           close(111)
        endif

        if (dumpmt) then
           open(111,file='/home/shiva/git/ela/r/denout.txt')
           do 116 ii = 1,szpgin*blklen
              write(111,*) v(jpgin + ii-1)
 116       continue
           close(111)
        endif


C     Computation loop. 

CC$omp parallel default(shared) 
CC$omp+private(ip, tid, ipgin, ipas, ipat, ipaa)
CC$omp+private(ipbs, ipbt, ipba, ibfmat, ibfntt)
CC$omp+private(iscr, iscrrs, iscrra, iscris, iscria)
CC$omp+private(izaa, izbb, izab, izba)
CC$omp+private(indrs, indra, indis, india)
CC$omp+private(iat, ifas, ifat, ifaa, ifbs, ifbt)
CC$omp+private(ifba, idgout, ibfgg, inddg, dg)

        tid = omp_get_thread_num()

C     Compute pointers to private scratch arrays per thread:
        
        ibfmat  = jbfmat + tid*szbfmat
        ibfgg   = jbfgg  + tid*szbfgg
        ibfntt  = jbfntt + tid*szbfntt
        iscr    = jscr   + tid*lscr
        
C     Compute pointers inside the bfntt array (scratch for lower triangular matrices)
        
        iscrrs = ibfntt
        iscrra = ibfntt + ntt
        iscris = ibfntt + 2*ntt
        iscria = ibfntt + 3*ntt
           
C     Compute pointers to various spin blocks of the array bfgg: 
        izaa = ibfmat
        
        if ( itype .eq. 2 ) then
           izbb = ibfmat + 2*nbsq
        elseif ( itype .eq. 3 ) then
           izab = ibfmat + 2*nbsq
           izba = ibfmat + 4*nbsq
           izbb = ibfmat + 6*nbsq
        endif
C     debug

C        write(iout,'(A,I5,I5,I5,I5,I5,I5)') 
C     & '1: ip, tid, ibfmat, ibfgg, ibfntt, iscr',
C     &  ip, tid, ibfmat, ibfgg, ibfntt, iscr

C        write(iout,'(A,I3,I5,I5,I5,I5)') 
C     & '1: tid, iscrrs, iscrra, iscris, iscria', tid,
C     &       iscrrs, iscrra, iscris, iscria

CC$omp do schedule(static,1)
        do 102 ip = 1, blklen

C     Compute pointers to the IO arrays
C     ..Decomposed density matrices..
           ipgin  = jpgin + (ip - 1)*szpgin
           ipas   = jpas  + (ip - 1)*szpas
           ipat   = jpat  + (ip - 1)*szpat
           ipaa   = jpaa  + (ip - 1)*szpaa
              
           if ( itype .eq. 2 .or. itype .eq. 3 ) then
              ipbs = jpbs + (ip - 1)*szpbs
              ipbt = jpbt + (ip - 1)*szpbt
              ipba = jpba + (ip - 1)*szpba
           endif

C     debug
C     omp single
C        write(iout,'(A,I3,I4,I5,I5,I5)')
C     &  '1: ip, tid, ipgin, ipas, ipat, ipaa', ip, tid, ipgin, 
C     &       ipas, ipat, ipaa 
C     omp end single
C        write(iout,'(A,I3,I5,I5,I5)')
C     &  '1: tid, ipbs, ipbt, ipba', tid, 
C     &       ipbs, ipbt, ipba 

C     Split P(g) at current grid point into spin-blocks.
C     Here we define pointers inside the array bfpg

          if ( itype .eq. 1 ) then
             call amove (2*nbsq, v(ipgin), v(izaa))

          elseif ( itype .eq. 2 ) then
             call amove (2*nbsq, v(ipgin), v(izaa))
             call amove (2*nbsq, v(ipgin+2*nbsq),v(izbb))
             
          elseif ( itype .eq. 3 ) then

C             .. split into spin blocks (in square format) ..
            call dmblck (iout, iprint, v(iscr), lscr, nbasis, v(ipgin),
     $           v(izaa), v(izab), v(izba), v(izbb), 2, -1)
          endif

C     Split each block into real-imaginary, symmetric-antisymmetric.
C     Place the result into corresponding places in the IO arrays
C     Each thread works on array with pointer ipXX
C           .. [alpha, alpha] block

          if ( itype .lt. 4 ) then
            indrs = ipas
            indis = indrs + ntt
            indra = ipaa
            india = indra + ntt

            call sptblk (iout, iprint, v(iscr), lscr, nbasis, v(izaa),
     $           v(iscrrs), v(iscrra), v(iscris), v(iscria), 1)

            call amove (ntt, v(iscrrs), v(indrs))
            call amove (ntt, v(iscrra), v(indra))
            call amove (ntt, v(iscris), v(indis))
            call amove (ntt, v(iscria), v(india))

          endif

C           .. [beta, beta] block

          if ( itype .gt. 1 ) then
            indrs = ipbs
            indis = indrs + ntt
            indra = ipba
            india = indra + ntt

            call sptblk (iout, iprint, v(iscr), lscr, nbasis, v(izbb),
     $           v(iscrrs), v(iscrra), v(iscris), v(iscria), 1)

            call amove (ntt, v(iscrrs), v(indrs))
            call amove (ntt, v(iscrra), v(indra))
            call amove (ntt, v(iscris), v(indis))
            call amove (ntt, v(iscria), v(india))
        
          endif

C           .. [alpha, beta] block

          if ( itype .eq. 3 ) then
            indrs = ipat
            indis = indrs + ntt
            indra = ipaa  + 2 * ntt
            india = indra + ntt

            call sptblk (iout, iprint, v(iscr), lscr, nbasis, v(izab),
     $           v(iscrrs), v(iscrra), v(iscris), v(iscria), 1)

            call amove (ntt, v(iscrrs), v(indrs))
            call amove (ntt, v(iscrra), v(indra))
            call amove (ntt, v(iscris), v(indis))
            call amove (ntt, v(iscria), v(india))

          endif

C           .. [beta, alpha] block

          if ( itype .eq. 3 ) then
            indrs = ipbt
            indis = indrs + ntt
            indra = ipba  + 2 * ntt
            india = indra + ntt

            call sptblk (iout, iprint, v(iscr), lscr, nbasis, v(izba),
     $           v(iscrrs), v(iscrra), v(iscris), v(iscria), 1)

            call amove (ntt, v(iscrrs), v(indrs))
            call amove (ntt, v(iscrra), v(indra))
            call amove (ntt, v(iscris), v(indis))
            call amove (ntt, v(iscria), v(india))
        
          endif

 102   continue

CC$omp end do
CC$omp end parallel
          
C     Contract P(g) against two-electron integrals in erictr.

       call derictr (iout, iprtf, v(jend), lenv-jend+1, iopcl,
     $ nmaxs*blklen, nmaxt*blklen, nmaxa*blklen, v(jpas),v(jpbs),
     $ v(jfas), v(jfbs), ipflag, fmm, fmflag, fmflg1, nfxflg, lseall,
     $ nomega, omega, allowp, jsym2e, nopuse, nop1, v(jneq), 
     $ v(jneqsh), v(jneqs2), rotall, neqall, v(jpbcsta), celvec,
     $ accdes)

CC$omp parallel default(shared) 
CC$omp+private(ip, tid, ipgin, ipas, ipat, ipaa)
CC$omp+private(ipbs, ipbt, ipba, ibfmat, ibfntt)
CC$omp+private(iscr, iscrrs, iscrra, iscris, iscria)
CC$omp+private(izaa, izbb, izab, izba)
CC$omp+private(indrs, indra, indis, india)
CC$omp+private(iat, ifas, ifat, ifaa, ifbs, ifbt)
CC$omp+private(ifba, idgout, ibfgg, inddg, dg)

        tid = omp_get_thread_num()

C     Compute pointers to private scratch arrays per thread:
        
        ibfmat  = jbfmat + tid*szbfmat
        ibfgg   = jbfgg  + tid*szbfgg
        ibfntt  = jbfntt + tid*szbfntt
        iscr    = jscr   + tid*lscr
        
C     Compute pointers inside the bfntt array (scratch for lower triangular matrices)
        
        iscrrs = ibfntt
        iscrra = ibfntt + ntt
        iscris = ibfntt + 2*ntt
        iscria = ibfntt + 3*ntt

C     debug
C        write(iout,'(A,I5,I5,I5,I5,I5,I5)')
C     & '2: tid, ip, ibfmat,ibfgg,ibfntt, iscr', tid, ip, 
C     &       ibfmat, ibfgg, ibfntt, iscr

C        write(iout,'(A,I5,I5,I5,I5,I5)')
C     & '2: tid, iscrrs, iscrra, iscris, iscria', tid,
C     &       iscrrs, iscrra, iscris, iscria
           
C     Compute pointers to various spin blocks of the array bfgg: 
        izaa = ibfmat
        
        if ( itype .eq. 2 ) then
           izbb = ibfmat + 2*nbsq
        elseif ( itype .eq. 3 ) then
           izab = ibfmat + 2*nbsq
           izba = ibfmat + 4*nbsq
           izbb = ibfmat + 6*nbsq
        endif

CC$omp do schedule(static,1)
C     Collect  Fock matrices and compute forces d(g).

        do 103 ip = 1, blklen

C     For given density matrix (lower triangular) PRISM returns  
C     3*natoms derivative Fock matrices consequitively, hence 
C     another loop here

            do 104 iat = 1, nat3

C     Compute pointers to the IO arrays
C     ..Density matrices loaded.. 
              ipgin  = jpgin + (ip - 1)*szpgin

C     ..Ponters to decomposed Fock matrices..
              ifas   = jfas + (ip - 1)*szfas + (iat - 1)*ntt
              ifat   = jfat + (ip - 1)*szfat + (iat - 1)*ntt
              ifaa   = jfaa + (ip - 1)*szfaa + (iat - 1)*ntt
           
              if ( itype .eq. 2 .or. itype .eq. 3 ) then
                 ifbs = jfbs + (ip - 1)*szfbs + (iat - 1)*ntt
                 ifbt = jfbt + (ip - 1)*szfbt + (iat - 1)*ntt
                 ifba = jfba + (ip - 1)*szfba + (iat - 1)*ntt
              endif

              idgout   = jdgout + (ip - 1)*szdgout + (iat - 1)*2
C     debug
C              write(iout,'(A,I3,I4,I5,I5,I5)')
C     & '2: tid, ipgin, ifas, ifat, ifaa',tid,ipgin, 
C     &             ifas, ifat, ifaa 
C              write(iout,'(A,I3,I5,I5,I5,I5)')
C     & '2: tid, ifbs, ifbt,ifba,idgout',tid,
C     &             ifbs, ifbt, ifba, idgout

C     Merge the real-imaginary, symmetric-antisymmetric blocks into a
C     single spin-block.
              
C     .. [alpha, alpha] block
              
              if ( itype .lt. 4 ) then
                 indrs = ifas
                 indis = indrs + ntt * nat3
                 indra = ifaa 
                 india = indra + ntt * nat3
                                  
                 call amove (ntt, v(indrs), v(iscrrs))
                 call amove (ntt, v(indra), v(iscrra))
                 call amove (ntt, v(indis), v(iscris))
                 call amove (ntt, v(india), v(iscria))

C     .. scale the antisymmetric blocks by a factor of -1 ..
                 
C                 call ascale (ntt, -1.0d0, v(iscrra), v(iscrra))
C                 call ascale (ntt, -1.0d0, v(iscria), v(iscria))
                 
                 call sptblk (iout,iprint,v(iscr),lscr,nbasis,v(izaa),
     $                v(iscrrs), v(iscrra), v(iscris), v(iscria), 2)
              
              endif

C     .. [beta, beta] block
              
              if ( itype .gt. 1 ) then
                 indrs = ifbs
                 indis = indrs + ntt * nat3
                 indra = ifba 
                 india = indra + ntt * nat3
                 
                 call amove (ntt, v(indrs), v(iscrrs))
                 call amove (ntt, v(indra), v(iscrra))
                 call amove (ntt, v(indis), v(iscris))
                 call amove (ntt, v(india), v(iscria))

C     .. scale the antisymmetric blocks by a factor of -1 ..
                 
C                 call ascale (ntt, -1.0d0, v(iscrra), v(iscrra))
C                 call ascale (ntt, -1.0d0, v(iscria), v(iscria))
                 
                 call sptblk (iout,iprint,v(iscr),lscr,nbasis,v(izbb),
     $                v(iscrrs), v(iscrra), v(iscris), v(iscria), 2)
                 
              endif

C     .. [alpha, beta] block

              if ( itype .eq. 3 ) then
                 indrs = ifat
                 indis = indrs + ntt * nat3
                 indra = ifaa  + 2 * ntt * nat3
                 india = indra + ntt * nat3
                                  
                 call amove (ntt, v(indrs), v(iscrrs))
                 call amove (ntt, v(indra), v(iscrra))
                 call amove (ntt, v(indis), v(iscris))
                 call amove (ntt, v(india), v(iscria))

C     .. scale the antisymmetric blocks by a factor of -1 ..
                 
C                 call ascale (ntt, -1.0d0, v(iscrra), v(iscrra))
C                 call ascale (ntt, -1.0d0, v(iscria), v(iscria))
                 
                 call sptblk (iout,iprint,v(iscr),lscr,nbasis,v(izab),
     $                v(iscrrs),v(iscrra),v(iscris),v(iscria), 2)
              endif

C     .. [beta, alpha] block

              if ( itype .eq. 3 ) then
                 indrs = ifbt
                 indis = indrs + ntt * nat3
                 indra = ifba  + 2 * ntt * nat3
                 india = indra + ntt * nat3

                 call amove (ntt, v(indrs), v(iscrrs))
                 call amove (ntt, v(indra), v(iscrra))
                 call amove (ntt, v(indis), v(iscris))
                 call amove (ntt, v(india), v(iscria))

C     .. scale the antisymmetric blocks by a factor of -1 ..
                 
C                 call ascale (ntt, -1.0d0, v(iscrra), v(iscrra))
C                 call ascale (ntt, -1.0d0, v(iscria), v(iscria))
        
                 call sptblk (iout,iprint,v(iscr),lscr,nbasis,v(izba),
     $                v(iscrrs), v(iscrra), v(iscris), v(iscria), 2)

              endif
              
C     Build dG(g)/da from it's spin blocks.

              if ( itype .eq. 1 ) then
                 call amove (2*nbsq, v(izaa), v(ibfgg))
              elseif ( itype .eq. 2 ) then
                 call amove (2*nbsq, v(izaa), v(ibfgg))
                 call amove (2*nbsq, v(izbb), v(ibfgg + 2*nbsq))
              elseif ( itype .eq. 3 ) then
                 call dmblck (iout, iprint, v(iscr), lscr, nbasis,
     $                v(ibfgg), v(izaa), v(izab), v(izba),
     $                v(izbb), 1, -1)
              endif

              dumpfk = .false.
              if (dumpfk) then
                 open(111,file='/home/shiva/git/ela/r/fkout.txt')
                 do 117 ii = 1,szmat
                    write(111,*) v(jbfgg + ii-1)
 117             continue
                 close(111)
              endif
              
C     Contract with P(g)
C     Warning! trcabc does tr[ A.B* ], not tr[ A.B ] and not tr[ A*B ] as
C     it is said in the code, so we need to get G*

              call vconjg (szmat,v(ibfgg))

              if (usemt) then
                 open(111,file='/home/shiva/git/ela/r/denin.txt')
                 do 125 ii = 1,2*nbsq
                    read(111,*) v(ipgin + ii-1)
 125             continue
                 close(111)

                 open(111,file='/home/shiva/git/ela/r/gamin.txt')
                 do 126 ii = 1,2*nbsq
                    read(111,*) v(ibfgg + ii-1)
 126             continue
                 close(111)
              endif

              if ( itype .eq. 1 ) then
                 dg  = trcabc(1, nbasis, v(ipgin), v(ibfgg) )

              elseif ( itype .eq. 2 ) then
                 dg  = half * trcabc(2, nbasis, v(ipgin), v(ibfgg))
                 
              elseif ( itype .eq. 3 ) then
                 dg = half * trcabc(1, 2*nbasis, v(ipgin), v(ibfgg))
                 
              endif

C     Save forces into array

              v(idgout)     =  dble (dg)
              v(idgout + 1) =  aimag (dg)
C     debug
C              write(iout,*)
C     & '2: tid, dg, ip, iat',tid, dg, ip, iat
 104     continue
 103    continue
CC$omp end do
CC$omp end parallel
C     Save forces d(g) into file

        inddg = (iblock - 1)*szdgout
        call fileio (1, -irwd, blklen*szdgout, v(jdgout), inddg)

 101  continue

      return
      end


