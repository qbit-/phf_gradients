


      subroutine evallg (iout, iprint, v, lenv, itype, ncik, nbasis,
     $     eref, ngrdt, irwh, irwg, irwt, irwl)

      implicit none


C +------------------------------------------------------+
C |  evallg   --   RSS,  06.2013                         |
C |                                                      |
C |                                                      |
C |  Compute the L-term needed for the computation of    |
C |  the one-electron part of the gradient.              |
C |                                                      |
C +------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       itype  - decides wich spn blocks need to be kept
C       ncik   - number of basis states for complex conj projection
C       nbasis - size of the basis
C       eref   - reference electronic energy
C       ngrdt  - total number of grid points
C       irwh   - RWF file with the one particle Hamiltonian
C       irwg   - RWF file with the two particle Hamiltonian G(g)
C       irwt   - RWF file with the transition density matrices P(g)
C       irwl   - RWF file containing the L-term [ out ]

      integer     iout, iprint, lenv, itype, ncik, nbasis, ngrdt
      integer     irwh, irwg, irwt, irwo, irwl

      real*8      v(*), eref

C     other variables

      integer    ntt, nbsq, szmat, szmatc
      integer    jh, jhf, jscr, jend, mdv 
      integer    szbf
      integer    lenp, lenio, npmax, npmem, npuse, mblklen, ipar
      integer    jpgin, jggin, jlgout
      integer    jl, jtmp
      integer    iblock, blklen, tid, ip, n1c, n2c, igrdt
      integer    il, itmp, iend, ilgout, ipgin, iggin, indl

      complex*16 eg

C     functions

      integer     omp_get_thread_num
      integer     nproc, lappar
      complex*16  trcabc, trcsq

C     constants

      complex*16  zero, one, mone, two

C     debug 

      complex*16  ts
      integer     sztst
      integer     jtst, itst
      integer     irwtst
      integer     indtst
      irwtst  = 824

      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      mone = cmplx (-1.0d0, 0.0d0)
      two  = cmplx (2.0d0, 0.0d0)

C     Construct useful variables.

      ntt  = nbasis * (nbasis + 1) / 2
      nbsq = nbasis * nbasis

      if ( itype .eq. 1 ) then
         szmatc = nbsq
         szmat  = 2 * szmatc
      elseif ( itype .eq. 2 ) then
         szmatc = 2 * nbsq
         szmat  = 2 * szmatc
      elseif ( itype .eq. 3 ) then
         szmatc = 4 * nbsq
         szmat  = 2 * szmatc
      endif

C     %%%%%%%%%%%%%%%%%%%%%%%%
C     %  Memory allocatioin  %
C     %%%%%%%%%%%%%%%%%%%%%%%%

C      Memory allocation
C     ===================

C       Allocate space for:
C         - alpha part of 1-particle Hamiltonian (jh)
C         - full hamiltonian
C         - scratch space

      jh   = 1 
      jhf  = jh  + 2*ntt
      jend = jhf + szmat

      if ( itype .eq. 3 ) then 
         jscr = jend 
         jend = jscr + 2*nbsq
      endif

C     Prepare memory for parallel execution
C     =====================================

C       - szbf   - size of the IO buffer at each grid point 
C                        for each of the matrices (P,G or L)
C
C       - lenp,  total length of scratch space per processor
C       - lenio, length of the total output per one grid point

C     IO:
      szbf   = ncik * ncik * szmat

C     debug
      sztst  = 2 * ncik * ncik

C     Scratch buffers:

      lenp   = szmat + szmat  
      lenio  = szbf + szbf + szbf + sztst

C     ...define how many blocks we can keep in memory

      npmem = (lenv - jend + 1) / (lenp + lenio)

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in evallg for a single CPU.')
      endif

C     ...define how many processes we can use
C     ...check if LAPACK is allowed to run in parallel

      ipar  = lappar (0)
      if ( ipar .eq. 1 ) then
         npmax = nproc (-1)
      else
         npmax = 1
      endif

      if ( npmem .gt. npmax ) then
         npuse   = npmax
         mblklen = min (ngrdt,(lenv - jend + 1 - lenp*npmem) / (lenio))
      else
         npuse   = npmem
         mblklen = npmem
      endif

C     Allocate memory for the I/O scratch per block
C     IO:

      jpgin  = jend
      jggin  = jpgin  + mblklen*szbf
      jlgout = jggin  + mblklen*szbf
      jend   = jlgout + mblklen*szbf

C     debug
      jtst   = jend
      jend   = jtst  + mblklen*sztst

C     Allocate scratch for threads
C     SCRATCH:
      
      jl   = jend
      jtmp = jl   + npuse*szmat
      jend = jtmp + npuse*szmat

      call tstcor (jend - 1, lenv, 'evallg')

      mdv = lenv - jend + 1

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %   Prepare calculation    %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     Load core Hamiltonian


      if ( itype .eq. 1 ) then
         call fileio (2, -irwh, 2*ntt, v(jh), 0)
         call zsquare (v(jh), v(jhf), nbasis, nbasis, 0)

      elseif ( itype .eq. 2 ) then
         call fileio (2, -irwh, 2*ntt, v(jh), 0)
         call zsquare (v(jh), v(jhf), nbasis, nbasis, 0)
         call fileio (2, -irwh, 2*ntt, v(jh), 0)
         call zsquare (v(jh), v(jhf+2*nbsq), nbasis, nbasis, 0)

      elseif ( itype .eq. 3 ) then
         call fileio (2, -irwh, 2*ntt, v(jh), 0)
         call aclear (2*nbsq, v(jscr))
         call dmblck (iout, iprint, v(jend), mdv, nbasis, v(jhf), 
     $        v(jh), v(jscr), v(jscr), v(jh), 1, -2)
      endif

C     Request correct size of RWF to store L(g)

      call conddf (irwl, ngrdt*szbf)

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Construct the terms L(g)                  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Loop over grid points.

      do 103 iblock = 1, ngrdt, mblklen

        blklen = min (mblklen, ngrdt - iblock + 1)

C     Load P(g) and G(g) in blocks
        call fileio (2, -irwt, blklen*szbf, v(jpgin),
     $       (iblock-1)*szbf)
        call fileio (2, -irwg, blklen*szbf, v(jggin),
     $       (iblock-1)*szbf)
 
C     Computation loop. 

C$omp   parallel default(shared)
C$omp+  private (ip, tid, n1c, n2c, igrdt)
C$omp+  private (il, itmp, ilgout, ipgin) 
C$omp+  private (eg,ts,itst)
        tid = omp_get_thread_num()

C$omp   do schedule(static,1) collapse (3)
        do 104 ip = 1, blklen
           do 105 n1c = 1, ncik
              do 106 n2c = 1, ncik

                 igrdt = iblock + ip - 1

C     Memory allocation per processor:
                 il     = jl   + tid*szmat
                 itmp   = jtmp + tid*szmat 

C     debug 
                 itst   = jtst + (ip - 1)*ncik*ncik*2 + 
     $                (n1c - 1)*ncik*2 + (n2c - 1)*2

C     Compute pointers for P(g) and L(g) in the IO arrays

                 ilgout = jlgout + (ip - 1)*ncik*ncik*szmat + 
     $                (n1c - 1)*ncik*szmat + (n2c - 1)*szmat 
                 ipgin = jpgin + (ip - 1)*ncik*ncik*szmat + 
     $                (n1c - 1)*ncik*szmat + (n2c - 1)*szmat
                 iggin = jggin + (ip - 1)*ncik*ncik*szmat + 
     $                (n1c - 1)*ncik*szmat + (n2c - 1)*szmat

C     Compute energy at point g: tr((h + 1/2 G(g)) . P(g))
C     Warning! trcabc computes tr(a . b*), so we conjugate G first

                 call amove (szmat, v(jhf), v(itmp))

                 call zaxpy (szmatc, one/two, v(iggin), 1, v(itmp), 1)
                 call vconjg (szmat, v(itmp))

                 if ( itype .eq. 1 ) then
                    eg = trcabc(1, nbasis, v(ipgin), v(itmp))
                    eg = eg * two - cmplx(eref)
                 elseif ( itype .eq. 2 ) then
                    eg = trcabc(2, nbasis, v(ipgin), v(itmp))
                    eg = eg - cmplx(eref)
                 elseif ( itype .eq. 3 ) then
                    eg = trcabc(1, 2*nbasis, v(ipgin), v(itmp))
                    eg = eg - cmplx(eref)
                 endif

                 call vconjg (szmat, v(itmp))

C     debug 
                 call vconjg (szmat, v(iggin))

                 if ( itype .eq. 1 ) then
                    ts = trcabc(1, nbasis, v(ipgin), v(iggin))
                 elseif ( itype .eq. 2 ) then
                    ts = trcabc(2, nbasis, v(ipgin), v(iggin))
                    ts = ts*cmplx(0.5d0,0.0d0)
                 elseif ( itype .eq. 3 ) then
                    ts = trcabc(1, 2*nbasis, v(ipgin), v(iggin))
                    ts = ts*cmplx(0.5d0,0.0d0)
                 endif

                 call vconjg (szmat, v(iggin))
                 v(itst)   = dble(ts)
                 v(itst+1) = aimag(ts)

C     Compute L(g) = eg . P(g) - P(g) . (h + G(g)) . P(g)

                 call zaxpy (szmatc, one/two, v(iggin), 1, v(itmp), 1)

                 if ( itype .eq. 1 ) then
                    call zgemm ('n', 'n', nbasis, nbasis, nbasis, one,
     $                   v(itmp), nbasis, v(ipgin), nbasis, zero,
     $                   v(il), nbasis)
                    call zgemm ('n', 'n', nbasis, nbasis, nbasis, one,
     $                   v(ipgin), nbasis, v(il), nbasis, zero,
     $                   v(itmp), nbasis)

                 elseif ( itype .eq. 2 ) then
                    call zgemm ('n', 'n', nbasis, nbasis, nbasis, one,
     $                   v(itmp), nbasis, v(ipgin), nbasis, zero,
     $                   v(il), nbasis)
                    call zgemm ('n', 'n', nbasis, nbasis, nbasis, one,
     $                   v(itmp + 2*nbsq), nbasis, v(ipgin + 2*nbsq),
     $                   nbasis, zero, v(il + 2*nbsq), nbasis)

                    call zgemm ('n', 'n', nbasis, nbasis, nbasis, one,
     $                   v(ipgin), nbasis, v(il), nbasis, zero,
     $                   v(itmp), nbasis)
                    call zgemm ('n', 'n', nbasis, nbasis, nbasis, one,
     $                   v(ipgin + 2*nbsq), nbasis, v(il + 2*nbsq),
     $                   nbasis, zero, v(itmp + 2*nbsq), nbasis)

                 elseif ( itype .eq. 3 ) then
                    call zgemm ('n', 'n', 2*nbasis, 2*nbasis, 2*nbasis,
     $                   one, v(itmp), 2*nbasis, v(ipgin),
     $                   2*nbasis, zero, v(il), 2*nbasis)
                    call zgemm ('n', 'n', 2*nbasis, 2*nbasis, 2*nbasis,
     $                   one, v(ipgin), 2*nbasis, v(il),
     $                   2*nbasis, zero, v(itmp), 2*nbasis)
                 endif

                 call aclear (szmat, v(ilgout))
                 call zaxpy (szmatc, eg, v(ipgin), 1, v(ilgout), 1)
                 call zaxpy (szmatc, mone, v(itmp), 1, v(ilgout), 1)

 106          continue
 105       continue
 104    continue
C$omp  end do
C$omp  end parallel

C     Save matrices
       
       indl  = (iblock - 1)*szbf
       
C     Save matrix L(g) into RWF file

       call fileio (1, -irwl, blklen*szbf, v(jlgout), indl)

       indtst = (iblock - 1)*sztst
       call fileio (1, -irwtst, blklen*sztst, v(jtst), indtst)

 103  continue

C     debug 
      npmax = nproc (0)


      return 
      end
