

      subroutine formden (iout, iprint, v, lenv, iphf, nbasis, 
     $     norbs, na, nb, ngrdt, indarr, grda, grdb, grdg,
     $     irwca, irwcb, irws, irwt, irwo)

      implicit none


C +------------------------------------------------------------+
C |  formden   --   RSS,  05.2013                              |
C |                                                            |
C |                                                            |
C |  Build the transition density  matrices, P(g) and the      |
C |  traces of the overlap N(g), needed for integration        |
C |  over the grid.                                            |
C +------------------------------------------------------------+

C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       nbasis - size of the basis
C       norbs  - number of orbitals
C       na     - number of alpha electrons
C       nb     - number of beta electrons
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grd?   - integration grid for each variable
C       irwca  - RWF file with alpha MO matrix in AO basis (Ca)
C       irwcb  - RWF file with beta MO matrix in AO basis  (Cb)
C       irws   - RWF file with the overlap matrix S
C       irwt   - RWF file with the transition density matrices P(g) [ out ]
C       irwo  - RWF file with the determinants x(g) of the overlap matrix N(g) [ out ]

      real*8      grda(*), grdb(*), grdg(*)
      integer     iout, iprint, lenv, iphf, nbasis, norbs, na, nb
      integer     ngrdt, indarr(ngrdt,*)
      integer     irwca, irwcb, irws, irwt, irwo

      real*8      v(*)

C     other variables

      integer     itype, nel, ncik, ntt, nbsq, szmo, szn, szp
      integer     js, jmo, jmot, jend
      integer     szbfx, szbfp, szbf, szngscr, szpgscr, szscr 

      integer     lenp, npmem, npmax, npuse, ipar
      integer     tid, ip,  iblock, blklen, n1c, n2c
      integer     lenio, mblklen, szpar
      integer     jpbuf, jxbuf 
      integer     ing, idtng, iscr, ipg, igrdt, iend, indp, inddt
      integer     ixgout, ipgout, indx 
      complex*16  xg

      logical     icmplx, ispuhf, ispghf, ispin

C     functions

      integer     omp_get_thread_num
      integer     intowp, nproc, lappar, igetrw

C     constants

      complex*16  zero, one

 1001 format (' formden: will use ', I2, ' processes out of ', I2)

      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)


C     %%%%%%%%%%%%%%%%%
C     %   Initialize  %
C     %%%%%%%%%%%%%%%%%

C     Construct useful variables.

      nel  = na + nb
      ntt  = nbasis * (nbasis + 1) / 2
      nbsq = nbasis * nbasis

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
        call gauerr ('Unrecognized iphf in formden.')
      endif

C       .. double the number of grid points if complex conjugation
C          restoration was performed

      ncik = 1
      if ( icmplx ) ncik = 2


C     Define useful array dimensions:
C       - size of N(g) matrices (szn)
C       - size of MO coefficients (szmo)

      if ( itype .eq. 1 ) then
        szn   = 2*na*na
        szmo  = 2*na*nbasis
        szp   = 2*nbsq

      elseif ( itype .eq. 2 ) then
        szn   = 2*(na*na + nb*nb)
        szmo  = 2*nel*nbasis
        szp   = 2*2*nbsq

      elseif ( itype .eq. 3 ) then
        szn   = 2*nel*nel
        szmo  = 2*nel*2*nbasis
        szp   = 2*2*2*nbsq

      endif


C     %%%%%%%%%%%%%%%%%%%%%%%%
C     %  Memory allocatioin  %
C     %%%%%%%%%%%%%%%%%%%%%%%%

C      Memory allocation
C     ===================

C       Allocate space for:
C         - rectangular block of MO coefficients C and C! (jmo,jmot)
C         - overlap matrix (js)

      js   = 1 
      jmo  = js   + 2*nbsq
      jmot = jmo  + szmo
      jend = jmot + szmo

      call tstcor (jend, lenv, 'formden')


C     Prepare memory for parallel execution
C     =====================================

C       - scratch space to call evalng for each processor (szngscr)
C         the memory requirements are described in evalng
C       - scratch space to call evalpg for each processor (szpgscr)
C       - szbfx - size of the buffer for x(g) at each grid point
C       - szbfp - size of the buffer for P(g) at each grid point
C       - lenp,  total length of scratch space per processor
C       - lenio, length of the total output per one grid point

      szbfx = ncik * ncik * 2
      szbfp = ncik * ncik * szp

      if ( itype .eq. 1 ) then
         szngscr = szmo + 2*szn + na
         szpgscr = szmo + szp

      elseif ( itype .eq. 2 ) then
         szngscr = szmo + 2*szn + nel
         szpgscr = szmo + szp

      elseif ( itype .eq. 3 ) then
         szngscr = szmo + 2*szn + nel + 40*nbsq
         szpgscr = szmo + szp   + 32*nbsq

      endif
      
      szscr  = max (szngscr, szpgscr)
      lenp   = szscr + szn
      lenio  = szp*ncik*ncik  + 2*ncik*ncik

C     ...define how many blocks we can keep in memory

      npmem = (lenv - jend + 1) / (lenp + lenio)

      if ( npmem .lt. 1 ) then
        call gauerr ('Not enough memory in formden for a single CPU.')
      endif

C     ...define how many processes we are allowed
C     ...check if LAPACK is allowed to run in parallel

      ipar  = lappar (0)
      if ( ipar .eq. 1 ) then
         npmax = nproc (0)
      else
         npmax = 1
      endif
      
      if ( npmem .gt. npmax ) then
         npuse  = npmax
         mblklen = min (ngrdt,(lenv - jend + 1 - lenp*npmem) / (lenio))
         szpar  = lenp*npmax
      else
         npuse  = npmem
         mblklen = npmem
         szpar  = mblklen*lenp
      endif

C     Allocate memory for the I/O scratch per block
      
      jpbuf  = jend
      jxbuf  = jpbuf  + mblklen*szbfp
      jend   = jxbuf  + mblklen*szbfx

      call tstcor (jend + szpar  - 1, lenv, 'formden')

C     %%%%%%%%%%%%%%%%%%%%%%%%
C     %  Prepare calculatoin %  
C     %%%%%%%%%%%%%%%%%%%%%%%%

C     Load overlap matrix

      call fileio (2, -irws, 2*ntt, v(js), 0)
      call zsquare (v(js), v(js), nbasis, nbasis, 0)

C     Load MO coefficients matrix C and construct C!

      if ( itype .eq. 1 ) then
         
         call fileio (2, -irwca, szmo, v(jmo), 0)

C        ...form a conjugate transpose

         call trprct (2,nbasis,na,v(jmo),v(jmot))
         call vconjg (nbasis*na,v(jmot)) 
         
      elseif ( itype .eq. 2 ) then

         call fileio (2, -irwca, 2*na*nbasis, v(jmo), 0)
         call fileio (2, -irwcb, 2*nb*nbasis, v(jmo+2*na*nbasis), 0)

C        ...arrange spin blocks properly in C!: (alpha)
C                                               (beta )
         call trprct (2,nbasis,na,v(jmo),v(jmot))
         call trprct (2,nbasis,nb,v(jmo+2*na*nbasis),
     $        v(jmot+2*na*nbasis))

         call vconjg (nbasis*nel,v(jmot)) 
         
      elseif ( itype .eq. 3 .and. .not. ispuhf ) then

C        ...na first alpha orbitals
         call fileio (2, -irwca, 2*nel*2*nbasis, v(jmo), 0)

C        ...form a conjugate transpose

         call trprct (2,2*nbasis,nel,v(jmo),v(jmot))
         call vconjg (2*nbasis*nel,v(jmot)) 

      elseif ( ispuhf ) then 

C        ...special care is needed if PSUHF is requested
C        ...we need to extend the orbitals to length 2*nbasis

C        ...aa block
         call fileio (2, -irwca, 2*na*nbasis, v(jmo), 0)

C        ...add a stride of zeros
         call aclear (2*nel*nbasis, v(jmo+2*na*nbasis))

C        ...bb block
         call fileio (2, -irwcb, 2*nb*nbasis, v(jmo +
     $        2*na*nbasis + 2*nel*nbasis), 0)
         
C        ...arrange blocks properly
         call trprct (2,nbasis,nel, v(jmo), v(jmot))
         call trprct (2,nbasis,nel, v(jmo + 2*nel*nbasis),
     $        v(jmot + 2*nel*nbasis))

C        ...form C
         call trprct (2,nel, 2*nbasis,v(jmot), v(jmo))
C        ...form C!
         call vconjg (2*nbasis*nel, v(jmot))
      endif

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Construct the matrices N(g) and then P(g) %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      do 103 iblock = 1, ngrdt, mblklen

        blklen = min (mblklen, ngrdt - iblock + 1)

C     Computation loop. 

C$omp   parallel default(shared)
C$omp+  private (tid, ip, n1c, n2c, igrdt)
C$omp+  private (ing, iscr, iend, ixgout, ipgout)
C$omp+  private (xg)
        tid = omp_get_thread_num()
C$omp   do schedule(static,1) collapse(3)
        do 104 ip = 1, blklen
           do 105 n1c = 1, ncik
              do 106 n2c = 1, ncik

                 igrdt = iblock + ip - 1

C     Memory allocation per processor:
C              ing     - current N(g) matrix
C              ingscr  - scratch space to run evalng or evalpg

                 ing     = jend  + tid * lenp
                 iscr    = ing   + szn  
                 iend    = iscr  + szscr

C     Compute pointers for P(g) and x(g) in the output array

                 ixgout = jxbuf + (ip - 1)*ncik*ncik*2 + 
     $                (n1c - 1)*ncik*2 + (n2c - 1)*2 
                 ipgout = jpbuf + (ip - 1)*ncik*ncik*szp + 
     $                (n1c - 1)*ncik*szp + (n2c - 1)*szp

C     Calculate N(g), x(g) and P(g) and save results in the output array

                 call evalng (iout, iprint, v(iscr), szscr, itype, 
     $                ispin, nbasis, na, nb, ngrdt, igrdt, indarr, grda,   
     $                grdb, grdg, v(jmo), v(jmot), v(js), n1c, n2c, 
     $                v(ing), xg)

                 call evalpg (iout, iprint, v(iscr), szscr, itype, 
     $                ispin, nbasis, na, nb, ngrdt, igrdt, indarr, grda,  
     $                grdb, grdg, v(jmo), v(jmot), v(ing), n1c, n2c, 
     $                v(ipgout))

                 v(ixgout)     = dble (xg)
                 v(ixgout + 1) = aimag (xg)

 106          continue
 105       continue
 104    continue
C$omp   end do
C$omp  end parallel
C     Save matrices
       
       indp  = (iblock - 1)*szbfp
       indx  = (iblock - 1)*szbfx
       
C     Save matrices P(g) into RWF file

       call fileio (1, -irwt, blklen*szbfp, v(jpbuf), indp)
       call fileio (1, -irwo, blklen*szbfx, v(jxbuf), indx)
       
 103  continue

      return
      end


