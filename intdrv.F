


      subroutine intdrv (iout, iprint, v, lenv, iphf, ncis, ncik,  
     $     ngrda, ngrdb, ngrdg, ngrdt, indarr, grda, wgta, grdb, wgtb,  
     $     grdg, wgtg, dmt, dind, civec, irwo, irwq, lenq, qint)

      implicit none

C +------------------------------------------------------------+
C |  intdrv   --   RSS, 06.2013                                |
C |                                                            |
C |                                                            |
C |                                                            |
C |  Integrates some quntity Q over the projection             |
C |  grid and performs a summation over CI coefficients:       |
C |                                                            |
C |        Q =                                                 |
C |  sum_{i,j,k,l} C*_{i,k} C_{j,l} int dg w(g;i,j,k,l) .      |
C |                                    . q(g;i,j) =            | 
C |          = sum_{i,j} int dg y(g,i,k) . q(g,i,k)            |
C |                                                            |
C |  The quantity q(g;i,j) is assumed to depend on CI          |
C |  indicies (i,j) corresponding to the complex conjugation   |
C |  projection and on the grid coordinate g.                  |
C |  The weight factor w(g;i,j,k,l) depends on the projection  |
C |  we do and on the overlap x(g,i,j) of state at point g and |
C |  a reference state.                                        |
C |  If it's only complex conjugation projection, then         |
C |  w = 1; if we also do spiin projection, then w is a        |
C |  Wigner D matrix.                                          |
C |  The CI coefficients are the ones we get from the small CI |
C |  problem associated with projection. (See C. Himenez-Hoyos |
C |  and G. Scuseria for details).  Also, as there is no       |
C |  analithycal expression for integration weights in case of |
C |  complex conjugation, this information is also obtained by |
C |  solving  CI problem                                       |
C |  (so in case of complex conjugation we  solve 4 times      |
C |  larger CI problem). The indicies i,j correspond to        |
C |  complex conjugation, and k,l to spin projection.          |
C |  To clarify what we say let's consider an example.         |
C |                                                            |
C |  Assume we did KSGHF projection. Then, we will have        |                     
C |  four different values of Q  at each grid point   :        |
C |       f(g;i,j) =                                           |
C |         <phi| I! F(g) I |phi>    i = 1, j = 1              |
C |         <phi| I! F(g) K |phi>    i = 1, j = 2              |
C |         <phi| K! F(g) I |phi>    i = 2, j = 1              |
C |         <phi| K! F(g) K |phi>    i = 1, j = 2              |
C |                                                            |
C |  As for the weight factors, if the S^2 = S, the total      |
C |  length of CI coefficients is ncis*ncik, ncis =  2S + 1.   |
C |   Then                                                     |
C |    w(g; i,j,k,l) =                                         |
C |      [D^S_{k,l} (Omega)]* . x(i,j) / N                     |
C |                                                            |
C |    k,l = -S...S                                            |
C |                                                            |
C |    where x(g,i,j) is the determinant of the overlap        |
C |    (see formden) and N is a normalization                  |
C |                                                            |
C |  N = sum_{i,j,k,l} C*(i,k) . C(l,j) . int dg w(g; i,j,k,l) |
C |                                                            |
C |  The final weight factor z can be thus computed as         |
C |  y(g,i,j) = sum_{k,l} . C*_{i,k} . w(g,i,j,k,l) . C_{l,j}  |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iphf   - type of PHF calculation (see phfdrv for details)
C       ncis   - number of CI vectors in spin projection
C       ncik   - number of CI vectors in complex conjugation projection
C       ngrda  - number of grid points in alpha integration
C       ngrdb  - number of grid points in beta integration
C       ngrdg  - number of grid points in gamma integration
C       ngrdt  - total number of grid points
C       indarr - mapping array igrdt => igrdX
C       grdX   - integration grid for each variable
C       ygtX   - array with quadrature weights for each variable
C       dmt    - array with small Wigner d matrix for each beta, see wigarr
C       dind   - indicies for Wigner d matrix
C       civec  - vector of CI coefficients associated with proojectors 
C       irwo   - RWF file with overlaps x(g)
C       irwq   - RWF file with values of Q at each grid point
C       lenq   - length of q (in real numbers)
C       qint   - integrated quantity[ out ]

      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)
      complex*16  civec(*)
      integer     ncis, ncik
      integer     iout, iprint, lenv, iphf, lenq
      integer     ngrda, ngrdb, ngrdg, ngrdt
      integer     indarr(ngrdt,*)
      integer     irwci, irwq, irwo

      real*8      v(*), dmt(ngrdb,ncis,ncis), dind(ncis), qint(*)

C     other variables
      
      real*8      facw 
      integer     jend, jqin, jxin, jqout, jqbuf, jyout
      integer     szqint, szqintc
      integer     lenp, lenio, npmem, npmax, npuse
      integer     szqloc, szxloc, szbfx, szbfq
      integer     mblklen, blklen, szpar, ipar
      integer     iblock, ip, iqin, ixin, tid, iyout, iqout, iqbuf
      integer     igrda, igrdb, igrdg, igrdt, n1c, n2c, n1s, n2s, n1, n2 

      logical     icmplx, ispuhf, ispghf, ispin
      complex*16  yg, facd, norm

C     functions

      integer     intowp, nproc, lappar
      integer     omp_get_thread_num

C     constants

      complex*16  zero, one, im


      zero = cmplx (0.0d0, 0.0d0)
      one  = cmplx (1.0d0, 0.0d0)
      im   = cmplx (0.0d0, 1.0d0)


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
        call gauerr ('Unrecognized iphf in intfrc.')
      endif

C     Define dimensions of quantities used:
C      - szqint - size of the Q at each grid point

      szqintc = lenq/2
      szqint  = lenq

C     Initialize pointer

      jend  = 1

C     Prepare memory for parallel execution
C     =====================================

C       - lenp,  total length of scratch space per processor
C       - lenio, length of the total output per one grid point

C     IO:
      szbfq  = ncik * ncik * szqint
      szbfx  = ncik * ncik * 2

C     Scratch buffers:
      szqloc  = 2*szqint   ! We need two buffers to store 
      szxloc  = 2          ! current value and sum

      lenp   = szqloc + szxloc 
      lenio  = szbfq + szbfx

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
         szpar   = lenp*npmax
      else
         npuse   = npmem
         mblklen = npmem
         szpar   = mblklen*lenp
      endif
      
C     Allocate memory for the I/O buffers
      
      jqin   = jend
      jxin   = jqin  + mblklen*szbfq
      jqout  = jxin  + mblklen*szbfx
      jqbuf  = jqout + npuse*szqint
      jyout  = jqbuf + npuse*szqint
      jend   = jyout + npuse*2

      call tstcor (jend + szpar - 1, lenv, 'intfrc')

C     Initialize arrays used to accumulate sums

      call aclear (npuse*szqint, v(jqout))
      call aclear (npuse*2, v(jyout))

C     Loop over grid points.

      do 103 iblock = 1, ngrdt, mblklen
        blklen = min (mblklen, ngrdt - iblock + 1)

C     Load integrants and weight factors.
        
        call fileio (2, -irwq, blklen*szbfq, v(jqin),
     $       (iblock-1)*szbfq)
        call fileio (2, -irwo, blklen*szbfx, v(jxin), 
     $       (iblock-1)*szbfx)

C     Computation loop.

C$omp   parallel do schedule(static,1)
C$omp+  private (ip, n1c, n2c, tid, igrdt, iqin)
C$omp+  private (iqout, iqbuf, ixin, iyout)
C$omp+  private (igrda, igrdb, igrdg, facw, yg)
C$omp+  private (n1s, n2s, n1, n2)  default(shared)
C$omp+  collapse(3) 
        do 105 ip = 1, blklen
           do 106 n1c = 1, ncik
              do 107 n2c = 1, ncik

C                Compute the pointers to x(g) and QINT(g)
C                iqin   - pointer to the next matrix QINT
C                ixin   - pointer to the next x(g)
C                iyout  - partial sum of y(g) in each process
C                iqout  - partial sum of Q(g) in each process
C                iqbuf  - buffer for QINT
C                igrdt  - current grid point number

                 iqin  = jqin + (ip - 1)*ncik*ncik*szqint + 
     $                (n1c - 1)*ncik*szqint + (n2c - 1)*szqint
                 ixin  = jxin + (ip - 1)*ncik*ncik*2 + 
     $                (n1c - 1)*ncik*2 + (n2c - 1)*2

                 tid   = omp_get_thread_num()
                 iyout = jyout + tid*2  
                 iqout = jqout + tid*szqint
                 iqbuf = jqbuf + tid*szqint

C               Compute factor x(g)
                 
                 if ( .not. ispin ) then
C                   Only complex conjugate projection
                    yg = conjg (civec(n1c)) * civec(n2c)                    
                 
                 else
C                   Spin projection

C                   Retrieve grid coordinates.
                    igrdt = iblock + ip - 1                 
                    
                    igrda = indarr(igrdt,1)
                    igrdb = indarr(igrdt,2)
                    igrdg = indarr(igrdt,3)

C                   Compute integration weight factor
                    facw = wgta(igrda) * wgtb(igrdb) * wgtg(igrdg)

C     Compute the factor sum_{k,l} C*_{i,k} . D^S_{k,l} . C_{l,j}
C     where D^S_{k,l} is a conjugate of  Wigner's rotation matrix
                   
C     Note that
C     ( D^S_{MK} (Omega) )*  =  exp (+i * alpha * M) *
C                               exp (+i * gamma * K) *
C                                   d^S_{MK} (beta),

C     where d^S_{MK} (beta) is Wigner's small d-matrix.
C     
C     Note also, that for UHF the integration w.r.t alpha and 
C     gamma vanishes, and the only Wigner's  matrix we need
C     is D^S_{Sz, Sz}

                   
                    if ( ispghf ) then
C                      Compute the w_{i,j} element.
C                      Get the shifts for "spin" indicies i and j in civect
                       n1 = (n1c - 1)*ncis
                       n2 = (n2c - 1)*ncis

                       yg = zero
                       do 111 n2s = 1, ncis
                          facd = zero
                          do 121 n1s = 1, ncis   
                             facd = facd + conjg (civec(n1 + n1s)) *
     $                     exp (+im * cmplx (grda(igrda) * dind(n1s))) *
     $                     exp (+im * cmplx (grdg(igrdg) * dind(n2s))) *
     $                     cmplx (dmt(igrdb, n1s, n2s))
                            
 121                      continue
                          yg = yg + facd*civec(n2 + n2s)
                         
 111                   continue
                   
                    elseif ( ispuhf ) then
                       yg = cmplx (dmt(igrdb,1,1)) * 
     $                      conjg (civec(n1c)) * civec(n2c)

                    endif
C                   Multiply by the integration weight 
                    yg = yg * facw

                 endif
                
C                Multiply by the overlap x(g)
                 yg = yg * cmplx (v(ixin), v(ixin+1))

C                Accumulate y(g) to compute normalization
                 v(iyout)   = v(iyout) + dble (yg)
                 v(iyout + 1) = v(iyout + 1) + aimag (yg)
             
C                Accumulate QINT
                 call amove (szqint, v(iqin), v(iqbuf))
                 call zscal (szqintc, yg, v(iqbuf), 1)
                 call aadd (szqint, v(iqout), v(iqbuf), v(iqout))
                                  
 107          continue
 106       continue
 105    continue
C$omp   end parallel do

 103  continue

C     Accumulate results from all threads

      call aclear (szqint, qint)
      norm = zero

      do 201 ip = 1, npuse

         iqout = jqout + (ip-1)*szqint
         iyout = jyout + (ip-1)*2  
         yg    = cmplx (v(iyout), v(iyout + 1))
         
         norm  = norm  + yg
         call aadd (szqint, v(iqout), qint, qint)

 201  continue

C     Do final normalization

      call zscal (szqintc, one/norm, qint, 1)
      
      return
      end
