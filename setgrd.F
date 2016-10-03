

      subroutine setgrd (iout, iprint, v, lenv, ngrdt,
     $     ngrda, ngrdb, ngrdg, indarr, grda, wgta,
     $     grdb, wgtb, grdg, wgtg)

      implicit none

C     input / output variables

      real*8      v(*)
      integer     iout, iprint, lenv
      integer     nbasis, norbs, nbas6d
      integer     ngrda, ngrdb, ngrdg, ngrdt
      integer     indarr(*)
      real*8      grda(*), grdb(*), grdg(*)
      real*8      wgta(*), wgtb(*), wgtg(*)

C     memory allocation variables

      integer     jgrda, jgrdb, jgrdg, jwgta, jwgtb, jwgtg

C     other variables

      real*8      pi
      integer     i

      call indgrd (iout, iprint, indarr, ngrdt, ngrda, ngrdb, ngrdg,
     $     1)

      pi = 4.0d0 * atan (1.0d0)

      if ( ngrda .gt. 1 ) then
        call trpgrd (0.0d0, 2*pi, grda, wgta, ngrda)
      else
        grda(1) = 0.0d0
        wgta(1) = 1.0d0
      endif

      if ( ngrdb .gt. 1 ) then
        call gauleg (0.0d0, pi, grdb, wgtb, ngrdb)

        do 10 i = 1, ngrdb
          wgtb(i) = wgtb(i) * sin (grdb(i))
 10     continue
      else
        grdb(1) = 0.0d0
        wgtb(1) = 1.0d0
      endif

      if ( ngrdg .gt. 1 ) then
        call trpgrd (0.0d0, 2*pi, grdg, wgtg, ngrdg)
      else
        grdg(1) = 0.0d0
        wgtg(1) = 1.0d0
      endif

      end
