

      subroutine fixgaufmt (iout, iprint, v, lenv, iopcl, nbasis, norbs,
     $     irwpa, irwca, irwpa1, irwca1)

      implicit none


C +------------------------------------------------------------+
C |  fixgaufmt   --   RSS, 05.2013                             |
C |                                                            |
C |  This routine rearranges GHF-type density matrix and MO    |
C |  coefficients to keep them as spin blocks.                 |
C |  That is, a general matrix S                               |
C |  will be stored as                                         |
C |                                                            |
C |    S  =  ( S_aa  S_ab )                                    |
C |          ( S_ba  S_bb )                                    |
C |                                                            |
C |  On the other hand, Gaussian stores GHF matrices in a      |
C |  different order. The basis functions are ordered          |
C |  according to                                              |
C |                                                            |
C |    PHF:       1a, 2a, ..., 1b, 2b, ...                     |
C |    Gaussian:  1a, 1b, 2a, 2b, ...                          |
C |                                                            |
C |  Thus, for GHF matrices, we need to reorganize rows and    |
C |  columns.                                                  |
C |                                                            |
C |  The variable iopcl is used to determine how Gaussian is   |
C |  storing these matrices.                                   |
C |                                                            |
C |     iopcl = 0,  real RHF                                   |
C |           = 1,  real UHF                                   |
C |           = 2,  complex RHF                                |
C |           = 3,  complex UHF                                |
C |           = 6,  complex GHF (there is no real GHF)         |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       v      - scratch array
C       lenv   - size of scratch array
C       iopcl  - describes type of HF wavefunction
C       nbasis - number of basis functions
C       norbs  - number of orbitals
C       irwpa  - RWF file with 1 PDM
C       irwca  - RWF file with MO coefficients
C       irwpa1  - working version of irwpa [ out ]
C       irwca1  - working version of irwca [ out ]

      real*8      v(*)
      integer     iout, iprint, lenv, iopcl, nbasis, norbs
      integer     irwpa, irwca, irwpa1, irwca1

C     other variables

      integer     nbsq, ntt2 
      integer     i, j, ind1, ind2
      integer     jorb, jp, jz, jscr
      integer     jend, mdv

C     functions

      integer     igetrw



C     Return if iopcl <= 3. (RHF of UHF)

      if ( iopcl .le. 3 ) then
        irwpa1 = irwpa
        irwca1 = irwca

        return
      endif


C     Construct useful variables.

      ntt2 = nbasis * (2*nbasis + 1)
      nbsq = nbasis * nbasis
      

C     Create scratch files for matrices of the correct dimension.

      irwpa1 = igetrw (2*ntt2, 0)
      irwca1 = igetrw (8*nbasis*norbs, 0)

C     Memory allocation.

C       Allocate space for:
C         - orbital coefficients (jorb)
C         - density matrix, square format (jp)
C         - scratch array for matrix reorganization (jscr)
C         - scratch arrays for lt to square transformation (jz)

      jorb = jend
      jp   = jorb + 8*nbsq
      jscr = jp   + 2*ntt2
      jz  = jscr  + 8*nbsq
      jend = jz   + 8*nbsq
   
      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'fixgaufmt')

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %      Orbital coefficients      %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C     Clear arrays of orbital coefficients.

      call aclear (8*nbsq, v(jorb))

C     Read gaussian orbitals from file.

      call fileio (2, -irwca, 8*nbasis*norbs, v(jorb), 0)
            
C     debugging...
C       print molecular orbitals

C$$$      write (iout, *) ' '
C$$$      call ltoutx (iout, 3, 2, 1, 0, 'MO coefficients', 2*nbasis,
C$$$     $     2*norbs, v(jorb), v(jorb), 0)


C     Reorganize the matrix for GHF-type MOs (see fixmat for more
C     details).

      call amove (8*nbsq, v(jorb), v(jscr))
      
C     .. now reorganize rows
      
      call trspnc (.false., 2*nbasis, 2*nbasis, v(jscr))
      
      do 11 i = 1, 2*nbasis
         if ( mod (i,2) .eq. 1 ) then
            j = (i-1)/2 + 1
         elseif ( mod (i,2) .eq. 0 ) then
            j = i/2 + nbasis
         endif
         
         ind1 = (i - 1) * 4*nbasis + jscr
         ind2 = (j - 1) * 4*nbasis + jorb
         
         call amove (4*nbasis, v(ind1), v(ind2))
 11   continue

      call trspnc (.false., 2*nbasis, 2*nbasis, v(jorb))


C     debugging...
C       print molecular orbitals

C$$$      write (iout, *) ' '
C$$$      call ltoutx (iout, 3, 2, 1, 0, 'MO coefficients', 2*nbasis,
C$$$     $     2*norbs, v(jorb), v(jorb), 0)

C     Save orbitals to RWF file.

      call fileio (1, -irwca1, 8*nbasis*norbs, v(jorb), 0)
      
C     %%%%%%%%%%%%%%%%%%%%%%%
C     %    Density matrix   %
C     %%%%%%%%%%%%%%%%%%%%%%%

C     Read density matrix from file.

      call fileio (2, -irwpa, 2*ntt2, v(jp), 0)

C     Transform density matrix to square format.

      call zsquare (v(jp), v(jz), 2*nbasis, 2*nbasis, 3)

C     debugging...
C       print density matrix

C$$$      write (iout, *) ' '
C$$$      call ltoutx (iout, 1, 2, 1, 0, 'Density matrix', 2*nbasis,
C$$$     $     2*nbasis, v(jz), v(jz), 0)

C     Clear the scratch array for transformation

      call aclear (8*nbsq, v(jscr))

C     Reorder rows and columns in GHF-type density matrix.
C       .. reorganize columns

      do 20 i = 1, 2*nbasis
         if ( mod (i,2) .eq. 1 ) then
            j = (i-1)/2 + 1
         elseif ( mod (i,2) .eq. 0 ) then
            j = i/2 + nbasis
         endif
         
         ind1 = (i - 1) * 4*nbasis + jz
         ind2 = (j - 1) * 4*nbasis + jscr
         
         call amove (4*nbasis, v(ind1), v(ind2))
 20   continue
      
C       .. now reorganize rows

      call trspnc (.false., 2*nbasis, 2*nbasis, v(jscr))
      
      do 21 i = 1, 2*nbasis
         if ( mod (i,2) .eq. 1 ) then
            j = (i-1)/2 + 1
         elseif ( mod (i,2) .eq. 0 ) then
            j = i/2 + nbasis
         endif
         
         ind1 = (i - 1) * 4*nbasis + jscr
         ind2 = (j - 1) * 4*nbasis + jz
         
         call amove (4*nbasis, v(ind1), v(ind2))
 21   continue
      
      call trspnc (.false., 2*nbasis, 2*nbasis, v(jz))

C     debugging...
C       print density matrix

C$$$      write (iout, *) ' '
C$$$      call ltoutx (iout, 1, 2, 1, 0, 'Density matrix', 2*nbasis,
C$$$     $     2*nbasis, v(jz), v(jz), 0)

C     Place the reordered density matrix to lower triangular
      
      call zlinear(v(jz),v(jp),2*nbasis,2*nbasis)

C     Save the reordered density matrix to temporary file

      call fileio (1, -irwpa1, 2*ntt2, v(jp), 0)

      return
      end


