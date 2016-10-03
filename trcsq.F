      Complex*16 Function TrcSQ(Nmat,N,A)
      Implicit Real*8(A-H,O-Z)
C
C     Trace(A) with complex square matrix.
C
      Complex*16 GCmplx
      Dimension A(2,N,N,NMat)
      Save Zero
      Data Zero/0.0d0/
C
      SR = Zero
      SI = Zero
      Do 10 IM = 1, NMat
        Do 10 I = 1, N
            SR = SR + A(1,I,I,IM)
   10       SI = SI + A(2,I,I,IM)
      TrcSQ = GCmplx(SR,SI)
      Return
      End
