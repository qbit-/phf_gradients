ela
===

This code provides analytic gradients:

Schutski, Roman, et al. "Analytic energy gradient for the projected Hartree–Fock method"
The Journal of chemical physics 140.20 (2014): 204101

for the Projected Hartree Fock method of
C.A. Jimenez-Hoyos and G.E. Scuseria:

Jimenez-Hoyos, Carlos A., et al. "Projected hartree–fock theory."
The Journal of chemical physics 136.16 (2012): 164109.

The code is written as a part ("link")
of Gaussian.

To compile:
   1. Get a development version of Gaussian (gdv). The last one tested with this code
   is gdv-h21.

   2. Get PGI fortran compiler. The last one tested with this code is PGI-13.3

   3. Adjust setenv.csh and setenv.sh to include correct locations of Gaussian and PGI

   4. Go to C-shell (csh), source setenv.csh and type "mk". You should get l535.exe
      if you succeed.

To use:
   1. This link replaces what is normally done by links 701-703. In order to compute
      gradients using this link, first use command "testrt "<job keywords>"" to obtain a numeric
      representation of the Gaussian job.

   2. For the input file as follows:
      
       %subst l535 /path/to/l535/
       #p nonstd
       !
       ! Optional description of the route 
       !
       <Output of testrt>

       <Normal Gaussian comment + molecule specification>

   3. Replace occurences of 5/options/1 (HF calculation) with 5/options/30 (PHF calculation)
      and 7/options/1,2 (gradient code for HF) with 5/options/35 (Gradient code for PHF)

   4. The options supported by the gradient code are listed in the beginning of the file.
      Basically it is only possible to select the integration grid and a PHF type
      (should match the one used in the PHF calculation).

      
