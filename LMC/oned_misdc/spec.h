      integer maxspec, maxspnml
      parameter (maxspec = 20)
      parameter (maxspnml = 16)

      integer eg_nodes, eg_IFLAG, eg_ITLS, egr, egi, LLINKMC
      parameter (eg_nodes = 1)
      parameter (eg_IFLAG = 7)
      parameter (eg_ITLS  = 3)
      parameter ( egr = 23 + 14*maxspec + 32*maxspec**2 + 13*eg_nodes
     &     + 30*eg_nodes*maxspec + 5*eg_nodes*maxspec**2)
      parameter (egi = maxspec)
      parameter (LLINKMC = 56)
      
      integer traninit, EGIWRK(egi), Nelt, Nspec, Nreac, Nfit, LeEQ1
      common / trani / EGIWRK, traninit, Nelt, Nspec, Nreac, Nfit, LeEQ1
      double precision EGRWRK(egr), invmwt(maxspec), P1ATM, TMIN_TRANS,
     &     Pr, Sc, thickFacTR
      common / tranr / EGRWRK, invmwt, P1ATM, TMIN_TRANS, Pr, Sc,
     &     thickFacTR
      character*256 tranfile
      common / tranc / tranfile
