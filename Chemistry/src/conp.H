      integer NP, NRHO, NWT, NWTI, NZ, iN2, NEQ, nchemdiag, negative_Y_test
      common /conpcom/ NP, NRHO, NWT, NWTI, NZ, iN2, NEQ, nchemdiag, negative_Y_test
      save /conpcom/

      REAL_T typVal_Density, typVal_Temp, typVal_RhoH, typVal_Trac, typVal_Y(maxspec), typVal_Vel,
     &     typVal_YMIN, typVal_YMAX
      common /typVals/ typVal_Density, typVal_Temp, typVal_RhoH, typVal_Trac, typVal_Y, typVal_Vel,
     &     typVal_YMIN, typVal_YMAX

      save /typVals/
