      REAL_T pamb, dpdt_factor
      common /htdata1/pamb, dpdt_factor
      
      integer Density, Temp, RhoH, Trac, FirstSpec, LastSpec
      common /scalID/ Density, Temp, RhoH, Trac, FirstSpec, LastSpec

      REAL_T typVal_Density, typVal_Temp, typVal_RhoH, typVal_Trac, typVal_Y(maxspec), typVal_Vel,
     &     typVal_YMIN, typVal_YMAX
      common /typVals/ typVal_Density, typVal_Temp, typVal_RhoH, typVal_Trac, typVal_Y, typVal_Vel,
     &     typVal_YMIN, typVal_YMAX

      save /htdata1/,/scalID/,/typVals/
