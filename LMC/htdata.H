      REAL_T pamb, dpdt_factor
      common /htdata1/pamb, dpdt_factor
      
      integer Density, Temp, RhoH, Trac, FirstSpec, LastSpec
      common /scalID/ Density, Temp, RhoH, Trac, FirstSpec, LastSpec

      save /htdata1/,/scalID/
