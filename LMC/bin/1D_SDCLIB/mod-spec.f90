module spec_module
  integer nx,nfine,nlevs,rr,nscal,Density,Temp, RhoH, Nspec, FirstSpec
  common / amri / nx,nfine,nlevs,rr
      common / speci / Nelt, Nspec, Nreac, Nfit, iH2, iO2, iCH4, &
           iCH3OCH3, iCO2, iH2O, iN2, Density, Temp, RhoH, &
          RhoRT, FirstSpec, LastSpec, nscal

  ! common / speci / nscal, Density, Temp, RhoH, Nspec, FirstSpec
  ! save /amri/
end module spec_module
