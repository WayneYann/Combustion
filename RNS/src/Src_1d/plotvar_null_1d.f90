subroutine rns_makeplotvar(lo, hi, dx, &
     prim, pm_l1, pm_h1,  &
     plot, pt_l1, pt_h1,  &
     npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
     icomp_X, icomp_omegadot, icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
     fuelID)
  use meth_params_module, only : NVAR, QRHO, QU, QTEMP, QFY, NSPEC
  use eos_module, only : eos_get_c
  implicit none

  double precision, intent(in) :: dx(1)
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: pm_l1, pm_h1
  integer, intent(in) :: pt_l1, pt_h1
  integer, intent(in) :: npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
       icomp_X, icomp_omegadot,	icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
       fuelID
  double precision, intent(in   ) :: prim(pm_l1:pm_h1,NVAR) 
  double precision, intent(inout) :: plot(pt_l1:pt_h1,0:npv-1) 
    
  integer :: i
  double precision :: Yt(NSPEC), c, dxinv(1)

  if (icomp_magvel .ge. 0) then
     do i=lo(1), hi(1)
        plot(i,icomp_magvel) = abs(prim(i,QU))
     end do
  end if

  if (icomp_Mach .ge. 0) then
     do i=lo(1), hi(1)
        Yt = prim(i,QFY:QFY+NSPEC-1)
        call eos_get_c(c,prim(i,QRHO),prim(i,QTEMP),Yt)
        plot(i,icomp_Mach) = abs(prim(i,QU) / c)
     end do
  end if

  if (icomp_divu .ge. 0) then
     dxinv = 1.d0/dx
     do i=lo(1),hi(1)
        plot(i,icomp_divu) = (prim(i+1,QU)-prim(i-1,QU))*dxinv(1)*0.5d0
     end do
  end if

end subroutine rns_makeplotvar
