subroutine rns_makeplotvar(lo, hi, dx, &
     prim, pm_l1, pm_h1,  &
     plot, pt_l1, pt_h1,  &
     npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
     icomp_X, icomp_omegadot, icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
     fuelID)
  use meth_params_module, only : NVAR, QRHO, QU, QTEMP, QFY, NSPEC
  use eos_module, only : eos_get_c
  use chemistry_module, only : molecular_weight, h0 => std_heat_formation
  implicit none

  double precision, intent(in) :: dx(1)
  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: pm_l1, pm_h1
  integer, intent(in) :: pt_l1, pt_h1
  integer, intent(in) :: npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
       icomp_X, icomp_omegadot, icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
       fuelID
  double precision, intent(in   ) :: prim(pm_l1:pm_h1,NVAR) 
  double precision, intent(inout) :: plot(pt_l1:pt_h1,0:npv-1) 
    
  integer :: i, n, iwrk, np
  double precision :: Yt(NSPEC), Xt(NSPEC), c, dxinv(1)
  double precision :: Y1d(lo(1):hi(1),nspec), wdot(lo(1):hi(1),nspec), rwrk

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

  if (icomp_X .ge. 0) then
     do i=lo(1),hi(1)
        Yt = prim(i,QFY:QFY+NSPEC-1)
        call ckytx (Yt, iwrk, rwrk, Xt)
        plot(i,icomp_X:icomp_X+NSPEC-1) = Xt
     end do
  end if

  if (icomp_omegadot.ge.0 .or. icomp_dYdt.ge.0 .or.  &
       icomp_heatRelease.ge.0 .or. icomp_fuelConsumption.ge.0) then

     np = hi(1) - lo(1) + 1

     do n=1, nspec
        do i=lo(1),hi(1)
           Y1d(i,n) = prim(i,QFY+n-1)
        end do
     end do
       
     call vckwyr(np, prim(lo(1),QRHO), prim(lo(1),QTEMP), Y1d, iwrk, rwrk, wdot)

     if (icomp_omegadot .ge. 0) then
        do n=1, nspec
           do i=lo(1),hi(1)
              plot(i,icomp_omegadot+n-1) = wdot(i,n) * molecular_weight(n)
           end do
        end do
     end if

     if (icomp_dYdt .ge. 0) then
        do n=1, nspec
           do i=lo(1),hi(1)
              plot(i,icomp_dYdt+n-1) = wdot(i,n) * molecular_weight(n) / prim(i,QRHO)
           end do
        end do
     end if

     if (icomp_heatRelease .ge. 0) then
        do i=lo(1),hi(1)
           plot(i,icomp_heatRelease) = 0.d0
        end do

        do n=1, nspec
           do i=lo(1),hi(1)
              plot(i,icomp_heatRelease) = plot(i,icomp_heatRelease) &
                   - h0(n) * wdot(i,n) * molecular_weight(n)
           end do
        end do
     end if

     if (icomp_fuelConsumption .ge. 0) then
        do i=lo(1),hi(1)
           plot(i,icomp_fuelConsumption) = -wdot(i,fuelID+1) * molecular_weight(fuelID+1)
        end do
     end if

  end if

end subroutine rns_makeplotvar


