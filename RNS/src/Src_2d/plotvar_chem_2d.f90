subroutine rns_makeplotvar(lo, hi, dx, &
     prim, pm_l1, pm_l2, pm_h1, pm_h2, &
     plot, pt_l1, pt_l2, pt_h1, pt_h2, &
     npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
     icomp_X, icomp_omegadot, icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
     fuelID)
  use meth_params_module, only : NVAR, QRHO, QU, QV, QTEMP, QFY, NSPEC
  use eos_module, only : eos_get_c
  use chemistry_module, only : molecular_weight, h0 => std_heat_formation
  implicit none

  double precision, intent(in) :: dx(2)
  integer, intent(in) :: lo(2), hi(2)
  integer, intent(in) :: pm_l1, pm_h1, pm_l2, pm_h2
  integer, intent(in) :: pt_l1, pt_h1, pt_l2, pt_h2
  integer, intent(in) :: npv, icomp_magvel, icomp_Mach, icomp_divu, icomp_magvort, &
       icomp_X, icomp_omegadot,	icomp_dYdt, icomp_heatRelease, icomp_fuelConsumption, &
       fuelID
  double precision, intent(in   ) :: prim(pm_l1:pm_h1,pm_l2:pm_h2,NVAR) 
  double precision, intent(inout) :: plot(pt_l1:pt_h1,pt_l2:pt_h2,0:npv-1) 
    
  integer :: i, j, n, iwrk, np
  double precision :: Yt(NSPEC), Xt(NSPEC), v, c, dxinv(2)
  double precision :: Y1d(lo(1):hi(1),nspec), wdot(lo(1):hi(1),nspec), rwrk

  if (icomp_magvel .ge. 0) then
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           plot(i,j,icomp_magvel) = sqrt(prim(i,j,QU)**2+prim(i,j,QV)**2)
        end do
     end do
  end if

  if (icomp_Mach .ge. 0) then
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           Yt = prim(i,j,QFY:QFY+NSPEC-1)
           call eos_get_c(c,prim(i,j,QRHO),prim(i,j,QTEMP),Yt)
           if (icomp_magvel .ge. 0) then
              v = plot(i,j,icomp_magvel)
           else
              v = sqrt(prim(i,j,QU)**2+prim(i,j,QV)**2)
           end if
           plot(i,j,icomp_Mach) = v / c
        end do
     end do
  end if

  if (icomp_divu .ge. 0) then
     dxinv = 1.d0/dx
     do j=lo(2), hi(2)
        do i=lo(1),hi(1)
           plot(i,j,icomp_divu) = 0.5d0 * &
                ( (prim(i+1,j,QU)-prim(i-1,j,QU))*dxinv(1) &
                + (prim(i,j+1,QV)-prim(i,j-1,QV))*dxinv(2) )
        end do
     end do
  end if

  if (icomp_magvort .ge. 0) then
     dxinv = 1.d0/dx
     do j=lo(2), hi(2)
        do i=lo(1),hi(1)
           plot(i,j,icomp_magvort) = 0.5d0 * &
                ( (prim(i+1,j,QV)-prim(i-1,j,QV))*dxinv(1) &
                + (prim(i,j+1,QU)-prim(i,j-1,QU))*dxinv(2) )
        end do
     end do
  end if

  if (icomp_X .ge. 0) then
     do j=lo(2), hi(2)
        do i=lo(1),hi(1)
           Yt = prim(i,j,QFY:QFY+NSPEC-1)
           call ckytx (Yt, iwrk, rwrk, Xt)
           plot(i,j,icomp_X:icomp_X+NSPEC-1) = Xt
        end do
     end do
  end if

  if (icomp_omegadot.ge.0 .or. icomp_dYdt.ge.0 .or.  &
       icomp_heatRelease.ge.0 .or. icomp_fuelConsumption.ge.0) then

     np = hi(1) - lo(1) + 1
     
     do j=lo(2),hi(2)

        do n=1, nspec
           do i=lo(1),hi(1)
              Y1d(i,n) = prim(i,j,QFY+n-1)
           end do
        end do
       
        call vckwyr(np, prim(lo(1),j,QRHO), prim(lo(1),j,QTEMP), Y1d, iwrk, rwrk, wdot)

        if (icomp_omegadot .ge. 0) then
           do n=1, nspec
              do i=lo(1),hi(1)
                 plot(i,j,icomp_omegadot+n-1) = wdot(i,n) * molecular_weight(n)
              end do
           end do
        end if

        if (icomp_dYdt .ge. 0) then
           do n=1, nspec
              do i=lo(1),hi(1)
                 plot(i,j,icomp_dYdt+n-1) = wdot(i,n) * molecular_weight(n) / prim(i,j,QRHO)
              end do
           end do
        end if

        if (icomp_heatRelease .ge. 0) then
           do i=lo(1),hi(1)
              plot(i,j,icomp_heatRelease) = 0.d0
           end do

           do n=1, nspec
              do i=lo(1),hi(1)
                 plot(i,j,icomp_heatRelease) = plot(i,j,icomp_heatRelease) &
                      - h0(n) * wdot(i,n) * molecular_weight(n)
              end do
           end do
        end if

        if (icomp_fuelConsumption .ge. 0) then
           do i=lo(1),hi(1)
              plot(i,j,icomp_fuelConsumption) = -wdot(i,fuelID+1) * molecular_weight(fuelID+1)
           end do
        end if
        
     end do
  end if

end subroutine rns_makeplotvar
