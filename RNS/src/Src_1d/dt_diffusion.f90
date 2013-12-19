subroutine rns_estdt_diff(u,u_l1,u_h1,lo,hi,dx,dt)
  use meth_params_module
  use egz_module, only : egzini, egzpar, egzvr1
  use chemistry_module, only : inv_mwt
  integer, intent(in) :: u_l1,u_h1
  integer, intent(in) :: lo(1), hi(1)
  double precision, intent(in) :: u(u_l1:u_h1,NVAR)
  double precision, intent(in) :: dx(1)
  double precision, intent(inout) :: dt

  integer :: i, n, np, iwrk
  double precision :: rwrk, Yt(NSPEC), Xt(NSPEC), rhoInv, Wbar, D, dx2
  double precision, allocatable :: DZ(:,:), XZ(:,:), CPZ(:,:)
  
  np = hi(1)-lo(1)+1
  
  dx2 = 0.5*dx(1)**2
  
  allocate(DZ (lo(1):hi(1),nspec))
  allocate(XZ (lo(1):hi(1),nspec))
  allocate(CPZ(lo(1):hi(1),nspec))
  
  CPZ = 0.d0
  
  call egzini(np)
  
  do i=lo(1),hi(1)
     rhoInv = 1.d0/u(i,URHO)
     Yt = u(i,UFS:UFS+nspec-1)*rhoInv
     call ckytx(Yt, iwrk, rwrk, Xt)
     XZ(i,:) = Xt
  end do

  call egzpar(u(lo(1):hi(1),UTEMP), XZ, CPZ)
  
  call egzvr1(u(lo(1):hi(1),UTEMP), DZ)
  
  do i=lo(1),hi(1)
     Xt = XZ(i,:)
     call ckmmwx(Xt, iwrk, rwrk, Wbar)
     D = 0.d0
     do n=1,nspec
        D = max(D, DZ(i,n)*Wbar*inv_mwt(n))
     end do
     dt = min(dt, dx2/D*u(i,URHO))
  end do
     
  deallocate(DZ,XZ,CPZ)
  
end subroutine rns_estdt_diff


