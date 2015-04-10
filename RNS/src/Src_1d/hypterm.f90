module hypterm_module

  use riemann_module, only : riemann

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,flx,dx)
    use meth_params_module, only : NVAR, URHO, UMX, UTEMP, difmag
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in) :: dx(1)
    double precision, intent(in ) ::   U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: flx( lo(1): hi(1)+1,NVAR)

    double precision, allocatable :: UL(:,:), UR(:,:)

    integer :: i, n
    double precision, allocatable :: divv(:), v(:)

    allocate(UL(lo(1):hi(1)+1,NVAR))
    allocate(UR(lo(1):hi(1)+1,NVAR))

    call reconstruct(lo(1), hi(1), U, Ulo(1), Uhi(1), UL, UR, lo(1), hi(1)+1)

    call riemann(lo(1),hi(1), UL, UR, lo(1),hi(1)+1, flx, lo(1),hi(1)+1)
    
    deallocate(UL,UR)

    if (difmag .gt. 0.d0) then

       allocate(v   (lo(1)-1:hi(1)+1))
       allocate(divv(lo(1)  :hi(1)+1))

       do i=lo(1)-1,hi(1)+1
          v(i) = U(i,UMX)/U(i,URHO)
       end do

       do i=lo(1),hi(1)+1
          divv(i) = difmag * min(0.d0, v(i)-v(i-1))
       end do

       do n=1,NVAR
          if (n.ne.UTEMP) then
             do i = lo(1),hi(1)+1
                flx(i,n) = flx(i,n) + divv(i)*(U(i,n) - U(i-1,n))
             end do
          end if
       end do

       deallocate(v,divv)
       
    end if

  end subroutine hypterm

  subroutine reconstruct(lo, hi, U, Ulo, Uhi, UL, UR, flo, fhi)
    use meth_params_module, only : NVAR, URHO, UMX, UTEMP, UFS, UEDEN, NSPEC, NCHARV, CFS
    use reconstruct_module, only : get_eigen_matrices_q
    use renorm_module, only : floor_species
    use eos_module, only : eos_get_eref
    use mdcd_module, only : mdcd
    integer, intent(in) :: lo, hi, Ulo, Uhi, flo, fhi
    double precision, intent(in) :: U(Ulo:Uhi,NVAR)
    double precision, intent(out) :: UL(flo:fhi,NVAR)
    double precision, intent(out) :: UR(flo:fhi,NVAR)

    integer :: i, n, ii, ivar, m
    double precision, allocatable :: Y(:,:), RoeW(:), v(:)
    double precision :: rhoInv, rho0, Y0(nspec), T0, v0(3), fac, eref
    double precision :: egv1(NCHARV,NCHARV), egv2(NCHARV,NCHARV), Uii(NCHARV), charv(-3:2,NCHARV)
    double precision :: cvl(NCHARV), cvr(NCHARV)

    allocate(Y(lo-3:hi+3,nspec))
    allocate(RoeW(lo-3:hi+3))
    allocate(v(lo-3:hi+3))

    v0 = 0.d0

    UL = 0.d0
    UR = 0.d0

    do i=lo-3,hi+3
       rhoInv = 1.d0/U(i,URHO)

       RoeW(i) = sqrt(U(i,URHO))

       do n=1,nspec
          Y(i,n) = U(i,UFS+n-1)*rhoInv*RoeW(i)
       end do

       v(i) = U(i,UMX)*rhoInv*RoeW(i)
    end do
    
    do i=lo, hi+1
       rho0 = RoeW(i-1)*RoeW(i)
       fac = 1.d0 / (RoeW(i-1)+RoeW(i))
       Y0 = (Y(i-1,:) + Y(i,:)) * fac
       T0 = (U(i-1,UTEMP)*RoeW(i-1) + U(i,UTEMP)*RoeW(i)) * fac
       v0(1) = (v(i-1) + v(i)) * fac

       call floor_species(nspec, Y0)

       ! egv1: left matrix;  egv2: right matrix
       call get_eigen_matrices_q(rho0, Y0, T0, v0, egv1, egv2)

       do ii=-3,2
          Uii(1) = U(i+ii,UMX)
          Uii(2) = 0.d0  ! because of 1d
          Uii(3) = 0.d0
          Uii(4) = U(i+ii,UEDEN)
          
          eref = eos_get_eref(Y(i+ii,:))
          Uii(4) = Uii(4) - U(i+ii,URHO) * eref

          do n=1,nspec
             Uii(CFS+n-1) = U(i+ii,UFS+n-1)
          end do

          do n=1,NCHARV
             charv(ii,n) = dot_product(egv1(:,n),Uii)
          end do
       end do

       do ivar=1,NCHARV
          call mdcd(charv(:,ivar), cvl(ivar), cvr(ivar))
       end do

       egv1 = transpose(egv2)  ! egv1 now holds transposed right matrix

       do n=1,NCHARV
          UL(i,UMX  ) = UL(i,UMX  ) + cvl(n)*egv2(1,n)
          UR(i,UMX  ) = UR(i,UMX  ) + cvr(n)*egv2(1,n)
          UL(i,UEDEN) = UL(i,UEDEN) + cvl(n)*egv2(4,n)
          UR(i,UEDEN) = UR(i,UEDEN) + cvr(n)*egv2(4,n)
       end do

       do m=1,nspec
          UL(i,UFS+m-1) = UL(i,UFS+m-1) + dot_product(cvl, egv1(:,CFS+m-1))
          UR(i,UFS+m-1) = UR(i,UFS+m-1) + dot_product(cvr, egv1(:,CFS+m-1))
          UL(i,URHO   ) = UL(i,URHO) + UL(i,UFS+m-1)
          UR(i,URHO   ) = UR(i,URHO) + UR(i,UFS+m-1)
       end do

       UL(i,UTEMP) = U(i-1,UTEMP)
       UR(i,UTEMP) = U(i  ,UTEMP)
          
       rhoInv = 1.d0/UL(i,URHO)
       do n=1,nspec
          Y0(n) = UL(i,UFS+n-1) * rhoInv
       end do
       call floor_species(nspec, Y0)
       eref = eos_get_eref(Y0)
       UL(i,UEDEN) = UL(i,UEDEN) + UL(i,URHO) * eref
       do n=1,nspec
          UL(i,UFS+n-1) = UL(i,URHO)*Y0(n)
       end do

       rhoInv = 1.d0/UR(i,URHO)
       do n=1,nspec
          Y0(n) = UR(i,UFS+n-1) * rhoInv
       end do
       call floor_species(nspec, Y0)
       eref = eos_get_eref(Y0)
       UR(i,UEDEN) = UR(i,UEDEN) + UR(i,URHO) * eref
       do n=1,nspec
          UR(i,UFS+n-1) = UR(i,URHO)*Y0(n)
       end do
    end do

  end subroutine reconstruct

end module hypterm_module

