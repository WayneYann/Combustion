module hypterm_module

  use meth_params_module, only : NVAR, URHO, UMX, UTEMP, difmag
  use reconstruct_module, only : reconstruct
  use riemann_module, only : riemann

  implicit none

  private

  public :: hypterm

contains

  subroutine hypterm(lo,hi,U,Ulo,Uhi,flx,dx)
    integer, intent(in) :: lo(1), hi(1), Ulo(1), Uhi(1)
    double precision, intent(in) :: dx(1)
    double precision, intent(in ) ::   U(Ulo(1):Uhi(1)  ,NVAR)
    double precision, intent(out) :: flx( lo(1): hi(1)+1,NVAR)

    double precision, allocatable :: UL(:,:), UR(:,:)

    integer :: i, n
    double precision, allocatable :: divv(:), v(:)

    allocate(UL(lo(1):hi(1)+1,NVAR))
    allocate(UR(lo(1):hi(1)+1,NVAR))

    call reconstruct(lo(1), hi(1), U, Ulo(1), Uhi(1), UL=UL, UR=UR)
    
    call riemann(lo(1), hi(1), UL, UR, flx)
    
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

end module hypterm_module

