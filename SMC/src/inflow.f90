  subroutine store_inflow(qin, U)
    type(physbndry_reg), intent(inout) :: qin
    type(multifab), intent(in) :: U

    integer :: n,ng
    integer ::  lo(U%dim),  hi(U%dim)
    integer :: qlo(U%dim), qhi(U%dim)
    double precision, pointer, dimension(:,:,:,:) :: qp, up
    
    ng = nghost(U)

    do n=1,nboxes(U)
       if (isValid(qin,n)) then
          lo = lwb(get_box(U,n))
          hi = upb(get_box(U,n))

          qlo = lwb(get_box(qin%data,n))
          qhi = upb(get_box(qin%data,n))

          qp => dataptr(qin%data, n)
          up => dataptr(U,n)

          call ctoprim_inflow(lo,hi,ng,up,qlo,qhi,qp)
       end if
    end do

  end subroutine store_inflow

