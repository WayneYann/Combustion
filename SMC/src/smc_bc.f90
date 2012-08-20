module smc_bc_module

  use bc_module

  implicit none

  ! boxbclo(ndim,nbox) & boxbchi contain boundray information for each box
  integer, allocatable, save :: boxbclo(:,:), boxbchi(:,:)

  ! datalo(ndim,nbox), datahi(ndim,nbox): the region with good data
  integer, allocatable, save :: datalo(:,:), datahi(:,:)

  private

  public :: smc_bc_init, smc_bc_close, get_boxbc, get_data_lo_hi

contains

  subroutine smc_bc_init(la, U)
    use layout_module
    use multifab_module
    use probin_module, only : bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi
    use derivative_stencil_module, only : stencil_ng
    type(layout), intent(in) :: la
    type(multifab), intent(in) :: U ! xxxxx we will need to save data for inflow

    integer :: ndm, nbx, i, j, bclo(3), bchi(3)
    integer ::  lo(la%lap%dim),  hi(la%lap%dim)
    integer :: plo(la%lap%dim), phi(la%lap%dim)
    type(box) :: pd, bx

    ndm = get_dim(la)
    nbx = nboxes(la)
    pd = get_pd(la)

    plo = lwb(pd)
    phi = upb(pd)
    
    bclo(1) = bcx_lo
    bclo(2) = bcy_lo
    bclo(3) = bcz_lo

    bchi(1) = bcx_hi
    bchi(2) = bcy_hi
    bchi(3) = bcz_hi

    do i=1,ndm
       if ( bclo(i) .ne. PERIODIC .and. &
            bclo(i) .ne. OUTLET ) then
          call bl_error("Unknown boundary")
       end if
       if ( bchi(i) .ne. PERIODIC .and. &
            bchi(i) .ne. OUTLET ) then
          call bl_error("Unknown boundary")
       end if
    end do

    allocate(boxbclo(ndm,nbx))
    allocate(boxbchi(ndm,nbx))
    allocate(datalo(ndm,nbx))
    allocate(datahi(ndm,nbx))

    do j=1,nbx
       bx = get_box(la,j)
       lo = lwb(bx)
       hi = upb(bx)

       do i=1,ndm
          if (lo(i) .eq. plo(i)) then
             boxbclo(i,j) = bclo(i)
          else
             boxbclo(i,j) = INTERIOR
          end if

          if (hi(i) .eq. phi(i)) then
             boxbchi(i,j) = bchi(i)
          else
             boxbchi(i,j) = INTERIOR
          end if

          if (boxbclo(i,j) .eq. PERIODIC .or. boxbclo(i,j) .eq. INTERIOR) then
             datalo(i,j) = lo(i) - stencil_ng
          else
             datalo(i,j) = lo(i)
          end if

          if (boxbchi(i,j) .eq. PERIODIC .or. boxbchi(i,j) .eq. INTERIOR) then
             datahi(i,j) = hi(i) + stencil_ng
          else
             datahi(i,j) = hi(i)
           end if

       end do
    end do    

  end subroutine smc_bc_init


  subroutine smc_bc_close()
    deallocate(boxbclo,boxbchi,datalo,datahi)
  end subroutine smc_bc_close


  subroutine get_boxbc(i,bbclo,bbchi)
    integer, intent(in) :: i
    integer, intent(out) :: bbclo(:), bbchi(:)
    bbclo = boxbclo(:,i)
    bbchi = boxbchi(:,i)
  end subroutine get_boxbc


  subroutine get_data_lo_hi(i, dlo, dhi)
    integer, intent(in) :: i
    integer, intent(out) :: dlo(:), dhi(:)
    dlo = datalo(:,i)
    dhi = datahi(:,i)
  end subroutine get_data_lo_hi
end module smc_bc_module
