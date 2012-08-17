module physbndry_reg_module

  use multifab_module
  use layout_module

  implicit none

  type :: physbndry_reg
     integer :: idim = -1
     integer :: iface = -1
     integer :: nc = -1
     type(multifab) :: data
     type(layout) :: la
  end type physbndry_reg

  interface isValid
     module procedure validFab
  end interface

  private

  public physbndry_reg_build, physbndry_reg_destroy, physbndry_reg, isValid

contains

  subroutine physbndry_reg_build(pbr, la0, nc_in, idim_in, iface_in)
    use bl_error_module
    type(layout), intent(in) :: la0
    type(physbndry_reg), intent(out) :: pbr
    integer, intent(in) :: nc_in, idim_in, iface_in

    type(box) :: rbox, pd
    type(box), allocatable :: bxs(:)
    type(boxarray)         :: baa
    integer :: ibox, nb
    integer ::  lo(la0%lap%dim),  hi(la0%lap%dim)
    integer :: plo(la0%lap%dim), phi(la0%lap%dim)
    logical :: pmask(la0%lap%dim)

    pbr%nc = nc_in
    pbr%idim = idim_in
    pbr%iface = iface_in

    nb = nboxes(la0)
    pd = get_pd(la0)
    pmask = get_pmask(la0)

    plo = lwb(pd)
    phi = upb(pd)

    allocate(bxs(nb))

    do ibox=1,nb
       rbox = get_box(la0,ibox)
       lo = lwb(rbox)
       hi = upb(rbox)

       if (pmask(idim_in)) then
          hi = lo
       else if (iface_in .eq. 0) then
          if (lo(idim_in) .ne. plo(idim_in)) then ! interior
             hi = lo
          else
             hi(idim_in) = lo(idim_in)
          end if
       else if (iface_in .eq. 1) then
          if (hi(idim_in) .ne. phi(idim_in)) then ! interior
             lo = hi
          else
             lo(idim_in) = hi(idim_in)
          end if
       else
          call bl_error("physbndry_reg_build: invalid iface")
       end if
       
       call build(bxs(ibox), lo, hi)
       
    end do

    call build(baa, bxs, sort=.false.)
    call build(pbr%la, baa, boxarray_bbox(baa), explicit_mapping=get_proc(la0))
    call build(pbr%data, pbr%la, nc=nc_in, ng=0)
    call destroy(baa)

    deallocate(bxs)

  end subroutine physbndry_reg_build


  subroutine physbndry_reg_destroy(pbr)
    type(physbndry_reg), intent(inout) :: pbr
    pbr%idim = -1
    pbr%iface = -1
    pbr%nc = -1
    call destroy(pbr%la)
    call destroy(pbr%data)

  end subroutine physbndry_reg_destroy

  function validFab(pbr, ibox) result(r)
    integer, intent(in) :: ibox
    type(physbndry_reg), intent(in) :: pbr
    logical :: r
    type(box) :: bx
    if (remote(pbr%data,ibox)) then
       r = .false.
    else
       bx = get_box(pbr%data,ibox)
       if (box_volume(bx) .eq. 1) then
          r = .false.
       else
          r = .true.
       end if
    end if
  end function validFab

end module physbndry_reg_module
