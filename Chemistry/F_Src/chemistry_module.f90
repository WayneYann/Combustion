module chemistry_module

  use bl_types
  use eglib_module

  implicit none

  integer, save :: nelements   ! number of elements
  integer, save :: nspecies    ! number of species
  integer, save :: nreactions  ! number of reactions

  logical, save :: chemistry_initialized = .false.

  character*2, allocatable, save :: elem_names(:)
  character*4, allocatable, save :: spec_names(:)

  double precision, allocatable, save :: molecular_weight(:)

contains

  subroutine chemistry_init(eglib_np)
    integer, intent(in) :: eglib_np

    integer :: iwrk, nfit, i, ic, ii
    double precision :: rwrk
    integer, allocatable :: names(:)

    call ckindx(iwrk, rwrk, nelements, nspecies, nreactions, nfit)

    allocate(elem_names(nelements))
    allocate(spec_names(nspecies))
    allocate(molecular_weight(nspecies))

    allocate(names(nspecies*4))  ! Each species name has at most 4 characters

    call cksyme(names, 2)  ! Two chars for element names

    ic = 1
    do i = 1, nelements
       do ii=1, 2
          elem_names(i)(ii:ii) = char(names(ic))
          ic = ic + 1
       end do
    end do

    call cksyms(names, 4) ! Four chars for species names

    ic = 1
    do i = 1, nspecies
       do ii=1, 4
          spec_names(i)(ii:ii) = char(names(ic))
          ic = ic+1
       end do
    end do

    deallocate(names)

    call ckwt(iwrk, rwrk, molecular_weight)

!    call eglib_init(nspecies, eglib_np)

    chemistry_initialized = .true.

  end subroutine chemistry_init


  subroutine chemistry_close()

    deallocate(elem_names,spec_names,molecular_weight)

!    call eglib_close()

  end subroutine chemistry_close

end module chemistry_module
