module chemistry_module

  use bl_types

  implicit none

  integer, save :: nelements   ! number of elements
  integer, save :: nspecies    ! number of species
  integer, save :: nreactions  ! number of reactions

  logical, save :: chemistry_initialized = .false.

  integer, private, parameter :: L_elem_name = 3 ! Each element name has at most 3 characters
  character*(L_elem_name), allocatable, save :: elem_names(:)

  integer, private, parameter :: L_spec_name = 8 ! Each species name has at most 8 characters
  character*(L_spec_name), allocatable, save :: spec_names(:)

  double precision, allocatable, save :: molecular_weight(:)

contains

  subroutine chemistry_init()
    integer :: iwrk, nfit, i, ic, ii
    double precision :: rwrk
    integer, allocatable :: names(:)

    call ckindx(iwrk, rwrk, nelements, nspecies, nreactions, nfit)

    allocate(elem_names(nelements))
    allocate(spec_names(nspecies))
    allocate(molecular_weight(nspecies))

    allocate(names(nspecies*L_spec_name))  

    call cksyme(names, L_elem_name) 

    ic = 1
    do i = 1, nelements
       do ii=1, L_elem_name
          elem_names(i)(ii:ii) = char(names(ic))
          ic = ic + 1
       end do
    end do

    call cksyms(names, L_spec_name) 

    ic = 1
    do i = 1, nspecies
       do ii=1, L_spec_name
          spec_names(i)(ii:ii) = char(names(ic))
          ic = ic+1
       end do
    end do

    deallocate(names)

    call ckwt(iwrk, rwrk, molecular_weight)

    chemistry_initialized = .true.

  end subroutine chemistry_init


  subroutine chemistry_close()
    deallocate(elem_names,spec_names,molecular_weight)
  end subroutine chemistry_close

end module chemistry_module
