module chemistry_module

  implicit none

  integer, save :: nelements   ! number of elements
  integer, save :: nspecies    ! number of species
  integer, save :: nreactions  ! number of reactions

  logical, save :: chemistry_initialized = .false.

  integer, private, parameter :: L_elem_name = 3 ! Each element name has at most 3 characters
  character*(L_elem_name), allocatable, save :: elem_names(:)

  integer, private, parameter :: L_spec_name = 16 ! Each species name has at most 8 characters
  character*(L_spec_name), allocatable, save :: spec_names(:)

  double precision, allocatable, save :: molecular_weight(:), inv_mwt(:)

  double precision, save :: Ru, Ruc, Patm

  double precision, allocatable, save :: std_heat_formation(:)

contains

  subroutine chemistry_init()
    integer :: iwrk, nfit, i, ic, ii
    double precision :: rwrk, T0
    integer, allocatable :: names(:)

    call ckindx(iwrk, rwrk, nelements, nspecies, nreactions, nfit)

    allocate(elem_names(nelements))
    allocate(spec_names(nspecies))
    allocate(molecular_weight(nspecies))
    allocate(inv_mwt(nspecies))

    allocate(names(nspecies*L_spec_name))  

    allocate(std_heat_formation(nspecies))

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
    inv_mwt = 1.d0 / molecular_weight

    call ckrp(iwrk, rwrk, Ru, Ruc, Patm)

    T0 = 298.15d0
    call ckhms(T0, iwrk, rwrk, std_heat_formation)

    chemistry_initialized = .true.

  end subroutine chemistry_init


  subroutine chemistry_close()
    deallocate(elem_names,spec_names,molecular_weight,inv_mwt,std_heat_formation)
  end subroutine chemistry_close


  function get_species_index(name) result (iname)
    character(len=*), intent(in) :: name
    integer :: iname
    integer :: i
    iname = -1
    do i = 1, nspecies
       if (trim(spec_names(i)) .eq. trim(name)) then
          iname = i
          exit
       end if
    end do
  end function get_species_index


end module chemistry_module
