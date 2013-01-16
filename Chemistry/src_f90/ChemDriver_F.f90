subroutine cd_initchem(len_en, len_sn, ne, ns)
  use chemistry_module
  implicit none
  integer, intent(out) :: len_en, len_sn, ne, ns
  call chemistry_init()
  len_en = L_elem_name
  len_sn = L_spec_name
  ne = nelements
  ns = nspecies
end subroutine cd_initchem


subroutine cd_closechem()

  use chemistry_module

  implicit none

  call chemistry_close()

end subroutine cd_closechem


subroutine cd_getelemname(ielem, lname, name)
  use chemistry_module
  implicit none
  integer, intent(in)  :: ielem
  integer, intent(out) :: lname
  integer, intent(inout) :: name(L_elem_name)

  integer :: i

  lname = len_trim(elem_names(ielem))

  do i = 1, lname
     name(i) = ichar(elem_names(ielem)(i:i))
  end do
end subroutine cd_getelemname


subroutine cd_getspecname(ispec, lname, name)
  use chemistry_module
  implicit none
  integer, intent(in)  :: ispec
  integer, intent(out) :: lname
  integer, intent(inout) :: name(L_spec_name)

  integer :: i

  lname = len_trim(spec_names(ispec))

  do i = 1, lname
     name(i) = ichar(spec_names(ispec)(i:i))
  end do
end subroutine cd_getspecname
