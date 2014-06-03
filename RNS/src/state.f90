module state

  double precision, save :: time
  integer, save :: level, iteration

end module state

subroutine rns_set_state(t,l,k)
  use state
  implicit none
  double precision, intent(in) :: t
  integer, intent(in) :: l, k
  time = t
  level = l
  iteration = k
end subroutine rns_set_state

