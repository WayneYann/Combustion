
subroutine initialize_pmf(filename)
  use pmf_module
  character (len=*) :: filename
  pmf_filename = filename
end subroutine initialize_pmf

subroutine pmf(xlo,xhi,y_vector)
  use pmf_module
  double precision :: xlo,xhi,y_vector(*)
  call interp_pmf(xlo,xhi,y_vector)
end subroutine pmf
