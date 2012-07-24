!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types

  implicit none

  integer, save :: irho, imx, imy, imz, iene, iry1
  integer, save :: qrho, qu, qv, qw, qpres, qtemp, qy1

  ! the total number of plot components
  integer, save :: n_plot_comps = 0

  integer, save :: ncons, nprim

contains

  function get_next_plot_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_plot_comps + 1
    n_plot_comps = n_plot_comps + num

    return
  end function get_next_plot_index

  subroutine init_variables()

    use probin_module, only: dm_in, nspec

    irho = 1
    imx = 2
    imy = 3
    if (dm_in .eq. 3) then
       imz = 4
       iene = 5
    else
       imz = -1
       iene = 4
       ncons = 5
    end if
    iry1 = iene+1

    ncons = iry1-1 + nspec

    qrho = 1
    qu = 2
    qv = 3
    if (dm_in .eq. 3) then
       qw = 4
       qpres = 5
       qtemp = 6
       qy1 = 7
    else
       qw = -1
       qpres = 4
       qtemp = 5
       qy1 = 7
    end if

    nprim = qy1-1 + nspec
    
  end subroutine init_variables

  subroutine init_plot_variables()

    return

  end subroutine init_plot_variables

end module variables
