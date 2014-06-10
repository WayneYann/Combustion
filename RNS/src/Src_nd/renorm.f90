module renorm_module

  implicit none

contains

  subroutine floor_species(nspec, Y)
    integer, intent(in) :: nspec
    double precision, intent(inout) :: Y(nspec)
    integer :: n
    do n=1,nspec
       Y(n) = max(Y(n), 0.d0)
    end do
    Y = Y / sum(Y)
  end subroutine floor_species

  subroutine JBBhack(nspec, Y, ierr)
    integer, intent(in) :: nspec
    double precision, intent(inout) :: Y(nspec)
    integer, intent(out), optional :: ierr

    integer          :: n, int_dom_spec
    logical          :: any_negative
    double precision :: dom_spec
    double precision, parameter :: eps = -1.d-16
    
    any_negative = .false.

    ! First deal with tiny undershoots by just setting them to zero
    do n=1,nspec
       if (Y(n) .lt. 0.d0) then
          if (Y(n) .gt. eps) then
             Y(n) = 0.d0
          else
             any_negative = .true.
          end if
       end if
    end do

    ! We know there are one or more undershoots needing correction 
    if (any_negative) then

       ! Find the dominant species
       dom_spec = 0.d0
       int_dom_spec = 0
       do n = 1, nspec
          if (Y(n) .gt. dom_spec) then
             dom_spec = Y(n)
             int_dom_spec = n
          end if
       end do

       ! Now take care of undershoots greater in magnitude than 1e-16.
       do n = 1, nspec
          if (Y(n) .lt. 0.d0) then
             ! Take enough from the dominant species to fill the negative one.
             Y(int_dom_spec) = Y(int_dom_spec) + Y(n)
             
             ! Test that we didn't make the dominant species negative
             if (Y(int_dom_spec) .lt. 0.d0) then 
                print *,'renorm: Just made dominant species negative ',int_dom_spec
                print *,'We were fixing species ',n,' which had value ',Y(n)
                print *,'Dominant species became ',Y(int_dom_spec)
                if (present(ierr)) then
                   ierr = 1
                   return
                else
                   call bl_error("renorm")
                end if
             end if
             
             ! Now set the negative species to zero
             Y(n) = 0.d0
          end if
       end do
       
    end if

    Y = Y / sum(Y)
    ierr = 0

  end subroutine JBBhack

end module renorm_module
