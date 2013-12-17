module debug
  implicit none
  integer, parameter :: un = 66
  logical, save      :: connected = .false.
contains

  subroutine plotX(q, window, options, wait)
    real(8),      intent(in   ) :: q(:,:)
    integer,      intent(in   ) :: window
    character(*), intent(in   ) :: options
    logical,      intent(in   ) :: wait

    integer :: i

#ifdef GNUPLOT
    if (.not. connected) then
       open(unit=un, file="/tmp/gp")
       connected = .true.
    end if

    write(un,*) "set term wxt ", window
    write(un,*) "plot '-' ", options
    do i = 1, size(q)
       write(un,*) i, q(:,i)
    end do
    write(un,"(a)"), "e"

    if (wait) then
       write (*,*) '==> paused'
       read  (*,*)
    end if
#endif
  end subroutine plotX

  subroutine plot1(q1, window, options, wait)
    real(8),      intent(in   ) :: q1(:)
    integer,      intent(in   ) :: window
    character(*), intent(in   ) :: options
    logical,      intent(in   ) :: wait
    real(8) :: q(1,size(q1))
    q(1,:) = q1
    call plotX(q, window, options, wait)
  end subroutine plot1

  subroutine plot2(q1, q2, window, options, wait)
    real(8),      intent(in   ) :: q1(:), q2(:)
    integer,      intent(in   ) :: window
    character(*), intent(in   ) :: options
    logical,      intent(in   ) :: wait
    real(8) :: q(2,size(q1))
    q(1,:) = q1
    q(2,:) = q2
    call plotX(q, window, options, wait)
  end subroutine plot2

end module debug
