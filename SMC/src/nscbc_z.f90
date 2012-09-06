
  subroutine outlet_zlo(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(2):hi(2))

    call bl_error("outlet_zlo not implemented")


  end subroutine outlet_zlo


  subroutine inlet_zlo(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,qin,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(2):hi(2))
    double precision,intent(in   )::qin(nqin,lo(1):hi(1),lo(2):hi(2))

    call bl_error("inlet_zlo not implemented")

  end subroutine inlet_zlo


  subroutine outlet_zhi(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(2):hi(2))

    call bl_error("outlet_zhi not implemented")

  end subroutine outlet_zhi


  subroutine inlet_zhi(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,qin,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(3):hi(2))
    double precision,intent(in   )::qin(nqin,lo(1):hi(1),lo(3):hi(2))

    call bl_error("inlet_zhi not implemented")

  end subroutine inlet_zhi

