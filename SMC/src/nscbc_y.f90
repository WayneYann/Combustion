
  subroutine outlet_ylo(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(3):hi(3))

    call bl_error("outlet_ylo not implemented")


  end subroutine outlet_ylo


  subroutine inlet_ylo(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,qin,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(3):hi(3))
    double precision,intent(in   )::qin(nqin,lo(1):hi(1),lo(3):hi(3))

    call bl_error("inlet_ylo not implemented")

  end subroutine inlet_ylo


  subroutine outlet_yhi(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(3):hi(3))

    call bl_error("outlet_yhi not implemented")

  end subroutine outlet_yhi


  subroutine inlet_yhi(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,qin,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(3):hi(3))
    double precision,intent(in   )::qin(nqin,lo(1):hi(1),lo(3):hi(3))

    call bl_error("inlet_yhi not implemented")

  end subroutine inlet_yhi

