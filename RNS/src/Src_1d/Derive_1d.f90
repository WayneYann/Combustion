!-----------------------------------------------------------------------

subroutine rns_derpres(p,p_l1,p_h1,ncomp_p, &
     &                 u,u_l1,u_h1,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)
  
  use eos_module, only : eos_get_p
  use meth_params_module, only : URHO, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in) :: p_l1,p_h1,ncomp_p
  integer,intent(in) :: u_l1,u_h1,ncomp_u
  integer,intent(in) :: lo(1), hi(1), domlo(1), domhi(1)
  double precision,intent(out) :: p(p_l1:p_h1,ncomp_p)
  double precision,intent(in ) :: u(u_l1:u_h1,ncomp_u)
  double precision,intent(in) :: dx(1), xlo(1), time, dt
  integer,intent(in) :: bc(1,2,ncomp_u), level, grid_no

  double precision :: Y(NSPEC), rhoInv
  integer          :: i,n

  do i = lo(1),hi(1)
     
     rhoInv = 1.d0/u(i,URHO)
     do n = 1,NSPEC
        Y(n)=u(i,UFS+n-1)*rhoInv
     enddo
     
     call eos_get_p(p(i,1), u(i,URHO), u(i,UTEMP), Y)
     
  enddo
  
end subroutine rns_derpres

!-----------------------------------------------------------------------

subroutine rns_dersoundspeed(c,c_l1,c_h1,ncomp_c, &
     &                       u,u_l1,u_h1,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use eos_module, only : eos_get_c
  use meth_params_module, only : URHO, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in):: c_l1,c_h1,ncomp_c
  integer,intent(in):: u_l1,u_h1,ncomp_u
  integer,intent(in):: lo(1), hi(1), domlo(1), domhi(1)
  double precision,intent(out):: c(c_l1:c_h1,ncomp_c)
  double precision,intent(in ):: u(u_l1:u_h1,ncomp_u)
  double precision,intent(in):: dx(1), xlo(1), time, dt
  integer,intent(in):: bc(1,2,ncomp_u), level, grid_no

  double precision :: Y(NSPEC), rhoInv
  integer          :: i,n

  do i = lo(1),hi(1)

     rhoInv = 1.d0/u(i,URHO)
     do n = 1,NSPEC
        Y(n)=u(i,UFS+n-1)*rhoInv
     enddo
     
     call eos_get_c(c(i,1), u(i,URHO), u(i,UTEMP), Y)

  enddo
  
end subroutine rns_dersoundspeed

!-----------------------------------------------------------------------

subroutine rns_dermachnumber(mach,mach_l1,mach_h1,ncomp_mach,&
     &                          u,   u_l1,   u_h1,ncomp_u,   &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use eos_module, only : eos_get_c
  use meth_params_module, only : URHO, UMX, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in):: mach_l1,mach_h1,ncomp_mach
  integer,intent(in)::    u_l1,   u_h1,ncomp_u
  integer,intent(in):: lo(1), hi(1), domlo(1), domhi(1)
  double precision,intent(out)::mach(mach_l1:mach_h1,ncomp_mach)
  double precision,intent(in )::   u(   u_l1:   u_h1,ncomp_u)
  double precision,intent(in) :: dx(1), xlo(1), time, dt
  integer,intent(in) :: bc(1,2,ncomp_u), level, grid_no

  double precision :: c, Y(NSPEC), rhoInv, ux
  integer          :: i,n

  do i = lo(1),hi(1)

     rhoInv = 1.d0/u(i,URHO)
     ux = u(i,UMX)*rhoInv
     do n = 1,NSPEC
        Y(n)=u(i,UFS+n-1)*rhoInv
     enddo

     call eos_get_c(c, u(i,URHO), u(i,UTEMP), Y)

     mach(i,1) = abs(ux) / c

  enddo

end subroutine rns_dermachnumber

!-----------------------------------------------------------------------

subroutine rns_derdivu(divu,div_l1,div_h1,nd, &
     &                  dat,dat_l1,dat_h1,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  implicit none

  integer,intent(in):: lo(1), hi(1)
  integer,intent(in):: div_l1,div_h1,nd
  integer,intent(in):: dat_l1,dat_h1,nc
  integer,intent(in):: domlo(1), domhi(1)
  integer,intent(in):: bc(1,2,nc)
  double precision,intent(in):: delta(1), xlo(1), time, dt
  double precision,intent(out):: divu(div_l1:div_h1,nd)
  double precision,intent(in )::  dat(dat_l1:dat_h1,nc)
  integer,intent(in):: level, grid_no

  integer          :: i
  double precision :: ulo,uhi

  do i = lo(1), hi(1)
     uhi = dat(i+1,2) / dat(i+1,1)
     ulo = dat(i-1,2) / dat(i-1,1)
     divu(i,1) = 0.5d0 * ( (uhi-ulo) / delta(1) )
  end do
  
end subroutine rns_derdivu

!-----------------------------------------------------------------------

subroutine rns_derspec(spec,spec_l1,spec_h1,ns, &
     &                  dat, dat_l1, dat_h1,nd, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  implicit none
  
  integer,intent(in):: lo(1), hi(1)
  integer,intent(in):: spec_l1,spec_h1,ns
  integer,intent(in)::  dat_l1, dat_h1,nd
  integer,intent(in):: domlo(1), domhi(1)
  integer,intent(in):: bc(1,2,nd)
  double precision,intent(in):: delta(1), xlo(1), time, dt
  double precision,intent(out):: spec(spec_l1:spec_h1,ns)
  double precision,intent(in )::  dat( dat_l1: dat_h1,nd)
  integer,intent(in):: level, grid_no
 
  integer i
  
  do i = lo(1), hi(1)
     spec(i,1) = dat(i,2) / dat(i,1)
  end do

end subroutine rns_derspec

!-----------------------------------------------------------------------

subroutine rns_dermolefrac(spec,spec_l1,spec_h1,ns, &
     &                      dat, dat_l1, dat_h1,nd, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  use meth_params_module, only : NSPEC
  use eos_module, only : eos_YtoX

  implicit none
  
  integer,intent(in):: lo(1), hi(1)
  integer,intent(in):: spec_l1,spec_h1,ns
  integer,intent(in)::  dat_l1, dat_h1,nd
  integer,intent(in):: domlo(1), domhi(1)
  integer,intent(in):: bc(1,2,nd)
  double precision,intent(in):: delta(1), xlo(1), time, dt
  double precision,intent(out):: spec(spec_l1:spec_h1,ns)
  double precision,intent(in )::  dat( dat_l1: dat_h1,nd)
  integer,intent(in):: level, grid_no
 
  integer :: i
  double precision :: Yt(NSPEC),Xt(NSPEC), rhoInv

  do i = lo(1), hi(1)
     rhoInv = 1.d0 / dat(i,1)
     Yt = dat(i,2:NSPEC+1) * rhoInv
     
     call eos_YtoX(Yt, Xt)
     
     spec(i,:) = Xt
  end do

end subroutine rns_dermolefrac


!-----------------------------------------------------------------------

subroutine rns_dervel(vel,vel_l1,vel_h1,nv, &
     &                dat,dat_l1,dat_h1,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will derive the velocity from the momentum.
  !
  implicit none

  integer,intent(in):: lo(1), hi(1)
  integer,intent(in):: vel_l1,vel_h1,nv
  integer,intent(in):: dat_l1,dat_h1,nc
  integer,intent(in):: domlo(1), domhi(1)
  integer,intent(in):: bc(1,2,nc)
  double precision,intent(in):: delta(1), xlo(1), time, dt
  double precision,intent(out):: vel(vel_l1:vel_h1,nv)
  double precision,intent(in ):: dat(dat_l1:dat_h1,nc)
  integer,intent(in):: level, grid_no
 
  integer i

  do i = lo(1), hi(1)
     vel(i,1) = dat(i,2) / dat(i,1)
  end do

end subroutine rns_dervel

