!-----------------------------------------------------------------------

subroutine rns_derpres(p,p_l1,p_l2,p_h1,p_h2,ncomp_p, &
     &                 u,u_l1,u_l2,u_h1,u_h2,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use eos_module, only : eos_get_p
  use meth_params_module, only : URHO, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in) :: p_l1,p_l2,p_h1,p_h2,ncomp_p
  integer,intent(in) :: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer,intent(in) :: lo(2), hi(2), domlo(2), domhi(2)
  double precision,intent(out) :: p(p_l1:p_h1,p_l2:p_h2,ncomp_p)
  double precision,intent(in ) :: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  double precision,intent(in) :: dx(2), xlo(2), time, dt
  integer,intent(in) :: bc(2,2,ncomp_u), level, grid_no

  double precision :: Y(NSPEC), rhoInv
  integer          :: i,j,n

  do j = lo(2),hi(2)
     do i = lo(1),hi(1)
        
        rhoInv = 1.d0/u(i,j,URHO)
        do n = 1,NSPEC
           Y(n)=u(i,j,UFS+n-1)*rhoInv
        enddo
        
        call eos_get_p(p(i,j,1), u(i,j,URHO), u(i,j,UTEMP), Y)
        
     enddo
  enddo
  
end subroutine rns_derpres

!-----------------------------------------------------------------------

subroutine rns_dersoundspeed(c,c_l1,c_l2,c_h1,c_h2,ncomp_c, &
     &                       u,u_l1,u_l2,u_h1,u_h2,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use eos_module, only : eos_get_c
  use meth_params_module, only : URHO, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in):: c_l1,c_l2,c_h1,c_h2,ncomp_c
  integer,intent(in):: u_l1,u_l2,u_h1,u_h2,ncomp_u
  integer,intent(in):: lo(2), hi(2), domlo(2), domhi(2)
  double precision,intent(out):: c(c_l1:c_h1,c_l2:c_h2,ncomp_c)
  double precision,intent(in ):: u(u_l1:u_h1,u_l2:u_h2,ncomp_u)
  double precision,intent(in):: dx(2), xlo(2), time, dt
  integer,intent(in):: bc(2,2,ncomp_u), level, grid_no

  double precision :: Y(NSPEC), rhoInv
  integer          :: i,j,n

  do j = lo(2),hi(2)
     do i = lo(1),hi(1)

        rhoInv = 1.d0/u(i,j,URHO)
        do n = 1,NSPEC
           Y(n)=u(i,j,UFS+n-1)*rhoInv
        enddo
        
        call eos_get_c(c(i,j,1), u(i,j,URHO), u(i,j,UTEMP), Y)
        
     enddo
  enddo
  
end subroutine rns_dersoundspeed

!-----------------------------------------------------------------------

subroutine rns_dermachnumber(mach,mach_l1,mach_l2,mach_h1,mach_h2,ncomp_mach,&
     &                          u,   u_l1,   u_l2,   u_h1,   u_h2,ncomp_u,   &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use eos_module, only : eos_get_c
  use meth_params_module, only : URHO, UMX, UMY, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in):: mach_l1,mach_l2,mach_h1,mach_h2,ncomp_mach
  integer,intent(in)::    u_l1,   u_l2,   u_h1,   u_h2,ncomp_u
  integer,intent(in):: lo(2), hi(2), domlo(2), domhi(2)
  double precision,intent(out)::mach(mach_l1:mach_h1,mach_l2:mach_h2,ncomp_mach)
  double precision,intent(in )::   u(   u_l1:   u_h1,   u_l2:   u_h2,ncomp_u)
  double precision,intent(in) :: dx(2), xlo(2), time, dt
  integer,intent(in) :: bc(2,2,ncomp_u), level, grid_no

  double precision :: c, Y(NSPEC), rhoInv, ux, uy
  integer          :: i,j,n

  do j = lo(2),hi(2)
     do i = lo(1),hi(1)
        
        rhoInv = 1.d0/u(i,j,URHO)
        ux = u(i,j,UMX)*rhoInv
        uy = u(i,j,UMY)*rhoInv
        do n = 1,NSPEC
           Y(n)=u(i,j,UFS+n-1)*rhoInv
        enddo
        
        call eos_get_c(c, u(i,j,URHO), u(i,j,UTEMP), Y)

        mach(i,j,1) = sqrt(ux**2 + uy**2) / c

     enddo
  enddo

end subroutine rns_dermachnumber

!-----------------------------------------------------------------------

subroutine rns_dermagvort(vort,  v_l1,  v_l2,  v_h1,  v_h2,nv, & 
     &                     dat,dat_l1,dat_l2,dat_h1,dat_h2,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will calculate vorticity
  !     
  implicit none

  integer,intent(in):: lo(2), hi(2)
  integer,intent(in)::  v_l1,  v_l2,  v_h1,  v_h2,nv
  integer,intent(in)::dat_l1,dat_l2,dat_h1,dat_h2,nc
  integer,intent(in):: domlo(2), domhi(2), level, grid_no
  integer,intent(in):: bc(2,2,nc)
  double precision,intent(in):: delta(2), xlo(2), time, dt
  double precision,intent(out):: vort(  v_l1:  v_h1,  v_l2:  v_h2,nv)
  double precision,intent(in )::  dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)

  integer          :: i,j
  double precision :: uy,vx, dxinv(2),rhoinv
  double precision, allocatable :: u(:,:), v(:,:)

  dxinv(1) = 1.d0/delta(1)
  dxinv(2) = 1.d0/delta(2)

  allocate(u(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1))
  allocate(v(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1))

  do j = lo(2)-1, hi(2)+1
     do i = lo(1)-1, hi(1)+1
        rhoinv = 1.d0/dat(i,j,1)
        u(i,j) = dat(i,j,2) * rhoinv
        v(i,j) = dat(i,j,3) * rhoinv
     end do
  end do

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        uy = 0.5d0 * (u(i,j+1) - u(i,j-1)) * dxinv(2) 
        vx = 0.5d0 * (v(i+1,j) - v(i-1,j)) * dxinv(1)
        vort(i,j,1) = vx - uy
     end do
  end do

  deallocate(u,v)
  
end subroutine rns_dermagvort

!-----------------------------------------------------------------------

subroutine rns_derdivu(divu,div_l1,div_l2,div_h1,div_h2,nd, &
     &                  dat,dat_l1,dat_l2,dat_h1,dat_h2,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  implicit none

  integer,intent(in):: lo(2), hi(2)
  integer,intent(in):: div_l1,div_l2,div_h1,div_h2,nd
  integer,intent(in):: dat_l1,dat_l2,dat_h1,dat_h2,nc
  integer,intent(in):: domlo(2), domhi(2)
  integer,intent(in):: bc(2,2,nc)
  double precision,intent(in):: delta(2), xlo(2), time, dt
  double precision,intent(out):: divu(div_l1:div_h1,div_l2:div_h2,nd)
  double precision,intent(in )::  dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
  integer,intent(in):: level, grid_no

  integer          :: i,j
  double precision :: ulo,uhi,vlo,vhi, dxinv(2)

  dxinv(1) = 1.d0/delta(1)
  dxinv(2) = 1.d0/delta(2)

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        uhi = dat(i+1,j,2) / dat(i+1,j,1)
        ulo = dat(i-1,j,2) / dat(i-1,j,1)
        vhi = dat(i,j+1,3) / dat(i,j+1,1)
        vlo = dat(i,j-1,3) / dat(i,j-1,1)
        divu(i,j,1) = 0.5d0 * ( (uhi-ulo) *dxinv(1) + (vhi-vlo) *dxinv(2) )
     end do
  end do
  
end subroutine rns_derdivu

!-----------------------------------------------------------------------

subroutine rns_derspec(spec,spec_l1,spec_l2,spec_h1,spec_h2,ns, &
     &                  dat, dat_l1, dat_l2, dat_h1, dat_h2,nd, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  implicit none
  
  integer,intent(in):: lo(2), hi(2)
  integer,intent(in):: spec_l1,spec_l2,spec_h1,spec_h2,ns
  integer,intent(in)::  dat_l1, dat_l2, dat_h1, dat_h2,nd
  integer,intent(in):: domlo(2), domhi(2)
  integer,intent(in):: bc(2,2,nd)
  double precision,intent(in):: delta(2), xlo(2), time, dt
  double precision,intent(out):: spec(spec_l1:spec_h1,spec_l2:spec_h2,ns)
  double precision,intent(in )::  dat( dat_l1: dat_h1, dat_l2: dat_h2,nd)
  integer,intent(in):: level, grid_no
 
  integer i,j
  
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        spec(i,j,1) = dat(i,j,2) / dat(i,j,1)
     end do
  end do

end subroutine rns_derspec

!-----------------------------------------------------------------------

subroutine rns_dermolefrac(spec,spec_l1,spec_l2,spec_h1,spec_h2,ns, &
     &                      dat, dat_l1, dat_l2, dat_h1, dat_h2,nd, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  use meth_params_module, only : NSPEC
  use eos_module, only : eos_YtoX

  implicit none
  
  integer,intent(in):: lo(2), hi(2)
  integer,intent(in):: spec_l1,spec_l2,spec_h1,spec_h2,ns
  integer,intent(in)::  dat_l1, dat_l2, dat_h1, dat_h2,nd
  integer,intent(in):: domlo(2), domhi(2)
  integer,intent(in):: bc(2,2,nd)
  double precision,intent(in):: delta(2), xlo(2), time, dt
  double precision,intent(out):: spec(spec_l1:spec_h1,spec_l2:spec_h2,ns)
  double precision,intent(in )::  dat( dat_l1: dat_h1, dat_l2: dat_h2,nd)
  integer,intent(in):: level, grid_no
 
  integer :: i,j
  double precision :: Yt(NSPEC),Xt(NSPEC), rhoInv

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        rhoInv = 1.d0 / dat(i,j,1)
        Yt = dat(i,j,2:NSPEC+1) * rhoInv
        
        call eos_YtoX(Yt, Xt)
        
        spec(i,j,:) = Xt
     end do
  end do

end subroutine rns_dermolefrac


!-----------------------------------------------------------------------

subroutine rns_dervel(vel,vel_l1,vel_l2,vel_h1,vel_h2,nv, &
     &                dat,dat_l1,dat_l2,dat_h1,dat_h2,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will derive the velocity from the momentum.
  !
  implicit none

  integer,intent(in):: lo(2), hi(2)
  integer,intent(in):: vel_l1,vel_l2,vel_h1,vel_h2,nv
  integer,intent(in):: dat_l1,dat_l2,dat_h1,dat_h2,nc
  integer,intent(in):: domlo(2), domhi(2)
  integer,intent(in):: bc(2,2,nc)
  double precision,intent(in):: delta(2), xlo(2), time, dt
  double precision,intent(out):: vel(vel_l1:vel_h1,vel_l2:vel_h2,nv)
  double precision,intent(in ):: dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
  integer,intent(in):: level, grid_no
 
  integer i,j

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        vel(i,j,1) = dat(i,j,2) / dat(i,j,1)
     end do
  end do
  
end subroutine rns_dervel

