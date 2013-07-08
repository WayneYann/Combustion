!-----------------------------------------------------------------------

subroutine rns_derpres(p,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,ncomp_p, &
     &                 u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use eos_module, only : eos_get_p
  use meth_params_module, only : URHO, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in) :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,ncomp_p
  integer,intent(in) :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
  integer,intent(in) :: lo(3), hi(3), domlo(3), domhi(3)
  double precision,intent(out) :: p(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,ncomp_p)
  double precision,intent(in ) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
  double precision,intent(in) :: dx(3), xlo(3), time, dt
  integer,intent(in) :: bc(3,2,ncomp_u), level, grid_no

  double precision :: Y(NSPEC), rhoInv
  integer          :: i,j,k,n
  !
  ! Compute pressure from the EOS
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k,n,Y,rhoInv)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rhoInv = 1.d0/u(i,j,k,URHO)
           do n = 1,NSPEC
              Y(n)=u(i,j,k,UFS+n-1)*rhoInv
           enddo

           call eos_get_p(p(i,j,k,1), u(i,j,k,URHO), u(i,j,k,UTEMP), Y)

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  
end subroutine rns_derpres

!-----------------------------------------------------------------------

subroutine rns_dersoundspeed(c,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c, &
     &                       u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use eos_module, only : eos_get_c
  use meth_params_module, only : URHO, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in):: c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c
  integer,intent(in):: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
  integer,intent(in):: lo(3), hi(3), domlo(3), domhi(3)
  double precision,intent(out):: c(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,ncomp_c)
  double precision,intent(in ):: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
  double precision,intent(in):: dx(3), xlo(3), time, dt
  integer,intent(in):: bc(3,2,ncomp_u), level, grid_no

  double precision :: Y(NSPEC), rhoInv
  integer          :: i,j,k,n

  !
  ! Compute soundspeed from the EOS.
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k,n,Y,rhoInv)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rhoInv = 1.d0/u(i,j,k,URHO)
           do n = 1,NSPEC
              Y(n)=u(i,j,k,UFS+n-1)*rhoInv
           enddo

           call eos_get_c(c(i,j,k,1), u(i,j,k,URHO), u(i,j,k,UTEMP), Y)

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  
end subroutine rns_dersoundspeed

!-----------------------------------------------------------------------

subroutine rns_dermachnumber(mach,mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach,&
     &                          u,   u_l1,   u_l2,   u_l3,   u_h1,   u_h2,   u_h3,ncomp_u,   &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use eos_module, only : eos_get_c
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UFS, NSPEC

  implicit none

  integer,intent(in):: mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach
  integer,intent(in)::    u_l1,   u_l2,   u_l3,   u_h1,   u_h2,   u_h3,ncomp_u
  integer,intent(in):: lo(3), hi(3), domlo(3), domhi(3)
  double precision,intent(out)::mach(mach_l1:mach_h1,mach_l2:mach_h2,mach_l3:mach_h3,ncomp_mach)
  double precision,intent(in )::   u(   u_l1:   u_h1,   u_l2:   u_h2,   u_l3:   u_h3,ncomp_u)
  double precision,intent(in) :: dx(3), xlo(3), time, dt
  integer,intent(in) :: bc(3,2,ncomp_u), level, grid_no

  double precision :: c, Y(NSPEC), rhoInv, ux, uy, uz
  integer          :: i,j,k,n

  !
  ! Compute Mach number of the flow.
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k,n,c,Y,rhoInv,ux,uy,uz)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rhoInv = 1.d0/u(i,j,k,URHO)
           ux = u(i,j,k,UMX)*rhoInv
           uy = u(i,j,k,UMY)*rhoInv
           uz = u(i,j,k,UMZ)*rhoInv
           do n = 1,NSPEC
              Y(n)=u(i,j,k,UFS+n-1)*rhoInv
           enddo

           call eos_get_c(c, u(i,j,k,URHO), u(i,j,k,UTEMP), Y)

           mach(i,j,k,1) = sqrt(ux**2 + uy**2 + uz**2) / c

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine rns_dermachnumber

!-----------------------------------------------------------------------

subroutine rns_dermagvort(vort,  v_l1,  v_l2,  v_l3,  v_h1,  v_h2,  v_h3,nv, & 
     &                     dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will calculate vorticity
  !     
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in)::  v_l1,  v_l2,  v_l3,  v_h1,  v_h2,  v_h3,nv
  integer,intent(in)::dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
  integer,intent(in):: domlo(3), domhi(3), level, grid_no
  integer,intent(in):: bc(3,2,nc)
  double precision,intent(in):: delta(3), xlo(3), time, dt
  double precision,intent(out):: vort(  v_l1:  v_h1,  v_l2:  v_h2,  v_l3:  v_h3,nv)
  double precision,intent(in )::  dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)

  integer          :: i,j,k
  double precision :: uy,uz,vx,vz,wx,wy, dxinv(3)
  double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)

  allocate(u(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1))
  allocate(v(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1))
  allocate(w(lo(1)-1:hi(1)+1, lo(2)-1:hi(2)+1, lo(3)-1:hi(3)+1))

  dxinv = 1.d0/delta

  !
  ! Convert momentum to velocity.
  !
  !$OMP PARALLEL PRIVATE(i,j,k,uy,uz,vx,vz,wx,wy)

  !$OMP DO 
  do k = lo(3)-1, hi(3)+1
     do j = lo(2)-1, hi(2)+1
        do i = lo(1)-1, hi(1)+1
           u(i,j,k) = dat(i,j,k,2) / dat(i,j,k,1)
           v(i,j,k) = dat(i,j,k,3) / dat(i,j,k,1)
           w(i,j,k) = dat(i,j,k,4) / dat(i,j,k,1)
        end do
     end do
  end do
  !$OMP END DO
  !
  ! Calculate vorticity.
  !
  !$OMP DO 
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           uy = 0.5d0 * (u(i,j+1,k) - u(i,j-1,k)) * dxinv(2)
           uz = 0.5d0 * (u(i,j,k+1) - u(i,j,k-1)) * dxinv(3)
           vx = 0.5d0 * (v(i+1,j,k) - v(i-1,j,k)) * dxinv(1)
           vz = 0.5d0 * (v(i,j,k+1) - v(i,j,k-1)) * dxinv(3)
           wx = 0.5d0 * (w(i+1,j,k) - w(i-1,j,k)) * dxinv(1)
           wy = 0.5d0 * (w(i,j+1,k) - w(i,j-1,k)) * dxinv(2)
           vort(i,j,k,1) = sqrt((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(u,v,w)

end subroutine rns_dermagvort

!-----------------------------------------------------------------------

subroutine rns_derdivu(divu,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd, &
     &                  dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will divergence of velocity.
  !
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd
  integer,intent(in):: dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
  integer,intent(in):: domlo(3), domhi(3)
  integer,intent(in):: bc(3,2,nc)
  double precision,intent(in):: delta(3), xlo(3), time, dt
  double precision,intent(out):: divu(div_l1:div_h1,div_l2:div_h2,div_l3:div_h3,nd)
  double precision,intent(in )::  dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
  integer,intent(in):: level, grid_no

  integer          :: i,j,k
  double precision :: ulo,uhi,vlo,vhi,wlo,whi, dxinv(3)

  dxinv = 1.d0/delta

  !$OMP PARALLEL DO PRIVATE(i,j,k,ulo,uhi,vlo,vhi,wlo,whi)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           uhi = dat(i+1,j,k,2) / dat(i+1,j,k,1)
           ulo = dat(i-1,j,k,2) / dat(i-1,j,k,1)
           vhi = dat(i,j+1,k,3) / dat(i,j+1,k,1)
           vlo = dat(i,j-1,k,3) / dat(i,j-1,k,1)
           whi = dat(i,j,k+1,4) / dat(i,j,k+1,1)
           wlo = dat(i,j,k-1,4) / dat(i,j,k-1,1)
           divu(i,j,k,1) = 0.5d0 * ( (uhi-ulo) * dxinv(1) + &
                (vhi-vlo) * dxinv(2) + &
                (whi-wlo) * dxinv(3) )
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine rns_derdivu

!-----------------------------------------------------------------------

subroutine rns_derspec(spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,ns, &
     &                  dat, dat_l1, dat_l2, dat_l3, dat_h1, dat_h2, dat_h3,nd, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  implicit none
  
  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,ns
  integer,intent(in)::  dat_l1, dat_l2, dat_l3, dat_h1, dat_h2, dat_h3,nd
  integer,intent(in):: domlo(3), domhi(3)
  integer,intent(in):: bc(3,2,nd)
  double precision,intent(in):: delta(3), xlo(3), time, dt
  double precision,intent(out):: spec(spec_l1:spec_h1,spec_l2:spec_h2,spec_l3:spec_h3,ns)
  double precision,intent(in )::  dat( dat_l1: dat_h1, dat_l2: dat_h2, dat_l3: dat_h3,nd)
  integer,intent(in):: level, grid_no
 
  integer i,j,k
  
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           spec(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine rns_derspec

!-----------------------------------------------------------------------

subroutine rns_dermolefrac(spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,ns, &
     &                      dat, dat_l1, dat_l2, dat_l3, dat_h1, dat_h2, dat_h3,nd, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  use meth_params_module, only : NSPEC
  use eos_module, only : eos_YtoX

  implicit none
  
  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,ns
  integer,intent(in)::  dat_l1, dat_l2, dat_l3, dat_h1, dat_h2, dat_h3,nd
  integer,intent(in):: domlo(3), domhi(3)
  integer,intent(in):: bc(3,2,nd)
  double precision,intent(in):: delta(3), xlo(3), time, dt
  double precision,intent(out):: spec(spec_l1:spec_h1,spec_l2:spec_h2,spec_l3:spec_h3,ns)
  double precision,intent(in )::  dat( dat_l1: dat_h1, dat_l2: dat_h2, dat_l3: dat_h3,nd)
  integer,intent(in):: level, grid_no
 
  integer :: i,j,k
  double precision :: Yt(NSPEC),Xt(NSPEC), rhoInv

  !$OMP PARALLEL DO PRIVATE(i,j,k,rhoInv,Xt,Yt)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rhoInv = 1.d0 / dat(i,j,k,1)
           Yt = dat(i,j,k,2:NSPEC+1) * rhoInv

           call eos_YtoX(Yt, Xt)

           spec(i,j,k,:) = Xt
        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine rns_dermolefrac


!-----------------------------------------------------------------------

subroutine rns_dervel(vel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
     &                dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will derive the velocity from the momentum.
  !
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
  integer,intent(in):: dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
  integer,intent(in):: domlo(3), domhi(3)
  integer,intent(in):: bc(3,2,nc)
  double precision,intent(in):: delta(3), xlo(3), time, dt
  double precision,intent(out):: vel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
  double precision,intent(in ):: dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
  integer,intent(in):: level, grid_no
 
  integer i,j,k

  !$OMP PARALLEL DO PRIVATE(i,j,k) 
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           vel(i,j,k,1) = dat(i,j,k,2) / dat(i,j,k,1)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine rns_dervel

