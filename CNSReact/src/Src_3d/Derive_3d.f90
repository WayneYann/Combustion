!-----------------------------------------------------------------------

subroutine cns_derpres(p,p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,ncomp_p, &
     &                 u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use chemistry_module, only : nspecies
  use eos_module, only : eos_get_p
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UFS

  implicit none

  integer,intent(in) :: p_l1,p_l2,p_l3,p_h1,p_h2,p_h3,ncomp_p
  integer,intent(in) :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
  integer,intent(in) :: lo(3), hi(3), domlo(3), domhi(3)
  double precision,intent(out) :: p(p_l1:p_h1,p_l2:p_h2,p_l3:p_h3,ncomp_p)
  double precision,intent(in ) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
  double precision,intent(in) :: dx(3), xlo(3), time, dt
  integer,intent(in) :: bc(3,2,ncomp_u), level, grid_no

  double precision :: Y(nspecies), rhoInv
  integer          :: i,j,k,n
  !
  ! Compute pressure from the EOS
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k,n,Y,rhoInv)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rhoInv = 1.d0/u(i,j,k,URHO)
           do n = 1,nspecies
              Y(n)=u(i,j,k,UFS+n-1)*rhoInv
           enddo

           call eos_get_p(p(i,j,k,1), u(i,j,k,URHO), u(i,j,k,UTEMP), Y)

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  
end subroutine cns_derpres

!-----------------------------------------------------------------------

subroutine cns_derkineng(kineng,ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk, &
     &                      dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will derive kinetic energy = 1/2 rho (u^2 + v^2)
  !
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: ken_l1,ken_l2,ken_l3,ken_h1,ken_h2,ken_h3,nk
  integer,intent(in):: dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
  integer,intent(in):: domlo(3), domhi(3)
  integer,intent(in):: bc(3,2,nc), level, grid_no
  double precision,intent(in):: delta(3), xlo(3), time, dt
  double precision,intent(out):: kineng(ken_l1:ken_h1,ken_l2:ken_h2,ken_l3:ken_h3,nk)
  double precision,intent(in )::    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)

  integer i,j,k

  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           kineng(i,j,k,1) = 0.5d0 / dat(i,j,k,1) * ( dat(i,j,k,2)**2 + &
                &                                     dat(i,j,k,3)**2 + &
                &                                     dat(i,j,k,4)**2 )
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine cns_derkineng


subroutine cns_dersoundspeed(c,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c, &
     &                       u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use chemistry_module, only : nspecies
  use eos_module, only : eos_get_c
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UFS

  implicit none

  integer,intent(in):: c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,ncomp_c
  integer,intent(in):: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
  integer,intent(in):: lo(3), hi(3), domlo(3), domhi(3)
  double precision,intent(out):: c(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,ncomp_c)
  double precision,intent(in ):: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
  double precision,intent(in):: dx(3), xlo(3), time, dt
  integer,intent(in):: bc(3,2,ncomp_u), level, grid_no

  double precision :: Y(nspecies), rhoInv
  integer          :: i,j,k,n

  !
  ! Compute soundspeed from the EOS.
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k,n,Y,rhoInv)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rhoInv = 1.d0/u(i,j,k,URHO)
           do n = 1,nspecies
              Y(n)=u(i,j,k,UFS+n-1)*rhoInv
           enddo

           call eos_get_c(c(i,j,k,1), u(i,j,k,URHO), u(i,j,k,UTEMP), Y)

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  
end subroutine cns_dersoundspeed


subroutine cns_dermachnumber(mach,mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach,&
     &                          u,   u_l1,   u_l2,   u_l3,   u_h1,   u_h2,   u_h3,ncomp_u,   &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use chemistry_module, only : nspecies
  use eos_module, only : eos_get_c
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UTEMP, UFS

  implicit none

  integer,intent(in):: mach_l1,mach_l2,mach_l3,mach_h1,mach_h2,mach_h3,ncomp_mach
  integer,intent(in)::    u_l1,   u_l2,   u_l3,   u_h1,   u_h2,   u_h3,ncomp_u
  integer,intent(in):: lo(3), hi(3), domlo(3), domhi(3)
  double precision,intent(out)::mach(mach_l1:mach_h1,mach_l2:mach_h2,mach_l3:mach_h3,ncomp_mach)
  double precision,intent(in )::   u(   u_l1:   u_h1,   u_l2:   u_h2,   u_l3:   u_h3,ncomp_u)
  double precision,intent(in) :: dx(3), xlo(3), time, dt
  integer,intent(in) :: bc(3,2,ncomp_u), level, grid_no

  double precision :: c, Y(nspecies), rhoInv, ux, uy, uz
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
           do n = 1,nspecies
              Y(n)=u(i,j,k,UFS+n-1)*rhoInv
           enddo

           call eos_get_c(c, u(i,j,k,URHO), u(i,j,k,UTEMP), Y)

           mach(i,j,k,1) = sqrt(ux**2 + uy**2 + uz**2) / c

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine cns_dermachnumber


!-----------------------------------------------------------------------

subroutine cns_dermagvort(vort,  v_l1,  v_l2,  v_l3,  v_h1,  v_h2,  v_h3,nv, & 
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
  double precision :: uy,uz,vx,vz,wx,wy
  double precision,dimension(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1) :: u,v,w
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
           uy = 0.5d0 * (u(i,j+1,k) - u(i,j-1,k)) / delta(2)
           uz = 0.5d0 * (u(i,j,k+1) - u(i,j,k-1)) / delta(3)
           vx = 0.5d0 * (v(i+1,j,k) - v(i-1,j,k)) / delta(1)
           vz = 0.5d0 * (v(i,j,k+1) - v(i,j,k-1)) / delta(3)
           wx = 0.5d0 * (w(i+1,j,k) - w(i-1,j,k)) / delta(1)
           wy = 0.5d0 * (w(i,j+1,k) - w(i,j-1,k)) / delta(2)
           vort(i,j,k,1) = sqrt((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
end subroutine cns_dermagvort


!-----------------------------------------------------------------------

subroutine cns_derdivu(divu,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3,nd, &
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
  double precision :: ulo,uhi,vlo,vhi,wlo,whi

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
           divu(i,j,k,1) = 0.5d0 * ( (uhi-ulo) / delta(1) + &
                (vhi-vlo) / delta(2) + &
                (whi-wlo) / delta(3) )
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine cns_derdivu


!-----------------------------------------------------------------------

subroutine cns_dereint1(e,e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e, &
     &                  u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEDEN 

  implicit none

  integer,intent(in):: e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e
  integer,intent(in):: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
  integer,intent(in):: lo(3), hi(3), domlo(3), domhi(3)
  double precision,intent(out):: e(e_l1:e_h1,e_l2:e_h2,e_l3:e_h3,ncomp_e)
  double precision,intent(in ):: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
  double precision,intent(in):: dx(3), xlo(3), time, dt
  integer,intent(in):: bc(3,2,ncomp_u), level, grid_no

  double precision :: rhoInv,ux,uy,uz
  integer          :: i,j,k
  !
  ! Compute internal energy from (rho E).
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k,rhoInv,ux,uy,uz)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           rhoInv = 1.d0/u(i,j,k,URHO)
           ux = u(i,j,k,UMX)*rhoInv
           uy = u(i,j,k,UMY)*rhoInv
           uz = u(i,j,k,UMZ)*rhoInv
           e(i,j,k,1) = u(i,j,k,UEDEN)*rhoInv-0.5d0*(ux**2+uy**2+uz**2)
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine cns_dereint1


!-----------------------------------------------------------------------

subroutine cns_dereint2(e,e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e, &
     &                  u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u, &
     lo,hi,domlo,domhi,dx,xlo,time,dt,bc,level,grid_no)

  use meth_params_module, only : URHO, UEINT 

  implicit none

  integer,intent(in):: e_l1,e_l2,e_l3,e_h1,e_h2,e_h3,ncomp_e
  integer,intent(in):: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,ncomp_u
  integer,intent(in):: lo(3), hi(3), domlo(3), domhi(3)
  double precision,intent(out):: e(e_l1:e_h1,e_l2:e_h2,e_l3:e_h3,ncomp_e)
  double precision,intent(in ):: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,ncomp_u)
  double precision,intent(in):: dx(3), xlo(3), time, dt
  integer,intent(in):: bc(3,2,ncomp_u), level, grid_no

  integer          :: i,j,k
  !
  ! Compute internal energy from (rho e).
  !
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           e(i,j,k,1) = u(i,j,k,UEINT) / u(i,j,k,URHO)
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO
  
end subroutine cns_dereint2


!-----------------------------------------------------------------------

subroutine cns_derspec(spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,ns, &
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

end subroutine cns_derspec

!-----------------------------------------------------------------------

subroutine cns_dermolefrac(spec,spec_l1,spec_l2,spec_l3,spec_h1,spec_h2,spec_h3,ns, &
     &                      dat, dat_l1, dat_l2, dat_l3, dat_h1, dat_h2, dat_h3,nd, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)

  use chemistry_module, only : nspecies

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
 
  integer :: i,j,k, iwrk
  double precision :: Yt(nspecies),Xt(nspecies), rhoInv, rwrk

  !$OMP PARALLEL DO PRIVATE(i,j,k,rhoInv,Xt,Yt,iwrk,rwrk)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rhoInv = 1.d0 / dat(i,j,k,1)
           Yt = dat(i,j,k,2:nspecies+1) * rhoInv

           call ckytx(Yt, iwrk, rwrk, Xt)

           spec(i,j,k,:) = Xt
        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine cns_dermolefrac


!-----------------------------------------------------------------------

subroutine cns_dervel(vel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
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
  
end subroutine cns_dervel


!-----------------------------------------------------------------------

subroutine cns_dermagvel(magvel,vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv, &
     &                      dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will derive magnitude of velocity.
  !
  implicit none
  
  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: vel_l1,vel_l2,vel_l3,vel_h1,vel_h2,vel_h3,nv
  integer,intent(in):: dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
  integer,intent(in):: domlo(3), domhi(3)
  integer,intent(in):: bc(3,2,nc)
  double precision,intent(in):: delta(3), xlo(3), time, dt
  double precision,intent(out):: magvel(vel_l1:vel_h1,vel_l2:vel_h2,vel_l3:vel_h3,nv)
  double precision,intent(in )::    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
  integer,intent(in):: level, grid_no

  integer i,j,k

  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           magvel(i,j,k,1) = sqrt( dat(i,j,k,2)**2 + &
                &                  dat(i,j,k,3)**2 + & 
                &                  dat(i,j,k,4)**2 ) / dat(i,j,k,1)
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine cns_dermagvel

!-----------------------------------------------------------------------

subroutine cns_dermagmom(magmom,mom_l1,mom_l2,mom_l3,mom_h1,mom_h2,mom_h3,nv, &
     &                      dat,dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc, &
     lo,hi,domlo,domhi,delta,xlo,time,dt,bc,level,grid_no)
  !
  ! This routine will derive magnitude of momentum.
  !
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in):: mom_l1,mom_l2,mom_l3,mom_h1,mom_h2,mom_h3,nv
  integer,intent(in):: dat_l1,dat_l2,dat_l3,dat_h1,dat_h2,dat_h3,nc
  integer,intent(in):: domlo(3), domhi(3)
  integer,intent(in):: bc(3,2,nc)
  double precision,intent(in):: delta(3), xlo(3), time, dt
  double precision,intent(out):: magmom(mom_l1:mom_h1,mom_l2:mom_h2,mom_l3:mom_h3,nv)
  double precision,intent(in )::    dat(dat_l1:dat_h1,dat_l2:dat_h2,dat_l3:dat_h3,nc)
  integer,intent(in):: level, grid_no

  integer i,j,k
  
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           magmom(i,j,k,1) = sqrt( dat(i,j,k,1)**2 + dat(i,j,k,2)**2 + dat(i,j,k,3)**2 )
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  
end subroutine cns_dermagmom

