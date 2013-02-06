
  subroutine outlet_zlo(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(2):hi(2))

    integer :: i,j,k,n
    double precision :: dxinv(3)
    double precision :: rho, u, v, w, T, pres, Y(nspecies), h(nspecies), rhoE
    double precision :: dpdn, dudn, dvdn, dwdn, drhodn, dYdn(nspecies)
    double precision :: L(5+nspecies), Ltr(5+nspecies), lhs(ncons)
    double precision :: S_p, S_Y(nspecies), d_u, d_v, d_w, d_p
    double precision :: hcal, cpWT, gam1

    double precision, dimension(lo(1):hi(1),lo(2):hi(2)) :: &
         dpdx, dpdy, dudx, dvdy, dwdx, dwdy

    k = lo(3)

    do n=1,3
       dxinv(n) = 1.0d0 / dx(n)
    end do

    call comp_trans_deriv_z(k,lo,hi,ngq,Q,dxinv,dlo,dhi,dpdx,dpdy,dudx,dvdy,dwdx,dwdy)

    !$omp parallel do private(i,j,n,rho,u,v,w,T,pres,Y,h,rhoE) &
    !$omp private(dpdn, dudn,dvdn,dwdn, drhodn, dYdn, L, Ltr, lhs) &
    !$omp private(S_p, S_Y, d_u, d_v, d_w, d_p, hcal, cpWT, gam1)
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          rho  = q  (i,j,k,qrho)
          u    = q  (i,j,k,qu)
          v    = q  (i,j,k,qv)
          w    = q  (i,j,k,qw)
          pres = q  (i,j,k,qpres)
          T    = q  (i,j,k,qtemp)
          Y    = q  (i,j,k,qy1:qy1+nspecies-1)
          h    = q  (i,j,k,qh1:qh1+nspecies-1)
          rhoE = con(i,j,k,iene)

          drhodn     = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qrho))
          dudn       = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qu))
          dvdn       = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qv))
          dwdn       = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qw))
          dpdn       = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qpres))
          do n=1,nspecies
             dYdn(n) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qy1+n-1))
          end do

          ! Simple 1D LODI 
          L(1) = (w-aux(ics,i,j))*0.5d0*(dpdn-rho*aux(ics,i,j)*dwdn)
          L(2) = w*(drhodn-dpdn/aux(ics,i,j)**2)
          L(3) = w*dudn
          L(4) = w*dvdn
          L(5) = sigma*aux(ics,i,j)*(1.d0-Ma2_zlo)/(2.d0*Lzdomain)*(pres-Pinfty)
          L(6:) = w*dYdn

          ! multi-D effects
          Ltr(5) = -0.5d0*(u*dpdx(i,j)+v*dpdy(i,j) &
               + aux(igamma,i,j)*pres*(dudx(i,j)+dvdy(i,j)) &
               + rho*aux(ics,i,j)*(u*dwdx(i,j)+v*dwdy(i,j)))
          L(5) = L(5) + min(0.99d0, (1.0d0-abs(w)/aux(ics,i,j))) * Ltr(5) 

          ! viscous and reaction effects
          d_u = fd(i,j,k,imx) / rho
          d_v = fd(i,j,k,imy) / rho
          d_w = fd(i,j,k,imz) / rho

          S_p = 0.d0
          d_p = fd(i,j,k,iene) - (d_u*con(i,j,k,imx)+d_v*con(i,j,k,imy)+d_w*con(i,j,k,imz))
          cpWT = aux(icp,i,j)*aux(iWbar,i,j)*T
          gam1 = aux(igamma,i,j) - 1.d0
          do n=1,nspecies
             hcal = h(n) - cpWT*inv_mwt(n)
             S_p    = S_p - hcal*aux(iwdot1+n-1,i,j)
             S_Y(n) = aux(iwdot1+n-1,i,j) / rho
             d_p    = d_p - hcal*fd(i,j,k,iry1+n-1)
          end do
          S_p = gam1 * S_p
          d_p = gam1 * d_p 

          ! S_Y seems to cause instabilities
          ! So set it to zero for now
          S_Y = 0.d0
          
          L(5) = L(5) + 0.5d0*(S_p + d_p + rho*aux(ics,i,j)*d_w)

          if (w > 0.d0) then
             L(2) = -S_p/aux(ics,i,j)**2
             L(3) = 0.d0
             L(4) = 0.d0
             L(5) = 0.5d0*S_p
             L(6:) = S_Y      
             
             L(5) = 0.5d0*(S_p + d_p + rho*aux(ics,i,j)*d_w) + Ltr(5) &
                  + outlet_eta*aux(igamma,i,j)*pres*(1.d0-Ma2_zlo)/(2.d0*Lzdomain)*w
          end if
          
          call LtoLHS(3, L, lhs, aux(:,i,j), rho, u, v, w, T, Y, h, rhoE)
          
          rhs(i,j,k,:) = rhs(i,j,k,:) - lhs

       end do
    end do
    !$omp end parallel do

  end subroutine outlet_zlo


!   subroutine inlet_zlo(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,qin,dlo,dhi)
!     integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
!     double precision,intent(in   )::dx(3)
!     double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
!     double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
!     double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
!     double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
!     double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(2):hi(2))
!     double precision,intent(in   )::qin(nqin,lo(1):hi(1),lo(2):hi(2))

!     call bl_error("inlet_zlo not implemented")

!   end subroutine inlet_zlo


  subroutine outlet_zhi(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,dlo,dhi)
    integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
    double precision,intent(in   )::dx(3)
    double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
    double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
    double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(2):hi(2))

    integer :: i,j,k,n
    double precision :: dxinv(3)
    double precision :: rho, u, v, w, T, pres, Y(nspecies), h(nspecies), rhoE
    double precision :: dpdn, dudn, dvdn, dwdn, drhodn, dYdn(nspecies)
    double precision :: L(5+nspecies), Ltr(5+nspecies), lhs(ncons)
    double precision :: S_p, S_Y(nspecies), d_u, d_v, d_w, d_p
    double precision :: hcal, cpWT, gam1

    double precision, dimension(lo(1):hi(1),lo(2):hi(2)) :: &
         dpdx, dpdy, dudx, dvdy, dwdx, dwdy

    k = hi(3)

    do n=1,3
       dxinv(n) = 1.0d0 / dx(n)
    end do

    call comp_trans_deriv_z(k,lo,hi,ngq,Q,dxinv,dlo,dhi,dpdx,dpdy,dudx,dvdy,dwdx,dwdy)

    !$omp parallel do private(i,j,n,rho,u,v,w,T,pres,Y,h,rhoE) &
    !$omp private(dpdn, dudn, dvdn, dwdn, drhodn, dYdn, L, Ltr, lhs) &
    !$omp private(S_p, S_Y, d_u, d_v, d_w, d_p, hcal, cpWT, gam1)
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          rho  = q  (i,j,k,qrho)
          u    = q  (i,j,k,qu)
          v    = q  (i,j,k,qv)
          w    = q  (i,j,k,qw)
          pres = q  (i,j,k,qpres)
          T    = q  (i,j,k,qtemp)
          Y    = q  (i,j,k,qy1:qy1+nspecies-1)
          h    = q  (i,j,k,qh1:qh1+nspecies-1)
          rhoE = con(i,j,k,iene)
          
          drhodn     = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qrho))
          dudn       = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qu))
          dvdn       = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qv))
          dwdn       = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qw))
          dpdn       = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qpres))
          do n=1,nspecies
             dYdn(n) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qy1+n-1))
          end do
          
          ! Simple 1D LODI 
          L(1) = sigma*aux(ics,i,j)*(1.d0-Ma2_zhi)/(2.d0*Lzdomain)*(pres-Pinfty)
          L(2) = w*(drhodn-dpdn/aux(ics,i,j)**2)
          L(3) = w*dudn
          L(4) = w*dvdn
          L(5) = (w+aux(ics,i,j))*0.5d0*(dpdn+rho*aux(ics,i,j)*dwdn)
          L(6:) = w*dYdn
          
          ! multi-D effects
          Ltr(1) = -0.5d0*(u*dpdx(i,j)+v*dpdy(i,j) &
               + aux(igamma,i,j)*pres*(dudx(i,j)+dvdy(i,j)) &
               - rho*aux(ics,i,j)*(u*dwdx(i,j)+v*dwdy(i,j)))
          L(1) = L(1) + min(0.99d0, (1.0d0-abs(w)/aux(ics,i,j))) * Ltr(1)
          
          ! viscous and reaction effects
          d_u = fd(i,j,k,imx) / rho
          d_v = fd(i,j,k,imy) / rho
          d_w = fd(i,j,k,imz) / rho

          S_p = 0.d0
          d_p = fd(i,j,k,iene) - (d_u*con(i,j,k,imx)+d_v*con(i,j,k,imy)+d_w*con(i,j,k,imz))
          cpWT = aux(icp,i,j)*aux(iWbar,i,j)*T
          gam1 = aux(igamma,i,j) - 1.d0
          do n=1,nspecies
             hcal = h(n) - cpWT*inv_mwt(n)
             S_p    = S_p - hcal*aux(iwdot1+n-1,i,j)
             S_Y(n) = aux(iwdot1+n-1,i,j) / rho
             d_p    = d_p - hcal*fd(i,j,k,iry1+n-1)
          end do
          S_p = gam1 * S_p
          d_p = gam1 * d_p 

          ! S_Y seems to cause instabilities
          ! So set it to zero for now
          S_Y = 0.d0
          
          L(1) = L(1) + 0.5d0*(S_p + d_p - rho*aux(ics,i,j)*d_w)
          
          if (w < 0.d0) then
             L(1) = 0.5d0*S_p
             L(2) = -S_p/aux(ics,i,j)**2
             L(3) = 0.d0
             L(4) = 0.d0
             L(6:) = S_Y      
             
             L(1) = 0.5d0*(S_p + d_p - rho*aux(ics,i,j)*d_w) + Ltr(1) &
                  - outlet_eta*aux(igamma,i,j)*pres*(1.d0-Ma2_zhi)/(2.d0*Lzdomain)*w
          end if
          
          call LtoLHS(3, L, lhs, aux(:,i,j), rho, u, v, w, T, Y, h, rhoE)
          
          rhs(i,j,k,:) = rhs(i,j,k,:) - lhs

       end do
    end do
    !$omp end parallel do

  end subroutine outlet_zhi


!   subroutine inlet_zhi(lo,hi,ngq,ngc,dx,Q,con,fd,rhs,aux,qin,dlo,dhi)
!     integer, intent(in) :: lo(3), hi(3), ngq, ngc, dlo(3), dhi(3)
!     double precision,intent(in   )::dx(3)
!     double precision,intent(in   )::Q  (-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
!     double precision,intent(in   )::con(-ngc+lo(1):hi(1)+ngc,-ngc+lo(2):hi(2)+ngc,-ngc+lo(3):hi(3)+ngc,ncons)
!     double precision,intent(in   )::fd (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
!     double precision,intent(inout)::rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)
!     double precision,intent(in   )::aux(naux,lo(1):hi(1),lo(2):hi(2))
!     double precision,intent(in   )::qin(nqin,lo(1):hi(1),lo(2):hi(2))
    
!     call bl_error("inlet_zhi not implemented")

!   end subroutine inlet_zhi


  subroutine comp_trans_deriv_z(k,lo,hi,ngq,Q,dxinv,dlo,dhi,dpdx,dpdy,dudx,dvdy,dwdx,dwdy)
    integer, intent(in) :: k, lo(3), hi(3), ngq,dlo(3),dhi(3)
    double precision, intent(in) :: dxinv(3)
    double precision, intent(in) :: Q(-ngq+lo(1):hi(1)+ngq,-ngq+lo(2):hi(2)+ngq,-ngq+lo(3):hi(3)+ngq,nprim)
    double precision, intent(out):: dpdx(lo(1):hi(1),lo(2):hi(2))
    double precision, intent(out):: dpdy(lo(1):hi(1),lo(2):hi(2))
    double precision, intent(out):: dudx(lo(1):hi(1),lo(2):hi(2))
    double precision, intent(out):: dvdy(lo(1):hi(1),lo(2):hi(2))
    double precision, intent(out):: dwdx(lo(1):hi(1),lo(2):hi(2))
    double precision, intent(out):: dwdy(lo(1):hi(1),lo(2):hi(2))

    integer :: i, j, slo(3), shi(3)

    slo = dlo + stencil_ng
    shi = dhi - stencil_ng

    !$omp parallel private(i,j)

    !d()/dx
    !$omp do
    do j=lo(2),hi(2)

       do i=slo(1),shi(1)
          dudx(i,j) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qpres))
       end do
       
       ! lo-x boundary
       if (dlo(1) .eq. lo(1)) then
          i = lo(1)
          ! use completely right-biased stencil
          dudx(i,j) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qpres))

          i = lo(1)+1
          ! use 3rd-order slightly right-biased stencil
          dudx(i,j) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qpres))

          i = lo(1)+2
          ! use 4th-order stencil
          dudx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qpres))

          i = lo(1)+3
          ! use 6th-order stencil
          dudx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qpres))
       end if

       ! hi-x boundary
       if (dhi(1) .eq. hi(1)) then
          i = hi(1)-3
          ! use 6th-order stencil
          dudx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qpres))

          i = hi(1)-2
          ! use 4th-order stencil
          dudx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qpres))

          i = hi(1)-1
          ! use 3rd-order slightly left-biased stencil
          dudx(i,j) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qpres))

          i = hi(1)
          ! use completely left-biased stencil
          dudx(i,j) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qu))
          dwdx(i,j) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qw))
          dpdx(i,j) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qpres))
       end if

    end do
    !$omp end do nowait

    ! d()/dy
    !$omp do
    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qpres))
       end do
    end do
    !$omp end do nowait

    !$omp master
 
    ! lo-y boundary
    if (dlo(2) .eq. lo(2)) then
       j = lo(2)
       ! use completely right-biased stencil
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qpres))
       end do

       j = lo(2)+1
       ! use 3rd-order slightly right-biased stencil
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qpres))
       end do

       j = lo(2)+2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qpres))
       end do

       j = lo(2)+3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qpres))
       end do
    end if

    ! hi-y boundary
    if (dhi(2) .eq. hi(2)) then
       j = hi(2)-3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qpres))
       end do

       j = hi(2)-2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qpres))
       end do

       j = hi(2)-1
       ! use 3rd-order slightly left-biased stencil
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qpres))
       end do

       j = hi(2)
       ! use completely left-biased stencil
       do i=lo(1),hi(1)
          dvdy(i,j) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qv))
          dwdy(i,j) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qw))
          dpdy(i,j) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qpres))
       end do
    end if

    !$omp end master

    !$omp end parallel
  end subroutine comp_trans_deriv_z
