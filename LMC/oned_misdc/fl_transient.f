c-----------------------------------------------------------------------
      subroutine sn_transient(nx,scal,const_src,
     $                        lin_src_old,lin_src_new,tstep)

      implicit none

c     Quantities passed in.
      integer nx
      real*8         scal(-1:nx  ,nscal)
      real*8    const_src(0 :nx-1,nscal)
      real*8  lin_src_old(0 :nx-1,nscal)
      real*8  lin_src_new(0 :nx-1,nscal)
      real*8  tstep
      
      include 'nums.fi'
      include 'sndata.fi'
      include 'nkbrn.fi'

c     Local variables
      real*8 rhoxold(nspec), aion(nspec), zion(nspec)
      real*8 rhoxnew(nspec)
      real*8 temperature, dens, rho_hstar
      real*8 src_const(nspec+1)
      real*8 src_old(nspec+1)
      real*8 src_new(nspec+1)
      integer i,n
      integer itfc

      call nb_get_az(aion, zion, nspec)
      
      do i = 0,nx-1

         if (scal(i,FirstSpec) .le. 0.d0) then
           do n = 1, nscal
             scal(i,n) = scal(i,n) + tstep* const_src(i,n) 
     $                     + 0.5d0 * tstep*(lin_src_old(i,n)
     $                                     +lin_src_new(i,n))
           enddo
         else
           temperature = scal(i,Temp)
           dens        = scal(i,Density)
           rho_hstar   = scal(i,RhoH)
           src_const(1) = const_src(i,RhoH)
           src_old(1)   = lin_src_old(i,RhoH)
           src_new(1)   = lin_src_new(i,RhoH)
           do n = 1, nspec
             rhoxold(n) = scal(i,FirstSpec-1+n)
             src_const(1+n) =   const_src(i,FirstSpec-1+n)
             src_old(1+n)   = lin_src_old(i,FirstSpec-1+n)
             src_new(1+n)   = lin_src_new(i,FirstSpec-1+n)
           end do
           call burner(i, tstep, temperature, dens, rho_hstar, nspec, 
     &                 rhoxold, aion, zion, rhoxnew,
     &                 nreac, itfc, 
     $                 src_const, src_old, src_new)

           scal(i,Density) = dens
           scal(i,Temp)    = temperature
           scal(i,RhoH)    = rho_hstar
  
           do n = 1, nspec
              scal(i,FirstSpec-1+n) = rhoxnew(n)
           end do
         endif

      end do
      end

C-----------------------------------------------------------------------

      subroutine nb_hgivenrty(
     &    H, H_l1,H_h1,
     &    R,  R_l1,R_h1,
     &    T, T_l1,T_h1,
     &    Y,    Y_l1,Y_h1,
     &    lo, hi)
      implicit none
      integer H_l1,H_h1
      integer R_l1,R_h1
      integer T_l1,T_h1
      integer Y_l1,Y_h1
      DOUBLE PRECISION H(H_l1:H_h1)
      DOUBLE PRECISION R(R_l1:R_h1)
      DOUBLE PRECISION T(T_l1:T_h1)
      DOUBLE PRECISION Y(Y_l1:Y_h1,*)
      integer lo,hi

      include 'nums.fi'
      include 'nkbrn.fi'
C
      integer i, n
      DOUBLE PRECISION aion(MAXSPEC), zion(MAXSPEC)
      DOUBLE PRECISION dens, temp, xmass(MAXSPEC)
      DOUBLE PRECISION pres, enthalpy, eint
      DOUBLE PRECISION c_v, c_p
      DOUBLE PRECISION ne, eta, pele, dpdt, dpdr, dedt, dedr
      integer ins(MAXREAC), outs(MAXREAC)
      DOUBLE PRECISION qs(MAXREAC)
      DOUBLE PRECISION qreac
C
      call nb_get_az(aion, zion, nspec)
      call nb_get_reac(ins, outs, qs, nreac)
      do i = lo,hi
          do n = 1, nspec
            xmass(n) = Y(i,n)
          end do
          call eos(1, R(i), T(i),
     &             nspec, xmass, aion, zion,
     &             pres, enthalpy, eint, c_v, c_p, ne, eta, pele,
     &             dpdt,dpdr,dedt,dedr)
          qreac = 0.0d0
          do n = 1, nreac
             qreac = qreac + Y(i,ins(n))*qs(n)
          end do
          H(i) = enthalpy + qreac
      end do
      end
C-----------------------------------------------------------------------
      subroutine nb_cpgivenrty(
     &     Cp,   Cp_l1,Cp_h1,
     &     R,    R_l1,R_h1,
     &     T,    T_l1,T_h1,
     &     Y,    Y_l1,Y_h1,
     &     lo, hi)
      implicit none
      integer Cp_l1,Cp_h1
      integer R_l1,R_h1
      integer T_l1,T_h1
      integer Y_l1,Y_h1
      DOUBLE PRECISION Cp(Cp_l1:Cp_h1)
      DOUBLE PRECISION R(R_l1:R_h1)
      DOUBLE PRECISION T(T_l1:T_h1)
      DOUBLE PRECISION Y(Y_l1:Y_h1,*)
      integer lo, hi
      include 'nums.fi'
      include 'nkbrn.fi'
C     
      integer i, j, n
      DOUBLE PRECISION aion(MAXSPEC), zion(MAXSPEC)
      DOUBLE PRECISION dens, temp, xmass(MAXSPEC)
      DOUBLE PRECISION pres, enthalpy, eint
      DOUBLE PRECISION c_v, c_p
      DOUBLE PRECISION ne, eta, pele, dpdt, dpdr, dedt, dedr
C     
      call nb_get_az(aion, zion, nspec)
      do i = lo,hi
            do n = 1, nspec
               xmass(n) = Y(i,n)
            end do
            call eos(1, R(i), T(i),
     &           nspec, xmass, aion, zion,
     &           pres, enthalpy, eint, c_v, Cp(i), ne, eta, pele,
     &           dpdt,dpdr,dedt,dedr)
      end do
      end
C-----------------------------------------------------------------------
      subroutine nb_pgivenrty(
     &     P,   P_l1,P_h1,
     &     R,   R_l1,R_h1,
     &     T,   T_l1,T_h1,
     &     Y,   Y_l1,Y_h1,
     &     lo, hi)
      implicit none
      integer P_l1,P_h1
      integer R_l1,R_h1
      integer T_l1,T_h1
      integer Y_l1,Y_h1
      DOUBLE PRECISION P(P_l1:P_h1)
      DOUBLE PRECISION R(R_l1:R_h1)
      DOUBLE PRECISION T(T_l1:T_h1)
      DOUBLE PRECISION Y(Y_l1:Y_h1,*)
      integer lo, hi
      include 'nums.fi'
      include 'nkbrn.fi'
C     
      integer i, n
      DOUBLE PRECISION aion(MAXSPEC), zion(MAXSPEC)
      DOUBLE PRECISION dens, temp, xmass(MAXSPEC)
      DOUBLE PRECISION pres, enthalpy, eint
      DOUBLE PRECISION c_v, c_p
      DOUBLE PRECISION ne, eta, pele, dpdt, dpdr, dedt, dedr
C     
      call nb_get_az(aion, zion, nspec)
      do i = lo,hi
            do n = 1, nspec
               xmass(n) = Y(i,n)
            end do
            call eos(1, R(i), T(i),
     &           nspec, xmass, aion, zion,
     &           P(i), enthalpy, eint, c_v, c_p, ne, eta, pele,
     &           dpdt,dpdr,dedt,dedr)
      end do
      end
C-----------------------------------------------------------------------
      subroutine nb_tgivenrhy(
     &     T, T_l1,T_h1,
     &     R, R_l1,R_h1,
     &     H, H_l1,H_h1,
     &     Y, Y_l1,Y_h1,
     &     lo, hi, cnt)
      implicit none
      integer T_l1,T_h1
      integer R_l1,R_h1
      integer H_l1,H_h1
      integer Y_l1,Y_h1
      DOUBLE PRECISION T(T_l1:T_h1)
      DOUBLE PRECISION R(R_l1:R_h1)
      DOUBLE PRECISION H(H_l1:H_h1)
      DOUBLE PRECISION Y(Y_l1:Y_h1,*)
      integer lo, hi
      integer cnt
      include 'nums.fi'
      include 'nkbrn.fi'
C     
      integer i, n
      DOUBLE PRECISION aion(MAXSPEC), zion(MAXSPEC)
      DOUBLE PRECISION xmass(MAXSPEC)
      DOUBLE PRECISION pres, enthalpy, eint
      DOUBLE PRECISION c_v, c_p
      DOUBLE PRECISION ne, eta, pele, dpdt, dpdr, dedt, dedr
      DOUBLE PRECISION ttt
      integer itemp
      integer ins(MAXREAC), outs(MAXREAC)
      DOUBLE PRECISION qs(MAXREAC)
      DOUBLE PRECISION qreac
C     
      call nb_get_az(aion, zion, nspec)
      call nb_get_reac(ins, outs, qs, nreac)
      cnt = 0
      do i = lo,hi
            do n = 1, nspec
               xmass(n) = Y(i,n)
            end do
            qreac = 0.0d0
            do n = 1, nreac
               qreac = qreac + Y(i,ins(n))*qs(n)
            end do
            enthalpy = H(i) - qreac
            ttt = T(i)
            itemp = 2
            call eos(itemp, R(i), T(i),
     &           nspec, xmass, aion, zion,
     &           pres, enthalpy, eint, c_v, c_p, ne, eta, pele,
     &           dpdt,dpdr,dedt,dedr)
            cnt = max(cnt,itemp)
      end do
      end
C-----------------------------------------------------------------------
      subroutine nb_rgivenpty(
     &     R,   R_l1,R_h1,
     &     P,   P_l1,P_h1,
     &     T,   T_l1,T_h1,
     &     Y,   Y_l1,Y_h1,
     &     lo, hi, cnt)
      implicit none
      integer R_l1,R_h1
      integer P_l1,P_h1
      integer T_l1,T_h1
      integer Y_l1,Y_h1
      DOUBLE PRECISION R(R_l1:R_h1)
      DOUBLE PRECISION P(P_l1:P_h1)
      DOUBLE PRECISION T(T_l1:T_h1)
      DOUBLE PRECISION Y(Y_l1:Y_h1,*)
      integer lo, hi
      integer cnt
      include 'nums.fi'
      include 'nkbrn.fi'
C     
      integer i, j, n
      DOUBLE PRECISION aion(MAXSPEC), zion(MAXSPEC)
      DOUBLE PRECISION temp, xmass(MAXSPEC)
      DOUBLE PRECISION pres, enthalpy, eint
      DOUBLE PRECISION c_v, c_p
      DOUBLE PRECISION ne, eta, pele, dpdt, dpdr, dedt, dedr
      integer itemp
C     
      call nb_get_az(aion, zion, nspec)
      cnt = 0
      do i = lo,hi
            do n = 1, nspec
               xmass(n) = Y(i,n)
            end do
            itemp = 3
            call eos(itemp, R(i), T(i),
     &           nspec,xmass, aion, zion,
     &           P(i), enthalpy, eint, c_v, c_p, ne, eta, pele,
     &           dpdt,dpdr,dedt,dedr)
            cnt = max(cnt,itemp)
      end do
      end

C-----------------------------------------------------------------------
      subroutine nb_get_az(aion,zion,nspecies)

      DOUBLE PRECISION aion(nspecies), zion(nspecies)

c     aion is species atomic weight
c     zion is species proton number

      if (nspecies .eq. 2) then

        aion(1) = 12
        aion(2) = 24

        zion(1) = 6
        zion(2) = 12

      else if (nspecies .eq. 3) then

        aion(1) = 12
        aion(2) = 24
        aion(3) = 16

        zion(1) = 6
        zion(2) = 12
        zion(3) = 8

      else if (nspecies .eq. 4) then

        aion(1) = 12
        aion(2) = 24
        aion(3) = 16
        aion(4) = 32

        zion(1) = 6
        zion(2) = 12
        zion(3) = 8
        zion(4) = 16

      endif

      end

C-----------------------------------------------------------------------

      subroutine nb_get_reac(ins, outs, qs, nnreac)

      integer nnreac
      integer ins(*), outs(*)
      DOUBLE PRECISION qs(*)

      include 'nums.fi'

      nnreac = nreac

      if (nnreac .ge. 1) then

          ins(1) = 1
         outs(1) = 2
           qs(1) = 5.57d17

      endif

      if (nnreac .eq. 2) then

          ins(2) = 3
         outs(2) = 4
           qs(2) = 4.98764d17

      endif

      end
