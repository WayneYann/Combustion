      subroutine normmass(lo, hi, xsID,Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, &
         Y_h3, Ynorm, YNORM_l1, YNORM_l2, YNORM_l3, YNORM_h1, YNORM_h2, &
         YNORM_h3)

      use cdwrk_module
      implicit none

      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer, intent(in) :: Ynorm_l1, Ynorm_l2, Ynorm_l3, Ynorm_h1,  &
            Ynorm_h2, Ynorm_h3
      integer, intent(in) :: xsID
      double precision, intent(in) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(out) ::  Ynorm(Ynorm_l1:Ynorm_h1,  &
           Ynorm_l2:Ynorm_h2, Ynorm_l3:Ynorm_h3,*)

      integer i, j, k, n
      double precision sum

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               sum = zero
               do n=1,Nspec
                  Ynorm(i,j,k,n) =  MAX( Y(i,j,k,n),zero)
                  sum = sum + Ynorm(i,j,k,n)
               end do
               Ynorm(i,j,k,xsID) = Y(i,j,k,xsID)+ one - sum
            end do
         end do
      end do
      end subroutine normmass

      subroutine frratextp(lo,hi,,X,X_l1, X_l2, X_l3, X_h1, X_h2, X_h3,  &
          T,T_l1, T_l2, T_l3, T_h1, T_h2, T_h3,  &
          FwdK,FwdK_l1, FwdK_l2, FwdK_l3, FwdK_h1, FwdK_h2, FwdK_h3,  &
          RevK,RevK_l1, RevK_l2, RevK_l3, RevK_h1, RevK_h2, RevK_h3,  &
          Patm,rxns,Nrxns)

      use cdwrk_module
      use conp_module

      implicit none
      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: X_l1, X_l2, X_l3, X_h1, X_h2, X_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: FwdK_l1, FwdK_l2, FwdK_l3, FwdK_h1, FwdK_h2, FwdK_h3
      integer, intent(in) :: RevK_l1, RevK_l2, RevK_l3, RevK_h1, RevK_h2, RevK_h3
      integer, intent(in) :: Nrxns
      integer rxns(Nrxns)
      double precision, intent(in) ::  X(X_l1:X_h1, X_l2:X_h2,  & 
           X_l3:X_h3,*)
      double precision, intent(in) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(out) :: FwdK(FwdK_l1:FwdK_h1,  & 
           FwdK_l2:FwdK_h2, FwdK_l3:FwdK_h3,*)
      double precision, intent(out) :: RevK(RevK_l1:RevK_h1,  & 
           RevK_l2:RevK_h2, RevK_l3:RevK_h3,*)
      double precision, intent(in) :: Patm
      double precision, parameter:: scale = 1.d6

      double precision Xt(maxspec),FwdKt(maxreac),RevKt(maxreac)
      double precision sum, Yt(maxspec)
      integer i,j,k,n
      double precision P1atm,RU,RUC,Pdyne
      logical, parameter :: JBB_HACK = .true.


      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      Pdyne = Patm * P1atm

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Xt(n) = X(i,j,k,n)
               end do
               if(JBB_HACK)then
                 CALL CKXTY(Xt,IWRK(ckbi),RWRK(ckbr),Yt)
                 sum = zero
                 do n=1,Nspec
                    Yt(n) =MAX( Yt(n),zero)
                    sum = sum+Yt(n)
                 end do
                 if (iN2 .gt. 0) then
                    Yt(iN2) = Yt(iN2)+one-sum
                 endif
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),Xt)
               endif
!  #ifdef MIKE1
!                 CALL CKKFKR(Pdyne,T(i,j,k),Xt,IWRK(ckbi),RWRK(ckbr),FwdKt,RevKt)
!  #else
               call bl_abort('FORT_FRrateXTP: CKKFKR not available')
!  #endif
               do n=1,Nrxns
                  FwdK(i,j,k,n) = FwdKt(rxns(n)+1)*scale
                  RevK(i,j,k,n) = RevKt(rxns(n)+1)*scale
               end do
            end do
         end do
      end do
      end subroutine frratextp

      subroutine htrls(lo,hi,Y,Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, & 
         T,T_l1, T_l2, T_l3, T_h1, T_h2, T_h3,  &
         Q,Q_l1, Q_l2, Q_l3, Q_h1, Q_h2, Q_h3,Patm)
      use cdwrk_module
      use conp_module
      implicit none
      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Q_l1, Q_l2, Q_l3, Q_h1, Q_h2, Q_h3
      double precision, intent(in)  :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(in)  :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(out) :: Q(Q_l1:Q_h1, Q_l2:Q_h2, Q_l3:Q_h3)
      double precision, intent(in)  :: Patm

      double precision :: Zt(maxspec+1),Zdott(maxspec+1)
      integer :: i,j,k,n
      integer :: ndummy
      double precision :: tdummy,P1atm,RU,RUC
      double precision :: RHO, CPB
      double precision, parameter ::  scal = 0.1d0
!
!     NOTE: scal converts result from assumed cgs to MKS (1 erg/s.cm^3 = .1 J/s.m^3)
!

      ndummy = Nspec
      tdummy = 0.0D0
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RWRK(NP) = Patm * P1atm

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Zt(1) = T(i,j,k)
               do n=1,Nspec
                  Zt(n+1) = Y(i,j,k,n)
               end do
               call conpFY(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               CALL CKRHOY(RWRK(NP),Zt(1),Zt(2),IWRK(ckbi),RWRK(ckbr),RHO)
               CALL CKCPBS(Zt(1),Zt(2),IWRK(ckbi),RWRK(ckbr),CPB)
               Q(i,j,k) = Zdott(1) * RHO * CPB * scal
            end do
         end do
      end do
      end subroutine htrls

      subroutine rratey(lo,hi,Y,Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3,  &
            T,T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, &
            Ydot,Ydot_l1, Ydot_l2, Ydot_l3, Ydot_h1, Ydot_h2, Ydot_h3, &
            Patm)
      use cdwrk_module
      use conp_module

      implicit none
      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Ydot_l1, Ydot_l2, Ydot_l3, Ydot_h1, Ydot_h2, Ydot_h3
      double precision, intent(in)  :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(in)  :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(out) :: Ydot(Ydot_l1:Ydot_h1, Ydot_l2:Ydot_h2, Ydot_l3:Ydot_h3,*)
      double precision, intent(in)  :: Patm

      double precision :: Zt(maxspec+1),Zdott(maxspec+1)
      integer :: i,j,k,n
      integer :: ndummy
      double precision :: tdummy,P1atm,RU,RUC

      ndummy = Nspec
      tdummy = 0.
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      RWRK(NP) = Patm * P1atm

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               Zt(1) = T(i,j,k)
               do n=1,Nspec+1
                  Zt(n+1) = Y(i,j,k,n)
               end do
               call conpFY(ndummy,tdummy,Zt,Zdott,RWRK,IWRK)
               do n=1,Nspec
                  Ydot(i,j,k,n) = Zdott(n+1)
               end do
            end do
         end do
      end do
      end subroutine rratey

      subroutine masstomole(lo, hi, Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, &
         X, X_l1, X_l2, X_l3, X_h1, X_h2, X_h3)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer, intent(in) :: X_l1, X_l2, X_l3, X_h1, X_h2, X_h3
      double precision, intent(in)  :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(out) :: X(X_l1:X_h1, X_l2:X_h2, X_l3:X_h3,*)

      double precision ::  Xt(maxspec), Yt(maxspec)
      integer :: i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),Xt)
               do n = 1,Nspec
                  X(i,j,k,n) = Xt(n)
               end do
            end do
         end do
      end do
      end subroutine masstomole
      
      subroutine moletomass(lo, hi, Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, &
         X, X_l1, X_l2, X_l3, X_h1, X_h2, X_h3)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer, intent(in) :: X_l1, X_l2, X_l3, X_h1, X_h2, X_h3
      double precision, intent(in)  :: X(X_l1:X_h1, X_l2:X_h2, X_l3:X_h3,*)
      double precision, intent(out) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)

      double precision ::  Xt(maxspec), Yt(maxspec)
      integer :: i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Xt(n) = X(i,j,k,n)
               end do
               CALL CKXTY(Xt,IWRK(ckbi),RWRK(ckbr),Yt)
               do n = 1,Nspec
                  Y(i,j,k,n) = Yt(n)
               end do
            end do
         end do
      end do
      end subroutine moletomass

      subroutine masstp_to_conc(lo, hi, Patm,Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, & 
         T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, C, C_l1, C_l2, C_l3, C_h1, C_h2, C_h3)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: C_l1, C_l2, C_l3, C_h1, C_h2, C_h3
      double precision, intent(in) :: Patm
      double precision, intent(in) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(in) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(out) :: C(C_l1:C_h1, C_l2:C_h2, C_l3:C_h3,*)

      
      double precision :: Yt(maxspec), Ct(maxspec), RU, RUC, P1ATM, Ptmp
      integer :: i,j,k,n

      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKYTCP(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),Ct)
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)*1.d6
               end do
            end do
         end do
      end do
      end subroutine masstp_to_conc

      subroutine massr_to_conc(lo, hi, Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, &
         T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3,  &
         RHO, RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3, &
         C, C_l1, C_l2, C_l3, C_h1, C_h2, C_h3)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: C_l1, C_l2, C_l3, C_h1, C_h2, C_h3
      integer, intent(in) :: RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3
      double precision, intent(in) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(in) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(out) :: C(C_l1:C_h1, C_l2:C_h2, C_l3:C_h3,*)
      double precision, intent(in) :: RHO(RHO_l1:RHO_h1, RHO_l2:RHO_h2, RHO_l3:RHO_h3)

      double precision :: Yt(maxspec), Ct(maxspec), rhoScl
      integer :: i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               rhoScl = RHO(i,j,k)*1.d-3
               CALL CKYTCR(rhoScl,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),Ct)
               do n = 1,Nspec
                  C(i,j,k,n) = Ct(n)*1.d6
               end do
            end do
         end do
      end do
      end subroutine massr_to_conc

      subroutine conc_to_mole(lo, hi,C, C_l1, C_l2, C_l3, C_h1, C_h2, C_h3, &
         X, X_l1, X_l2, X_l3, X_h1, X_h2, X_h3)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3)
      integer, intent(in) :: hi(3)
      integer, intent(in) :: C_l1, C_l2, C_l3, C_h1, C_h2, C_h3
      integer, intent(in) :: X_l1, X_l2, X_l3, X_h1, X_h2, X_h3
      double precision, intent(in) :: C(C_l1:C_h1, C_l2:C_h2, C_l3:C_h3,*)
      double precision, intent(out) ::  X(X_l1:X_h1, X_l2:X_h2, X_l3:X_h3,*)


      double precision :: Ct(maxspec), Xt(maxspec)
      integer :: i,j,k,n

      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Ct(n) = C(i,j,k,n)*1.d-6
               end do
               CALL CKCTX(Ct,IWRK(ckbi),RWRK(ckbr),Xt)
               do n = 1,Nspec
                  X(i,j,k,n) = Xt(n)
               end do
            end do
         end do
      end do
      end subroutine conc_to_mole

      subroutine molprod(lo, hi, id, Q, Q_l1, Q_l2, Q_l3, Q_h1, Q_h2, Q_h3, &
           C, C_l1, C_l2, C_l3, C_h1, C_h2, C_h3, T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3 )
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3), hi(3), id
      integer, intent(in) :: Q_l1, Q_l2, Q_l3, Q_h1, Q_h2, Q_h3
      integer, intent(in) :: C_l1, C_l2, C_l3, C_h1, C_h2, C_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      double precision, intent(out) :: Q(Q_l1:Q_h1, Q_l2:Q_h2, Q_l3:Q_h3,*)
      double precision, intent(in) :: C(C_l1:C_h1, C_l2:C_h2, C_l3:C_h3,*)
      double precision, intent(in) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)



      double preciison :: Ct(maxspec), Qt(maxreac), Qkt(maxreac), millionth
      integer  :: i,j,k,n

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n = 1,Nspec
                  Ct(n) = C(i,j,k,n)*1.d-6
               end do
               CALL CKQC(T(i,j,k),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
!  #ifdef MIKE
!                 CALL CKCONT(id,Qt,IWRK(ckbi),RWRK(ckbr),Qkt)
!  #else
               call bl_abort('FORT_MOLPROD: CKCONT not available')
!  #endif
               do n = 1,Nreac
                  Q(i,j,k,n) = Qkt(n)*1.d6
               end do
            end do
         end do
      end do
      end subroutine molprod

      subroutine geteltmoles(namenc, namlen, lo, hi,  &
         Celt, Celt_l1, Celt_l2, Celt_l3, Celt_h1, Celt_h2, Celt_h3, &
         C, C_l1, C_l2, C_l3, C_h1,C_h2, C_h3)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: namlen, maxlen
      integer, intent(in) :: namenc(namlen)
      integer, intent(in) :: Celt_l1, Celt_l2, Celt_l3, Celt_h1, Celt_h2, Celt_h3
      integer, intent(in) :: C_l1, C_l2, C_l3, C_h1, C_h2, C_h3
      double precision, intent(inout) :: Celt(Celt_l1:Celt_h1, Celt_l2:Celt_h2, Celt_l3:Celt_h3)
      double precision, intent(in) :: C(C_l1:C_h1, C_l2:C_h2, C_l3:C_h3,*)

      integer :: thenames(maxelts*2)
      logical :: match
      integer :: i, j, k, theidx, n, lout
      integer :: NCF(Nelt,Nspec)
!     Find index of desired element
      CALL CKSYME(thenames,2)
      theidx = -1
      do i=1,Nelt
         match = .true.
         do j=1,namlen               
            if (namenc(j) .NE. thenames((i-1)*2+j)) match = .false.
         enddo
         if (match .eqv. .true.) theidx = i
      end do
      if (theidx.lt.0) then
         call bl_pd_abort()
      endif
!     Get the matrix of elements versus species
      call CKNCF(Nelt,IWRK,RWRK,NCF)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               Celt(i,j,k) = zero
               do n = 1,Nspec
                  Celt(i,j,k) = Celt(i,j,k) +C(i,j,k,n)*NCF(theidx,n)
               end do
            end do
         end do
      end do
      end subroutine geteltmoles

      subroutine conpsolv(lo, hi, &
        Ynew, Ynew_l1, Ynew_l2, Ynew_l3, Ynew_h1 , Ynew_h2, Ynew_h3, &
        Tnew, Tnew_l1, Tnew_l2, Tnew_l3, Tnew_h1, Tnew_h2, Tnew_h3, &
        Yold, Yold_l1, Yold_l2, Yold_l3, Yold_h1, Yold_h2, Yold_h3, &
        Told, Told_l1, Told_l2, Told_l3, Told_h1, Told_h2, Told_h3, & 
        FuncCount, FuncCount_l1, FuncCount_l2, FuncCount_l3, &
        FuncCount_h1, FuncCount_h2, FuncCount_h3,Patm,dt,diag,do_diag)

      use cdwrk_module
      use conp_module
      implicit none

      include "vode.H"

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: Yold_l1, Yold_l2, Yold_l3, Yold_h1, Yold_h2, Yold_h3
      integer, intent(in) :: Told_l1, Told_l2, Told_l3, Told_h1, Told_h2, Told_h3
      integer, intent(in) :: Ynew_l1, Ynew_l2, Ynew_l3, Ynew_h1, Ynew_h2, Ynew_h3
      integer, intent(in) :: Tnew_l1, Tnew_l2, Tnew_l3, Tnew_h1, Tnew_h2, Tnew_h3
      integer, intent(in) :: FuncCount_l1, FuncCount_l2, FuncCount_l3, FuncCount_h1, &
           FuncCount_h2, FuncCount_h3
      integer, intent(in) do_diag
      double precision, intent(inout) :: Yold(Yold_l1:Yold_h1, Yold_l2:Yold_h2, Yold_l3:Yold_h3,*)
      double precision, intent(inout) :: Told(Told_l1:Told_h1, Told_l2:Told_h2, Told_l3:Told_h3)
      double precision, intent(inout) :: Ynew(Ynew_l1:Ynew_h1, Ynew_l2:Ynew_h2, Ynew_l3:Ynew_h3,*)
      double precision, intent(inout) :: Tnew(Tnew_l1:Tnew_h1, Tnew_l2:Tnew_h2, Tnew_l3:Tnew_h3)
      double precision, intent(inout) :: FuncCount(FuncCount_l1:FuncCount_h1,  & 
           FuncCount_l2:FuncCount_h2, FuncCount_l3:FuncCount_h3)
      double precision, intent(inout) :: Patm, dt
      double precision, intent(inout) :: diag(FuncCount_l1:FuncCount_h1, &
           FuncCount_l2:FuncCount_h2, FuncCount_l3:FuncCount_h3,*)

      integer, parameter ::  ITOL = 1, IOPT = 1, ITASK = 1
      integer :: open_vode_failure_file
      double precision, parameter :: RTOL =1.0D-8, ATOLEPS=1.0D-8
      double precision ::  ATOL(maxspec+1)
      double precision, paramter :: spec_scalT=2000.d0, NEWJ_TOL=0.01d0
      external conpFY, conpJY, open_vode_failure_file
      double precision :: TT1, TT2, RU, RUC, P1atm, sum, atoln
      integer :: i, j, k, m, MF, ISTATE, LOUTCK, lout

      integer :: nsubchem, nsub, node, tid, nthrds
      double precision ::  dtloc, weight, TT1save, scale, tspecies(0:maxspec), nwdot(0:maxspec)
      double precision ::  Ct(maxspec),Qt(maxreac),Ytemp(maxspec),Yres(maxspec)

      character*(maxspnml) name
      integer, parameter :: LOUTCK=6

      logical  :: newJ_triggered, bad_soln

#ifdef BL_USE_OMP
      include "omp_lib.h"
      REAL_T,  pointer, save :: dvrwrk(:,:)
      integer, pointer, save :: dviwrk(:,:)
      integer,          save :: dvrlen
      integer,          save :: dvilen
#endif
!      
!     Set IOPT=1 parameter settings for VODE (only really need to be set once).
!
      RWRK(dvbr+4) = 0
      RWRK(dvbr+5) = 0
      RWRK(dvbr+6) = 1.d-19
      IWRK(dvbi+4) = 0
      IWRK(dvbi+5) = max_vode_subcycles
      IWRK(dvbi+6) = 0
      
      if (do_diag.eq.1) nsubchem = nchemdiag
!
!     Set molecular weights and pressure in area accessible by conpF
!
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)

      RWRK(NP) = Patm * P1atm

      if (nstiff .eq. 1) then
         MF = 22  ! finite difference jacobian
      else
         MF = 10
      endif
!      
!     Set up ATOL
!
      if (ITOL.eq.2) then
         ATOL(1) = spec_scalT*ATOLEPS
         do m=1,Nspec
            ATOL(m+1) = ATOLEPS*spec_scalY(m)
         end do
      else
         ATOL(1) = ATOLEPS
      endif

      if (do_diag.eq.1) then
         nsub  = nsubchem
         dtloc = dt/nsubchem
      else
         nsub  = 1
         dtloc = dt
      endif

!$omp parallel copyin(/VHACK/)

#ifdef BL_USE_OMP

!$omp critical
      if (.not.associated(dvrwrk)) then
         dvrlen = 22 + 9*(Nspec+1) + 2*(Nspec+1)**2
         dvilen = 30 + Nspec + 1
         nthrds = omp_get_num_threads()

         allocate(dvrwrk(dvrlen,nthrds),dviwrk(dvilen,nthrds))

         do i = 1,nthrds
            dvrwrk(5,i) = 0
            dvrwrk(6,i) = 0
            dvrwrk(7,i) = ten2minus19
            dviwrk(5,i) = 0
            dviwrk(6,i) = max_vode_subcycles
            dviwrk(7,i) = 0
         end do
      endif
!$omp end critical

#define DVRWK   dvrwrk(1,tid)
#define DVIWK   dviwrk(1,tid)
#define DVIFCNT dviwrk(12,tid)
#define DVRLEN  dvrlen
#define DVILEN  dvilen

#else

#define DVRWK   RWRK(dvbr)
#define DVIWK   IWRK(dvbi)
#define DVIFCNT IWRK(dvbi+11)
#define DVRLEN  dvr
#define DVILEN  dvi

#endif /*BL_USE_OMP*/

!$omp do private(i,j,k,m,TT1,ISTATE,sum,atoln,Ytemp,Yres,newJ_triggered)
!$omp&private(scale,bad_soln,Ct,Qt,node,weight,TT1save,TT2,tid,tspecies)
!$omp&private(lout,name) schedule(dynamic,1)

      do k=lo(3),hi(3)

#ifdef BL_USE_OMP
         tid = omp_get_thread_num() + 1
#endif
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               TT1    = zero
               ISTATE = 1

#ifdef SOLN_IS_1D
            if (i.ne.lo(1) .and. j.ne.lo(2)) then
               do m=1,Nspec
                  Ynew(i,j,k,m) = Ynew(lo(1),lo(2),k,m)
               end do
               Tnew(i,j,k) = Tnew(lo(1),lo(2),k)
            else
#endif
               tspecies(0) = Told(i,j,k)

#ifdef DO_JBB_HACK
#ifdef DO_JBB_HACK_TEMP

               tspecies(0) = MAX(HACK_TEMP_MIN, MIN(HACK_TEMP_MAX, tspecies(0)))
#endif
            sum = zero
            do m=1,Nspec
               Ytemp(m) = Yold(i,j,k,m)
               atoln = ATOLEPS
               if (ITOL.eq.2) atoln = atol(m+1)
               Ytemp(m) = MAX(Yold(i,j,k,m),zero)
               sum = sum+Ytemp(m)
            end do
            if (iN2 .gt. 0) Ytemp(iN2) = Ytemp(iN2)+one-sum
#else
            do m=1,Nspec
               Ytemp(m) = Yold(i,j,k,m)
            end do
#endif
            do m = 1,Nspec
               tspecies(m) = Ytemp(m)
            end do

#ifdef TRIGGER_NEW_J
            newJ_triggered = .FALSE.
            sum = zero
            do m=1,NEQ
               scale = spec_scalT
               if (m.ne.1) scale = spec_scalY(m-1)
               sum = sum + ABS(tspecies(m-1)-YJ_SAVE(m))/scale
            end do
            if (sum .gt. NEWJ_TOL) then
               FIRST = .TRUE.
               newJ_triggered = .TRUE.
            end if
#endif

#ifdef ALWAYS_NEW_J
            FIRST = .TRUE.
#endif
            if (do_diag.eq.1) then
               FuncCount(i,j,k) = 0
               CALL CKYTCP(RWRK(NP),tspecies(0),tspecies(1),IWRK(ckbi),RWRK(ckbr),Ct)
               CALL CKQC(tspecies(0),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
               do m=1,Nreac
                  diag(i,j,k,m) = diag(i,j,k,m)+half*dtloc*Qt(m)*million
               enddo
            endif

            do node = 1,nsub
               if (node.lt.nsub) then
                  weight = one
               else
                  weight = half
               endif

               TT1save = TT1
               TT2     = TT1 + dtloc

#if !defined(BL_USE_DOUBLE) || defined(BL_T3E)
               CALL SVODE
#else 
               CALL DVODE
#endif
     &              (CONPF_FILE, NEQ, tspecies, TT1, TT2, ITOL, RTOL, ATOL,
     &              ITASK, ISTATE, IOPT, DVRWK, DVRLEN, DVIWK, DVILEN,
     &              CONPJ_FILE, MF, RWRK, IWRK)
c
c     If the step was bad, and we reused an old Jacobian, try again from scratch
c               
#if defined(TRIGGER_NEW_J) && defined(DO_JBB_HACK)
            if (newJ_triggered .EQV. .FALSE.) then 
               bad_soln = .FALSE.
               do m=1,Nspec
                  if (tspecies(m) .lt. -1.D-6*spec_scalY(m))
     &                 bad_soln = .TRUE.
               end do
               if (bad_soln .EQV. .TRUE.) then
                  FIRST = .TRUE. 
                  tspecies(0) = Told(i,j,k)
                  do m=1,Nspec
                     tspecies(m) = Ytemp(m)
                  end do
                  TT1    = zero
                  ISTATE = 1
#if !defined(BL_USE_DOUBLE) || defined(BL_T3E)
                  CALL SVODE
#else 
                  CALL DVODE
#endif
     &                 (CONPF_FILE, NEQ, tspecies, TT1, TT2, ITOL, RTOL, ATOL,
     &                 ITASK, ISTATE, IOPT, DVRWK, DVRLEN, DVIWK, DVILEN,
     &                 CONPJ_FILE, MF, RWRK, IWRK)
               end if
            end if
#endif
               TT1 = TT2

               if (do_diag.eq.1) then
                  CALL CKYTCP(RWRK(NP),tspecies(0),tspecies(1),IWRK(ckbi),RWRK(ckbr),Ct)
                  CALL CKQC(tspecies(0),Ct,IWRK(ckbi),RWRK(ckbr),Qt)
                  do m=1,Nreac
                     diag(i,j,k,m) = diag(i,j,k,m)+weight*dtloc*Qt(m)*million
                  enddo
                  FuncCount(i,j,k) = FuncCount(i,j,k) + DVIFCNT
               else
                  FuncCount(i,j,k) = DVIFCNT
               endif

               if (verbose_vode .eq. 1) then
#ifdef BL_USE_OMP
                  write(6,*) '......dvode done:'
                  write(6,*) ' last successful step size = ',dvrwrk(11,tid)
                  write(6,*) '          next step to try = ',dvrwrk(12,tid)
                  write(6,*) '   integrated time reached = ',dvrwrk(13,tid)
                  write(6,*) '      number of time steps = ',dviwrk(11,tid)
                  write(6,*) '              number of fs = ',dviwrk(12,tid)
                  write(6,*) '              number of Js = ',dviwrk(13,tid)
                  write(6,*) '    method order last used = ',dviwrk(14,tid)
                  write(6,*) '   method order to be used = ',dviwrk(15,tid)
                  write(6,*) '            number of LUDs = ',dviwrk(19,tid)
                  write(6,*) ' number of Newton iterations ',dviwrk(20,tid)
                  write(6,*) ' number of Newton failures = ',dviwrk(21,tid)
                  if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
                     call get_spec_name(name,dviwrk(16,tid))
                     write(6,*) '   spec with largest error = ',name
                  end if
#else
                  write(6,*) '......dvode done:'
                  write(6,*) ' last successful step size = ',RWRK(dvbr+10)
                  write(6,*) '          next step to try = ',RWRK(dvbr+11)
                  write(6,*) '   integrated time reached = ',RWRK(dvbr+12)
                  write(6,*) '      number of time steps = ',IWRK(dvbi+10)
                  write(6,*) '              number of fs = ',IWRK(dvbi+11)
                  write(6,*) '              number of Js = ',IWRK(dvbi+12)
                  write(6,*) '    method order last used = ',IWRK(dvbi+13)
                  write(6,*) '   method order to be used = ',IWRK(dvbi+14)
                  write(6,*) '            number of LUDs = ',IWRK(dvbi+18)
                  write(6,*) ' number of Newton iterations ',IWRK(dvbi+19)
                  write(6,*) ' number of Newton failures = ',IWRK(dvbi+20)
                  if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
                     call get_spec_name(name,IWRK(dvbi+15))
                     write(6,*) '   spec with largest error = ',name
                  end if
#endif
               end if

               if (ISTATE .LE. -1) then
                  call CONPF_FILE(NEQ, TT1, tspecies, nwdot, RWRK, IWRK)
                  lout = open_vode_failure_file()                  
                  write(lout,*)
                  write(lout,995) 'VODE Failed at (i,j,k) = (',i,',',j,
     &                 ',',k,'),   Return code = ',ISTATE
                  write(lout,996) 'time(T2,Tl,dt)  ',dt, TT1, dt-TT1
                write(lout,995) 'State ID, old, last, dY/dt, dY/dt*(dt)'
                  write(lout,996) 'T               ',
     &                 Told(i,j,k),tspecies(0),nwdot(0),nwdot(0)*(dt-TT1)
                  do m=1,Nspec
                     call get_spec_name(name,m)
                     write(lout,996) name,Yold(i,j,k,m),
     &                    tspecies(m),nwdot(m),nwdot(m)*(dt-TT1)
                  end do

 995              format(a,4(i4,a))
 996              format(a16,1x,4e30.22)
                  close(lout)
                  call bl_abort('VODE failed, see drop file, exiting')
               end if
             enddo

#ifdef DO_JBB_HACK_POST

               tspecies(0) = MAX(HACK_TEMP_MIN, MIN(HACK_TEMP_MAX, tspecies(0)))
#endif
               Tnew(i,j,k) = tspecies(0)

               do m= 1,Nspec
                  Yres(m) = tspecies(m)
               end do

#if defined(DO_JBB_HACK)
               do m=1,Nspec
                  Ynew(i,j,k,m) = Yold(i,j,k,m)+Yres(m)-Ytemp(m)
               end do
#elif  defined(DO_JBB_HACK_POST)
            sum = zero
            do m=1,Nspec
               Ytemp(m) = Yres(i,j,k,m)
               atoln = ATOLEPS
               if (ITOL.eq.2) atoln = atol(m+1)
               Ytemp(m) = MAX(Yres(i,j,k,m),zero)
               sum = sum+Ytemp(m)
            end do
            if (iN2 .gt. 0) then
               Ytemp(iN2) = Ytemp(iN2)+one-sum
            endif
#else
               do m=1,Nspec
                  Ynew(i,j,k,m) = Yres(m)
               end do
#endif

#ifdef SOLN_IS_1D
               endif
#endif
            end do
         end do
      end do
!$omp end do

!$omp end parallel


#undef DVSRWK
#undef DVSIWK
#undef DVRLEN
#undef DVILEN

      end subroutine conpsolv
 
      subroutine mixavg_rhodiff_temp(lo, hi, & 
        RD, RD_l1, RD_l2, RD_l3, RD_h1, RD_h2, RD_h3, &
        T,T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, & 
        Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, Patm, do_temp, do_VelVisc)

      use cdwrk_module
      implicit none

      integer, intent(in) :: lo(3), hi(3), do_temp, do_VelVisc
      integer, intent(in) :: RD_l1, RD_l2, RD_l3, RD_h1, RD_h2, RD_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(inout) :: RD(RD_l1:RD_h1, RD_l2:RD_h2, RD_l3:RD_h3,*)
      double precision, intent(inout) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(inout) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(in) :: Patm

      integer :: i, j, k, n
      double precision :: RU, RUC, P1ATM, Ptmp, Yt(maxspec), Dt(maxspec)
      double precision :: SCAL, TSCAL, Wavg, RHO, Tt, invmwt(maxspec)
      double precision :: alpha, l1, l2, X(maxspec), CPMS(maxspec)

#ifdef BL_USE_OMP
      include "omp_lib.h"
      REAL_T,  pointer :: egrwrk(:,:)
      integer, pointer :: egiwrk(:,:)
      common /egwork/ egrwrk,egiwrk
      save /egwork/
#endif
      integer egrlen, egilen, tid, nthrds

      double precision, parameter :: SCAL = 0.1d0, TSCAL = 1.0d-5

      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM
      call CKWT(IWRK(ckbi),RWRK(ckbr),invmwt)

      do n=1,Nspec
         invmwt(n) = one / invmwt(n)
      end do

!$omp parallel

#ifdef BL_USE_OMP

!$omp critical
      if (.not.associated(egrwrk)) then
         egrlen = 23 + 14*Nspec + 32*Nspec**2 + 13*eg_nodes
     &        + 30*eg_nodes*Nspec + 5*eg_nodes*Nspec**2
         egilen = Nspec
         nthrds = omp_get_num_threads()

         allocate(egrwrk(egrlen,nthrds),egiwrk(egilen,nthrds))

         do i = 1,nthrds
            do n=1,egrlen
               egrwrk(n,i) = RWRK(egbr+n-1)
            enddo
         enddo
         do i = 1,nthrds
            do n=1,egilen
               egiwrk(n,i) = IWRK(egbi+n-1)
            enddo
         end do
      endif
!$omp end critical

#define EGSRWK egrwrk(1,tid)
#define EGSIWK egiwrk(1,tid)

#else

#define EGSRWK RWRK(egbr)
#define EGSIWK IWRK(egbi)

#endif /*BL_USE_OMP*/

!$omp do private(i,j,k,n,Yt,Tt,alpha,Dt,Wavg)
!$omp&private(CPMS,X,RHO,l1,l2,tid)

      do k=lo(3),hi(3)

#ifdef BL_USE_OMP
         tid = omp_get_thread_num() + 1
#endif
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               Tt = MAX(T(i,j,k),TMIN_TRANS)
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
               CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL EGSPAR(Tt,X,Yt,CPMS,EGSRWK,EGSIWK)
               CALL EGSV1(Ptmp,Tt,Yt,Wavg,EGSRWK,Dt)
               CALL CKRHOY(Ptmp,Tt,Yt,IWRK(ckbi),RWRK(ckbr),RHO)
               do n=1,Nspec
                  RD(i,j,k,n) = RHO * Wavg * invmwt(n) * Dt(n) * SCAL
               end do

               if (do_temp .ne. 0) then
                  alpha = 1
                  CALL EGSL1(alpha,Tt,X,EGSRWK,l1)
                  alpha = -1
                  CALL EGSL1(alpha,Tt,X,EGSRWK,l2)
                  RD(i,j,k,Nspec+1) = half * (l1 + l2) * TSCAL
               endif

               if (do_VelVisc .ne. 0) then
                  CALL EGSE3(Tt,Yt,EGSRWK,RD(i,j,k,Nspec+2))
                  RD(i,j,k,Nspec+2) = RD(i,j,k,Nspec+2) * SCAL
               endif

            end do
         end do
      end do
!$omp end do

!$omp end parallel

#undef EGSRWK
#undef EGSIWK
      end subroutine maxavg_rhodiff_temp

      subroutine mix_shear_visc(lo, hi, &
          eta, eta_l1, eta_l2, eta_l3, eta_h1, eta_h2, eta_h3, &
          T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, &
          Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3)
      use cwdrk_module
      implicit none


      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: eta_l1, eta_l2, eta_l3, eta_h1, eta_h2, eta_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(inout) :: eta(eta_l1:eta_h1, eta_l2:eta_h2, eta_l3:eta_h3)
      double precision, intent(inout) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(inout) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)

      double precision, parameter ::  SCAL = 0.1d0

      
      integer :: i, j, k, n
      double precision :: X(maxspec), Yt(maxspec), CPMS(maxspec), Tt

#ifdef BL_USE_OMP
      include "omp_lib.h"
      REAL_T,  pointer :: egrwrk(:,:)
      integer, pointer :: egiwrk(:,:)
      common /egwork/ egrwrk,egiwrk
      save /egwork/
#endif
      integer egrlen, egilen, tid, nthrds

      parameter(SCAL = tenth)

!$omp parallel

#ifdef BL_USE_OMP

!$omp critical
      if (.not.associated(egrwrk)) then
         egrlen = 23 + 14*Nspec + 32*Nspec**2 + 13*eg_nodes &
            + 30*eg_nodes*Nspec + 5*eg_nodes*Nspec**2
         egilen = Nspec
         nthrds = omp_get_num_threads()

         allocate(egrwrk(egrlen,nthrds),egiwrk(egilen,nthrds))

         do i = 1,nthrds
            do n=1,egrlen
               egrwrk(n,i) = RWRK(egbr+n-1)
            enddo
         enddo
         do i = 1,nthrds
            do n=1,egilen
               egiwrk(n,i) = IWRK(egbi+n-1)
            enddo
         end do
      endif
!$omp end critical

#define EGSRWK egrwrk(1,tid)
#define EGSIWK egiwrk(1,tid)

#else

#define EGSRWK RWRK(egbr)
#define EGSIWK IWRK(egbi)

#endif /*BL_USE_OMP*/
!
! The following computes the mixture averaged shear viscosity using EGLib
! Note that SCAL converts assumed cgs units to MKS (1 g/cm.s = .1 kg/m.s)
!
!$omp do private(i,j,k,n,Yt,Tt,CPMS,X,tid)
      do k=lo(3),hi(3)

#ifdef BL_USE_OMP
         tid = omp_get_thread_num() + 1
#endif
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               Tt = MAX(T(i,j,k),TMIN_TRANS) 
               CALL CKCPMS(Tt,IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL EGSPAR(Tt,X,Yt,CPMS,EGSRWK,EGSIWK)
               CALL EGSE3(Tt,Yt,EGSRWK,eta(i,j,k))
               eta(i,j,k) = eta(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end do

!$omp end parallel

#undef EGSRWK
#undef EGSIWK
      end subroutine mix_shear_visc

      subroutine rhofrompty(lo, hi, &
         RHO, RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3, &
         T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, &
         Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, Patm)
      use cdwrk_module
      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(out) :: RHO(RHO_l1:RHO_h1, RHO_l2:RHO_h2, RHO_l3:RHO_h3)
      double precision, intent(in) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(in) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(in) :: Patm

      
      integer :: i, j, k, n
      double precision :: RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL
!
!     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
!
      double precision, parameter :: SCAL = 1000.0D0)
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKRHOY(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),RHO(i,j,k))
               RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end subroutine rhofrompty
      
      subroutine RHOfromPvTY(lo, hi, & 
         RHO, RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3,  & 
         T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, & 
         Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, &
         P, P_l1, P_l2, P_l3, P_h1, P_h2, P_h3)
      use cdwrk_module
      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      integer, intent(in) :: P_l1, P_l2, P_l3, P_h1, P_h2, P_h3
      double precision, intent(out) ::  RHO(RHO_l1:RHO_h1, RHO_l2:RHO_h2, RHO_l3:RHO_h3)
      double precision, intent(in) ::  T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(in) ::  Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(in) ::  P(P_l1:P_h1, P_l2:P_h2, P_l3:P_h3)


      
      integer :: i, j, k, n
      double precision RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL
!
!     NOTE: SCAL converts result from assumed cgs to MKS (1 g/cm^3 = 1.e3 kg/m^3)
!
      double precision, parameter :: SCAL = 1000.d0
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               Ptmp = P(i,j,k) * P1ATM
               CALL CKRHOY(Ptmp,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),RHO(i,j,k))
               RHO(i,j,k) = RHO(i,j,k) * SCAL
            end do
         end do
      end do
      end subroutine RHOfromPvTY
      
      subroutine PfromRTY(lo, hi, &
          P, P_l1, P_l2, P_l3, P_h1, P_h2, P_h3, & 
          RHO, RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3, &
          T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, &
          Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3)
      use cdwrk_module
      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: P_l1, P_l2, P_l3, P_h1, P_h2, P_h3
      integer, intent(in) :: RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(inout) :: P(P_l1:P_h1, P_l2:P_h2, P_l3:P_h3)
      double precision, intent(inout) :: RHO(RHO_l1:RHO_h1, RHO_l2:RHO_h2, RHO_l3:RHO_h3)
      double precision, intent(inout) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(inout) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      
      integer :: i, j, k, n
      double precision :: Yt(maxspec), RHOt, SCAL, SCAL1
!
!     NOTE: SCAL converts result from assumed cgs to MKS (1 dyne/cm^2 = .1 Pa)
!           SCAL1 converts density (1 kg/m^3 = 1.e-3 g/cm^3)
!
      double precision, parameter :: SCAL = 0.1d0, SCAL1 = 1.e-3

!$omp parallel do private(i,j,k,n,Yt,RHOt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               RHOt = RHO(i,j,k) * SCAL1
               CALL CKPY(RHOt,T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),P(i,j,k))
               P(i,j,k) = P(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end subroutine PfromRTY
      
      subroutine TfromPRY(lo, hi, &
           T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3,& 
           RHO, RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3, &
           Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3, Patm)

      use cdwrk_module

      implicit none
      integer,intent(in) :: lo(3), hi(3)
      integer,intent(in) :: RHO_l1, RHO_l2, RHO_l3, RHO_h1, RHO_h2, RHO_h3
      integer,intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer,intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(inout) :: RHO(RHO_l1:RHO_h1, RHO_l2:RHO_h2, RHO_l3:RHO_h3)
      double precision, intent(inout) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(inout) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      double precision, intent(inout) :: Patm

      
      integer :: i, j, k, n
      double precision :: RU, RUC, P1ATM, Ptmp, Yt(maxspec), SCAL, Wavg, RHOt
!
!     NOTE: SCAL converts density (1 kg/m^3 = 1.e-3 g/cm^3)
!
      double precision, parameter :: SCAL = 1.d-3
      
      CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
      Ptmp = Patm * P1ATM

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),Wavg)
               RHOt = RHO(i,j,k) * SCAL
               T(i,j,k) = Ptmp / (RHOt * RU / Wavg)
            end do
         end do
      end do
      end subroutine TfromPRY
      
      subroutine CPMIXfromTY(lo, hi, &
        CPMIX, CPMIX_l1, CPMIX_l2, CPMIX_l3, CPMIX_h1, CPMIX_h2, CPMIX_h3, &
        T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3 ,&
        Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: CPMIX_l1, CPMIX_l2, CPMIX_l3, CPMIX_h1, CPMIX_h2, CPMIX_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(inout) :: CPMIX(CPMIX_l1:CPMIX_h1, &
                  CPMIX_l2:CPMIX_h2, CPMIX_l3:CPMIX_h3)
      double precision, intent(in) :: T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(in) :: Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)

      
      integer  :: i, j, k, n
      double precision ::  Yt(maxspec)
!
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
!
      double precision, parameter :: SCAL = 1.d-4

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKCPBS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),CPMIX(i,j,k))
               CPMIX(i,j,k) = CPMIX(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do
      end subroutine CPMIXfromTY
      
      subroutine CVMIXfromTY(lo, hi, &
         CVMIX, CVMIX_l1, CVMIX_l2, CVMIX_l3, CVMIX_h1, CVMIX_h2, CVMIX_h3, &
         T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, &
         Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: CVMIX_l1, CVMIX_l2, CVMIX_l3, CVMIX_h1, CVMIX_h2, CVMIX_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(out) :: CVMIX(CVMIX_l1:CVMIX_h1, &
                            CVMIX_l2:CVMIX_h2, CVMIX_l3:CVMIX_h3)
      double precision, intent(in) ::  T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(in) ::  Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)

      
      integer  :: i, j, k, n
      double precision Yt(maxspec)
!
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
!
      double precision, parameter :: SCAL = 1.d-4

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKCVBS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),CVMIX(i,j,k))
               CVMIX(i,j,k) = CVMIX(i,j,k) * SCAL
            end do
         end do
      end do
      end subroutine CVMIXfromTY
      
      subroutine HMIXfromTY(lo, hi,  & 
        HMIX, HMIX_l1, HMIX_l2, HMIX_l3, HMIX_h1, HMIX_h2, HMIX_h3,  & 
        T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, & 
        Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3)
      use cdwrk_module
      implicit none


      integer, intent(in) ::  lo(3), hi(3)
      integer, intent(in) ::  HMIX_l1, HMIX_l2, HMIX_l3, HMIX_h1, HMIX_h2, HMIX_h3
      integer, intent(in) ::  T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) ::  Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(out) ::  HMIX(HMIX_l1:HMIX_h1, &
                    HMIX_l2:HMIX_h2, HMIX_l3:HMIX_h3)
      double precision, intent(in) ::  T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(in) ::  Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)
      
      integer :: i, j, k, n
      double precision ::  Yt(maxspec)
!
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
!
      double preicison, parameter :: SCAL = 1.d-4

!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKHBMS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),HMIX(i,j,k))
               HMIX(i,j,k) = HMIX(i,j,k) * SCAL
            end do
         end do
      end do
!$omp end parallel do

      end subroutine HMIXfromTY
      
      subroutine MWMIXfromY(lo, hi, &
          MWMIX, MWMIX_l1, MWMIX_l2, MWMIX_l3, MWMIX_h1, MWMIX_h2, MWMIX_h3,  &
          Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3)
      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: MWMIX_l1, MWMIX_l2, MWMIX_l3, MWMIX_h1, MWMIX_h2, MWMIX_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(inout) ::  MWMIX(MWMIX_l1:MWMIX_h1, &
                             MWMIX_l2:MWMIX_h2, MWMIX_l3:MWMIX_h3)
      double precision, intent(inout) ::  Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)

      
      integer ::  i, j, k,n
      double precision :: Yt(maxspec)
!
!     Returns mean molecular weight in kg/kmole
!
!$omp parallel do private(i,j,k,n,Yt)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do
               CALL CKMMWY(Yt,IWRK(ckbi),RWRK(ckbr),MWMIX(i,j,k))
            end do
         end do
      end do
!$omp end parallel do

      end subroutine MWMIXfromY
      
      subroutine CPfromT(lo, hi,P, CP_l1, CP_l2, CP_l3, CP_h1, CP_h2, CP_h3, &
            T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3)

      use cdwrk_module
      implicit none

      integer, intent(in) ::  lo(3), hi(3)
      integer, intent(in) ::  CP_l1, CP_l2, CP_l3, CP_h1, CP_h2, CP_h3
      integer, intent(in) ::  T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      double precision, intent(out) ::  CP(CP_l1:CP_h1, CP_l2:CP_h2, CP_l3:CP_h3,*)
      double precision, intent(in)  ::  T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      
      integer :: i, j, k, n
      double precision ::  CPt(maxspec)
!
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g.K = 1.e-4 J/kg.K)
!
       double precision, parameter :: SCAL = 1.d-4
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKCPMS(T(i,j,k),IWRK(ckbi),RWRK(ckbr),CPt)
               do n=1,Nspec
                  CP(i,j,k,n) = CPt(n) * SCAL
               end do
            end do
         end do
      end do
      end subroutine CPfromT
      
      subroutine HfromT(lo, hi, H, H_l1, H_l2, H_l3, H_h1, H_h2, H_h3, &
         T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3)
      use cdwrk_module
      implicit none
      integer,intent(in) ::  lo(3), hi(3)
      integer,intent(in) ::  H_l1, H_l2, H_l3, H_h1, H_h2, H_h3
      integer,intent(in) ::  T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      double precision, intent(out) ::  H(H_l1:H_h1, H_l2:H_h2, H_l3:H_h3,*)
      double precision, intent(in) ::  T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)

      
      integer :: i, j, k, n
      double precision :: Ht(maxspec)
!
!     NOTE: SCAL converts result from assumed cgs to MKS (1 erg/g = 1.e-4 J/kg)
!
      double precision, parameter :: SCAL = 1.d-4

!$omp parallel do private(i,j,k,n,Ht)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               CALL CKHMS(T(i,j,k),IWRK(ckbi),RWRK(ckbr),Ht)
               do n=1,Nspec
                  H(i,j,k,n) = Ht(n) * SCAL
               end do
            end do
         end do
      end do
!$omp end parallel do

      end subroutine HfromT

      integer function TfromHY(lo, hi,  & 
         T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3, & 
         HMIX, HMIX_l1, HMIX_l2, HMIX_l3, HMIX_h1, HMIX_h2, HMIX_h3,  & 
         Y, Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3,errMax, NiterMAX, res)
      use cdwrk_module
      implicit none


      integer, intent(inout) :: NiterMAX
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: HMIX_l1, HMIX_l2, HMIX_l3, HMIX_h1, HMIX_h2, HMIX_h3
      integer, intent(in) :: Y_l1, Y_l2, Y_l3, Y_h1, Y_h2, Y_h3
      double precision, intent(out) ::  T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(in) ::  HMIX(HMIX_l1:HMIX_h1, HMIX_l2:HMIX_h2, HMIX_l3:HMIX_h3)
      double precision, intent(in) ::  Y(Y_l1:Y_h1, Y_l2:Y_h2, Y_l3:Y_h3,*)

      double precision, intent(inout) :: errMAX
      double precision, intent(inout) :: res(0:NiterMAX-1)
      double precision, intent(inout) :: Yt(maxspec)
      integer :: i,j,k,n,Niter,MAXiters

      MAXiters = 0

!$omp parallel do private(i,j,k,n,Yt,res,Niter) reduction(max:MAXiters)
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               do n=1,Nspec
                  Yt(n) = Y(i,j,k,n)
               end do

               call FORT_TfromHYpt(T(i,j,k),HMIX(i,j,k),Yt,errMax,NiterMAX,res,Niter)

               if (Niter .lt. 0) then
                  write(6,*) 'T from h,y solve in FORT_TfromHY failed',Niter
                  call bl_abort(" ")
               end if
               
               if (Niter .gt. MAXiters) MAXiters = Niter

            end do
         end do
      end do
!$omp end parallel do

!     Set max iters taken during this solve, and exit
      TfromHY = MAXiters
      return
      end function TfromHY
!
!     Optically thin radiation model, specified at
!            http://www.ca.sandia.gov/tdf/Workshop/Submodels.html
!     
!     Q(T,species) = 4*sigma*SUM{pi*aP,i} *(T4-Tb4) 
!     
!     sigma=5.669e-08 W/m2K4 is the Steffan-Boltzmann constant, 
!     SUM{ } represents a summation over the species in the radiation calculation, 
!     pi is partial pressure of species i in atm (Xi times local pressure)
!     aP,i is the Planck mean absorption coefficient of species i, 1/[m.atm]
!     T is the local flame temperature (K)
!     Tb is the background temperature (300K or as spec. in expt)
!
!     For H2O and CO2,
!         aP = exp{c0 + c1*ln(T) + c2*{ln(T)}2 + c3*{ln(T)}3 + c4*{ln(T)}4} 
!     
!                            H2O                  CO2
!              c0       0.278713E+03         0.96986E+03
!              c1      -0.153240E+03        -0.58838E+03
!              c2       0.321971E+02         0.13289E+03
!              c3      -0.300870E+01        -0.13182E+02
!              c4       0.104055E+00         0.48396E+00
!     For CH4:
!     
!     aP,ch4 = 6.6334 - 0.0035686*T + 1.6682e-08*T2 + 2.5611e-10*T3 - 2.6558e-14*T4
!
!     For CO:   aP,co = c0+T*(c1 + T*(c2 + T*(c3 + T*c4)))
!
!           T <= 750                 else
!      
!         c0   4.7869              10.09       
!         c1  -0.06953             -0.01183    
!         c2   2.95775e-4          4.7753e-6   
!         c3  -4.25732e-7          -5.87209e-10
!         c4   2.02894e-10         -2.5334e-14 
!      
      
      subroutine OTrad_TDF(lo, hi,  & 
         Qloss, Qloss_l1, Qloss_l2, Qloss_l3, Qloss_h1, Qloss_h2, Qloss_h3, & 
         T, T_l1, T_l2, T_l3, T_h1, T_h2, T_h3   ,&
         X, X_l1, X_l2, X_l3, X_h1, X_h2, X_h3, Patm, T_bg)
      use cdwrk_module
      implicit none
      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) :: Qloss_l1, Qloss_l2, Qloss_l3, Qloss_h1, Qloss_h2, Qloss_h3
      integer, intent(in) :: T_l1, T_l2, T_l3, T_h1, T_h2, T_h3
      integer, intent(in) :: X_l1, X_l2, X_l3, X_h1, X_h2, X_h3
      double precision, intent(inout) ::  Qloss(Qloss_l1:Qloss_h1, &
                       Qloss_l2:Qloss_h2, Qloss_l3:Qloss_h3)
      double precision, intent(inout) ::  T(T_l1:T_h1, T_l2:T_h2, T_l3:T_h3)
      double precision, intent(inout) ::  X(X_l1:X_h1, X_l2:X_h2, X_l3:X_h3,*)

      double precision, intent(in) :: Patm, T_bg
      
      character*(maxspnml) name      
      integer :: n, i, j, k, iH2O, iCO2, iCH4, iCO
      double precision :: aP, c0, c1, c2, c3, c4
      double precision :: T1,T2,T3,T4,lnT1,lnT2,lnT3,lnT4,Tb4

      double precision, parameter :: sigma = 5.669D-08
      
      iH2O = 0
      iCO2 = 0
      iCH4 = 0
      iCO  = 0
      
      do n = 1,Nspec
         call get_spec_name(name, n)
         if (name .EQ. 'H20') iH2O = n
         if (name .EQ. 'CO2') iCO2 = n
         if (name .EQ. 'CH4') iCH4 = n
         if (name .EQ. 'CO')  iCO  = n
      end do
      
      Tb4 = T_bg**4

!$omp parallel do private(i,j,k,T1,T2,T3,T4,lnT1,lnT2,lnT3,lnT4)
!$omp&private(aP,c0,c1,c2,c3,c4)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               
               T1 = T(i,j,k)
               T2 = T1*T1
               T3 = T2*T1
               T4 = T3*T1
               
               if ( (iH2O.gt.0) .or. (iCO2.gt.0) ) then
                  lnT1 = LOG(T1)
                  lnT2 = lnT1*lnT1
                  lnT3 = lnT2*lnT1
                  lnT4 = lnT3*lnT1
               end if
               
               aP = zero
               
               if ((iH2O.gt.0).and.(X(i,j,k,iH2O).gt.zero)) then
                  aP = aP + X(i,j,k,iH2O)*EXP(
     &                 + 0.278713D+03
     &                 - 0.153240D+03*lnT1
     &                 + 0.321971D+02*lnT2
     &                 - 0.300870D+01*lnT3
     &                 + 0.104055D+00*lnT4 )
               end if
               
               if ((iCO2.gt.0).and.(X(i,j,k,iCO2).gt.zero)) then            
                  aP = aP + X(i,j,k,iCO2)*EXP(
     &                 + 0.96986D+03
     &                 - 0.58838D+03*lnT1
     &                 + 0.13289D+03*lnT2
     &                 - 0.13182D+02*lnT3
     &                 + 0.48396D+00*lnT4 )
               end if

               if ((iCH4.gt.0).and.(X(i,j,k,iCH4).gt.zero)) then
                  aP = aP + X(i,j,k,iCH4)*
     &                 ( 6.6334D0
     &                 - 0.0035686D0 *T1
     &                 + 1.6682D-08*T2
     &                 + 2.5611D-10*T3
     &                 - 2.6558D-14*T4 )         
               end if
               
               if ((iCO.gt.0).and.(X(i,j,k,iCO).gt.zero)) then
                  if ( T1 .le. 750.0D0 ) then
                     c0 =  4.7869D0
                     c1 = -0.06953D0
                     c2 =  2.95775D-4
                     c3 = -4.25732D-7
                     c4 =  2.02894D-10
                  else
                     c0 =  10.09D0
                     c1 = -0.01183D0
                     c2 =  4.7753D-6
                     c3 = -5.87209D-10
                     c4 = -2.5334D-14
                  endif
                  aP = aP + X(i,j,k,iCO)*(c0 + c1*T1 + c2*T2 + c3*T3 + c4*T4)
               end if

               Qloss(i,j,k) = four*sigma*Patm*(T4-Tb4)*aP
               
            end do
         end do
      end do
!$omp end parallel do

      end subroutine OTrad_TDF
