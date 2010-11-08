
      block data chemdat
      include 'spec.h'
      data traninit / -1 /
      data tranfile / 'tran.asc.chem-H' /
c      data tranfile / 'tran.asc.CH4-2step' /
      data TMIN_TRANS / 0.d0 /
      data LeEQ1 / 1 /
      data max_vode_subcycles / 15000 /
      data min_vode_timestep / 1.d-19 /
      data Pcgs / -1.d0 /
      data iH2  / -1 /
      data iO2  / -1 /
      data iN2  / -1 /
      data iCH4 / -1 /
      data dvd_debug / 0 /
      end

      subroutine initchem
      implicit none
      include 'spec.h'
      integer lout, n, len, offset,i
      integer kname(maxspnml*maxspec)

      if (traninit.lt.0) then
         open(unit=51,status='old',form='formatted',file=tranfile)
         call CNVTTR(51)
         lout = 6
         call EGINI(eg_nodes, lout, eg_IFLAG, eg_ITLS,
     &        EGRWRK, egr, EGIWRK, egi)

c     Other useful things
         call CKRP(IWRK, RWRK, RU, RUC, P1ATM)
         call CKWT(IWRK, RWRK, mwt)         
         do n=1,Nspec
            invmwt(n) = 1.d0 / mwt(n)
         end do

         call cksyms(kname,maxspnml)
         specNameLen = 0
         do n=1,Nspec
            offset = (n-1)*maxspnml+1
            call convStr(kname(offset),maxspnml,specNames(n),len)
            specNamelen = MAX(specNameLen,len)
            if (specNames(n).eq.'H2')  iH2=n
            if (specNames(n).eq.'O2')  iO2=n
            if (specNames(n).eq.'N2')  iN2=n
            if (specNames(n).eq.'CH4') iCH4=n
         enddo
         traninit = 1
      endif
      end

      subroutine CNVTTR (LINKMC)      
      implicit none
      include 'spec.h'
      character*(maxspnml) CDUMMY
      character*16 CFMT, IFMT, LFMT, RFMT
      PARAMETER
     1     (CFMT='(8A16)', IFMT='(10I12)', LFMT='(L8)',
     2     RFMT='(1P,5E24.16)')
      integer MaxOrder, NOrd, NS, LINKEG, LINKMC, IDUMMY, K, NONS,
     &     NONSNS, N
      parameter (MaxOrder = 4, LINKEG=30)
      double precision RDUMMY, PATM
      logical LDUMMY

      double precision WT(maxspec), POLA(maxspec), DIPO(maxspec),
     &     SIGM(maxspec),
     &     EPSI(maxspec), ZROT(maxspec), COFE(maxspec*MaxOrder),
     &     COFL(maxspec*MaxOrder), COFD(maxspec*maxspec*MaxOrder)
      integer LINA(maxspec)
C-----------------------------------------------------------------------
C     CHEMKIN-III TRANSPORT LINKING FILE (Warning: CK I and II different)
C-----------------------------------------------------------------------
      READ (LINKMC,CFMT) CDUMMY
      READ (LINKMC,CFMT) CDUMMY
      READ (LINKMC,CFMT) CDUMMY
      READ (LINKMC,LFMT) LDUMMY
      READ (LINKMC,IFMT) IDUMMY, IDUMMY, NOrd, NS, IDUMMY


      if (NOrd .GT. MaxOrder) then
c         call bl_abort('Polyfit too large for EGLib!')
      end if
c
c     Make sure the transport & chemistry files are consistent.
c
      CALL CKINDX(IDUMMY,RDUMMY,Nelt,Nspec,Nreac,Nfit)

      if (Nelt.gt.maxelts) then
         print *,'Too many elements, increase maxelts'
         print *, 'Nelt = ',Nelt
         print *, 'maxelts = ',maxelts
      endif
      if (Nspec.gt.maxspec) then
         print *,'Too many species, increase maxspec'
      endif
      if (Nreac.gt.maxreac) then
         print *,'Too many reactions, increase maxreac'
      endif

      if (NS .NE. Nspec) then
         print*, 'transport database thinks Nspec = ', NS
         print*, 'chemistry database thinks Nspec = ', Nspec
c         call bl_abort("chemistry & transport are not consistent")
      endif

      NONS = NOrd * NS
      NONSNS = NONS * NS
      READ (LINKMC,RFMT) PATM
      READ (LINKMC,RFMT) (WT(K),   K=1, NS)
      READ (LINKMC,RFMT) (EPSI(K), K=1, NS) 
      READ (LINKMC,RFMT) (SIGM(K), K=1, NS)
      READ (LINKMC,RFMT) (DIPO(K), K=1, NS)
      READ (LINKMC,RFMT) (POLA(K), K=1, NS)
      READ (LINKMC,RFMT) (ZROT(K), K=1, NS) 
      READ (LINKMC,IFMT) (LINA(K), K=1, NS)
      READ (LINKMC,RFMT) (COFL(N), N=1, NONS)
      READ (LINKMC,RFMT) (COFE(N), N=1, NONS)
      READ (LINKMC,RFMT) (COFD(N), N=1, NONSNS)
      rewind(unit=LLINKMC)
C-----------------------------------------------------------------------
C     WRITE EGLIB TRANSPORT LINKING FILE
C-----------------------------------------------------------------------
      OPEN (UNIT=LINKEG,STATUS='unknown',FORM='UNFORMATTED',
     &     FILE='Linkeg')
      WRITE (LINKEG) NS, NOrd, (WT(K), K=1, NS), (EPSI(K), K=1, NS), 
     &                (SIGM(K), K=1, NS), (DIPO(K), K=1, NS), 
     &                (POLA(K), K=1, NS), (ZROT(K), K=1, NS), 
     &                (LINA(K), K=1, NS), (COFE(N), N=1, NONS), 
     &                (COFL(N), N=1, NONS), (COFD(N), N=1, NONSNS)
      CLOSE(UNIT=LINKEG)
      end

      subroutine convStr(codedString,maxLen,string,strLen)
      implicit none
      integer codedString(*), maxLen, k, strLen
      character*(*) string
      do k=1,maxLen
         string(k:k) = ' '
         if (codedString(k) .ne. ICHAR(' ')) then
            string(k:k) = CHAR(codedString(k))
            strLen = k
         endif
      enddo
      end

      integer function FORT_GETCKSPECNAME(i, coded)
      include 'spec.h'
      integer i
      integer coded(*)
      integer names(maxspec*maxspnml)
      integer j, str_len
      call CKSYMS(names, maxspnml)
      do j = 1, maxspnml
         coded(j) = names(maxspnml*(i-1)+j)
      end do
      str_len = 1
      do j = 1, maxspnml
         if (coded(j).eq.ICHAR(' ')) then
            str_len = j
            exit
         endif 
      end do
      FORT_GETCKSPECNAME = str_len - 1
      end

      integer function get_spec_name(name, j)
      include 'spec.h'
      integer i, j, FORT_GETCKSPECNAME
      integer coded(maxspnml)
      character*(maxspnml) name
      get_spec_name = FORT_GETCKSPECNAME(j, coded)
      do i = 1, maxspnml
         name(i:i) = ' '
      end do
      do i = 1, get_spec_name
         name(i:i) = char(coded(i))
      end do
      end

      subroutine FORT_TfromHYpt(T,Hin,Y,Nspec,errMax,NiterMAX,res,Niter)
      implicit none
      double precision T,Y(*),H,Hin
      double precision TMIN,TMAX,errMAX
      integer Nspec,NiterMAX,Niter,n,NiterDAMP, Discont_NiterMAX
      parameter (TMIN=250.d0, TMAX=5000.d0, Discont_NiterMAX=100)
      double precision  T0,cp,cv,dH,temp,RoverWbar,Wbar,RU,RUC,P1ATM
      double precision res(0:NiterMAX-1),dT, Htarg
      logical out_of_bounds, converged, soln_bad, stalled, discont
      double precision HMIN,cpMIN,HMAX,cpMAX
      integer ihitlo,ihithi,j,IWRK
      double precision old_T, old_H, Tsec, Hsec, RWRK

      out_of_bounds(temp) = (temp.lt.TMIN-1.d0) .or. (temp.gt.TMAX)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = 0.5d0*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      soln_bad = .FALSE.
c     Hin in MKS, convert to CGS
c      Htarg = Hin * 1.d4
      Htarg = Hin
      ihitlo = 0
      ihithi = 0

      CALL CKHBMS(T,Y,IWRK,RWRK,H)

      old_T = T
      old_H = H

      Htarg = Hin

      dH = 2.d0*ABS(H - Htarg)/(1.d0 + ABS(H) + ABS(Htarg))
      res(Niter) = dH
      converged = dH.le.errMAX
      stalled = .false.

      do while ((.not.converged) .and. 
     &     (.not.stalled) .and.
     &     (.not.soln_bad))

         CALL CKCPBS(T,Y,IWRK,RWRK,cp)
         dT = (Htarg - H)/cp
         old_T = T
         if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
            T = TMAX
            ihithi = 1
         else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
            T = TMIN
            ihitlo = 1
         else
            T = T + dT
         end if
         soln_bad = out_of_bounds(T)
         if (soln_bad) then
            Niter = -1
            goto 100
         else
            old_H = H
            CALL CKHBMS(T,Y,IWRK,RWRK,H)
            dH = 2.d0*ABS(H - Htarg)/(1.d0 + ABS(H) + ABS(Htarg))
            res(Niter) = min(dH,abs(dT))
            Niter = Niter + 1
         end if
         converged = (dH.le.errMAX)
         stalled = (ABS(dT).le.errMAX)
         if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
            T = TMIN
            CALL CKHBMS(T,Y,IWRK,RWRK,HMIN)
            CALL CKCPBS(T,Y,IWRK,RWRK,cpMIN)
            T=TMIN+(Htarg-HMIN)/cpMIN
            converged = .true.
         endif
         if ((ihithi.eq.1).and.(H.lt.Htarg)) then
            T = TMAX
            CALL CKHBMS(T,Y,IWRK,RWRK,HMAX)
            CALL CKCPBS(T,Y,IWRK,RWRK,cpMAX)
            T=TMAX+(Htarg-HMAX)/cpMAX
            converged = .true.
         endif

c     If the iterations are failing, perhaps it is because the fits are discontinuous
c     The following implements a modified secant method to hone in on a discontinity in h
c     with T.  Discont_NiterMAX is fairly large because this process can be particularly
c     slow to converge if the Htarg value happens to lay between the discontinuous function
c     values.  
         if (Niter .ge. NiterMAX) then
            do while (.not. stalled)
               dT = - (H - Htarg) * (old_T - T)/(old_H - H)
               Tsec = T + dT
               soln_bad = out_of_bounds(Tsec)
               if (soln_bad) then
                  Niter = -3
                  goto 100
               endif
               CALL CKHBMS(Tsec,Y,IWRK,RWRK,Hsec)
               if ( (Hsec-Htarg)*(Htarg-H) .gt. 0.d0 ) then
                  old_H = H
                  old_T = T
               endif
               H = Hsec
               T = Tsec
               stalled = (ABS(dT).le.errMAX)
               Niter = Niter + 1
               if (Niter.gt.NiterMAX+Discont_NiterMAX) then
                  Niter = -2
                  goto 100
               endif
            enddo
            converged = .true.
         endif
      end do

      return
c
c     Error condition....dump state and bail out
c
 100  continue

      write(6,997) 'T from (H,Y): failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,998) 'previous T = ',old_T
      write(6,998) 'target H = ',Htarg
      write(6,998) 'species mass fracs:'
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual:'
      do n = 0,NiterMAX-1
         write(6,998) '  ',res(n)
      end do

 997  format(a,3(i4,a))
 998  format(a,d21.12)
      end
  
      subroutine dscal (n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if ( n.le.0 .or. incx.le.0 ) return
      if (incx.eq.1) goto 20
c
c     code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        dx(i) = da*dx(i)
      enddo
      return
c
c     code for increment equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dx(i) = da*dx(i)
      enddo
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
      enddo
      end

      double precision function vddot (n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      vddot = 0.0d0
      dtemp = 0.0d0
      if (n.le.0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      enddo
      vddot = dtemp
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,5)
      if ( m .eq. 0 ) goto 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      enddo
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      enddo
   60 vddot = dtemp
      end

      subroutine dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if (n.le.0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
      enddo
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) goto 40
      do i = 1,m
        dy(i) = dx(i)
      enddo
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
      enddo
      end

