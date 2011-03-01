c *************************************************************************
c ** TRIDIAG **
c ** Do a tridiagonal solve 
c *************************************************************************

      subroutine tridiag(a,b,c,r,u,gam,n)

      implicit none

      integer n
      real*8 a(n),b(n),c(n),r(n),u(n),gam(n)
      integer j
      real*8 bet 
      if (b(1) .eq. 0) print *,'CANT HAVE B(1) = ZERO'

      bet = b(1)
      u(1) = r(1)/bet

      do j = 2,n
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) then
          print *,'TRIDIAG FAILED '
          stop
        endif
        u(j) = (r(j)-a(j)*u(j-1))/bet
      enddo

      do j = n-1,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      enddo

      end

      subroutine print_soln(step,time,scal,filename,dx)
      include 'spec.h'
      integer step
      real*8 time, scal(maxscal,0:nx+1), dx
      real*8 Peos(1:nx), Y(maxspec)
      character*(*) filename
      character*(100) fstr
      character*(100) fname
      character*(maxspnml) names(maxspec)
      integer n,i,j, get_spec_name, nlen(maxspec)

      do n=1,Nspec
         nlen = get_spec_name(names(n), n)
      enddo
      do i=1,nx
         do n=1,Nspec
            Y(n) = scal(FirstSpec+n-1,i)/scal(Density,i)
         enddo
         call CKPY(scal(Density,i),scal(Temp,i),Y,IWRK,RWRK,Peos(i))
      enddo
      write(fname,'(a,a,I0.6,a)') trim(filename),'_',step,'.dat'
      open(unit=12,file=trim(fname))
      write(fstr,'(a1,i2,a2)') '(',Nspec+2,'a)'
      write(12,fstr) 'VARIABLES=X Rho ',(trim(names(n)),' ',n=1,Nspec),
     &     ' RhoH Temp Peos'
      write(12,'(a,I0.6,a,I0.6,5(a,g20.12))') 'ZONE I=',nx,' T= "STEP=',step,
     &     ' time=',time,'" DATAPACKING=POINT STRANDID=1 SOLUTIONTIME=',
     &     time
      do i=1,nx
         write(12,'(50g20.12)') (i+0.5d0)*dx + problo, scal(1,i),
     &        (scal(1+n,i)/scal(1,i),n=1,Nspec),scal(Nspec+2,i),
     &        scal(Nspec+3,i),(Peos(i)-Pcgs)/Pcgs
      enddo
      close(12)
      end

      subroutine print_update(step,time,scal,Peos,rho,filename,dx)
      include 'spec.h'
      integer step
      real*8 time, scal(maxscal,1:nx), Peos(1:nx), rho(1:nx), dx
      character*(*) filename
      character*(100) fstr
      character*(100) fname
      character*(maxspnml) names(maxspec)
      integer n,i,j, get_spec_name, nlen(maxspec)
      do n=1,Nspec
         nlen = get_spec_name(names(n), n)
      enddo
      write(fname,'(a,a,I0.6,a)') trim(filename),'_',step,'.dat'
      open(unit=12,file=trim(fname))
      write(fstr,'(a1,i2,a2)') '(',Nspec+2,'a)'
      write(12,fstr) 'VARIABLES=X Rho ',(trim(names(n)),' ',n=1,Nspec),
     &     ' RhoH Temp Peos'
      write(12,'(a,I0.6,a,I0.6,5(a,g20.12))') 'ZONE I=',nx,' T= "STEP=',step,
     &     ' time=',time,'" DATAPACKING=POINT STRANDID=1 SOLUTIONTIME=',
     &     time
      do i=1,nx
         write(12,'(50g20.12)') (i+0.5d0)*dx + problo, scal(Density,i),
     &        (scal(FirstSpec+n-1,i)/rho(i),n=1,Nspec),scal(RhoH,i),
     &        scal(Nspec+3,i),(Peos(i)-Pcgs)/Pcgs
      enddo
      close(12)
      end

      subroutine print_cp_prime(S,fname,dx)
      include 'spec.h'
      real*8 S(maxscal,0:nx+1), dx, mass(maxspec)
      real*8 cpip(maxspec), cpim(maxspec), cpbp, cpbm, dcpdt, T, dT, deltaT
      character*(100) fstr
      character*(100) fname
      character*(maxspnml) names(maxspec)
      integer n,i,j, get_spec_name, nlen(maxspec)

      real*8 mindiff, thisdiff
      integer mindiffloc


      do n=1,Nspec
         nlen = get_spec_name(names(n), n)
      enddo
      open(unit=12,file=trim(fname))
      write(fstr,'(a1,i2,a2)') '(',Nspec+2,'a)'
      write(12,fstr) 'VARIABLES=X T Cp ',(trim(names(n)),' ',n=1,Nspec)
      write(12,'(a,I0.6,a)') 'ZONE I=',nx,' DATAPACKING=POINT'


      mindiffloc = -1
      mindiff = 1.d8
      do i=0,nx+1
         thisdiff = ABS(S(Temp,i)-1000.d0)
         if (thisdiff.lt.mindiff) then
            mindiff = thisdiff
            mindiffloc = i
         endif
      enddo

      dT = 1.d0
      do i=1,nx
         do n=1,Nspec
            mass(n) = S(FirstSpec+n-1,i)/S(Density,i)
         enddo
         T = S(Temp,i) + dT
         call CKCPMS(T,IWRK,RWRK,cpip)
         call CKCPBS(T,mass,IWRK,RWRK,cpbp) 
         T = S(Temp,i) - dT
         call CKCPMS(T,IWRK,RWRK,cpim)
         call CKCPBS(T,mass,IWRK,RWRK,cpbm) 
         T = S(Temp,i)
         deltaT = 2*dT

         write(12,'(50g20.12)') (i+0.5d0)*dx + problo, T,
     &        (cpbp-cpbm)/deltaT,((cpip(n)-cpim(n))/deltaT,n=1,Nspec)
      enddo
      close(12)
      end

      subroutine print_cp(S,fname,dx)
      include 'spec.h'
      real*8 S(maxscal,0:nx+1), dx, mass(maxspec)
      real*8 cpi(maxspec), cpb
      character*(100) fstr
      character*(100) fname
      character*(maxspnml) names(maxspec)
      integer n,i,j, get_spec_name, nlen(maxspec)

      real*8 mindiff, thisdiff
      integer mindiffloc


      do n=1,Nspec
         nlen = get_spec_name(names(n), n)
      enddo
      open(unit=12,file=trim(fname))
      write(fstr,'(a1,i2,a2)') '(',Nspec+2,'a)'
      write(12,fstr) 'VARIABLES=X T Cp ',(trim(names(n)),' ',n=1,Nspec)
      write(12,'(a,I0.6,a)') 'ZONE I=',nx,' DATAPACKING=POINT'


      mindiffloc = -1
      mindiff = 1.d8
      do i=0,nx+1
         thisdiff = ABS(S(Temp,i)-1000.d0)
         if (thisdiff.lt.mindiff) then
            mindiff = thisdiff
            mindiffloc = i
         endif
      enddo

      do i=1,nx
         do n=1,Nspec
            mass(n) = S(FirstSpec+n-1,i)/S(Density,i)
         enddo
         call CKCPMS(S(Temp,i),IWRK,RWRK,cpi)
         call CKCPBS(S(Temp,i),mass,IWRK,RWRK,cpb) 
         write(12,'(50g20.12)') (i+0.5d0)*dx + problo, S(Temp,i),
     &        cpb,(cpi(n),n=1,Nspec)
      enddo
      close(12)
      end

