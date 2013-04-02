
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"
#include "ArrayLim.H"
#include "SLABSTAT_NS_F.H"
#include "SLABSTAT_HT_F.H"

#define SDIM 3

c
c ::: -----------------------------------------------------------
c ::: This is a general routine to calculate the basic running statistics
c ::: based on the velocity, density, tracer, pressure, temperature, 
c ::: enthalpy, RhoRT and the species mass fractions.  The tracer
c ::: passed in to this routine is not a conserved quantity.  The data
c ::: saved by this routine are sufficient to calculate a core set of
c ::: statistics on these fields using both Reynolds and Favre averages.
c ::: The state should be passed in to this routine in the order,
c :::   Rho, U, V, W, Tr, P, T, RhoH, RhoRT, RhoY(...)
c ::: and the saved statistics are those for FORT_NS_BASICSTATS_NCTRAC 
c ::: with these tacked on at the end,
c :::   <T> <Rho T> <T^2> <Rho T^2> 
c :::   <U T> <V T> <W T> <Rho U T> <Rho V T> <Rho W T>
c :::   <Rho H> <RhoRT> <Rho H^2> <RhoRT^2> <Rho U H> <Rho V H> <Rho W H>
c :::   <Rho Y1> <Rho Y1^2> <Rho U Y1> <Rho V Y1> <Rho W Y1>
c :::   <Rho Y2> <Rho Y2^2> <Rho U Y2> <Rho V Y2> <Rho W Y2>
c :::   <Rho Y3> <Rho Y3^2> <Rho U Y3> <Rho V Y3> <Rho W Y3>
c :::   <Rho Y4> <Rho Y4^2> <Rho U Y4> <Rho V Y4> <Rho W Y4>
c :::      .                    .                      .
c :::      .                    .                      .
c :::      .                    .                      .
c :::      .                    .                      .
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: runStats      <=  Array to contain the running statistics
c ::: DIMS(runStats) => Dimensions of the running statistics array
c ::: nStats         => Number of components in the statistics array
c ::: state          =) Array containing the state data
c ::: DIMS(state)    => Dimensions of the state data array
c ::: nStateComp     => Number of components in the state array
c ::: dt             => Time step
c ::: delta          => Grid spacing
c ::: -----------------------------------------------------------
c
      subroutine FORT_HT_BASICSTATS_NCTRAC(state, DIMS(state), nStateComp,
     $                                     runStats, DIMS(runStats), nStats,
     $                                     dt, delta)

      implicit none

c
c     ::::: passed variables
c
      integer nStats, nStateComp
      integer DIMDEC(state)
      integer DIMDEC(runStats)
      REAL_T dt
      REAL_T delta(SDIM)
      REAL_T runStats(DIMV(runStats),nStats)
      REAL_T state(DIMV(state),nStateComp)

c
c     ::::: local variables
c
      integer i, j, k, n, Rho, xVel, yVel, zVel, Trac, Pres, Temp, RhoH, 
     $        RhoRT, RhoYlo, RhoYhi, nStatsExpect, nGhostExpect, nSpecies,
     $        num, sp, nNSstats
      integer loStats(SDIM), hiStats(SDIM)
      integer loState(SDIM), hiState(SDIM)
      REAL_T  enthalpy, massfrac

      PARAMETER (nGhostExpect = 0)

c
c     ===================================
c     ::: Set the Values of LO and HI :::
c     ===================================
c
      call SET_LOHI(DIMS(runStats), loStats, hiStats)
      call SET_LOHI(DIMS(state), loState, hiState)

      nNSstats = 35
      nSpecies = nStateComp - 9
      nStatsExpect = nNSstats + 17 + 5 * nSpecies

#ifndef NDEBUG
      if (nStats .NE. nStatsExpect) then
        write(*,1000) nStatsExpect, nStats, nStateComp, nSpecies
 1000   format('Incorrect number of statistic and/or state components',
     $        /'FORT_HT_BASICSTATS_NCTRAC.  Should have nStats = ',I5,
     $        /'   nStats = ', I5, 5x, 'nStateComp = ', I5, 'nSpecies = ', I5)
        call bl_abort(" ")
      endif

      do n = 1, SDIM
        if (loState(n) .GT. loStats(n) - nGhostExpect .OR.
     $      hiState(n) .LT. hiStats(n) + nGhostExpect) then
          write(*,1010) n, nGhostExpect, loStats, hiStats, loState, hiState
 1010     format('Incorrect number of ghost cells in the state date in',
     $          /'FORT_HT_BASICSTATS_NCTRAC.',
     $          /'   Direction = ', I2, 5x, 'nGhostExpect = ', I2,
     $          /'   loStats = ', SDIM (I2,1x), 5x, 'hiStats = ', SDIM (I2,1x),
     $          /'   loState = ', SDIM (I2,1x), 5x, 'hiState = ', SDIM (I2,1x))
          call bl_abort(" ")
        endif
      enddo
#endif

      num = 0


c
c     =========================
c     ::: Set State Indices :::
c     =========================
c
      Rho  = 1
      xVel = 2
      yVel = 3
      zVel = 4
      Trac = 5
      Pres = 6
      Temp = 7
      RhoH = 8
      RhoRT = 9
      RhoYlo = 10
      RhoYhi = RhoYlo + nSpecies - 1

c
c     ====================================
c     ::: Calculate Running Statistics :::
c     ====================================
c
c     === Call NS Statistics Routine ===
c
      CALL FORT_NS_BASICSTATS_NCTRAC(state, DIMS(state), Pres,
     $                               runStats, DIMS(runStats), nNSstats,
     $                               dt, delta)

c
c     === Add HT Specific Statistics ===
c
      do k = loStats(3), hiStats(3)
        do j = loStats(2), hiStats(2)
          do i = loStats(1), hiStats(1)
            num = nNSstats + 1

            do n = 1, 2
              runStats(i,j,k,num) = runStats(i,j,k,num) 
     $                                              + dt * state(i,j,k,Temp)**n
              runStats(i,j,k,num+1) = runStats(i,j,k,num+1) 
     $                           + dt * state(i,j,k,Rho) * state(i,j,k,Temp)**n
              num = num + 2
            enddo

            runStats(i,j,k,num) = runStats(i,j,k,num) 
     $                         + dt * state(i,j,k,xVel) * state(i,j,k,Temp)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) 
     $                         + dt * state(i,j,k,yVel) * state(i,j,k,Temp)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) 
     $                         + dt * state(i,j,k,zVel) * state(i,j,k,Temp)
            runStats(i,j,k,num+3) = runStats(i,j,k,num+3) 
     $                         + dt * state(i,j,k,Rho) * state(i,j,k,xVel) 
     $                                                 * state(i,j,k,Temp)
            runStats(i,j,k,num+4) = runStats(i,j,k,num+4) 
     $                         + dt * state(i,j,k,Rho) * state(i,j,k,yVel) 
     $                                                 * state(i,j,k,Temp)
            runStats(i,j,k,num+5) = runStats(i,j,k,num+5) 
     $                         + dt * state(i,j,k,Rho) * state(i,j,k,zVel) 
     $                                                 * state(i,j,k,Temp)
            num = num + 6

            enthalpy = state(i,j,k,RhoH) / state(i,j,k,Rho)
            do n = 1, 2
              runStats(i,j,k,num) = runStats(i,j,k,num) 
     $                                    + dt * state(i,j,k,Rho) * enthalpy**n
              runStats(i,j,k,num+1) = runStats(i,j,k,num+1) 
     $                                             + dt * state(i,j,k,RhoRT)**n
              num = num + 2
            enddo

            runStats(i,j,k,num) = runStats(i,j,k,num) 
     $                         + dt * state(i,j,k,RhoH) * state(i,j,k,xVel) 
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1) 
     $                         + dt * state(i,j,k,RhoH) * state(i,j,k,yVel) 
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2) 
     $                         + dt * state(i,j,k,RhoH) * state(i,j,k,zVel) 
            num = num + 3

            do sp = RhoYlo, RHoYhi
              massfrac = state(i,j,k,sp) / state(i,j,k,Rho)
              runStats(i,j,k,num) = runStats(i,j,k,num)
     $                                    + dt * state(i,j,k,sp)
              runStats(i,j,k,num+1) = runStats(i,j,k,num+1)
     $                                    + dt * state(i,j,k,sp) * massfrac
              runStats(i,j,k,num+2) = runStats(i,j,k,num+2)
     $                               + dt * state(i,j,k,sp) * state(i,j,k,xVel) 
              runStats(i,j,k,num+3) = runStats(i,j,k,num+3)
     $                               + dt * state(i,j,k,sp) * state(i,j,k,yVel) 
              runStats(i,j,k,num+4) = runStats(i,j,k,num+4)
     $                               + dt * state(i,j,k,sp) * state(i,j,k,zVel) 
              num = num + 5
            enddo
          enddo
        enddo
      enddo


      if (num-1 .NE. nStatsExpect) then
        write(*,1020) num-1, nStatsExpect
 1020   format('ERROR: Num-1 = ',I5,'not equal to nStatsExpect = ',I5)
        call bl_abort(" ")
      endif

c
c
      return
      end



c
c ::: -----------------------------------------------------------
c ::: This is a general routine to calculate the basic running statistics
c ::: based on the velocity, density, tracer, pressure, temperature,
c ::: enthalpy, RhoRT and the species mass fractions.  The tracer
c ::: passed in to this routine is a conserved quantity.  The data
c ::: saved by this routine are sufficient to calculate a core set of
c ::: statistics on these fields using both Reynolds and Favre averages.
c ::: The state should be passed in to this routine in the order,
c :::   Rho, U, V, W, Tr, P, T, RhoH, RhoRT, RhoY(...)
c ::: and the saved statistics are those for FORT_NS_BASICSTATS_NCTRAC
c ::: with these tacked on at the end,
c :::   <T> <Rho T> <T^2> <Rho T^2>
c :::   <U T> <V T> <W T> <Rho U T> <Rho V T> <Rho W T>
c :::   <Rho H> <RhoRT> <Rho H^2> <RhoRT^2> <Rho U H> <Rho V H> <Rho W H>
c :::   <Rho Y1> <Rho Y1^2> <Rho U Y1> <Rho V Y1> <Rho W Y1>
c :::   <Rho Y2> <Rho Y2^2> <Rho U Y2> <Rho V Y2> <Rho W Y2>
c :::   <Rho Y3> <Rho Y3^2> <Rho U Y3> <Rho V Y3> <Rho W Y3>
c :::   <Rho Y4> <Rho Y4^2> <Rho U Y4> <Rho V Y4> <Rho W Y4>
c :::      .                    .                      .
c :::      .                    .                      .
c :::      .                    .                      .
c :::      .                    .                      .
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: runStats      <=  Array to contain the running statistics
c ::: DIMS(runStats) => Dimensions of the running statistics array
c ::: nStats         => Number of components in the statistics array
c ::: state          =) Array containing the state data
c ::: DIMS(state)    => Dimensions of the state data array
c ::: nStateComp     => Number of components in the state array
c ::: dt             => Time step
c ::: delta          => Grid spacing
c ::: -----------------------------------------------------------
c
      subroutine FORT_HT_BASICSTATS_CTRAC(state, DIMS(state), nStateComp,
     $                                    runStats, DIMS(runStats), nStats,
     $                                    dt, delta)

      implicit none

c
c     ::::: passed variables
c
      integer nStats, nStateComp
      integer DIMDEC(state)
      integer DIMDEC(runStats)
      REAL_T dt
      REAL_T delta(SDIM)
      REAL_T runStats(DIMV(runStats),nStats)
      REAL_T state(DIMV(state),nStateComp)

c
c     ::::: local variables
c
      integer i, j, k, n, Rho, xVel, yVel, zVel, RhoTrac, Pres, Temp, RhoH,
     $        RhoRT, RhoYlo, RhoYhi, nStatsExpect, nGhostExpect, nSpecies,
     $        num, sp, nNSstats
      integer loStats(SDIM), hiStats(SDIM)
      integer loState(SDIM), hiState(SDIM)
      REAL_T  enthalpy, massfrac

      PARAMETER (nGhostExpect = 0)

c
c     ===================================
c     ::: Set the Values of LO and HI :::
c     ===================================
c
      call SET_LOHI(DIMS(runStats), loStats, hiStats)
      call SET_LOHI(DIMS(state), loState, hiState)

      nNSstats = 35
      nSpecies = nStateComp - 9
      nStatsExpect = nNSstats + 17 + 5 * nSpecies
      num = 0

#ifndef NDEBUG
      if (nStats .NE. nStatsExpect) then
        write(*,1000) nStatsExpect, nStats, nStateComp, nSpecies
 1000   format('Incorrect number of statistic and/or state components',
     $        /'FORT_HT_BASICSTATS_CTRAC.  Should have nStats = ',I5,
     $        /'   nStats = ', I5, 5x, 'nStateComp = ', I5, 'nSpecies = ', I5)
        call bl_abort(" ")
      endif

      do n = 1, SDIM
        if (loState(n) .GT. loStats(n) - nGhostExpect .OR.
     $      hiState(n) .LT. hiStats(n) + nGhostExpect) then
          write(*,1010) n, nGhostExpect, loStats, hiStats, loState, hiState
 1010     format('Incorrect number of ghost cells in the state date in',
     $          /'FORT_HT_BASICSTATS_CTRAC.',
     $          /'   Direction = ', I2, 5x, 'nGhostExpect = ', I2,
     $          /'   loStats = ', SDIM (I2,1x), 5x, 'hiStats = ', SDIM (I2,1x),
     $          /'   loState = ', SDIM (I2,1x), 5x, 'hiState = ', SDIM (I2,1x))
          call bl_abort(" ")
        endif
      enddo
#endif


c
c     =========================
c     ::: Set State Indices :::
c     =========================
c
      Rho  = 1
      xVel = 2
      yVel = 3
      zVel = 4
      RhoTrac = 5
      Pres = 6
      Temp = 7
      RhoH = 8
      RhoRT = 9
      RhoYlo = 10
      RhoYhi = RhoYlo + nSpecies - 1


c
c     ====================================
c     ::: Calculate Running Statistics :::
c     ====================================
c
c     === Call NS Statistics Routine ===
c
      CALL FORT_NS_BASICSTATS_CTRAC(state, DIMS(state), Pres,
     $                              runStats, DIMS(runStats), nNSstats,
     $                              dt, delta)

c
c     === Add HT Specific Statistics ===
c
      do k = loStats(3), hiStats(3)
        do j = loStats(2), hiStats(2)
          do i = loStats(1), hiStats(1)
            num = nNSstats + 1

            do n = 1, 2
              runStats(i,j,k,num) = runStats(i,j,k,num)
     $                                              + dt * state(i,j,k,Temp)**n
              runStats(i,j,k,num+1) = runStats(i,j,k,num+1)
     $                           + dt * state(i,j,k,Rho) * state(i,j,k,Temp)**n
              num = num + 2
            enddo

            runStats(i,j,k,num) = runStats(i,j,k,num)
     $                         + dt * state(i,j,k,xVel) * state(i,j,k,Temp)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1)
     $                         + dt * state(i,j,k,yVel) * state(i,j,k,Temp)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2)
     $                         + dt * state(i,j,k,zVel) * state(i,j,k,Temp)
            runStats(i,j,k,num+3) = runStats(i,j,k,num+3)
     $                         + dt * state(i,j,k,Rho) * state(i,j,k,xVel)
     $                                                 * state(i,j,k,Temp)
            runStats(i,j,k,num+4) = runStats(i,j,k,num+4)
     $                         + dt * state(i,j,k,Rho) * state(i,j,k,yVel)
     $                                                 * state(i,j,k,Temp)
            runStats(i,j,k,num+5) = runStats(i,j,k,num+5)
     $                         + dt * state(i,j,k,Rho) * state(i,j,k,zVel)
     $                                                 * state(i,j,k,Temp)
            num = num + 6

            enthalpy = state(i,j,k,RhoH) / state(i,j,k,Rho)
            do n = 1, 2
              runStats(i,j,k,num) = runStats(i,j,k,num)
     $                                    + dt * state(i,j,k,Rho) * enthalpy**n
              runStats(i,j,k,num+1) = runStats(i,j,k,num+1)
     $                                             + dt * state(i,j,k,RhoRT)**n
              num = num + 2
            enddo

            runStats(i,j,k,num) = runStats(i,j,k,num)
     $                         + dt * state(i,j,k,RhoH) * state(i,j,k,xVel)
            runStats(i,j,k,num+1) = runStats(i,j,k,num+1)
     $                         + dt * state(i,j,k,RhoH) * state(i,j,k,yVel)
            runStats(i,j,k,num+2) = runStats(i,j,k,num+2)
     $                         + dt * state(i,j,k,RhoH) * state(i,j,k,zVel)
            num = num + 3

            do sp = RhoYlo, RHoYhi
              massfrac = state(i,j,k,sp) / state(i,j,k,Rho)
              runStats(i,j,k,num) = runStats(i,j,k,num)
     $                                    + dt * state(i,j,k,sp)
              runStats(i,j,k,num+1) = runStats(i,j,k,num+1)
     $                                    + dt * state(i,j,k,sp) * massfrac
              runStats(i,j,k,num+2) = runStats(i,j,k,num+2)
     $                               + dt * state(i,j,k,sp) * state(i,j,k,xVel)
              runStats(i,j,k,num+3) = runStats(i,j,k,num+3)
     $                               + dt * state(i,j,k,sp) * state(i,j,k,yVel)
              runStats(i,j,k,num+4) = runStats(i,j,k,num+4)
     $                               + dt * state(i,j,k,sp) * state(i,j,k,zVel)
              num = num + 5
            enddo
          enddo
        enddo
      enddo


      if (num-1 .NE. nStatsExpect) then
        write(*,1020) num-1, nStatsExpect
 1020   format('ERROR: Num-1 = ',I5,'not equal to nStatsExpect = ',I5)
        call bl_abort(" ")
      endif

      end

