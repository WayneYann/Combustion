        program readplt

        implicit none

        integer nx
        integer n
        integer nsteps
        real*8 x
        real*8 time

        real*8 spec1,spec2,spec3
        real*8 rho,temperature,rhoh,vel
        real*8 rhort

        character pltfile*(8)
        character velfile*(8)
        character rhofile*(8)
        character tmpfile*(9)
        character rhohfile*(9)
        character specfile*(9)
        character rhortfile*(10)
        character char_of_int*(5)

        integer i
        pltfile(1:3) = 'plt'
        velfile(1:3) = 'VEL'
        tmpfile(1:4) = 'TEMP'
        rhofile(1:3) = 'RHO'
        rhohfile(1:4) = 'RHOH'
        specfile(1:4) = 'SPEC'
        rhortfile(1:5) = 'RHORT'

        read(5,*) nsteps

        write(char_of_int,1005) nsteps
        pltfile(4:8) = char_of_int
        velfile(4:8) = char_of_int
        tmpfile(5:9) = char_of_int
        rhofile(4:8) = char_of_int
        rhohfile(5:9) = char_of_int
        specfile(5:9) = char_of_int
        rhortfile(6:10) = char_of_int
1005    format(i5.5)

        open(11,file=velfile,form='formatted')
        open(12,file=tmpfile,form='formatted')
        open(13,file=rhofile,form='formatted')
        open(14,file=rhohfile,form='formatted')
        open(15,file=specfile,form='formatted')
        open(16,file=rhortfile,form='formatted')

        open(10,file=pltfile,form='formatted')
        print *,'...reading data from ',pltfile
        read(10,*) n
        read(10,*) nx
        read(10,*) time
        print *,'NSTEPS ',n,' AT TIME ',time
        do i = 0,nx
          read(10,*)          x,spec1,spec2,spec3,
     $                        rho,temperature,rhoh,vel,rhort
          write(11,*) x,vel
          write(12,*) x,temperature
          write(13,*) x,rho
          write(14,*) x,rhoh
          write(15,*) x,spec1
          write(16,*) x,rhort
        enddo

        end
