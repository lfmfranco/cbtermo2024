!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code to execute a NVE MD simulation with velocity-Verlet algorithm  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! University of Campinas                                              !
! School of Chemical Engineering                                      !
! Prof. Luis Fernando Mercier Franco                                  !
! Date: May 7th, 2024                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disclaimer:                                                         !
! The author does not accept any liability for the use of this code   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globalvar
        integer*8           :: n_particles
        integer*8           :: max_steps
        integer*8           :: n_save
        real*8, parameter   :: avogadro  = 6.0221409d23
        real*8, parameter   :: boltzmann = 1.38064852d-23
        real*8, parameter   :: r         = avogadro*boltzmann
        real*8, parameter   :: pi        = 4.d0*datan(1.d0)
        real*8              :: molar_mass
        real*8              :: density
        real*8              :: temperature
        real*8              :: v
        real*8              :: w
        real*8              :: vij
        real*8              :: wij
        real*8              :: rijsq
        real*8              :: box_length
        real*8              :: time_step
        real*8              :: eps
        real*8              :: eps24
        real*8              :: sigma
        real*8              :: sigmasq
        real*8              :: rcut
        real*8              :: rcutsq
        real*8              :: volume
        real*8              :: rho
        real*8              :: ircut
        real*8              :: ircut3
        real*8              :: ircut9
        real*8              :: vlrc
        real*8              :: plrc
        real*8, allocatable :: rx(:)
        real*8, allocatable :: ry(:)
        real*8, allocatable :: rz(:)
        real*8, allocatable :: vx(:)
        real*8, allocatable :: vy(:)
        real*8, allocatable :: vz(:)
        real*8, allocatable :: ax(:)
        real*8, allocatable :: ay(:)
        real*8, allocatable :: az(:)

        contains

        subroutine init_var()
               implicit none
               character :: comment*18

               open(1,file="file.inp")
               read(1,'(a18,i6)') comment,n_particles
               read(1,'(a18,f12.3)') comment,temperature
               read(1,'(a18,f12.3)') comment,density
               read(1,'(a18,f12.3)') comment,molar_mass
               read(1,'(a18,f12.3)') comment,eps
               read(1,'(a18,f12.3)') comment,sigma
               read(1,'(a18,f12.3)') comment,rcut
               read(1,'(a18,f12.3)') comment,time_step
               read(1,'(a18,i10)') comment,max_steps
               read(1,'(a18,i10)') comment,n_save
               close(1)

               allocate(rx(n_particles))
               allocate(ry(n_particles))
               allocate(rz(n_particles))
               allocate(vx(n_particles))
               allocate(vy(n_particles))
               allocate(vz(n_particles))
               allocate(ax(n_particles))
               allocate(ay(n_particles))
               allocate(az(n_particles))
        
               ! Converting from g/mol to kg/mol
               molar_mass = molar_mass*1.d-3

               box_length = 1.d10*(molar_mass*dble(n_particles)        &
                            /density/avogadro)**(1.d0/3.d0)

               ! Epsilon in J/mol
               eps        = eps*r
               eps24      = 24.d-10*eps/molar_mass

               ! Sigma in Angstroms
               sigmasq    = sigma*sigma

               ! Cut-off radius in Angstroms
               rcutsq     = rcut*rcut
        
               volume     = 1.d-30*box_length*box_length*box_length
               rho        = dble(n_particles)/volume
               ircut      = sigma/rcut
               ircut3     = ircut*ircut*ircut
               ircut9     = ircut3*ircut3*ircut3
               vlrc       = 8.d-30*pi*rho*sigma**3.d0*eps              &
                            *(ircut9/3.d0-ircut3)/3.d0
               plrc       = 16.d-30*pi*rho*rho*sigma**3.d0*eps         &
                            *(2.d0*ircut9/3.d0-ircut3)/avogadro/3.d0

        end subroutine init_var

end module globalvar   
      
      
program md
         use globalvar
         implicit none
         integer*8 :: i
         integer*8 :: steps
         real*8    :: pressure
         real*8    :: kinetic
         real*8    :: potential
         real*8    :: total_energy
         character :: atom*1

         call init_var()

         open(1,file="conf.xyz")
         read(1,*) n_particles
         do i=1,n_particles
            read(1,*) atom,rx(i),ry(i),rz(i),vx(i),vy(i),vz(i)
         end do
         close(1)
    
         kinetic     = 0.5d10*molar_mass*sum(vx(:)**2.d0+vy(:)**2.d0   &
                       +vz(:)**2.d0)/dble(n_particles)
         temperature = 2.d0*kinetic/r/3.d0

         call compute_acceleration()
    
         potential    = 4.d0*eps*v/dble(n_particles)+vlrc  
         pressure     = 1.d-6*((dble(n_particles)*boltzmann            &
                        *temperature/volume)                           &
                        +1.d10*eps24*w*molar_mass/volume/avogadro/3.d0 &
                        +plrc)
         total_energy = kinetic+potential

         open(2,file="traj.xyz")
         open(3,file="thermo.dat")
         write(*,*) '    Steps       T(K)        p(MPa)    K(J/mol)',&
                    '     U(J/mol)   E(J/mol)'
         do steps=1,max_steps
    
           write(*,'(i10,5f12.2)') steps,temperature,pressure,kinetic, &
                                   potential,total_energy
           write(3,'(i10,5e15.7)') steps,temperature,pressure,kinetic, &
                                   potential,total_energy

           ! Velocity-Verlet algorithm
           vx(:) = vx(:)+0.5d0*ax(:)*time_step
           vy(:) = vy(:)+0.5d0*ay(:)*time_step
           vz(:) = vz(:)+0.5d0*az(:)*time_step
           rx(:) = rx(:)+vx(:)*time_step
           ry(:) = ry(:)+vy(:)*time_step
           rz(:) = rz(:)+vz(:)*time_step
           call compute_acceleration()
           vx(:) = vx(:)+0.5d0*ax(:)*time_step
           vy(:) = vy(:)+0.5d0*ay(:)*time_step
           vz(:) = vz(:)+0.5d0*az(:)*time_step
       
           kinetic      = 0.5d10*molar_mass*sum(vx(:)**2.d0+vy(:)**2.d0&
                          +vz(:)**2.d0)/dble(n_particles)
           temperature  = 2.d0*kinetic/r/3.d0
           potential    = 4.d0*eps*v/dble(n_particles)+vlrc
           pressure     = 1.d-6*((dble(n_particles)*boltzmann          &
                          *temperature/volume)                         &
                          +1.d10*eps24*w*molar_mass/volume/avogadro/3.d0&
                          +plrc)
           total_energy = kinetic+potential

           if (mod(steps,n_save) == 0) then
              write(2,*) n_particles
              write(2,*) ''
              do i=1,n_particles
                 write(2,*) 'C',rx(i),ry(i),rz(i)
              end do
           end if
          
    end do
    close(2)
    close(3)

end program md

subroutine compute_acceleration()
        use globalvar
        implicit none
        integer*8 :: i
        integer*8 :: j
        real*8    :: rxi
        real*8    :: ryi
        real*8    :: rzi
        real*8    :: rxj
        real*8    :: ryj
        real*8    :: rzj
        real*8    :: rxij
        real*8    :: ryij
        real*8    :: rzij
        real*8    :: aij
        real*8    :: axi
        real*8    :: ayi
        real*8    :: azi
    
        ax(:) = 0.d0
        ay(:) = 0.d0
        az(:) = 0.d0
        v     = 0.d0
        w     = 0.d0
    
        do i=1,n_particles-1
           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)
           axi = ax(i)
           ayi = ay(i)
           azi = az(i)
           do j=i+1,n_particles
              rxj   = rx(j)
              ryj   = ry(j)
              rzj   = rz(j)
              rxij  = rxi-rxj
              ryij  = ryi-ryj
              rzij  = rzi-rzj
              ! Minimum image convention
              rxij  = rxij-box_length*dnint(rxij/box_length)
              ryij  = ryij-box_length*dnint(ryij/box_length)
              rzij  = rzij-box_length*dnint(rzij/box_length)
              rijsq = rxij*rxij+ryij*ryij+rzij*rzij

              if (rijsq <= rcutsq) then
                 call compute_potential()
                 v     = v+vij
                 w     = w+wij
                 aij   = wij/rijsq
                 axi   = axi+aij*rxij
                 ayi   = ayi+aij*ryij
                 azi   = azi+aij*rzij
                 ax(j) = ax(j)-aij*rxij
                 ay(j) = ay(j)-aij*ryij
                 az(j) = az(j)-aij*rzij
              end if
           end do
           ax(i) = axi
           ay(i) = ayi
           az(i) = azi
        end do

        ax(:) = eps24*ax(:)
        ay(:) = eps24*ay(:)
        az(:) = eps24*az(:)

end subroutine compute_acceleration

subroutine compute_potential()
        use globalvar
        implicit none
        real*8 :: sr2
        real*8 :: sr6
        real*8 :: sr12

        sr2  = sigmasq/rijsq
        sr6  = sr2*sr2*sr2
        sr12 = sr6*sr6
        vij  = sr12-sr6
        wij  = 2.d0*sr12-sr6

end subroutine compute_potential
