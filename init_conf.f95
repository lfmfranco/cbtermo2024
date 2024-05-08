!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code to generate an initial configuration of spherical particles    !
! in a FCC lattice with initial velocities distributed according to   !
! Maxwell-Boltzmann                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! University of Campinas                                              !
! School of Chemical Engineering                                      !
! Prof. Luis Fernando Mercier Franco                                  !
! Date: May 7th, 2024                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disclaimer:                                                         !
! The author does not accept any liability for the use of this code   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program init_conf
        implicit none
        integer*8           :: i
        integer*8           :: j
        integer*8           :: k
        integer*8           :: particle
        integer*8           :: n_particles
        integer*8           :: unit_cells
        real*8              :: volume
        real*8              :: density
        real*8              :: temperature
        real*8              :: molar_mass
        real*8              :: box_length
        real*8              :: side
        real*8              :: sumv
        real*8              :: norm
        real*8              :: avg_vx
        real*8              :: avg_vy
        real*8              :: avg_vz
        real*8              :: gauss
        real*8, parameter   :: mtoang    = 1.d10
        real*8, parameter   :: avogadro  = 6.0221409d23
        real*8, parameter   :: boltzmann = 1.38064852d-23
        real*8, parameter   :: r         = avogadro*boltzmann
        real*8, allocatable :: rx(:)
        real*8, allocatable :: ry(:)
        real*8, allocatable :: rz(:)
        real*8, allocatable :: vx(:)
        real*8, allocatable :: vy(:)
        real*8, allocatable :: vz(:)
        character           :: comment*18

        open(1,file="file.inp")
        read(1,'(a18,i6)') comment,n_particles
        read(1,'(a18,f12.3)') comment,temperature
        read(1,'(a18,f12.3)') comment,density
        read(1,'(a18,f12.3)') comment,molar_mass
        close(1)

        allocate(rx(n_particles))
        allocate(ry(n_particles))
        allocate(rz(n_particles))
        allocate(vx(n_particles))
        allocate(vy(n_particles))
        allocate(vz(n_particles))

        ! Converting molar mass from g/mol to kg/mol
        molar_mass = 1.d-3*molar_mass

        ! Volume in m^3
        volume     = dble(n_particles)*molar_mass/density/avogadro

        box_length = mtoang*volume**(1.d0/3.d0)
                
        ! Unit cell length in Angstroms
        side       = mtoang*(volume/dble(n_particles))**(1.d0/3.d0)

        ! Number of unit cells in one direction
        unit_cells = nint(dble(n_particles)**(1.d0/3.d0))

        ! Generating positions in a FCC lattice
        particle = 1
        do i=1,unit_cells
           do j=1,unit_cells
              do k=1,unit_cells
                 rx(particle) = dble(i-1)*side
                 ry(particle) = dble(j-1)*side
                 rz(particle) = dble(k-1)*side
                 particle     = particle+1
              end do
           end do
        end do
        ! Centralizing the box
        rx(:)  = rx(:)-0.5d0*box_length
        ry(:)  = ry(:)-0.5d0*box_length
        rz(:)  = rz(:)-0.5d0*box_length

        ! Generating velocities according Maxwell-Boltzmann distribution
        do particle=1,n_particles
           vx(particle)  = gauss()
           vy(particle)  = gauss()
           vz(particle)  = gauss()
        end do
        sumv   = dsqrt(sum(vx(:)**2.d0+vy(:)**2.d0+vz(:)**2.d0))
        norm   = 1d-5*dsqrt(3.d0*dble(n_particles)*r*temperature       &
                 /molar_mass)
        vx(:)  = vx(:)*norm/sumv
        vy(:)  = vy(:)*norm/sumv
        vz(:)  = vz(:)*norm/sumv
        avg_vx = sum(vx)/dble(n_particles)
        avg_vy = sum(vy)/dble(n_particles)
        avg_vz = sum(vz)/dble(n_particles)
        vx(:)  = vx(:)-avg_vx
        vy(:)  = vy(:)-avg_vy
        vz(:)  = vz(:)-avg_vz
       
        ! Writing initial configuration in conf.xyz 
        open(1,file="conf.xyz")
        write(1,*) n_particles
        write(1,*) ''
        do i=1,n_particles
           write(1,*) 'C',rx(i),ry(i),rz(i),vx(i),vy(i),vz(i)
        end do
        close(1)

end program init_conf

function gauss() result (fun)
        implicit none
        integer*4         :: i
        real*8, parameter :: a1  = 3.949846138d0
        real*8, parameter :: a3  = 0.252408784d0
        real*8, parameter :: a5  = 0.076542912d0
        real*8, parameter :: a7  = 0.008355968d0
        real*8, parameter :: a9  = 0.029899776d0
        real*8            :: summ 
        real*8            :: r
        real*8            :: r2
        real*8            :: fun
        real*8            :: random_n

        summ = 0.d0
        do i=1,12
           call random_number(random_n)
           summ = summ+random_n
        end do

        r   = 0.25d0*(summ-6.d0)
        r2  = r*r
        fun = ((((a9*r2+a7)*r2+a5)*r2+a3)*r2+a1)*r
        
end function gauss
