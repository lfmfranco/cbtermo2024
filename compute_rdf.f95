program compute_rdf
        implicit none
        integer*8                            :: i,j
        integer*8                            :: bin
        integer*8                            :: step
        integer*8                            :: n_particles
        integer*8                            :: max_bin 
        integer*8                            :: max_steps
        integer*8                            :: n_save
        integer*8                            :: n_equil
        integer*8, dimension(:), allocatable :: hist
        real*8, parameter                    :: avogadro  = 6.0221409d23
        real*8, parameter                    :: mtoang    = 1.d10
        real*8                               :: time_step
        real*8                               :: density
        real*8                               :: temperature
        real*8                               :: molar_mass
        real*8                               :: eps
        real*8                               :: sigma
        real*8                               :: rcut
        real*8                               :: box_length
        real*8                               :: volume
        real*8                               :: rho
        real*8                               :: delr
        real*8                               :: pi
        real*8, dimension(:), allocatable    :: rx
        real*8, dimension(:), allocatable    :: ry
        real*8, dimension(:), allocatable    :: rz
        real*8                               :: vx
        real*8                               :: vy
        real*8                               :: vz
        real*8                               :: rxi,ryi,rzi,rxj,ryj,rzj
        real*8                               :: rxij,ryij,rzij
        real*8                               :: rijsq,rij
        real*8                               :: cons
        real*8                               :: rlower
        real*8                               :: rupper
        real*8                               :: nideal
        real*8                               :: rdf
        real*8                               :: dummy
        character                            :: atom*1
        character                            :: comment*18

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
        read(1,'(a18,i10)') comment,n_equil
        read(1,'(a18,i10)') comment,max_bin
        close(1)

        allocate(hist(max_bin))

        n_equil    = n_equil/n_save
        max_steps  = max_steps/n_save
        molar_mass = 1.d-3*molar_mass
        box_length = mtoang*(dble(n_particles)*molar_mass              &
                     /density/avogadro)**(1.d0/3.d0)
        volume     = box_length*box_length*box_length
        rho        = dble(n_particles)/volume
        delr       = 0.5d0*box_length/dble(max_bin)
        pi         = 4.d0*datan(1.d0)

        hist(:)    = 0

        allocate(rx(n_particles))
        allocate(ry(n_particles))
        allocate(rz(n_particles))

        open(1,file="traj.xyz")
        do step=1,max_steps
           if (step <= n_equil) then
              read(1,*) n_particles
              do i=1,n_particles
                 read(1,*) atom,dummy,dummy,dummy
              end do
           else 
              read(1,*) n_particles
              do i=1,n_particles
                 read(1,*) atom,rx(i),ry(i),rz(i)
              end do

              do i=1,n_particles-1
                 rxi = rx(i)
                 ryi = ry(i)
                 rzi = rz(i)
                 do j=i+1,n_particles
                    rxj = rx(j)
                    ryj = ry(j)
                    rzj = rz(j)

                    rxij = rxi-rxj
                    ryij = ryi-ryj
                    rzij = rzi-rzj

                    rxij = rxij-box_length*dnint(rxij/box_length)
                    ryij = ryij-box_length*dnint(ryij/box_length)
                    rzij = rzij-box_length*dnint(rzij/box_length)

                    rijsq = rxij*rxij+ryij*ryij+rzij*rzij
                    rij   = dsqrt(rijsq)
                    bin   = dnint(rij/delr)
                    if (bin <= max_bin) then
                        hist(bin) = hist(bin)+2
                    end if
                 end do
              end do
           end if
        end do
        close(1)


        deallocate(rx,ry,rz)

        cons = 4.d0*pi*rho/3.d0

        open(1,file="rdf.dat")
        do bin=1,max_bin
           rlower = delr*dble(bin)
           rupper = rlower+delr
           nideal = cons*(rupper**3.d0-rlower**3.d0)
           rdf    = dble(hist(bin))/nideal/dble(max_steps-n_equil)/dble(n_particles)
           write(1,*) rlower+0.5d0*delr,rdf
        end do
        close(1)

end program compute_rdf
