program compute_cv
        implicit none
        integer*8         :: i,step
        integer*8         :: n_particles
        integer*8         :: max_steps,n_equil,n_save
        real*8, parameter :: avogadro  = 6.0221409d23
        real*8, parameter :: boltzmann = 1.38064852d-23
        real*8, parameter :: r         = avogadro*boltzmann
        real*8            :: eps
        real*8            :: sigma
        real*8            :: rcut
        real*8            :: time_step
        real*8            :: density
        real*8            :: molar_mass
        real*8            :: temperature
        real*8            :: pressure
        real*8            :: kinetic
        real*8            :: potential
        real*8            :: total_energy
        real*8            :: avg_temp
        real*8            :: avg_pot
        real*8            :: avg_pot_sq
        real*8            :: sig
        real*8            :: cvres
        character         :: comment*18

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
        close(1)

        avg_temp   = 0.d0
        avg_pot    = 0.d0
        avg_pot_sq = 0.d0
        open(1,file="thermo.dat")
        do i=1,max_steps
           if (i <= n_equil) then
               read(1,*) step,temperature,pressure,kinetic,potential,  &
                         total_energy
           else
               read(1,*) step,temperature,pressure,kinetic,potential,  &
                         total_energy
               avg_temp    = avg_temp+temperature
               avg_pot     = avg_pot+potential
               avg_pot_sq  = avg_pot_sq+potential*potential
           end if
         end do
         close(1)
         avg_temp   = avg_temp/dble(max_steps-n_equil)
         avg_pot    = avg_pot/dble(max_steps-n_equil)
         avg_pot_sq = avg_pot_sq/dble(max_steps-n_equil)
         sig        = avg_pot_sq-avg_pot*avg_pot

         cvres      = 1.5d0*r*sig*dble(n_particles)/(1.5d0*r*r         &
                      *avg_temp*avg_temp-sig*dble(n_particles))

         write(*,*) 'Temperature = ',avg_temp,' K'
         write(*,*) 'cv res      = ',cvres,' J/mol/K'

end program compute_cv

