program amatter
   !import Fortran 2008 kinds
   use, intrinsic :: iso_fortran_env
   use active_mod !use active matter module
   use fluid_mod
#ifdef _OPENMP
   use omp_lib
#endif

   implicit none
   ! ACTIVE MATTER PARAMETERS (DEFAULT PARAMS)

   integer(int64) :: np = 80
   integer(int64) :: nworms = 1125
   integer(int64) :: nsteps = 250000
   real(real64) :: fdogic = 0.06d0
   real(real64) :: fdogicwall = 0.0000d0
   real(real64) :: fdep = 1.0d0
   real(real64) :: fdepwall = 0.0d0
   real(real64) :: diss = 0.08d0
   real(real64) :: dt = 0.02d0 ! TODO longer timestep?
   real(real64) :: kspring = 57.146436d0
   real(real64) :: kbend = 40.0d0
   real(real64) :: length0 = 0.8d0
   real(real64) :: rcut = 2.5d0
   real(real64) :: rcutsmall = 2.0d0**(1.0/6.0)
   logical :: fluid_cpl = .false.
   integer(int64) :: save_interval = 1000


   ! LATTICE BOLTZMANN FLUID PARAMETERS
   integer(int64), parameter :: dr = 200
   !future setting for if the fluid needs to have smaller decompositon
   ! e.g. hxcell*3 = lx
   integer(int64), parameter :: resolution = 2
   real(real64), parameter :: rho0 = 1.0d0
   real(real64), parameter :: niu = 0.01!viscosity
   real(real64), parameter :: epsilon = 0.1d0 ! bead-fluid coupling

   ! VARIABLES
   integer(int64):: itime = 0
   integer(int64):: i, tid!,nx,ny
   real(real64) :: t, calc_ftime, calc_stime, io_ftime, io_stime
   real(real64) :: elapsed_calc = 0.0d0
   real(real64) :: elapsed_io = 0.0d0
   
   ! INITALIZATION
   !call setup_args(np,nworms,nsteps,save_interval,fdogic,fdogicwall,fdep,&
   !                fdepwall,diss,dt,kspring,kbend,length0,rcut,rcutsmall,fluid_cpl)
   call init_worms(nsteps, dt, nworms, np, rcut, rcutsmall, length0, fluid_cpl,resolution)
   if (fluid_cpl) then
      call init_fluid(nxcell*resolution, nycell*resolution, rho0, niu,resolution)
      print *, ubound(rho), lbound(rho)
   end if
   call write_xyz(itime,nworms, np)

#ifdef _OPENMP
!$omp parallel
   tid = omp_get_thread_num()
   if (tid.eq.0) print*,"OpenMP Threads: ",omp_get_num_threads()
!$omp end parallel
#endif

   ! MAIN LOOP
   main_loop: do itime = 1, nsteps
      t = float(itime)*dt
      if (mod(itime, 100) .eq. 0) then
         print *, itime, elapsed_calc + elapsed_io, "s Iter/s:", itime/(elapsed_calc + elapsed_io)
      end if
      call cpu_time(calc_stime)
      !first update positions and store old forces
      call update_pos(nworms, np, dt)

      !calculate forces
      call calc_forces(nworms, np, kspring, kbend, length0)

      !worm-wall interactions here, and dissipation force proportional to velocity
      call worm_wall(nworms, np, fdogic, fdogicwall, diss)

      !add calculation of vxave and vyave
      do i = 1, ncells
         hhead(i) = -1
      end do
      !assign particles to their cells
      call cell_sort(nworms, np, fdogic, fdep)

      ! evolving the fluid substrate
      if (fluid_cpl) then
         call fluid_step(itime, save_interval)
      end if
      !update velocities and normalize average velocities (they are one time step behind)
      call update_vel(nworms, np)
      if (fluid_cpl) then
         call bead_fluid_interaction(nworms, np,resolution,epsilon)
      endif

      call cpu_time(calc_ftime)
      elapsed_calc = elapsed_calc + (calc_ftime - calc_stime)

      call cpu_time(io_stime)
      if (mod(itime, save_interval) .eq. 0) then
         print *, "writing data at:", itime
         !if (fluid_cpl) call write_out(itime,nxcell*resolution,nycell*resolution)
         call write_xyz(itime,nworms, np)
         
         !call write_active_vel(itime,nxcell*resolution,nycell*resolution)
      end if
      call cpu_time(io_ftime)
      elapsed_io = elapsed_io + (io_ftime - io_stime)
   end do main_loop

   call finalize(fluid_cpl,elapsed_calc, elapsed_io)
   

   contains
           subroutine setup_args(np,nworms,nsteps,save_interval,fdogic,fdogicwall,fdep,&
                                 fdepwall,diss,dt,kspring,kbend,length0,rcut,rcutsmall,&
                                 fluid_cpl)
                   integer(int64), intent(in out) :: np,nworms,nsteps,save_interval
                   real(real64), intent(in out) :: fdogic,fdogicwall,fdep,fdepwall,diss,dt,&
                                                   kspring,kbend,length0,rcut,rcutsmall
                   logical, intent(in out) :: fluid_cpl
                   !subroutine variables
                   character(len=100) :: buffer, label
                   character(len=100) :: arg
                   integer :: pos
                   integer, parameter :: fh = 15
                   integer :: ios = 0
                   integer :: line = 0
                  
                   ! this assumes that the parameter file is the second argument
                   call get_command_argument(1,arg)
                   if (len_trim(arg) .eq. 0) then
                        print*,"NO INPUT FILE"
                        stop
                   endif
                   print*,trim(arg)
                   open(fh, file=trim(arg))

                   ! ios is negative if an end of record condition is encountered
                   ! or if an endfile condition was detected. It is positive if
                   ! an error was detected. ios is zero otherwise.

                   do while (ios .eq. 0)
                        read(fh,'(A)',iostat=ios) buffer
                        if (ios .eq. 0) then
                                line = line + 1
                                ! find the first instance of whilespace, split label and data
                                pos = scan(buffer,'=')
                                label = buffer(1:pos-1)
                                buffer = buffer(pos+1:)
                                !print*,label,buffer
                                select case(label)
                                case('np')
                                        read(buffer, *, iostat=ios) np
                                        print*,"np:",np
                                case('nworms')
                                        read(buffer, *, iostat=ios) nworms
                                        print*,"nworms:",nworms
                                case('nsteps')
                                        read(buffer, *, iostat=ios) nsteps
                                        print*,"nsteps",nsteps
                                case('fdogic')
                                        read(buffer, *, iostat=ios) fdogic
                                        print*,"fdogic:",fdogic
                                case('fdogicwall')
                                        read(buffer,*,iostat=ios) fdogicwall
                                        print*,"fdogicwall",fdogicwall
                                case('fdep')
                                        read(buffer, *, iostat=ios) fdep
                                        print*,"fdep:",fdep
                                case('fdepwall')
                                        read(buffer, *, iostat=ios) fdepwall
                                        print*,"fdepwall:",fdepwall
                                case('diss')
                                        read(buffer, *, iostat=ios) diss
                                        print*,"diss:",diss
                                case('dt')
                                        read(buffer, *, iostat=ios) dt
                                        print*,"dt:",dt
                                case('kspring')
                                        read(buffer, *, iostat=ios) kspring
                                        print*,"kspring:",kspring
                                case('kbend')
                                        read(buffer,*,iostat=ios) kbend
                                        print*,"kbend:",kbend
                                case("length0")
                                        read(buffer,*,iostat=ios) length0
                                        print*,"length0:",length0
                                case("rcut")
                                        read(buffer,*,iostat=ios) rcut
                                        print*,"rcut:",rcut
                                case("rcutsmall")
                                        read(buffer,*,iostat=ios) rcutsmall
                                        print*,"rcutsmall:",rcutsmall
                                case("fluid_cpl")
                                        read(buffer,*,iostat=ios) fluid_cpl
                                        print*,"fluid_cpl:",fluid_cpl
                                case("save_interval")
                                        read(buffer,*,iostat=ios) save_interval
                                        print*,"save_interval:",save_interval
                                case default
                                        print*,"Skipping invalid label at line,",line
                                end select
                        end if
                end do
                return
           end subroutine setup_args
end program amatter
