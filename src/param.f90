! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_param
use mod_types
#if defined(_OPENACC) && !defined(_USE_DIEZDECOMP)
#if defined (_USE_HIP)
use hipdecomp
#else
use cudecomp
#endif
#endif
implicit none
public
!
! parameters
!
real(rp), parameter :: pi = acos(-1._rp)
#if !defined(_EPS_EXACT_ZERO) /* recommended */
real(rp), parameter :: eps = epsilon(1._rp)
#else
real(rp), parameter :: eps = 0._rp
#endif
real(rp), parameter :: small = epsilon(1._rp)*10**(precision(1._rp)/2)
character(len=100), parameter :: datadir = 'data/'
!
! variables to be determined from the input file
!
integer , protected, dimension(3) :: ng
real(rp), protected, dimension(3) :: l
integer , protected :: gtype
real(rp), protected :: gr
real(rp), protected :: cfl,dtmax,dt_f
logical , protected :: is_solve_ns,is_track_interface
!
character(len=100), protected :: inivel
logical, protected :: is_wallturb,is_forced_hit
!
integer , protected :: nstep
real(rp), protected :: time_max,tw_max
logical , protected, dimension(3) :: stop_type
logical , protected :: restart,is_overwrite_save
integer , protected :: nsaves_max
integer , protected :: icheck,iout0d,iout1d,iout2d,iout3d,isave
real(rp), protected :: tout0d,tout1d,tout2d,tout3d,tsave
!
integer , dimension(2) :: dims ! not protected -> overwritten when autotuning
integer , protected :: ipencil_axis = 1
!
integer, dimension(0:1,3) :: nb
logical, dimension(0:1,3) :: is_bound
character(len=1), protected, dimension(0:1,3,3) ::  cbcvel
real(rp)        , protected, dimension(0:1,3,3) ::   bcvel
character(len=1), protected, dimension(0:1,3)   ::  cbcpre
real(rp)        , protected, dimension(0:1,3)   ::   bcpre
!
real(rp), protected, dimension(3) :: bforce
logical , protected, dimension(3) :: is_forced
!
real(rp), protected, dimension(3) :: dl,dli
!
! scalar input parameters
!
integer , protected :: nscal ! number of transported scalars (0 or 1)
real(rp), protected, dimension(3) :: gacc
character(len=100), protected                 :: inisca
character(len=1)  , protected, dimension(0:1,3) :: cbcsca
real(rp)          , protected, dimension(0:1,3) ::  bcsca
real(rp), protected :: ssource
logical , protected :: is_boussinesq_buoyancy = .false.
!
! two-fluid input parameters
!
character(len=100), protected     :: inipsi
character(len=1), protected, dimension(0:1,3)   ::  cbcpsi
real(rp)        , protected, dimension(0:1,3)   ::   bcpsi
character(len=1), protected, dimension(0:1,3)   ::  cbccur
real(rp)        , protected, dimension(0:1,3)   ::   bccur
character(len=1), protected, dimension(0:1,3,3) ::  cbcnor
real(rp)        , protected, dimension(0:1,3,3) ::   bcnor
real(rp), protected               :: sigma
real(rp), protected, dimension(2) :: rho12,mu12
real(rp), protected, dimension(2) :: ka12,cp12,beta12
real(rp), protected               :: acdi_gam_factor,acdi_gam_min,vof_thinc_beta
real(rp), protected               :: psi_thickness_factor
real(rp), protected               :: rho0 ! not an input
!
#if defined(_OPENACC)
!
! cuDecomp input parameters
!
integer, protected :: cudecomp_t_comm_backend,cudecomp_h_comm_backend
logical, protected :: cudecomp_is_t_comm_autotune ,cudecomp_is_h_comm_autotune , &
                      cudecomp_is_t_enable_nccl   ,cudecomp_is_h_enable_nccl   , &
                      cudecomp_is_t_enable_nvshmem,cudecomp_is_h_enable_nvshmem, &
                      cudecomp_is_t_in_place
logical, protected :: is_use_diezdecomp = .false.,is_diezdecomp_x2z_z2x_transposes = .false.
#endif
!
! other options: debugging/benchmarking
!
logical, protected :: is_debug = .true., is_debug_poisson = .false., &
                      is_timing = .true., &
#if defined(_CONSTANT_COEFFS_POISSON)
                      is_mask_divergence_check = .false.
#else
                      is_mask_divergence_check = .true. ! the divergence follows the tolerance of the iterative solver
#endif
!
! other options: numerics
!
logical, protected :: is_impdiff = .false., is_impdiff_1d = .false., &
                      is_poisson_dtdma = .false., &
                      is_gridpoint_natural_channel = .false.
!
! HYPRE solver options for the variable-coefficient problems; each input has two entries:
! the first for the pressure, the second for the (implicit) diffusion
!
integer , protected, dimension(2) :: hypre_solver_i  = 5, &      ! 1: SMG, 2: PFMG, 3: GMRES, 4: BiCGSTAB, 5: FlexGMRES
                                     hypre_maxiter   = 50, &
                                     hypre_precond_i = 1, &      ! 0: SMG, 1: PFMG, 9: none
                                     hypre_precond_maxiter = 1, &
                                     hypre_gmres_k_dim     = 0, & ! 0: HYPRE's default
                                     hypre_pfmg_relax_type = -1, & ! < 0: red/black Gauss-Seidel
                                     hypre_pfmg_num_pre_relax  = 1, &
                                     hypre_pfmg_num_post_relax = 1, &
                                     hypre_pfmg_max_levels = -1    ! 0: HYPRE's default, < 0: see `configure_pfmg`
real(rp), protected, dimension(2) :: hypre_tol = 1.e-6
!
! halo widths of the velocity, pressure, scalar, and interface fields; the grid arrays take the largest
! (one layer suffices for CaNS; extensions with wider stencils may increase them here)
!
integer, parameter :: nh_v = 1, nh_p = 1, nh_s = 1, nh_psi = 1, nh_g = max(nh_v,nh_p,nh_s,nh_psi)
!
! i/o backend parameters
!
character(len=16), protected :: io_backend = 'mpiio'
character(len=6) , protected :: io_ext = '.bin'
logical          , protected :: is_use_compression = .false.
logical          , protected :: is_full_restart = .false.

contains
  subroutine read_input(myid)
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use mpi
    implicit none
    character(len=255) :: input_file,filename
    integer, intent(in) :: myid
    integer :: iunit,ierr,numarg
    character(len=1024) :: c_iomsg
    logical :: is_io_fallback
    namelist /dns/ &
                  ng, &
                  l, &
                  gtype,gr, &
                  cfl,dtmax,dt_f, &
                  is_solve_ns,is_track_interface, &
                  inivel, &
                  is_wallturb,is_forced_hit, &
                  nstep,time_max,tw_max, &
                  stop_type, &
                  restart,is_overwrite_save,nsaves_max, &
                  icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                  tout0d,tout1d,tout2d,tout3d,tsave, &
                  cbcvel,cbcpre,bcvel,bcpre, &
                  bforce, &
                  is_forced, &
                  gacc, &
                  nscal, &
                  dims,ipencil_axis
    namelist /io/  &
                  io_backend, &
                  is_use_compression, &
                  is_full_restart
#if defined(_OPENACC) && !defined(_USE_DIEZDECOMP)
    namelist /cudecomp/ &
                       cudecomp_t_comm_backend,cudecomp_is_t_enable_nccl,cudecomp_is_t_enable_nvshmem, &
                       cudecomp_h_comm_backend,cudecomp_is_h_enable_nccl,cudecomp_is_h_enable_nvshmem
#endif
    namelist /scalar/ &
                     inisca, &
                     cbcsca,bcsca, &
                     ssource, &
                     is_boussinesq_buoyancy
    namelist /two_fluid/ &
                        inipsi, &
                        cbcpsi,cbcnor,bcnor,bcpsi, &
                        sigma,rho12,mu12, &
                        ka12,cp12,beta12, &
                        psi_thickness_factor
    namelist /numerics/ &
                       is_impdiff,is_impdiff_1d, &
                       is_poisson_dtdma, &
                       is_gridpoint_natural_channel
    namelist /hypre/ &
                    hypre_solver_i,hypre_tol,hypre_maxiter, &
                    hypre_precond_i,hypre_precond_maxiter,hypre_gmres_k_dim, &
                    hypre_pfmg_relax_type,hypre_pfmg_num_pre_relax,hypre_pfmg_num_post_relax, &
                    hypre_pfmg_max_levels
    namelist /other_options/ &
                            is_debug,is_debug_poisson, &
                            is_timing, &
                            is_mask_divergence_check
    !
    ! defaults
    !
    ng(:) = [128,128,128]
    l(:)  = [1. ,1. ,1. ]
    gtype = 1; gr = 0.
    cfl = 0.95; dtmax = 1.e9; dt_f = -1.
    is_solve_ns = .true.; is_track_interface = .true.
    inivel = 'zer'
    is_wallturb = .false.; is_forced_hit = .false.
    nstep = 1000; time_max = 1.; tw_max = 0.5
    stop_type = [.true.,.false.,.false.]
    restart = .false.; is_overwrite_save = .true.; nsaves_max = 0
    icheck = 10; iout0d = 10; iout1d = 100; iout2d = 1000; iout3d = 500; isave = 1000
    tout0d = 0.; tout1d = 0.; tout2d = 0.; tout3d = 0.; tsave = 0.
    cbcvel(:,:,:) = 'P'; cbcpre(:,:) = 'P'; bcvel(:,:,:) = 0.; bcpre(:,:) = 0.
    bforce(:) = 0.; is_forced(:) = .false.
    gacc(:) = 0.
    nscal = 0
    dims(:) = 0
    !
    inisca = 'zer'
    cbcsca(:,:) = 'P'; bcsca(:,:) = 0.
    ssource = 0.
    !
    inipsi = 'uni'
    cbcpsi(:,:) = 'P'; cbcnor(:,:,:) = 'P'; bcpsi(:,:) = 0.; bcnor(:,:,:) = 0.
    cbccur(:,:) = 'P'; bccur(:,:) = 0.
    sigma = 0.; rho12(:) = 1.; mu12(:) = 0.01
    ka12(:) = 0.01; cp12(:) = 1.; beta12(:) = 1.
    psi_thickness_factor = 0.51; acdi_gam_factor = 1.; acdi_gam_min = 1.e-12
#if defined(_INTERFACE_CAPTURING_VOF)
    psi_thickness_factor = 0.50 ! 0.25?
#endif
    !
    ! use a custom input file if one was provided on the command line
    !
    input_file = 'input.nml'
    numarg = command_argument_count()
    if(numarg == 1) then
      call get_command_argument(1,filename)
      input_file = trim(filename)
    end if
    open(newunit=iunit,file=input_file,status='old',action='read',iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'ERROR: reading the input file: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      ! read `dns` namelist
      !
      rewind(iunit)
      read(iunit,nml=dns,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'ERROR: reading `dns` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      dl(:) = l(:)/(1.*ng(:))
      dli(:) = dl(:)**(-1)
      if(all([1,2,3] /= ipencil_axis)) then
        ipencil_axis = 1 ! default to one
        if(myid == 0) print*, 'Warning: prescribed value of `ipencil_axis` is different from 1/2/3.'
        if(myid == 0) print*, 'Defaulting to `ipencil_axis = 1` (x-aligned pencils)...'
      end if
#if defined(_OPENACC)
#if !defined(_USE_DIEZDECOMP)
      !
      ! reading cuDecomp parameters, if these are set
      !
      ! defaults
      !
      cudecomp_is_t_comm_autotune  = .true.
      cudecomp_is_h_comm_autotune  = .true.
      cudecomp_is_t_enable_nccl    = .true.
      cudecomp_is_h_enable_nccl    = .true.
      cudecomp_is_t_enable_nvshmem = .true.
      cudecomp_is_h_enable_nvshmem = .true.
      rewind(iunit)
      read(iunit,nml=cudecomp,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0) then
        if(myid == 0) print*, 'ERROR: reading `cudecomp` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      if(cudecomp_t_comm_backend >= 1 .and. cudecomp_t_comm_backend <= 8) then
        cudecomp_is_t_comm_autotune = .false. ! do not autotune if backend is prescribed
        select case(cudecomp_t_comm_backend)
        case(1)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
        case(2)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P_PL
        case(3)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_A2A
        case(4)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL
        case(5)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NCCL_PL
        case(6)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM
        case(7)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL
#if !defined(_USE_HIP)
        case(8)
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_NVSHMEM_SM
#endif
        case default
          cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
        end select
      end if
      if(cudecomp_h_comm_backend >= 1 .and. cudecomp_h_comm_backend <= 5) then
        cudecomp_is_h_comm_autotune = .false. ! do not autotune if backend is prescribed
        select case(cudecomp_h_comm_backend)
        case(1)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
        case(2)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI_BLOCKING
        case(3)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NCCL
        case(4)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM
        case(5)
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_NVSHMEM_BLOCKING
        case default
          cudecomp_h_comm_backend = CUDECOMP_HALO_COMM_MPI
        end select
      end if
      !
      ! manually set cuDecomp out-of-place transposes by default
      !
      cudecomp_is_t_in_place = .false.
#else
      is_use_diezdecomp = .true.
      is_diezdecomp_x2z_z2x_transposes = .false.
#endif
#endif
      !
      ! reading scalar transport parameters, if these are set
      !
      if(nscal > 0) then
        !
        ! read `scalar` namelist
        !
        rewind(iunit)
        read(iunit,nml=scalar,iostat=ierr,iomsg=c_iomsg)
        if(ierr /= 0) then
          if(myid == 0) print*, 'ERROR: reading `scalar` namelist: ', trim(c_iomsg)
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          close(iunit)
          error stop
        end if
      else
        nscal = 0 ! negative values equivalent to nscal = 0
      end if
      if(is_boussinesq_buoyancy) then
        if(nscal == 0) then
          if(myid == 0) print*, 'ERROR: reading the input file: `is_boussinesq_buoyancy = T` requires `nscal > 0`.'
          if(myid == 0) print*, 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        end if
      end if
      !
      ! reading two-fluid parameters, if these are set
      !
      rewind(iunit)
      read(iunit,nml=two_fluid,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `two_fluid` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      cbccur(:,:) = cbcpsi(:,:)
      bccur(:,:) = bcpre(:,:)
      rho0 = 1.
#if defined(_CONSTANT_COEFFS_POISSON)
      rho0 = minval(rho12(:))
#endif
      vof_thinc_beta = 1._rp/(2._rp*psi_thickness_factor)
      !
      ! read `numerics` namelist
      !
      rewind(iunit)
      read(iunit,nml=numerics,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `numerics` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      if(is_impdiff_1d) is_impdiff = .true.
      !
      ! read `hypre` namelist
      !
      rewind(iunit)
      read(iunit,nml=hypre,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `hypre` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      ! read `other_options` namelist
      !
      rewind(iunit)
      read(iunit,nml=other_options,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `other_options` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      ! read `io` namelist
      !
      rewind(iunit)
      read(iunit,nml=io,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `io` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      is_io_fallback = .false.
      select case(trim(io_backend))
      case('mpiio')
        io_ext = '.bin'
      case('hdf5')
#if defined(_USE_HDF5)
        io_ext = '.h5'
#else
        if(myid == 0) print*, 'Warning: `io_backend = hdf5` requires an HDF5 build.'
        is_io_fallback = .true.
#endif
      case('adios2')
#if defined(_USE_ADIOS2)
        io_ext = '.bp'
#else
        if(myid == 0) print*, 'Warning: `io_backend = adios2` requires an ADIOS2 build.'
        is_io_fallback = .true.
#endif
      case default
        if(myid == 0) print*, 'Warning: unknown `io_backend = ', trim(io_backend), '`.'
        is_io_fallback = .true.
      end select
      if(is_io_fallback) then
        if(myid == 0) print*, 'Defaulting to `io_backend = mpiio`...'
        io_backend = 'mpiio'
        io_ext = '.bin'
      end if
      if(is_use_compression .and. trim(io_backend) == 'mpiio') then
        if(myid == 0) print*, 'Warning: compression is ignored for `io_backend = mpiio`.'
      end if
    close(iunit)
  end subroutine read_input
end module mod_param
