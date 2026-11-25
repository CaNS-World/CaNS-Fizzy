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
real(rp), parameter :: eps = epsilon(1._rp)
real(rp), parameter :: small = epsilon(1._rp)*10**(precision(1._rp)/2)
real(rp), parameter :: psi_eps = max(1.e-12_rp,1.e3_rp*epsilon(1._rp)) ! cells with psi_eps < psi < 1-psi_eps hold an
                                                                       ! interface (1.2e-4 in single precision)
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
character(len=3), protected :: forced_hit_type ! 'tgv' (Taylor-Green) or 'abc' (Arnold-Beltrami-Childress)
real(rp), protected :: forced_hit_amp
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
real(rp)        , protected, dimension(0:1,3)   ::  slip_length ! Navier slip length of the walls (zero: no slip)
character(len=1), protected, dimension(0:1,3,3) ::  cbcvel_halo ! (not inputs: conditions imposed on the halo cells of the
real(rp)        , protected, dimension(0:1,3,3) ::   bcvel_halo !  velocity, with the Navier slip of the walls that have it)
character(len=1), protected, dimension(0:1,3)   ::  cbcpre
real(rp)        , protected, dimension(0:1,3)   ::   bcpre
!
real(rp), protected, dimension(3) :: bforce
logical , protected, dimension(3) :: is_forced
real(rp), protected, dimension(3) :: velf ! bulk (mass-averaged) velocity of the initial condition along a forced direction
!
real(rp), protected, dimension(3) :: dl,dli
!
! scalar input parameters
!
integer , protected :: nscal ! number of transported scalars
real(rp), protected, dimension(3) :: gacc
character(len=100), protected, allocatable, dimension(:)     :: iniscal
character(len=1)  , protected, allocatable, dimension(:,:,:) :: cbcscal ! size (0:1,3,nscal)
real(rp)          , protected, allocatable, dimension(:,:,:) ::  bcscal ! size (0:1,3,nscal)
real(rp), protected, allocatable, dimension(:) :: ssource
logical , protected, allocatable, dimension(:) :: is_sforced ! zero-net-heating forcing of each scalar
real(rp), protected, allocatable, dimension(:) :: scalf      ! reference value of the initial scalar field, if forced
logical , protected :: is_boussinesq_buoyancy = .false.
!
! two-fluid input parameters
!
character(len=100), protected     :: inipsi
character(len=1), protected, dimension(0:1,3)   ::  cbcpsi
real(rp)        , protected, dimension(0:1,3)   ::   bcpsi
character(len=1), protected, dimension(0:1,3)   ::  cbccur ! (not inputs: the conditions of the curvature and of the
real(rp)        , protected, dimension(0:1,3)   ::   bccur !  interface normals follow from the others)
character(len=1), protected, dimension(0:1,3,3) ::  cbcnor
real(rp)        , protected, dimension(0:1,3,3) ::   bcnor
real(rp)        , protected, dimension(0:1,3)   :: contact_angle ! static contact angle at walls, in degrees through phase 1
real(rp), protected               :: sigma
real(rp), protected               :: dsigma_ds,s_ref ! surface tension linear in the first scalar: sigma+dsigma_ds*(s-s_ref)
real(rp), protected, dimension(2) :: rho12,mu12
real(rp), protected, allocatable, dimension(:,:) :: ka12,cp12 ! size (2,nscal): conductivity and heat capacity of each scalar
real(rp), protected, dimension(2) :: beta12
real(rp), protected               :: vof_thinc_beta
real(rp), parameter               :: acdi_gam_factor = 1._rp, acdi_gam_min = 1.e-12_rp ! regularization velocity of ACDI: factor
                                                                                       ! of the largest velocity, and lower bound
real(rp), protected               :: psi_thickness_factor
!
! interface method, selected with three inputs; blank entries take the default of the reconstruction:
!   ifc_reconstruction: 'acdi' (ACDI phase field, default), 'thinc' (THINC/QQ volume of fluid) or
!                       'plic' (geometric volume of fluid, piecewise-linear interface)
!   ifc_advection     : 'ab2' (default of 'acdi' and 'thinc') or
!                       'split' (conservative direction-split advection; default of 'plic')
!   ifc_transport     : transport of momentum and scalar: 'ab2' (default with the 'ab2' advection), or, with the
!                       'split' advection, 'one-step' (default) or 'sweeps' (along with the sweeps)
! the valid combinations are checked in `mod_sanity`
!
character(len=10), protected :: ifc_reconstruction,ifc_advection,ifc_transport
integer , protected :: nmarker ! number of markers (phase fields) that the dispersed phase is distributed over
logical , protected :: is_sdf_normals ! normals (and THINC/QQ reconstruction) from a signed-distance function
logical , protected :: is_brackbill_normalization ! density-based scaling of the surface tension force (Brackbill et al.)
logical , protected :: is_hf_normals  ! normals of the geometric reconstruction from height functions
logical , protected :: is_ifc_acdi,is_ifc_thinc,is_ifc_plic,is_ifc_split,is_vof_sweep_transport ! not inputs
real(rp), protected               :: rho0 ! not an input
real(rp), protected               :: visc0 ! reference diffusivities of the constant-coefficient implicit diffusion
real(rp), protected, allocatable, dimension(:) :: alpha0 ! (momentum, and each scalar); not inputs
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
                      is_gridpoint_natural_channel = .false., &
                      is_harmonic_interpolation = .false., & ! harmonic (instead of arithmetic) interpolation of the viscosity
                                                             ! and thermal conductivity between the two fluids
                      is_visc_limiter = .false. ! (true by default with the sharp interface of `ifc_reconstruction = 'plic'`)
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
real(rp), protected, dimension(2) :: hypre_tol = 1.e-6_rp
!
! halo widths of the velocity, pressure, scalar, and interface fields; the grid arrays take the largest
! (one layer suffices for CaNS; extensions with wider stencils may increase them here); that of the phase
! indicator follows the interface method, so it is set, with the largest one, when the input is read
!
integer, parameter :: nh_v = 1, nh_p = 1, nh_s = 1, nh_ifc = 1
integer, protected :: nh_psi = 1, nh_g = max(nh_v,nh_p,nh_s,nh_ifc)
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
    integer :: iunit,ierr,numarg,idir,ib,is
    character(len=1024) :: c_iomsg
    logical :: is_io_fallback
    namelist /dns/ &
                  ng, &
                  l, &
                  gtype,gr, &
                  cfl,dtmax,dt_f, &
                  is_solve_ns,is_track_interface, &
                  inivel, &
                  is_wallturb,is_forced_hit,forced_hit_type,forced_hit_amp, &
                  nstep,time_max,tw_max, &
                  stop_type, &
                  restart,is_overwrite_save,nsaves_max, &
                  icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                  tout0d,tout1d,tout2d,tout3d,tsave, &
                  cbcvel,cbcpre,bcvel,bcpre, &
                  slip_length, &
                  bforce, &
                  is_forced,velf, &
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
                     iniscal, &
                     cbcscal,bcscal, &
                     ssource, &
                     is_sforced,scalf, &
                     is_boussinesq_buoyancy
    namelist /two_fluid/ &
                        inipsi, &
                        cbcpsi,bcpsi,contact_angle, &
                        sigma,dsigma_ds,s_ref,rho12,mu12, &
                        ka12,cp12,beta12, &
                        psi_thickness_factor, &
                        ifc_reconstruction,ifc_advection,ifc_transport, &
                        nmarker, &
                        is_sdf_normals,is_hf_normals, &
                        is_brackbill_normalization
    namelist /numerics/ &
                       is_impdiff,is_impdiff_1d, &
                       is_poisson_dtdma, &
                       is_gridpoint_natural_channel, &
                       is_harmonic_interpolation, &
                       is_visc_limiter
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
    cfl = 0.95_rp; dtmax = 1.e9; dt_f = -1.
    is_solve_ns = .true.; is_track_interface = .true.
    inivel = 'zer'
    is_wallturb = .false.; is_forced_hit = .false.; forced_hit_type = 'tgv'; forced_hit_amp = 0.5
    nstep = 1000; time_max = 1.; tw_max = 0.5
    stop_type = [.true.,.false.,.false.]
    restart = .false.; is_overwrite_save = .true.; nsaves_max = 0
    icheck = 10; iout0d = 10; iout1d = 100; iout2d = 1000; iout3d = 500; isave = 1000
    tout0d = 0.; tout1d = 0.; tout2d = 0.; tout3d = 0.; tsave = 0.
    cbcvel(:,:,:) = 'P'; cbcpre(:,:) = 'P'; bcvel(:,:,:) = 0.; bcpre(:,:) = 0.
    slip_length(:,:) = 0.
    bforce(:) = 0.; is_forced(:) = .false.; velf(:) = 1.
    gacc(:) = 0.
    nscal = 0
    dims(:) = 0
    !
    inipsi = 'uni'
    cbcpsi(:,:) = 'P'; cbcnor(:,:,:) = 'P'; bcpsi(:,:) = 0.; bcnor(:,:,:) = 0.
    cbccur(:,:) = 'P'; bccur(:,:) = 0.
    contact_angle(:,:) = 90.
    sigma = 0.; dsigma_ds = 0.; s_ref = 0.; rho12(:) = 1.; mu12(:) = 0.01_rp
    beta12(:) = 1.
    psi_thickness_factor = -1. ! (default set below)
    ifc_reconstruction = 'acdi'; ifc_advection = ''; ifc_transport = ''
    nmarker = 1
    is_sdf_normals = .true.; is_hf_normals = .false.
    is_brackbill_normalization = .true.
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
      !
      ! conditions imposed on the halo cells of the velocity: Navier slip (Robin condition, type 'R', with the
      ! slip length as its value) for the wall-tangential components at the walls with a positive slip length
      ! (see `mod_sanity` for the requirements)
      !
      cbcvel_halo(:,:,:) = cbcvel(:,:,:)
      bcvel_halo(:,:,:)  = bcvel(:,:,:)
      do idir=1,3
        do ib=0,1
          if(slip_length(ib,idir) > 0.) then
            cbcvel_halo(ib,idir,:) = 'R'
            bcvel_halo(ib,idir,:)  = slip_length(ib,idir)
            cbcvel_halo(ib,idir,idir) = cbcvel(ib,idir,idir)
            bcvel_halo(ib,idir,idir)  = bcvel(ib,idir,idir)
          end if
        end do
      end do
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
      cudecomp_t_comm_backend      = 0
      cudecomp_h_comm_backend      = 0
      cudecomp_is_t_comm_autotune  = .true.
      cudecomp_is_h_comm_autotune  = .true.
      cudecomp_is_t_enable_nccl    = .true.
      cudecomp_is_h_enable_nccl    = .true.
      cudecomp_is_t_enable_nvshmem = .true.
      cudecomp_is_h_enable_nvshmem = .true.
      rewind(iunit)
      read(iunit,nml=cudecomp,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
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
        ! allocate memory
        !
        allocate(iniscal(nscal), &
                 cbcscal(0:1,3,nscal),bcscal(0:1,3,nscal), &
                 ssource(nscal),is_sforced(nscal),scalf(nscal))
        !
        ! set default values
        !
        iniscal(:)     = 'zer'
        cbcscal(:,:,:) = 'P'
        bcscal(:,:,:)  = 0.
        ssource(:)     = 0.
        is_sforced(:)  = .false.
        scalf(:)       = 0.
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
      ! reading two-fluid parameters, if these are set; the properties of the scalars are per scalar (and
      ! allocated for one scalar if there is none, such that the input files may always set them)
      !
      allocate(ka12(2,max(nscal,1)),cp12(2,max(nscal,1)))
      ka12(:,:) = 0.01_rp; cp12(:,:) = 1.
      rewind(iunit)
      read(iunit,nml=two_fluid,iostat=ierr,iomsg=c_iomsg)
      if(ierr /= 0 .and. ierr /= iostat_end) then
        if(myid == 0) print*, 'ERROR: reading `two_fluid` namelist: ', trim(c_iomsg)
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      if(any(rho12(:) <= 0.) .or. any(mu12(:) < 0.) .or. (nscal > 0 .and. any(cp12(:,:) <= 0.)) .or. &
         psi_thickness_factor == 0.) then
        if(myid == 0) print*, 'ERROR: reading the input file: `rho12` and, with `nscal > 0`, `cp12` must be positive, &
                                     &`mu12` must not be negative, and `psi_thickness_factor` must not be zero.'
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        close(iunit)
        error stop
      end if
      !
      ! interface method: defaults of the blank entries, and what follows from it
      !
      is_ifc_acdi  = ifc_reconstruction == 'acdi'
      is_ifc_thinc = ifc_reconstruction == 'thinc'
      is_ifc_plic  = ifc_reconstruction == 'plic'
      if(ifc_advection == '') ifc_advection = merge('split','ab2  ',is_ifc_plic)
      is_ifc_split = ifc_advection == 'split'
      if(ifc_transport == '') ifc_transport = merge('one-step','ab2     ',is_ifc_split)
      is_vof_sweep_transport = ifc_transport == 'sweeps'
      if(is_ifc_plic) is_sdf_normals = .false. ! no signed-distance function with the geometric reconstruction
      if(psi_thickness_factor < 0.) psi_thickness_factor = merge(0.51_rp,0.50_rp,is_ifc_acdi) ! 0.25?
      if(is_ifc_split) nh_psi = 4 ! nine-cell columns of the height functions (geometric reconstruction), and
                                  ! interface band of the transport along with the sweeps
      !
      ! height-function normals: the columns of the donor cells in the halo, one layer, or two with the
      ! transport along with the sweeps
      !
      if(is_ifc_split .and. is_hf_normals) nh_psi = nh_psi + merge(2,1,is_vof_sweep_transport)
      nh_g = max(nh_v,nh_p,nh_s,nh_psi,nh_ifc)
      is_visc_limiter = is_ifc_plic
      !
      cbccur(:,:) = cbcpsi(:,:)
      bccur(:,:) = 0._rp
      !
      ! conditions of the interface normal, from those of the phase indicator: at a wall (zero gradient of
      ! the phase indicator) the interface is normal to it, i.e., the wall-normal component of the interface
      ! normal vanishes and the other two have zero gradient
      !
      do idir=1,3
        do ib=0,1
          select case(cbcpsi(ib,idir))
          case('P')
            cbcnor(ib,idir,:) = 'P'
          case default
            cbcnor(ib,idir,:) = 'N'
            if(cbcpsi(ib,idir) == 'N') cbcnor(ib,idir,idir) = 'D'
          end select
        end do
      end do
      bcnor(:,:,:) = 0.
      rho0 = 1. ! reference density of the pressure splitting; one with the variable-coefficient solver, where it has no role
#if defined(_CONSTANT_COEFFS_POISSON)
      rho0 = minval(rho12(:))
#endif
      visc0  = maxval(mu12(:)/rho12(:))
      allocate(alpha0(max(nscal,1)))
      alpha0(:) = 0.
      do is=1,nscal
        alpha0(is) = maxval(ka12(:,is)/(rho12(:)*cp12(:,is)))
      end do
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
  !
  pure elemental logical function is_contact_angle(theta)
    !
    ! true for a contact angle (in degrees) other than 90 degrees, to the tolerance `psi_eps` on its cosine
    !
    implicit none
    real(rp), intent(in) :: theta
    is_contact_angle = abs(cos(theta*pi/180._rp)) > psi_eps
  end function is_contact_angle
end module mod_param
