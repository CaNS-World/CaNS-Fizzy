! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_param
use mod_types
!@acc use cudecomp
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
real(rp), parameter, dimension(2,3) :: rkcoeff = reshape([32._rp/60._rp,  0._rp        , &
                                                          25._rp/60._rp, -17._rp/60._rp, &
                                                          45._rp/60._rp, -25._rp/60._rp], shape(rkcoeff))
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
character(len=100), protected :: inivel,inisca
logical, protected :: is_wallturb,is_forced_hit
!
integer , protected :: nstep
real(rp), protected :: time_max,tw_max
logical , protected, dimension(3) :: stop_type
logical , protected :: restart,is_overwrite_save
integer , protected :: nsaves_max
integer , protected :: icheck,iout0d,iout1d,iout2d,iout3d,isave
!
integer , dimension(2) :: dims
!
integer, dimension(0:1,3) :: nb
logical, dimension(0:1,3) :: is_bound
character(len=1), protected, dimension(0:1,3,3) ::  cbcvel
real(rp)        , protected, dimension(0:1,3,3) ::   bcvel
character(len=1), protected, dimension(0:1,3)   ::  cbcpre
real(rp)        , protected, dimension(0:1,3)   ::   bcpre
character(len=1), protected, dimension(0:1,3)   ::  cbcpsi
real(rp)        , protected, dimension(0:1,3)   ::   bcpsi
character(len=1), protected, dimension(0:1,3)   ::  cbcsca
real(rp)        , protected, dimension(0:1,3)   ::   bcsca
character(len=1), protected, dimension(0:1,3,3) ::  cbcnor
real(rp)        , protected, dimension(0:1,3,3) ::   bcnor
!
real(rp), protected, dimension(3) :: bforce,gacc
real(rp), protected               :: ssource
!
real(rp), protected, dimension(3) :: dl,dli
!
! two-fluid input parameters
!
character(len=100), protected     :: inipsi
real(rp), protected               :: sigma
real(rp), protected, dimension(2) :: rho12,mu12
real(rp), protected, dimension(2) :: ka12,cp12,beta12
real(rp), protected               :: acdi_gam_factor,acdi_gam_min,vof_thinc_beta
real(rp), protected               :: psi_thickness_factor
real(rp), protected               :: rho0 ! not an input
#if defined(_OPENACC)
!
! cuDecomp input parameters
!
integer, protected :: cudecomp_t_comm_backend,cudecomp_h_comm_backend
logical, protected :: cudecomp_is_t_comm_autotune ,cudecomp_is_h_comm_autotune , &
                      cudecomp_is_t_enable_nccl   ,cudecomp_is_h_enable_nccl   , &
                      cudecomp_is_t_enable_nvshmem,cudecomp_is_h_enable_nvshmem, &
                      cudecomp_is_t_in_place
logical :: exists
#endif
contains
  subroutine read_input(myid)
    use mpi
    implicit none
    integer, intent(in) :: myid
    integer :: iunit,ierr
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
                  cbcvel,cbcpre,bcvel,bcpre, &
                  bforce,gacc, &
                  dims
    namelist /scalar/ &
                  inisca, &
                  cbcsca,bcsca, &
                  ssource
    namelist /two_fluid/ &
                  inipsi, &
                  cbcpsi,cbcnor,bcnor,bcpsi, &
                  sigma,rho12,mu12, &
                  ka12,cp12,beta12, &
                  psi_thickness_factor
#if defined(_OPENACC)
    namelist /cudecomp/ &
                       cudecomp_t_comm_backend,cudecomp_is_t_enable_nccl,cudecomp_is_t_enable_nvshmem, &
                       cudecomp_h_comm_backend,cudecomp_is_h_enable_nccl,cudecomp_is_h_enable_nvshmem
#endif
    !
    ! set-up default parameters
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
    cbcvel(:,:,:) = 'P'; cbcpre(:,:) = 'P'; bcvel(:,:,:) = 0.; bcpre(:,:) = 0.
    bforce(:) = 0.; gacc(:) = 0
    dims(:) = 0
    !
    inisca = 'zer'
    cbcsca(:,:) = 'P'; bcsca(:,:) = 0.
    ssource = 0.
    !
    inipsi = 'uni'
    cbcpsi(:,:) = 'P'; cbcnor(:,:,:) = 'P'; bcpsi(:,:) = 0.; bcnor(:,:,:) = 0.
    sigma = 0.; rho12(:) = 1.; mu12(:) = 0.01
    ka12(:) = 0.01; cp12(:) = 1.; beta12(:) = 1.
    psi_thickness_factor = 0.51; acdi_gam_factor = 1.; acdi_gam_min = 1.e-12
#if defined(_INTERFACE_CAPTURING_VOF)
    psi_thickness_factor = 0.50 ! 0.25?
#endif
    !
    ! read input file
    !
    open(newunit=iunit,file='input.nml',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,nml=dns,iostat=ierr)
        rewind(iunit)
        read(iunit,nml=scalar,iostat=ierr)
        rewind(iunit)
        read(iunit,nml=two_fluid,iostat=ierr)
      else
        if(myid == 0) print*, 'Error reading the input file'
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        error stop
      end if
    close(iunit)
    !
    dl(:) = l(:)/(1.*ng(:))
    dli(:) = dl(:)**(-1)
    rho0 = 1.
#if defined(_CONSTANT_COEFFS_POISSON)
    rho0 = minval(rho12(:))
#endif
    vof_thinc_beta = 1._rp/(2._rp*psi_thickness_factor)
    !
#if defined(_OPENACC)
    !
    ! read cuDecomp parameter file cudecomp.in, if it exists
    !
    ! defaults
    !
    cudecomp_is_t_comm_autotune  = .true.
    cudecomp_is_h_comm_autotune  = .true.
    cudecomp_is_t_enable_nccl    = .true.
    cudecomp_is_h_enable_nccl    = .true.
    cudecomp_is_t_enable_nvshmem = .true.
    cudecomp_is_h_enable_nvshmem = .true.
    open(newunit=iunit,file='input.nml',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,nml=cudecomp,iostat=ierr)
      else
        if(myid == 0) print*, 'Error reading the input file'
        if(myid == 0) print*, 'Aborting...'
        call MPI_FINALIZE(ierr)
        error stop
      end if
    close(iunit)
    if(cudecomp_t_comm_backend >= 1 .and. cudecomp_t_comm_backend <= 7) then
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
      case default
        cudecomp_t_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
      end select
    end if
    if(cudecomp_h_comm_backend >= 1 .and. cudecomp_h_comm_backend <= 4) then
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
#endif
  end subroutine read_input
end module mod_param
