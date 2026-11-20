! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
!
!        CCCCCCCCCCCCC                    NNNNNNNN        NNNNNNNN    SSSSSSSSSSSSSSS
!     CCC::::::::::::C                    N:::::::N       N::::::N  SS:::::::::::::::S
!   CC:::::::::::::::C                    N::::::::N      N::::::N S:::::SSSSSS::::::S
!  C:::::CCCCCCCC::::C                    N:::::::::N     N::::::N S:::::S     SSSSSSS
! C:::::C       CCCCCC   aaaaaaaaaaaaa    N::::::::::N    N::::::N S:::::S
!C:::::C                 a::::::::::::a   N:::::::::::N   N::::::N S:::::S
!C:::::C                 aaaaaaaaa:::::a  N:::::::N::::N  N::::::N  S::::SSSS
!C:::::C                          a::::a  N::::::N N::::N N::::::N   SS::::::SSSSS
!C:::::C                   aaaaaaa:::::a  N::::::N  N::::N:::::::N     SSS::::::::SS
!C:::::C                 aa::::::::::::a  N::::::N   N:::::::::::N        SSSSSS::::S
!C:::::C                a::::aaaa::::::a  N::::::N    N::::::::::N             S:::::S
! C:::::C       CCCCCC a::::a    a:::::a  N::::::N     N:::::::::N             S:::::S
!  C:::::CCCCCCCC::::C a::::a    a:::::a  N::::::N      N::::::::N SSSSSSS     S:::::S
!   CC:::::::::::::::C a:::::aaaa::::::a  N::::::N       N:::::::N S::::::SSSSSS:::::S
!     CCC::::::::::::C  a::::::::::aa:::a N::::::N        N::::::N S:::::::::::::::SS
!        CCCCCCCCCCCCC   aaaaaaaaaa  aaaa NNNNNNNN         NNNNNNN  SSSSSSSSSSSSSSS
!-------------------------------------------------------------------------------------
! CaNS -- Canonical Navier-Stokes Solver
!-------------------------------------------------------------------------------------
program cans
  use, intrinsic :: iso_fortran_env, only: compiler_version,compiler_options
  use, intrinsic :: iso_c_binding  , only: C_PTR
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan
  use mpi
  use decomp_2d
  use mod_bound          , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,ierr,dinfo_dtdma
  use mod_correc         , only: correc
  use mod_fft            , only: fftini,fftend
  use mod_fillps         , only: fillps
  use mod_forcing        , only: lscale_forcing
  use mod_initflow       , only: initflow,initscal
  use mod_initgrid       , only: initgrid,save_grid
  use mod_initmpi        , only: initmpi
  use mod_initsolver     , only: initsolver
  use mod_load           , only: load_one
  use mod_mom            , only: cmpt_momentum
  use mod_rk             , only: tm => rk,tm_scal => rk_scal,tm_2fl => rk_2fl
  use mod_output         , only: is_output_due, &
                                 out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  use mod_param          , only: ng,l,dl,dli, &
                                 gtype,gr, &
                                 cfl,dtmax,dt_f, &
                                 is_solve_ns,is_track_interface, &
                                 inivel,inisca,inipsi, &
                                 is_wallturb,is_forced_hit, &
                                 nstep,time_max,tw_max,stop_type, &
                                 restart,is_overwrite_save,nsaves_max, &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 tout0d,tout1d,tout2d,tout3d,tsave, &
                                 cbcvel,bcvel,cbcpre,bcpre, &
                                 cbcsca,bcsca,cbcpsi,bcpsi,cbcnor,bcnor,cbccur,bccur, &
                                 is_forced,bforce,ssource, &
                                 gacc,nscal, &
                                 rho0,rho12,mu12,sigma,ka12,cp12,beta12, &
                                 psi_thickness_factor, &
                                 acdi_gam_factor,acdi_gam_min, &
                                 vof_thinc_beta, &
                                 dims, &
                                 nb,is_bound, &
                                 small, &
                                 datadir,io_ext, &
                                 read_input, &
                                 is_debug,is_debug_poisson, &
                                 is_timing, &
                                 is_impdiff, &
                                 is_poisson_dtdma, &
                                 is_mask_divergence_check, &
                                 nh_v,nh_p,nh_s,nh_psi,nh_g
  use mod_sanity         , only: test_sanity_input,test_sanity_solver
  use mod_solve_helmholtz, only: rhs_bound
  use mod_acdi           , only: acdi_set_epsilon
#if !defined(_INTERFACE_CAPTURING_VOF)
  use mod_acdi           , only: acdi_set_gamma,acdi_cmpt_phi
#elif defined(_SDF_NORMALS)
  use mod_vof_thinc_qq   , only: vof_thinc_cmpt_phi
#endif
  use mod_two_fluid      , only: init2fl,cmpt_norm_curv => cmpt_norm_curv_youngs,bulk_mean_12_stag
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mod_solver_vc      , only: solver_vc
#endif
#if !defined(_OPENACC)
  use mod_solver         , only: solver
#else
  use mod_solver_gpu     , only: solver => solver_gpu
  use mod_workspaces     , only: init_wspace_arrays,set_cufft_wspace,cudecomp_finalize
  use mod_common_cudecomp, only: istream_acc_queue_1,ap_z_dtdma
#endif
  use mod_updatep        , only: updatep
  use mod_posix          , only: rename_file,read_link
  use mod_utils          , only: bulk_mean
#if defined(_OPENACC)
  use mod_utils          , only: device_memory_footprint
#endif
  use mod_types
  use omp_lib
  implicit none
  integer , dimension(3) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
  real(rp), target, allocatable, dimension(:,:,:) :: u,v,w,p,pp
  real(rp), allocatable, dimension(:,:,:) :: pn,po
  real(rp), dimension(3) :: rho_av
  logical , dimension(3) :: is_cmpt_rho_av
#if !defined(_OPENACC) || defined(_USE_HIP)
  type(C_PTR), dimension(2,2) :: arrplanp
#else
  integer    , dimension(2,2) :: arrplanp
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  integer , dimension(3) :: n_z_d
  real(rp), allocatable, dimension(:,:,:) :: ap_d,cp_d
  logical :: is_dtdma_update_p
  real(rp) :: normfftp
  type(rhs_bound) :: rhsbp
  real(rp) :: alpha
  !
  real(rp) :: dt,dto,dt_r,dti,dt_cfl,time,time_start,dtrk,dtrki,divtot,divmax
  real(rp) :: gam,seps
  integer :: istep,iunit
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g, &
                                         xc_g,xf_g,yc_g,yf_g, &
                                         grid_vol_ratio_c,grid_vol_ratio_f
  real(rp), dimension(3) :: dpdl,meanvel,momentum
  real(rp), dimension(3,6) :: force_balance
  real(rp), dimension(2) :: tm_coeff
  !
  ! scalar field
  !
  real(rp), target, allocatable, dimension(:,:,:) :: s
  integer :: is,nio
  !
  ! two-fluid solver specific
  !
  real(rp), target, allocatable, dimension(:,:,:) :: psi
  real(rp), allocatable, dimension(:,:,:) :: psio,phi,kappa,normx,normy,normz, &
                                             psiflx_x,psiflx_y,psiflx_z
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
  real(rp), allocatable, dimension(:,:,:) :: surfx_n,surfy_n,surfz_n, &
                                             surfx_o,surfy_o,surfz_o
#endif
  !
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(42) :: var
  !
  real(rp) :: dt12,dt12av,dt12min,dt12max
  real(dp) :: twi,tw
  !
  integer  :: savecounter
  character(len=7  ) :: fldnum
  character(len=4  ) :: chkptnum
  character(len=100) :: filename
  character(len=256) :: fname_save
  logical :: is_save_tmp
  integer :: istat
  type :: arr_ptr
    real(rp), pointer, contiguous , dimension(:,:,:) :: arr
  end type arr_ptr
  type(arr_ptr)    , allocatable, dimension(:) ::   io_vars
  integer          , allocatable, dimension(:,:) ::    nh_io
  character(len=10), allocatable, dimension(:) :: c_io_vars
  integer :: i,j,k,kk
  logical :: is_done,kill
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! initialize MPI/OpenMP
  !
  !$ call omp_set_num_threads(omp_get_max_threads())
  call initmpi(ng,dims,cbcpre,nh_g,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,nb,is_bound)
  twi = MPI_WTIME()
  savecounter = 0
  !
  ! allocate variables
  !
  allocate(u( 1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v), &
           v( 1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v), &
           w( 1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v), &
           p( 1-nh_p:n(1)+nh_p,1-nh_p:n(2)+nh_p,1-nh_p:n(3)+nh_p), &
           pp(1-nh_p:n(1)+nh_p,1-nh_p:n(2)+nh_p,1-nh_p:n(3)+nh_p))
#if defined(_CONSTANT_COEFFS_POISSON)
  allocate(pn,mold=pp)
#endif
  allocate(po,mold=pp)
  if(nscal > 0) allocate(s(1-nh_s:n(1)+nh_s,1-nh_s:n(2)+nh_s,1-nh_s:n(3)+nh_s))
  nio = 5+nscal
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  if(is_poisson_dtdma) then
#if defined(_OPENACC)
    n_z_d(:) = ap_z_dtdma%shape(:)
#else
    n_z_d(:) = dinfo_dtdma%zsz(:)
#endif
    allocate(ap_d(n_z_d(1),n_z_d(2),n_z_d(3)), &
             cp_d(n_z_d(1),n_z_d(2),n_z_d(3)))
  end if
  allocate(dzc( 1-nh_g:n(3)+nh_g), &
           dzf( 1-nh_g:n(3)+nh_g), &
           zc(  1-nh_g:n(3)+nh_g), &
           zf(  1-nh_g:n(3)+nh_g), &
           dzci(1-nh_g:n(3)+nh_g), &
           dzfi(1-nh_g:n(3)+nh_g))
  allocate(dzc_g( 1-nh_g:ng(3)+nh_g), &
           dzf_g( 1-nh_g:ng(3)+nh_g), &
           zc_g(  1-nh_g:ng(3)+nh_g), &
           zf_g(  1-nh_g:ng(3)+nh_g), &
           dzci_g(1-nh_g:ng(3)+nh_g), &
           dzfi_g(1-nh_g:ng(3)+nh_g))
  allocate(xc_g(1-nh_g:ng(1)+nh_g),xf_g(1-nh_g:ng(1)+nh_g), &
           yc_g(1-nh_g:ng(2)+nh_g),yf_g(1-nh_g:ng(2)+nh_g))
  allocate(grid_vol_ratio_c,mold=dzc)
  allocate(grid_vol_ratio_f,mold=dzf)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
  allocate(psi(1-nh_psi:n(1)+nh_psi,1-nh_psi:n(2)+nh_psi,1-nh_psi:n(3)+nh_psi))
  allocate(kappa,normx,normy,normz,mold=psi)
  allocate(psio,mold=psi)
#if defined(_SDF_NORMALS)
  allocate(phi,mold=psi)
#endif
  allocate(psiflx_x,psiflx_y,psiflx_z,mold=psi)
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
  allocate(surfx_n(n(1),n(2),n(3)),surfy_n(n(1),n(2),n(3)),surfz_n(n(1),n(2),n(3)), &
           surfx_o(n(1),n(2),n(3)),surfy_o(n(1),n(2),n(3)),surfz_o(n(1),n(2),n(3)))
  surfx_n(:,:,:) = 0._rp
  surfy_n(:,:,:) = 0._rp
  surfz_n(:,:,:) = 0._rp
  surfx_o(:,:,:) = 0._rp
  surfy_o(:,:,:) = 0._rp
  surfz_o(:,:,:) = 0._rp
#endif
  !
  if(is_debug) then
    if(myid == 0) print*, '*** CaNS build information ***'
    if(myid == 0) print*, 'Compiler version: ', compiler_version()
    if(myid == 0) print*, 'Compiler options: ', compiler_options()
    block
      character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version
      integer :: ilen
      call MPI_GET_LIBRARY_VERSION(mpi_version,ilen,ierr)
      if(myid == 0) print*, 'MPI version: ', trim(mpi_version)
    end block
    if(myid == 0) print*, ''
  end if
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, '*** Beginning of simulation ***'
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, ''
  call initgrid(gtype,ng(3),nh_g,gr,l(3),dzc_g,dzf_g,zc_g,zf_g,cbcpre(0,3)//cbcpre(1,3) == 'PP')
  do kk=1-nh_g,ng(1)+nh_g
    xc_g(kk) = (kk-0.5_rp)*dl(1)
    xf_g(kk) = kk*dl(1)
  end do
  do kk=1-nh_g,ng(2)+nh_g
    yc_g(kk) = (kk-0.5_rp)*dl(2)
    yf_g(kk) = kk*dl(2)
  end do
  if(myid == 0) then
    open(newunit=iunit,file=trim(datadir)//'geometry.out',status='replace')
    write(iunit,'(*(1x,i0))') ng(1),ng(2),ng(3)
    write(iunit,'(*(1x,es24.16e3))') l(1),l(2),l(3)
    close(iunit)
    call save_grid(datadir,'grid',ng,nh_g,zc_g,zf_g,dzc_g,dzf_g)
  end if
  !$acc enter data copyin(lo,hi,n) async
  !$acc enter data copyin(bforce,dl,dli,l) async
  !$acc enter data copyin(gacc) async
  !$acc enter data copyin(rho12,mu12,ka12,cp12,beta12) async
  !$acc enter data copyin(zc_g,zf_g,dzc_g,dzf_g) async
  !$acc enter data create(zc,zf,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g) async
  !
  !$acc parallel loop default(present) private(k) async
  do kk=lo(3)-nh_g,hi(3)+nh_g
    k = kk-(lo(3)-1)
    zc( k) = zc_g(kk)
    zf( k) = zf_g(kk)
    dzc(k) = dzc_g(kk)
    dzf(k) = dzf_g(kk)
    dzci(k) = dzc(k)**(-1)
    dzfi(k) = dzf(k)**(-1)
  end do
  !$acc data copy(ng) async
  !$acc parallel loop default(present) async
  do k=1-nh_g,ng(3)+nh_g
    dzci_g(k) = dzc_g(k)**(-1)
    dzfi_g(k) = dzf_g(k)**(-1)
  end do
  !$acc end data
  !$acc enter data create(grid_vol_ratio_c,grid_vol_ratio_f) async
  !$acc parallel loop default(present) async
  do k=1-nh_g,n(3)+nh_g
    grid_vol_ratio_c(k) = dl(1)*dl(2)*dzc(k)/(l(1)*l(2)*l(3))
    grid_vol_ratio_f(k) = dl(1)*dl(2)*dzf(k)/(l(1)*l(2)*l(3))
  end do
  !$acc update self(zc,zf,dzc,dzf,dzci,dzfi) async
  !$acc exit data copyout(zc_g,zf_g,dzc_g,dzf_g,dzci_g,dzfi_g) async ! not needed on the device
  !$acc wait
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_forced)
  !
  ! initialize Poisson solver
  !
  call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
                  lambdaxyp,['c','c','c'],ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
  !$acc enter data copyin(lambdaxyp,ap,bp,cp) async
  !$acc enter data copyin(rhsbp,rhsbp%x,rhsbp%y,rhsbp%z) async
  if(is_poisson_dtdma) then
    !$acc enter data create(ap_d,cp_d) async
  end if
  !$acc wait
#if defined(_OPENACC)
  !
  ! determine workspace sizes and allocate the memory
  !
  call init_wspace_arrays()
  call set_cufft_wspace(pack(arrplanp,.true.),istream_acc_queue_1)
  if(myid == 0) print*,'*** Device memory footprint (Gb): ', &
                  device_memory_footprint(n,nh_v,nh_p,nh_s,nh_g,n_z,nscal)/(1._sp*1024**3), ' ***'
#endif
  if(is_debug_poisson) then
    call test_sanity_solver(ng,lo,hi,n,nh_v,nh_p,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,n_z,l,dl,dli, &
                            dzc,dzf,dzci,dzfi,dzci_g,dzfi_g,nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
  end if
  !
  call acdi_set_epsilon(dl,nh_g,dzfi,psi_thickness_factor,seps) ! interface width, used by the initial conditions of all methods
  !
  allocate(io_vars(nio),c_io_vars(nio),nh_io(3,nio)) ! u,v,w,p,psi,s
  nh_io(:,1:3) = nh_v; nh_io(:,4) = nh_p; nh_io(:,5) = nh_psi
  io_vars(1)%arr => u  ; c_io_vars(1) = 'u'
  io_vars(2)%arr => v  ; c_io_vars(2) = 'v'
  io_vars(3)%arr => w  ; c_io_vars(3) = 'w'
  io_vars(4)%arr => p  ; c_io_vars(4) = 'p'
  io_vars(5)%arr => psi; c_io_vars(5) = 'psi'
  if(nscal > 0) then
    nh_io(:,6) = nh_s
    io_vars(6)%arr => s; c_io_vars(6) = 's'
  end if
  !
  is_dtdma_update_p = .true.
  !
  u(:,:,:)  = 0.
  v(:,:,:)  = 0.
  w(:,:,:)  = 0.
  p(:,:,:)  = 0.
  pp(:,:,:) = 0.
  po(:,:,:) = 0.
#if defined(_CONSTANT_COEFFS_POISSON)
  pn(:,:,:) = 0.
#endif
  if(nscal > 0) s(:,:,:) = 0.
  psi(:,:,:) = 0.
  if(.not.restart) then
    istep = 0
    time = 0.
    call initflow(inivel,cbcvel,bcvel,ng,lo,nh_v,nh_p,nh_g,l,dl,zc,zf,dzc,dzf,rho12(2),mu12(2), &
                  bforce,is_wallturb,time,u,v,w,p)
    if(nscal > 0) call initscal(inisca,bcsca,ng,lo,nh_s,nh_g,l,dl,zc,dzf,s)
    call init2fl(inipsi,cbcpsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,psi)
    if(myid == 0) print*, '*** Initial condition successfully set ***'
  else
    do is=1,nio
      call load_one('r',trim(datadir)//'fld_'//trim(c_io_vars(is))//trim(io_ext), &
                    MPI_COMM_WORLD,ng,nh_io(:,is),nh_g,lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                    is_check_time=is > 1)
    end do
    if(myid == 0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    if(is_overwrite_save.and.nsaves_max > 0) then
      if(myid == 0) then
        istat = read_link(trim(datadir)//'fld_'//trim(c_io_vars(1))//trim(io_ext),fname_save)
        if(istat == 0) read(fname_save(5:8),'(i4)',iostat=istat) savecounter
        if(istat /= 0) savecounter = 0
      end if
      call MPI_BCAST(savecounter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    end if
  end if
  !$acc enter data copyin(u,v,w,p,pp,pn,po)
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
  !$acc enter data copyin(surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o)
#endif
  call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,p)
#if defined(_CONSTANT_COEFFS_POISSON)
  !$acc parallel loop collapse(3) default(present) async(1)
  !$OMP parallel do   collapse(3) DEFAULT(shared)
  do k=1-nh_p,n(3)+nh_p
    do j=1-nh_p,n(2)+nh_p
      do i=1-nh_p,n(1)+nh_p
        pn(i,j,k) = p(i,j,k)
        po(i,j,k) = p(i,j,k)
      end do
    end do
  end do
#endif
  if(nscal > 0) then
    !$acc enter data copyin(s) async(1)
    call boundp(cbcsca,n,nh_s,nh_g,bcsca,nb,is_bound,dl,dzc,s)
  end if
  !$acc enter data copyin(psi) create(kappa,normx,normy,normz) async(1)
#if defined(_SDF_NORMALS)
  !$acc enter data create(phi) async(1)
#endif
  !$acc enter data create(psio,psiflx_x,psiflx_y,psiflx_z) async(1)
  call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
  !$acc wait
  !
#if defined(_SDF_NORMALS)
#if !defined(_INTERFACE_CAPTURING_VOF)
  call acdi_cmpt_phi(n,nh_psi,seps,psi,phi)
#else
  call vof_thinc_cmpt_phi(n,nh_psi,vof_thinc_beta,psi,phi)
#endif
  call cmpt_norm_curv(n,nh_psi,nh_g,dli,dzci,dzfi,phi,normx,normy,normz,kappa)
#else
  call cmpt_norm_curv(n,nh_psi,nh_g,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
#endif
  call boundp(cbccur,n,nh_psi,nh_g,bccur,nb,is_bound,dl,dzc,kappa)
  call boundp(cbcnor(:,:,1),n,nh_psi,nh_g,bcnor(:,:,1),nb,is_bound,dl,dzc,normx)
  call boundp(cbcnor(:,:,2),n,nh_psi,nh_g,bcnor(:,:,2),nb,is_bound,dl,dzc,normy)
  call boundp(cbcnor(:,:,3),n,nh_psi,nh_g,bcnor(:,:,3),nb,is_bound,dl,dzc,normz)
  !
#if !defined(_INTERFACE_CAPTURING_VOF)
  call acdi_set_gamma(n,nh_v,acdi_gam_factor,u,v,w,gam)
  gam = max(gam,acdi_gam_min)
  if(myid == 0) print*, 'ACDI parameters. Gamma: ', gam, 'Epsilon: ', seps
#endif
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !$acc wait ! not needed but to prevent possible future issues
  !$acc update self(u,v,w,p,psi,kappa)
  if(nscal > 0) then
    !$acc update self(s)
  end if
  if(is_output_due(istep,iout1d,tout1d,time,time,restart)) then
    include 'out1d.h90'
  end if
  if(is_output_due(istep,iout2d,tout2d,time,time,restart)) then
    include 'out2d.h90'
  end if
  if(is_output_due(istep,iout3d,tout3d,time,time,restart)) then
    include 'out3d.h90'
  end if
  !
  call chkdt(n,nh_v,nh_g,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dt_cfl,gam,seps,ka12,cp12)
  dt = min(cfl*dt_cfl,dtmax)
  if(dt_f > 0.) then
    if(dt_f > dt) then
      if(myid == 0) print*, 'WARNING: fixed time step exceeds estimated stability limit.'
      if(myid == 0) print*, 'dt_f = ', dt_f, 'Estimated stable dt = ', dt
    end if
    dt = dt_f
  end if
  if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
  dto = dt
  dti = 1./dt
  kill = .false.
  force_balance(:,:) = 0._rp
  force_balance(:,5) = bforce(:)
  !
  ! main loop
  !
  if(myid == 0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)
    if(is_timing) then
      !$acc wait(1)
      dt12 = MPI_WTIME()
    end if
    istep = istep + 1
    time_start = time
    time = time + dt
    if(myid == 0) print*, 'Time step #', istep, 'Time = ', time
    dt_r = dt/dto
    tm_coeff(1) = 1. + 0.5*dt_r
    tm_coeff(2) =    - 0.5*dt_r
    dtrk = dt
    dtrki = dtrk**(-1)
    !
    ! phase field update
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1-nh_psi,n(3)+nh_psi
      do j=1-nh_psi,n(2)+nh_psi
        do i=1-nh_psi,n(1)+nh_psi
          psio(i,j,k) = psi(i,j,k)
        end do
      end do
    end do
    if(is_track_interface) then
      call tm_2fl(tm_coeff,n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,dt,gam,seps,vof_thinc_beta,u,v,w,normx,normy,normz,phi,psi, &
                  psiflx_x,psiflx_y,psiflx_z)
      call bounduvw(cbcvel,n,nh_psi,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
      call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
#if defined(_SDF_NORMALS)
#if !defined(_INTERFACE_CAPTURING_VOF)
      call acdi_cmpt_phi(n,nh_psi,seps,psi,phi)
#else
      call vof_thinc_cmpt_phi(n,nh_psi,vof_thinc_beta,psi,phi)
#endif
      call cmpt_norm_curv(n,nh_psi,nh_g,dli,dzci,dzfi,phi,normx,normy,normz,kappa)
#else
      call cmpt_norm_curv(n,nh_psi,nh_g,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
#endif
      call boundp(cbccur,n,nh_psi,nh_g,bccur,nb,is_bound,dl,dzc,kappa)
      call boundp(cbcnor(:,:,1),n,nh_psi,nh_g,bcnor(:,:,1),nb,is_bound,dl,dzc,normx)
      call boundp(cbcnor(:,:,2),n,nh_psi,nh_g,bcnor(:,:,2),nb,is_bound,dl,dzc,normy)
      call boundp(cbcnor(:,:,3),n,nh_psi,nh_g,bcnor(:,:,3),nb,is_bound,dl,dzc,normz)
    else
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1-nh_psi,n(3)+nh_psi
        do j=1-nh_psi,n(2)+nh_psi
          do i=1-nh_psi,n(1)+nh_psi
            psiflx_x(i,j,k) = 0.
            psiflx_y(i,j,k) = 0.
            psiflx_z(i,j,k) = 0.
          end do
        end do
      end do
    end if
    if(nscal > 0) then
      call tm_scal(tm_coeff,n,nh_v,nh_s,nh_psi,nh_g,dli,dzci,dzfi,dt,ssource,rho12,ka12,cp12,psi,u,v,w, &
                   psio,psiflx_x,psiflx_y,psiflx_z,s)
      call boundp(cbcsca,n,nh_s,nh_g,bcsca,nb,is_bound,dl,dzc,s)
    end if
    if(.not.is_solve_ns) then
      call initflow(inivel,cbcvel,bcvel,ng,lo,nh_v,nh_p,nh_g,l,dl,zc,zf,dzc,dzf,rho12(2),mu12(2), &
                    bforce,is_wallturb,time,u,v,w,p)
      !$acc wait(1)
      !$acc update device(u,v,w,p) async(1)
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w)
    else
      rho_av(:) = 0.
      is_cmpt_rho_av(:) = (abs(gacc(:)) > 0.) .and. cbcpre(0,:)//cbcpre(1,:) == 'PP'
      if(is_cmpt_rho_av(1)) &
        call bulk_mean_12_stag(n,nh_psi,nh_g,1,grid_vol_ratio_f,psi,rho12,rho_av(1))
      if(is_cmpt_rho_av(2)) &
        call bulk_mean_12_stag(n,nh_psi,nh_g,2,grid_vol_ratio_f,psi,rho12,rho_av(2))
      if(is_cmpt_rho_av(3)) &
        call bulk_mean_12_stag(n,nh_psi,nh_g,3,grid_vol_ratio_c,psi,rho12,rho_av(3))
      call tm(tm_coeff,n,nh_v,nh_p,nh_s,nh_psi,nh_g,is_bound,is_forced,l,dli,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f, &
              dt,dt_r,bforce,gacc,sigma,rho_av,rho12,mu12,beta12,rho0,psi,kappa,p,pn,po,s, &
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
              surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o, &
#endif
              psio,psiflx_x,psiflx_y,psiflx_z,u,v,w,force_balance)
      if(is_forced_hit) then
        call lscale_forcing(2,lo,hi,nh_v,nh_g,0.5_rp,dtrk,l,dl,zc,zf,u,v,w)
      end if
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w,.false.,-1)
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1-nh_p,n(3)+nh_p
        do j=1-nh_p,n(2)+nh_p
          do i=1-nh_p,n(1)+nh_p
            pp(i,j,k) = p(i,j,k)
          end do
        end do
      end do
      call fillps(n,nh_v,nh_p,nh_g,dli,dzfi,dtrki*rho0,u,v,w,p)
#if defined(_CONSTANT_COEFFS_POISSON)
      call updt_rhs_b(['c','c','c'],cbcpre,n,nh_p,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,p)
      call solver(n,nh_p,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],p,is_dtdma_update_p,ap_d,cp_d, &
                  .not.is_impdiff)
#else
      call solver_vc(ng,lo,hi,nh_p,nh_psi,nh_g,cbcpre,bcpre,dli,dzci,dzfi,is_bound,rho12,psi,p,po)
#endif
      call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,p)
#if defined(_CONSTANT_COEFFS_POISSON)
      call correc(n,nh_v,nh_p,nh_psi,nh_g,dli,dzci,rho12,dtrk/rho0,p,psi,u,v,w)
#else
      call correc(n,nh_v,nh_p,nh_psi,nh_g,dli,dzci,rho12,dtrk,p,psi,u,v,w)
#endif
      call bounduvw(cbcvel,n,nh_v,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,u,v,w,.true.,1,0)
      if(any(is_forced(:))) then
        call cmpt_momentum(n,nh_v,nh_psi,nh_g,is_forced,grid_vol_ratio_c,grid_vol_ratio_f,rho12,psi,u,v,w,momentum)
        where(is_forced(:)) force_balance(:,6) = (momentum(:)-force_balance(:,6))/dt
      end if
      call updatep(n,nh_p,nh_g,dli,dzci,dzfi,alpha,pp,p)
      call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,p)
    end if
#if defined(_CONSTANT_COEFFS_POISSON)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1-nh_p,n(3)+nh_p
      do j=1-nh_p,n(2)+nh_p
        do i=1-nh_p,n(1)+nh_p
          po(i,j,k) = pn(i,j,k)
          pn(i,j,k) =  p(i,j,k)
        end do
      end do
    end do
#endif
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep >= nstep   ) is_done = is_done.or..true.
    end if
    if(stop_type(2)) then ! maximum simulation time reached
      if(time  >= time_max) is_done = is_done.or..true.
    end if
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600.
      call MPI_ALLREDUCE(MPI_IN_PLACE,tw,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    end if
    dto = dt
#if !defined(_INTERFACE_CAPTURING_VOF)
    call acdi_set_gamma(n,nh_v,acdi_gam_factor,u,v,w,gam)
    gam = max(gam,acdi_gam_min)
#endif
    if(icheck > 0.and.mod(istep,max(icheck,1)) == 0) then
      if(myid == 0) print*, 'Checking stability and divergence...'
      call chkdt(n,nh_v,nh_g,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dt_cfl,gam,seps,ka12,cp12)
      dt = min(cfl*dt_cfl,dtmax)
      if(dt_f > 0.) then
        if(dt_f > dt) then
          if(myid == 0) print*, 'ERROR: fixed time step exceeds estimated stability limit.'
          if(myid == 0) print*, 'dt_f = ', dt_f, 'Estimated stable dt = ', dt
          if(myid == 0) print*, 'Aborting...'
          is_done = .true.
          kill = .true.
        else
          dt = dt_f
        end if
      end if
      if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
      dti = 1./dt
      call chkdiv(lo,hi,nh_v,nh_g,l,dli,dzfi,u,v,w,divtot,divmax)
      if(myid == 0) print*, 'Velocity divergence norms: Mean = ', divtot, '| Maximum = ', divmax
      if((divmax > small.and..not.is_mask_divergence_check).or.is_nan(divtot)) then
        if(myid == 0) print*, 'ERROR: maximum divergence is too large.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
    end if
    !
    ! output routines below
    !
    if(is_output_due(istep,iout0d,tout0d,time_start,time)) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      dpdl(:) = 0._rp
      where(is_forced(:)) dpdl(:) = -(force_balance(:,5)-bforce(:))
      if(any(abs(bforce(:)) > 0.).or.any(abs(dpdl(:)) > 0.)) then
        meanvel(:) = 0._rp
        if(abs(bforce(1)) > 0..or.abs(dpdl(1)) > 0.) then
          call bulk_mean(n,nh_v,nh_g,grid_vol_ratio_f,u,meanvel(1))
        end if
        if(abs(bforce(2)) > 0..or.abs(dpdl(2)) > 0.) then
          call bulk_mean(n,nh_v,nh_g,grid_vol_ratio_f,v,meanvel(2))
        end if
        if(abs(bforce(3)) > 0..or.abs(dpdl(3)) > 0.) then
          call bulk_mean(n,nh_v,nh_g,grid_vol_ratio_c,w,meanvel(3))
        end if
        var(1)   = time
        var(2:4) = bforce(1:3)
        var(5:7) = dpdl(1:3)
        var(8:10) = meanvel(1:3)
        call out0d(trim(datadir)//'forcing.out',10,var)
      end if
      if(any(is_forced(:))) then
        var(1) = time
        var(2:4) = force_balance(:,1)
        var(5:7) = force_balance(:,2)
        var(8:10) = force_balance(:,3)
        var(11:13) = force_balance(:,4)
        var(14:16) = force_balance(:,5)
        var(17:19) = force_balance(:,6)
        call out0d(trim(datadir)//'forcing_balance.out',19,var)
      end if
      !
#if !defined(_INTERFACE_CAPTURING_VOF)
      var(1) = 1.*istep
      var(2) = time
      var(3) = gam
      var(4) = seps
      call out0d(trim(datadir)//'log_acdi.out',4,var)
#endif
    end if
    write(fldnum,'(i7.7)') istep
    if(is_output_due(istep,iout1d,tout1d,time_start,time)) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
      if(nscal > 0) then
        !$acc update self(s)
      end if
      include 'out1d.h90'
    end if
    if(is_output_due(istep,iout2d,tout2d,time_start,time)) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
      if(nscal > 0) then
        !$acc update self(s)
      end if
      include 'out2d.h90'
    end if
    if(is_output_due(istep,iout3d,tout3d,time_start,time)) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
      if(nscal > 0) then
        !$acc update self(s)
      end if
      include 'out3d.h90'
    end if
    if((isave > 0.or.tsave > 0.).and.(.not.kill).and.(is_output_due(istep,isave,tsave,time_start,time).or.is_done)) then
      if(is_overwrite_save.and.nsaves_max <= 0) then
        filename = 'fld'
      else
        filename = 'fld_'//fldnum
        if(is_overwrite_save) then
          if(savecounter >= nsaves_max) savecounter = 0
          savecounter = savecounter + 1
          write(chkptnum,'(i4.4)') savecounter
          filename = 'fld_'//chkptnum
        end if
      end if
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
      if(nscal > 0) then
        !$acc update self(s)
      end if
      is_save_tmp = trim(io_ext) /= '.bp'
      do is=1,nio
        fname_save = trim(datadir)//trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext)
        if(is_save_tmp) fname_save = trim(datadir)//'fld_'//trim(c_io_vars(is))//'.tmp'//trim(io_ext)
        select case(trim(c_io_vars(is)))
        case('u')
          call load_one('w',trim(fname_save), &
                        MPI_COMM_WORLD,ng,nh_io(:,is),nh_g,lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        xf_g,yc_g,zc_g)
        case('v')
          call load_one('w',trim(fname_save), &
                        MPI_COMM_WORLD,ng,nh_io(:,is),nh_g,lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        xc_g,yf_g,zc_g)
        case('w')
          call load_one('w',trim(fname_save), &
                        MPI_COMM_WORLD,ng,nh_io(:,is),nh_g,lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        xc_g,yc_g,zf_g)
        case default
          call load_one('w',trim(fname_save), &
                        MPI_COMM_WORLD,ng,nh_io(:,is),nh_g,lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        xc_g,yc_g,zc_g)
        end select
      end do
      if(is_save_tmp) then
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid == 0) then
          do is=1,nio
            istat = rename_file(trim(datadir)//'fld_'//trim(c_io_vars(is))//'.tmp'//trim(io_ext), &
                                trim(datadir)//trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext))
            if(istat /= 0) print*, 'WARNING: could not rename the temporary checkpoint file of ', trim(c_io_vars(is)), '.'
          end do
        end if
      end if
      if(.not.is_overwrite_save.or.nsaves_max > 0) then
        !
        ! fld_*.<ext> -> last checkpoint file (symbolic link)
        !
        do is=1,nio
          call gen_alias(myid,trim(datadir),trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext), &
                         'fld_'//trim(c_io_vars(is))//trim(io_ext))
        end do
      end if
      if(is_overwrite_save.and.nsaves_max > 0) then
        var(1) = 1.*istep
        var(2) = time
        var(3) = 1.*savecounter
        call out0d(trim(datadir)//'log_checkpoints.out',3,var)
      end if
      if(myid == 0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    end if
    if(is_timing) then
      !$acc wait(1)
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) print*, 'Average, minimum & maximum elapsed time: '
      if(myid == 0) print*, dt12av/(1.*product(dims)),dt12min,dt12max
    end if
  end do
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
#if defined(_OPENACC)
  call cudecomp_finalize
#endif
  call MPI_FINALIZE(ierr)
  if(kill) error stop
end program cans
