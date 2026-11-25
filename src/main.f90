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
  use mod_bound              , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv             , only: chkdiv
  use mod_chkdt              , only: chkdt
  use mod_common_mpi         , only: myid,ierr,dinfo_dtdma
  use mod_correc             , only: correc
  use mod_fft                , only: fftini,fftend
  use mod_fillps             , only: fillps
  use mod_forcing            , only: lscale_forcing
  use mod_initflow           , only: initflow,initscal,set_bulk_velocity
  use mod_initgrid           , only: initgrid,save_grid
  use mod_initmpi            , only: initmpi
  use mod_impdiff            , only: impdiff_z_mom_solve,impdiff_z_scal_solve,impdiff_updt_forcing,impdiff_updt_forcing_s, &
                                     impdiff_cc_mom_updt,impdiff_cc_scal_updt
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mod_impdiff            , only: impdiff_vc_mom_updt,impdiff_vc_mom_solve,impdiff_vc_scal_updt,impdiff_vc_scal_solve
#endif
  use mod_initsolver         , only: initsolver
  use mod_interface_capturing, only: ifc_init_params,ifc_init_field,ifc_init,ifc_advance,ifc_advecting_velocity, &
                                     ifc_update_params,ifc_write_log, &
                                     marker,ifc_markers_init,ifc_markers_advance,ifc_markers_force
  use mod_load               , only: load_one
  use mod_mom                , only: cmpt_momentum,marker_fields
  use mod_rk                 , only: tm => rk,tm_scal => rk_scal
  use mod_output             , only: is_output_due, &
                                     out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  use mod_param              , only: ng,l,dl,dli, &
                                     gtype,gr, &
                                     cfl,dtmax,dt_f, &
                                     is_solve_ns,is_track_interface, &
                                     inivel,iniscal,inipsi, &
                                     is_wallturb,is_forced_hit,forced_hit_type,forced_hit_amp, &
                                     nstep,time_max,tw_max,stop_type, &
                                     restart,is_overwrite_save,nsaves_max, &
                                     icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                     tout0d,tout1d,tout2d,tout3d,tsave, &
                                     cbcvel,bcvel,cbcvel_halo,bcvel_halo,cbcpre,bcpre, &
                                     cbcscal,bcscal, &
                                     is_forced,velf,bforce,ssource,is_sforced,scalf,cbcpsi,bcpsi, &
                                     gacc,nscal,nmarker, &
                                     rho0,visc0,alpha0,rho12,mu12,sigma,ka12,cp12,beta12, &
                                     psi_thickness_factor, &
                                     dims, &
                                     nb,is_bound, &
                                     small, &
                                     datadir,io_ext,is_full_restart, &
                                     read_input, &
                                     is_debug,is_debug_poisson, &
                                     is_timing, &
                                     is_impdiff,is_impdiff_1d,is_ifc_split,is_ifc_plic,is_sdf_normals,is_vof_sweep_transport, &
                                     is_poisson_dtdma, &
                                     is_mask_divergence_check, &
                                     nh_v,nh_p,nh_s,nh_psi,nh_ifc,nh_g
  use mod_sanity             , only: test_sanity_input,test_sanity_solver
  use mod_scal               , only: scalar,cmpt_scal_mean
  use mod_solve_helmholtz    , only: solve_helmholtz,rhs_bound
  use mod_two_fluid          , only: bulk_mean_12_stag
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mod_solver_vc          , only: solver_vc,solver_vc_init,solver_vc_end
#endif
#if !defined(_OPENACC)
  use mod_solver             , only: solver
#else
  use mod_solver_gpu         , only: solver => solver_gpu
  use mod_workspaces         , only: init_wspace_arrays,set_cufft_wspace,cudecomp_finalize
  use mod_common_cudecomp    , only: istream_acc_queue_1,ap_z_dtdma
#endif
  use mod_updatep            , only: updatep
  use mod_posix              , only: rename_file,read_link
  use mod_utils              , only: bulk_mean
#if defined(_OPENACC)
  use mod_utils              , only: device_memory_footprint
#endif
  use mod_types
  use omp_lib
  implicit none
  integer , dimension(3) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
  real(rp), target, allocatable, dimension(:,:,:) :: u,v,w,p,pp
  real(rp), target, allocatable, dimension(:,:,:) :: po
  real(rp), target, allocatable, dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko,dpsidtrko
  real(rp), target, allocatable, dimension(:,:,:) :: uo,vo,wo ! previous velocity (direction-split advection)
  real(rp), dimension(3) :: wall_force_old
  logical , dimension(3) :: is_velf ! direction whose initial mean velocity is set from `velf` (see `initflow`)
  logical :: is_first
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
#if !defined(_OPENACC) || defined(_USE_HIP)
  type(C_PTR), dimension(2,2) :: arrplanu,arrplanv,arrplanw
#else
  integer    , dimension(2,2) :: arrplanu,arrplanv,arrplanw
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyu,lambdaxyv,lambdaxyw
  real(rp), allocatable, dimension(:) :: au,av,aw,bu,bv,bw,cu,cv,cw
  real(rp) :: normfftu,normfftv,normfftw
  type(rhs_bound) :: rhsbu,rhsbv,rhsbw
  logical :: is_impdiff_cc
  !
  real(rp) :: dt,dto,dt_r,dt_cfl,time,time_start,dtrk,dtrki,divtot,divmax
  real(rp), dimension(2) :: dt_restart
  real(rp) :: gam,seps
  integer :: istep,iunit
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g, &
                                         xc_g,xf_g,yc_g,yf_g, &
                                         grid_vol_ratio_c,grid_vol_ratio_f,wzc
  real(rp), dimension(3) :: dpdl,meanvel,momentum,momentum_o,wall_force_imp
  !
  ! zero-net-heating forcing of the scalars (the counterpart of the zero-net-acceleration forcing): heat lost through
  ! the boundaries in the previous step, balance of each scalar, and the part of the implicit diffusion
  !
  real(rp), allocatable, dimension(:)   :: wall_heat_old
  real(rp), allocatable, dimension(:,:) :: sforce_balance
  real(rp) :: wall_heat_imp,stot,stot_o,rhocp_mean
  real(rp), dimension(3,7) :: force_balance
  real(rp), dimension(2) :: tm_coeff
  !
  ! scalar fields (see `mod_scal`)
  !
  type(scalar), target, allocatable, dimension(:) :: scalars
  type(scalar), pointer :: s
  integer :: iscal
  character(len=3) :: scalnum
  integer :: is,nio,nio_s
  !
  ! two-fluid solver specific
  !
  real(rp), target, allocatable, dimension(:,:,:) :: psi
  real(rp), allocatable, dimension(:,:,:) :: psio,phi,kappa,normx,normy,normz, &
                                             psiflx_x,psiflx_y,psiflx_z
  !
  ! several markers (see `mod_interface_capturing`): psi is then the total phase indicator
  !
  type(marker), target, allocatable, dimension(:) :: markers
  type(marker_fields) , allocatable, dimension(:) :: mfields
  real(rp), allocatable, dimension(:,:,:) :: fsx,fsy,fsz ! surface tension force summed over the markers
  integer :: imark,nmk
  character(len=3) :: marknum
  real(rp), target, allocatable, dimension(:,:,:) :: uadv,vadv,wadv ! advecting (or transported) velocity of the geometric
                                                                    ! VoF method
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
  real(rp), target, allocatable, dimension(:,:,:) :: surfx_n,surfy_n,surfz_n, &
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
  character(len=16), allocatable, dimension(:) :: c_io_vars
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
#if !defined(_CONSTANT_COEFFS_POISSON)
  call solver_vc_init
#endif
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
  allocate(po,mold=pp)
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  if(.not.is_ifc_split) allocate(dpsidtrko(n(1),n(2),n(3))) ! (only the AB2 transport of the phase indicator uses it)
  !
  ! scalars (one element is allocated if there is none: `scalars(1)%val` is then not allocated, and can be
  ! passed as an absent argument)
  !
  allocate(scalars(max(nscal,1)))
  do iscal=1,nscal
    s => scalars(iscal)
    allocate(s%val(1-nh_s:n(1)+nh_s,1-nh_s:n(2)+nh_s,1-nh_s:n(3)+nh_s),s%dsdtrko(n(1),n(2),n(3)))
  end do
  !$acc enter data copyin(scalars(:))
  nmk = merge(nmarker,0,nmarker > 1) ! (a single marker is the phase indicator itself)
  nio = 5+nscal+nmk
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
  allocate(wzc,mold=dzc)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
  !
  ! implicit diffusion in all directions: with constant coefficients and the FFT-based solver, as the
  ! pressure, or with variable coefficients and HYPRE, as the pressure
  !
#if defined(_CONSTANT_COEFFS_POISSON)
  is_impdiff_cc = is_impdiff .and. .not.is_impdiff_1d
#else
  is_impdiff_cc = .false.
#endif
  if(is_impdiff_cc) then
    allocate(lambdaxyu(n_z(1),n_z(2)), &
             lambdaxyv(n_z(1),n_z(2)), &
             lambdaxyw(n_z(1),n_z(2)))
    allocate(au(n_z(3)),bu(n_z(3)),cu(n_z(3)), &
             av(n_z(3)),bv(n_z(3)),cv(n_z(3)), &
             aw(n_z(3)),bw(n_z(3)),cw(n_z(3)))
    allocate(rhsbu%x(n(2),n(3),0:1), &
             rhsbu%y(n(1),n(3),0:1), &
             rhsbu%z(n(1),n(2),0:1), &
             rhsbv%x(n(2),n(3),0:1), &
             rhsbv%y(n(1),n(3),0:1), &
             rhsbv%z(n(1),n(2),0:1), &
             rhsbw%x(n(2),n(3),0:1), &
             rhsbw%y(n(1),n(3),0:1), &
             rhsbw%z(n(1),n(2),0:1))
    do iscal=1,nscal
      s => scalars(iscal)
      allocate(s%lambdaxy(n_z(1),n_z(2)))
      allocate(s%a(n_z(3)),s%b(n_z(3)),s%c(n_z(3)))
      allocate(s%rhsb%x(n(2),n(3),0:1), &
               s%rhsb%y(n(1),n(3),0:1), &
               s%rhsb%z(n(1),n(2),0:1))
    end do
  end if
  allocate(psi(1-nh_psi:n(1)+nh_psi,1-nh_psi:n(2)+nh_psi,1-nh_psi:n(3)+nh_psi))
  allocate(psio,mold=psi)
  allocate(kappa(1-nh_ifc:n(1)+nh_ifc,1-nh_ifc:n(2)+nh_ifc,1-nh_ifc:n(3)+nh_ifc))
  if(.not.is_ifc_plic) allocate(normx,normy,normz,mold=kappa) ! (these three are otherwise absent arguments of the
  if(is_sdf_normals) allocate(phi,mold=kappa)                 !  routines that take them)
  allocate(markers(nmk))
  do imark=1,nmk
    allocate(markers(imark)%psi,mold=psi)
    if(.not.is_ifc_split) allocate(markers(imark)%dpsidtrko(n(1),n(2),n(3)))
    allocate(markers(imark)%kappa,mold=kappa)
    if(.not.is_ifc_plic) allocate(markers(imark)%normx,markers(imark)%normy,markers(imark)%normz,mold=kappa)
    if(is_sdf_normals) allocate(markers(imark)%phi,mold=kappa)
  end do
  if(nmk > 0) then
    allocate(fsx(n(1),n(2),n(3)),fsy(n(1),n(2),n(3)),fsz(n(1),n(2),n(3)))
    allocate(mfields(nmk))
    do imark=1,nmk
      mfields(imark)%psi   => markers(imark)%psi
      mfields(imark)%kappa => markers(imark)%kappa
    end do
  end if
  allocate(psiflx_x,psiflx_y,psiflx_z,mold=kappa)
  if(is_ifc_split) then
    allocate(uadv,vadv,wadv,uo,vo,wo,mold=u)
    do iscal=1,nscal
      s => scalars(iscal)
      allocate(s%so,s%sadv,mold=s%val)
    end do
  end if
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
  !
  ! weight of cell `k` in the volume average of a cell-centered field at the face between cells `k` and `k+1`
  ! (the control volume of the vertical velocity), which is the average that conserves mass
  !
  !$acc enter data create(wzc) async
  !$acc parallel loop default(present) async
  do k=1-nh_g,n(3)+nh_g
    wzc(k) = 0.5
    if(k < n(3)+nh_g) wzc(k) = dzf(k)/(dzf(k)+dzf(k+1))
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
  if(is_impdiff_cc) then
    call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,1),bcvel(:,:,1), &
                    lambdaxyu,['f','c','c'],au,bu,cu,arrplanu,normfftu,rhsbu%x,rhsbu%y,rhsbu%z)
    call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,2),bcvel(:,:,2), &
                    lambdaxyv,['c','f','c'],av,bv,cv,arrplanv,normfftv,rhsbv%x,rhsbv%y,rhsbv%z)
    call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcvel(:,:,3),bcvel(:,:,3), &
                    lambdaxyw,['c','c','f'],aw,bw,cw,arrplanw,normfftw,rhsbw%x,rhsbw%y,rhsbw%z)
    !$acc enter data copyin(lambdaxyu,au,bu,cu,lambdaxyv,av,bv,cv,lambdaxyw,aw,bw,cw) async
    !$acc enter data copyin(rhsbu,rhsbu%x,rhsbu%y,rhsbu%z) async
    !$acc enter data copyin(rhsbv,rhsbv%x,rhsbv%y,rhsbv%z) async
    !$acc enter data copyin(rhsbw,rhsbw%x,rhsbw%y,rhsbw%z) async
    do iscal=1,nscal
      s => scalars(iscal)
      call initsolver(ng,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcscal(:,:,iscal),bcscal(:,:,iscal), &
                      s%lambdaxy,['c','c','c'],s%a,s%b,s%c,s%arrplan,s%normfft,s%rhsb%x,s%rhsb%y,s%rhsb%z)
      !$acc enter data copyin(s%lambdaxy,s%a,s%b,s%c) async
      !$acc enter data copyin(s%rhsb,s%rhsb%x,s%rhsb%y,s%rhsb%z) async
    end do
  end if
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
  if(is_impdiff_cc) then
    call set_cufft_wspace(pack(arrplanu,.true.),istream_acc_queue_1)
    call set_cufft_wspace(pack(arrplanv,.true.),istream_acc_queue_1)
    call set_cufft_wspace(pack(arrplanw,.true.),istream_acc_queue_1)
    do iscal=1,nscal
      call set_cufft_wspace(pack(scalars(iscal)%arrplan,.true.),istream_acc_queue_1)
    end do
  end if
  if(myid == 0) print*,'*** Device memory footprint (Gb): ', &
                  device_memory_footprint(n,nh_v,nh_p,nh_s,nh_g,n_z,nscal)/(1._sp*1024**3), ' ***'
#endif
  if(is_debug_poisson) then
    call test_sanity_solver(ng,lo,hi,n,nh_v,nh_p,nh_g,n_x_fft,n_y_fft,lo_z,hi_z,n_z,l,dl,dli, &
                            dzc,dzf,dzci,dzfi,dzci_g,dzfi_g,nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
  end if
  !
  call ifc_init_params(dl,nh_g,dzfi,psi_thickness_factor,seps)
  !
  ! checkpoint variables: the fields, followed by the terms of the previous time step needed to
  ! restart the time integration exactly, saved and loaded with `is_full_restart = T`
  !
  nio_s = nio
  if(is_full_restart) then
    nio_s = nio + 4 + nscal
    if(.not.is_ifc_split) nio_s = nio_s + 1 + nmk
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
    nio_s = nio_s + 6
#endif
    if(is_ifc_split) then
      nio_s = nio_s + 3 + nscal
      if(is_vof_sweep_transport) nio_s = nio_s + 3
    end if
  end if
  allocate(io_vars(nio_s),c_io_vars(nio_s),nh_io(3,nio_s)) ! u,v,w,p,psi,scalars,markers
  nh_io(:,1:3) = nh_v; nh_io(:,4) = nh_p; nh_io(:,5) = nh_psi
  io_vars(1)%arr => u  ; c_io_vars(1) = 'u'
  io_vars(2)%arr => v  ; c_io_vars(2) = 'v'
  io_vars(3)%arr => w  ; c_io_vars(3) = 'w'
  io_vars(4)%arr => p  ; c_io_vars(4) = 'p'
  io_vars(5)%arr => psi; c_io_vars(5) = 'psi'
  do iscal=1,nscal
    write(scalnum,'(i3.3)') iscal
    nh_io(:,5+iscal) = nh_s
    io_vars(5+iscal)%arr => scalars(iscal)%val; c_io_vars(5+iscal) = 's_'//scalnum
  end do
  do imark=1,nmk
    write(marknum,'(i3.3)') imark
    nh_io(:,5+nscal+imark) = nh_psi
    io_vars(5+nscal+imark)%arr => markers(imark)%psi; c_io_vars(5+nscal+imark) = 'psi_'//marknum
  end do
  if(is_full_restart) then
    is = nio
    nh_io(:,nio+1:) = 0; nh_io(:,nio+1) = nh_p
    io_vars(is+1)%arr => po       ; c_io_vars(is+1) = 'po'
    io_vars(is+2)%arr => dudtrko  ; c_io_vars(is+2) = 'dudtrko'
    io_vars(is+3)%arr => dvdtrko  ; c_io_vars(is+3) = 'dvdtrko'
    io_vars(is+4)%arr => dwdtrko  ; c_io_vars(is+4) = 'dwdtrko'
    is = is + 4
    if(.not.is_ifc_split) then
      io_vars(is+1)%arr => dpsidtrko; c_io_vars(is+1) = 'dpsidtrko'
      is = is + 1
    end if
    do iscal=1,nscal
      write(scalnum,'(i3.3)') iscal
      io_vars(is+1)%arr => scalars(iscal)%dsdtrko; c_io_vars(is+1) = 'dsdtrko_'//scalnum
      is = is + 1
    end do
    if(.not.is_ifc_split) then
      do imark=1,nmk
        write(marknum,'(i3.3)') imark
        io_vars(is+1)%arr => markers(imark)%dpsidtrko; c_io_vars(is+1) = 'dpsidtrko_'//marknum
        is = is + 1
      end do
    end if
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
    io_vars(is+1)%arr => surfx_n; c_io_vars(is+1) = 'surfx_n'
    io_vars(is+2)%arr => surfy_n; c_io_vars(is+2) = 'surfy_n'
    io_vars(is+3)%arr => surfz_n; c_io_vars(is+3) = 'surfz_n'
    io_vars(is+4)%arr => surfx_o; c_io_vars(is+4) = 'surfx_o'
    io_vars(is+5)%arr => surfy_o; c_io_vars(is+5) = 'surfy_o'
    io_vars(is+6)%arr => surfz_o; c_io_vars(is+6) = 'surfz_o'
    is = is + 6
#endif
    if(is_ifc_split) then
      nh_io(:,is+1:is+3) = nh_v
      io_vars(is+1)%arr => uo; c_io_vars(is+1) = 'uo'
      io_vars(is+2)%arr => vo; c_io_vars(is+2) = 'vo'
      io_vars(is+3)%arr => wo; c_io_vars(is+3) = 'wo'
      is = is + 3
      if(is_vof_sweep_transport) then
        !
        ! (the transported velocity and scalars of the previous time step replace the previous scalars)
        !
        nh_io(:,is+1:is+3) = nh_v
        io_vars(is+1)%arr => uadv; c_io_vars(is+1) = 'uadv'
        io_vars(is+2)%arr => vadv; c_io_vars(is+2) = 'vadv'
        io_vars(is+3)%arr => wadv; c_io_vars(is+3) = 'wadv'
        is = is + 3
        do iscal=1,nscal
          write(scalnum,'(i3.3)') iscal
          nh_io(:,is+1) = nh_s
          io_vars(is+1)%arr => scalars(iscal)%sadv; c_io_vars(is+1) = 'sadv_'//scalnum
          is = is + 1
        end do
      else
        do iscal=1,nscal
          write(scalnum,'(i3.3)') iscal
          nh_io(:,is+1) = nh_s
          io_vars(is+1)%arr => scalars(iscal)%so; c_io_vars(is+1) = 'so_'//scalnum
          is = is + 1
        end do
      end if
    end if
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
  dudtrko(:,:,:) = 0.
  dvdtrko(:,:,:) = 0.
  dwdtrko(:,:,:) = 0.
  if(.not.is_ifc_split) dpsidtrko(:,:,:) = 0.
  do iscal=1,nscal
    scalars(iscal)%dsdtrko(:,:,:) = 0.
  end do
  if(is_ifc_split) then
    uo(:,:,:) = 0.
    vo(:,:,:) = 0.
    wo(:,:,:) = 0.
    uadv(:,:,:) = 0.
    vadv(:,:,:) = 0.
    wadv(:,:,:) = 0.
    do iscal=1,nscal
      scalars(iscal)%so(:,:,:) = 0.
      scalars(iscal)%sadv(:,:,:) = 0.
    end do
  end if
  allocate(wall_heat_old(nscal),sforce_balance(4,nscal))
  wall_force_old(:) = -huge(1._rp) ! (dummy values, saved also where nothing is forced; see `fld_dt.out`)
  wall_heat_old(:)  = -huge(1._rp)
  sforce_balance(:,:) = 0.
  is_first = .true.
  do iscal=1,nscal
    scalars(iscal)%val(:,:,:) = 0.
  end do
  psi(:,:,:) = 0.
  do imark=1,nmk
    markers(imark)%psi(:,:,:) = 0.
    if(.not.is_ifc_split) markers(imark)%dpsidtrko(:,:,:) = 0.
    markers(imark)%kappa(:,:,:) = 0.
    if(.not.is_ifc_plic) then
      markers(imark)%normx(:,:,:) = 0.
      markers(imark)%normy(:,:,:) = 0.
      markers(imark)%normz(:,:,:) = 0.
    end if
    if(is_sdf_normals) markers(imark)%phi(:,:,:) = 0.
  end do
  psio(:,:,:) = 0.
  kappa(:,:,:) = 0.
  if(.not.is_ifc_plic) then
    normx(:,:,:) = 0.
    normy(:,:,:) = 0.
    normz(:,:,:) = 0.
  end if
  if(is_sdf_normals) phi(:,:,:) = 0.
  psiflx_x(:,:,:) = 0.
  psiflx_y(:,:,:) = 0.
  psiflx_z(:,:,:) = 0.
  if(.not.restart) then
    istep = 0
    time = 0.
    call initflow(inivel,cbcvel,bcvel,ng,lo,nh_v,nh_p,nh_g,l,dl,zc,zf,dzc,dzf,rho12(2),mu12(2), &
                  is_forced,velf,bforce,is_wallturb,time,u,v,w,p,is_velf)
    do iscal=1,nscal
      call initscal(iniscal(iscal),bcscal(:,:,iscal),ng,lo,nh_s,nh_g,l,dl,zc,dzf,is_sforced(iscal),scalf(iscal), &
                    scalars(iscal)%val)
    end do
    call ifc_init_field(inipsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,psi)
    do imark=1,nmk
      call ifc_init_field(inipsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,markers(imark)%psi,imark)
    end do
    if(myid == 0) print*, '*** Initial condition successfully set ***'
  else
    do is=1,nio
      call load_one('r',trim(datadir)//'fld_'//trim(c_io_vars(is))//trim(io_ext), &
                    MPI_COMM_WORLD,ng,nh_io(:,is),nh_g,lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                    is_check_time=is > 1)
    end do
    if(myid == 0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
    if(is_full_restart) then
      !
      ! terms of the previous time step, if they were all saved; otherwise the time integration restarts
      ! as in a new run, with a first-order step
      !
      is_first = .false.
      if(myid == 0) then
        do is=nio+1,nio_s
          inquire(file=trim(datadir)//'fld_'//trim(c_io_vars(is))//trim(io_ext),exist=is_done)
          is_first = is_first.or..not.is_done
        end do
        open(newunit=iunit,file=trim(datadir)//'fld_dt.out',status='old',action='read',iostat=istat)
        if(istat == 0) then
          read(iunit,*,iostat=istat) dt_restart(1:2),wall_force_old(1:3),wall_heat_old(1:nscal)
          if(istat /= 0 .and. nscal > 0) then
            !
            ! a file written before the heat losses of the scalars were saved: the values before them, with the
            ! heat losses left as dummies
            !
            rewind(iunit)
            read(iunit,*,iostat=istat) dt_restart(1:2),wall_force_old(1:3)
            wall_heat_old(:) = -huge(1._rp)
          end if
          close(iunit)
        end if
        is_first = is_first.or.istat /= 0
      end if
      call MPI_BCAST(is_first      ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dt_restart    ,2,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(wall_force_old,3,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(wall_heat_old ,nscal,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
      if(.not.is_first) then
        do is=nio+1,nio_s
          call load_one('r',trim(datadir)//'fld_'//trim(c_io_vars(is))//trim(io_ext), &
                        MPI_COMM_WORLD,ng,nh_io(:,is),nh_g,lo,hi,io_vars(is)%arr,trim(c_io_vars(is)),time,istep, &
                        is_check_time=.true.)
        end do
        if(myid == 0) print*, '*** Terms of the previous time step loaded. ***'
      else
        if(myid == 0) print*, 'WARNING: terms of the previous time step not found; restarting with a first-order step.'
      end if
    end if
    if(is_overwrite_save.and.nsaves_max > 0) then
      if(myid == 0) then
        istat = read_link(trim(datadir)//'fld_'//trim(c_io_vars(1))//trim(io_ext),fname_save)
        if(istat == 0) read(fname_save(5:8),'(i4)',iostat=istat) savecounter
        if(istat /= 0) savecounter = 0
      end if
      call MPI_BCAST(savecounter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    end if
  end if
  !$acc enter data copyin(u,v,w,p,pp,po)
  !$acc enter data copyin(dudtrko,dvdtrko,dwdtrko)
  if(.not.is_ifc_split) then
    !$acc enter data copyin(dpsidtrko)
  end if
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
  !$acc enter data copyin(surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o)
#endif
  call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,u,v,w,restart)
  call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,p)
  if(is_ifc_split) then
    !$acc enter data copyin(uo,vo,wo,uadv,vadv,wadv)
    if(.not.is_first) then
      call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,uo,vo,wo,restart)
      if(is_vof_sweep_transport) call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,uadv,vadv,wadv,.false.)
    end if
  end if
#if defined(_CONSTANT_COEFFS_POISSON)
  !
  ! pressure of the previous time step, for the extrapolation of the pressure splitting
  !
  if(is_first) then
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1-nh_p,n(3)+nh_p
      do j=1-nh_p,n(2)+nh_p
        do i=1-nh_p,n(1)+nh_p
          po(i,j,k) = p(i,j,k)
        end do
      end do
    end do
  else
    call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,po)
  end if
#endif
  do iscal=1,nscal
    s => scalars(iscal)
    !$acc enter data copyin(s%val,s%dsdtrko) async(1)
    call boundp(cbcscal(:,:,iscal),n,nh_s,nh_g,bcscal(:,:,iscal),nb,is_bound,dl,dzc,s%val)
    if(is_ifc_split) then
      !$acc enter data copyin(s%so,s%sadv) async(1)
      if(.not.is_first) then
        if(is_vof_sweep_transport) then
          call boundp(cbcscal(:,:,iscal),n,nh_s,nh_g,bcscal(:,:,iscal),nb,is_bound,dl,dzc,s%sadv)
        else
          call boundp(cbcscal(:,:,iscal),n,nh_s,nh_g,bcscal(:,:,iscal),nb,is_bound,dl,dzc,s%so)
        end if
      end if
    end if
  end do
  !$acc enter data copyin(psi) create(kappa) async(1)
  if(.not.is_ifc_plic) then
    !$acc enter data create(normx,normy,normz) async(1)
  end if
  if(is_sdf_normals) then
    !$acc enter data create(phi) async(1)
  end if
  !$acc enter data create(psio,psiflx_x,psiflx_y,psiflx_z) async(1)
  if(.not.restart .and. is_solve_ns) then
    if(any(is_velf(:))) then
      !
      ! bulk velocity of the initial condition along a forced direction, before the interface method uses the
      ! velocity (not with a prescribed velocity field, which stays the one of `initflow`, as at every time step)
      !
      call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
      call set_bulk_velocity(n,nh_v,nh_psi,nh_g,is_velf,velf,grid_vol_ratio_c,grid_vol_ratio_f,wzc,rho12,psi,u,v,w)
      call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,u,v,w,restart)
    end if
  end if
  call ifc_init(n,nh_v,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,u,v,w,psi,phi,normx,normy,normz,kappa,gam)
  if(nmk > 0) then
    !$acc enter data copyin(markers(:))
    do imark=1,nmk
      !$acc enter data copyin(markers(imark)%psi) create(markers(imark)%kappa) async(1)
      if(.not.is_ifc_split) then
        !$acc enter data copyin(markers(imark)%dpsidtrko) async(1)
      end if
      if(.not.is_ifc_plic) then
        !$acc enter data create(markers(imark)%normx,markers(imark)%normy,markers(imark)%normz) async(1)
      end if
      if(is_sdf_normals) then
        !$acc enter data create(markers(imark)%phi) async(1)
      end if
    end do
    !$acc enter data create(fsx,fsy,fsz) async(1)
    call ifc_markers_init(n,nh_v,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,restart,markers,psi)
  end if
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !$acc wait ! not needed but to prevent possible future issues
  !$acc update self(u,v,w,p,psi,kappa)
  do imark=1,nmk
    !$acc update self(markers(imark)%psi,markers(imark)%kappa)
  end do
  do iscal=1,nscal
    !$acc update self(scalars(iscal)%val)
  end do
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
  call chkdt(n,nh_v,nh_g,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dt_cfl,gam,seps, &
             ka12,cp12)
  dt = min(cfl*dt_cfl,dtmax)
  if(dt_f > 0.) then
    if(dt_f > dt_cfl) then
      if(myid == 0) print*, 'WARNING: fixed time step exceeds estimated stability limit.'
      if(myid == 0) print*, 'dt_f = ', dt_f, 'Estimated stable dt = ', dt_cfl
    end if
    dt = dt_f
  end if
  dto = dt
  if(restart.and..not.is_first) then
    !
    ! time steps of the run that is restarted
    !
    if(dt_f <= 0.) dt = dt_restart(1)
    dto = dt_restart(2)
  end if
  if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
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
    alpha = 0. ! no correction of the pressure update for the implicit Z diffusion, whose coefficients vary (see `updatep`)
    if(is_impdiff_cc) alpha = -.5*visc0*dtrk
    dtrki = dtrk**(-1)
    !
    ! phase field update
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1-nh_psi,n(3)+nh_psi
      do j=1-nh_psi,n(2)+nh_psi
        do i=1-nh_psi,n(1)+nh_psi
          psio(i,j,k) = psi(i,j,k)
        end do
      end do
    end do
    if(is_track_interface .and. nmk > 0) then
      call ifc_markers_advance(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,wzc,dt,dto,gam,seps,u,v,w, &
                               is_first,istep,markers,psi,psiflx_x,psiflx_y,psiflx_z,uo,vo,wo,uadv,vadv,wadv, &
                               scalars(1:nscal))
    else if(is_track_interface) then
      call ifc_advance(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dli,dzc,dzf,dzci,dzfi,wzc,dt,dto,gam,seps,u,v,w, &
                       normx,normy,normz,phi,kappa,psi,psiflx_x,psiflx_y,psiflx_z,is_first,dpsidtrko, &
                       istep,uo,vo,wo,uadv,vadv,wadv,scalars(1:nscal))
    else
      call ifc_advecting_velocity(tm_coeff,is_first,istep,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dl,dzc,dzf,dzci,dzfi,wzc,dt,dto, &
                                  u,v,w,psi,psiflx_x,psiflx_y,psiflx_z,uo,vo,wo,uadv,vadv,wadv,scalars(1:nscal))
    end if
    do iscal=1,nscal
      s => scalars(iscal)
      call tm_scal(tm_coeff,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,is_bound,is_sforced(iscal),l,dli,dzci,dzfi,grid_vol_ratio_f,dt, &
                   ssource(iscal),rho12,ka12(:,iscal),cp12(:,iscal),alpha0(iscal),psi,u,v,w, &
                   psio,psiflx_x,psiflx_y,psiflx_z,is_first,wall_heat_old(iscal),s%dsdtrko,s%val,sforce_balance(:,iscal), &
                   wall_heat_imp,s%so,s%sadv,uadv,vadv,wadv)
      if(is_impdiff) then
        if(is_sforced(iscal)) then
          !
          ! zero-net-heating forcing: the heat lost through the boundaries by the terms integrated implicitly is
          ! balanced here, first with its value at the old time level (returned by `tm_scal`), before the solve
          !
          call cmpt_scal_mean(n,nh_s,nh_psi,nh_g,grid_vol_ratio_f,rho12*cp12(:,iscal),psi,s%val,stot_o,rhocp_mean)
          wall_heat_imp = wall_heat_imp*dt
          call impdiff_updt_forcing_s(n,nh_s,nh_psi,rho12*cp12(:,iscal),psi,wall_heat_imp,s%val)
        end if
        if(is_impdiff_1d) then
          call impdiff_z_scal_solve(n,nh_s,nh_psi,nh_g,cbcscal(:,:,iscal),bcscal(:,:,iscal),dzc,dzci,dzfi,0.5_rp*dt, &
                                    rho12*cp12(:,iscal),ka12(:,iscal),psi,s%val)
        else if(is_impdiff_cc) then
          call impdiff_cc_scal_updt(n,nh_s,0.5_rp*dt,s%val)
          call solve_helmholtz(n,nh_s,ng,hi,s%arrplan,s%normfft,-0.5*alpha0(iscal)*dt, &
                               s%lambdaxy,s%a,s%b,s%c,s%rhsb%x,s%rhsb%y,s%rhsb%z,is_bound,cbcscal(:,:,iscal),['c','c','c'],s%val)
#if !defined(_CONSTANT_COEFFS_POISSON)
        else
          call impdiff_vc_scal_updt(n,nh_s,nh_psi,0.5_rp*dt,rho12*cp12(:,iscal),psi,s%val)
          call boundp(cbcscal(:,:,iscal),n,nh_s,nh_g,bcscal(:,:,iscal),nb,is_bound,dl,dzc,s%val)
          call impdiff_vc_scal_solve(ng,lo,hi,n,nh_s,nh_psi,nh_g,cbcscal(:,:,iscal),cbcpre,is_bound,dli,dzc,dzf,dzci,dzfi, &
                                     0.5_rp*dt,rho12*cp12(:,iscal),ka12(:,iscal),psi,s%val)
#endif
        end if
        if(is_sforced(iscal)) then
          !
          ! and then with the small remainder that keeps the total of rho*cp*s unchanged by the implicit integration
          !
          call cmpt_scal_mean(n,nh_s,nh_psi,nh_g,grid_vol_ratio_f,rho12*cp12(:,iscal),psi,s%val,stot,rhocp_mean)
          stot = stot_o-stot
          call impdiff_updt_forcing_s(n,nh_s,nh_psi,rho12*cp12(:,iscal),psi,stot,s%val)
          sforce_balance(1,iscal) = sforce_balance(1,iscal) + (wall_heat_imp+stot)/dt
        end if
      end if
      call boundp(cbcscal(:,:,iscal),n,nh_s,nh_g,bcscal(:,:,iscal),nb,is_bound,dl,dzc,s%val)
      if(is_sforced(iscal)) then
        !
        ! measured rate of change of the mean of rho*cp*s, from its value at the start of the step (returned by
        ! `tm_scal`), and the means
        !
        call cmpt_scal_mean(n,nh_s,nh_psi,nh_g,grid_vol_ratio_f,rho12*cp12(:,iscal),psi,s%val,stot,rhocp_mean)
        sforce_balance(2,iscal) = (stot-sforce_balance(2,iscal))/dt
        sforce_balance(3,iscal) = stot
        sforce_balance(4,iscal) = stot/rhocp_mean
      end if
    end do
    if(.not.is_solve_ns) then
      call initflow(inivel,cbcvel,bcvel,ng,lo,nh_v,nh_p,nh_g,l,dl,zc,zf,dzc,dzf,rho12(2),mu12(2), &
                    is_forced,velf,bforce,is_wallturb,time,u,v,w,p,is_velf)
      !$acc wait(1)
      !$acc update device(u,v,w,p) async(1)
      call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,u,v,w)
    else
      !
      ! the weight of the mean density is left out of the gravity term along the periodic directions, where the
      ! net weight has nothing to balance it, and along those between walls (zero pressure gradient at both
      ! ends), where it only changes the pressure by the hydrostatic one of the mean density: the pressure that
      ! is solved for then varies less in space and time, which the pressure splitting benefits from
      !
      rho_av(:) = 0.
      is_cmpt_rho_av(:) = (abs(gacc(:)) > 0.) .and. (cbcpre(0,:)//cbcpre(1,:) == 'PP' .or. cbcpre(0,:)//cbcpre(1,:) == 'NN')
      if(is_cmpt_rho_av(1)) &
        call bulk_mean_12_stag(n,nh_psi,nh_g,1,grid_vol_ratio_f,wzc,psi,rho12,rho_av(1))
      if(is_cmpt_rho_av(2)) &
        call bulk_mean_12_stag(n,nh_psi,nh_g,2,grid_vol_ratio_f,wzc,psi,rho12,rho_av(2))
      if(is_cmpt_rho_av(3)) &
        call bulk_mean_12_stag(n,nh_psi,nh_g,3,grid_vol_ratio_c,wzc,psi,rho12,rho_av(3))
      !
      ! (between walls any uniform density would do: the mean one is rounded to single precision there, such
      ! that the round-off of the global sum, which depends on the domain decomposition, does not reach the
      ! flow, and runs with different numbers of tasks stay identical)
      !
      where(is_cmpt_rho_av(:) .and. cbcpre(0,:)//cbcpre(1,:) == 'NN') rho_av(:) = real(real(rho_av(:),sp),rp)
      if(nmk > 0) call ifc_markers_force(n,nh_psi,nh_ifc,nh_g,dli,dzci,wzc,markers,psi,fsx,fsy,fsz)
      call tm(tm_coeff,n,nh_v,nh_p,nh_s,nh_psi,nh_ifc,nh_g,is_bound,is_forced,l,dli,dzci,dzfi,wzc, &
              grid_vol_ratio_c,grid_vol_ratio_f, &
              dt,dt_r,bforce,gacc,sigma,rho_av,rho12,mu12,beta12,rho0,psi,kappa,p,po,scalars(1)%val, &
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
              surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o, &
#endif
              psio,psiflx_x,psiflx_y,psiflx_z,is_first,wall_force_old,dudtrko,dvdtrko,dwdtrko,u,v,w,force_balance, &
              wall_force_imp,uadv,vadv,wadv,fsx,fsy,fsz,mfields)
      if(is_forced_hit) then
        call lscale_forcing(forced_hit_type,lo,hi,nh_v,nh_g,forced_hit_amp,dtrk,l,dl,zc,u,v,w)
      end if
      if(is_impdiff) then
        !
        ! implicit diffusion (Crank-Nicolson): adds the explicit half and solves for the implicit one
        !
        if(any(is_forced(:))) then
          !
          ! zero-net-acceleration forcing: the wall shear of the terms integrated implicitly is balanced
          ! here, first with its value at the old time level (returned by `tm`), before the solve
          !
          call cmpt_momentum(n,nh_v,nh_psi,nh_g,is_forced,grid_vol_ratio_c,grid_vol_ratio_f,wzc,rho12,psi,u,v,w,momentum_o)
          wall_force_imp(:) = wall_force_imp(:)*dtrk
          call impdiff_updt_forcing(n,nh_v,nh_psi,nh_g,wzc,rho12,psi,wall_force_imp,u,v,w)
        end if
        if(is_impdiff_1d) then
          call impdiff_z_mom_solve(n,nh_v,nh_psi,nh_g,cbcvel,bcvel,dzc,dzf,dzci,dzfi,wzc,0.5_rp*dtrk,rho12,mu12,psi,u,v,w)
        else if(is_impdiff_cc) then
          call impdiff_cc_mom_updt(n,nh_v,0.5_rp*dtrk,u,v,w)
          call solve_helmholtz(n,nh_v,ng,hi,arrplanu,normfftu,alpha, &
                               lambdaxyu,au,bu,cu,rhsbu%x,rhsbu%y,rhsbu%z,is_bound,cbcvel(:,:,1),['f','c','c'],u)
          call solve_helmholtz(n,nh_v,ng,hi,arrplanv,normfftv,alpha, &
                               lambdaxyv,av,bv,cv,rhsbv%x,rhsbv%y,rhsbv%z,is_bound,cbcvel(:,:,2),['c','f','c'],v)
          call solve_helmholtz(n,nh_v,ng,hi,arrplanw,normfftw,alpha, &
                               lambdaxyw,aw,bw,cw,rhsbw%x,rhsbw%y,rhsbw%z,is_bound,cbcvel(:,:,3),['c','c','f'],w)
#if !defined(_CONSTANT_COEFFS_POISSON)
        else
          call impdiff_vc_mom_updt(n,nh_v,nh_psi,nh_g,wzc,0.5_rp*dtrk,rho12,psi,u,v,w)
          call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,u,v,w)
          call impdiff_vc_mom_solve(ng,lo,hi,n,nh_v,nh_psi,nh_g,cbcvel,cbcpre,is_bound,dli,dzc,dzf,dzci,dzfi,wzc,0.5_rp*dtrk, &
                                    rho12,mu12,psi,u,v,w)
#endif
        end if
        if(any(is_forced(:))) then
          !
          ! and then with the small remainder that keeps the momentum unchanged by the implicit integration
          !
          call cmpt_momentum(n,nh_v,nh_psi,nh_g,is_forced,grid_vol_ratio_c,grid_vol_ratio_f,wzc,rho12,psi,u,v,w,momentum)
          momentum(:) = momentum_o(:)-momentum(:)
          call impdiff_updt_forcing(n,nh_v,nh_psi,nh_g,wzc,rho12,psi,momentum,u,v,w)
          force_balance(:,1) = force_balance(:,1) + (wall_force_imp(:)+momentum(:))*dtrki
          force_balance(:,5) = force_balance(:,5) + (wall_force_imp(:)+momentum(:))*dtrki
        end if
      end if
      call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,u,v,w,.false.,-1)
      call fillps(n,nh_v,nh_p,nh_g,dli,dzfi,dtrki*rho0,u,v,w,pp)
#if defined(_CONSTANT_COEFFS_POISSON)
      call updt_rhs_b(['c','c','c'],cbcpre,n,nh_p,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      call solver(n,nh_p,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],pp,is_dtdma_update_p,ap_d,cp_d, &
                  .not.is_impdiff)
#else
      call solver_vc(ng,lo,hi,nh_p,nh_psi,nh_g,cbcpre,bcpre,dli,dzci,dzf,wzc,is_bound,rho12,psi,pp,po)
#endif
      call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,pp)
      call correc(n,nh_v,nh_p,nh_psi,nh_g,dli,dzci,wzc,rho12,dtrk/rho0,pp,psi,u,v,w)
      call bounduvw(cbcvel_halo,n,nh_v,nh_g,bcvel_halo,nb,is_bound,dl,dzc,dzf,u,v,w,.true.,1,0)
      if(any(is_forced(:))) then
        call cmpt_momentum(n,nh_v,nh_psi,nh_g,is_forced,grid_vol_ratio_c,grid_vol_ratio_f,wzc,rho12,psi,u,v,w,momentum)
        where(is_forced(:)) force_balance(:,6) = (momentum(:)-force_balance(:,6))/dt
      end if
#if defined(_CONSTANT_COEFFS_POISSON)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1-nh_p,n(3)+nh_p
        do j=1-nh_p,n(2)+nh_p
          do i=1-nh_p,n(1)+nh_p
            po(i,j,k) = p(i,j,k)
          end do
        end do
      end do
#endif
      call updatep(n,nh_p,nh_g,dli,dzci,dzfi,alpha,pp,p)
      call boundp(cbcpre,n,nh_p,nh_g,bcpre,nb,is_bound,dl,dzc,p)
    end if
    is_first = .false.
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
    call ifc_update_params(n,nh_v,u,v,w,gam)
    if(icheck > 0.and.mod(istep,max(icheck,1)) == 0) then
      if(myid == 0) print*, 'Checking stability and divergence...'
      call chkdt(n,nh_v,nh_g,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dt_cfl,gam,seps, &
                 ka12,cp12)
      dt = min(cfl*dt_cfl,dtmax)
      if(dt_f > 0.) dt = dt_f
      if(myid == 0) print*, 'dt_cfl = ', dt_cfl, 'dt = ', dt
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
      var(1) = 1._rp*istep
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
        var(20:22) = force_balance(:,7)
        call out0d(trim(datadir)//'forcing_balance.out',22,var)
      end if
      do iscal=1,nscal
        if(is_sforced(iscal)) then
          write(scalnum,'(i3.3)') iscal
          var(1) = time
          var(2:5) = sforce_balance(1:4,iscal)
          call out0d(trim(datadir)//'forcing_scalar_'//scalnum//'.out',5,var)
        end if
      end do
      !
      call ifc_write_log(istep,time,n,nh_psi,nh_g,dl,dzf,gam,seps,psi)
    end if
    write(fldnum,'(i7.7)') istep
    if(is_output_due(istep,iout1d,tout1d,time_start,time)) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
      do imark=1,nmk
        !$acc update self(markers(imark)%psi)
      end do
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out1d.h90'
    end if
    if(is_output_due(istep,iout2d,tout2d,time_start,time)) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
      do imark=1,nmk
        !$acc update self(markers(imark)%psi)
      end do
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      include 'out2d.h90'
    end if
    if(is_output_due(istep,iout3d,tout3d,time_start,time)) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
      do imark=1,nmk
        !$acc update self(markers(imark)%psi,markers(imark)%kappa)
      end do
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
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
      do imark=1,nmk
        !$acc update self(markers(imark)%psi)
      end do
      do iscal=1,nscal
        !$acc update self(scalars(iscal)%val)
      end do
      if(is_full_restart) then
        !$acc update self(po,dudtrko,dvdtrko,dwdtrko)
        if(.not.is_ifc_split) then
          !$acc update self(dpsidtrko)
          do imark=1,nmk
            !$acc update self(markers(imark)%dpsidtrko)
          end do
        end if
        do iscal=1,nscal
          !$acc update self(scalars(iscal)%dsdtrko)
        end do
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
        !$acc update self(surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o)
#endif
        if(is_ifc_split) then
          !$acc update self(uo,vo,wo)
          if(is_vof_sweep_transport) then
            !$acc update self(uadv,vadv,wadv)
          end if
          do iscal=1,nscal
            !$acc update self(scalars(iscal)%so,scalars(iscal)%sadv)
          end do
        end if
      end if
      is_save_tmp = trim(io_ext) /= '.bp'
      do is=1,nio_s
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
          do is=1,nio_s
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
        do is=1,nio_s
          call gen_alias(myid,trim(datadir),trim(filename)//'_'//trim(c_io_vars(is))//trim(io_ext), &
                         'fld_'//trim(c_io_vars(is))//trim(io_ext))
        end do
      end if
      if(is_full_restart) then
        !
        ! time steps, and wall force and heat loss (of each scalar) of the previous time step
        !
        if(myid == 0) then
          open(newunit=iunit,file=trim(datadir)//trim(filename)//'_dt.out',status='replace')
          write(iunit,'(*(es24.16e3,1x))') dt,dto,wall_force_old(1:3),wall_heat_old(1:nscal)
          close(iunit)
        end if
        if(.not.is_overwrite_save.or.nsaves_max > 0) then
          call gen_alias(myid,trim(datadir),trim(filename)//'_dt.out','fld_dt.out')
        end if
      end if
      if(is_overwrite_save.and.nsaves_max > 0) then
        var(1) = 1._rp*istep
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
  if(is_impdiff_cc) then
    call fftend(arrplanu)
    call fftend(arrplanv)
    call fftend(arrplanw)
    do iscal=1,nscal
      call fftend(scalars(iscal)%arrplan)
    end do
  end if
#if !defined(_CONSTANT_COEFFS_POISSON)
  call solver_vc_end
#endif
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
#if defined(_OPENACC)
  call cudecomp_finalize
#endif
  call MPI_FINALIZE(ierr)
  if(kill) error stop
end program cans
