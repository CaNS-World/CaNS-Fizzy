! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
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
#if defined(_DEBUG)
  use, intrinsic :: iso_fortran_env, only: compiler_version,compiler_options
#endif
  use, intrinsic :: iso_c_binding  , only: C_PTR
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan
  use mpi
  use decomp_2d
  use mod_bound          , only: boundp,bounduvw,boundnor,updt_rhs_b
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,ierr
  use mod_correc         , only: correc
  use mod_fft            , only: fftini,fftend
  use mod_fillps         , only: fillps
  use mod_forcing        , only: lscale_forcing
  use mod_initflow       , only: initflow,initscal
  use mod_initgrid       , only: initgrid
  use mod_initmpi        , only: initmpi
  use mod_initsolver     , only: initsolver
  use mod_load           , only: load_one
  use mod_rk             , only: tm => rk,tm_scal => rk_scal,tm_2fl => rk_2fl
  use mod_output         , only: out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d
  use mod_param          , only: rkcoeff,small, &
                                 nb,is_bound,cbcvel,bcvel,cbcpre,bcpre,cbcsca,bcsca,cbcpsi,bcpsi,cbcnor,bcnor, &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 nstep,time_max,tw_max,stop_type,restart,is_overwrite_save,nsaves_max, &
                                 datadir,   &
                                 is_solve_ns,is_track_interface, &
                                 cfl,dtmin,dt_f, &
                                 inivel,inisca,inipsi, &
                                 is_wallturb,is_forced_hit, &
                                 dims, &
                                 gtype,gr, &
                                 bforce,ssource, &
                                 ng,l,dl,dli, &
                                 read_input, &
                                 rho0,rho12,mu12,sigma,gacc,ka12,cp12,beta12, &
                                 acdi_gam_factor,acdi_gam_min,acdi_eps_factor
#if 1
  use mod_sanity         , only: test_sanity_input
#endif
  use mod_acdi           , only: acdi_set_epsilon,acdi_set_gamma,acdi_cmpt_phi
  use mod_two_fluid      , only: init2fl,cmpt_norm_curv => cmpt_norm_curv_youngs
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mod_solver_vc      , only: solver_vc
#endif
#if !defined(_OPENACC)
  use mod_solver         , only: solver
#else
  use mod_solver_gpu     , only: solver => solver_gpu
  use mod_workspaces     , only: init_wspace_arrays,set_cufft_wspace
  use mod_common_cudecomp, only: istream_acc_queue_1
#endif
  use mod_timer          , only: timer_tic,timer_toc,timer_print
  use mod_updatep        , only: updatep
  use mod_utils          , only: bulk_mean,bulk_mean_12
  !@acc use mod_utils    , only: device_memory_footprint
  use mod_types
  implicit none
  integer , dimension(3) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,pp,pn,po
  real(rp) :: rho_av
#if !defined(_OPENACC)
  type(C_PTR), dimension(2,2) :: arrplanp
#else
  integer    , dimension(2,2) :: arrplanp
#endif
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  real(rp) :: normfftp
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  type(rhs_bound) :: rhsbp
  real(rp) :: alpha
  real(rp) :: dt,dto,dt_r,dti,dtmax,dtrk,dtrki,time,divtot,divmax
  real(rp) :: gam,seps
  integer :: irk,istep
  real(rp), allocatable, dimension(:) :: dzc  ,dzf  ,zc  ,zf  ,dzci  ,dzfi, &
                                         dzc_g,dzf_g,zc_g,zf_g,dzci_g,dzfi_g, &
                                         grid_vol_ratio_c,grid_vol_ratio_f
  real(rp), dimension(3) :: dpdl
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(42) :: var
#if defined(_TIMING)
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  real(rp) :: twi,tw
  !
  integer  :: savecounter
  character(len=7  ) :: fldnum
  character(len=4  ) :: chkptnum
  character(len=100) :: filename,fexts(6)
  integer :: k,kk
  logical :: is_done,kill
  real(rp), dimension(2) :: tm_coeff
  !
  ! scalar field
  !
  real(rp), allocatable, dimension(:,:,:) :: s
  !
  ! two-fluid solver specific
  !
  real(rp), allocatable, dimension(:,:,:) :: psi,phi,kappa,normx,normy,normz, &
                                             acdi_rgx,acdi_rgy,acdi_rgz
  real(rp), allocatable, dimension(:,:,:) :: psio
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
  call initmpi(ng,dims,cbcpre,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,nb,is_bound)
  twi = MPI_WTIME()
  savecounter = 0
  !
  ! allocate variables
  !
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
#if !defined(_CONSTANT_COEFFS_POISSON)
  allocate(po,mold=pp)
  po(:,:,:) = 0._rp
#else
  pp(:,:,:) = 0._rp
  allocate(pn,mold=pp)
  allocate(po,mold=pp)
  pn(:,:,:) = pp(:,:,:)
  po(:,:,:) = pn(:,:,:)
#endif
#if defined(_SCALAR)
  allocate(s,mold=pp)
#endif
  allocate(lambdaxyp(n_z(1),n_z(2)))
  allocate(ap(n_z(3)),bp(n_z(3)),cp(n_z(3)))
  allocate(dzc( 0:n(3)+1), &
           dzf( 0:n(3)+1), &
           zc(  0:n(3)+1), &
           zf(  0:n(3)+1), &
           dzci(0:n(3)+1), &
           dzfi(0:n(3)+1))
  allocate(dzc_g( 0:ng(3)+1), &
           dzf_g( 0:ng(3)+1), &
           zc_g(  0:ng(3)+1), &
           zf_g(  0:ng(3)+1), &
           dzci_g(0:ng(3)+1), &
           dzfi_g(0:ng(3)+1))
  allocate(grid_vol_ratio_c,mold=dzc)
  allocate(grid_vol_ratio_f,mold=dzf)
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
  allocate(psi,phi,kappa,normx,normy,normz,mold=pp)
#if defined(_CONSERVATIVE_MOMENTUM)
  allocate(acdi_rgx,acdi_rgy,acdi_rgz,mold=pp)
  if(.not. allocated(psio)) allocate(psio(0:n(1)+1,0:n(2)+1,0:n(3)+1))
#endif
#if defined(_DEBUG)
  if(myid == 0) print*, 'This executable of CaNS was built with compiler: ', compiler_version()
  if(myid == 0) print*, 'Using the options: ', compiler_options()
  block
    character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version
    integer :: ilen
    call MPI_GET_LIBRARY_VERSION(mpi_version,ilen,ierr)
    if(myid == 0) print*, 'MPI Version: ', trim(mpi_version)
  end block
  if(myid == 0) print*, ''
#endif
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, '*** Beginning of simulation ***'
  if(myid == 0) print*, '*******************************'
  if(myid == 0) print*, ''
  call initgrid(gtype,ng(3),gr,l(3),dzc_g,dzf_g,zc_g,zf_g)
  if(myid == 0) then
    open(99,file=trim(datadir)//'grid.bin',action='write',form='unformatted',access='stream',status='replace')
    write(99) dzc_g(1:ng(3)),dzf_g(1:ng(3)),zc_g(1:ng(3)),zf_g(1:ng(3))
    close(99)
    open(99,file=trim(datadir)//'grid.out')
    do kk=0,ng(3)+1
      write(99,*) 0.,zf_g(kk),zc_g(kk),dzf_g(kk),dzc_g(kk)
    end do
    close(99)
    open(99,file=trim(datadir)//'geometry.out')
      write(99,*) ng(1),ng(2),ng(3)
      write(99,*) l(1),l(2),l(3)
    close(99)
  end if
  !$acc enter data copyin(lo,hi,n) async
  !$acc enter data copyin(bforce,gacc,dl,dli,l) async
  !$acc enter data copyin(zc_g,zf_g,dzc_g,dzf_g) async
  !$acc enter data create(zc,zf,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g) async
  !
  !$acc enter data copyin(rho12,mu12,ka12,cp12,beta12) async
  !
  !$acc parallel loop default(present) private(k) async
  do kk=lo(3)-1,hi(3)+1
    k = kk-(lo(3)-1)
    zc( k) = zc_g(kk)
    zf( k) = zf_g(kk)
    dzc(k) = dzc_g(kk)
    dzf(k) = dzf_g(kk)
    dzci(k) = dzc(k)**(-1)
    dzfi(k) = dzf(k)**(-1)
  end do
  !$acc kernels default(present) async
  dzci_g(:) = dzc_g(:)**(-1)
  dzfi_g(:) = dzf_g(:)**(-1)
  !$acc end kernels
  !$acc enter data create(grid_vol_ratio_c,grid_vol_ratio_f) async
  !$acc kernels default(present) async
  grid_vol_ratio_c(:) = dl(1)*dl(2)*dzc(:)/(l(1)*l(2)*l(3))
  grid_vol_ratio_f(:) = dl(1)*dl(2)*dzf(:)/(l(1)*l(2)*l(3))
  !$acc end kernels
  !$acc update self(zc,zf,dzc,dzf,dzci,dzfi) async
  !$acc exit data copyout(zc_g,zf_g,dzc_g,dzf_g,dzci_g,dzfi_g) async ! not needed on the device
  !$acc wait
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity_input(ng,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre)
  !
  ! initialize Poisson solver
  !
  call initsolver(ng,n_x_fft,n_y_fft,lo_z,hi_z,dli,dzci_g,dzfi_g,cbcpre,bcpre(:,:), &
                  lambdaxyp,['c','c','c'],ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
  !$acc enter data copyin(lambdaxyp,ap,bp,cp) async
  !$acc enter data copyin(rhsbp,rhsbp%x,rhsbp%y,rhsbp%z) async
  !$acc wait
#if defined(_OPENACC)
  !
  ! determine workspace sizes and allocate the memory
  !
  call init_wspace_arrays()
  call set_cufft_wspace(pack(arrplanp,.true.),istream_acc_queue_1)
  if(myid == 0) print*,'*** Device memory footprint (Gb): ', &
                  device_memory_footprint(n,n_z)/(1._sp*1024**3), ' ***'
#endif
#if defined(_DEBUG_SOLVER)
  call test_sanity_solver(ng,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,dli,dzc,dzf,dzci,dzfi,dzci_g,dzfi_g, &
                          nb,is_bound,cbcvel,cbcpre,bcvel,bcpre)
#endif
  !
  call acdi_set_epsilon(dl,dzfi,acdi_eps_factor,seps)
  !
  fexts(1) = 'u'
  fexts(2) = 'v'
  fexts(3) = 'w'
  fexts(4) = 'p'
  fexts(5) = 'psi'
#if defined(_SCALAR)
  fexts(6) = 's'
#endif
  if(.not.restart) then
    istep = 0
    time = 0.
    !$acc update self(zc,dzc,dzf)
    call initflow(inivel,bcvel,ng,lo,l,dl,zc,zf,dzc,dzf,rho12(2),mu12(2),bforce,is_wallturb,time,u,v,w,p)
#if defined(_SCALAR)
    call initscal(inisca,bcsca,ng,lo,l,dl,dzf,zc,s)
#endif
    call init2fl(inipsi,cbcpsi,seps,lo,hi,l,dl,zc_g,psi)
    if(myid == 0) print*, '*** Initial condition succesfully set ***'
  else
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(1))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,u,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(2))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,v,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(3))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,w,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(4))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,p,time,istep)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(5))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,psi,time,istep)
#if defined(_SCALAR)
    call load_one('r',trim(datadir)//'fld_'//trim(fexts(6))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,s,time,istep)
#endif
    if(myid == 0) print*, '*** Checkpoints loaded at time = ', time, 'time step = ', istep, '. ***'
  end if
  !$acc enter data copyin(u,v,w,p,pp,pn,po) async
  !$acc wait
  call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
  !$acc kernels default(present) async(1)
  pn(:,:,:) =  p(:,:,:)
  po(:,:,:) = pn(:,:,:)
  !$acc end kernels
#if defined(_SCALAR)
  !$acc enter data copyin(s)
  call boundp(cbcsca,n,bcsca,nb,is_bound,dl,dzc,s)
#endif
  !$acc enter data copyin(psi) create(phi,kappa,normx,normy,normz)
  !$acc enter data create(acdi_rgx,acdi_rgy,acdi_rgz)
  !$acc enter data create(psio)
  call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,psi)
  !
  call acdi_cmpt_phi(n,seps,psi,phi)
  call cmpt_norm_curv(n,dli,dzci,dzfi,phi,normx,normy,normz,kappa)
  call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,kappa)
  call boundnor(cbcnor,n,bcnor,nb,is_bound,dl,dzc,normx,normy,normz)
  !
#if defined(_CONSERVATIVE_MOMENTUM)
  !$acc kernels default(present) async(1)
  psio(:,:,:) = psi(:,:,:)
  !$acc end kernels
#endif
  !
  call acdi_set_gamma(n,acdi_gam_factor,u,v,w,gam)
  gam = max(gam,acdi_gam_min)
  if(myid == 0) print*, 'ACDI parameters. Gamma: ', gam, 'Epsilon: ', seps
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !$acc wait
  !$acc update self(u,v,w,p,psi,kappa)
#include "out1d.h90"
#include "out2d.h90"
#include "out3d.h90"
  !
  call chkdt(n,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dtmax,gam,seps,ka12,cp12)
  dt = min(cfl*dtmax,dtmin); if(dt_f > 0.) dt = dt_f
  if(myid == 0) print*, 'dtmax = ', dtmax, 'dt = ', dt
  dto = dt
  dti = 1./dt
  kill = .false.
  !
  ! main loop
  !
  if(myid == 0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)
#if defined(_TIMING)
    !$acc wait(1)
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time = time + dt
    if(myid == 0) print*, 'Time step #', istep, 'Time = ', time
    do irk=1,3
      tm_coeff(:) = rkcoeff(:,irk)
      dtrk = sum(rkcoeff(:,1:irk))*dt
      dtrki = dtrk**(-1)
      dt_r = dtrk/dto
      !
      ! phase field update
      !
#if defined(_CONSERVATIVE_MOMENTUM)
      !$acc kernels default(present) async(1)
      psio(:,:,:)   = psi(:,:,:)
      !$acc end kernels
#endif
      if(is_track_interface) then
        call tm_2fl(tm_coeff,n,dli,dzci,dzfi,dt,gam,seps,u,v,w,normx,normy,normz,phi,psi,acdi_rgx,acdi_rgy,acdi_rgz)
#if defined(_CONSERVATIVE_MOMENTUM)
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,acdi_rgx,acdi_rgy,acdi_rgz)
#endif
        call boundp(cbcpsi,n,bcpsi,nb,is_bound,dl,dzc,psi)
        call acdi_cmpt_phi(n,seps,psi,phi)
        call cmpt_norm_curv(n,dli,dzci,dzfi,phi,normx,normy,normz,kappa)
        call boundp(cbcpsi,n,bcpre,nb,is_bound,dl,dzc,kappa)
        call boundnor(cbcnor,n,bcnor,nb,is_bound,dl,dzc,normx,normy,normz)
      end if
#if defined(_SCALAR)
      call tm_scal(tm_coeff,n,dli,dzci,dzfi,dt,ssource,rho12,ka12,cp12,psi,u,v,w,s)
      call boundp(cbcsca,n,bcsca,nb,is_bound,dl,dzc,s)
#endif
      if(.not.is_solve_ns) then
        call initflow(inivel,bcvel,ng,lo,l,dl,zc,zf,dzc,dzf,rho12(2),mu12(2),bforce,is_wallturb,time,u,v,w,p)
        !$acc wait(1)
        !$acc update device(u,v,w,p) async(1)
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
      else
        rho_av = 0.
        if(any(abs(gacc(:))>0. .and. cbcpre(0,:)//cbcpre(1,:) == 'PP')) then
          call bulk_mean_12(n,grid_vol_ratio_c,psi,rho12,rho_av)
        end if
        call tm(tm_coeff,n,dli,dzci,dzfi,dt,dt_r, &
                bforce,gacc,sigma,rho_av,rho12,mu12,beta12,rho0,psi,kappa,s,p,pn,po,psio, &
                acdi_rgx,acdi_rgy,acdi_rgz,u,v,w)
        if(is_forced_hit) then
          call lscale_forcing(2,lo,hi,0.5_rp,dtrk,l,dl,zc,zf,u,v,w)
        end if
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.false.,dl,dzc,dzf,u,v,w)
        !$acc kernels default(present) async(1)
        pp(:,:,:) = p(:,:,:)
        !$acc end kernels
        call fillps(n,dli,dzfi,dtrki,rho0,u,v,w,p)
#if defined(_CONSTANT_COEFFS_POISSON)
        call updt_rhs_b(['c','c','c'],cbcpre,n,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,p)
        call solver(n,ng,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre,['c','c','c'],p)
#else
        call solver_vc(ng,lo,hi,cbcpre,bcpre,dli,dzci,dzfi,is_bound,rho12,psi,p,po)
#endif
        call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
        call correc(n,dli,dzci,rho0,rho12,dtrk,p,psi,u,v,w)
        call bounduvw(cbcvel,n,bcvel,nb,is_bound,.true.,dl,dzc,dzf,u,v,w)
        call updatep(pp,p)
        call boundp(cbcpre,n,bcpre,nb,is_bound,dl,dzc,p)
      end if
    end do
    !$acc kernels default(present) async(1)
    po(:,:,:) = pn(:,:,:)
    pn(:,:,:) =  p(:,:,:)
    !$acc end kernels
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
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    end if
    dto = dt
    if(mod(istep,1) == 0) then
      call acdi_set_gamma(n,acdi_gam_factor,u,v,w,gam)
      gam = max(gam,acdi_gam_min)
    end if
    if(mod(istep,icheck) == 0) then
      if(myid == 0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,is_solve_ns,is_track_interface,mu12,rho12,sigma,gacc,u,v,w,dtmax,gam,seps)
      dt = min(cfl*dtmax,dtmin); if(dt_f > 0.) dt = dt_f
      if(myid == 0) print*, 'dtmax = ', dtmax, 'dt = ', dt
      if(dtmax < small) then
        if(myid == 0) print*, 'ERROR: time step is too small.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
      dti = 1./dt
      call chkdiv(lo,hi,dli,dzfi,u,v,w,divtot,divmax)
      if(myid == 0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
#if !defined(_MASK_DIVERGENCE_CHECK)
      if(divmax > small.or.is_nan(divtot)) then
        if(myid == 0) print*, 'ERROR: maximum divergence is too large.'
        if(myid == 0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
#endif
    end if
    !
    ! output routines below
    !
    if(mod(istep,iout0d) == 0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
      !
      var(1) = 1.*istep
      var(2) = time
      var(3) = gam
      var(4) = seps
      call out0d(trim(datadir)//'log_acdi.out',4,var)
    end if
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
#include "out1d.h90"
    end if
    if(mod(istep,iout2d) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
#include "out2d.h90"
    end if
    if(mod(istep,iout3d) == 0) then
      !$acc wait
      !$acc update self(u,v,w,p,psi,kappa)
#include "out3d.h90"
    end if
    if(mod(istep,isave ) == 0.or.(is_done.and..not.kill)) then
      if(is_overwrite_save) then
        filename = 'fld'
      else
        filename = 'fld_'//fldnum
        if(nsaves_max > 0) then
          if(savecounter >= nsaves_max) savecounter = 0
          savecounter = savecounter + 1
          write(chkptnum,'(i4.4)') savecounter
          filename = 'fld_'//chkptnum
          var(1) = 1.*istep
          var(2) = time
          var(3) = 1.*savecounter
          call out0d(trim(datadir)//'log_checkpoints.out',3,var)
        end if
      end if
      !$acc wait
      !$acc update self(u,v,w,p,psi)
      call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(1))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,u,time,istep)
      call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(2))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,v,time,istep)
      call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(3))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,w,time,istep)
      call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(4))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,p,time,istep)
      call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(5))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,psi,time,istep)
#if defined(_SCALAR)
      call load_one('w',trim(datadir)//trim(filename)//'_'//trim(fexts(6))//'.bin',MPI_COMM_WORLD,ng,[1,1,1],lo,hi,s,time,istep)
#endif
      if(.not.is_overwrite_save) then
        !
        ! fld_?.bin -> last checkpoint file (symbolic link)
        !
        do k = 1,5
          call gen_alias(myid,trim(datadir),trim(filename)//'_'//trim(fexts(k))//'.bin','fld_'//trim(fexts(k))//'.bin')
        end do
#if defined(_SCALAR)
        call gen_alias(myid,trim(datadir),trim(filename)//'_'//trim(fexts(k))//'.bin','fld_'//trim(fexts(k))//'.bin') ! k = 6 now
#endif
      end if
      if(myid == 0) print*, '*** Checkpoints saved at time = ', time, 'time step = ', istep, '. ***'
    end if
#if defined(_TIMING)
      !$acc wait(1)
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) print*, 'Avrg, min & max elapsed time: '
      if(myid == 0) print*, dt12av/(1.*product(dims)),dt12min,dt12max
#endif
  end do
  !
  ! clear ffts
  !
  call fftend(arrplanp)
  if(myid == 0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
end program cans
