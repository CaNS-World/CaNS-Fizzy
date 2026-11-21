! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_interface_capturing
  !
  ! common interface to the interface-capturing methods, selected at compile time:
  !   - ACDI diffuse-interface method (default)
  !   - THINC/QQ volume-of-fluid method (`INTERFACE_CAPTURING_VOF=1`)
  !   - geometric (PLIC) volume-of-fluid method with conservative direction-split advection
  !     (`INTERFACE_CAPTURING_SPLIT_VOF=1`)
  !
  ! the main program calls, in this order:
  !   ifc_init_params   -> method parameters needed by the initial condition
  !   ifc_init_field    -> initial phase indicator
  !   ifc_init          -> initialization from the initial (or loaded) phase indicator
  !   ifc_advance       -> one time step of the phase indicator, returning the (volumetric) phase fluxes
  !                        used for the mass-consistent transport of momentum and scalars
  !   ifc_update_params -> update of the method parameters after each time step
  !   ifc_write_log     -> log file of the method
  !
  use mod_types
  use mod_common_mpi, only: myid
  use mod_param     , only: datadir,ng,nb,is_bound,cbcvel,bcvel,cbcpsi,bcpsi,cbcnor,bcnor,cbccur,bccur, &
                            acdi_gam_factor,acdi_gam_min,vof_thinc_beta
  use mod_bound     , only: boundp,bounduvw
  use mod_output    , only: out0d
  use mod_two_fluid , only: init2fl,init2fl_sharp,cmpt_norm_curv => cmpt_norm_curv_youngs
  use mod_rk        , only: rk_2fl
  use mod_acdi      , only: acdi_set_epsilon
#if defined(_INTERFACE_CAPTURING_SPLIT_VOF)
  use mod_vof_split , only: vof_split_init,vof_split_advance,vof_split_diagnostics
#elif !defined(_INTERFACE_CAPTURING_VOF)
  use mod_acdi      , only: acdi_set_gamma,acdi_cmpt_phi
#elif defined(_SDF_NORMALS)
  use mod_vof_thinc_qq, only: vof_thinc_cmpt_phi
#endif
  implicit none
  private
  public ifc_init_params,ifc_init_field,ifc_init,ifc_advance,ifc_update_params,ifc_write_log
#if defined(_INTERFACE_CAPTURING_SPLIT_VOF)
  real(rp) :: vol0 ! initial volume of phase 1
#endif
  contains
  subroutine ifc_init_params(dl,nh_g,dzfi,psi_thickness_factor,seps)
    !
    ! method parameters needed by the initial condition: the interface width, used by all methods
    !
    implicit none
    real(rp), intent(in   ), dimension(3)       :: dl
    integer , intent(in   )                     :: nh_g
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzfi
    real(rp), intent(in   ) :: psi_thickness_factor
    real(rp), intent(inout) :: seps
    call acdi_set_epsilon(dl,nh_g,dzfi,psi_thickness_factor,seps)
  end subroutine ifc_init_params
  !
  subroutine ifc_init_field(inipsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,psi)
    !
    ! initial condition for the phase indicator: a smooth profile of width `seps`, or
    ! the exact volume fractions for the geometric volume-of-fluid method
    !
    implicit none
    character(len=*), intent(in) :: inipsi
    real(rp), intent(in ) :: seps
    integer , intent(in ), dimension(3)        :: lo,hi
    integer , intent(in )                      :: nh_psi,nh_g
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(1-nh_g:)  :: zc_g,dzf_g
    real(rp), intent(out), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:) :: psi
#if !defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    call init2fl(inipsi,cbcpsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,psi)
#else
    call init2fl_sharp(inipsi,cbcpsi,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,psi)
#endif
  end subroutine ifc_init_field
  !
  subroutine ifc_init(n,nh_v,nh_psi,nh_g,dl,dli,dzc,dzf,dzci,dzfi,seps,u,v,w,psi,phi,normx,normy,normz,kappa,gam)
    !
    ! initialization from the initial (or loaded) phase indicator, already present on the device
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi
    real(rp), intent(in   ) :: seps
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: u,v,w
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,phi,normx,normy,normz,kappa
    real(rp), intent(inout) :: gam
#if defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    real(rp) :: clip,fmin,fmax,cfl(3)
    integer  :: nsub,i,j,k
#endif
    !
    call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
    !$acc wait
#if defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    call vof_split_init(ng,n,nh_psi)
    call vof_split_diagnostics(n,nh_psi,nh_g,dl,dzf,psi,vol0,clip,fmin,fmax,cfl,nsub)
    if(myid == 0) print*, 'Geometric VoF. Initial volume of phase 1: ', vol0
    !
    ! the interface normals and curvature are not computed yet with this method, which for now only
    ! transports the interface (see `ifc_geometry`)
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1-nh_psi,n(3)+nh_psi
      do j=1-nh_psi,n(2)+nh_psi
        do i=1-nh_psi,n(1)+nh_psi
          normx(i,j,k) = 0.
          normy(i,j,k) = 0.
          normz(i,j,k) = 0.
          kappa(i,j,k) = 0.
        end do
      end do
    end do
#endif
    call ifc_geometry(n,nh_psi,nh_g,dl,dli,dzc,dzci,dzfi,seps,psi,phi,normx,normy,normz,kappa)
#if !defined(_INTERFACE_CAPTURING_VOF) && !defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    call acdi_set_gamma(n,nh_v,acdi_gam_factor,u,v,w,gam)
    gam = max(gam,acdi_gam_min)
    if(myid == 0) print*, 'ACDI parameters. Gamma: ', gam, 'Epsilon: ', seps
#endif
  end subroutine ifc_init
  !
  subroutine ifc_advance(tm_coeff,n,nh_v,nh_psi,nh_g,dl,dli,dzc,dzf,dzci,dzfi,dt,gam,seps,u,v,w, &
                         normx,normy,normz,phi,kappa,psi,psiflx_x,psiflx_y,psiflx_z,is_first,dpsidtrko)
    !
    ! advances the phase indicator over one time step, returns the phase fluxes (with their boundary
    ! conditions), and updates the interface normals and curvature
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: tm_coeff
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi
    real(rp), intent(in   ) :: dt,gam,seps
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:)       :: u,v,w
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: normx,normy,normz,phi,kappa,psi
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psiflx_x,psiflx_y,psiflx_z
    logical , intent(in   ) :: is_first
    real(rp), intent(inout), dimension(:,:,:) :: dpsidtrko
    !
#if !defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    call rk_2fl(tm_coeff,n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,dt,gam,seps,vof_thinc_beta,u,v,w,normx,normy,normz,phi,psi, &
                psiflx_x,psiflx_y,psiflx_z,is_first,dpsidtrko)
    call bounduvw(cbcvel,n,nh_psi,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
    call boundp(cbcpsi,n,nh_psi,nh_g,bcpsi,nb,is_bound,dl,dzc,psi)
#else
    call vof_split_advance(tm_coeff,n,nh_v,nh_psi,nh_g,dl,dzc,dzf,dzfi,dt,cbcpsi,bcpsi,nb,is_bound,u,v,w, &
                           psi,psiflx_x,psiflx_y,psiflx_z)
    call bounduvw(cbcvel,n,nh_psi,nh_g,bcvel,nb,is_bound,dl,dzc,dzf,psiflx_x,psiflx_y,psiflx_z)
#endif
    call ifc_geometry(n,nh_psi,nh_g,dl,dli,dzc,dzci,dzfi,seps,psi,phi,normx,normy,normz,kappa)
  end subroutine ifc_advance
  !
  subroutine ifc_geometry(n,nh_psi,nh_g,dl,dli,dzc,dzci,dzfi,seps,psi,phi,normx,normy,normz,kappa)
    !
    ! interface normals and curvature, with their boundary conditions
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzci,dzfi
    real(rp), intent(in   ) :: seps
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: phi,normx,normy,normz,kappa
    !
#if !defined(_INTERFACE_CAPTURING_SPLIT_VOF)
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
#endif
  end subroutine ifc_geometry
  !
  subroutine ifc_update_params(n,nh_v,u,v,w,gam)
    !
    ! update of the method parameters after a time step
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(inout) :: gam
#if !defined(_INTERFACE_CAPTURING_VOF) && !defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    call acdi_set_gamma(n,nh_v,acdi_gam_factor,u,v,w,gam)
    gam = max(gam,acdi_gam_min)
#endif
  end subroutine ifc_update_params
  !
  subroutine ifc_write_log(istep,time,n,nh_psi,nh_g,dl,dzf,gam,seps,psi)
    !
    ! log file of the method
    !
    implicit none
    integer , intent(in) :: istep
    real(rp), intent(in) :: time
    integer , intent(in), dimension(3) :: n
    integer , intent(in)               :: nh_psi,nh_g
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(1-nh_g:) :: dzf
    real(rp), intent(in) :: gam,seps
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
#if defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    real(rp), dimension(11) :: var
    real(rp) :: vol,clip,fmin,fmax,cfl(3)
    integer  :: nsub
    !
    ! volume of phase 1 and its relative change, accumulated volume added by clipping, extrema of the
    ! volume fraction before clipping, CFL numbers and number of sub-steps of the last time step
    !
    call vof_split_diagnostics(n,nh_psi,nh_g,dl,dzf,psi,vol,clip,fmin,fmax,cfl,nsub)
    var(1)    = 1.*istep
    var(2)    = time
    var(3)    = vol
    var(4)    = (vol-vol0)/max(vol0,tiny(1._rp))
    var(5)    = clip
    var(6)    = fmin
    var(7)    = fmax
    var(8:10) = cfl(:)
    var(11)   = 1.*nsub
    call out0d(trim(datadir)//'log_vof.out',11,var)
#elif !defined(_INTERFACE_CAPTURING_VOF)
    real(rp), dimension(4) :: var
    var(1) = 1.*istep
    var(2) = time
    var(3) = gam
    var(4) = seps
    call out0d(trim(datadir)//'log_acdi.out',4,var)
#endif
  end subroutine ifc_write_log
end module mod_interface_capturing
