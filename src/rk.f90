! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
#define _FAST_MOM_KERNELS
module mod_rk
  use mod_impdiff, only: impdiff_z_mom_rhs,impdiff_z_scal_rhs,impdiff_cc_mom_rhs,impdiff_cc_scal_rhs,impdiff_wall_force
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mod_impdiff, only: impdiff_vc_mom_rhs,impdiff_vc_scal_rhs
#endif
  use mod_mom    , only: mom_xyz_ad,mom_xyz_a_cons,mom_xyz_oth,cmpt_wallshear,cmpt_bodyforce
  use mod_param  , only: is_impdiff,is_impdiff_1d,is_vof_sweep_transport,visc0,alpha0
  use mod_types
  implicit none
  private
  public rk,rk_scal,rk_2fl
  contains
  subroutine rk(rkpar,n,nh_v,nh_p,nh_s,nh_psi,nh_ifc,nh_g,is_bound,is_forced,l,dli,dzci,dzfi,wzc, &
                grid_vol_ratio_c,grid_vol_ratio_f,dt,dt_r, &
                bforce,gacc,sigma,rho_av,rho12,mu12,beta12,rho0,psi,kappa,p,pn,po,s, &
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
                surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o, &
#endif
                psio,psiflx_x,psiflx_y,psiflx_z,is_first,wall_force_old,dudtrko,dvdtrko,dwdtrko,u,v,w,force_balance, &
                uadv,vadv,wadv)
    !
    ! AB2 scheme for time integration of the momentum equations
    !
    ! `is_first` is true in the first time step of a run without the terms of the previous step
    ! (`d?dtrko` and `wall_force_old`), which are then set to those of the present one
    !
    ! with the geometric volume-of-fluid method, the advection terms are integrated in one step, consistently
    ! with the transport of the volume fraction: they use its phase fluxes and its advecting velocity `uadv`,
    ! `vadv`, `wadv` (extrapolated to the middle of the time step), which is also the transported velocity;
    ! with `is_vof_sweep_transport`, the velocity has been transported along with the sweeps of the volume
    ! fraction instead, and `uadv`, `vadv`, `wadv` is the transported velocity
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_p,nh_s,nh_psi,nh_ifc,nh_g
    logical , intent(in), dimension(0:1,3) :: is_bound
    logical , intent(in), dimension(3) :: is_forced
    real(rp), intent(in), dimension(3) :: l,dli
    real(rp), intent(in) :: dt,dt_r
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi,wzc,grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(in   ), dimension(3)  :: bforce,gacc
    real(rp), intent(in   )                :: sigma,rho_av(3)
    real(rp), intent(in   ), dimension(2)  :: rho12,mu12,beta12
    real(rp), intent(in   )                :: rho0
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:)           :: psi
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:)           :: kappa
    real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:)           :: p
    real(rp), intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:), optional :: pn,po
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: s
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
    real(rp), intent(inout), dimension(:,:,:) :: surfx_n,surfy_n,surfz_n, &
                                                 surfx_o,surfy_o,surfz_o
#endif
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:), optional :: psio
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: psiflx_x,psiflx_y,psiflx_z
    logical , intent(in   )                                 :: is_first
    real(rp), intent(inout), dimension(3)                   :: wall_force_old
    real(rp), intent(inout), dimension(:,:,:)               :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:)           :: u,v,w
    real(rp), intent(out  ), dimension(3,6)                 :: force_balance
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:), optional :: uadv,vadv,wadv
    real(rp), allocatable, dimension(:,:,:), save :: dudtrk,dvdtrk,dwdtrk
    real(rp) :: rho,drho,rhox_n,rhox_p,rhoy_n,rhoy_p,rhoz_n,rhoz_p
    real(rp), save, dimension(3) :: bforce_eff
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(3) :: taux,tauy,tauz,wall_force,wall_force_imp,dwall_force_old
    integer :: i,j,k
    !
    factor1  = rkpar(1)*dt
    factor2  = rkpar(2)*dt
    factor12 = factor1+factor2
    !
    ! initialization
    !
    if(.not.allocated(dudtrk)) then ! leverage save attribute to allocate these arrays on the device only once
      allocate(dudtrk(n(1),n(2),n(3)),dvdtrk(n(1),n(2),n(3)),dwdtrk(n(1),n(2),n(3)))
      !$acc enter data create(dudtrk,dvdtrk,dwdtrk) async(1)
      !$acc enter data create(bforce_eff) async(1)
    end if
    !
    ! zero-net-acceleration forcing
    !
    force_balance(:,:) = 0._rp
    bforce_eff(:) = bforce(:)
    if(any(is_forced(:))) then
      call cmpt_wallshear(n,nh_v,nh_psi,nh_g,is_forced,is_bound,l,dli,dzci,dzfi,mu12,psi,u,v,w,taux,tauy,tauz)
      wall_force(1) = sum(taux(:)/l(:))
      wall_force(2) = sum(tauy(:)/l(:))
      wall_force(3) = sum(tauz(:)/l(:))
      if(is_first) wall_force_old(:) = wall_force(:)
      force_balance(:,1) = rkpar(1)*wall_force(:)+rkpar(2)*wall_force_old(:)
      call cmpt_bodyforce(n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,is_forced,dli,dzci,wzc,grid_vol_ratio_c,grid_vol_ratio_f,dt_r, &
                          rho12,beta12,gacc,sigma,rho0,rho_av,psi,psio,kappa,s,u,v,w, &
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
                          surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o, &
#endif
                          force_balance(:,2),force_balance(:,3),force_balance(:,4),force_balance(:,6))
      where(is_forced(:))
        bforce_eff(:) = force_balance(:,1)-force_balance(:,2)-force_balance(:,3)-force_balance(:,4)
      end where
      wall_force_old(:) = wall_force(:)
    end if
    force_balance(:,5) = bforce_eff(:)
    !$acc update device(bforce_eff) async(1)
    !
    ! advection, diffusion and regularization terms
    !
    call mom_xyz_ad(n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,wzc,rho12,mu12,psiflx_x,psiflx_y,psiflx_z,u,v,w,psi, &
                    dudtrk,dvdtrk,dwdtrk)
    if(is_impdiff) then
      !
      ! the diffusion along z that is integrated implicitly leaves the explicit terms (see `mod_impdiff`)
      !
      if(is_impdiff_1d) then
        call impdiff_z_mom_rhs(n,nh_v,nh_psi,nh_g,dzci,dzfi,mu12,psi,u,v,w,dudtrk,dvdtrk,dwdtrk)
      else
#if defined(_CONSTANT_COEFFS_POISSON)
        call impdiff_cc_mom_rhs(n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,wzc,visc0,rho12,psi,psio,u,v,w,dudtrk,dvdtrk,dwdtrk, &
                                dudtrko,dvdtrko,dwdtrko)
#else
        call impdiff_vc_mom_rhs(n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,mu12,psi,u,v,w,dudtrk,dvdtrk,dwdtrk)
#endif
      end if
      if(any(is_forced(:))) then
        !
        ! the forcing leaves out the wall shear of the terms integrated implicitly, which is balanced
        ! together with the implicit solve (see the main program)
        !
        call impdiff_wall_force(n,nh_psi,nh_g,is_forced,grid_vol_ratio_c,grid_vol_ratio_f,wzc,rho12,psi,psio, &
                                wall_force_imp,dwall_force_old)
        wall_force_old(:)  = wall_force_old(:) - wall_force_imp(:)
        wall_force_imp(:)  = merge(rkpar(1)+rkpar(2),rkpar(1),is_first)*wall_force_imp(:)
        if(.not.is_first) wall_force_imp(:) = wall_force_imp(:) - rkpar(2)*dwall_force_old(:)
        force_balance(:,1) = force_balance(:,1) - wall_force_imp(:)
        bforce_eff(:)      = bforce_eff(:)      - wall_force_imp(:)
        force_balance(:,5) = bforce_eff(:)
        !$acc update device(bforce_eff) async(1)
      end if
    end if
    !
    if(is_first) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dudtrko(i,j,k) = dudtrk(i,j,k)
            dvdtrko(i,j,k) = dvdtrk(i,j,k)
            dwdtrko(i,j,k) = dwdtrk(i,j,k)
          end do
        end do
      end do
    end if
    !
#if !defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          rhox_n = rho + drho*0.5*(psio(i,j,k)+psio(i+1,j,k))
          rhoy_n = rho + drho*0.5*(psio(i,j,k)+psio(i,j+1,k))
          rhoz_n = rho + drho*(wzc(k)*psio(i,j,k)+(1.-wzc(k))*psio(i,j,k+1)) ! volume average
          rhox_p = rho + drho*0.5*(psi( i,j,k)+psi( i+1,j,k))
          rhoy_p = rho + drho*0.5*(psi( i,j,k)+psi( i,j+1,k))
          rhoz_p = rho + drho*(wzc(k)*psi( i,j,k)+(1.-wzc(k))*psi( i,j,k+1))
          u(i,j,k) = (u(i,j,k)*rhox_n + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k))/rhox_p
          v(i,j,k) = (v(i,j,k)*rhoy_n + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k))/rhoy_p
          w(i,j,k) = (w(i,j,k)*rhoz_n + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k))/rhoz_p
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        end do
      end do
    end do
#else
    if(.not.is_vof_sweep_transport) then
      !
      ! the same update, together with the advection terms, in a single kernel
      !
      call mom_xyz_a_cons(n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,wzc,rho12,psi,psio,psiflx_x,psiflx_y,psiflx_z, &
                          uadv,vadv,wadv,factor1,factor2,factor12,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko,u,v,w)
    else
      !
      ! the other terms are added to the velocity transported along with the volume fraction
      !
      rho = rho12(2); drho = rho12(1)-rho12(2)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            rhox_p = rho + drho*0.5*(psi(i,j,k)+psi(i+1,j,k))
            rhoy_p = rho + drho*0.5*(psi(i,j,k)+psi(i,j+1,k))
            rhoz_p = rho + drho*(wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1))
            u(i,j,k) = uadv(i,j,k) + (factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k))/rhox_p
            v(i,j,k) = vadv(i,j,k) + (factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k))/rhoy_p
            w(i,j,k) = wadv(i,j,k) + (factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k))/rhoz_p
            dudtrko(i,j,k) = dudtrk(i,j,k)
            dvdtrko(i,j,k) = dvdtrk(i,j,k)
            dwdtrko(i,j,k) = dwdtrk(i,j,k)
          end do
        end do
      end do
    end if
#endif
    !
    ! pressure, surface tension, and buoyancy terms
    !
    call mom_xyz_oth(n,nh_v,nh_p,nh_s,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,wzc,dt_r,rho12,beta12,bforce_eff,gacc,sigma,rho0,rho_av, &
                     p,psi,kappa,s,pn,po, &
#if defined(_BALANCED_CAPILLARY_PRESSURE_SPLIT)
                     surfx_n,surfy_n,surfz_n,surfx_o,surfy_o,surfz_o, &
#endif
                     factor12,u,v,w)
  end subroutine rk
  !
  subroutine rk_scal(rkpar,n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,dt,ssource,rho12,ka12,cp12,psi,u,v,w, &
                     psio,psiflx_x,psiflx_y,psiflx_z,is_first,dsdtrko,s,so,sadv,uadv,vadv,wadv)
    use mod_scal, only: scal_ad,scal_a_cons
    !
    ! AB2 scheme for time integration of the scalar field
    !
    ! with the geometric volume-of-fluid method, the advection term is integrated in one step, consistently
    ! with the transport of the volume fraction (see `rk`), with the transported scalar `sadv` extrapolated
    ! to the middle of the time step from the present scalar and that of the previous time step, `so`;
    ! with `is_vof_sweep_transport`, the scalar has been transported along with the sweeps of the volume
    ! fraction instead, and `sadv` is the transported scalar
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in   )                :: dt
    real(rp), intent(in   )                :: ssource
    real(rp), intent(in   ), dimension(2)  :: rho12,ka12,cp12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psio
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    logical , intent(in   )                   :: is_first
    real(rp), intent(inout), dimension(:,:,:) :: dsdtrko
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:), optional :: so,sadv
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:), optional :: uadv,vadv,wadv
#if defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    real(rp) :: c1,c2
#endif
    real(rp), allocatable, dimension(:,:,:), save :: dsdtrk
    real(rp), dimension(2), save :: rhocp12
    real(rp) :: rhocp,drhocp,rhocp_n,rhocp_p
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    if(.not.allocated(dsdtrk)) then ! leverage save attribute to allocate these arrays on the device only once
      allocate(dsdtrk(n(1),n(2),n(3)))
      !$acc enter data create(dsdtrk) async(1)
      rhocp12(:) = rho12(:)*cp12(:)
      !$acc wait(1)
    end if
#if defined(_INTERFACE_CAPTURING_SPLIT_VOF)
    if(.not.is_vof_sweep_transport) then
      c1 = rkpar(1)
      c2 = rkpar(2)
      if(is_first) then
        c1 = 1.
        c2 = 0.
      end if
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1-nh_s,n(3)+nh_s
        do j=1-nh_s,n(2)+nh_s
          do i=1-nh_s,n(1)+nh_s
            sadv(i,j,k) = c1*s(i,j,k)+c2*so(i,j,k)
            so(i,j,k) = s(i,j,k)
          end do
        end do
      end do
    end if
#endif
    call scal_ad(n(1),n(2),n(3),nh_v,nh_s,nh_psi,nh_ifc,nh_g,dli(1),dli(2),dzci,dzfi,ka12,rhocp12,psiflx_x,psiflx_y,psiflx_z, &
                 psi,u,v,w,s,dsdtrk)
    if(is_impdiff) then
      !
      ! the diffusion along z that is integrated implicitly leaves the explicit terms (see `mod_impdiff`)
      !
      if(is_impdiff_1d) then
        call impdiff_z_scal_rhs(n,nh_s,nh_psi,nh_g,dzci,dzfi,ka12,psi,s,dsdtrk)
      else
#if defined(_CONSTANT_COEFFS_POISSON)
        call impdiff_cc_scal_rhs(n,nh_s,nh_psi,nh_g,dli,dzci,dzfi,alpha0,rhocp12,psi,psio,s,dsdtrk,dsdtrko)
#else
        call impdiff_vc_scal_rhs(n,nh_s,nh_psi,nh_g,dli,dzci,dzfi,ka12,psi,s,dsdtrk)
#endif
      end if
    end if
    if(is_first) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dsdtrko(i,j,k) = dsdtrk(i,j,k)
          end do
        end do
      end do
    end if
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    if(.not.is_vof_sweep_transport) then
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            rhocp_p = rhocp + drhocp*psi( i,j,k)
            rhocp_n = rhocp + drhocp*psio(i,j,k)
            s(i,j,k) = s(i,j,k)*rhocp_n/rhocp_p + (factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource)/rhocp_p
            dsdtrko(i,j,k) = dsdtrk(i,j,k)
          end do
        end do
      end do
#if defined(_INTERFACE_CAPTURING_SPLIT_VOF)
      call scal_a_cons(n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dli,dzfi,rhocp12,psi,psiflx_x,psiflx_y,psiflx_z,uadv,vadv,wadv,sadv, &
                       factor12,s)
#endif
    else
      !
      ! the other terms are added to the scalar transported along with the volume fraction
      !
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            rhocp_p = rhocp + drhocp*psi(i,j,k)
            s(i,j,k) = sadv(i,j,k) + (factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k) + factor12*ssource)/rhocp_p
            dsdtrko(i,j,k) = dsdtrk(i,j,k)
          end do
        end do
      end do
    end if
  end subroutine rk_scal
  !
  subroutine rk_2fl(rkpar,n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,dt,gam,seps,beta,u,v,w,normx,normy,normz,phi,psi, &
                    psiflx_x,psiflx_y,psiflx_z,is_first,dpsidtrko)
#if !defined(_INTERFACE_CAPTURING_VOF)
    use mod_acdi        , only: acdi_transport_pf
#else
    use mod_vof_thinc_qq, only: vof_thinc_transport_psi
#endif
    use mod_two_fluid   , only: clip_field
    !
    ! AB2 scheme for time integration of the phase field
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in) :: nh_v,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in   ) :: gam,seps,beta,dt
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: normx,normy,normz
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: phi
    real(rp), intent(inout), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(out  ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    logical , intent(in   )                   :: is_first
    real(rp), intent(inout), dimension(:,:,:) :: dpsidtrko
    real(rp), allocatable, dimension(:,:,:), save :: dpsidtrk
    real(rp) :: factor1,factor2,factor12
    integer :: i,j,k
    !
    factor1  = rkpar(1)*dt
    factor2  = rkpar(2)*dt
    factor12 = factor1 + factor2 ! unused but left here in case a source term is needed in the future
    if(.not.allocated(dpsidtrk)) then ! leverage save attribute to allocate these arrays on the device only once
      allocate(dpsidtrk(n(1),n(2),n(3)))
      !$acc enter data create(dpsidtrk) async(1)
    end if
#if !defined(_INTERFACE_CAPTURING_VOF)
    call acdi_transport_pf(n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,gam,seps,u,v,w,normx,normy,normz,phi,psi,dpsidtrk, &
                           psiflx_x,psiflx_y,psiflx_z)
#else
    call vof_thinc_transport_psi(n,nh_v,nh_psi,nh_ifc,nh_g,dli,dzfi,beta,u,v,w,normx,normy,normz,phi,psi,dpsidtrk, &
                                 psiflx_x,psiflx_y,psiflx_z)
#endif
    if(is_first) then
      !$acc parallel loop collapse(3) default(present) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            dpsidtrko(i,j,k) = dpsidtrk(i,j,k)
          end do
        end do
      end do
    end if
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          psi(i,j,k) = psi(i,j,k) + factor1*dpsidtrk(i,j,k) + factor2*dpsidtrko(i,j,k)
          dpsidtrko(i,j,k) = dpsidtrk(i,j,k)
        end do
      end do
    end do
    call clip_field([nh_psi,nh_psi,nh_psi],[0._rp,1._rp],psi)
  end subroutine rk_2fl
end module mod_rk
