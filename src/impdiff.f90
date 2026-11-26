! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_impdiff
  !
  ! implicit integration of the diffusion terms (Crank-Nicolson), in three variants:
  !
  ! - along z (`is_impdiff_1d = T`, routines `impdiff_z_*`): the part of the viscous term that is integrated
  !   implicitly is the diffusion of each velocity component along z due to its own gradient,
  !     d/dz(mu du/dz), d/dz(mu dv/dz) and d/dz(2 mu dw/dz),
  !   with the local viscosity and density, i.e., one tridiagonal system with variable coefficients per line
  !   along z; the remaining viscous terms are integrated explicitly; likewise for the scalar, d/dz(ka ds/dz)
  ! - in all directions with constant coefficients (`is_impdiff = T` with the FFT-based pressure solver, routines
  !   `impdiff_cc_*`): the Laplacian of each velocity component, and of the scalar, times a reference diffusivity,
  !   solved with the Helmholtz solver of CaNS; the remainder of the diffusion terms is integrated explicitly
  ! - in all directions with variable coefficients (`is_impdiff = T` with the HYPRE solver, routines `impdiff_vc_*`):
  !   the diffusion of each velocity component due to its own gradient, and that of the scalar, with the local
  !   properties
  !
  ! the viscosity limiter (`is_visc_limiter`) only acts on the terms that are integrated explicitly
  !
  use mpi
  use mod_param    , only: is_impdiff_1d,is_harmonic_interpolation
#if !defined(_CONSTANT_COEFFS_POISSON)
  use mod_solver_vc, only: solver_vc_helmholtz
#endif
  use mod_types
  implicit none
  private
  public impdiff_z_mom_rhs,impdiff_z_mom_solve,impdiff_z_scal_rhs,impdiff_z_scal_solve, &
         impdiff_cc_mom_rhs,impdiff_cc_mom_updt,impdiff_cc_scal_rhs,impdiff_cc_scal_updt, &
#if !defined(_CONSTANT_COEFFS_POISSON)
         impdiff_vc_mom_rhs,impdiff_vc_mom_updt,impdiff_vc_mom_solve, &
         impdiff_vc_scal_rhs,impdiff_vc_scal_updt,impdiff_vc_scal_solve, &
#endif
         impdiff_wall_force,impdiff_updt_forcing,impdiff_wall_heat,impdiff_updt_forcing_s
  real(rp), allocatable, dimension(:,:,:), save :: cw            ! work array of the line solver
  real(rp), allocatable, dimension(:,:,:), save :: dud,dvd,dwd,dsd ! implicit terms at the old time level
  contains
  subroutine impdiff_z_mom_rhs(n,nh_v,nh_psi,nh_g,dzci,dzfi,mu12,psi,u,v,w,dudt,dvdt,dwdt)
    !
    ! computes the terms that are integrated implicitly, discretized as in `mom_xyz_ad` (without its viscosity
    ! limiter),
    ! takes them out of the terms that are integrated explicitly (`dudt`, `dvdt`, `dwdt`),
    ! and keeps them for the explicit half of the Crank-Nicolson scheme (see `impdiff_z_mom_solve`)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(2) :: mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(inout), dimension(:,:,:) :: dudt,dvdt,dwdt
    real(rp) :: mu,dmu,muzp,muzm
    logical  :: is_harm
    integer  :: i,j,k
    !
    if(.not.allocated(dud)) then
      allocate(dud(n(1),n(2),n(3)),dvd(n(1),n(2),n(3)),dwd(n(1),n(2),n(3)))
      !$acc        enter data create(dud,dvd,dwd) async(1)
      !$omp target enter data map(alloc:dud,dvd,dwd)
    end if
    call visc_coeffs(mu12,mu,dmu,is_harm)
    !$acc parallel     loop collapse(3) default(present) private(muzp,muzm) async(1)
    !$omp target teams loop collapse(3)                  private(muzp,muzm)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          muzp = visc(is_harm,mu,dmu,0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)))
          muzm = visc(is_harm,mu,dmu,0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)))
          dud(i,j,k) = dzfi(k)*((u(i,j,k+1)-u(i,j,k))*dzci(k)*muzp-(u(i,j,k)-u(i,j,k-1))*dzci(k-1)*muzm)
          muzp = visc(is_harm,mu,dmu,0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1)))
          muzm = visc(is_harm,mu,dmu,0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1)))
          dvd(i,j,k) = dzfi(k)*((v(i,j,k+1)-v(i,j,k))*dzci(k)*muzp-(v(i,j,k)-v(i,j,k-1))*dzci(k-1)*muzm)
          muzp = visc(is_harm,mu,dmu,psi(i,j,k+1))
          muzm = visc(is_harm,mu,dmu,psi(i,j,k  ))
          dwd(i,j,k) = dzci(k)*(2.*(w(i,j,k+1)-w(i,j,k))*dzfi(k+1)*muzp-2.*(w(i,j,k)-w(i,j,k-1))*dzfi(k)*muzm)
          dudt(i,j,k) = dudt(i,j,k) - dud(i,j,k)
          dvdt(i,j,k) = dvdt(i,j,k) - dvd(i,j,k)
          dwdt(i,j,k) = dwdt(i,j,k) - dwd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_z_mom_rhs
  !
  subroutine impdiff_z_mom_solve(n,nh_v,nh_psi,nh_g,cbc,bc,dzc,dzf,dzci,dzfi,wzc,theta,rho12,mu12,psi,u,v,w)
    !
    ! solves, for each velocity component u_q and each line along z,
    !   rho_q u_q - theta*d/dz(c_q mu du_q/dz) = rho_q u_q* + theta*d/dz(c_q mu du_q/dz)^n,
    ! where u_q* is the input velocity, rho_q the density at the location of u_q, c_q = 2 for w and 1
    ! otherwise, theta half of the time step (Crank-Nicolson), and the last term the one that
    ! `impdiff_z_mom_rhs` computed at the old time level;
    ! the boundary conditions along z are those that `bounduvw` imposes
    !
    ! each line is solved with the Thomas algorithm; the domain must not be decomposed along z
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    character(len=1), intent(in   ), dimension(0:1,3,3) :: cbc
    real(rp), intent(in   ), dimension(0:1,3,3)      :: bc
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi,wzc
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    !
    if(.not.allocated(cw)) then
      allocate(cw(n(1),n(2),n(3)))
      !$acc        enter data create(cw) async(1)
      !$omp target enter data map(alloc:cw)
    end if
    call solve_lines_c(n,nh_v,nh_psi,nh_g,1,0,cbc(:,3,1),bc(:,3,1),dzc,dzci,dzfi,theta,rho12,mu12,psi,dud,cw,u)
    call solve_lines_c(n,nh_v,nh_psi,nh_g,0,1,cbc(:,3,2),bc(:,3,2),dzc,dzci,dzfi,theta,rho12,mu12,psi,dvd,cw,v)
    call solve_lines_f(n,nh_v,nh_psi,nh_g    ,cbc(:,3,3),bc(:,3,3),dzf,dzci,dzfi,wzc,theta,rho12,mu12,psi,dwd,cw,w)
  end subroutine impdiff_z_mom_solve
  !
  subroutine impdiff_cc_mom_rhs(n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,wzc,visc0,rho12,psi,psio,u,v,w,dudt,dvdt,dwdt, &
                                dudto,dvdto,dwdto)
    !
    ! constant-coefficient implicit diffusion in all directions: the term that is integrated implicitly is
    ! the Laplacian of each velocity component times a reference kinematic viscosity, visc0*lap(u_q);
    ! it is computed here, taken out of the terms that are integrated explicitly (times the local density,
    ! since `dudt`, `dvdt` and `dwdt` are per unit volume), and kept for the explicit half of the
    ! Crank-Nicolson scheme (see `impdiff_cc_mom_updt`);
    ! the explicit terms of the previous step (`dudto`, `dvdto`, `dwdto`) hold this term times the density
    ! of that step: they are updated to the present density, which keeps the scheme second-order accurate
    ! where the density changes in time
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi,wzc
    real(rp), intent(in   )               :: visc0
    real(rp), intent(in   ), dimension(2) :: rho12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,psio
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(inout), dimension(:,:,:) :: dudt,dvdt,dwdt,dudto,dvdto,dwdto
    real(rp) :: rho,drho,dxi2,dyi2
    integer  :: i,j,k
    !
    if(.not.allocated(dud)) then
      allocate(dud(n(1),n(2),n(3)),dvd(n(1),n(2),n(3)),dwd(n(1),n(2),n(3)))
      !$acc        enter data create(dud,dvd,dwd) async(1)
      !$omp target enter data map(alloc:dud,dvd,dwd)
    end if
    rho = rho12(2); drho = rho12(1)-rho12(2)
    dxi2 = dli(1)**2; dyi2 = dli(2)**2
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dud(i,j,k) = visc0*((u(i+1,j,k)-2.*u(i,j,k)+u(i-1,j,k))*dxi2 + &
                              (u(i,j+1,k)-2.*u(i,j,k)+u(i,j-1,k))*dyi2 + &
                              ((u(i,j,k+1)-u(i,j,k))*dzci(k)-(u(i,j,k)-u(i,j,k-1))*dzci(k-1))*dzfi(k))
          dvd(i,j,k) = visc0*((v(i+1,j,k)-2.*v(i,j,k)+v(i-1,j,k))*dxi2 + &
                              (v(i,j+1,k)-2.*v(i,j,k)+v(i,j-1,k))*dyi2 + &
                              ((v(i,j,k+1)-v(i,j,k))*dzci(k)-(v(i,j,k)-v(i,j,k-1))*dzci(k-1))*dzfi(k))
          dwd(i,j,k) = visc0*((w(i+1,j,k)-2.*w(i,j,k)+w(i-1,j,k))*dxi2 + &
                              (w(i,j+1,k)-2.*w(i,j,k)+w(i,j-1,k))*dyi2 + &
                              ((w(i,j,k+1)-w(i,j,k))*dzfi(k+1)-(w(i,j,k)-w(i,j,k-1))*dzfi(k))*dzci(k))
          dudt(i,j,k) = dudt(i,j,k) - (rho + drho*0.5*(psi(i,j,k)+psi(i+1,j,k)))*dud(i,j,k)
          dvdt(i,j,k) = dvdt(i,j,k) - (rho + drho*0.5*(psi(i,j,k)+psi(i,j+1,k)))*dvd(i,j,k)
          dwdt(i,j,k) = dwdt(i,j,k) - (rho + drho*(wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1)))*dwd(i,j,k)
          dudto(i,j,k) = dudto(i,j,k) + drho*0.5*(psio(i,j,k)+psio(i+1,j,k)-psi(i,j,k)-psi(i+1,j,k))*dud(i,j,k)
          dvdto(i,j,k) = dvdto(i,j,k) + drho*0.5*(psio(i,j,k)+psio(i,j+1,k)-psi(i,j,k)-psi(i,j+1,k))*dvd(i,j,k)
          dwdto(i,j,k) = dwdto(i,j,k) + drho*(wzc(k)*psio(i,j,k)+(1.-wzc(k))*psio(i,j,k+1) - &
                                              wzc(k)*psi( i,j,k)-(1.-wzc(k))*psi( i,j,k+1))*dwd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_cc_mom_rhs
  !
  subroutine impdiff_cc_mom_updt(n,nh_v,theta,u,v,w)
    !
    ! adds the explicit half of the Crank-Nicolson integration of the constant-coefficient diffusion,
    ! with `theta` half of the time step; the implicit half follows, with the Helmholtz solver
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v
    real(rp), intent(in   )               :: theta
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    integer  :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + theta*dud(i,j,k)
          v(i,j,k) = v(i,j,k) + theta*dvd(i,j,k)
          w(i,j,k) = w(i,j,k) + theta*dwd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_cc_mom_updt
  !
  subroutine impdiff_cc_scal_rhs(n,nh_s,nh_psi,nh_g,dli,dzci,dzfi,alpha0,rhocp12,psi,psio,s,dsdt,dsdto)
    !
    ! constant-coefficient implicit diffusion of the scalar, alpha0*lap(s), with a reference diffusivity
    ! (see `impdiff_cc_mom_rhs`; the role of the density is played by the volumetric heat capacity)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_s,nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in   )               :: alpha0
    real(rp), intent(in   ), dimension(2) :: rhocp12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,psio
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(inout), dimension(:,:,:) :: dsdt,dsdto
    real(rp) :: rhocp,drhocp,dxi2,dyi2
    integer  :: i,j,k
    !
    if(.not.allocated(dsd)) then
      allocate(dsd(n(1),n(2),n(3)))
      !$acc        enter data create(dsd) async(1)
      !$omp target enter data map(alloc:dsd)
    end if
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    dxi2 = dli(1)**2; dyi2 = dli(2)**2
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dsd(i,j,k) = alpha0*((s(i+1,j,k)-2.*s(i,j,k)+s(i-1,j,k))*dxi2 + &
                               (s(i,j+1,k)-2.*s(i,j,k)+s(i,j-1,k))*dyi2 + &
                               ((s(i,j,k+1)-s(i,j,k))*dzci(k)-(s(i,j,k)-s(i,j,k-1))*dzci(k-1))*dzfi(k))
          dsdt(i,j,k)  = dsdt(i,j,k)  - (rhocp + drhocp*psi(i,j,k))*dsd(i,j,k)
          dsdto(i,j,k) = dsdto(i,j,k) + drhocp*(psio(i,j,k)-psi(i,j,k))*dsd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_cc_scal_rhs
  !
  subroutine impdiff_cc_scal_updt(n,nh_s,theta,s)
    !
    ! adds the explicit half of the Crank-Nicolson integration of the constant-coefficient diffusion
    ! of the scalar (see `impdiff_cc_mom_updt`)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_s
    real(rp), intent(in   )               :: theta
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    integer  :: i,j,k
    !
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + theta*dsd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_cc_scal_updt
  !
#if !defined(_CONSTANT_COEFFS_POISSON)
  subroutine impdiff_vc_mom_rhs(n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,mu12,psi,u,v,w,dudt,dvdt,dwdt)
    !
    ! variable-coefficient implicit diffusion in all directions: the terms that are integrated implicitly
    ! are the diffusion of each velocity component due to its own gradient,
    !   div(mu grad u_q) + d/dx_q(mu du_q/dx_q),
    ! discretized as in `mom_xyz_ad` (without its viscosity limiter); they are computed here, taken out of the terms that are
    ! integrated explicitly, and kept for the explicit half of the Crank-Nicolson scheme
    ! (see `impdiff_vc_mom_updt`)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(2) :: mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(inout), dimension(:,:,:) :: dudt,dvdt,dwdt
    integer  :: i,j,k
    !
    if(.not.allocated(dud)) allocate(dud(n(1),n(2),n(3)),dvd(n(1),n(2),n(3)),dwd(n(1),n(2),n(3)))
    call impdiff_vc_mom_cmpt(n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,mu12,psi,u,v,w,dud,dvd,dwd)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dudt(i,j,k) = dudt(i,j,k) - dud(i,j,k)
          dvdt(i,j,k) = dvdt(i,j,k) - dvd(i,j,k)
          dwdt(i,j,k) = dwdt(i,j,k) - dwd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_vc_mom_rhs
  !
  subroutine impdiff_vc_mom_cmpt(n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,mu12,psi,u,v,w,bu,bv,bw)
    !
    ! computes the diffusion terms of `impdiff_vc_mom_rhs` for the velocity field `u`, `v`, `w`
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: mu12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(out), dimension(:,:,:) :: bu,bv,bw
    real(rp) :: mu,dmu,dxi2,dyi2,muxp,muxm,muyp,muym,muzp,muzm
    logical  :: is_harm
    real(rp) :: c_ccc,c_pcc,c_mcc,c_cpc,c_cmc,c_ccp,c_ccm,c_ppc,c_pmc,c_mpc,c_pcp,c_pcm,c_mcp,c_cpp,c_cpm,c_cmp
    integer  :: i,j,k
    !
    call visc_coeffs(mu12,mu,dmu,is_harm)
    dxi2 = dli(1)**2; dyi2 = dli(2)**2
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          c_ccc = psi(i  ,j  ,k  )
          c_pcc = psi(i+1,j  ,k  )
          c_mcc = psi(i-1,j  ,k  )
          c_cpc = psi(i  ,j+1,k  )
          c_cmc = psi(i  ,j-1,k  )
          c_ccp = psi(i  ,j  ,k+1)
          c_ccm = psi(i  ,j  ,k-1)
          c_ppc = psi(i+1,j+1,k  )
          c_pmc = psi(i+1,j-1,k  )
          c_mpc = psi(i-1,j+1,k  )
          c_pcp = psi(i+1,j  ,k+1)
          c_pcm = psi(i+1,j  ,k-1)
          c_mcp = psi(i-1,j  ,k+1)
          c_cpp = psi(i  ,j+1,k+1)
          c_cpm = psi(i  ,j+1,k-1)
          c_cmp = psi(i  ,j-1,k+1)
          !
          muxp = visc(is_harm,mu,dmu,c_pcc)
          muxm = visc(is_harm,mu,dmu,c_ccc)
          muyp = visc(is_harm,mu,dmu,0.25*(c_ccc+c_cpc+c_ppc+c_pcc))
          muym = visc(is_harm,mu,dmu,0.25*(c_ccc+c_cmc+c_pmc+c_pcc))
          muzp = visc(is_harm,mu,dmu,0.25*(c_ccc+c_pcc+c_ccp+c_pcp))
          muzm = visc(is_harm,mu,dmu,0.25*(c_ccc+c_pcc+c_ccm+c_pcm))
          bu(i,j,k) = 2.*dxi2*((u(i+1,j,k)-u(i,j,k))*muxp-(u(i,j,k)-u(i-1,j,k))*muxm) + &
                         dyi2*((u(i,j+1,k)-u(i,j,k))*muyp-(u(i,j,k)-u(i,j-1,k))*muym) + &
                      dzfi(k)*((u(i,j,k+1)-u(i,j,k))*dzci(k)*muzp-(u(i,j,k)-u(i,j,k-1))*dzci(k-1)*muzm)
          !
          muxp = visc(is_harm,mu,dmu,0.25*(c_ccc+c_pcc+c_ppc+c_cpc))
          muxm = visc(is_harm,mu,dmu,0.25*(c_ccc+c_mcc+c_mpc+c_cpc))
          muyp = visc(is_harm,mu,dmu,c_cpc)
          muym = visc(is_harm,mu,dmu,c_ccc)
          muzp = visc(is_harm,mu,dmu,0.25*(c_ccc+c_cpc+c_cpp+c_ccp))
          muzm = visc(is_harm,mu,dmu,0.25*(c_ccc+c_cpc+c_cpm+c_ccm))
          bv(i,j,k) =    dxi2*((v(i+1,j,k)-v(i,j,k))*muxp-(v(i,j,k)-v(i-1,j,k))*muxm) + &
                      2.*dyi2*((v(i,j+1,k)-v(i,j,k))*muyp-(v(i,j,k)-v(i,j-1,k))*muym) + &
                      dzfi(k)*((v(i,j,k+1)-v(i,j,k))*dzci(k)*muzp-(v(i,j,k)-v(i,j,k-1))*dzci(k-1)*muzm)
          !
          muxp = visc(is_harm,mu,dmu,0.25*(c_ccc+c_pcc+c_ccp+c_pcp))
          muxm = visc(is_harm,mu,dmu,0.25*(c_ccc+c_mcc+c_ccp+c_mcp))
          muyp = visc(is_harm,mu,dmu,0.25*(c_ccc+c_cpc+c_ccp+c_cpp))
          muym = visc(is_harm,mu,dmu,0.25*(c_ccc+c_cmc+c_ccp+c_cmp))
          muzp = visc(is_harm,mu,dmu,c_ccp)
          muzm = visc(is_harm,mu,dmu,c_ccc)
          bw(i,j,k) =    dxi2*((w(i+1,j,k)-w(i,j,k))*muxp-(w(i,j,k)-w(i-1,j,k))*muxm) + &
                         dyi2*((w(i,j+1,k)-w(i,j,k))*muyp-(w(i,j,k)-w(i,j-1,k))*muym) + &
                      2.*dzci(k)*((w(i,j,k+1)-w(i,j,k))*dzfi(k+1)*muzp-(w(i,j,k)-w(i,j,k-1))*dzfi(k)*muzm)
        end do
      end do
    end do
  end subroutine impdiff_vc_mom_cmpt
  !
  subroutine impdiff_vc_mom_updt(n,nh_v,nh_psi,nh_g,wzc,theta,rho12,psi,u,v,w)
    !
    ! adds the explicit half of the Crank-Nicolson integration of the terms of `impdiff_vc_mom_rhs`,
    ! with `theta` half of the time step
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in   ), dimension(1-nh_g:) :: wzc
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp) :: rho,drho
    integer  :: i,j,k
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + theta*dud(i,j,k)/(rho + drho*0.5*(psi(i,j,k)+psi(i+1,j,k)))
          v(i,j,k) = v(i,j,k) + theta*dvd(i,j,k)/(rho + drho*0.5*(psi(i,j,k)+psi(i,j+1,k)))
          w(i,j,k) = w(i,j,k) + theta*dwd(i,j,k)/(rho + drho*(wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1)))
        end do
      end do
    end do
  end subroutine impdiff_vc_mom_updt
  !
  subroutine impdiff_vc_mom_solve(ng,lo,hi,n,nh_v,nh_psi,nh_g,cbc,cbcpre,is_bound,dli,dzci,dzfi,wzc,theta,rho12,mu12,psi, &
                                  u,v,w)
    !
    ! solves for the implicit half of the Crank-Nicolson integration, in correction form: with u_q the input
    ! velocity (which includes the explicit half, and whose halos and boundary values must be up to date),
    !   rho_q du_q - theta*B_q(du_q) = theta*B_q(u_q),   u_q <- u_q + du_q,
    ! where B_q is the diffusion operator of `impdiff_vc_mom_rhs`, and the correction du_q satisfies the
    ! homogeneous counterpart of the boundary conditions of u_q; the three systems are solved with HYPRE
    ! (second entry of the `hypre_*` inputs), with a tolerance relative to the size of the viscous update
    !
    implicit none
    integer , intent(in   ), dimension(3) :: ng,lo,hi,n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    character(len=1), intent(in   ), dimension(0:1,3,3) :: cbc
    character(len=1), intent(in   ), dimension(0:1,3)   :: cbcpre
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi,wzc
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), allocatable, dimension(:,:,:,:) :: coeffs
    integer :: i,j,k
    !
    call impdiff_vc_mom_cmpt(n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,mu12,psi,u,v,w,dud,dvd,dwd)
    allocate(coeffs(0:6,n(1),n(2),n(3)))
    call helmholtz_coeffs(n,nh_psi,nh_g,1,cbc(:,:,1),is_bound,dli,dzci,dzfi,wzc,theta,rho12,mu12,psi,coeffs,dud)
    call solver_vc_helmholtz(ng,lo,hi,cbcpre,coeffs,dud)
    call helmholtz_coeffs(n,nh_psi,nh_g,2,cbc(:,:,2),is_bound,dli,dzci,dzfi,wzc,theta,rho12,mu12,psi,coeffs,dvd)
    call solver_vc_helmholtz(ng,lo,hi,cbcpre,coeffs,dvd)
    call helmholtz_coeffs(n,nh_psi,nh_g,3,cbc(:,:,3),is_bound,dli,dzci,dzfi,wzc,theta,rho12,mu12,psi,coeffs,dwd)
    call solver_vc_helmholtz(ng,lo,hi,cbcpre,coeffs,dwd)
    deallocate(coeffs)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + dud(i,j,k)
          v(i,j,k) = v(i,j,k) + dvd(i,j,k)
          w(i,j,k) = w(i,j,k) + dwd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_vc_mom_solve
  !
  subroutine impdiff_vc_scal_rhs(n,nh_s,nh_psi,nh_g,dli,dzci,dzfi,ka12,psi,s,dsdt)
    !
    ! variable-coefficient implicit diffusion of the scalar in all directions, div(ka grad s), discretized
    ! exactly as in `scal_ad` (see `impdiff_vc_mom_rhs`)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_s,nh_psi,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(2) :: ka12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(inout), dimension(:,:,:) :: dsdt
    integer  :: i,j,k
    !
    if(.not.allocated(dsd)) allocate(dsd(n(1),n(2),n(3)))
    call impdiff_vc_scal_cmpt(n,nh_s,nh_psi,nh_g,dli,dzci,dzfi,ka12,psi,s,dsd)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dsdt(i,j,k) = dsdt(i,j,k) - dsd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_vc_scal_rhs
  !
  subroutine impdiff_vc_scal_cmpt(n,nh_s,nh_psi,nh_g,dli,dzci,dzfi,ka12,psi,s,bs)
    !
    ! computes the diffusion term of `impdiff_vc_scal_rhs` for the scalar field `s`
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_s,nh_psi,nh_g
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: ka12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(out), dimension(:,:,:) :: bs
    real(rp) :: ka,dka,dxi2,dyi2,kaxp,kaxm,kayp,kaym,kazp,kazm
    logical  :: is_harm
    integer  :: i,j,k
    !
    call visc_coeffs(ka12,ka,dka,is_harm)
    dxi2 = dli(1)**2; dyi2 = dli(2)**2
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          kaxp = visc(is_harm,ka,dka,0.5*(psi(i+1,j,k)+psi(i,j,k)))
          kaxm = visc(is_harm,ka,dka,0.5*(psi(i-1,j,k)+psi(i,j,k)))
          kayp = visc(is_harm,ka,dka,0.5*(psi(i,j+1,k)+psi(i,j,k)))
          kaym = visc(is_harm,ka,dka,0.5*(psi(i,j-1,k)+psi(i,j,k)))
          kazp = visc(is_harm,ka,dka,0.5*(psi(i,j,k+1)+psi(i,j,k)))
          kazm = visc(is_harm,ka,dka,0.5*(psi(i,j,k-1)+psi(i,j,k)))
          bs(i,j,k) = dxi2*((s(i+1,j,k)-s(i,j,k))*kaxp-(s(i,j,k)-s(i-1,j,k))*kaxm) + &
                      dyi2*((s(i,j+1,k)-s(i,j,k))*kayp-(s(i,j,k)-s(i,j-1,k))*kaym) + &
                      dzfi(k)*((s(i,j,k+1)-s(i,j,k))*dzci(k)*kazp-(s(i,j,k)-s(i,j,k-1))*dzci(k-1)*kazm)
        end do
      end do
    end do
  end subroutine impdiff_vc_scal_cmpt
  !
  subroutine impdiff_vc_scal_updt(n,nh_s,nh_psi,theta,rhocp12,psi,s)
    !
    ! adds the explicit half of the Crank-Nicolson integration of the term of `impdiff_vc_scal_rhs`
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_s,nh_psi
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rhocp12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp) :: rhocp,drhocp
    integer  :: i,j,k
    !
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + theta*dsd(i,j,k)/(rhocp + drhocp*psi(i,j,k))
        end do
      end do
    end do
  end subroutine impdiff_vc_scal_updt
  !
  subroutine impdiff_vc_scal_solve(ng,lo,hi,n,nh_s,nh_psi,nh_g,cbc,cbcpre,is_bound,dli,dzci,dzfi,theta,rhocp12,ka12,psi,s)
    !
    ! solves for the implicit half of the Crank-Nicolson integration of the scalar diffusion, in correction
    ! form (see `impdiff_vc_mom_solve`); the halos and boundary values of `s` must be up to date
    !
    implicit none
    integer , intent(in   ), dimension(3) :: ng,lo,hi,n
    integer , intent(in   )               :: nh_s,nh_psi,nh_g
    character(len=1), intent(in   ), dimension(0:1,3) :: cbc,cbcpre
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rhocp12,ka12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), allocatable, dimension(:,:,:,:) :: coeffs
    integer :: i,j,k
    !
    call impdiff_vc_scal_cmpt(n,nh_s,nh_psi,nh_g,dli,dzci,dzfi,ka12,psi,s,dsd)
    allocate(coeffs(0:6,n(1),n(2),n(3)))
    call helmholtz_coeffs(n,nh_psi,nh_g,0,cbc,is_bound,dli,dzci,dzfi,dzci,theta,rhocp12,ka12,psi,coeffs,dsd) ! (no weights needed)
    call solver_vc_helmholtz(ng,lo,hi,cbcpre,coeffs,dsd)
    deallocate(coeffs)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + dsd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_vc_scal_solve
  !
  subroutine helmholtz_coeffs(n,nh_psi,nh_g,iq,cbc,is_bound,dli,dzci,dzfi,wzc,theta,rho12,mu12,psi,coeffs,rhs)
    !
    ! matrix and right-hand side of the system of `impdiff_vc_mom_solve` for a variable staggered along
    ! direction `iq` (0: cell-centered, with `rho12` and `mu12` its storage and diffusion coefficients):
    !   coeffs(0)*dp_c + sum_nb coeffs(nb)*dp_nb = rhs,
    ! with the lower and upper neighbours along direction `idir` at entries 2*idir-1 and 2*idir;
    ! each equation is multiplied by the height of its control volume, which makes the matrix symmetric;
    ! on input `rhs` holds the diffusion term of the present field; `wzc` are the weights of the volume average
    ! at the z faces, used for the density of a variable staggered along z
    !
    ! boundaries: the correction of a variable centered along the direction mirrors its interior value in the
    ! ghost cell, with a minus (plus) sign for a Dirichlet (Neumann) condition; a boundary node of a
    ! staggered variable is not an unknown: it vanishes for a Dirichlet condition and equals its interior
    ! neighbour for a Neumann one, and the node on the upper boundary, which is part of the grid,
    ! gets a trivial equation
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_psi,nh_g,iq
    character(len=1), intent(in   ), dimension(0:1,3) :: cbc
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi,wzc
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(out  ), dimension(0:,:,:,:) :: coeffs
    real(rp), intent(inout), dimension(:,:,:)    :: rhs
    real(rp) :: rho,drho,mu,dmu,vol,rhop,cc,psic,psif,sgn
    logical  :: is_harm
    real(rp), dimension(0:1,3) :: cs
    integer , dimension(3) :: ijk,ijkn,e
    logical :: is_node
    integer :: i,j,k,idir,ib
    !
    call visc_coeffs(mu12,mu,dmu,is_harm)
    rho = rho12(2); drho = rho12(1)-rho12(2)
    e(:) = 0; if(iq > 0) e(iq) = 1
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ijk(:) = [i,j,k]
          !
          ! phase indicator at the location of the variable, and diffusion coefficient at the faces of
          ! its control volume
          !
          psic = 0.5*(psi(i,j,k)+psi(i+e(1),j+e(2),k+e(3)))
          do idir=1,3
            do ib=0,1
              ijkn(:) = ijk(:); ijkn(idir) = ijkn(idir) + 2*ib-1
              if(idir == iq) then
                psif = psi(i+ib*e(1),j+ib*e(2),k+ib*e(3))
                cs(ib,idir) = 2.*visc(is_harm,mu,dmu,psif)
              else
                psif = 0.5*(psi(ijkn(1),ijkn(2),ijkn(3))+psi(ijkn(1)+e(1),ijkn(2)+e(2),ijkn(3)+e(3)))
                cs(ib,idir) = visc(is_harm,mu,dmu,0.5*(psic+psif))
              end if
            end do
          end do
          !
          ! metric terms, times the height of the control volume
          !
          if(iq /= 3) then
            vol = 1./dzfi(k)
            cs(0,3) = cs(0,3)*dzci(k-1)
            cs(1,3) = cs(1,3)*dzci(k  )
          else
            vol = 1./dzci(k)
            cs(0,3) = cs(0,3)*dzfi(k  )
            cs(1,3) = cs(1,3)*dzfi(k+1)
          end if
          cs(:,1) = cs(:,1)*dli(1)**2*vol
          cs(:,2) = cs(:,2)*dli(2)**2*vol
          cs(:,:) = cs(:,:)*theta
          if(iq == 3) psic = wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1) ! volume average, for the density
          rhop = (rho + drho*psic)*vol
          cc   = rhop + sum(cs(:,:))
          !
          ! boundaries
          !
          is_node = .false.
          do idir=1,3
            do ib=0,1
              if(.not.(is_bound(ib,idir) .and. ijk(idir) == merge(n(idir),1,ib == 1))) cycle
              if(cbc(ib,idir) == 'P') cycle
              sgn = merge(-1._rp,1._rp,cbc(ib,idir) == 'D')
              if(idir /= iq) then
                cc = cc - sgn*cs(ib,idir)
                cs(ib,idir) = 0.
              else if(ib == 0) then
                if(sgn > 0.) cc = cc - cs(ib,idir)
                cs(ib,idir) = 0.
              else
                is_node = .true.
              end if
            end do
            if(idir == iq .and. is_bound(1,idir) .and. ijk(idir) == n(idir)-1 .and. cbc(1,idir) /= 'P') then
              if(cbc(1,idir) /= 'D') cc = cc - cs(1,idir)
              cs(1,idir) = 0.
            end if
          end do
          if(is_node) then
            coeffs(0  ,i,j,k) = rhop
            coeffs(1:6,i,j,k) = 0.
            rhs(i,j,k) = 0.
          else
            coeffs(0,i,j,k) = cc
            do idir=1,3
              coeffs(2*idir-1,i,j,k) = -cs(0,idir)
              coeffs(2*idir  ,i,j,k) = -cs(1,idir)
            end do
            rhs(i,j,k) = theta*rhs(i,j,k)*vol
          end if
        end do
      end do
    end do
  end subroutine helmholtz_coeffs
#endif
  !
  subroutine impdiff_wall_force(n,nh_psi,nh_g,is_forced,grid_vol_ratio_c,grid_vol_ratio_f,wzc,rho12,psi,psio, &
                                wall_force,dwall_force_old)
    !
    ! computes the wall shear force per unit volume carried by the terms that are integrated implicitly,
    ! i.e., minus their volume average, at the old time level (see `impdiff_z_mom_rhs`, `impdiff_cc_mom_rhs`);
    ! the zero-net-acceleration forcing leaves it out of the wall shear that it balances explicitly,
    ! and balances what the implicit integration applies together with the solve (see the main program);
    ! with constant coefficients (in all directions, with the constant-coefficient pressure solver) the terms
    ! are per unit mass, and `dwall_force_old` is the change of the force of the previous step due to the
    ! update of its density (zero otherwise)
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_psi,nh_g
    logical , intent(in ), dimension(3) :: is_forced
    real(rp), intent(in ), dimension(1-nh_g:) :: grid_vol_ratio_c,grid_vol_ratio_f,wzc
    real(rp), intent(in ), dimension(2) :: rho12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,psio
    real(rp), intent(out), dimension(3) :: wall_force,dwall_force_old
    real(rp) :: rho,drho,fx,fy,fz,gx,gy,gz,psixp,psiyp,psizp
    real(rp), dimension(6) :: force
    logical  :: is_cc
    integer  :: i,j,k,ierr
    !
#if defined(_CONSTANT_COEFFS_POISSON)
    is_cc = .not.is_impdiff_1d
#else
    is_cc = .false.
#endif
    rho = rho12(2); drho = rho12(1)-rho12(2)
    fx = 0._rp; fy = 0._rp; fz = 0._rp
    gx = 0._rp; gy = 0._rp; gz = 0._rp
    !$acc        data copy(fx,fy,fz,gx,gy,gz) async(1)
    !$omp target data map(tofrom:fx,fy,fz,gx,gy,gz)
    if(.not.is_cc) then
      !$acc parallel     loop collapse(3) default(present) reduction(+:fx,fy,fz) async(1)
      !$omp target teams loop collapse(3)                  reduction(+:fx,fy,fz)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            fx = fx - dud(i,j,k)*grid_vol_ratio_f(k)
            fy = fy - dvd(i,j,k)*grid_vol_ratio_f(k)
            fz = fz - dwd(i,j,k)*grid_vol_ratio_c(k)
          end do
        end do
      end do
    else
      !$acc parallel     loop collapse(3) default(present) private(psixp,psiyp,psizp) reduction(+:fx,fy,fz,gx,gy,gz) async(1)
      !$omp target teams loop collapse(3)                  private(psixp,psiyp,psizp) reduction(+:fx,fy,fz,gx,gy,gz)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psixp = 0.5*(psi(i,j,k)+psi(i+1,j,k))
            psiyp = 0.5*(psi(i,j,k)+psi(i,j+1,k))
            psizp = wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1)
            fx = fx - (rho+drho*psixp)*dud(i,j,k)*grid_vol_ratio_f(k)
            fy = fy - (rho+drho*psiyp)*dvd(i,j,k)*grid_vol_ratio_f(k)
            fz = fz - (rho+drho*psizp)*dwd(i,j,k)*grid_vol_ratio_c(k)
            gx = gx - drho*(0.5*(psio(i,j,k)+psio(i+1,j,k))-psixp)*dud(i,j,k)*grid_vol_ratio_f(k)
            gy = gy - drho*(0.5*(psio(i,j,k)+psio(i,j+1,k))-psiyp)*dvd(i,j,k)*grid_vol_ratio_f(k)
            gz = gz - drho*(wzc(k)*psio(i,j,k)+(1.-wzc(k))*psio(i,j,k+1)-psizp)*dwd(i,j,k)*grid_vol_ratio_c(k)
          end do
        end do
      end do
    end if
    !$omp end target data
    !$acc end        data
    !$acc wait(1)
    force(:) = [fx,fy,fz,gx,gy,gz]
    call MPI_ALLREDUCE(MPI_IN_PLACE,force,6,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    wall_force(:)      = force(1:3)
    dwall_force_old(:) = force(4:6)
    where(.not.is_forced(:))
      wall_force(:)      = 0._rp
      dwall_force_old(:) = 0._rp
    end where
  end subroutine impdiff_wall_force
  !
  subroutine impdiff_updt_forcing(n,nh_v,nh_psi,nh_g,wzc,rho12,psi,dtf,u,v,w)
    !
    ! adds a uniform force per unit volume, times the time step (`dtf`), to the velocity; it is used by
    ! the zero-net-acceleration forcing to balance the wall shear of the terms integrated implicitly, and
    ! the force of a variable surface tension (see `rk`)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in   ), dimension(1-nh_g:) :: wzc
    real(rp), intent(in   ), dimension(2) :: rho12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(3) :: dtf
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp) :: rho,drho,dtfx,dtfy,dtfz
    integer  :: i,j,k
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    dtfx = dtf(1); dtfy = dtf(2); dtfz = dtf(3)
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + dtfx/(rho + drho*0.5*(psi(i,j,k)+psi(i+1,j,k)))
          v(i,j,k) = v(i,j,k) + dtfy/(rho + drho*0.5*(psi(i,j,k)+psi(i,j+1,k)))
          w(i,j,k) = w(i,j,k) + dtfz/(rho + drho*(wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1)))
        end do
      end do
    end do
  end subroutine impdiff_updt_forcing
  !
  subroutine impdiff_wall_heat(n,nh_psi,nh_g,grid_vol_ratio,rhocp12,psi,psio,wall_heat,dwall_heat_old)
    !
    ! the counterpart of `impdiff_wall_force` for a scalar: computes the heat lost through the boundaries per unit
    ! volume by the terms that are integrated implicitly, i.e., minus their volume average, at the old time level
    ! (see `impdiff_z_scal_rhs`, `impdiff_cc_scal_rhs`, `impdiff_vc_scal_rhs`); the zero-net-heating forcing leaves
    ! it out of the heat loss that it balances explicitly, and balances what the implicit integration applies
    ! together with the solve (see the main program); with constant coefficients the terms are per unit rho*cp,
    ! and `dwall_heat_old` is the change of the heat loss of the previous step due to the update of its rho*cp
    ! (zero otherwise)
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_psi,nh_g
    real(rp), intent(in ), dimension(1-nh_g:) :: grid_vol_ratio
    real(rp), intent(in ), dimension(2) :: rhocp12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,psio
    real(rp), intent(out) :: wall_heat,dwall_heat_old
    real(rp) :: rhocp,drhocp,fs,gs
    real(rp), dimension(2) :: heat
    logical  :: is_cc
    integer  :: i,j,k,ierr
    !
#if defined(_CONSTANT_COEFFS_POISSON)
    is_cc = .not.is_impdiff_1d
#else
    is_cc = .false.
#endif
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    fs = 0._rp; gs = 0._rp
    !$acc        data copy(fs,gs) async(1)
    !$omp target data map(tofrom:fs,gs)
    if(.not.is_cc) then
      !$acc parallel     loop collapse(3) default(present) reduction(+:fs) async(1)
      !$omp target teams loop collapse(3)                  reduction(+:fs)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            fs = fs - dsd(i,j,k)*grid_vol_ratio(k)
          end do
        end do
      end do
    else
      !$acc parallel     loop collapse(3) default(present) reduction(+:fs,gs) async(1)
      !$omp target teams loop collapse(3)                  reduction(+:fs,gs)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            fs = fs - (rhocp+drhocp*psi(i,j,k))*dsd(i,j,k)*grid_vol_ratio(k)
            gs = gs - drhocp*(psio(i,j,k)-psi(i,j,k))*dsd(i,j,k)*grid_vol_ratio(k)
          end do
        end do
      end do
    end if
    !$omp end target data
    !$acc end        data
    !$acc wait(1)
    heat(:) = [fs,gs]
    call MPI_ALLREDUCE(MPI_IN_PLACE,heat,2,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    wall_heat      = heat(1)
    dwall_heat_old = heat(2)
  end subroutine impdiff_wall_heat
  !
  subroutine impdiff_updt_forcing_s(n,nh_s,nh_psi,rhocp12,psi,dtf,s)
    !
    ! the counterpart of `impdiff_updt_forcing` for a scalar: adds a uniform source per unit volume, times the time
    ! step (`dtf`), to the scalar; it is used by the zero-net-heating forcing to balance the heat lost through the
    ! boundaries by the terms integrated implicitly (see the main program)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_s,nh_psi
    real(rp), intent(in   ), dimension(2) :: rhocp12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ) :: dtf
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp) :: rhocp,drhocp
    integer  :: i,j,k
    !
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    !$acc parallel     loop collapse(3) default(present) async(1)
    !$omp target teams loop collapse(3)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + dtf/(rhocp + drhocp*psi(i,j,k))
        end do
      end do
    end do
  end subroutine impdiff_updt_forcing_s
  !
  subroutine impdiff_z_scal_rhs(n,nh_s,nh_psi,nh_g,dzci,dzfi,ka12,psi,s,dsdt)
    !
    ! computes the term of the scalar equation that is integrated implicitly, discretized exactly as in `scal_ad`,
    ! takes it out of the terms that are integrated explicitly (`dsdt`), and keeps it for the explicit half
    ! of the Crank-Nicolson scheme (see `impdiff_z_scal_solve`)
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_s,nh_psi,nh_g
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in   ), dimension(2) :: ka12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(inout), dimension(:,:,:) :: dsdt
    real(rp) :: ka,dka,kazp,kazm
    logical  :: is_harm
    integer  :: i,j,k
    !
    if(.not.allocated(dsd)) then
      allocate(dsd(n(1),n(2),n(3)))
      !$acc        enter data create(dsd) async(1)
      !$omp target enter data map(alloc:dsd)
    end if
    call visc_coeffs(ka12,ka,dka,is_harm)
    !$acc parallel     loop collapse(3) default(present) private(kazp,kazm) async(1)
    !$omp target teams loop collapse(3)                  private(kazp,kazm)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          kazp = visc(is_harm,ka,dka,0.5*(psi(i,j,k+1)+psi(i,j,k  )))
          kazm = visc(is_harm,ka,dka,0.5*(psi(i,j,k  )+psi(i,j,k-1)))
          dsd(i,j,k) = (kazp*(s(i,j,k+1)-s(i,j,k))*dzci(k)-kazm*(s(i,j,k)-s(i,j,k-1))*dzci(k-1))*dzfi(k)
          dsdt(i,j,k) = dsdt(i,j,k) - dsd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_z_scal_rhs
  !
  subroutine impdiff_z_scal_solve(n,nh_s,nh_psi,nh_g,cbc,bc,dzc,dzci,dzfi,theta,rhocp12,ka12,psi,s)
    !
    ! solves, for each line along z,
    !   rhocp s - theta*d/dz(ka ds/dz) = rhocp s* + theta*d/dz(ka ds/dz)^n,
    ! where s* is the input scalar and theta half of the time step (see `impdiff_z_mom_solve`);
    ! the boundary conditions along z are those that `boundp` imposes
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_s,nh_psi,nh_g
    character(len=1), intent(in   ), dimension(0:1,3) :: cbc
    real(rp), intent(in   ), dimension(0:1,3)      :: bc
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzci,dzfi
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rhocp12,ka12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    !
    if(.not.allocated(cw)) then
      allocate(cw(n(1),n(2),n(3)))
      !$acc        enter data create(cw) async(1)
      !$omp target enter data map(alloc:cw)
    end if
    call solve_lines_c(n,nh_s,nh_psi,nh_g,0,0,cbc(:,3),bc(:,3),dzc,dzci,dzfi,theta,rhocp12,ka12,psi,dsd,cw,s)
  end subroutine impdiff_z_scal_solve
  !
  subroutine solve_lines_c(n,nh,nh_psi,nh_g,di,dj,cbc,bc,dzc,dzci,dzfi,theta,rho12,mu12,psi,dpd,cw,p)
    !
    ! lines of a variable `p` that is cell-centered along z, and staggered along x (di = 1), along y (dj = 1),
    ! or cell-centered (di = dj = 0); `rho12` and `mu12` are the storage and diffusion coefficients of the two fluids,
    ! and `dpd` the diffusion term at the old time level;
    ! system: -a_k p_(k-1) + b_k p_k - c_k p_(k+1) = rho_k p*_k + theta*dpd_k, with b_k = rho_k + a_k + c_k
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh,nh_psi,nh_g,di,dj
    character(len=1), intent(in   ), dimension(0:1) :: cbc
    real(rp), intent(in   ), dimension(0:1)      :: bc
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzci,dzfi
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(:,:,:) :: dpd
    real(rp), intent(inout), dimension(:,:,:) :: cw
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: rho,drho,mu,dmu,a,b,c,d,rhop,z,cwm,pm,psim,psic,psip,muzm,muzp
    logical  :: is_harm
    real(rp) :: sgn0,sgn1,fac0,fac1
    integer  :: i,j,k,nz
    !
    nz = n(3)
    call visc_coeffs(mu12,mu,dmu,is_harm)
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    ! ghost values at the lower (0) and upper (1) boundaries: p(ghost) = sgn*p(interior) + fac (see `set_bc`)
    !
    call bc_ghost(cbc(0),2.*bc(0),-dzc(0 )*bc(0),sgn0,fac0)
    call bc_ghost(cbc(1),2.*bc(1), dzc(nz)*bc(1),sgn1,fac1)
    !$acc parallel     loop collapse(2) default(present) private(a,b,c,d,rhop,z,cwm,pm,psim,psic,psip,muzm,muzp) async(1)
    !$omp target teams loop collapse(2)                  private(a,b,c,d,rhop,z,cwm,pm,psim,psic,psip,muzm,muzp)
    do j=1,n(2)
      do i=1,n(1)
        !
        ! `psic` is the phase indicator at the location of p(i,j,k), `psim` and `psip` those below and above,
        ! and `muzm` and `muzp` the diffusion coefficients at the lower and upper faces
        !
        psim = 0.5*(psi(i,j,0)+psi(i+di,j+dj,0))
        psic = 0.5*(psi(i,j,1)+psi(i+di,j+dj,1))
        muzm = visc(is_harm,mu,dmu,0.5*(psim+psic))
        cwm = 0.; pm = 0.
        !$acc loop seq
        do k=1,nz
          psip = 0.5*(psi(i,j,k+1)+psi(i+di,j+dj,k+1))
          muzp = visc(is_harm,mu,dmu,0.5*(psic+psip))
          rhop = rho + drho*psic
          a = theta*dzfi(k)*dzci(k-1)*muzm
          c = theta*dzfi(k)*dzci(k  )*muzp
          b = rhop + a + c
          d = rhop*p(i,j,k) + theta*dpd(i,j,k)
          psic = psip
          muzm = muzp
          if(k == 1) then
            b = b - a*sgn0
            d = d + a*fac0
            a = 0.
          end if
          if(k == nz) then
            b = b - c*sgn1
            d = d + c*fac1
            c = 0.
          end if
          !
          ! forward elimination (`cwm` and `pm` are the values of the previous point, and `a` vanishes for the first)
          !
          z   = 1./(b - a*cwm)
          pm  = (d + a*pm)*z
          cwm = c*z
          p(i,j,k)  = pm
          cw(i,j,k) = cwm
        end do
        !
        ! backward substitution
        !
        !$acc loop seq
        do k=nz-1,1,-1
          pm = p(i,j,k) + cw(i,j,k)*pm
          p(i,j,k) = pm
        end do
      end do
    end do
  end subroutine solve_lines_c
  !
  subroutine solve_lines_f(n,nh,nh_psi,nh_g,cbc,bc,dzf,dzci,dzfi,wzc,theta,rho12,mu12,psi,dpd,cw,p)
    !
    ! lines of a variable `p` that is staggered along z (the vertical velocity), with twice the viscosity
    ! (see `solve_lines_c` for the other arguments);
    ! the node `nz` lies on the upper boundary, and is prescribed by the boundary condition
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh,nh_psi,nh_g
    character(len=1), intent(in   ), dimension(0:1) :: cbc
    real(rp), intent(in   ), dimension(0:1)      :: bc
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzf,dzci,dzfi,wzc
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(:,:,:) :: dpd
    real(rp), intent(inout), dimension(:,:,:) :: cw
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: rho,drho,mu,dmu,a,b,c,d,rhop,z,cwm,pm,psim,psip,muzm,muzp
    logical  :: is_harm
    real(rp) :: sgn0,sgn1,fac0,fac1
    integer  :: i,j,k,nz
    !
    nz = n(3)
    call visc_coeffs(mu12,mu,dmu,is_harm)
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    ! lower boundary node:  p(0 ) = sgn0*p(1   ) + fac0
    ! upper boundary node:  p(nz) = sgn1*p(nz-1) + fac1 (see `set_bc`)
    !
    call bc_ghost(cbc(0),bc(0),-dzf(0 )*bc(0),sgn0,fac0); if(cbc(0) == 'D') sgn0 = 0.
    call bc_ghost(cbc(1),bc(1), dzf(nz)*bc(1),sgn1,fac1); if(cbc(1) == 'D') sgn1 = 0.
    !$acc parallel     loop collapse(2) default(present) private(a,b,c,d,rhop,z,cwm,pm,psim,psip,muzm,muzp) async(1)
    !$omp target teams loop collapse(2)                  private(a,b,c,d,rhop,z,cwm,pm,psim,psip,muzm,muzp)
    do j=1,n(2)
      do i=1,n(1)
        !
        ! `psim` and `psip` are the phase indicator at the cells below and above p(i,j,k),
        ! and `muzm` and `muzp` the viscosities there
        !
        psim = psi(i,j,1)
        muzm = visc(is_harm,mu,dmu,psim)
        cwm = 0.; pm = 0.
        !$acc loop seq
        do k=1,nz
          if(k < nz) then
            psip = psi(i,j,k+1)
            muzp = visc(is_harm,mu,dmu,psip)
            rhop = rho + drho*(wzc(k)*psim+(1.-wzc(k))*psip)
            a = theta*dzci(k)*2.*dzfi(k  )*muzm
            c = theta*dzci(k)*2.*dzfi(k+1)*muzp
            b = rhop + a + c
            d = rhop*p(i,j,k) + theta*dpd(i,j,k)
            psim = psip
            muzm = muzp
            if(k == 1) then
              b = b - a*sgn0
              d = d + a*fac0
              a = 0.
            end if
          else
            a = sgn1
            b = 1.
            c = 0.
            d = fac1
          end if
          z   = 1./(b - a*cwm)
          pm  = (d + a*pm)*z
          cwm = c*z
          p(i,j,k)  = pm
          cw(i,j,k) = cwm
        end do
        !$acc loop seq
        do k=nz-1,1,-1
          pm = p(i,j,k) + cw(i,j,k)*pm
          p(i,j,k) = pm
        end do
      end do
    end do
  end subroutine solve_lines_f
  !
  subroutine bc_ghost(ctype,fac_d,fac_n,sgn,fac)
    !
    ! coefficients of the ghost (or boundary-node) value, p(ghost) = sgn*p(interior) + fac,
    ! of a Dirichlet (fac_d) or Neumann (fac_n) boundary condition
    !
    implicit none
    character(len=1), intent(in ) :: ctype
    real(rp)        , intent(in ) :: fac_d,fac_n
    real(rp)        , intent(out) :: sgn,fac
    select case(ctype)
    case('D')
      sgn = -1.; fac = fac_d
    case default
      sgn =  1.; fac = fac_n
    end select
  end subroutine bc_ghost
  !
  subroutine visc_coeffs(mu12,mu,dmu,is_harm)
    !
    ! coefficients of the interpolation of the viscosity between the two fluids (see `visc`),
    ! and whether it is harmonic
    !
    implicit none
    real(rp), intent(in ), dimension(2) :: mu12
    real(rp), intent(out) :: mu,dmu
    logical , intent(out) :: is_harm
    is_harm = is_harmonic_interpolation
    if(is_harm) then
      mu = 1._rp/mu12(2); dmu = 1._rp/mu12(1)-1._rp/mu12(2)
    else
      mu = mu12(2); dmu = mu12(1)-mu12(2)
    end if
  end subroutine visc_coeffs
  !
  pure real(rp) function visc(is_harm,mu,dmu,psi)
    !
    ! viscosity for a value `psi` of the phase indicator: arithmetic or harmonic interpolation
    !
    !$acc routine seq
    !$omp declare target
    implicit none
    logical , intent(in) :: is_harm
    real(rp), intent(in) :: mu,dmu,psi
    visc = mu+dmu*psi
    if(is_harm) visc = 1._rp/visc
  end function visc
end module mod_impdiff
