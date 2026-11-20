! -
!
! SPDX-FileCopyrightText: Copyright (c) 2024 The CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_impdiff
  !
  ! implicit integration of the diffusion along z (`is_impdiff_1d = T`)
  !
  ! the part of the viscous term that is integrated implicitly is the diffusion of each velocity component
  ! along z due to its own gradient,
  !   d/dz(mu du/dz), d/dz(mu dv/dz) and d/dz(2 mu dw/dz),
  ! with the local viscosity and density, i.e., one tridiagonal system with variable coefficients per line
  ! along z; the remaining viscous terms are integrated explicitly; likewise for the scalar, d/dz(ka ds/dz)
  !
  use mpi
  use mod_types
  implicit none
  private
  public impdiff_z_mom_rhs,impdiff_z_mom_solve,impdiff_z_scal_rhs,impdiff_z_scal_solve
  public impdiff_cc_mom_rhs,impdiff_cc_mom_updt,impdiff_cc_scal_rhs,impdiff_cc_scal_updt
  public impdiff_wall_force,impdiff_updt_forcing
  real(rp), allocatable, dimension(:,:,:), save :: cw            ! work array of the line solver
  real(rp), allocatable, dimension(:,:,:), save :: dud,dvd,dwd,dsd ! implicit terms at the old time level
  contains
  subroutine impdiff_z_mom_rhs(n,nh_v,nh_psi,nh_g,dzci,dzfi,mu12,psi,u,v,w,dudt,dvdt,dwdt)
    !
    ! computes the terms that are integrated implicitly, discretized exactly as in `mom_xyz_ad`,
    ! takes them out of the terms that are integrated explicitly (`dudt`, `dvdt`, `dwdt`),
    ! and keeps them for the explicit half of the Crank-Nicolson scheme (see `impdiff_z_mom_solve`)
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: mu12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(inout), dimension(:,:,:) :: dudt,dvdt,dwdt
    real(rp) :: mu,dmu,muzp,muzm
    integer  :: i,j,k
    !
    if(.not.allocated(dud)) then
      allocate(dud(n(1),n(2),n(3)),dvd(n(1),n(2),n(3)),dwd(n(1),n(2),n(3)))
      !$acc enter data create(dud,dvd,dwd) async(1)
    end if
    call visc_coeffs(mu12,mu,dmu)
    !$acc parallel loop collapse(3) default(present) private(muzp,muzm) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)  private(muzp,muzm)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          muzp = visc(mu,dmu,0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)))
          muzm = visc(mu,dmu,0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)))
          dud(i,j,k) = dzfi(k)*((u(i,j,k+1)-u(i,j,k))*dzci(k)*muzp-(u(i,j,k)-u(i,j,k-1))*dzci(k-1)*muzm)
          muzp = visc(mu,dmu,0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1)))
          muzm = visc(mu,dmu,0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1)))
          dvd(i,j,k) = dzfi(k)*((v(i,j,k+1)-v(i,j,k))*dzci(k)*muzp-(v(i,j,k)-v(i,j,k-1))*dzci(k-1)*muzm)
          muzp = visc(mu,dmu,psi(i,j,k+1))
          muzm = visc(mu,dmu,psi(i,j,k  ))
          dwd(i,j,k) = dzci(k)*(2.*(w(i,j,k+1)-w(i,j,k))*dzfi(k+1)*muzp-2.*(w(i,j,k)-w(i,j,k-1))*dzfi(k)*muzm)
          dudt(i,j,k) = dudt(i,j,k) - dud(i,j,k)
          dvdt(i,j,k) = dvdt(i,j,k) - dvd(i,j,k)
          dwdt(i,j,k) = dwdt(i,j,k) - dwd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_z_mom_rhs
  !
  subroutine impdiff_z_mom_solve(n,nh_v,nh_psi,nh_g,cbc,bc,dzc,dzf,dzci,dzfi,theta,rho12,mu12,psi,u,v,w)
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
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    real(rp), intent(in   ), dimension(0:1,3,3)      :: bc
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzf,dzci,dzfi
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    !
    if(.not.allocated(cw)) then
      allocate(cw(n(1),n(2),n(3)))
      !$acc enter data create(cw) async(1)
    end if
    call solve_lines_c(n,nh_v,nh_psi,nh_g,1,0,cbc(:,3,1),bc(:,3,1),dzc,dzci,dzfi,theta,rho12,mu12,psi,dud,cw,u)
    call solve_lines_c(n,nh_v,nh_psi,nh_g,0,1,cbc(:,3,2),bc(:,3,2),dzc,dzci,dzfi,theta,rho12,mu12,psi,dvd,cw,v)
    call solve_lines_f(n,nh_v,nh_psi,nh_g    ,cbc(:,3,3),bc(:,3,3),dzf,dzci,dzfi,theta,rho12,mu12,psi,dwd,cw,w)
  end subroutine impdiff_z_mom_solve
  !
  subroutine impdiff_cc_mom_rhs(n,nh_v,nh_psi,nh_g,dli,dzci,dzfi,visc0,rho12,psi,psio,u,v,w,dudt,dvdt,dwdt,dudto,dvdto,dwdto)
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
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzci,dzfi
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
      !$acc enter data create(dud,dvd,dwd) async(1)
    end if
    rho = rho12(2); drho = rho12(1)-rho12(2)
    dxi2 = dli(1)**2; dyi2 = dli(2)**2
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
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
          dwdt(i,j,k) = dwdt(i,j,k) - (rho + drho*0.5*(psi(i,j,k)+psi(i,j,k+1)))*dwd(i,j,k)
          dudto(i,j,k) = dudto(i,j,k) + drho*0.5*(psio(i,j,k)+psio(i+1,j,k)-psi(i,j,k)-psi(i+1,j,k))*dud(i,j,k)
          dvdto(i,j,k) = dvdto(i,j,k) + drho*0.5*(psio(i,j,k)+psio(i,j+1,k)-psi(i,j,k)-psi(i,j+1,k))*dvd(i,j,k)
          dwdto(i,j,k) = dwdto(i,j,k) + drho*0.5*(psio(i,j,k)+psio(i,j,k+1)-psi(i,j,k)-psi(i,j,k+1))*dwd(i,j,k)
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
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
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
      !$acc enter data create(dsd) async(1)
    end if
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    dxi2 = dli(1)**2; dyi2 = dli(2)**2
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
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
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s(i,j,k) = s(i,j,k) + theta*dsd(i,j,k)
        end do
      end do
    end do
  end subroutine impdiff_cc_scal_updt
  !
  subroutine impdiff_wall_force(n,nh_psi,nh_g,is_forced,is_cc,grid_vol_ratio_c,grid_vol_ratio_f,rho12,psi,psio, &
                                wall_force,dwall_force_old)
    !
    ! computes the wall shear force per unit volume carried by the terms that are integrated implicitly,
    ! i.e., minus their volume average, at the old time level (see `impdiff_z_mom_rhs`, `impdiff_cc_mom_rhs`);
    ! the zero-net-acceleration forcing leaves it out of the wall shear that it balances explicitly,
    ! and balances what the implicit integration applies together with the solve (see the main program);
    ! with constant coefficients (`is_cc`) the terms are per unit mass, and `dwall_force_old` is the change of
    ! the force of the previous step due to the update of its density (zero otherwise)
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_psi,nh_g
    logical , intent(in ), dimension(3) :: is_forced
    logical , intent(in )               :: is_cc
    real(rp), intent(in ), dimension(1-nh_g:) :: grid_vol_ratio_c,grid_vol_ratio_f
    real(rp), intent(in ), dimension(2) :: rho12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi,psio
    real(rp), intent(out), dimension(3) :: wall_force,dwall_force_old
    real(rp) :: rho,drho,fx,fy,fz,gx,gy,gz,psixp,psiyp,psizp
    real(rp), dimension(6) :: force
    integer  :: i,j,k,ierr
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    fx = 0._rp; fy = 0._rp; fz = 0._rp
    gx = 0._rp; gy = 0._rp; gz = 0._rp
    !$acc data copy(fx,fy,fz,gx,gy,gz) async(1)
    if(.not.is_cc) then
      !$acc parallel loop collapse(3) default(present) reduction(+:fx,fy,fz) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)  reduction(+:fx,fy,fz)
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
      !$acc parallel loop collapse(3) default(present) private(psixp,psiyp,psizp) reduction(+:fx,fy,fz,gx,gy,gz) async(1)
      !$OMP parallel do   collapse(3) DEFAULT(shared)  private(psixp,psiyp,psizp) reduction(+:fx,fy,fz,gx,gy,gz)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psixp = 0.5*(psi(i,j,k)+psi(i+1,j,k))
            psiyp = 0.5*(psi(i,j,k)+psi(i,j+1,k))
            psizp = 0.5*(psi(i,j,k)+psi(i,j,k+1))
            fx = fx - (rho+drho*psixp)*dud(i,j,k)*grid_vol_ratio_f(k)
            fy = fy - (rho+drho*psiyp)*dvd(i,j,k)*grid_vol_ratio_f(k)
            fz = fz - (rho+drho*psizp)*dwd(i,j,k)*grid_vol_ratio_c(k)
            gx = gx - drho*(0.5*(psio(i,j,k)+psio(i+1,j,k))-psixp)*dud(i,j,k)*grid_vol_ratio_f(k)
            gy = gy - drho*(0.5*(psio(i,j,k)+psio(i,j+1,k))-psiyp)*dvd(i,j,k)*grid_vol_ratio_f(k)
            gz = gz - drho*(0.5*(psio(i,j,k)+psio(i,j,k+1))-psizp)*dwd(i,j,k)*grid_vol_ratio_c(k)
          end do
        end do
      end do
    end if
    !$acc end data
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
  subroutine impdiff_updt_forcing(n,nh_v,nh_psi,rho12,psi,dtf,u,v,w)
    !
    ! adds a uniform force per unit volume, times the time step (`dtf`), to the velocity; it is used by
    ! the zero-net-acceleration forcing to balance the wall shear of the terms integrated implicitly
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh_v,nh_psi
    real(rp), intent(in   ), dimension(2) :: rho12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(3) :: dtf
    real(rp), intent(inout), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp) :: rho,drho,dtfx,dtfy,dtfz
    integer  :: i,j,k
    !
    rho = rho12(2); drho = rho12(1)-rho12(2)
    dtfx = dtf(1); dtfy = dtf(2); dtfz = dtf(3)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          u(i,j,k) = u(i,j,k) + dtfx/(rho + drho*0.5*(psi(i,j,k)+psi(i+1,j,k)))
          v(i,j,k) = v(i,j,k) + dtfy/(rho + drho*0.5*(psi(i,j,k)+psi(i,j+1,k)))
          w(i,j,k) = w(i,j,k) + dtfz/(rho + drho*0.5*(psi(i,j,k)+psi(i,j,k+1)))
        end do
      end do
    end do
  end subroutine impdiff_updt_forcing
  !
  subroutine impdiff_z_scal_rhs(n,nh_s,nh_psi,nh_g,dzci,dzfi,ka12,psi,s,dsdt)
    !
    ! computes the term of the scalar equation that is integrated implicitly, discretized exactly as in `scal_ad`,
    ! takes it out of the terms that are integrated explicitly (`dsdt`), and keeps it for the explicit half
    ! of the Crank-Nicolson scheme (see `impdiff_z_scal_solve`)
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_s,nh_psi,nh_g
    real(rp), intent(in ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: ka12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(inout), dimension(:,:,:) :: dsdt
    real(rp) :: ka,dka,kazp,kazm
    integer  :: i,j,k
    !
    if(.not.allocated(dsd)) then
      allocate(dsd(n(1),n(2),n(3)))
      !$acc enter data create(dsd) async(1)
    end if
    call visc_coeffs(ka12,ka,dka)
    !$acc parallel loop collapse(3) default(present) private(kazp,kazm) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)  private(kazp,kazm)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          kazp = visc(ka,dka,0.5*(psi(i,j,k+1)+psi(i,j,k  )))
          kazm = visc(ka,dka,0.5*(psi(i,j,k  )+psi(i,j,k-1)))
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
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp), intent(in   ), dimension(0:1,3)      :: bc
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzci,dzfi
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rhocp12,ka12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    !
    if(.not.allocated(cw)) then
      allocate(cw(n(1),n(2),n(3)))
      !$acc enter data create(cw) async(1)
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
    character(len=1), intent(in), dimension(0:1) :: cbc
    real(rp), intent(in   ), dimension(0:1)      :: bc
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzc,dzci,dzfi
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(:,:,:) :: dpd
    real(rp), intent(inout), dimension(:,:,:) :: cw
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: rho,drho,mu,dmu,a,b,c,d,rhop,z,cwm,pm,psim,psic,psip,muzm,muzp
    real(rp) :: sgn0,sgn1,fac0,fac1
    integer  :: i,j,k,nz
    !
    nz = n(3)
    call visc_coeffs(mu12,mu,dmu)
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    ! ghost values at the lower (0) and upper (1) boundaries: p(ghost) = sgn*p(interior) + fac (see `set_bc`)
    !
    call bc_ghost(cbc(0),2.*bc(0),-dzc(0 )*bc(0),sgn0,fac0)
    call bc_ghost(cbc(1),2.*bc(1), dzc(nz)*bc(1),sgn1,fac1)
    !$acc parallel loop collapse(2) default(present) private(a,b,c,d,rhop,z,cwm,pm,psim,psic,psip,muzm,muzp) async(1)
    !$OMP parallel do   collapse(2) DEFAULT(shared)  private(a,b,c,d,rhop,z,cwm,pm,psim,psic,psip,muzm,muzp)
    do j=1,n(2)
      do i=1,n(1)
        !
        ! `psic` is the phase indicator at the location of p(i,j,k), `psim` and `psip` those below and above,
        ! and `muzm` and `muzp` the diffusion coefficients at the lower and upper faces
        !
        psim = 0.5*(psi(i,j,0)+psi(i+di,j+dj,0))
        psic = 0.5*(psi(i,j,1)+psi(i+di,j+dj,1))
        muzm = visc(mu,dmu,0.5*(psim+psic))
        cwm = 0.; pm = 0.
        !$acc loop seq
        do k=1,nz
          psip = 0.5*(psi(i,j,k+1)+psi(i+di,j+dj,k+1))
          muzp = visc(mu,dmu,0.5*(psic+psip))
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
  subroutine solve_lines_f(n,nh,nh_psi,nh_g,cbc,bc,dzf,dzci,dzfi,theta,rho12,mu12,psi,dpd,cw,p)
    !
    ! lines of a variable `p` that is staggered along z (the vertical velocity), with twice the viscosity
    ! (see `solve_lines_c` for the other arguments);
    ! the node `nz` lies on the upper boundary, and is prescribed by the boundary condition
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   )               :: nh,nh_psi,nh_g
    character(len=1), intent(in), dimension(0:1) :: cbc
    real(rp), intent(in   ), dimension(0:1)      :: bc
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzf,dzci,dzfi
    real(rp), intent(in   )               :: theta
    real(rp), intent(in   ), dimension(2) :: rho12,mu12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(:,:,:) :: dpd
    real(rp), intent(inout), dimension(:,:,:) :: cw
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: rho,drho,mu,dmu,a,b,c,d,rhop,z,cwm,pm,psim,psip,muzm,muzp
    real(rp) :: sgn0,sgn1,fac0,fac1
    integer  :: i,j,k,nz
    !
    nz = n(3)
    call visc_coeffs(mu12,mu,dmu)
    rho = rho12(2); drho = rho12(1)-rho12(2)
    !
    ! lower boundary node:  p(0 ) = sgn0*p(1   ) + fac0
    ! upper boundary node:  p(nz) = sgn1*p(nz-1) + fac1 (see `set_bc`)
    !
    call bc_ghost(cbc(0),bc(0),-dzf(0 )*bc(0),sgn0,fac0); if(cbc(0) == 'D') sgn0 = 0.
    call bc_ghost(cbc(1),bc(1), dzf(nz)*bc(1),sgn1,fac1); if(cbc(1) == 'D') sgn1 = 0.
    !$acc parallel loop collapse(2) default(present) private(a,b,c,d,rhop,z,cwm,pm,psim,psip,muzm,muzp) async(1)
    !$OMP parallel do   collapse(2) DEFAULT(shared)  private(a,b,c,d,rhop,z,cwm,pm,psim,psip,muzm,muzp)
    do j=1,n(2)
      do i=1,n(1)
        !
        ! `psim` and `psip` are the phase indicator at the cells below and above p(i,j,k),
        ! and `muzm` and `muzp` the viscosities there
        !
        psim = psi(i,j,1)
        muzm = visc(mu,dmu,psim)
        cwm = 0.; pm = 0.
        !$acc loop seq
        do k=1,nz
          if(k < nz) then
            psip = psi(i,j,k+1)
            muzp = visc(mu,dmu,psip)
            rhop = rho + drho*0.5*(psim+psip)
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
  subroutine visc_coeffs(mu12,mu,dmu)
    !
    ! coefficients of the interpolation of the viscosity between the two fluids (see `visc`)
    !
    implicit none
    real(rp), intent(in ), dimension(2) :: mu12
    real(rp), intent(out) :: mu,dmu
#if defined(_VISC_HARMONIC_INTERPOLATION)
    mu = 1._rp/mu12(2); dmu = 1._rp/mu12(1)-1._rp/mu12(2)
#else
    mu = mu12(2); dmu = mu12(1)-mu12(2)
#endif
  end subroutine visc_coeffs
  !
  pure real(rp) function visc(mu,dmu,psi)
    !
    ! viscosity for a value `psi` of the phase indicator: arithmetic or harmonic interpolation
    !
    !$acc routine seq
    implicit none
    real(rp), intent(in) :: mu,dmu,psi
#if defined(_VISC_HARMONIC_INTERPOLATION)
    visc = 1._rp/(mu+dmu*psi)
#else
    visc = mu+dmu*psi
#endif
  end function visc
end module mod_impdiff
