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
  use mod_types
  implicit none
  private
  public impdiff_z_mom_rhs,impdiff_z_mom_solve,impdiff_z_scal_rhs,impdiff_z_scal_solve,impdiff_updt_forcing
  real(rp), allocatable, dimension(:,:,:), save :: cw ! work array of the line solver
  contains
  subroutine impdiff_z_mom_rhs(n,nh_v,nh_psi,nh_g,dzci,dzfi,mu12,psi,u,v,w,dudtd,dvdtd,dwdtd)
    !
    ! computes the terms that are integrated implicitly, discretized exactly as in `mom_xyz_ad`
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_v,nh_psi,nh_g
    real(rp), intent(in ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: mu12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(out), dimension(:,:,:) :: dudtd,dvdtd,dwdtd
    real(rp) :: mu,dmu,muzp,muzm
    integer  :: i,j,k
    !
    call visc_coeffs(mu12,mu,dmu)
    !$acc parallel loop collapse(3) default(present) private(muzp,muzm) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)  private(muzp,muzm)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          muzp = visc(mu,dmu,0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i,j,k+1)+psi(i+1,j,k+1)))
          muzm = visc(mu,dmu,0.25*(psi(i,j,k)+psi(i+1,j,k)+psi(i,j,k-1)+psi(i+1,j,k-1)))
          dudtd(i,j,k) = dzfi(k)*((u(i,j,k+1)-u(i,j,k))*dzci(k)*muzp-(u(i,j,k)-u(i,j,k-1))*dzci(k-1)*muzm)
          muzp = visc(mu,dmu,0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k+1)+psi(i,j,k+1)))
          muzm = visc(mu,dmu,0.25*(psi(i,j,k)+psi(i,j+1,k)+psi(i,j+1,k-1)+psi(i,j,k-1)))
          dvdtd(i,j,k) = dzfi(k)*((v(i,j,k+1)-v(i,j,k))*dzci(k)*muzp-(v(i,j,k)-v(i,j,k-1))*dzci(k-1)*muzm)
          muzp = visc(mu,dmu,psi(i,j,k+1))
          muzm = visc(mu,dmu,psi(i,j,k  ))
          dwdtd(i,j,k) = dzci(k)*(2.*(w(i,j,k+1)-w(i,j,k))*dzfi(k+1)*muzp-2.*(w(i,j,k)-w(i,j,k-1))*dzfi(k)*muzm)
        end do
      end do
    end do
  end subroutine impdiff_z_mom_rhs
  !
  subroutine impdiff_z_mom_solve(n,nh_v,nh_psi,nh_g,cbc,bc,dzc,dzf,dzci,dzfi,theta,rho12,mu12,psi,u,v,w)
    !
    ! solves, for each velocity component u_q and each line along z,
    !   rho_q u_q - theta*d/dz(c_q mu du_q/dz) = rho_q u_q*,
    ! where u_q* is the input velocity, rho_q the density at the location of u_q, c_q = 2 for w and 1
    ! otherwise, and theta the implicit fraction of the time step (dt/2 for Crank-Nicolson);
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
    call solve_lines_c(n,nh_v,nh_psi,nh_g,1,0,cbc(:,3,1),bc(:,3,1),dzc,dzci,dzfi,theta,rho12,mu12,psi,cw,u)
    call solve_lines_c(n,nh_v,nh_psi,nh_g,0,1,cbc(:,3,2),bc(:,3,2),dzc,dzci,dzfi,theta,rho12,mu12,psi,cw,v)
    call solve_lines_f(n,nh_v,nh_psi,nh_g    ,cbc(:,3,3),bc(:,3,3),dzf,dzci,dzfi,theta,rho12,mu12,psi,cw,w)
  end subroutine impdiff_z_mom_solve
  !
  subroutine impdiff_updt_forcing(n,nh_v,nh_psi,rho12,psi,dtf,u,v,w)
    !
    ! adds a uniform force per unit volume, times the time step (`dtf`), to the velocity;
    ! it corrects the zero-net-acceleration forcing, since the wall shear of the implicit half of the
    ! diffusion is only known after the implicit solve
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
  subroutine impdiff_z_scal_rhs(n,nh_s,nh_psi,nh_g,dzci,dzfi,ka12,psi,s,dsdtd)
    !
    ! computes the term of the scalar equation that is integrated implicitly, discretized exactly as in `scal_ad`
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in )               :: nh_s,nh_psi,nh_g
    real(rp), intent(in ), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: ka12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(out), dimension(:,:,:) :: dsdtd
    real(rp) :: ka,dka,kazp,kazm
    integer  :: i,j,k
    !
    call visc_coeffs(ka12,ka,dka)
    !$acc parallel loop collapse(3) default(present) private(kazp,kazm) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)  private(kazp,kazm)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          kazp = visc(ka,dka,0.5*(psi(i,j,k+1)+psi(i,j,k  )))
          kazm = visc(ka,dka,0.5*(psi(i,j,k  )+psi(i,j,k-1)))
          dsdtd(i,j,k) = (kazp*(s(i,j,k+1)-s(i,j,k))*dzci(k)-kazm*(s(i,j,k)-s(i,j,k-1))*dzci(k-1))*dzfi(k)
        end do
      end do
    end do
  end subroutine impdiff_z_scal_rhs
  !
  subroutine impdiff_z_scal_solve(n,nh_s,nh_psi,nh_g,cbc,bc,dzc,dzci,dzfi,theta,rhocp12,ka12,psi,s)
    !
    ! solves, for each line along z,
    !   rhocp s - theta*d/dz(ka ds/dz) = rhocp s*,
    ! where s* is the input scalar and theta the implicit fraction of the time step (see `impdiff_z_mom_solve`);
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
    call solve_lines_c(n,nh_s,nh_psi,nh_g,0,0,cbc(:,3),bc(:,3),dzc,dzci,dzfi,theta,rhocp12,ka12,psi,cw,s)
  end subroutine impdiff_z_scal_solve
  !
  subroutine solve_lines_c(n,nh,nh_psi,nh_g,di,dj,cbc,bc,dzc,dzci,dzfi,theta,rho12,mu12,psi,cw,p)
    !
    ! lines of a variable `p` that is cell-centered along z, and staggered along x (di = 1), along y (dj = 1),
    ! or cell-centered (di = dj = 0); `rho12` and `mu12` are the storage and diffusion coefficients of the two fluids;
    ! system: -a_k p_(k-1) + b_k p_k - c_k p_(k+1) = rho_k p*_k, with b_k = rho_k + a_k + c_k
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
    real(rp), intent(inout), dimension(:,:,:) :: cw
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: rho,drho,mu,dmu,a,b,c,d,rhop,z
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
    !$acc parallel loop collapse(2) default(present) private(a,b,c,d,rhop,z) async(1)
    !$OMP parallel do   collapse(2) DEFAULT(shared)  private(a,b,c,d,rhop,z)
    do j=1,n(2)
      do i=1,n(1)
        !$acc loop seq
        do k=1,nz
          rhop = rho + drho*0.5*(psi(i,j,k)+psi(i+di,j+dj,k))
          a = theta*dzfi(k)*dzci(k-1)*visc(mu,dmu,0.25*(psi(i,j,k)+psi(i+di,j+dj,k)+psi(i,j,k-1)+psi(i+di,j+dj,k-1)))
          c = theta*dzfi(k)*dzci(k  )*visc(mu,dmu,0.25*(psi(i,j,k)+psi(i+di,j+dj,k)+psi(i,j,k+1)+psi(i+di,j+dj,k+1)))
          b = rhop + a + c
          d = rhop*p(i,j,k)
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
          ! forward elimination
          !
          if(k == 1) then
            z = 1./b
            p(i,j,k) = d*z
          else
            z = 1./(b - a*cw(i,j,k-1))
            p(i,j,k) = (d + a*p(i,j,k-1))*z
          end if
          cw(i,j,k) = c*z
        end do
        !
        ! backward substitution
        !
        !$acc loop seq
        do k=nz-1,1,-1
          p(i,j,k) = p(i,j,k) + cw(i,j,k)*p(i,j,k+1)
        end do
      end do
    end do
  end subroutine solve_lines_c
  !
  subroutine solve_lines_f(n,nh,nh_psi,nh_g,cbc,bc,dzf,dzci,dzfi,theta,rho12,mu12,psi,cw,p)
    !
    ! lines of a variable `p` that is staggered along z (the vertical velocity), with twice the viscosity;
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
    real(rp), intent(inout), dimension(:,:,:) :: cw
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: rho,drho,mu,dmu,a,b,c,d,rhop,z
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
    !$acc parallel loop collapse(2) default(present) private(a,b,c,d,rhop,z) async(1)
    !$OMP parallel do   collapse(2) DEFAULT(shared)  private(a,b,c,d,rhop,z)
    do j=1,n(2)
      do i=1,n(1)
        !$acc loop seq
        do k=1,nz
          if(k < nz) then
            rhop = rho + drho*0.5*(psi(i,j,k)+psi(i,j,k+1))
            a = theta*dzci(k)*2.*dzfi(k  )*visc(mu,dmu,psi(i,j,k  ))
            c = theta*dzci(k)*2.*dzfi(k+1)*visc(mu,dmu,psi(i,j,k+1))
            b = rhop + a + c
            d = rhop*p(i,j,k)
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
          if(k == 1) then
            z = 1./b
            p(i,j,k) = d*z
          else
            z = 1./(b - a*cw(i,j,k-1))
            p(i,j,k) = (d + a*p(i,j,k-1))*z
          end if
          cw(i,j,k) = c*z
        end do
        !$acc loop seq
        do k=nz-1,1,-1
          p(i,j,k) = p(i,j,k) + cw(i,j,k)*p(i,j,k+1)
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
