! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_scal
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mod_param, only: is_harmonic_interpolation,is_ifc_split
  use mod_solve_helmholtz, only: rhs_bound
  use mod_types
  implicit none
  private
  public scal_ad,scal_a_cons,cmpt_scalflux,scalar
  !
  ! scalar derived type, as in CaNS: the field, its right-hand side of the previous time step, the previous
  ! and transported fields of the direction-split advection, and the data of the constant-coefficient
  ! implicit diffusion; the input parameters of each scalar are those of `mod_param`
  !
  type scalar
    real(rp), allocatable, dimension(:,:,:) :: val
    real(rp), allocatable, dimension(:,:,:) :: dsdtrko
    real(rp), allocatable, dimension(:,:,:) :: so,sadv
#if !defined(_OPENACC) || defined(_USE_HIP)
    type(C_PTR), dimension(2,2) :: arrplan
#else
    integer    , dimension(2,2) :: arrplan
#endif
    real(rp), allocatable, dimension(:,:) :: lambdaxy
    real(rp), allocatable, dimension(:) :: a,b,c
    real(rp) :: normfft
    type(rhs_bound) :: rhsb
  end type scalar
  !
  contains
  subroutine scal_ad(nx,ny,nz,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dxi,dyi,dzci,dzfi,ka12,rhocp12,psiflx_x,psiflx_y,psiflx_z,psi,u,v,w, &
                     s,dsdt)
    !
    implicit none
    integer , intent(in) :: nx,ny,nz
    integer , intent(in) :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in) :: dxi,dyi
    real(rp), intent(in), dimension(1-nh_g:) :: dzci,dzfi
    real(rp), intent(in), dimension(2) :: ka12,rhocp12
    real(rp), intent(in), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:), optional :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: u,v,w
    real(rp), intent(in), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), dimension(:,:,:), intent(out) :: dsdt
    integer :: i,j,k
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: kaxp,kaxm,kayp,kaym,kazp,kazm,rhocp_c
    real(rp) :: psixp,psixm,psiyp,psiym,psizp,psizm
    real(rp) :: ka,dka,rhocp,drhocp
    logical  :: is_adv,is_harm
    !
    is_adv = .not.is_ifc_split ! (otherwise, advection is computed in `scal_a_cons` or along with the phase indicator)
    is_harm = is_harmonic_interpolation
    if(is_harm) then
      ! N.b.: care should be taken if one of the conductivities is zero.
      ka = 1._rp/ka12(2); dka = 1._rp/ka12(1)-1._rp/ka12(2)
    else
      ka    = ka12(2)   ; dka    = ka12(1)-ka12(2)
    end if
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          if(is_adv) then
            usip = 0.5*(rhocp*u(i  ,j,k)+drhocp*psiflx_x(i  ,j,k))*(s(i+1,j,k)+s(i,j,k))
            usim = 0.5*(rhocp*u(i-1,j,k)+drhocp*psiflx_x(i-1,j,k))*(s(i-1,j,k)+s(i,j,k))
            vsjp = 0.5*(rhocp*v(i,j  ,k)+drhocp*psiflx_y(i,j  ,k))*(s(i,j+1,k)+s(i,j,k))
            vsjm = 0.5*(rhocp*v(i,j-1,k)+drhocp*psiflx_y(i,j-1,k))*(s(i,j-1,k)+s(i,j,k))
            wskp = 0.5*(rhocp*w(i,j,k  )+drhocp*psiflx_z(i,j,k  ))*(s(i,j,k+1)+s(i,j,k))
            wskm = 0.5*(rhocp*w(i,j,k-1)+drhocp*psiflx_z(i,j,k-1))*(s(i,j,k-1)+s(i,j,k))
          else
            usip = 0.; usim = 0.
            vsjp = 0.; vsjm = 0.
            wskp = 0.; wskm = 0.
          end if
          !
          psixp = 0.5*(psi(i+1,j,k)+psi(i  ,j,k))
          psixm = 0.5*(psi(i  ,j,k)+psi(i-1,j,k))
          psiyp = 0.5*(psi(i,j+1,k)+psi(i,j  ,k))
          psiym = 0.5*(psi(i,j  ,k)+psi(i,j-1,k))
          psizp = 0.5*(psi(i,j,k+1)+psi(i,j,k  ))
          psizm = 0.5*(psi(i,j,k  )+psi(i,j,k-1))
          kaxp = ka+dka*psixp
          kaxm = ka+dka*psixm
          kayp = ka+dka*psiyp
          kaym = ka+dka*psiym
          kazp = ka+dka*psizp
          kazm = ka+dka*psizm
          if(is_harm) then
            kaxp = 1._rp/kaxp
            kaxm = 1._rp/kaxm
            kayp = 1._rp/kayp
            kaym = 1._rp/kaym
            kazp = 1._rp/kazp
            kazm = 1._rp/kazm
          end if
          !
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxi
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxi
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyi
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyi
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          !
          dsdt(i,j,k) = dxi*(     -usip + usim ) + &
                        dyi*(     -vsjp + vsjm ) + &
                        dzfi(k)*( -wskp + wskm ) + &
                        ( (kaxp*dsdxp-kaxm*dsdxm)*dxi + &
                          (kayp*dsdyp-kaym*dsdym)*dyi + &
                          (kazp*dsdzp-kazm*dsdzm)*dzfi(k) )
        end do
      end do
    end do
  end subroutine scal_ad
  !
  subroutine scal_a_cons(n,nh_v,nh_s,nh_psi,nh_ifc,nh_g,dli,dzfi,rhocp12,psi,psiflx_x,psiflx_y,psiflx_z,ua,va,wa,sa,factor,s)
    !
    ! advection term of the scalar in conservative form, consistent with the transport of the geometric
    ! volume-of-fluid method: integrated over the time `factor` and added to the scalar, divided by rho*cp
    !
    ! the face fluxes of rho*cp are rhocp_2*ua + (rhocp_1-rhocp_2)*psiflx, with psiflx the (time-averaged)
    ! phase-1 volume fluxes of the volume-of-fluid step and ua its advecting velocity; sa is the transported
    ! scalar, interpolated centrally to the faces
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: nh_v,nh_s,nh_psi,nh_ifc,nh_g
    real(rp), intent(in   ), dimension(3) :: dli
    real(rp), intent(in   ), dimension(1-nh_g:) :: dzfi
    real(rp), intent(in   ), dimension(2) :: rhocp12
    real(rp), intent(in   ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in   ), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:) :: psiflx_x,psiflx_y,psiflx_z
    real(rp), intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: ua,va,wa
    real(rp), intent(in   ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: sa
    real(rp), intent(in   ) :: factor
    real(rp), intent(inout), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm,rhocp,drhocp,dxi,dyi
    integer :: i,j,k
    !
    rhocp = rhocp12(2); drhocp = rhocp12(1)-rhocp12(2)
    dxi = dli(1)
    dyi = dli(2)
    !$acc parallel loop collapse(3) default(present) private(usip,usim,vsjp,vsjm,wskp,wskm) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          usip = 0.5*(rhocp*ua(i  ,j,k)+drhocp*psiflx_x(i  ,j,k))*(sa(i+1,j,k)+sa(i,j,k))
          usim = 0.5*(rhocp*ua(i-1,j,k)+drhocp*psiflx_x(i-1,j,k))*(sa(i-1,j,k)+sa(i,j,k))
          vsjp = 0.5*(rhocp*va(i,j  ,k)+drhocp*psiflx_y(i,j  ,k))*(sa(i,j+1,k)+sa(i,j,k))
          vsjm = 0.5*(rhocp*va(i,j-1,k)+drhocp*psiflx_y(i,j-1,k))*(sa(i,j-1,k)+sa(i,j,k))
          wskp = 0.5*(rhocp*wa(i,j,k  )+drhocp*psiflx_z(i,j,k  ))*(sa(i,j,k+1)+sa(i,j,k))
          wskm = 0.5*(rhocp*wa(i,j,k-1)+drhocp*psiflx_z(i,j,k-1))*(sa(i,j,k-1)+sa(i,j,k))
          s(i,j,k) = s(i,j,k) + factor*(dxi*(-usip+usim) + dyi*(-vsjp+vsjm) + dzfi(k)*(-wskp+wskm))/(rhocp+drhocp*psi(i,j,k))
        end do
      end do
    end do
  end subroutine scal_a_cons
  !
  subroutine cmpt_scalflux(n,nh_s,nh_psi,nh_g,is_bound,l,dli,dzci,dzfi,ka12,psi,s,flux)
    use mpi
    implicit none
    integer , intent(in ), dimension(3) :: n
    integer , intent(in) :: nh_s,nh_psi,nh_g
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(1-nh_g:)    :: dzci,dzfi
    real(rp), intent(in ), dimension(2) :: ka12
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in ), dimension(1-nh_s:,1-nh_s:,1-nh_s:) :: s
    real(rp), intent(out), dimension(3) :: flux
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(rp) :: kaxp,kaxm,kayp,kaym,kazp,kazm
    real(rp) :: flux_x,flux_y,flux_z
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz
    real(rp) :: ka,dka,psip
    logical  :: is_harm
    integer :: ierr
    !
    is_harm = is_harmonic_interpolation
    if(is_harm) then
      ! N.b.: care should be taken if one of the conductivities is zero.
      ka = 1._rp/ka12(2); dka = 1._rp/ka12(1)-1._rp/ka12(2)
    else
      ka = ka12(2); dka = ka12(1)-ka12(2)
    end if
    !
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    flux_x = 0._rp
    !$acc data copy(flux_x) async(1)
    if(is_bound(0,1)) then
      i = 0
      !$acc parallel loop collapse(2) default(present) private(dsdxp,kaxp,psip) reduction(+:flux_x) async(1)
      do k=1,nz
        do j=1,ny
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))*dxi
          psip = 0.5*(psi(i+1,j,k)+psi(i,j,k))
          kaxp = ka+dka*psip
          if(is_harm) kaxp = 1._rp/kaxp
          flux_x = flux_x + kaxp*dsdxp/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    if(is_bound(1,1)) then
      i = nx+1
      !$acc parallel loop collapse(2) default(present) private(dsdxm,kaxm,psip) reduction(+:flux_x) async(1)
      do k=1,nz
        do j=1,ny
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))*dxi
          psip = 0.5*(psi(i,j,k)+psi(i-1,j,k))
          kaxm = ka+dka*psip
          if(is_harm) kaxm = 1._rp/kaxm
          flux_x = flux_x - kaxm*dsdxm/(dyi*dzfi(k)*ly*lz)
        end do
      end do
    end if
    !$acc end data
    flux_y = 0._rp
    !$acc data copy(flux_y) async(1)
    if(is_bound(0,2)) then
      j = 0
      !$acc parallel loop collapse(2) default(present) private(dsdyp,kayp,psip) reduction(+:flux_y) async(1)
      do k=1,nz
        do i=1,nx
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))*dyi
          psip = 0.5*(psi(i,j+1,k)+psi(i,j,k))
          kayp = ka+dka*psip
          if(is_harm) kayp = 1._rp/kayp
          flux_y = flux_y + kayp*dsdyp/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    if(is_bound(1,2)) then
      j = ny+1
      !$acc parallel loop collapse(2) default(present) private(dsdym,kaym,psip) reduction(+:flux_y) async(1)
      do k=1,nz
        do i=1,nx
          dsdym = (s(i,j  ,k)-s(i,j-1,k))*dyi
          psip = 0.5*(psi(i,j,k)+psi(i,j-1,k))
          kaym = ka+dka*psip
          if(is_harm) kaym = 1._rp/kaym
          flux_y = flux_y - kaym*dsdym/(dxi*dzfi(k)*lx*lz)
        end do
      end do
    end if
    !$acc end data
    flux_z = 0._rp
    !$acc data copy(flux_z) async(1)
    if(is_bound(0,3)) then
      k = 0
      !$acc parallel loop collapse(2) default(present) private(dsdzp,kazp,psip) reduction(+:flux_z) async(1)
      do j=1,ny
        do i=1,nx
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))*dzci(k  )
          psip = 0.5*(psi(i,j,k+1)+psi(i,j,k))
          kazp = ka+dka*psip
          if(is_harm) kazp = 1._rp/kazp
          flux_z = flux_z + kazp*dsdzp/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    if(is_bound(1,3)) then
      k = nz+1
      !$acc parallel loop collapse(2) default(present) private(dsdzm,kazm,psip) reduction(+:flux_z) async(1)
      do j=1,ny
        do i=1,nx
          dsdzm = (s(i,j,k  )-s(i,j,k-1))*dzci(k-1)
          psip = 0.5*(psi(i,j,k)+psi(i,j,k-1))
          kazm = ka+dka*psip
          if(is_harm) kazm = 1._rp/kazm
          flux_z = flux_z - kazm*dsdzm/(dxi*dyi*lx*ly)
        end do
      end do
    end if
    !$acc end data
    !$acc wait(1)
    flux(:) = [flux_x,flux_y,flux_z]
    call MPI_ALLREDUCE(MPI_IN_PLACE,flux,3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine cmpt_scalflux
end module mod_scal
