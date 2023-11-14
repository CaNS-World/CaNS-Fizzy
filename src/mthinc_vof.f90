! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_two_fluid_common
  use mod_types
  !
  ! TODO: SPLIT MTHINC-SPECIFIC FROM COMMON ROUTINES
  !
  implicit none
  private
  public update_property
  real(rp), parameter :: limit = 1.e-8, beta = 2._rp
  contains
  subroutine update_property(nh,prop12,psi,p)
    !
    ! updates a certain transport property based on a one-fluid formulation
    !
    implicit none
    integer , intent(in ), dimension(3)        :: nh
    real(rp), intent(in ), dimension(2)        :: prop12
    real(rp), intent(in ), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: psi
    real(rp), intent(out), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: p
    real(rp) :: prop1,prop2
    prop1 = prop12(1)
    prop2 = prop12(2)
    !$acc kernels default(present) async(1)
       p(:,:,:) = psi(:,:,:)*prop1+(1.-psi(:,:,:))*prop2
    !$acc end kernels
  end subroutine update_property
  !
  subroutine clip_field(nh,minmax,p)
    !
    ! clips a scalar field to certain maximum and minimum values
    !
    implicit none
    integer , intent(in   ), dimension(3) :: nh
    real(rp), intent(in   ), dimension(2) :: minmax
    real(rp), intent(inout), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: p
    real(rp) :: rmin,rmax
    integer :: n1,n2,n3,i,j,k
    !
    n1 = size(p,1)-2*nh(1)
    n2 = size(p,2)-2*nh(2)
    n3 = size(p,3)-2*nh(3)
    rmin = minmax(1); rmax = minmax(2)
    !
    !$acc parallel loop collapse(3) default(present) private(rmin,rmax) async(1)
    do concurrent(k=0:n3-1,j=0:n2-1,i=0:n1-1)
      p(i,j,k) = min(max(rmin,p(i,j,k)),rmax)
    end do
  end subroutine clip_field
  !
#if 1
  subroutine mthinc_updt_vof(cbc,bc,idir,n,dt,dl,dli,dzc,dzf,dzci,dzfi,vel,psi,normx,normy,normz,kappa,lij,d)
    !
    ! updates mthinc parameters after advecting the vof fiel:
    !   - clips the vof to force its boundness and employs boundary conditions
    !   - updates normals, curvature, mthinc lij tensor, and the d parameter
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    real(rp), intent(in   ), dimension(0:1,3,3) :: bc
    integer , intent(in   )               :: idir
    integer , intent(in   ), dimension(3) :: n
    real(rp), intent(in   )               :: dt
    real(rp), intent(in   ), dimension(3) :: dl,dli
    real(rp), intent(in   ), dimension(0:) :: dzc,dzf,dzci,dzfi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: vel
    real(rp), intent(inout), dimension(0:,0:,0:) :: psi
    real(rp), intent(out  ), dimension(0:,0:,0:) :: normx,normy,normz,kappa
    real(rp), intent(out  ), dimension(0:,0:,0:,1:) :: lij
    real(rp), intent(out  ), dimension(0:,0:,0:) :: d
    real(rp) :: rmin,rmax
    integer :: q
    !
    call boundp(cbc,n,bc,dl,dzc,dzf,psi)
    call clip_field([1,1,1],[0._rp,1._rp],psi)
    call mthinc_cmpt_norm_and_curv(n,dl,dli,dzc,dzf,dzci,dzfi,psi,normx,normy,normz,kappa,lij)
    call boundp(cbc,n,bc,dl,dzc,dzf,normx)
    call boundp(cbc,n,bc,dl,dzc,dzf,normy)
    call boundp(cbc,n,bc,dl,dzc,dzf,normz)
    call boundp(cbc,n,bc,dl,dzc,dzf,kappa)
    do q=1,6
      call boundp(cbc,n,bc,dl,dzc,dzf,lij(:,:,:,q))
    enddo
    !
    call mthinc_cmpt_d_param(idir,n,dt,dli,dzci,vel,psi,normx,normy,normz,lij,d)
    call boundp(cbc,n,bc,dl,dzc,dzf,d)
  end subroutine mthinc_updt_vof
  !
#endif
  subroutine mthinc_avect_vof(cbc,bc,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,v,w,psi,normx,normy,normz,kappa,lij,d)
    !
    ! updates mthinc parameters after advecting the vof fiel:
    !   - clips the vof to force its boundness and employs boundary conditions
    !   - updates normals, curvature, mthinc lij tensor, and the d parameter
    !
    ! TODO: allocate here lij with save attribute? since we do not need it elsewhere?
    integer, save :: iter = 0
    !
    iter = iter + 1
    idir = mod(iter-1,3)+1
    !
    select case(idir)
    case(1)
      call mthinc_cmpt_vof_flux(1,n,dt,dli,dzci,u,psi  ,normx,normy,normz,lij,d,flux,psi_1)
      call mthinc_updt_vof(cbc,bc,idir,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,psi_1,normx,normy,normz,kappa,lij,d)
      call mthinc_cmpt_vof_flux(2,n,dt,dli,dzci,v,psi_1,normx,normy,normz,lij,d,flux,psi_2)
      call mthinc_updt_vof(cbc,bc,idir,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,psi,normx,normy,normz,kappa,lij,d)
      call mthinc_cmpt_vof_flux(3,n,dt,dli,dzci,w,psi_2,normx,normy,normz,lij,d,flux,psi  )
      iidir = 3
    case(2)
      call mthinc_cmpt_vof_flux(2,n,dt,dli,dzci,v,psi  ,normx,normy,normz,lij,d,flux,psi_1)
      call mthinc_updt_vof(cbc,bc,2,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,psi,normx,normy,normz,kappa,lij,d)
      call mthinc_cmpt_vof_flux(3,n,dt,dli,dzci,w,psi_1,normx,normy,normz,lij,d,flux,psi_2)
      call mthinc_updt_vof(cbc,bc,3,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,psi,normx,normy,normz,kappa,lij,d)
      call mthinc_cmpt_vof_flux(1,n,dt,dli,dzci,u,psi_2,normx,normy,normz,lij,d,flux,psi  )
      iidir = 1
    case(3)
      call mthinc_cmpt_vof_flux(3,n,dt,dli,dzci,w,psi  ,normx,normy,normz,lij,d,flux,psi_1)
      call mthinc_updt_vof(cbc,bc,3,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,psi,normx,normy,normz,kappa,lij,d)
      call mthinc_cmpt_vof_flux(1,n,dt,dli,dzci,u,psi_1,normx,normy,normz,lij,d,flux,psi_2)
      call mthinc_updt_vof(cbc,bc,1,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,psi,normx,normy,normz,kappa,lij,d)
      call mthinc_cmpt_vof_flux(2,n,dt,dli,dzci,v,psi_2,normx,normy,normz,lij,d,flux,psi  )
      iidir = 2
    end select
      !
      ! divergence correction step
      !
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psi(i,j,k) = psi(i,j,k) - dt*(psi_1(i,j,k)*(u(i,j,k)-u(i-1,j,k))*dli(1)  + &
                                          psi_2(i,j,k)*(v(i,j,k)-v(i,j-1,k))*dli(2)  + &
                                            psi(i,j,k)*(w(i,j,k)-w(i,j,k-1))*dzfi(k) )
#if 1
            !
            ! volume deflation/inflation step for problems with non-zero velocity divergence
            !
            div = (u(i,j,k)-u(i-1,j,k))*dli(1) + &
                  (v(i,j,k)-v(i,j-1,k))*dli(2) + &
                  (w(i,j,k)-w(i,j,k-1))*dzfi(k)
            vof(i,j,k) = vof(i,j,k)/(1.-dt*div)
#endif
          enddo
        enddo
      enddo
      call mthinc_updt_vof(cbc,bc,iidir,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,psi,normx,normy,normz,kappa,lij,d)
    end subroutine mthinc_avect_vof(cbc,bc,n,dt,dl,dli,dzc,dzf,dzci,dzfi,u,v,w,psi,normx,normy,normz,kappa,lij,d)
#endif
  subroutine mthinc_cmpt_norm_and_curv(n,dl,dli,dzc,dzf,dzci,dzfi,psi,normx,normy,normz,kappa,lij)
    !
    ! computes the normals, curvature and lij tensor using finite-differences based on Youngs method
    ! see Li et al., JCP 231 2328–2358
    !
    implicit none
    integer , parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(3)           :: dl,dli
    real(rp), intent(in ), dimension(0:)          :: dzc,dzf,dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:)    :: psi
    real(rp), intent(out), dimension(0:,0:,0:)    :: normx,normy,normz,kappa
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: lij
    real(rp) :: norm
    logical , save :: is_first = .true.
    real(rp), allocatable, dimension(:), save :: mx,my,mz
    integer  :: i,j,k,q
    !
    if(is_first) then
      is_first = .false.
      allocate(mx(8),my(8),mz(8))
      !$acc enter data create(mx,my,mz)
    end if
    !
    !$acc parallel loop collapse(3) default(present) private(mx,my,mz,norm) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          mx(1) = 0.25*((psi(i+1,j  ,k  )+psi(i+1,j+1,k  )+psi(i+1,j  ,k+1)+psi(i+1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j+1,k+1)))*dli(1)
          mx(2) = 0.25*((psi(i+1,j  ,k  )+psi(i+1,j-1,k  )+psi(i+1,j  ,k+1)+psi(i+1,j-1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j-1,k+1)))*dli(1)
          mx(3) = 0.25*((psi(i+1,j  ,k  )+psi(i+1,j+1,k  )+psi(i+1,j  ,k-1)+psi(i+1,j+1,k-1)) - &
                        (psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j+1,k-1)))*dli(1)
          mx(4) = 0.25*((psi(i+1,j  ,k  )+psi(i+1,j-1,k  )+psi(i+1,j  ,k-1)+psi(i+1,j-1,k-1)) - &
                        (psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j-1,k-1)))*dli(1)
          mx(5) = 0.25*((psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j+1,k+1)) - &
                        (psi(i-1,j  ,k  )+psi(i-1,j+1,k  )+psi(i-1,j  ,k+1)+psi(i-1,j+1,k+1)))*dli(1)
          mx(6) = 0.25*((psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j-1,k+1)) - &
                        (psi(i-1,j  ,k  )+psi(i-1,j-1,k  )+psi(i-1,j  ,k+1)+psi(i-1,j-1,k+1)))*dli(1)
          mx(7) = 0.25*((psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j+1,k-1)) - &
                        (psi(i-1,j  ,k  )+psi(i-1,j+1,k  )+psi(i-1,j  ,k-1)+psi(i-1,j+1,k-1)))*dli(1)
          mx(8) = 0.25*((psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j-1,k-1)) - &
                        (psi(i-1,j  ,k  )+psi(i-1,j-1,k  )+psi(i-1,j  ,k-1)+psi(i-1,j-1,k-1)))*dli(1)
          !
          my(1) = 0.25*((psi(i  ,j+1,k  )+psi(i+1,j+1,k  )+psi(i  ,j+1,k+1)+psi(i+1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)))*dli(2)
          my(2) = 0.25*((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)) - &
                        (psi(i  ,j-1,k  )+psi(i+1,j-1,k  )+psi(i  ,j-1,k+1)+psi(i+1,j-1,k+1)))*dli(2)
          my(3) = 0.25*((psi(i  ,j+1,k  )+psi(i+1,j+1,k  )+psi(i  ,j+1,k-1)+psi(i+1,j+1,k-1)) - &
                        (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)))*dli(2)
          my(4) = 0.25*((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)) - &
                        (psi(i  ,j-1,k  )+psi(i+1,j-1,k  )+psi(i  ,j-1,k-1)+psi(i+1,j-1,k-1)))*dli(2)
          my(5) = 0.25*((psi(i  ,j+1,k  )+psi(i-1,j+1,k  )+psi(i  ,j+1,k+1)+psi(i-1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)))*dli(2)
          my(6) = 0.25*((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)) - &
                        (psi(i  ,j-1,k  )+psi(i-1,j-1,k  )+psi(i  ,j-1,k+1)+psi(i-1,j-1,k+1)))*dli(2)
          my(7) = 0.25*((psi(i  ,j+1,k  )+psi(i-1,j+1,k  )+psi(i  ,j+1,k-1)+psi(i-1,j+1,k-1)) - &
                        (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)))*dli(2)
          my(8) = 0.25*((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)) - &
                        (psi(i  ,j-1,k  )+psi(i-1,j-1,k  )+psi(i  ,j-1,k-1)+psi(i-1,j-1,k-1)))*dli(2)
          !
          mz(1) = 0.25*((psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)+psi(i  ,j+1,k+1)+psi(i+1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j+1,k  )+psi(i+1,j+1,k  )))*dzci(k  )
          mz(2) = 0.25*((psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)+psi(i  ,j-1,k+1)+psi(i+1,j-1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j-1,k  )+psi(i+1,j-1,k  )))*dzci(k  )
          mz(3) = 0.25*((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j+1,k  )+psi(i+1,j+1,k  )) - &
                        (psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)+psi(i  ,j+1,k-1)+psi(i+1,j+1,k-1)))*dzci(k  )
          mz(4) = 0.25*((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j-1,k  )+psi(i+1,j-1,k  )) - &
                        (psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)+psi(i  ,j-1,k-1)+psi(i+1,j-1,k-1)))*dzci(k  )
          mz(5) = 0.25*((psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)+psi(i  ,j+1,k+1)+psi(i-1,j+1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j+1,k  )+psi(i-1,j+1,k  )))*dzci(k+1)
          mz(6) = 0.25*((psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)+psi(i  ,j-1,k+1)+psi(i-1,j-1,k+1)) - &
                        (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j-1,k  )+psi(i-1,j-1,k  )))*dzci(k+1)
          mz(7) = 0.25*((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j+1,k  )+psi(i-1,j+1,k  )) - &
                        (psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)+psi(i  ,j+1,k-1)+psi(i-1,j+1,k-1)))*dzci(k  )
          mz(8) = 0.25*((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j-1,k  )+psi(i-1,j-1,k  )) - &
                        (psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)+psi(i  ,j-1,k-1)+psi(i-1,j-1,k-1)))*dzci(k  )
          !$acc loop seq
          do q=1,8
            norm = sqrt(mx(q)**2+my(q)**2+mz(q)**2)+eps
            mx(q) = mx(q)/norm
            my(q) = my(q)/norm
            mz(q) = mz(q)/norm
          end do
          !
          ! compute the normal vector
          !
          normx(i,j,k) = .125*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          normy(i,j,k) = .125*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          normz(i,j,k) = .125*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
#if 0
          norm = sqrt(normx(i,j,k)**2+normy(i,j,k)**2+normz(i,j,k)**2)+eps
          normx(i,j,k) = normx(i,j,k)/norm
          normy(i,j,k) = normy(i,j,k)/norm
          normz(i,j,k) = normz(i,j,k)/norm
#endif
          !
          ! compute the curvature tensor
          !
          lij(i,j,k,1) = ((mx(1)+mx(2)+mx(3)+mx(4))-(mx(5)+mx(6)+mx(7)+mx(8)))*dl( 1)*0.25
          lij(i,j,k,2) = ((my(1)+my(3)+my(5)+my(7))-(my(2)+my(4)+my(6)+my(8)))*dl( 2)*0.25
          lij(i,j,k,3) = ((mz(1)+mz(2)+mz(5)+mz(6))-(mz(3)+mz(4)+mz(7)+mz(8)))*dzf(k)*0.25
          kappa(i,j,k) = -(lij(i,j,k,1)*dli(1)**2+lij(i,j,k,2)*dli(2)**2+lij(i,j,k,3)*dzfi(k)**2) ! curvature
          !
          lij(i,j,k,4) = (((my(1)+my(2)+my(3)+my(4))-(my(5)+my(6)+my(7)+my(8)))*dl( 2) + &
                          ((mx(1)+mx(3)+mx(5)+mx(7))-(mx(2)+mx(4)+mx(6)+mx(8)))*dl( 1))*0.125
          lij(i,j,k,5) = (((mz(1)+mz(2)+mz(3)+mz(4))-(mz(5)+mz(6)+mz(7)+mz(8)))*dzf(k) + &
                          ((mx(1)+mx(2)+mx(5)+mx(6))-(mx(3)+mx(4)+mx(7)+mx(8)))*dl( 1))*0.125
          lij(i,j,k,6) = (((mz(1)+mz(3)+mz(5)+mz(7))-(mz(2)+mz(4)+mz(6)+mz(8)))*dzf(k) + &
                          ((my(1)+my(2)+my(5)+my(6))-(my(3)+my(4)+my(7)+my(8)))*dl( 2))*0.125
        enddo
      enddo
    enddo
  end subroutine mthinc_cmpt_norm_and_curv
  !
#if 1
  subroutine mthinc_cmpt_vof_flux(idir,n,dt,dli,dzci,vel,psi,normx,normy,normz,lij,d,flux,psi_out)
    !
    ! computes the mthinc vof fluxes
    ! see Li et al., JCP 231 2328–2358
    !
    use mod_common_cudecomp, only: buf => work
    implicit none
    real(rp), parameter :: gq = 1._rp/sqrt(3._rp)
    integer , parameter :: eps = epsilon(1._rp)
    !
    integer , intent(in )                         :: idir
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in )                         :: dt
    real(rp), intent(in ), dimension(3)           :: dli
    real(rp), intent(in ), dimension(0:)          :: dzci
    real(rp), intent(in ), dimension(0:,0:,0:)    :: vel,psi,normx,normy,normz
    real(rp), intent(in ), dimension(0:,0:,0:,1:) :: lij  ! curvature tensor
    real(rp), intent(in ), dimension(0:,0:,0:)    :: d
    real(rp), intent(out), dimension(0:,0:,0:)    :: psi_out ! flux,psi
    real(rp), pointer, contiguous, dimension(:,:,:) :: flux
    integer  :: i,j,k,ii,jj,kk,q
    integer  :: is_i,is_j,is_k
    real(rp) :: rvel,rsgn
    real(rp) :: dx_int,dy_int,dz_int,xl,xu,yl,yu,zl,zu
    real(rp) :: rnormx,rnormy,rnormz,rnorm_max
    integer  :: is_nx,is_ny,is_nz,cx,cy,cz
    real(rp) :: lxx,lxy,lxz,lyy,lyz,lzz,a001,a010,a100,a011,a101,a110,a002,a020,a200
    real(rp) :: xl_int_a,xu_int_a,x_int_nm,x_int_np, &
                yl_int_a,yu_int_a,y_int_nm,y_int_np, &
                zl_int_a,zu_int_a,z_int_nm,z_int_np
    real(rp) :: dd
    real(rp) :: p_xy_mm,p_xy_mm_l,p_xy_mm_u, &
                p_xy_mp,p_xy_mp_l,p_xy_mp_u, &
                p_xy_pm,p_xy_pm_l,p_xy_pm_u, &
                p_xy_pp,p_xy_pp_l,p_xy_pp_u, &
                p_yz_mm,p_yz_mm_l,p_yz_mm_u, &
                p_yz_mp,p_yz_mp_l,p_yz_mp_u, &
                p_yz_pm,p_yz_pm_l,p_yz_pm_u, &
                p_yz_pp,p_yz_pp_l,p_yz_pp_u, &
                p_zx_mm,p_zx_mm_l,p_zx_mm_u, &
                p_zx_mp,p_zx_mp_l,p_zx_mp_u, &
                p_zx_pm,p_zx_pm_l,p_zx_pm_u, &
                p_zx_pp,p_zx_pp_l,p_zx_pp_u
    real(rp) :: dl_a,rnorm,factor
    real(rp) :: rsum,is_x,is_y,is_z
    flux(0:n(1),0:n(2),0:n(3)) => buf(1:product(n(:)+1))
    !
    is_i = 0; if(idir == 1) is_i = 1
    is_j = 0; if(idir == 2) is_j = 1
    is_k = 0; if(idir == 3) is_k = 1
    do k=0,n(3)
      do j=0,n(2)
        do i=0,n(1)
          !
          ! check velocity sign
          !
          rvel = vel(i,j,k)
          rsgn = sign(1._rp,rvel)
          q = nint(0.5*(-rsgn+1.))
          ii = i + is_i*q
          jj = j + is_j*q
          kk = k + is_k*q
          !
          ! compute integration interval
          !
          dx_int = is_i*rvel*dt*dli(1)
          dy_int = is_j*rvel*dt*dli(2)
          dz_int = is_k*rvel*dt*dzci(kk)
          xl = q-dx_int
          xu = q
          yl = q-dy_int
          yu = q
          zl = q-dz_int
          zu = q
          !
          ! determine integration order
          !
          rnormx = abs(normx(ii,jj,kk))+eps; rnormy = abs(normy(ii,jj,kk))+eps; rnormz = abs(normz(ii,jj,kk))+eps
          rnorm_max = max(rnormx,rnormy,rnormz)
          is_nx = int(rnormx/rnorm_max+eps)
          is_ny = int(rnormy/rnorm_max+eps)
          is_nz = int(rnormz/rnorm_max+eps)
          cx = 1-is_nx
          cy = 1-is_ny
          cz = 1-is_nz
          rnormx = normx(ii,jj,kk); rnormy = normy(ii,jj,kk); rnormz = normz(ii,jj,kk)
          lxx  = lij(ii,jj,kk,1)
          lyy  = lij(ii,jj,kk,2)
          lzz  = lij(ii,jj,kk,3)
          lxy  = lij(ii,jj,kk,4)
          lxz  = lij(ii,jj,kk,5)
          lyz  = lij(ii,jj,kk,6)
          a200 = 0.5*cx*lxx
          a020 = 0.5*cy*lyy
          a002 = 0.5*cz*lzz
          a110 = cx*cy*lxy
          a101 = cx*cz*lxz
          a011 = cy*cz*lyz
          a100 = rnormx-(a200+0.5*(a110+a101))
          a010 = rnormy-(a020+0.5*(a110+a011))
          a001 = rnormz-(a002+0.5*(a101+a011))
          if( (psi(ii,jj,kk) <= limit).or.(psi(ii,jj,kk) >= 1.-limit) ) then
            flux(i,j,k) = psi(ii,jj,kk)*(xu-xl)*(yu-yl)*(zu-zl)
          else
            !
            xl_int_a  = is_nx*xl
            xu_int_a  = is_nx*xu
            x_int_nm = (1-is_nx)*0.5*((xl+xu)-gq*(xu-xl))
            x_int_np = (1-is_nx)*0.5*((xl+xu)+gq*(xu-xl))
            yl_int_a  = is_ny*yl
            yu_int_a  = is_ny*yu
            y_int_nm = (1-is_ny)*0.5*((yl+yu)-gq*(yu-yl))
            y_int_np = (1-is_ny)*0.5*((yl+yu)+gq*(yu-yl))
            zl_int_a  = is_nz*zl
            zu_int_a  = is_nz*zu
            z_int_nm = (1-is_nz)*0.5*((zl+zu)-gq*(zu-zl))
            z_int_np = (1-is_nz)*0.5*((zl+zu)+gq*(zu-zl))
            dd = d(ii,jj,kk)
            !
            ! -1,-1 quadrature point
            !
            p_xy_mm = a200*x_int_nm**2+a020*y_int_nm**2+a110*x_int_nm*y_int_nm+a100*x_int_nm+a010*y_int_nm
            p_xy_mm_l = p_xy_mm_l + a001*zl_int_a
            p_xy_mm_u = p_xy_mm_l + a001*zu_int_a
            p_yz_mm = a020*y_int_nm**2+a002*z_int_nm**2+a011*y_int_nm*z_int_nm+a010*y_int_nm+a001*z_int_nm
            p_yz_mm_l = p_yz_mm_l + a100*xl_int_a
            p_yz_mm_u = p_yz_mm_l + a100*xu_int_a
            p_zx_mm = a002*z_int_nm**2+a200*x_int_nm**2+a101*z_int_nm*x_int_nm+a001*z_int_nm+a100*x_int_nm
            p_zx_mm_l = p_zx_mm_l + a010*yl_int_a
            p_zx_mm_u = p_zx_mm_l + a010*yu_int_a
            !
            ! +1,-1 quadrature point
            !
            p_xy_pm = a200*x_int_np**2+a020*y_int_nm**2+a110*x_int_np*y_int_nm+a100*x_int_np+a010*y_int_nm
            p_xy_pm_l = p_xy_pm_l + a001*zl_int_a
            p_xy_pm_u = p_xy_pm_l + a001*zu_int_a
            p_yz_pm = a020*y_int_np**2+a002*z_int_nm**2+a011*y_int_np*z_int_nm+a010*y_int_np+a001*z_int_nm
            p_yz_pm_l = p_yz_pm_l + a100*xl_int_a
            p_yz_pm_u = p_yz_pm_l + a100*xu_int_a
            p_zx_pm = a002*z_int_np**2+a200*x_int_nm**2+a101*z_int_np*x_int_nm+a001*z_int_np+a100*x_int_nm
            p_zx_pm_l = p_zx_pm_l + a010*yl_int_a
            p_zx_pm_u = p_zx_pm_l + a010*yu_int_a
            !
            ! -1,+1 quadrature point
            !
            p_xy_mp = a200*x_int_nm**2+a020*y_int_np**2+a110*x_int_nm*y_int_np+a100*x_int_nm+a010*y_int_np
            p_xy_mp_l = p_xy_mp_l + a001*zl_int_a
            p_xy_mp_u = p_xy_mp_l + a001*zu_int_a
            p_yz_mp = a020*y_int_nm**2+a002*z_int_np**2+a011*y_int_nm*z_int_np+a010*y_int_nm+a001*z_int_np
            p_yz_mp_l = p_yz_mp_l + a100*xl_int_a
            p_yz_mp_u = p_yz_mp_l + a100*xu_int_a
            p_zx_mp = a002*z_int_nm**2+a200*x_int_np**2+a101*z_int_nm*x_int_np+a001*z_int_nm+a100*x_int_np
            p_zx_mp_l = p_zx_mp_l + a010*yl_int_a
            p_zx_mp_u = p_zx_mp_l + a010*yu_int_a
            !
            ! +1,+1 quadrature point
            !
            p_xy_pp = a200*x_int_np**2+a020*y_int_np**2+a110*x_int_np*y_int_np+a100*x_int_np+a010*y_int_np
            p_xy_pp_l = p_xy_pp_l + a001*zl_int_a
            p_xy_pp_u = p_xy_pp_l + a001*zu_int_a
            p_yz_pp = a020*y_int_np**2+a002*z_int_np**2+a011*y_int_np*z_int_np+a010*y_int_np+a001*z_int_np
            p_yz_pp_l = p_yz_pp_l + a100*xl_int_a
            p_yz_pp_u = p_yz_pp_l + a100*xu_int_a
            p_zx_pp = a002*z_int_np**2+a200*x_int_np**2+a101*z_int_np*x_int_np+a001*z_int_np+a100*x_int_np
            p_zx_pp_l = p_zx_pp_l + a010*yl_int_a
            p_zx_pp_u = p_zx_pp_l + a010*yu_int_a
            !
            ! time to compute the flux
            !
            dl_a = ( xu_int_a-xl_int_a + yu_int_a-yl_int_a + zu_int_a-zl_int_a )
            rnorm = rnormx*is_x+rnormy*is_y+rnormz*is_z
            factor = 1./((rnorm+eps)*beta)
            rsum =  dl_a + factor*(log(cosh(beta*(p_yz_mm_u+dd))/cosh(beta*(p_yz_mm_l+dd)))*is_x + &
                                   log(cosh(beta*(p_yz_pm_u+dd))/cosh(beta*(p_yz_pm_l+dd)))*is_x + &
                                   log(cosh(beta*(p_yz_mp_u+dd))/cosh(beta*(p_yz_mp_l+dd)))*is_x + &
                                   log(cosh(beta*(p_yz_pp_u+dd))/cosh(beta*(p_yz_pp_l+dd)))*is_x + &
                                   log(cosh(beta*(p_zx_mm_u+dd))/cosh(beta*(p_zx_mm_l+dd)))*is_y + &
                                   log(cosh(beta*(p_zx_pm_u+dd))/cosh(beta*(p_zx_pm_l+dd)))*is_y + &
                                   log(cosh(beta*(p_zx_mp_u+dd))/cosh(beta*(p_zx_mp_l+dd)))*is_y + &
                                   log(cosh(beta*(p_zx_pp_u+dd))/cosh(beta*(p_zx_pp_l+dd)))*is_y + &
                                   log(cosh(beta*(p_xy_mm_u+dd))/cosh(beta*(p_xy_mm_l+dd)))*is_z + &
                                   log(cosh(beta*(p_xy_pm_u+dd))/cosh(beta*(p_xy_pm_l+dd)))*is_z + &
                                   log(cosh(beta*(p_xy_mp_u+dd))/cosh(beta*(p_xy_mp_l+dd)))*is_z + &
                                   log(cosh(beta*(p_xy_pp_u+dd))/cosh(beta*(p_xy_pp_l+dd)))*is_z)
            flux(i,j,k) = (xu-xl)*(yu-yl)*(zu-zl)/dl_a*rsum
          end if
        end do
      end do
    end do
    select case(idir)
    case(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psi_out(i,j,k) = (psi(i,j,k)-(flux(i,j,k)-flux(i-1,j,k))*dli(1))/(1.-dt*dli(1)*(vel(i,j,k)-vel(i-1,j,k)))
          enddo
        enddo
      enddo
    case(2)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psi_out(i,j,k) = (psi(i,j,k)-(flux(i,j,k)-flux(i,j-1,k))*dli(2))/(1.-dt*dli(2)*(vel(i,j,k)-vel(i,j-1,k)))
          enddo
        enddo
      enddo
    case(3)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psi_out(i,j,k) = (psi(i,j,k)-(flux(i,j,k)-flux(i,j,k-1))*dzfi(k))/(1.-dt*dzfi(k)*(vel(i,j,k)-vel(i,j,k-1))) ! CHECK GRID
          enddo
        enddo
      enddo
    end select
  end subroutine mthinc_cmpt_vof_flux
#endif
  !
  subroutine mthinc_cmpt_d_param(idir,n,dt,dli,dzci,vel,psi,normx,normy,normz,lij,d)
    !
    ! computes the mthinc d paramater
    ! see Li et al., JCP 231 2328–2358
    !
    implicit none
    real(rp), parameter :: rrm = 0.5*(1-1._rp/sqrt(3._rp)), &
                           rrp = 0.5*(1+1._rp/sqrt(3._rp)), &
                           eps = epsilon(1._rp)
    integer , intent(in )                         :: idir
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in )                         :: dt
    real(rp), intent(in ), dimension(3)           :: dli
    real(rp), intent(in ), dimension(0:)          :: dzci
    real(rp), intent(in ), dimension(0:,0:,0:)    :: vel,psi,normx,normy,normz
    real(rp), intent(in ), dimension(0:,0:,0:,1:) :: lij ! curvature tensor
    real(rp), intent(out), dimension(0:,0:,0:)    :: d
    integer  :: i,j,k,ii,jj,kk,q
    integer  :: is_i,is_j,is_k
    real(rp) :: rvel,rsgn
    real(rp) :: dx_int,dy_int,dz_int,xl,xu,yl,yu,zl,zu
    real(rp) :: rnormx,rnormy,rnormz,rnorm_max
    integer  :: is_nx,is_ny,is_nz,cx,cy,cz
    real(rp) :: lxx,lxy,lxz,lyy,lyz,lzz,a001,a010,a100,a011,a101,a110,a002,a020,a200
    real(rp) :: a,b,aa,bb,bbm,bbp,cc,ccm,dd,ccp,qq,b0,b1,b2,b3,b4,c0,c1,c2,factor,z1
    !
    is_i = 0; if(idir == 1) is_i = 1
    is_j = 0; if(idir == 2) is_j = 1
    is_k = 0; if(idir == 3) is_k = 1
    do k=1-1,n(3)+1
      do j=1-1,n(2)+1
        do i=1-1,n(1)+1
          !
          ! check velocity sign
          !
          rvel = vel(i,j,k)
          rsgn = sign(1._rp,rvel)
          q = nint(0.5*(-rsgn+1.))
          ii = i + is_i*q
          jj = j + is_j*q
          kk = k + is_k*q
          !
          ! compute integration interval
          !
          dx_int = is_i*rvel*dt*dli(1)
          dy_int = is_j*rvel*dt*dli(2)
          dz_int = is_k*rvel*dt*dzci(kk) ! TODO: I should do something different for W?
          xl = q-dx_int
          xu = q
          yl = q-dy_int
          yu = q
          zl = q-dz_int
          zu = q
          !
          ! determine integration order
          !
          rnormx = abs(normx(ii,jj,kk))+eps; rnormy = abs(normy(ii,jj,kk))+eps; rnormz = abs(normz(ii,jj,kk))+eps
          rnorm_max = max(rnormx,rnormy,rnormz)
          is_nx = int(rnormx/rnorm_max+eps)
          is_ny = int(rnormy/rnorm_max+eps)
          is_nz = int(rnormz/rnorm_max+eps)
          cx = 1-is_nx
          cy = 1-is_ny
          cz = 1-is_nz
          rnormx = normx(ii,jj,kk); rnormy = normy(ii,jj,kk); rnormz = normz(ii,jj,kk)
          lxx  = lij(ii,jj,kk,1)
          lyy  = lij(ii,jj,kk,2)
          lzz  = lij(ii,jj,kk,3)
          lxy  = lij(ii,jj,kk,4)
          lxz  = lij(ii,jj,kk,5)
          lyz  = lij(ii,jj,kk,6)
          a200 = 0.5*cx*lxx
          a020 = 0.5*cy*lyy
          a002 = 0.5*cz*lzz
          a110 = cx*cy*lxy
          a101 = cx*cz*lxz
          a011 = cy*cz*lyz
          a100 = rnormx-(a200+0.5*(a110+a101))
          a010 = rnormy-(a020+0.5*(a110+a011))
          a001 = rnormz-(a002+0.5*(a101+a011))
          if( (psi(ii,jj,kk) <= limit).or.(psi(ii,jj,kk) >= 1.-limit) ) then
            dd = -huge(1._rp)
          else
            aa = is_nx*a100 + is_ny*a010 + is_nz*a001
            bbm = is_nx * a020*rrm**2+a002*rrm**2+a011*rrm*rrm+a010*rrm+a001*rrm + &
                  is_ny * a200*rrm**2+a002*rrm**2+a101*rrm*rrm+a100*rrm+a001*rrm + &
                  is_nz * a200*rrm**2+a020*rrm**2+a110*rrm*rrm+a100*rrm+a010*rrm
            bbp = is_nx * a020*rrp**2+a002*rrm**2+a011*rrp*rrm+a010*rrp+a001*rrm + &
                  is_ny * a200*rrp**2+a002*rrm**2+a101*rrp*rrm+a100*rrp+a001*rrm + &
                  is_nz * a200*rrp**2+a020*rrm**2+a110*rrp*rrm+a100*rrp+a010*rrm
            ccm = bbm
            ccp = is_nx * a020*rrm**2+a002*rrp**2+a011*rrm*rrp+a010*rrm+a001*rrp + &
                  is_ny * a200*rrm**2+a002*rrp**2+a101*rrm*rrp+a100*rrm+a001*rrp + &
                  is_nz * a200*rrm**2+a020*rrp**2+a110*rrm*rrp+a100*rrm+a010*rrp
            qq  = 2*aa*(2.*psi(i,j,k)-1.)
            aa  = exp(2*beta*aa)
            bbm = exp(2*beta*bbm)
            bbp = exp(2*beta*bbp)
            ccm = exp(2*beta*ccm)
            ccp = exp(2*beta*ccp)
            qq  = exp(2*beta*qq)
            !
            b4 = aa**2*bbm*bbp*ccm*ccp*(aa**2-qq)
            b3 = aa**2*(bbm*bbp*(ccm+ccp)+ccm*ccp*(bbm+bbp))*(aa-qq)/b4
            b2 = aa**2*((bbm+bbp)*(ccm+ccp)+bbm*bbp+ccm*ccp)*(1.-qq)/b4
            b1 = aa**((bbm+bbp)+(ccm+ccp))*(1.-aa*qq)/b4
            b0 = 1.-aa**2*qq/b4
            !
            ! solve quartic equation
            !
            c2 = -b2
            c1 = b1*b3-4._rp*b0
            c0 = b0*(4._rp*b2-b3**2)-b1**2
            !
            a  = -c2**2/9.+c1/3.
            b  = 2.*c2**3/27.-c1*c2/3.+c0
            factor = b**2+4.*a**3
            !
            if(factor >= 0._rp) then
              z1 = sign(1._rp,(-b+sqrt(factor)))*abs(0.5*(-b+sqrt(factor)))**(1._rp/3._rp) + &
                   sign(1._rp,(-b-sqrt(factor)))*abs(0.5*(-b-sqrt(factor)))**(1._rp/3._rp)
            else
              z1 = 2.*sqrt(-a)*cos(atan(sqrt(-factor)/(-b))/3.)
            endif
            z1 = z1 - c2/3.
            !
            aa = 0.5*b3
            bb = 0.5*z1
            dd = sqrt(bb**2-b0)
            cc = (-0.5_rp*b1+aa*bb)/dd
            !
            d(i,j,k) = 0.5*(-(aa-cc) + sqrt((aa-cc)**2-4.*(bb-dd)))
          end if
        end do
      end do
    end do
  end subroutine mthinc_cmpt_d_param
#if 0
  subroutine initvof(n,dli,vof)
    use mod_param, only: inivof,lx,ly,lz,cbcvof,xc,yc,zc,r
    use decomp_2d
    use mod_common_mpi, only: myid,ierr,coord
    !
    ! computes initial conditions for the velocity field
    !
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(3) :: dli
    real(8), dimension(0:,0:,0:), intent(out) :: vof
    integer :: i,j,k,q,ii,jj,kk
    real(8) :: x,y,z,xl,yl,zl,xx,yy,zz
    real(8) :: sdist,sdistmin
    real(8) :: zfilm_top,zfilm_bot,sdist1,sdist2,dfilm
    integer :: nbox
    real(8) :: eps,epsbox
    real(8), dimension(3) :: dl,dlbox
    integer, dimension(3) :: iperiod
    real(8) :: grid_vol_ratio
    !
    nbox = 100
    dl(:) = dli(:)**(-1)
    dlbox(:) = dl(:)/(1.d0*nbox)
    eps    = sqrt(dl(1)**2    + dl(2)**2    + dl(3)**2   )/2.d0
    epsbox = sqrt(dlbox(1)**2 + dlbox(2)**2 + dlbox(3)**2)/2.d0
    grid_vol_ratio = product(dlbox(:))/product(dl(:))
    iperiod(:) = 0
    do q=1,3
      if(cbcvof(0,q)//cbcvof(1,q).eq.'PP') iperiod(q) = 1
    enddo
    vof(:,:,:) = 0.d0
    select case(trim(inivof))
    case('bub')
      do k=1,n(3)
        z = (k-.5d0)*dl(3)
        do j=1,n(2)
          y = (j+coord(2)*n(2)-.5d0)*dl(2)
          do i=1,n(1)
            x = (i+coord(1)*n(1)-.5d0)*dl(1)
            sdistmin = max(lx,ly,lz)*2.d0
            do kk = -1,1
              do jj = -1,1
                do ii = -1,1
                  sdist = sqrt( (x+ii*iperiod(1)*lx-xc)**2 + &
                                (y+jj*iperiod(2)*ly-yc)**2 + &
                                (z+kk*iperiod(3)*lz-zc)**2 ) - r
                  if(abs(sdist).lt.sdistmin) sdistmin = sdist
                enddo
              enddo
            enddo
            sdist = sdistmin
            if(     sdist.lt.-eps ) then
              vof(i,j,k) = 1.d0
            elseif( sdist.gt. eps ) then
              vof(i,j,k) = 0.d0
            else
              zl = z-dl(3)/2.d0
              yl = y-dl(2)/2.d0
              xl = x-dl(1)/2.d0
              do kk=1,nbox
                zz = zl + kk*dlbox(3)
                do jj=1,nbox
                  yy = yl + jj*dlbox(2)
                  do ii=1,nbox
                    xx = xl + ii*dlbox(1)
                    sdist = sqrt((xx-xc)**2 + (yy-yc)**2 + (zz-zc)**2) - r
                    if(sdist.lt.-epsbox) vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
    case('flm')
      eps    = dl(3)
      epsbox = dlbox(3)
      grid_vol_ratio = product(dlbox(:))/product(dl(:))
      dfilm = 1./12.*.5
      zfilm_top = lz/2.+dfilm/2.
      zfilm_bot = lz/2.-dfilm/2.
      do k=1,n(3)
        z = (k-.5d0)*dl(3)
        do j=1,n(2)
          do i=1,n(1)
            sdist1 =  (z - zfilm_top)
            sdist2 = -(z - zfilm_bot)
            if(     all((/sdist1,sdist2/) .lt.-eps) ) then
              vof(i,j,k) = 1.d0
            elseif( all((/sdist1,sdist2/) .gt. eps) ) then
              vof(i,j,k) = 0.d0
            else
              zl = z-dl(3)/2.d0
              do kk=1,nbox
                zz = zl + kk*dlbox(3)
                do jj=1,nbox
                  do ii=1,nbox
                    sdist1 =  (zz - zfilm_top)
                    sdist2 = -(zz - zfilm_bot)
                    if( all((/sdist1,sdist2/).lt.-epsbox) ) vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
    case('uni')
      !$acc kernels
      do k=1,n3
        do j=1,n2
          do i=1,n1
            vof(i,j,k) = 1.0_rp
           enddo
         enddo
       enddo
       !$acc end kernels
       !
       !
    case('zer')
      !$acc kernels
      do k=1,n3
        do j=1,n2
          do i=1,n1
            vof(i,j,k) = 0.0_rp
          enddo
        enddo
      enddo
      !$acc end kernels
    case('twz')
      eps    = dl(3)
      epsbox = dlbox(3)
      grid_vol_ratio = product(dlbox(:))/product(dl(:))
      zfilm_bot = lz/2.0_rp
      !$acc kernels
      do k=1,n3
        kk = ijk_start(3) + k
        z = (kk-0.5_rp)*dl(3)
        do j=1,n2
          do i=1,n1
            sdist1 =  (z - zfilm_bot)
            if(     sdist1 .gt. +eps)  then
              vof(i,j,k) = 1.0_rp
            elseif( sdist1 .lt. -eps)  then
              vof(i,j,k) = 0.0_rp
            else
              zl = z-dl(3)/2.0_rp
              do kk=1,nbox
                zz = zl + kk*dlbox(3)
                do jj=1,nbox
                  do ii=1,nbox
                    sdist1 =  (zz - zfilm_bot)
                    if( sdist1.gt.+epsbox)  vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
      !
    case default
      if(myid.eq.0) print*, 'ERROR: invalid name for initial VoF field'
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation abortited due to errors in the case file ***'
      if(myid.eq.0) print*, '    check setup.h90'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    end select
  end subroutine initvof
#endif
end module mod_two_fluid_common
