! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_two_fluid_common
  use mod_types
  implicit none
  private
  public update_property
  type, public :: vof_t
    ! vof type containing volume fraction, {x,y,z}-normals, and curvature fields
    real(rp), allocatable, dimension(:,:,:) :: vof,nx,ny,nz,kappa
  end type
  contains
  subroutine update_property(nh,prop12,vof,p)
    !
    ! updates a certain transport property based on a one-fluid formulation
    !
    implicit none
    integer , intent(in ), dimension(3)        :: nh
    real(rp), intent(in ), dimension(2)        :: prop12
    real(rp), intent(in ), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: vof
    real(rp), intent(out), dimension(1-nh(1):,1-nh(2):,1-nh(3):) :: p
    real(rp) :: prop1,prop2
    prop1 = prop12(1)
    prop2 = prop12(2)
    !$acc kernels default(present) async(1)
       p(:,:,:) = vof(:,:,:)*prop1+(1._rp-vof(:,:,:))*prop2
    !$acc end kernels
  end subroutine update_property
  !
  subroutine clip_field(nh,minmax,p)
    !
    ! clips a scalar field to a certain maximum and minimum values
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
  subroutine diff_geo_surface_youngs_mthinc(n,dl,dli,dzc,dzf,dzci,dzfi,vof,lij)
    !
    ! see Li et al., JCP 231 2328–2358
    !
    implicit none
    integer , parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(3)           :: dl,dli
    real(rp), intent(in ), dimension(0:)          :: dzc,dzf,dzci,dzfi
    type(vof_t   ), intent(inout) :: vof
    real(rp), intent(inout), dimension(0:,0:,0:,1:) :: lij ! curvature tensor
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
    associate(vof => vof%vof, normx => vof%nx, normy => vof%ny, normz => vof%nz, kappa => vof%kappa)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            !i+1/2 j+1/2 k+1/2
            mx(1) = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j+1,k+1)) - &
                     (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)))*dli(1)
            !i+1/2 j-1/2 k+1/2
            mx(2) = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j-1,k+1)) - &
                     (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)))*dli(1)
            !i+1/2 j+1/2 k-1/2
            mx(3) = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j+1,k-1)) - &
                     (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)))*dli(1)
            !i+1/2 j-1/2 k-1/2
            mx(4) = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j-1,k-1)) - &
                     (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)))*dli(1)
            !i-1/2 j+1/2 k+1/2
            mx(5) = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)) - &
                     (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j+1,k+1)))*dli(1)
            !i-1/2 j-1/2 k+1/2
            mx(6) = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)) - &
                     (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j-1,k+1)))*dli(1)
            !i-1/2 j+1/2 k-1/2
            mx(7) = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)) - &
                     (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j+1,k-1)))*dli(1)
            !i-1/2 j-1/2 k-1/2
            mx(8) = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)) - &
                     (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j-1,k-1)))*dli(1)
            !
            !i+1/2 j+1/2 k+1/2
            my(1) = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                     (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)))*dli(2)
            !i+1/2 j-1/2 k+1/2
            my(2) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)) - &
                     (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)))*dli(2)
            !i+1/2 j+1/2 k-1/2
            my(3) = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)) - &
                     (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)))*dli(2)
            !i+1/2 j-1/2 k-1/2
            my(4) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)) - &
                     (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(2)
            !i-1/2 j+1/2 k+1/2
            my(5) = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                     (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)))*dli(2)
            !i-1/2 j-1/2 k+1/2
            my(6) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)) - &
                     (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)))*dli(2)
            !i-1/2 j+1/2 k-1/2
            my(7) = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)) - &
                     (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)))*dli(2)
            !i-1/2 j-1/2 k-1/2
            my(8) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)) - &
                     (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(2)
            !
            !i+1/2 j+1/2 k+1/2
            mz(1) = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                     (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )))*dzci(k  )
            !i+1/2 j-1/2 k+1/2
            mz(2) = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)) - &
                     (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )))*dzci(k  )
            !i+1/2 j+1/2 k-1/2
            mz(3) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )) - &
                     (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)))*dzci(k  )
            !i+1/2 j-1/2 k-1/2
            mz(4) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )) - &
                     (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dzci(k  )
            !i-1/2 j+1/2 k+1/2
            mz(5) = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                     (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )))*dzci(k+1)
            !i-1/2 j-1/2 k+1/2
            mz(6) = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)) - &
                     (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )))*dzci(k+1)
            !i-1/2 j+1/2 k-1/2
            mz(7) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )) - &
                     (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)))*dzci(k  )
            !i-1/2 j-1/2 k-1/2
            mz(8) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )) - &
                     (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dzci(k  )
            !$acc loop seq
            do q=1,8
              norm = sqrt(mx(q)**2+my(q)**2+mz(q)**2+eps)**(-1)
              mx(q) = mx(q)/norm
              my(q) = my(q)/norm
              mz(q) = mz(q)/norm
            end do
            !
            ! compute the normal vector
            !
            normx(i,j,k) = (mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))/8
            normy(i,j,k) = (my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))/8
            normz(i,j,k) = (mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))/8
            !
            ! compute the curvature tensor
            !
            lij(i,j,k,1) = ((mx(1)+mx(2)+mx(3)+mx(4))-(mx(5)+mx(6)+mx(7)+mx(8)))*dl( 1)*0.25_rp
            lij(i,j,k,2) = ((my(1)+my(3)+my(5)+my(7))-(my(2)+my(4)+my(6)+my(8)))*dl( 2)*0.25_rp
            lij(i,j,k,3) = ((mz(1)+mz(2)+mz(5)+mz(6))-(mz(3)+mz(4)+mz(7)+mz(8)))*dzf(k)*0.25_rp
            kappa(i,j,k) = -(lij(i,j,k,1)*dli(1)**2+lij(i,j,k,2)*dli(2)**2+lij(i,j,k,3)*dzfi(k)**2) ! curvature
            !
            lij(i,j,k,4) = (((my(1)+my(2)+my(3)+my(4))-(my(5)+my(6)+my(7)+my(8)))*dl( 2) + &
                            ((mx(1)+mx(3)+mx(5)+mx(7))-(mx(2)+mx(4)+mx(6)+mx(8)))*dl( 1))*0.125_rp
            lij(i,j,k,5) = (((mz(1)+mz(2)+mz(3)+mz(4))-(mz(5)+mz(6)+mz(7)+mz(8)))*dzf(k) + &
                            ((mx(1)+mx(2)+mx(5)+mx(6))-(mx(3)+mx(4)+mx(7)+mx(8)))*dl( 1))*0.125_rp
            lij(i,j,k,6) = (((mz(1)+mz(3)+mz(5)+mz(7))-(mz(2)+mz(4)+mz(6)+mz(8)))*dzf(k) + &
                            ((my(1)+my(2)+my(5)+my(6))-(my(3)+my(4)+my(7)+my(8)))*dl( 2))*0.125_rp
          enddo
        enddo
      enddo
    end associate
  end subroutine diff_geo_surface_youngs_mthinc
end module mod_two_fluid_common
