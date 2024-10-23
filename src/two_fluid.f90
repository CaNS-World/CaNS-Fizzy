! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_two_fluid
  use mpi
  use mod_param, only: pi,nh
  use mod_types
  implicit none
  private
  public init2fl,clip_field,cmpt_norm_curv_youngs
  type sphere
    real(rp) :: xyz(3),r
  end type sphere
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
  subroutine clip_field(nhalo,minmax,p)
    !
    ! clips a scalar field to certain maximum and minimum values
    !
    implicit none
    integer , intent(in   ), dimension(3) :: nhalo
    real(rp), intent(in   ), dimension(2) :: minmax
    real(rp), intent(inout), dimension(1-nhalo(1):,1-nhalo(2):,1-nhalo(3):) :: p
    real(rp) :: rmin,rmax
    integer :: n1,n2,n3,nh1,nh2,nh3,i,j,k
    !
    n1 = size(p,1)-2*nhalo(1)
    n2 = size(p,2)-2*nhalo(2)
    n3 = size(p,3)-2*nhalo(3)
    nh1 = nhalo(1)
    nh2 = nhalo(2)
    nh3 = nhalo(3)
    rmin = minmax(1); rmax = minmax(2)
    !
    !$acc parallel loop collapse(3) default(present) firstprivate(rmin,rmax) async(1)
    do concurrent(k=1-nh3:n3+nh3,j=1-nh2:n2+nh2,i=1-nh1:n1+nh1)
      p(i,j,k) = min(max(rmin,p(i,j,k)),rmax)
    end do
  end subroutine clip_field
  !
  subroutine cmpt_norm_curv_youngs(n,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
    !
    ! computes the normals and curvature based on a volume-of-fluid
    ! or level-set field using finite-differences based on Youngs method
    !
    implicit none
    real(rp), parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh:)          :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:)    :: psi
    real(rp), intent(out), dimension(1-nh:,1-nh:,1-nh:)    :: normx,normy,normz,kappa
    real(rp) :: psimmm,psimcm,psimpm,psicmm,psiccm,psicpm,psipmm,psipcm,psippm, &
                psimmc,psimcc,psimpc,psicmc,psiccc,psicpc,psipmc,psipcc,psippc, &
                psimmp,psimcp,psimpp,psicmp,psiccp,psicpp,psipmp,psipcp,psippp
    real(rp) :: mx_1,mx_2,mx_3,mx_4,mx_5,mx_6,mx_7,mx_8, &
                my_1,my_2,my_3,my_4,my_5,my_6,my_7,my_8, &
                mz_1,mz_2,mz_3,mz_4,mz_5,mz_6,mz_7,mz_8
    real(rp) :: norm,norm_x,norm_y,norm_z
    integer  :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          psimmm = psi(i-1,j-1,k-1)
          psimcm = psi(i-1,j  ,k-1)
          psimpm = psi(i-1,j+1,k-1)
          psicmm = psi(i  ,j-1,k-1)
          psiccm = psi(i  ,j  ,k-1)
          psicpm = psi(i  ,j+1,k-1)
          psipmm = psi(i+1,j-1,k-1)
          psipcm = psi(i+1,j  ,k-1)
          psippm = psi(i+1,j+1,k-1)
          psimmc = psi(i-1,j-1,k  )
          psimcc = psi(i-1,j  ,k  )
          psimpc = psi(i-1,j+1,k  )
          psicmc = psi(i  ,j-1,k  )
          psiccc = psi(i  ,j  ,k  )
          psicpc = psi(i  ,j+1,k  )
          psipmc = psi(i+1,j-1,k  )
          psipcc = psi(i+1,j  ,k  )
          psippc = psi(i+1,j+1,k  )
          psimmp = psi(i-1,j-1,k+1)
          psimcp = psi(i-1,j  ,k+1)
          psimpp = psi(i-1,j+1,k+1)
          psicmp = psi(i  ,j-1,k+1)
          psiccp = psi(i  ,j  ,k+1)
          psicpp = psi(i  ,j+1,k+1)
          psipmp = psi(i+1,j-1,k+1)
          psipcp = psi(i+1,j  ,k+1)
          psippp = psi(i+1,j+1,k+1)
          !
          mx_1 = 0.25*((psipcc+psippc+psipcp+psippp)-(psiccc+psicpc+psiccp+psicpp))*dli(1)
          mx_2 = 0.25*((psipcc+psipmc+psipcp+psipmp)-(psiccc+psicmc+psiccp+psicmp))*dli(1)
          mx_3 = 0.25*((psipcc+psippc+psipcm+psippm)-(psiccc+psicpc+psiccm+psicpm))*dli(1)
          mx_4 = 0.25*((psipcc+psipmc+psipcm+psipmm)-(psiccc+psicmc+psiccm+psicmm))*dli(1)
          mx_5 = 0.25*((psiccc+psicpc+psiccp+psicpp)-(psimcc+psimpc+psimcp+psimpp))*dli(1)
          mx_6 = 0.25*((psiccc+psicmc+psiccp+psicmp)-(psimcc+psimmc+psimcp+psimmp))*dli(1)
          mx_7 = 0.25*((psiccc+psicpc+psiccm+psicpm)-(psimcc+psimpc+psimcm+psimpm))*dli(1)
          mx_8 = 0.25*((psiccc+psicmc+psiccm+psicmm)-(psimcc+psimmc+psimcm+psimmm))*dli(1)
          !
          my_1 = 0.25*((psicpc+psippc+psicpp+psippp)-(psiccc+psipcc+psiccp+psipcp))*dli(2)
          my_2 = 0.25*((psiccc+psipcc+psiccp+psipcp)-(psicmc+psipmc+psicmp+psipmp))*dli(2)
          my_3 = 0.25*((psicpc+psippc+psicpm+psippm)-(psiccc+psipcc+psiccm+psipcm))*dli(2)
          my_4 = 0.25*((psiccc+psipcc+psiccm+psipcm)-(psicmc+psipmc+psicmm+psipmm))*dli(2)
          my_5 = 0.25*((psicpc+psimpc+psicpp+psimpp)-(psiccc+psimcc+psiccp+psimcp))*dli(2)
          my_6 = 0.25*((psiccc+psimcc+psiccp+psimcp)-(psicmc+psimmc+psicmp+psimmp))*dli(2)
          my_7 = 0.25*((psicpc+psimpc+psicpm+psimpm)-(psiccc+psimcc+psiccm+psimcm))*dli(2)
          my_8 = 0.25*((psiccc+psimcc+psiccm+psimcm)-(psicmc+psimmc+psicmm+psimmm))*dli(2)
          !
          mz_1 = 0.25*((psiccp+psipcp+psicpp+psippp)-(psiccc+psipcc+psicpc+psippc))*dzci(k  )
          mz_2 = 0.25*((psiccp+psipcp+psicmp+psipmp)-(psiccc+psipcc+psicmc+psipmc))*dzci(k  )
          mz_3 = 0.25*((psiccc+psipcc+psicpc+psippc)-(psiccm+psipcm+psicpm+psippm))*dzci(k-1)
          mz_4 = 0.25*((psiccc+psipcc+psicmc+psipmc)-(psiccm+psipcm+psicmm+psipmm))*dzci(k-1)
          mz_5 = 0.25*((psiccp+psimcp+psicpp+psimpp)-(psiccc+psimcc+psicpc+psimpc))*dzci(k  )
          mz_6 = 0.25*((psiccp+psimcp+psicmp+psimmp)-(psiccc+psimcc+psicmc+psimmc))*dzci(k  )
          mz_7 = 0.25*((psiccc+psimcc+psicpc+psimpc)-(psiccm+psimcm+psicpm+psimpm))*dzci(k-1)
          mz_8 = 0.25*((psiccc+psimcc+psicmc+psimmc)-(psiccm+psimcm+psicmm+psimmm))*dzci(k-1)
          !
          norm = sqrt(mx_1**2+my_1**2+mz_1**2)+eps
          mx_1 = mx_1/norm; my_1 = my_1/norm; mz_1 = mz_1/norm
          norm = sqrt(mx_2**2+my_2**2+mz_2**2)+eps
          mx_2 = mx_2/norm; my_2 = my_2/norm; mz_2 = mz_2/norm
          norm = sqrt(mx_3**2+my_3**2+mz_3**2)+eps
          mx_3 = mx_3/norm; my_3 = my_3/norm; mz_3 = mz_3/norm
          norm = sqrt(mx_4**2+my_4**2+mz_4**2)+eps
          mx_4 = mx_4/norm; my_4 = my_4/norm; mz_4 = mz_4/norm
          norm = sqrt(mx_5**2+my_5**2+mz_5**2)+eps
          mx_5 = mx_5/norm; my_5 = my_5/norm; mz_5 = mz_5/norm
          norm = sqrt(mx_6**2+my_6**2+mz_6**2)+eps
          mx_6 = mx_6/norm; my_6 = my_6/norm; mz_6 = mz_6/norm
          norm = sqrt(mx_7**2+my_7**2+mz_7**2)+eps
          mx_7 = mx_7/norm; my_7 = my_7/norm; mz_7 = mz_7/norm
          norm = sqrt(mx_8**2+my_8**2+mz_8**2)+eps
          mx_8 = mx_8/norm; my_8 = my_8/norm; mz_8 = mz_8/norm
          !
          ! compute the normal vector
          !
          norm_x = .125*(mx_1+mx_2+mx_3+mx_4+mx_5+mx_6+mx_7+mx_8)
          norm_y = .125*(my_1+my_2+my_3+my_4+my_5+my_6+my_7+my_8)
          norm_z = .125*(mz_1+mz_2+mz_3+mz_4+mz_5+mz_6+mz_7+mz_8)
          norm = sqrt(norm_x**2+norm_y**2+norm_z**2)+eps
          normx(i,j,k) = norm_x/norm
          normy(i,j,k) = norm_y/norm
          normz(i,j,k) = norm_z/norm
          !
          ! compute the curvature
          !
          kappa(i,j,k) = - ( 0.25*((mx_1+mx_2+mx_3+mx_4)-(mx_5+mx_6+mx_7+mx_8))*dli(1) + &
                             0.25*((my_1+my_3+my_5+my_7)-(my_2+my_4+my_6+my_8))*dli(2) + &
                             0.25*((mz_1+mz_2+mz_5+mz_6)-(mz_3+mz_4+mz_7+mz_8))*dzfi(k) )
        end do
      end do
    end do
  end subroutine cmpt_norm_curv_youngs
  !
  subroutine cmpt_norm_fd2(n,dli,dzfi,psi,normx,normy,normz)
    !
    ! computes the normals based on a volume-of-fluid or level-set field
    ! using second-order finite-differences
    !
    implicit none
    real(rp), parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh:)       :: dzfi
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:)    :: psi
    real(rp), intent(out), dimension(1-nh:,1-nh:,1-nh:)    :: normx,normy,normz
    real(rp) :: psixp,psixm,psiyp,psiym,psizp,psizm,dpsidx,dpsidy,dpsidz,norm
    integer  :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          psixp = 0.5*(psi(i+1,j,k)+psi(i  ,j,k))
          psixm = 0.5*(psi(i  ,j,k)+psi(i-1,j,k))
          psiyp = 0.5*(psi(i,j+1,k)+psi(i,j  ,k))
          psiym = 0.5*(psi(i,j  ,k)+psi(i,j-1,k))
          psizp = 0.5*(psi(i,j,k+1)+psi(i,j,k  ))
          psizm = 0.5*(psi(i,j,k  )+psi(i,j,k-1))
          dpsidx = (psixp-psixm)*dli(1)
          dpsidy = (psiyp-psiym)*dli(2)
          dpsidz = (psizp-psizm)*dzfi(k)
          norm   = sqrt(dpsidx**2+dpsidy**2+dpsidz**2)
          dpsidx = dpsidx/norm
          dpsidy = dpsidy/norm
          dpsidz = dpsidz/norm
        end do
      end do
    end do
  end subroutine cmpt_norm_fd2
  !
  subroutine cmpt_curv_fd2(n,dli,dzfi,psi,normx,normy,normz,kappa)
    !
    ! computes the curvature using second-order finite-differences
    !
    implicit none
    real(rp), parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh:)       :: dzfi
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:)    :: psi
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:)    :: normx,normy,normz
    real(rp), intent(out), dimension(1-nh:,1-nh:,1-nh:)    :: kappa
    real(rp) :: normxp,normxm,normyp,normym,normzp,normzm
    integer  :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          normxp = 0.5*(psi(i+1,j,k)+psi(i  ,j,k))
          normxm = 0.5*(psi(i  ,j,k)+psi(i-1,j,k))
          normyp = 0.5*(psi(i,j+1,k)+psi(i,j  ,k))
          normym = 0.5*(psi(i,j  ,k)+psi(i,j-1,k))
          normzp = 0.5*(psi(i,j,k+1)+psi(i,j,k  ))
          normzm = 0.5*(psi(i,j,k  )+psi(i,j,k-1))
          kappa(i,j,k) = - ( (normxp-normxm)*dli(1) + &
                             (normyp-normym)*dli(2) + &
                             (normzp-normzm)*dzfi(k) )
        end do
      end do
    end do
  end subroutine cmpt_curv_fd2
  !
  subroutine init2fl(inipsi,cbcpsi,seps,lo,hi,l,dl,zc_g,psi)
    use mod_common_mpi, only: myid
    !
    ! computes initial conditions for the volume fraction field psi
    !
    implicit none
    character(len=*), intent(in) :: inipsi
    character(len=*), intent(in), dimension(0:1,3) :: cbcpsi
    real(rp), intent(in)               :: seps
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), dimension(1-nh:), intent(in) :: zc_g
    real(rp), dimension(lo(1)-nh:,lo(2)-nh:,lo(3)-nh:), intent(out) :: psi
    integer  :: i,j,k,ii,jj,kk,q
    real(rp) :: x,y,z,xl,yl,zl,xx,yy,zz,xxc,yyc,zzc,r
    real(rp) :: sdist,sdistmin
    real(rp) :: dfilm,zfilm,zfilm_top,zfilm_bot,zfilm_max,sdist1,sdist2
    integer , dimension(3) :: iperiod,idir
    logical , dimension(3) :: is_dim
    real(rp) :: psi_aux
    logical :: is_sphere,is_swap_pf
    type(sphere), allocatable, dimension(:) :: spheres
    integer :: nspheres
    integer :: ierr
    !
    is_sphere  = .false.
    is_swap_pf = .false.
    psi(:,:,:) = 1._rp
    select case(trim(inipsi))
    case('uni')
    case('zer')
      psi(:,:,:) = 0._rp
    case('bub3')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.true. ,.true. ,.true.] ! sphere
      is_sphere = .true.
    case('bub2')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.true. ,.false.,.true. ] ! cylinder in xz plane
      is_sphere = .true.
      !is_dim(:) = [.true. ,.true. ,.false.] ! cylinder in xy plane
      !is_dim(:) = [.false.,.true. ,.true. ] ! cylinder in yz plane
    case('bub1')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.false.,.false.,.true.] ! planar film
      is_sphere = .true.
    case('drp3')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.true. ,.true. ,.true.] ! sphere
      is_sphere = .true.
      is_swap_pf = .true.
    case('drp2')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.true. ,.false.,.true. ] ! cylinder in xz plane
      !is_dim(:) = [.true. ,.true. ,.false.] ! cylinder in xy plane
      !is_dim(:) = [.false.,.true. ,.true. ] ! cylinder in yz plane
      is_sphere = .true.
      is_swap_pf = .true.
    case('drp1')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.false.,.false.,.true.] ! planar film
      is_sphere = .true.
      is_swap_pf = .true.
    case('flm')
      do k=lo(3),hi(3)
        z = zc_g(k)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            sdist = z-l(3)/2.
            psi(i,j,k) = smooth_step_tanh(sdist,seps)
          end do
        end do
      end do
    case('cap-wav-1d')
      !
      ! initial condition for the capillary wave benchmark
      ! (Prosperetti, Phys. Fluids 24 (1981) 1217)
      !
      ! interface at z = l(3)/2, with a sinusoidal bump
      ! with amplitude zfilm_max and wavelength = lx
      !
      zfilm_max = 0.01*l(1)
      do k=lo(3),hi(3)
        z = zc_g(k) - l(3)/2
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            x = (i-0.5)*dl(1)/l(1)
            sdist = z - zfilm_max*cos(2*pi*x)
            psi(i,j,k) = smooth_step_tanh(sdist,seps)
          end do
        end do
      end do
    case('zalesak-disk')
      !
      ! Zalesak's disk (WIP, not working yet)
      !
      block
        real(rp) :: sw,sl,shift,xxc_slot,zzc_slot
        r   = 0.15 ! disk radius
        xxc = 0.5
        zzc = 0.75
        sw  = 0.05/2 ! slot half width
        sl  = 0.25/2 ! slot half length
        shift = r/2
        xxc_slot = xxc
        zzc_slot = zzc-r+sl-shift
        sl = sl + shift
        !
        do k=lo(3),hi(3)
          z = zc_g(k)/l(3)
          do j=lo(2),hi(2)
            do i=lo(1),hi(1)
              x = (i-0.5)*dl(1)/l(1)
              !
              ! sdf of the disk
              !
              sdist1 = sqrt((x-xxc)**2+(z-zzc)**2) - r
              !
              ! sdf of the slot
              !
              xx = abs(x-xxc_slot)
              zz = abs(z-zzc_slot)
              sdist2 = sqrt(max(xx-sw,0._rp)**2 + max(zz-sl,0._rp)**2) + &
                       min(max(xx-sw,zz-sl),0._rp)
              !
              ! subtract sdfs
              !
              sdist  = -max(sdist1,-sdist2)
              !
              ! compute psi
              !
              psi(i,j,k) = smooth_step_tanh(sdist,seps)
            end do
          end do
        end do
      end block
    case default
      if(myid == 0) print*, 'ERROR: invalid name for initial psi field'
      if(myid == 0) print*, ''
      if(myid == 0) print*, '*** Simulation aborted due to errors in the case file ***'
      if(myid == 0) print*, '    check INFO_INPUT.md'
      call MPI_FINALIZE(ierr)
      error stop
    end select
    !
    if(is_sphere) then
      if(nspheres == 0) then
        if(myid == 0) print*, 'NOTE: `spheres.in` file not found.'
        if(myid == 0) print*, 'Initializing sphere/cylinder/plane in the domain center with radius r=0.25*minval(l(:)).'
        nspheres = nspheres + 1
        allocate(spheres(nspheres))
        q = nspheres
        spheres(q)%xyz(:) = l(:)/2
        spheres(q)%r = 0.25*minval(l(:),mask=is_dim)
      end if
      iperiod(:) = 0
      where(cbcpsi(0,:)//cbcpsi(1,:) == 'PP') iperiod(:) = 1
      idir(:) = 1
      where(.not.is_dim(:)) idir(:) = 0
      do q = 1,nspheres
        xxc = spheres(q)%xyz(1)
        yyc = spheres(q)%xyz(2)
        zzc = spheres(q)%xyz(3)
        r   = spheres(q)%r
        do k=lo(3),hi(3)
          z = zc_g(k)
          do j=lo(2),hi(2)
            y = (j-0.5)*dl(2)
            do i=lo(1),hi(1)
              x = (i-0.5)*dl(1)
              sdistmin = maxval(l(:)*idir(:))*2.0
              do kk = -1,1
                do jj = -1,1
                  do ii = -1,1
                    sdist = sqrt( (x+ii*iperiod(1)*l(1)-xxc)**2*idir(1) + &
                                  (y+jj*iperiod(2)*l(2)-yyc)**2*idir(2) + &
                                  (z+kk*iperiod(3)*l(3)-zzc)**2*idir(3) ) - r
                    if(abs(sdist) <= sdistmin) sdistmin = sdist
                  end do
                end do
              end do
              sdist = sdistmin
              psi_aux = smooth_step_tanh(sdist,seps)
              psi(i,j,k) = min(psi(i,j,k),psi_aux)
            end do
          end do
        end do
      end do
    end if
    if(is_swap_pf) psi(:,:,:) = 1.-psi(:,:,:)
  end subroutine init2fl
  !
  subroutine read_sphere_file(fname,spheres,nspheres)
    implicit none
    character(len=*), intent(in) :: fname
    type(sphere), allocatable, intent(out), dimension(:) :: spheres
    integer, intent(out) :: nspheres
    character(len=1024) :: dummy
    integer :: q,iunit,ierr
    logical :: is_bubble_file
    !
    inquire(file=fname,exist=is_bubble_file)
    if(.not.is_bubble_file) then
      nspheres = 0
    else
      open(newunit=iunit,file=fname,action='read',iostat=ierr)
      if (ierr /= 0) then
        error stop 'Error reading input file '//trim(fname)//'.'
      end if
      nspheres = -1
      do while(ierr == 0)
        read(iunit, '(A)', iostat=ierr) dummy
        if(ierr > 0) then
          error stop 'Error reading input file '//trim(fname)//'.'
        end if
        nspheres = nspheres + 1
      end do
      allocate(spheres(nspheres))
      rewind(iunit)
      do q = 1,nspheres
        read(iunit,*) spheres(q)%xyz(1),spheres(q)%xyz(2),spheres(q)%xyz(3),spheres(q)%r
      end do
      close(iunit)
    end if
  end subroutine read_sphere_file
  !
  pure elemental real(rp) function smooth_step_sin(r,eps) result(res)
    use mod_param, only:pi
    !$acc routine seq
    !
    ! smooth step function based on trigonometric functions
    !
    implicit none
    real(rp), intent(in) :: r,eps
    !
    if(r <= -eps) then
      res = 0.
    else if(r <= eps) then
      res = .5 + .5*r/eps + .5/pi*sin(pi*r/eps)
    else
      res = 1.
    end if
  end function smooth_step_sin
  !
  pure elemental real(rp) function smooth_step_erf(r,eps) result(res)
    !$acc routine seq
    !
    ! smooth step function based on the error function
    !
    implicit none
    !
    real(rp), intent(in) :: r,eps
    !
    res = .5*(1.+erf(r/eps))
  end function smooth_step_erf
  !
  pure elemental real(rp) function smooth_step_tanh(r,eps) result(res)
    !$acc routine seq
    !
    ! smooth step function based on the error function
    !
    implicit none
    !
    real(rp), intent(in) :: r,eps
    !
    res = .5*(1.+tanh(r/(2.*eps)))
  end function smooth_step_tanh
  !
  pure elemental real(rp) function smooth_sign(delta,phi) result(res)
    !$acc routine seq
    !
    ! smooth sign function
    !
    implicit none
    !
    real(rp), intent(in) :: delta,phi
    res = sign(1._rp,phi)
    if(abs(phi) <= delta) then
      res = phi/(sqrt(phi**2+delta**2))
    end if
  end function smooth_sign
  !
  pure elemental real(rp) function smooth_impulse(r,eps) result(res)
    use mod_param, only:pi
    !$acc routine seq
    !
    ! smooth impulse Dirac delta function using trigonometric functions
    !
    implicit none
    !
    real(rp), intent(in) :: r,eps
    if(abs(r) >= eps) then
      res = 0.
    else
      res = .5/eps + .5/eps*cos(pi*r/eps)
    end if
  end function smooth_impulse
#if 0
  subroutine cmpt_norm_curv_old(n,dli,dzci,dzfi,psi,normx,normy,normz,kappa)
    !
    ! computes the normals and curvature based on a VoF field
    ! using finite-differences based on Youngs method
    ! (replaced by `cmpt_norm_curv` above, which is slightly faster on GPUs)
    !
    implicit none
    real(rp), parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(3)           :: dli
    real(rp), intent(in ), dimension(1-nh:)       :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh:,1-nh:,1-nh:)    :: psi
    real(rp), intent(out), dimension(1-nh:,1-nh:,1-nh:)    :: normx,normy,normz,kappa
    real(rp) :: psimmm,psimcm,psimpm,psicmm,psiccm,psicpm,psipmm,psipcm,psippm, &
                psimmc,psimcc,psimpc,psicmc,psiccc,psicpc,psipmc,psipcc,psippc, &
                psimmp,psimcp,psimpp,psicmp,psiccp,psicpp,psipmp,psipcp,psippp
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
          psimmm = psi(i-1,j-1,k-1)
          psimcm = psi(i-1,j  ,k-1)
          psimpm = psi(i-1,j+1,k-1)
          psicmm = psi(i  ,j-1,k-1)
          psiccm = psi(i  ,j  ,k-1)
          psicpm = psi(i  ,j+1,k-1)
          psipmm = psi(i+1,j-1,k-1)
          psipcm = psi(i+1,j  ,k-1)
          psippm = psi(i+1,j+1,k-1)
          psimmc = psi(i-1,j-1,k  )
          psimcc = psi(i-1,j  ,k  )
          psimpc = psi(i-1,j+1,k  )
          psicmc = psi(i  ,j-1,k  )
          psiccc = psi(i  ,j  ,k  )
          psicpc = psi(i  ,j+1,k  )
          psipmc = psi(i+1,j-1,k  )
          psipcc = psi(i+1,j  ,k  )
          psippc = psi(i+1,j+1,k  )
          psimmp = psi(i-1,j-1,k+1)
          psimcp = psi(i-1,j  ,k+1)
          psimpp = psi(i-1,j+1,k+1)
          psicmp = psi(i  ,j-1,k+1)
          psiccp = psi(i  ,j  ,k+1)
          psicpp = psi(i  ,j+1,k+1)
          psipmp = psi(i+1,j-1,k+1)
          psipcp = psi(i+1,j  ,k+1)
          psippp = psi(i+1,j+1,k+1)
          !
          mx(1) = 0.25*((psipcc+psippc+psipcp+psippp)-(psiccc+psicpc+psiccp+psicpp))*dli(1)
          mx(2) = 0.25*((psipcc+psipmc+psipcp+psipmp)-(psiccc+psicmc+psiccp+psicmp))*dli(1)
          mx(3) = 0.25*((psipcc+psippc+psipcm+psippm)-(psiccc+psicpc+psiccm+psicpm))*dli(1)
          mx(4) = 0.25*((psipcc+psipmc+psipcm+psipmm)-(psiccc+psicmc+psiccm+psicmm))*dli(1)
          mx(5) = 0.25*((psiccc+psicpc+psiccp+psicpp)-(psimcc+psimpc+psimcp+psimpp))*dli(1)
          mx(6) = 0.25*((psiccc+psicmc+psiccp+psicmp)-(psimcc+psimmc+psimcp+psimmp))*dli(1)
          mx(7) = 0.25*((psiccc+psicpc+psiccm+psicpm)-(psimcc+psimpc+psimcm+psimpm))*dli(1)
          mx(8) = 0.25*((psiccc+psicmc+psiccm+psicmm)-(psimcc+psimmc+psimcm+psimmm))*dli(1)
          !
          my(1) = 0.25*((psicpc+psippc+psicpp+psippp)-(psiccc+psipcc+psiccp+psipcp))*dli(2)
          my(2) = 0.25*((psiccc+psipcc+psiccp+psipcp)-(psicmc+psipmc+psicmp+psipmp))*dli(2)
          my(3) = 0.25*((psicpc+psippc+psicpm+psippm)-(psiccc+psipcc+psiccm+psipcm))*dli(2)
          my(4) = 0.25*((psiccc+psipcc+psiccm+psipcm)-(psicmc+psipmc+psicmm+psipmm))*dli(2)
          my(5) = 0.25*((psicpc+psimpc+psicpp+psimpp)-(psiccc+psimcc+psiccp+psimcp))*dli(2)
          my(6) = 0.25*((psiccc+psimcc+psiccp+psimcp)-(psicmc+psimmc+psicmp+psimmp))*dli(2)
          my(7) = 0.25*((psicpc+psimpc+psicpm+psimpm)-(psiccc+psimcc+psiccm+psimcm))*dli(2)
          my(8) = 0.25*((psiccc+psimcc+psiccm+psimcm)-(psicmc+psimmc+psicmm+psimmm))*dli(2)
          !
          mz(1) = 0.25*((psiccp+psipcp+psicpp+psippp)-(psiccc+psipcc+psicpc+psippc))*dzci(k  )
          mz(2) = 0.25*((psiccp+psipcp+psicmp+psipmp)-(psiccc+psipcc+psicmc+psipmc))*dzci(k  )
          mz(3) = 0.25*((psiccc+psipcc+psicpc+psippc)-(psiccm+psipcm+psicpm+psippm))*dzci(k-1)
          mz(4) = 0.25*((psiccc+psipcc+psicmc+psipmc)-(psiccm+psipcm+psicmm+psipmm))*dzci(k-1)
          mz(5) = 0.25*((psiccp+psimcp+psicpp+psimpp)-(psiccc+psimcc+psicpc+psimpc))*dzci(k  )
          mz(6) = 0.25*((psiccp+psimcp+psicmp+psimmp)-(psiccc+psimcc+psicmc+psimmc))*dzci(k  )
          mz(7) = 0.25*((psiccc+psimcc+psicpc+psimpc)-(psiccm+psimcm+psicpm+psimpm))*dzci(k-1)
          mz(8) = 0.25*((psiccc+psimcc+psicmc+psimmc)-(psiccm+psimcm+psicmm+psimmm))*dzci(k-1)
          !
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
          norm = sqrt(normx(i,j,k)**2+normy(i,j,k)**2+normz(i,j,k)**2)+eps
          normx(i,j,k) = normx(i,j,k)/norm
          normy(i,j,k) = normy(i,j,k)/norm
          normz(i,j,k) = normz(i,j,k)/norm
          !
          ! compute the curvature
          !
          kappa(i,j,k) = - ( 0.25*((mx(1)+mx(2)+mx(3)+mx(4))-(mx(5)+mx(6)+mx(7)+mx(8)))*dli(1) + &
                             0.25*((my(1)+my(3)+my(5)+my(7))-(my(2)+my(4)+my(6)+my(8)))*dli(2) + &
                             0.25*((mz(1)+mz(2)+mz(5)+mz(6))-(mz(3)+mz(4)+mz(7)+mz(8)))*dzfi(k) )
        end do
      end do
    end do
  end subroutine cmpt_norm_curv_old
#endif
end module mod_two_fluid
