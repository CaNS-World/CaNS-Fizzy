! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_two_fluid
  use mpi
  use mod_param, only: pi
  use mod_types
  implicit none
  private
  public initvof,cmpt_norm_curv
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
  subroutine cmpt_norm_curv(n,dl,dli,dzc,dzf,dzci,dzfi,psi,kappa)
    !
    ! computes the normals and curvature based on a VoF field
    ! using finite-differences based on Youngs method
    !
    implicit none
    real(rp), parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(3)           :: dl,dli
    real(rp), intent(in ), dimension(0:)          :: dzc,dzf,dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:)    :: psi
    !real(rp), intent(out), dimension(0:,0:,0:)    :: normx,normy,normz,kappa
    real(rp), intent(out), dimension(0:,0:,0:)    :: kappa
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
#if 0
          normx(i,j,k) = .125*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          normy(i,j,k) = .125*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          normz(i,j,k) = .125*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          norm = sqrt(normx(i,j,k)**2+normy(i,j,k)**2+normz(i,j,k)**2)+eps
          normx(i,j,k) = normx(i,j,k)/norm
          normy(i,j,k) = normy(i,j,k)/norm
          normz(i,j,k) = normz(i,j,k)/norm
#endif
          !
          ! compute the curvature
          !
          kappa(i,j,k) = - ( 0.25*((mx(1)+mx(2)+mx(3)+mx(4))-(mx(5)+mx(6)+mx(7)+mx(8)))*dli(1) + &
                             0.25*((my(1)+my(3)+my(5)+my(7))-(my(2)+my(4)+my(6)+my(8)))*dli(2) + &
                             0.25*((mz(1)+mz(2)+mz(5)+mz(6))-(mz(3)+mz(4)+mz(7)+mz(8)))*dzfi(k) )
        end do
      end do
    end do
  end subroutine cmpt_norm_curv
  !
  subroutine initvof(inipsi,cbcpsi,seps,lo,hi,l,dl,dzf_g,zc_g,psi)
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
    real(rp), dimension(0:), intent(out) :: dzf_g,zc_g
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(out) :: psi
    integer  :: i,j,k,ii,jj,kk,q
    real(rp) :: x,y,z,xl,yl,zl,xx,yy,zz,xxc,yyc,zzc,r
    real(rp) :: sdist,sdistmin
    real(rp) :: dfilm,zfilm,zfilm_top,zfilm_bot,zfilm_max,sdist1,sdist2
    integer , dimension(3) :: iperiod,iexp
    logical , dimension(3) :: is_dim
    real(rp) :: psi_aux
    logical :: is_sphere
    type(sphere), allocatable, dimension(:) :: spheres
    integer :: nspheres
    integer :: ierr
    !
    is_sphere = .false.
    psi(:,:,:) = 0.
    select case(trim(inipsi))
    case('uni')
      psi(:,:,:) = 1._rp
    case('zer')
    case('bu3')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.true. ,.true. ,.true.] ! sphere
      is_sphere = .true.
    case('bu2')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.true. ,.false.,.true.] ! cylinder
      is_sphere = .true.
    case('bu1')
      call read_sphere_file('spheres.in',spheres,nspheres)
      is_dim(:) = [.false.,.false.,.true.] ! planar film
      is_sphere = .true.
    case('flm')
      do k=lo(3),hi(3)
        z = zc_g(k)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            sdist = z-l(3)/2.
            psi(i,j,k) = smooth_step_tanh(sdistmin,seps)
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
      zfilm_max = 0.05
      do k=lo(3),hi(3)
        z = zc_g(k)/l(3) - 0.5
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            x = (i-0.5)*dl(1)/l(1)
            sdist = z - zfilm_max*cos(2*pi*x)
            psi(i,j,k) = smooth_step_tanh(sdistmin,seps)
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
              psi(i,j,k) = smooth_step_tanh(sdistmin,seps)
            end do
          end do
        end do
      end block
    case default
      if(myid == 0) print*, 'ERROR: invalid name for initial VoF field'
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
      iexp(:) = 1
      where(.not.is_dim(:)) iexp(:) = 0
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
              sdistmin = maxval(l(:)*iexp(:))*2.0
              do kk = -1,1
                do jj = -1,1
                  do ii = -1,1
                    sdist = sqrt( (x+ii*iperiod(1)*l(1)-xxc)**2*iexp(1) + &
                                  (y+jj*iperiod(2)*l(2)-yyc)**2*iexp(2) + &
                                  (z+kk*iperiod(3)*l(3)-zzc)**2*iexp(3) ) - r
                    if(abs(sdist) <= sdistmin) sdistmin = sdist
                  end do
                end do
              end do
              sdist = sdistmin
              psi_aux = smooth_step_tanh(sdist,seps)
              psi(i,j,k) = max(psi(i,j,k),psi_aux)
            end do
          end do
        end do
      end do
    end if
  end subroutine initvof
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
end module mod_two_fluid
