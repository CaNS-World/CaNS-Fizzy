! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_two_fluid
  use mod_types
  implicit none
  private
  public initvof,cmpt_norm_curv
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
  subroutine cmpt_norm_curv(n,dl,dli,dzc,dzf,dzci,dzfi,psi,normx,normy,normz,kappa)
    !
    ! computes the normals and curvature based on a VoF field
    ! using finite-differences based on Youngs method
    !
    implicit none
    integer , parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(3)           :: dl,dli
    real(rp), intent(in ), dimension(0:)          :: dzc,dzf,dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:)    :: psi
    real(rp), intent(out), dimension(0:,0:,0:)    :: normx,normy,normz,kappa
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
#if 1
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
  subroutine initvof(inipsi,cbcpsi,lo,hi,l,dl,dzf,zc,psi)
    use mod_common_mpi, only: myid
    use mod_param, only: datadir
    !
    ! computes initial conditions for the velocity field
    !
    implicit none
    type sphere
      real(rp) :: x,y,z,r
    end type sphere
    character(len=*), intent(in) :: inipsi
    character(len=*), intent(in), dimension(0:1,3) :: cbcpsi
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), dimension(0:), intent(out) :: dzf,zc
    real(rp), dimension(0:,0:,0:), intent(out) :: psi
    integer  :: i,j,k,ii,jj,kk,q
    real(rp) :: x,y,z,xl,yl,zl,xx,yy,zz,xxc,yyc,zzc,r
    real(rp) :: sdist,sdistmin
    real(rp) :: dfilm,zfilm,zfilm_top,zfilm_bot,sdist1,sdist2
    integer  :: nbox
    real(rp), dimension(3) :: dlbox
    real(rp) :: eps,epsbox,grid_vol_ratio
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
      call read_sphere_file(trim(datadir)//'spheres.in',spheres,nspheres)
      is_dim(:) = [.true. ,.true. ,.true.] ! sphere
      is_sphere = .true.
    case('bu2')
      is_dim(:) = [.true. ,.false.,.true.] ! cylinder
      is_sphere = .true.
    case('bu1')
      is_dim(:) = [.false.,.false.,.true.] ! planar film
      is_sphere = .true.
    case('flm')
      nbox = 100
      grid_vol_ratio = 1./(1.*nbox**3)
      iperiod(:) = 0
      where(cbcpsi(0,:)//cbcpsi(1,:) == 'PP') iperiod(:) = 1
      do k=lo(3),hi(3)
        z = zc(k)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            sdist = z-l(3)/2.
            if(      sdist <= -eps ) then
              psi(i,j,k) = 1.
            else if( sdist <=  eps ) then
              psi(i,j,k) = 0.
            else
              dlbox(3) = dzf(k)/(1.*nbox)
              epsbox = dlbox(3)/2.
              zl = z-dzf(k)/2.
              yl = y-dl(2)/2.
              xl = x-dl(1)/2.
              do kk=1,nbox
                zz = zl + kk*dlbox(3)
                do jj=1,nbox
                  do ii=1,nbox
                    sdist = zz-l(3)/2.
                    if(sdist <= -epsbox) psi(i,j,k) = psi(i,j,k) + grid_vol_ratio
                  end do
                end do
              end do
            end if
          end do
        end do
      end do
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
      nbox = 100
      grid_vol_ratio = 1./(1.*nbox**3)
      iperiod(:) = 0
      where(cbcpsi(0,:)//cbcpsi(1,:) == 'PP') iperiod(:) = 1
      iexp(:) = 1
      where(.not.is_dim(:)) iexp(:) = 0
      do q = 1,nspheres
        xxc = spheres(q)%x
        yyc = spheres(q)%y
        zzc = spheres(q)%z
        r   = spheres(q)%r
        psi_aux = 0.
        do k=lo(3),hi(3)
          z = zc(k)
          do j=lo(2),hi(2)
            y = (j-0.5)*dl(2)
            do i=lo(1),hi(1)
              x = (i-0.5)*dl(1)
              sdistmin = maxval(l(:)*iexp(:))*2.0
              do kk = -1,1
                do jj = -1,1
                  do ii = -1,1
                    sdist = sqrt( (x+ii*iperiod(1)*l(1)-xxc)**(2*iexp(1)) + &
                                  (y+jj*iperiod(2)*l(2)-yyc)**(2*iexp(2)) + &
                                  (z+kk*iperiod(3)*l(3)-zzc)**(2*iexp(3)) ) - r
                    if(abs(sdist) <= sdistmin) sdistmin = sdist
                  end do
                end do
              end do
              sdist = sdistmin
              if(      sdist <= -eps ) then
                psi_aux = 1.
              else if( sdist <=  eps ) then
                psi_aux = 0.
              else
                dlbox(:) = [dl(1),dl(2),dzf(k)]/(1.*nbox)
                epsbox = sqrt(dlbox(1)**(2*iexp(1)) + dlbox(2)**(2*iexp(2)) + dlbox(3)**(2*iexp(3)))/2.
                zl = z-dzf(k)/2.
                yl = y-dl(2)/2.
                xl = x-dl(1)/2.
                do kk=1,nbox
                  zz = zl + kk*dlbox(3)
                  do jj=1,nbox
                    yy = yl + jj*dlbox(2)
                    do ii=1,nbox
                      xx = xl + ii*dlbox(1)
                      sdist = sqrt((xx-xxc)**(2*iexp(1)) + (yy-yyc)**(2*iexp(2)) + (zz-zzc)**(2*iexp(3))) - r
                      if(sdist <= -epsbox) psi_aux = psi_aux + grid_vol_ratio
                    end do
                  end do
                end do
              end if
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
    type sphere
      real(rp) :: x,y,z,r
    end type sphere 
    character(len=*), intent(in) :: fname
    type(sphere), allocatable, intent(out), dimension(:) :: spheres
    integer, intent(out) :: nspheres
    character(len=1024) :: dummy
    integer :: q,iunit,ierr
    !
    open(newunit=iunit,file=fname,action='read',iostat=ierr)
    if (ierr /= 0) then
      error stop 'Error reading input file '//trim(fname)//'.'
    end if
    nspheres = 0
    do while(ierr == 0)
      read(10, '(A)', iostat=ierr) dummy
      if(ierr > 0) then
        error stop 'Error reading input file '//trim(fname)//'.'
      end if
      nspheres = nspheres + 1
    end do
    allocate(spheres(nspheres))
    rewind(iunit)
    do q = 1,nspheres
      read(iunit,*) spheres(q)%x,spheres(q)%y,spheres(q)%z,spheres(q)%r
    end do
  end subroutine read_sphere_file
end module mod_two_fluid
