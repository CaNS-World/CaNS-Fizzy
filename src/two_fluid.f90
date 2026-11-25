! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_two_fluid
  use mpi
  use mod_param, only: pi
  use mod_types
  implicit none
  private
  public init2fl,init2fl_sharp,clip_field,update_property,cmpt_norm_curv_youngs,bulk_mean_12_stag
  type sphere
    real(rp) :: xyz(3),r
    integer  :: imarker ! marker of the sphere (fifth, optional column of `spheres.in`; 1 if absent)
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
    integer, dimension(3) :: lo,hi
    integer :: i,j,k
    !
    prop1 = prop12(1)
    prop2 = prop12(2)
    lo(:) = lbound(p)
    hi(:) = ubound(p)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          p(i,j,k) = psi(i,j,k)*prop1+(1.-psi(i,j,k))*prop2
        end do
      end do
    end do
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
    integer :: n1,n2,n3,nh1,nh2,nh3,i,j,k
    !
    n1 = size(p,1)-2*nh(1)
    n2 = size(p,2)-2*nh(2)
    n3 = size(p,3)-2*nh(3)
    nh1 = nh(1)
    nh2 = nh(2)
    nh3 = nh(3)
    rmin = minmax(1); rmax = minmax(2)
    !
    !$acc parallel loop collapse(3) default(present) firstprivate(rmin,rmax) async(1)
    do concurrent(k=1-nh3:n3+nh3,j=1-nh2:n2+nh2,i=1-nh1:n1+nh1)
      p(i,j,k) = min(max(rmin,p(i,j,k)),rmax)
    end do
  end subroutine clip_field
  !
  subroutine bulk_mean_12_stag(n,nh_psi,nh_g,idir,grid_vol_ratio,wzc,psi,p12,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    use mpi
    use mod_types
    implicit none
    integer , intent(in), dimension(3) :: n
    integer , intent(in) :: nh_psi,nh_g
    integer , intent(in)               :: idir
    real(rp), intent(in), dimension(1-nh_g:) :: grid_vol_ratio,wzc
    real(rp), intent(in), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:) :: psi
    real(rp), intent(in), dimension(2)        :: p12
    real(rp), intent(out) :: mean
    real(rp) :: psi_loc
    integer  :: i,j,k
    integer  :: ierr
    mean = 0.
    select case(idir)
    case(1)
      !$acc data copy(mean) async(1)
      !$acc parallel loop collapse(3) default(present) private(psi_loc) reduction(+:mean) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psi_loc = 0.5*(psi(i+1,j,k)+psi(i,j,k))
            mean = mean + (psi_loc*p12(1)+(1.-psi_loc)*p12(2))*grid_vol_ratio(k)
          end do
        end do
      end do
      !$acc end data
    case(2)
      !$acc data copy(mean) async(1)
      !$acc parallel loop collapse(3) default(present) private(psi_loc) reduction(+:mean) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psi_loc = 0.5*(psi(i,j+1,k)+psi(i,j,k))
            mean = mean + (psi_loc*p12(1)+(1.-psi_loc)*p12(2))*grid_vol_ratio(k)
          end do
        end do
      end do
      !$acc end data
    case(3)
      !$acc data copy(mean) async(1)
      !$acc parallel loop collapse(3) default(present) private(psi_loc) reduction(+:mean) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            psi_loc = wzc(k)*psi(i,j,k)+(1.-wzc(k))*psi(i,j,k+1)
            mean = mean + (psi_loc*p12(1)+(1.-psi_loc)*p12(2))*grid_vol_ratio(k)
          end do
        end do
      end do
      !$acc end data
    case default 
      !$acc data copy(mean) async(1)
      !$acc parallel loop collapse(3) default(present) reduction(+:mean) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            mean = mean + (psi(i,j,k)*p12(1)+(1.-psi(i,j,k))*p12(2))*grid_vol_ratio(k)
          end do
        end do
      end do
      !$acc end data
    end select
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine bulk_mean_12_stag
  !
  subroutine cmpt_norm_curv_youngs(n,nh_psi,nh_ifc,nh_g,dli,dzci,dzfi,psi,normx,normy,normz,kappa,is_bound,theta_c)
    !
    ! computes the normals and curvature based on a volume-of-fluid
    ! or level-set field using finite-differences based on Youngs method
    !
    ! the normalized normals at the cell vertices are computed once and stored, and then averaged to the
    ! cell centers and differentiated for the curvature (rather than computing the eight vertex normals of
    ! each cell, i.e., each vertex normal eight times)
    !
    ! contact angles `theta_c` (in degrees, measured through phase 1) other than 90 degrees are imposed at
    ! the walls held by the task (`is_bound`) as in the wall adhesion of Brackbill et al., J. Comput. Phys. 100
    ! (1992) 335--354: the normal at the vertices on the wall, which points into phase 1 and is tangential
    ! to the wall for a field with zero normal gradient, is rotated to n = sin(theta)*t - cos(theta)*e,
    ! with t its unit tangential part and e the unit normal to the wall pointing into the domain; the
    ! normals and the curvature of the cells next to the wall then see the contact angle, while the field
    ! itself keeps its condition at the wall
    !
    implicit none
    real(rp), parameter :: eps = epsilon(1._rp)
    integer , intent(in ), dimension(3) :: n
    integer , intent(in) :: nh_psi,nh_ifc,nh_g
    real(rp), intent(in ), dimension(3) :: dli
    real(rp), intent(in ), dimension(1-nh_g:)          :: dzci,dzfi
    real(rp), intent(in ), dimension(1-nh_psi:,1-nh_psi:,1-nh_psi:)    :: psi
    real(rp), intent(out), dimension(1-nh_ifc:,1-nh_ifc:,1-nh_ifc:)    :: normx,normy,normz,kappa
    logical , intent(in ), dimension(0:1,3), optional :: is_bound
    real(rp), intent(in ), dimension(0:1,3), optional :: theta_c
    real(rp), allocatable, dimension(:,:,:), save :: vnx,vny,vnz ! normals at the vertices (i+1/2,j+1/2,k+1/2)
    real(rp) :: cs,sn,cn,tn
    integer  :: idir,ib,iw
    real(rp) :: psiccc,psicpc,psiccp,psicpp,psipcc,psippc,psipcp,psippp
    real(rp) :: mx_1,mx_2,mx_3,mx_4,mx_5,mx_6,mx_7,mx_8, &
                my_1,my_2,my_3,my_4,my_5,my_6,my_7,my_8, &
                mz_1,mz_2,mz_3,mz_4,mz_5,mz_6,mz_7,mz_8
    real(rp) :: normi,norm_x,norm_y,norm_z
    integer  :: i,j,k
    !
    if(.not.allocated(vnx)) then
      allocate(vnx(0:n(1),0:n(2),0:n(3)),vny(0:n(1),0:n(2),0:n(3)),vnz(0:n(1),0:n(2),0:n(3)))
      !$acc enter data create(vnx,vny,vnz) async(1)
    end if
    !
    ! normals at the vertices, from the eight cells around them
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(psiccc,psicpc,psiccp,psicpp,psipcc,psippc,psipcp,psippp,mx_1,my_1,mz_1,normi) async(1)
    do k=0,n(3)
      do j=0,n(2)
        do i=0,n(1)
          psiccc = psi(i  ,j  ,k  )
          psicpc = psi(i  ,j+1,k  )
          psiccp = psi(i  ,j  ,k+1)
          psicpp = psi(i  ,j+1,k+1)
          psipcc = psi(i+1,j  ,k  )
          psippc = psi(i+1,j+1,k  )
          psipcp = psi(i+1,j  ,k+1)
          psippp = psi(i+1,j+1,k+1)
          mx_1 = 0.25*((psipcc+psippc+psipcp+psippp)-(psiccc+psicpc+psiccp+psicpp))*dli(1)
          my_1 = 0.25*((psicpc+psippc+psicpp+psippp)-(psiccc+psipcc+psiccp+psipcp))*dli(2)
          mz_1 = 0.25*((psiccp+psipcp+psicpp+psippp)-(psiccc+psipcc+psicpc+psippc))*dzci(k)
          !
          ! the vectors are normalized with the reciprocal of their magnitude: away from the interface they vanish,
          ! and a division with a zero numerator is slow on GPUs
          !
          normi = 1._rp/(sqrt(mx_1**2+my_1**2+mz_1**2)+eps)
          vnx(i,j,k) = mx_1*normi
          vny(i,j,k) = my_1*normi
          vnz(i,j,k) = mz_1*normi
        end do
      end do
    end do
    !
    ! contact angles: normals at the vertices on the walls
    !
    if(present(theta_c)) then
      do idir=1,3
        do ib=0,1
          cs = cos(theta_c(ib,idir)*acos(-1._rp)/180._rp)
          if(.not.is_bound(ib,idir) .or. abs(cs) <= max(1.e-12_rp,1.e3_rp*epsilon(1._rp))) cycle
          sn = sin(theta_c(ib,idir)*acos(-1._rp)/180._rp)
          cn = -cs*merge(1._rp,-1._rp,ib == 0) ! component along the direction, for a unit tangential part
          iw = ib*n(idir)                      ! index of the vertices on the wall
          select case(idir)
          case(1)
            !$acc parallel loop collapse(2) default(present) private(tn) async(1)
            do k=0,n(3)
              do j=0,n(2)
                tn = sqrt(vny(iw,j,k)**2+vnz(iw,j,k)**2)
                vnx(iw,j,k) = cn*tn
                vny(iw,j,k) = sn*vny(iw,j,k)
                vnz(iw,j,k) = sn*vnz(iw,j,k)
              end do
            end do
          case(2)
            !$acc parallel loop collapse(2) default(present) private(tn) async(1)
            do k=0,n(3)
              do i=0,n(1)
                tn = sqrt(vnx(i,iw,k)**2+vnz(i,iw,k)**2)
                vnx(i,iw,k) = sn*vnx(i,iw,k)
                vny(i,iw,k) = cn*tn
                vnz(i,iw,k) = sn*vnz(i,iw,k)
              end do
            end do
          case(3)
            !$acc parallel loop collapse(2) default(present) private(tn) async(1)
            do j=0,n(2)
              do i=0,n(1)
                tn = sqrt(vnx(i,j,iw)**2+vny(i,j,iw)**2)
                vnx(i,j,iw) = sn*vnx(i,j,iw)
                vny(i,j,iw) = sn*vny(i,j,iw)
                vnz(i,j,iw) = cn*tn
              end do
            end do
          end select
        end do
      end do
    end if
    !
    ! normals and curvature at the cell centers
    !
    !$acc parallel loop collapse(3) default(present) &
    !$acc private(mx_1,my_1,mz_1,mx_2,my_2,mz_2,mx_3,my_3,mz_3,mx_4,my_4,mz_4,mx_5,my_5,mz_5,mx_6,my_6,mz_6,mx_7,my_7) &
    !$acc private(mz_7,mx_8,my_8,mz_8,norm_x,norm_y,norm_z,normi) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          mx_1 = vnx(i  ,j  ,k  ); my_1 = vny(i  ,j  ,k  ); mz_1 = vnz(i  ,j  ,k  )
          mx_2 = vnx(i  ,j-1,k  ); my_2 = vny(i  ,j-1,k  ); mz_2 = vnz(i  ,j-1,k  )
          mx_3 = vnx(i  ,j  ,k-1); my_3 = vny(i  ,j  ,k-1); mz_3 = vnz(i  ,j  ,k-1)
          mx_4 = vnx(i  ,j-1,k-1); my_4 = vny(i  ,j-1,k-1); mz_4 = vnz(i  ,j-1,k-1)
          mx_5 = vnx(i-1,j  ,k  ); my_5 = vny(i-1,j  ,k  ); mz_5 = vnz(i-1,j  ,k  )
          mx_6 = vnx(i-1,j-1,k  ); my_6 = vny(i-1,j-1,k  ); mz_6 = vnz(i-1,j-1,k  )
          mx_7 = vnx(i-1,j  ,k-1); my_7 = vny(i-1,j  ,k-1); mz_7 = vnz(i-1,j  ,k-1)
          mx_8 = vnx(i-1,j-1,k-1); my_8 = vny(i-1,j-1,k-1); mz_8 = vnz(i-1,j-1,k-1)
          !
          ! compute the normal vector
          !
          norm_x = .125*(mx_1+mx_2+mx_3+mx_4+mx_5+mx_6+mx_7+mx_8)
          norm_y = .125*(my_1+my_2+my_3+my_4+my_5+my_6+my_7+my_8)
          norm_z = .125*(mz_1+mz_2+mz_3+mz_4+mz_5+mz_6+mz_7+mz_8)
          normi = 1._rp/(sqrt(norm_x**2+norm_y**2+norm_z**2)+eps)
          normx(i,j,k) = norm_x*normi
          normy(i,j,k) = norm_y*normi
          normz(i,j,k) = norm_z*normi
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
  subroutine init2fl(inipsi,cbcpsi,seps,lo,hi,nh_psi,nh_g,l,dl,zc_g,psi,imarker)
    use mod_common_mpi, only: myid
    !
    ! computes initial conditions for the volume fraction field psi; with `imarker`, only the spheres of that
    ! marker are initialized (none, if the file holds no sphere of that marker)
    !
    implicit none
    character(len=*), intent(in) :: inipsi
    character(len=*), intent(in), dimension(0:1,3) :: cbcpsi
    real(rp), intent(in)                :: seps
    integer , intent(in), dimension(3)  :: lo,hi
    integer , intent(in) :: nh_psi,nh_g
    real(rp), intent(in), dimension(3)  :: l,dl
    real(rp), intent(in), dimension(1-nh_g:) :: zc_g
    real(rp), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:), intent(out) :: psi
    integer , intent(in), optional :: imarker
    integer  :: i,j,k,ii,jj,kk,q
    real(rp) :: x,y,z,xx,zz,xxc,yyc,zzc,r
    real(rp) :: sdist,sdistmin
    real(rp) :: zfilm_max,sdist1,sdist2
    integer , dimension(3) :: iperiod,idir
    logical , dimension(3) :: is_dim
    real(rp) :: psi_aux
    logical :: is_sphere,is_swap_pf
    type(sphere), allocatable, dimension(:) :: spheres
    integer :: nspheres
    integer :: ierr
    character(len=len(inipsi)) :: ini
    integer :: isuffix,ih,ia,ib
    real(rp), dimension(3) :: xyz
    !
    ! two-dimensional initial conditions take a suffix `-x`, `-y` or `-z` with the direction
    ! normal to their plane; without a suffix, they lie in the x-z plane (`-y`)
    !
    ini = inipsi
    ih = 2
    isuffix = len_trim(ini)
    if(isuffix > 2) then
      if(ini(isuffix-1:isuffix-1) == '-' .and. index('xyz',ini(isuffix:isuffix)) > 0) then
        ih = index('xyz',ini(isuffix:isuffix))
        ini = ini(:isuffix-2)
      end if
    end if
    ia = minval(pack([1,2,3],[1,2,3] /= ih)) ! in-plane directions
    ib = maxval(pack([1,2,3],[1,2,3] /= ih))
    !
    is_sphere  = .false.
    is_swap_pf = .false.
    psi(:,:,:) = 1._rp
    select case(trim(ini))
    case('uni')
    case('zer')
      psi(:,:,:) = 0._rp
    case('bub3')
      call read_sphere_file('spheres.in',spheres,nspheres,imarker)
      is_dim(:) = [.true. ,.true. ,.true.] ! sphere
      is_sphere = .true.
    case('bub2')
      call read_sphere_file('spheres.in',spheres,nspheres,imarker)
      is_dim(:) = .true.; is_dim(ih) = .false. ! cylinder along the direction normal to the plane
      is_sphere = .true.
    case('bub1')
      call read_sphere_file('spheres.in',spheres,nspheres,imarker)
      is_dim(:) = [.false.,.false.,.true.] ! planar film
      is_sphere = .true.
    case('drp3','dis3')
      call read_sphere_file('spheres.in',spheres,nspheres,imarker)
      is_dim(:) = [.true. ,.true. ,.true.] ! sphere
      is_sphere = .true.
      is_swap_pf = .true.
    case('drp2','dis2')
      call read_sphere_file('spheres.in',spheres,nspheres,imarker)
      is_dim(:) = .true.; is_dim(ih) = .false. ! cylinder along the direction normal to the plane
      is_sphere = .true.
      is_swap_pf = .true.
    case('drp1','dis1')
      call read_sphere_file('spheres.in',spheres,nspheres,imarker)
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
      ! interface at the mid height of the plane, with a sinusoidal bump
      ! with amplitude zfilm_max and wavelength equal to the width of the plane
      !
      zfilm_max = 0.01_rp*l(ia)
      do k=lo(3),hi(3)
        xyz(3) = zc_g(k)
        do j=lo(2),hi(2)
          xyz(2) = (j-0.5)*dl(2)
          do i=lo(1),hi(1)
            xyz(1) = (i-0.5)*dl(1)
            z = xyz(ib) - l(ib)/2
            x = xyz(ia)/l(ia)
            sdist = z - zfilm_max*cos(2*pi*x)
            psi(i,j,k) = smooth_step_tanh(sdist,seps)
          end do
        end do
      end do
    case('zalesak-disk')
      !
      ! Zalesak's disk
      !
      block
        real(rp) :: sw,sl,shift,xxc_slot,zzc_slot
        r   = 0.15_rp ! disk radius
        xxc = 0.5
        zzc = 0.75
        sw  = 0.05_rp/2 ! slot half width
        sl  = 0.25/2 ! slot half length
        shift = r/2
        xxc_slot = xxc
        zzc_slot = zzc-r+sl-shift
        sl = sl + shift
        !
        do k=lo(3),hi(3)
          xyz(3) = zc_g(k)
          do j=lo(2),hi(2)
            xyz(2) = (j-0.5)*dl(2)
            do i=lo(1),hi(1)
              xyz(1) = (i-0.5)*dl(1)
              z = xyz(ib)/l(ib)
              x = xyz(ia)/l(ia)
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
      if(nspheres == 0 .and. .not.present(imarker)) then
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
  subroutine init2fl_sharp(inipsi,cbcpsi,lo,hi,nh_psi,nh_g,l,dl,zc_g,dzf_g,psi,imarker)
    use mod_common_mpi, only: myid
    use mod_plic      , only: plic_volume
    !
    ! computes a sharp initial condition for the volume fraction field psi, i.e., the exact volume
    ! fraction of phase 1 in each cell, for the same initial geometries as in `init2fl`; the interface
    ! is described by a signed-distance(-like) function (positive in phase 1), and cells cut by the
    ! interface are refined recursively, with a tangent plane approximation of the interface at the
    ! finest level (cell size/2**maxlevel)
    !
    implicit none
    integer, parameter :: maxlevel = 5
    character(len=*), intent(in) :: inipsi
    character(len=*), intent(in), dimension(0:1,3) :: cbcpsi
    integer , intent(in), dimension(3)  :: lo,hi
    integer , intent(in) :: nh_psi,nh_g
    real(rp), intent(in), dimension(3)  :: l,dl
    real(rp), intent(in), dimension(1-nh_g:) :: zc_g,dzf_g
    real(rp), dimension(lo(1)-nh_psi:,lo(2)-nh_psi:,lo(3)-nh_psi:), intent(out) :: psi
    integer , intent(in), optional :: imarker
    type(sphere), allocatable, dimension(:) :: spheres
    integer , dimension(3) :: iperiod,idir
    logical , dimension(3) :: is_dim
    logical  :: is_sphere,is_swap_pf
    integer  :: nspheres,i,j,k,ierr
    character(len=len(inipsi)) :: ini
    integer :: isuffix,ih,ia,ib,ishape
    !
    ! two-dimensional initial conditions take a suffix `-x`, `-y` or `-z` with the direction
    ! normal to their plane; without a suffix, they lie in the x-z plane (`-y`)
    !
    ini = inipsi
    ih = 2
    isuffix = len_trim(ini)
    if(isuffix > 2) then
      if(ini(isuffix-1:isuffix-1) == '-' .and. index('xyz',ini(isuffix:isuffix)) > 0) then
        ih = index('xyz',ini(isuffix:isuffix))
        ini = ini(:isuffix-2)
      end if
    end if
    ia = minval(pack([1,2,3],[1,2,3] /= ih)) ! in-plane directions
    ib = maxval(pack([1,2,3],[1,2,3] /= ih))
    !
    is_sphere  = .false.
    is_swap_pf = .false.
    ishape = 0 ! spheres, cylinders and films from `spheres.in`
    select case(trim(ini))
    case('uni')
      ishape = 1
    case('zer')
      ishape = 2
    case('flm')
      ishape = 3
    case('cap-wav-1d')
      ishape = 4
    case('zalesak-disk')
      ishape = 5
    case('bub3','drp3','dis3')
      is_dim(:) = [.true. ,.true. ,.true.] ! sphere
      is_sphere = .true.
    case('bub2','drp2','dis2')
      is_dim(:) = .true.; is_dim(ih) = .false. ! cylinder along the direction normal to the plane
      is_sphere = .true.
    case('bub1','drp1','dis1')
      is_dim(:) = [.false.,.false.,.true.] ! planar film
      is_sphere = .true.
    case default
      if(myid == 0) print*, 'ERROR: invalid name for initial psi field'
      if(myid == 0) print*, ''
      if(myid == 0) print*, '*** Simulation aborted due to errors in the case file ***'
      if(myid == 0) print*, '    check INFO_INPUT.md'
      call MPI_FINALIZE(ierr)
      error stop
    end select
    if(is_sphere) then
      is_swap_pf = ini(1:3) == 'drp' .or. ini(1:3) == 'dis'
      call read_sphere_file('spheres.in',spheres,nspheres,imarker)
      if(nspheres == 0 .and. .not.present(imarker)) then
        if(myid == 0) print*, 'NOTE: `spheres.in` file not found.'
        if(myid == 0) print*, 'Initializing sphere/cylinder/plane in the domain center with radius r=0.25*minval(l(:)).'
        nspheres = 1
        allocate(spheres(nspheres))
        spheres(1)%xyz(:) = l(:)/2
        spheres(1)%r = 0.25*minval(l(:),mask=is_dim)
      end if
      iperiod(:) = 0
      where(cbcpsi(0,:)//cbcpsi(1,:) == 'PP') iperiod(:) = 1
      idir(:) = 1
      where(.not.is_dim(:)) idir(:) = 0
    end if
    !
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          psi(i,j,k) = frac_box((i-1)*dl(1),(j-1)*dl(2),zc_g(k)-0.5_rp*dzf_g(k),dl(1),dl(2),dzf_g(k),0)
        end do
      end do
    end do
  contains
    recursive real(rp) function frac_box(x0,y0,z0,hx,hy,hz,level) result(fr)
      !
      ! volume fraction of phase 1 in the box [x0,x0+hx]x[y0,y0+hy]x[z0,z0+hz]
      !
      real(rp), intent(in) :: x0,y0,z0,hx,hy,hz
      integer , intent(in) :: level
      real(rp) :: xc,yc,zc,d,diag,eg,gx,gy,gz,gn,nx,ny,nz
      integer  :: q1,q2,q3
      xc = x0+0.5_rp*hx
      yc = y0+0.5_rp*hy
      zc = z0+0.5_rp*hz
      d = sdf(xc,yc,zc)
      diag = 0.5_rp*sqrt(hx**2+hy**2+hz**2)
      if(d >= diag) then
        fr = 1._rp
      else if(d <= -diag) then
        fr = 0._rp
      else if(level >= maxlevel) then
        eg = 0.5_rp*min(hx,hy,hz) ! (differences across the box: a much smaller step is lost in single precision)
        gx = (sdf(xc+eg,yc,zc)-sdf(xc-eg,yc,zc))/(2._rp*eg)
        gy = (sdf(xc,yc+eg,zc)-sdf(xc,yc-eg,zc))/(2._rp*eg)
        gz = (sdf(xc,yc,zc+eg)-sdf(xc,yc,zc-eg))/(2._rp*eg)
        gn = sqrt(gx**2+gy**2+gz**2)
        if(gn > 0._rp) then
          nx = -gx/gn; ny = -gy/gn; nz = -gz/gn ! pointing out of phase 1
          fr = plic_volume(nx*hx,ny*hy,nz*hz,d/gn+0.5_rp*(nx*hx+ny*hy+nz*hz))
        else
          fr = merge(1._rp,0._rp,d > 0._rp)
        end if
      else
        fr = 0._rp
        do q3=0,1
          do q2=0,1
            do q1=0,1
              fr = fr+frac_box(x0+q1*0.5_rp*hx,y0+q2*0.5_rp*hy,z0+q3*0.5_rp*hz, &
                               0.5_rp*hx,0.5_rp*hy,0.5_rp*hz,level+1)
            end do
          end do
        end do
        fr = 0.125_rp*fr
      end if
    end function frac_box
    !
    real(rp) function sdf(x,y,z)
      !
      ! signed-distance(-like) function of the initial interface, positive in phase 1
      !
      real(rp), intent(in) :: x,y,z
      real(rp), dimension(3) :: xyz,dd
      real(rp) :: xx,zz,d1,d2,dmin,sw,sl,shift,xxc,zzc,r,eta,deta,zfilm_max
      integer  :: q,ii
      xyz(:) = [x,y,z]
      select case(ishape)
      case(1) ! 'uni'
        sdf = huge(1._rp)
      case(2) ! 'zer'
        sdf = -huge(1._rp)
      case(3) ! 'flm'
        sdf = z-l(3)/2
      case(4) ! 'cap-wav-1d'
        zfilm_max = 0.01_rp*l(ia)
        eta  =  zfilm_max*cos(2*pi*xyz(ia)/l(ia))
        deta = -zfilm_max*2*pi/l(ia)*sin(2*pi*xyz(ia)/l(ia))
        sdf  = (xyz(ib)-l(ib)/2-eta)/sqrt(1._rp+deta**2)
      case(5) ! 'zalesak-disk'
        r   = 0.15_rp ! (same constants as in `init2fl`)
        xxc = 0.5
        zzc = 0.75
        sw  = 0.05_rp/2
        sl  = 0.25/2
        shift = r/2
        xx = xyz(ia)/l(ia)
        zz = xyz(ib)/l(ib)
        d1 = sqrt((xx-xxc)**2+(zz-zzc)**2) - r
        xx = abs(xx-xxc)
        zz = abs(zz-(zzc-r+sl-shift))
        sl = sl + shift
        d2 = sqrt(max(xx-sw,0._rp)**2 + max(zz-sl,0._rp)**2) + min(max(xx-sw,zz-sl),0._rp)
        sdf = -max(d1,-d2)
      case default ! spheres, cylinders and films from `spheres.in`
        !
        ! the nearest periodic image is found along each direction separately
        !
        dmin = huge(1._rp)
        do q=1,nspheres
          dd(:) = huge(1._rp)
          do ii=-1,1
            dd(:) = min(dd(:),(xyz(:)+ii*iperiod(:)*l(:)-spheres(q)%xyz(:))**2*idir(:))
          end do
          dmin = min(dmin,sqrt(dd(1)+dd(2)+dd(3))-spheres(q)%r)
        end do
        sdf = merge(-dmin,dmin,is_swap_pf)
      end select
    end function sdf
  end subroutine init2fl_sharp
  !
  subroutine read_sphere_file(fname,spheres,nspheres,imarker)
    !
    ! reads the spheres of the file: centre, radius and, optionally, marker (1 if absent); with `imarker`,
    ! only the spheres of that marker are returned, and the file must exist and hold markers between 1 and
    ! `nmarker` only; blank lines are skipped
    !
    use mod_common_mpi, only: myid
    use mod_param     , only: nmarker
    implicit none
    character(len=*), intent(in) :: fname
    type(sphere), allocatable, intent(out), dimension(:) :: spheres
    integer, intent(out) :: nspheres
    integer, intent(in ), optional :: imarker
    type(sphere), allocatable, dimension(:) :: spheres_all
    character(len=1024) :: dummy
    integer :: q,iunit,ierr
    logical :: is_bubble_file,is_valid
    !
    is_valid = .true.
    inquire(file=fname,exist=is_bubble_file)
    if(.not.is_bubble_file) then
      nspheres = 0
      if(present(imarker)) then
        if(myid == 0) print*, 'ERROR: several markers (`nmarker > 1`) need the file `'//trim(fname)//'`.'
        is_valid = .false.
      end if
    else
      open(newunit=iunit,file=fname,action='read',iostat=ierr)
      if(ierr /= 0) then
        if(myid == 0) print*, 'ERROR: reading the file `'//trim(fname)//'`.'
        is_valid = .false.
      end if
      nspheres = 0
      do while(ierr == 0)
        read(iunit,'(A)',iostat=ierr) dummy
        if(ierr == 0 .and. len_trim(dummy) > 0) nspheres = nspheres + 1
      end do
      allocate(spheres(nspheres))
      if(is_valid) rewind(iunit)
      q = 0
      do while(q < nspheres .and. is_valid)
        read(iunit,'(A)') dummy
        if(len_trim(dummy) == 0) cycle
        q = q + 1
        read(dummy,*,iostat=ierr) spheres(q)%xyz(1),spheres(q)%xyz(2),spheres(q)%xyz(3),spheres(q)%r,spheres(q)%imarker
        if(ierr /= 0) then
          read(dummy,*,iostat=ierr) spheres(q)%xyz(1),spheres(q)%xyz(2),spheres(q)%xyz(3),spheres(q)%r
          spheres(q)%imarker = 1
        end if
        if(ierr /= 0) then
          if(myid == 0) print*, 'ERROR: reading the file `'//trim(fname)//'`: invalid line: '//trim(dummy)
          is_valid = .false.
        end if
      end do
      close(iunit)
      if(present(imarker) .and. is_valid) then
        if(any(spheres(:)%imarker < 1 .or. spheres(:)%imarker > nmarker)) then
          if(myid == 0) print*, 'ERROR: the markers of the spheres of `'//trim(fname)//'` (fifth column) must be &
                                        &between 1 and `nmarker`.'
          is_valid = .false.
        end if
        call move_alloc(spheres,spheres_all)
        nspheres = count(spheres_all(:)%imarker == imarker)
        allocate(spheres(nspheres))
        spheres(:) = pack(spheres_all(:),spheres_all(:)%imarker == imarker)
      end if
    end if
    if(.not.is_valid) then
      if(myid == 0) print*, ''
      if(myid == 0) print*, '*** Simulation aborted due to errors in the case file ***'
      if(myid == 0) print*, '    check INFO_INPUT.md'
      call MPI_FINALIZE(ierr)
      error stop
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
    ! smooth step function based on the hyperbolic tangent function
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
