! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
program test_plic
  !
  ! unit tests of the PLIC geometric kernels (`src/plic.f90`):
  !   1. plane volume vs a quadruple-precision inclusion--exclusion reference
  !   2. inverse (plane constant from the volume) round trip
  !   3. additivity of sub-boxes and complement symmetry
  !   4. MYC, Youngs and centred-columns normals for exact planes (cubic, anisotropic, and stretched cells)
  !   5. normals for spheres and circles under grid refinement
  !   6. height-function curvature and normals for spheres and circles under grid refinement
  !   7. (GPU builds) device results vs host results
  !
  use mod_types
  use mod_plic
  implicit none
  integer , parameter :: qp = max(selected_real_kind(30),selected_real_kind(15)) ! quadruple precision, if available
  real(rp), parameter :: pi = acos(-1._rp)
  logical :: passed
  !
  passed = .true.
  call random_seed_fixed()
  call test_volume(passed)
  call test_alpha(passed)
  call test_additivity(passed)
  call test_normals_planes(passed)
  call test_normals_spheres(passed)
  call test_hf_spheres(passed)
#if defined(_OPENACC)
  call test_device(passed)
#endif
  if(passed) then
    print*, 'ALL PLIC UNIT TESTS PASSED'
  else
    print*, 'SOME PLIC UNIT TESTS FAILED'
    error stop 1
  end if
contains
  subroutine random_seed_fixed()
    integer :: n
    integer, allocatable :: seed(:)
    call random_seed(size=n)
    allocate(seed(n))
    seed(:) = 12345
    call random_seed(put=seed)
  end subroutine random_seed_fixed
  !
  real(rp) function urand(a,b)
    real(rp), intent(in) :: a,b
    real(rp) :: r
    call random_number(r)
    urand = a+(b-a)*r
  end function urand
  !
  subroutine random_normal(kind,mx,my,mz)
    !
    ! random normal of a given kind: 1 generic, 2 one tiny component, 3 two tiny components,
    ! 4 one zero component, 5 two zero components, 6 ties between components
    !
    integer , intent(in ) :: kind
    real(rp), intent(out) :: mx,my,mz
    real(rp) :: m(3),tiny
    integer :: q,p
    do q=1,3
      m(q) = urand(-1._rp,1._rp)
    end do
    call random_number(tiny)
    tiny = 10._rp**(-3-int(13*tiny)) ! 1e-3 ... 1e-15
    call random_number(mx); q = 1+int(3*mx)
    call random_number(mx); p = 1+mod(q+int(2*mx),3)
    select case(kind)
    case(2)
      m(q) = m(q)*tiny
    case(3)
      m(q) = m(q)*tiny
      if(p /= q) m(p) = m(p)*tiny
    case(4)
      m(q) = 0._rp
    case(5)
      m(q) = 0._rp
      if(p /= q) m(p) = 0._rp
    case(6)
      if(p /= q) m(p) = sign(abs(m(q)),m(p))
      call random_number(mx)
      if(mx < 0.3_rp) m(:) = sign(abs(m(1)),m(:))
      if(mx > 0.7_rp) then ! b1+b2 = b3
        m(1) = sign(0.25_rp,m(1)); m(2) = sign(0.25_rp,m(2)); m(3) = sign(0.5_rp,m(3))
      end if
    end select
    mx = m(1); my = m(2); mz = m(3)
  end subroutine random_normal
  !
  logical function is_well_conditioned(mx,my,mz)
    !
    ! the quadruple-precision inclusion--exclusion reference loses accuracy when the product of the
    ! nonzero normalized components is very small
    !
    real(rp), intent(in) :: mx,my,mz
    real(rp) :: b(3),s,p
    integer :: q
    b(:) = abs([mx,my,mz])
    s = sum(b)
    p = 1._rp
    do q=1,3
      if(b(q) > 0._rp) p = p*b(q)/s
    end do
    is_well_conditioned = p > 1.e-16_rp
  end function is_well_conditioned
  !
  real(rp) function vol_ref(mx,my,mz,alpha)
    !
    ! reference volume by inclusion--exclusion in quadruple precision
    !
    real(rp), intent(in) :: mx,my,mz,alpha
    real(qp) :: b(3),a,v,bb(3)
    integer :: q,nb,s1,s2,s3
    b(:) = abs(real([mx,my,mz],qp))
    a = real(alpha,qp)-min(real(mx,qp),0._qp)-min(real(my,qp),0._qp)-min(real(mz,qp),0._qp)
    nb = 0
    do q=1,3
      if(b(q) > 0._qp) then
        nb = nb+1
        bb(nb) = b(q)
      end if
    end do
    select case(nb)
    case(0)
      v = merge(1._qp,0._qp,a > 0._qp)
    case(1)
      v = min(max(a/bb(1),0._qp),1._qp)
    case(2)
      v = 0._qp
      do s1=0,1
        do s2=0,1
          v = v + (-1)**(s1+s2)*max(a-s1*bb(1)-s2*bb(2),0._qp)**2
        end do
      end do
      v = v/(2._qp*bb(1)*bb(2))
    case(3)
      v = 0._qp
      do s1=0,1
        do s2=0,1
          do s3=0,1
            v = v + (-1)**(s1+s2+s3)*max(a-s1*bb(1)-s2*bb(2)-s3*bb(3),0._qp)**3
          end do
        end do
      end do
      v = v/(6._qp*bb(1)*bb(2)*bb(3))
    end select
    vol_ref = real(min(max(v,0._qp),1._qp),rp)
  end function vol_ref
  !
  subroutine test_volume(passed)
    logical, intent(inout) :: passed
    integer , parameter :: nsamples = 200000
    real(rp), parameter :: tol = 1.e-14_rp
    integer , parameter :: ndump = 20000
    real(rp) :: mx,my,mz,alpha,s,alo,ahi,err,errmax(6),r,bp(6)
    integer :: n,kind,iunit
    errmax(:) = 0._rp
    !
    ! samples are also written to a file, checked in exact rational arithmetic by `check_exact.py`
    !
    open(newunit=iunit,file='volume_samples.txt',status='replace')
    do kind=1,6
      do n=1,nsamples
        call random_normal(kind,mx,my,mz)
        s = abs(mx)+abs(my)+abs(mz)
        alo = min(mx,0._rp)+min(my,0._rp)+min(mz,0._rp)
        ahi = max(mx,0._rp)+max(my,0._rp)+max(mz,0._rp)
        call random_number(r)
        if(r < 0.8_rp) then
          alpha = urand(alo-0.05_rp*s,ahi+0.05_rp*s)
        else ! close to the breakpoints of the piecewise volume function
          call random_number(r)
          bp(:) = [abs(mx),abs(my),abs(mz),abs(mx)+abs(my),abs(my)+abs(mz),abs(mx)+abs(mz)]
          alpha = alo + bp(1+int(6*r)) + urand(-1._rp,1._rp)*10._rp**(-int(urand(3._rp,16._rp)))*s
        end if
        if(n <= ndump) write(iunit,'(5es26.17e3)') mx,my,mz,alpha,plic_volume(mx,my,mz,alpha)
        if(precision(1._qp) > 2*precision(1._rp) .and. is_well_conditioned(mx,my,mz)) then
          err = abs(plic_volume(mx,my,mz,alpha)-vol_ref(mx,my,mz,alpha))
          errmax(kind) = max(errmax(kind),err)
        end if
      end do
    end do
    close(iunit)
    print '(a,6es10.2)', ' volume: max abs error per normal kind (generic,tiny1,tiny2,zero1,zero2,ties):', errmax(:)
    if(precision(1._qp) > 2*precision(1._rp)) then
      print '(a)', '         (vs quad-precision reference, where well conditioned; see also check_exact.py)'
    else
      print '(a)', '         (quad precision unavailable: reference check skipped; see check_exact.py)'
    end if
    if(any(errmax(:) > tol)) then
      print*, '  FAILED: volume error above ', tol
      passed = .false.
    end if
  end subroutine test_volume
  !
  subroutine test_alpha(passed)
    logical, intent(inout) :: passed
    integer , parameter :: nsamples = 200000
    real(rp), parameter :: tol = 1.e-14_rp
    real(rp), parameter, dimension(9) :: vspecial = [0._rp,1.e-300_rp,1.e-15_rp,1.e-8_rp,0.5_rp, &
                                                     1._rp-1.e-8_rp,1._rp-1.e-15_rp,1._rp-epsilon(1._rp),1._rp]
    real(rp) :: mx,my,mz,v,alpha,err,errmax(6),errend
    integer :: n,kind,q
    errmax(:) = 0._rp
    errend = 0._rp
    do kind=1,6
      do n=1,nsamples
        call random_normal(kind,mx,my,mz)
        if(n <= size(vspecial)) then
          v = vspecial(n)
        else
          call random_number(v)
        end if
        alpha = plic_alpha(mx,my,mz,v)
        err = abs(plic_volume(mx,my,mz,alpha)-v)
        errmax(kind) = max(errmax(kind),err)
      end do
    end do
    !
    ! empty and full cells map to the cube corners
    !
    do q=1,1000
      call random_normal(1,mx,my,mz)
      errend = max(errend,abs(plic_alpha(mx,my,mz,0._rp)-(min(mx,0._rp)+min(my,0._rp)+min(mz,0._rp))))
      errend = max(errend,abs(plic_alpha(mx,my,mz,1._rp)-(max(mx,0._rp)+max(my,0._rp)+max(mz,0._rp))))
    end do
    print '(a,6es10.2)', ' alpha : max |V(alpha(V))-V| per normal kind:                         ', errmax(:)
    print '(a,es10.2)' , ' alpha : max error of alpha for empty/full cells:', errend
    if(any(errmax(:) > tol) .or. errend > tol) then
      print*, '  FAILED: inverse error above ', tol
      passed = .false.
    end if
  end subroutine test_alpha
  !
  subroutine test_additivity(passed)
    logical, intent(inout) :: passed
    integer , parameter :: nsamples = 200000
    real(rp), parameter :: tol = 4.e-16_rp
    real(rp) :: mx,my,mz,alpha,c,v,va,vb,vc,err1,err2,r
    integer :: n,kind,idir
    err1 = 0._rp
    err2 = 0._rp
    do n=1,nsamples
      kind = 1+mod(n,6)
      call random_normal(kind,mx,my,mz)
      call random_number(v)
      alpha = plic_alpha(mx,my,mz,v)
      vc = plic_volume(mx,my,mz,alpha)
      call random_number(c)
      call random_number(r); idir = 1+int(3*r)
      select case(idir)
      case(1)
        va = plic_volume_box(mx,my,mz,alpha,0._rp,0._rp,0._rp,c,1._rp,1._rp)
        vb = plic_volume_box(mx,my,mz,alpha,c    ,0._rp,0._rp,1._rp-c,1._rp,1._rp)
      case(2)
        va = plic_volume_box(mx,my,mz,alpha,0._rp,0._rp,0._rp,1._rp,c,1._rp)
        vb = plic_volume_box(mx,my,mz,alpha,0._rp,c    ,0._rp,1._rp,1._rp-c,1._rp)
      case(3)
        va = plic_volume_box(mx,my,mz,alpha,0._rp,0._rp,0._rp,1._rp,1._rp,c)
        vb = plic_volume_box(mx,my,mz,alpha,0._rp,0._rp,c    ,1._rp,1._rp,1._rp-c)
      end select
      err1 = max(err1,abs(c*va+(1._rp-c)*vb-vc))
      err2 = max(err2,abs(plic_volume(-mx,-my,-mz,-alpha)-(1._rp-vc)))
      !
      ! the full box must reproduce the cell volume exactly (donor slab of unit width)
      !
      if(plic_volume_box(mx,my,mz,alpha,0._rp,0._rp,0._rp,1._rp,1._rp,1._rp) /= vc) err1 = huge(1._rp)
    end do
    print '(a,2es10.2)', ' boxes : max additivity error, max complement error:', err1,err2
    if(err1 > tol*4 .or. err2 > tol*4) then
      print*, '  FAILED: additivity/complement error above ', 4*tol
      passed = .false.
    end if
  end subroutine test_additivity
  !
  subroutine block_plane(nx,ny,nz,d,dx,dy,dz,c)
    !
    ! exact volume fractions of the region n.x <= d in a 3x3x3 block; x = 0 at the lower corner of the central cell
    !
    real(rp), intent(in ) :: nx,ny,nz,d,dx,dy
    real(rp), intent(in ), dimension(-1:1) :: dz
    real(rp), intent(out), dimension(-1:1,-1:1,-1:1) :: c
    real(rp) :: x0,y0,z0
    integer :: i,j,k
    do k=-1,1
      z0 = merge(-dz(-1),merge(0._rp,dz(0),k == 0),k == -1)
      do j=-1,1
        y0 = j*dy
        do i=-1,1
          x0 = i*dx
          c(i,j,k) = plic_volume(nx*dx,ny*dy,nz*dz(k),d-nx*x0-ny*y0-nz*z0)
        end do
      end do
    end do
  end subroutine block_plane
  !
  real(rp) function angle_err(mx,my,mz,dx,dy,dz0,nx,ny,nz)
    !
    ! angle between the cell-scaled estimate m (converted to physical space) and the exact normal n
    !
    real(rp), intent(in) :: mx,my,mz,dx,dy,dz0,nx,ny,nz
    real(rp) :: ex,ey,ez,s
    ex = mx/dx; ey = my/dy; ez = mz/dz0
    s = sqrt(ex**2+ey**2+ez**2)
    angle_err = acos(min(max((ex*nx+ey*ny+ez*nz)/s,-1._rp),1._rp))
  end function angle_err
  !
  subroutine stats(e,n,emean,e95,emax)
    real(rp), intent(inout), dimension(:) :: e
    integer , intent(in ) :: n
    real(rp), intent(out) :: emean,e95,emax
    real(rp) :: tmp
    integer :: i,j
    emean = sum(e(1:n))/max(n,1)
    emax  = maxval(e(1:n))
    ! simple selection of the 95th percentile (insertion sort on a copy is fine for these sizes)
    do i=2,n
      tmp = e(i)
      j = i-1
      do while(j >= 1)
        if(e(j) <= tmp) exit
        e(j+1) = e(j)
        j = j-1
      end do
      e(j+1) = tmp
    end do
    e95 = e(max(1,int(0.95*n)))
  end subroutine stats
  !
  subroutine test_normals_planes(passed)
    logical, intent(inout) :: passed
    integer , parameter :: nsamples = 20000
    real(rp) :: c(-1:1,-1:1,-1:1),dz(-1:1),dx,dy,nx,ny,nz,d,s,mx,my,mz,th,ph
    real(rp), allocatable :: e(:,:)
    real(rp) :: emean(3),e95(3),emax(3)
    character(len=*), parameter :: grids(3) = ['cubic      ','anisotropic','stretched  ']
    integer :: igrid,n,nn,q
    allocate(e(nsamples,3))
    do igrid=1,3
      select case(igrid)
      case(1)
        dx = 1._rp; dy = 1._rp; dz(:) = 1._rp
      case(2)
        dx = 1._rp; dy = 0.6_rp; dz(:) = 1.7_rp
      case(3)
        dx = 1._rp; dy = 1._rp; dz(:) = [0.8_rp,1._rp,1.25_rp]
      end select
      nn = 0
      do n=1,nsamples*4
        call random_number(th); th = acos(2._rp*th-1._rp)
        call random_number(ph); ph = 2._rp*pi*ph
        nx = sin(th)*cos(ph); ny = sin(th)*sin(ph); nz = cos(th)
        d = nx*urand(0._rp,dx)+ny*urand(0._rp,dy)+nz*urand(0._rp,dz(0)) ! plane through the central cell
        call block_plane(nx,ny,nz,d,dx,dy,dz,c)
        if(c(0,0,0) < 1.e-6_rp .or. c(0,0,0) > 1._rp-1.e-6_rp) cycle
        nn = nn+1
        call plic_normal_myc(c,dx,dy,dz,[.true.,.true.,.true.],mx,my,mz)
        e(nn,1) = angle_err(mx,my,mz,dx,dy,dz(0),nx,ny,nz)
        call plic_normal_youngs(c,dx,dy,dz,mx,my,mz)
        e(nn,2) = angle_err(mx,my,mz,dx,dy,dz(0),nx,ny,nz)
        call plic_normal_cc(c,dx,dy,dz,[.true.,.true.,.true.],mx,my,mz)
        e(nn,3) = angle_err(mx,my,mz,dx,dy,dz(0),nx,ny,nz)
        if(nn == nsamples) exit
      end do
      do q=1,3
        call stats(e(:,q),nn,emean(q),e95(q),emax(q))
      end do
      emean = emean*180/pi; e95 = e95*180/pi; emax = emax*180/pi
      print '(a,a11,a,i6,a)', ' planes (',grids(igrid),'), ',nn,' samples; angle error [deg] mean / 95% / max:'
      print '(a,3f9.4)', '    MYC   ', emean(1),e95(1),emax(1)
      print '(a,3f9.4)', '    Youngs', emean(2),e95(2),emax(2)
      print '(a,3f9.4)', '    CC    ', emean(3),e95(3),emax(3)
      if(emax(1) > 5._rp .or. emean(1) > min(emean(2),emean(3))) then
        print*, '  FAILED: MYC error above 5 degrees, or MYC worse than Youngs or CC on average'
        passed = .false.
      end if
    end do
  end subroutine test_normals_planes
  !
  real(rp) function frac_sphere(xc,yc,zc,r,x0,y0,z0,dx,dy,dz,is_2d,nref)
    !
    ! volume fraction of the box [x0,x0+dx]x[y0,y0+dy]x[z0,z0+dz] inside the sphere (or cylinder along y,
    ! if is_2d) of radius r centred at (xc,yc,zc): sub-boxes cut by the surface use the tangent plane;
    ! nref returns the (unnormalized) average of the surface normal over the sub-boxes cut by the surface
    !
    real(rp), intent(in) :: xc,yc,zc,r,x0,y0,z0,dx,dy,dz
    logical , intent(in) :: is_2d
    real(rp), intent(out) :: nref(3)
    integer , parameter :: ns = 16
    real(rp) :: hx,hy,hz,px,py,pz,rr,dist,diag,nx,ny,nz,f,ax,ay,az
    integer :: i,j,k
    hx = dx/ns; hy = dy/ns; hz = dz/ns
    f = 0._rp
    nref(:) = 0._rp
    do k=1,ns
      pz = z0+(k-0.5_rp)*hz
      do j=1,ns
        py = y0+(j-0.5_rp)*hy
        do i=1,ns
          px = x0+(i-0.5_rp)*hx
          ax = px-xc; ay = merge(0._rp,py-yc,is_2d); az = pz-zc
          rr = sqrt(ax**2+ay**2+az**2)
          dist = rr-r
          diag = 0.5_rp*sqrt(hx**2+merge(0._rp,hy**2,is_2d)+hz**2)
          if(dist <= -diag) then
            f = f+1._rp
          else if(dist < diag) then
            nx = ax/rr; ny = ay/rr; nz = az/rr ! outward normal, i.e., pointing out of phase 1 (inside)
            f = f+plic_volume(nx*hx,ny*hy,nz*hz,-dist+0.5_rp*(nx*hx+ny*hy+nz*hz))
            nref(:) = nref(:)+[nx,ny,nz]
          end if
        end do
      end do
    end do
    frac_sphere = f/ns**3
  end function frac_sphere
  !
  subroutine test_normals_spheres(passed)
    logical, intent(inout) :: passed
    integer , parameter :: nres = 4
    integer , parameter, dimension(nres) :: rres = [4,8,16,32]
    real(rp) :: r,xc,yc,zc,dzu(-1:1),c(-1:1,-1:1,-1:1),mx,my,mz,nx,ny,nz,rr,h
    real(rp) :: emean(3,nres,2),emax(3,nres,2)
    real(rp), allocatable :: f(:,:,:),fn(:,:,:,:)
    integer :: ires,ng,i,j,k,nn,idim,q
    logical :: is_2d,act(3)
    do idim=1,2
      is_2d = idim == 2
      act(:) = [.true.,.not.is_2d,.true.]
      do ires=1,nres
        h  = 1._rp
        r  = rres(ires)*h
        ng = 2*rres(ires)+8
        xc = 0.5_rp*ng*h+0.123_rp*h; yc = 0.5_rp*ng*h+0.311_rp*h; zc = 0.5_rp*ng*h-0.207_rp*h
        allocate(f(0:ng+1,0:ng+1,0:ng+1),fn(3,0:ng+1,0:ng+1,0:ng+1))
        do k=0,ng+1
          do j=0,ng+1
            if(is_2d .and. j > 1) then
              f(:,j,k) = f(:,1,k)
              fn(:,:,j,k) = fn(:,:,1,k)
              cycle
            end if
            do i=0,ng+1
              f(i,j,k) = frac_sphere(xc,yc,zc,r,(i-1)*h,(j-1)*h,(k-1)*h,h,h,h,is_2d,fn(:,i,j,k))
            end do
          end do
        end do
        if(is_2d) f(:,0,:) = f(:,1,:)
        if(is_2d) fn(:,:,0,:) = fn(:,:,1,:)
        dzu(:) = h
        emean(:,ires,idim) = 0._rp
        emax(:,ires,idim) = 0._rp
        nn = 0
        do k=1,ng
          do j=1,merge(1,ng,is_2d)
            do i=1,ng
              if(f(i,j,k) < 1.e-6_rp .or. f(i,j,k) > 1._rp-1.e-6_rp) cycle
              c(:,:,:) = f(i-1:i+1,j-1:j+1,k-1:k+1)
              nx = fn(1,i,j,k); ny = fn(2,i,j,k); nz = fn(3,i,j,k)
              rr = sqrt(nx**2+ny**2+nz**2); nx = nx/rr; ny = ny/rr; nz = nz/rr
              nn = nn+1
              do q=1,3
                select case(q)
                case(1)
                  call plic_normal_myc(c,h,h,dzu,act,mx,my,mz)
                case(2)
                  call plic_normal_youngs(c,h,h,dzu,mx,my,mz)
                case(3)
                  call plic_normal_cc(c,h,h,dzu,act,mx,my,mz)
                end select
                rr = angle_err(mx,my,mz,h,h,h,nx,ny,nz)
                emean(q,ires,idim) = emean(q,ires,idim)+rr
                emax(q,ires,idim) = max(emax(q,ires,idim),rr)
              end do
            end do
          end do
        end do
        emean(:,ires,idim) = emean(:,ires,idim)/nn
        deallocate(f,fn)
      end do
      print '(a,a)', ' ',merge('circles (2D, y-invariant)','spheres                  ',is_2d)
      print '(a)', '    R/h   mean angle error [rad]: MYC / Youngs / CC        max: MYC / Youngs / CC'
      do ires=1,nres
        print '(i7,3x,3es10.2,5x,3es10.2)', rres(ires),emean(:,ires,idim),emax(:,ires,idim)
      end do
      print '(a,3f6.2)', '    observed order of the mean error (last two resolutions): ', &
            log(emean(:,nres-1,idim)/emean(:,nres,idim))/log(2._rp)
      !
      ! MYC is not exact for planes, so its error tends to a small plateau (about 2e-3 rad, see the plane tests)
      ! rather than converging at a fixed order; require a monotone decrease and a smaller error than Youngs and CC
      !
      if(any(emean(1,2:nres,idim) >= emean(1,1:nres-1,idim)) .or. &
         any(emean(1,2:nres,idim) >= min(emean(2,2:nres,idim),emean(3,2:nres,idim)))) then
        print*, '  FAILED: MYC error does not decrease with resolution, or is not the smallest'
        passed = .false.
      end if
    end do
  end subroutine test_normals_spheres
  subroutine hf_cell(nq,f,i,j,k,h,act,kappa,nrm,is_valid)
    !
    ! height-function curvature and normal of cell (i,j,k) on a uniform grid of spacing h, trying the
    ! directions in decreasing order of the MYC normal components (as in the code)
    !
    integer , intent(in ) :: nq
    real(rp), intent(in ), dimension(0:,0:,0:) :: f
    integer , intent(in ) :: i,j,k
    real(rp), intent(in ) :: h
    logical , intent(in ) :: act(3)
    real(rp), intent(out) :: kappa,nrm(3)
    logical , intent(out) :: is_valid
    real(rp) :: c(-1:1,-1:1,-nq:nq),c3(-1:1,-1:1,-1:1),dq(-nq:nq),hgt(-1:1,-1:1),sgn,m(3),na,nb,nn3,dz3(-1:1)
    integer :: idir(3),d,a,b,q,o(3),ea(3),eb(3),ed(3),r,ia,ib
    dz3(:) = h
    c3(:,:,:) = f(i-1:i+1,j-1:j+1,k-1:k+1)
    call plic_normal_myc(c3,h,h,dz3,act,m(1),m(2),m(3))
    idir = [1,2,3]
    ! order by decreasing |m|
    do r=1,2
      do q=1,3-r
        if(abs(m(idir(q))) < abs(m(idir(q+1)))) idir([q,q+1]) = idir([q+1,q])
      end do
    end do
    is_valid = .false.
    kappa = 0._rp
    nrm(:) = 0._rp
    dq(:) = h
    do r=1,3
      d = idir(r)
      if(.not.act(d)) cycle
      ia = merge(2,1,d == 1); ib = merge(2,3,d == 3)
      ed = 0; ed(d) = 1; ea = 0; ea(ia) = 1; eb = 0; eb(ib) = 1
      do q=-nq,nq
        do b=-1,1
          do a=-1,1
            o = [i,j,k]+q*ed+a*ea+b*eb
            c(a,b,q) = f(o(1),o(2),o(3))
          end do
        end do
      end do
      call hf_heights(nq,c,dq,1.e-6_rp,hgt,sgn,is_valid)
      if(is_valid) then
        kappa = hf_curvature(hgt,sgn,h,h,h,h,act(ia),act(ib))
        call hf_normal(hgt,sgn,h,h,h,h,act(ia),act(ib),na,nb,nn3)
        nrm(ia) = na; nrm(ib) = nb; nrm(d) = nn3
        exit
      end if
    end do
  end subroutine hf_cell
  !
  subroutine test_hf_spheres(passed)
    !
    ! height-function curvature and normals for spheres and circles (exact volume fractions)
    !
    logical, intent(inout) :: passed
    integer , parameter :: nres = 4
    integer , parameter, dimension(nres) :: rres = [4,8,16,32]
    real(rp) :: r,xc,yc,zc,h,kex,kappa,nrm(3),nr(3),rr,e2,einf,frac,ean,eanmax
    real(rp) :: res(5,nres)
    real(rp), allocatable :: f(:,:,:),fn(:,:,:,:)
    integer :: ires,ng,i,j,k,nn,nv,idim,nq
    logical :: is_2d,act(3),is_valid
    do nq=3,4
    do idim=1,2
      is_2d = idim == 2
      act(:) = [.true.,.not.is_2d,.true.]
      do ires=1,nres
        h  = 1._rp
        r  = rres(ires)*h
        ng = 2*rres(ires)+12
        xc = 0.5_rp*ng*h+0.123_rp*h; yc = 0.5_rp*ng*h+0.311_rp*h; zc = 0.5_rp*ng*h-0.207_rp*h
        kex = merge(1._rp/r,2._rp/r,is_2d)
        allocate(f(0:ng+1,0:ng+1,0:ng+1),fn(3,0:ng+1,0:ng+1,0:ng+1))
        do k=0,ng+1
          do j=0,ng+1
            if(is_2d .and. j > 1) then
              f(:,j,k) = f(:,1,k); fn(:,:,j,k) = fn(:,:,1,k)
              cycle
            end if
            do i=0,ng+1
              f(i,j,k) = frac_sphere(xc,yc,zc,r,(i-1)*h,(j-1)*h,(k-1)*h,h,h,h,is_2d,fn(:,i,j,k))
            end do
          end do
        end do
        if(is_2d) then
          f(:,0,:) = f(:,1,:); fn(:,:,0,:) = fn(:,:,1,:)
        end if
        nn = 0; nv = 0; e2 = 0._rp; einf = 0._rp; ean = 0._rp; eanmax = 0._rp
        do k=5,ng-4
          do j=merge(1,5,is_2d),merge(1,ng-4,is_2d)
            do i=5,ng-4
              if(f(i,j,k) < 1.e-6_rp .or. f(i,j,k) > 1._rp-1.e-6_rp) cycle
              nn = nn+1
              ! in 2D, the y columns are periodic copies: use f with periodic wrap along y
              call hf_cell(nq,f,i,j,k,h,act,kappa,nrm,is_valid)
              if(.not.is_valid) cycle
              nv = nv+1
              e2 = e2+((kappa-kex)/kex)**2
              einf = max(einf,abs(kappa-kex)/kex)
              nr = fn(:,i,j,k)/sqrt(sum(fn(:,i,j,k)**2))
              rr = acos(min(max(sum(nrm*nr),-1._rp),1._rp))
              ean = ean+rr; eanmax = max(eanmax,rr)
            end do
          end do
        end do
        frac = 1._rp*nv/max(nn,1)
        res(:,ires) = [sqrt(e2/max(nv,1)),einf,frac,ean/max(nv,1),eanmax]
        deallocate(f,fn)
      end do
      print '(a,i2,a,a)', ' height functions, columns of ',2*nq+1,' cells, ',merge('circles (2D)','spheres     ',is_2d)
      print '(a)', '    R/h   curvature rel. error L2 / Linf    valid fraction    normal angle error mean / max [rad]'
      do ires=1,nres
        print '(i7,3x,2es11.2,6x,f8.3,6x,2es11.2)', rres(ires),res(1,ires),res(2,ires),res(3,ires),res(4,ires),res(5,ires)
      end do
      print '(a,2f6.2)', '    observed order (last two resolutions), curvature L2 and mean normal error: ', &
            log(res(1,nres-1)/res(1,nres))/log(2._rp),log(res(4,nres-1)/res(4,nres))/log(2._rp)
      if(log(res(1,nres-1)/res(1,nres))/log(2._rp) < 1.5_rp) then
        print*, '  FAILED: height-function curvature not second order'
        passed = .false.
      end if
    end do
    end do
  end subroutine test_hf_spheres
#if defined(_OPENACC)
  subroutine test_device(passed)
    !
    ! the kernels evaluated on the device must reproduce the host results
    !
    logical, intent(inout) :: passed
    integer, parameter :: n = 100000
    real(rp), allocatable, dimension(:) :: mx,my,mz,v,alpha_h,alpha_d,vol_h,vol_d
    real(rp), allocatable, dimension(:,:,:,:) :: cb
    real(rp), allocatable, dimension(:,:) :: nrm_h,nrm_d
    real(rp) :: c(-1:1,-1:1,-1:1),dz(-1:1),a,b,e,err1,err2,err3
    integer :: q
    allocate(mx(n),my(n),mz(n),v(n),alpha_h(n),alpha_d(n),vol_h(n),vol_d(n),cb(-1:1,-1:1,-1:1,n), &
             nrm_h(3,n),nrm_d(3,n))
    dz(:) = [0.8_rp,1._rp,1.25_rp]
    do q=1,n
      call random_normal(1+mod(q,6),mx(q),my(q),mz(q))
      call random_number(v(q))
      alpha_h(q) = plic_alpha(mx(q),my(q),mz(q),v(q))
      vol_h(q) = plic_volume(mx(q),my(q),mz(q),alpha_h(q)+0.1_rp*mx(q))
      call block_plane(mx(q),my(q),mz(q),alpha_h(q),1._rp,1._rp,dz,cb(:,:,:,q))
      c(:,:,:) = cb(:,:,:,q)
      call plic_normal_myc(c,1._rp,1._rp,dz,[.true.,.true.,.true.],nrm_h(1,q),nrm_h(2,q),nrm_h(3,q))
    end do
    !$acc data copyin(mx,my,mz,v,cb,dz) copyout(alpha_d,vol_d,nrm_d)
    !$acc parallel loop private(c,a,b,e) default(present)
    do q=1,n
      alpha_d(q) = plic_alpha(mx(q),my(q),mz(q),v(q))
      vol_d(q) = plic_volume(mx(q),my(q),mz(q),alpha_d(q)+0.1_rp*mx(q))
      c(:,:,:) = cb(:,:,:,q)
      call plic_normal_myc(c,1._rp,1._rp,dz,[.true.,.true.,.true.],a,b,e)
      nrm_d(1,q) = a; nrm_d(2,q) = b; nrm_d(3,q) = e
    end do
    !$acc end data
    err1 = maxval(abs(alpha_d-alpha_h))
    err2 = maxval(abs(vol_d-vol_h))
    err3 = maxval(abs(nrm_d-nrm_h))
    print '(a,3es10.2)', ' device: max |device-host| for alpha, volume, MYC normal:', err1,err2,err3
    if(max(err1,err2,err3) > 1.e-13_rp) then
      print*, '  FAILED: device results differ from host results'
      passed = .false.
    end if
  end subroutine test_device
#endif
end program test_plic
