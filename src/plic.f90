! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_plic
  !
  ! geometric kernels of the piecewise-linear interface calculation (PLIC) in a rectangular cell
  !
  ! all kernels are formulated in cell-scaled coordinates, i.e., in the unit cube xi = (x-x0)/dx in [0,1]^3;
  ! the interface is the plane m.xi = alpha, with phase 1 occupying the region m.xi <= alpha, i.e., m points
  ! out of phase 1; a physical normal n maps to cell-scaled coordinates as m = (n_x*dx,n_y*dy,n_z*dz)
  !
  ! references:
  !   R. Scardovelli & S. Zaleski, J. Comput. Phys. 164 (2000) 228--237 (plane--volume relations)
  !   E. Aulisa, S. Manservisi, R. Scardovelli & S. Zaleski, J. Comput. Phys. 225 (2007) 2301--2319 (MYC normals)
  !
  use mod_types
  implicit none
  private
  public plic_volume,plic_alpha,plic_volume_box,plic_normal_myc,plic_normal_youngs,plic_normal_cc, &
         hf_heights,hf_curvature,hf_normal
  contains
  pure real(rp) function plic_volume(mx,my,mz,alpha) result(vol)
    !
    ! volume fraction of the unit cube in the region m.xi <= alpha
    !
    implicit none
    real(rp), intent(in) :: mx,my,mz,alpha
    real(rp) :: s,a,b1,b2,b3
    !$acc routine seq
    s = abs(mx)+abs(my)+abs(mz)
    a = alpha-min(mx,0._rp)-min(my,0._rp)-min(mz,0._rp) ! mirror the cube such that m >= 0
    if(     a <= 0._rp) then
      vol = 0._rp
    else if(a >= s    ) then
      vol = 1._rp
    else
      a = a/s
      call sort3(abs(mx)/s,abs(my)/s,abs(mz)/s,b1,b2,b3)
      if(a <= 0.5_rp) then
        vol = vol_half(b1,b2,b3,a)
      else
        vol = 1._rp-vol_half(b1,b2,b3,1._rp-a)
      end if
    end if
  end function plic_volume
  !
  pure real(rp) function plic_alpha(mx,my,mz,vol) result(alpha)
    !
    ! plane constant alpha such that the region m.xi <= alpha has volume fraction vol
    !
    implicit none
    real(rp), intent(in) :: mx,my,mz,vol
    real(rp) :: s,v,a,b1,b2,b3
    !$acc routine seq
    s = abs(mx)+abs(my)+abs(mz)
    v = min(max(vol,0._rp),1._rp)
    a = 0._rp
    if(s > 0._rp) then
      call sort3(abs(mx)/s,abs(my)/s,abs(mz)/s,b1,b2,b3)
      if(v <= 0.5_rp) then
        a = alpha_half(b1,b2,b3,v)
      else
        a = 1._rp-alpha_half(b1,b2,b3,1._rp-v)
      end if
    end if
    alpha = a*s+min(mx,0._rp)+min(my,0._rp)+min(mz,0._rp)
  end function plic_alpha
  !
  pure real(rp) function plic_volume_box(mx,my,mz,alpha,x0,y0,z0,dx,dy,dz) result(vol)
    !
    ! volume fraction of the box [x0,x0+dx]x[y0,y0+dy]x[z0,z0+dz] (in cell-scaled coordinates)
    ! occupied by phase 1, relative to the volume of the box
    !
    implicit none
    real(rp), intent(in) :: mx,my,mz,alpha,x0,y0,z0,dx,dy,dz
    !$acc routine seq
    vol = plic_volume(mx*dx,my*dy,mz*dz,alpha-mx*x0-my*y0-mz*z0)
  end function plic_volume_box
  !
  pure subroutine plic_normal_myc(c,dx,dy,dz,is_active,mx,my,mz)
    !
    ! mixed Youngs-centred (MYC) estimate of the interface normal (Aulisa et al. 2007)
    !
    ! c(-1:1,-1:1,-1:1) are the volume fractions in the 3x3x3 block of cells centred on the cell of interest;
    ! dx and dy are the (uniform) grid spacings, and dz(-1:1) the sizes of the three cell layers along z;
    ! is_active(:) is false along directions in which the problem is invariant (e.g., 2D runs);
    ! the normal is returned in the cell-scaled coordinates of the central cell, pointing out of phase 1,
    ! and normalized such that |mx|+|my|+|mz| = 1
    !
    implicit none
    real(rp), intent(in ), dimension(-1:1,-1:1,-1:1) :: c
    real(rp), intent(in ) :: dx,dy
    real(rp), intent(in ), dimension(-1:1) :: dz
    logical , intent(in ), dimension(3) :: is_active
    real(rp), intent(out) :: mx,my,mz
    real(rp) :: mx_y,my_y,mz_y,t,tcc
    !$acc routine seq
    call normal_cc(c,dx,dy,dz,is_active,mx,my,mz,tcc)
    call plic_normal_youngs(c,dx,dy,dz,mx_y,my_y,mz_y)
    !
    ! a centred-columns normal that is more aligned with its column direction than the Youngs normal
    ! is with any grid direction is likely affected by columns that do not contain the interface;
    ! use Youngs' estimate instead
    !
    t = max(abs(mx_y),abs(my_y),abs(mz_y))
    if(tcc > t) then
      mx = mx_y
      my = my_y
      mz = mz_y
    end if
  end subroutine plic_normal_myc
  !
  pure subroutine plic_normal_youngs(c,dx,dy,dz,mx,my,mz)
    !
    ! Youngs' estimate of the interface normal: minus the gradient of the volume fraction,
    ! computed with the (1,2,1)-weighted 3x3x3 stencil; same conventions as `plic_normal_myc`
    !
    implicit none
    real(rp), intent(in ), dimension(-1:1,-1:1,-1:1) :: c
    real(rp), intent(in ) :: dx,dy
    real(rp), intent(in ), dimension(-1:1) :: dz
    real(rp), intent(out) :: mx,my,mz
    real(rp), parameter, dimension(-1:1) :: w = [1._rp,2._rp,1._rp]
    real(rp) :: gx,gy,gz,dzm,dzp,s
    integer :: p,q
    !$acc routine seq
    dzm = 0.5_rp*(dz(-1)+dz(0))
    dzp = 0.5_rp*(dz( 0)+dz(1))
    gx = 0._rp
    gy = 0._rp
    gz = 0._rp
    do q=-1,1
      do p=-1,1
        gx = gx + w(p)*w(q)*(c(1,p,q)-c(-1,p,q))
        gy = gy + w(p)*w(q)*(c(p,1,q)-c(p,-1,q))
        gz = gz + w(p)*w(q)*dcentred(c(p,q,-1),c(p,q,0),c(p,q,1),dzm,dzp)
      end do
    end do
    mx = -gx/(32._rp*dx)*dx
    my = -gy/(32._rp*dy)*dy
    mz = -gz/16._rp*dz(0)
    s = abs(mx)+abs(my)+abs(mz)
    if(s > 0._rp) then
      mx = mx/s
      my = my/s
      mz = mz/s
    else
      mx = 0._rp
      my = 0._rp
      mz = 1._rp
    end if
  end subroutine plic_normal_youngs
  !
  pure subroutine plic_normal_cc(c,dx,dy,dz,is_active,mx,my,mz)
    !
    ! centred-columns estimate of the interface normal: for each active direction, the interface is
    ! written as a height function of the other two coordinates, with heights given by column sums of
    ! the volume fractions and slopes by centred differences; among the candidates, the one with the
    ! largest component along its own column direction is returned; same conventions as `plic_normal_myc`
    !
    implicit none
    real(rp), intent(in ), dimension(-1:1,-1:1,-1:1) :: c
    real(rp), intent(in ) :: dx,dy
    real(rp), intent(in ), dimension(-1:1) :: dz
    logical , intent(in ), dimension(3) :: is_active
    real(rp), intent(out) :: mx,my,mz
    real(rp) :: tmax
    !$acc routine seq
    call normal_cc(c,dx,dy,dz,is_active,mx,my,mz,tmax)
  end subroutine plic_normal_cc
  !
  pure subroutine normal_cc(c,dx,dy,dz,is_active,mx,my,mz,tmax)
    !
    ! same as `plic_normal_cc`, also returning the component of the chosen normal along its column direction
    !
    implicit none
    real(rp), intent(in ), dimension(-1:1,-1:1,-1:1) :: c
    real(rp), intent(in ) :: dx,dy
    real(rp), intent(in ), dimension(-1:1) :: dz
    logical , intent(in ), dimension(3) :: is_active
    real(rp), intent(out) :: mx,my,mz,tmax
    real(rp) :: dzm,dzp,hm,h0,hp,hmy,hpy,sgn,ax,ay,az,s,t
    integer :: p
    !$acc routine seq
    dzm = 0.5_rp*(dz(-1)+dz(0))
    dzp = 0.5_rp*(dz( 0)+dz(1))
    tmax = -1._rp
    mx = 0._rp
    my = 0._rp
    mz = 1._rp
    if(is_active(1)) then
      !
      ! columns along x: x = h(y,z)
      !
      hmy = sum(c(:,-1,0))*dx
      hpy = sum(c(:, 1,0))*dx
      hm  = sum(c(:,0,-1))*dx
      h0  = sum(c(:,0, 0))*dx
      hp  = sum(c(:,0, 1))*dx
      sgn = sign(1._rp,sum(c(-1,:,:))-sum(c(1,:,:)))
      ax = sgn*dx
      ay = -(hpy-hmy)/(2._rp*dy)*dy
      az = -dcentred(hm,h0,hp,dzm,dzp)*dz(0)
      s = abs(ax)+abs(ay)+abs(az)
      t = abs(ax)/s
      if(t > tmax) then
        tmax = t
        mx = ax/s; my = ay/s; mz = az/s
      end if
    end if
    if(is_active(2)) then
      !
      ! columns along y: y = h(x,z)
      !
      hmy = sum(c(-1,:,0))*dy
      hpy = sum(c( 1,:,0))*dy
      hm  = sum(c(0,:,-1))*dy
      h0  = sum(c(0,:, 0))*dy
      hp  = sum(c(0,:, 1))*dy
      sgn = sign(1._rp,sum(c(:,-1,:))-sum(c(:,1,:)))
      ax = -(hpy-hmy)/(2._rp*dx)*dx
      ay = sgn*dy
      az = -dcentred(hm,h0,hp,dzm,dzp)*dz(0)
      s = abs(ax)+abs(ay)+abs(az)
      t = abs(ay)/s
      if(t > tmax) then
        tmax = t
        mx = ax/s; my = ay/s; mz = az/s
      end if
    end if
    if(is_active(3)) then
      !
      ! columns along z: z = h(x,y)
      !
      hm  = 0._rp
      hp  = 0._rp
      hmy = 0._rp
      hpy = 0._rp
      do p=-1,1
        hm  = hm  + c(-1,0,p)*dz(p)
        hp  = hp  + c( 1,0,p)*dz(p)
        hmy = hmy + c(0,-1,p)*dz(p)
        hpy = hpy + c(0, 1,p)*dz(p)
      end do
      sgn = sign(1._rp,sum(c(:,:,-1))-sum(c(:,:,1)))
      ax = -(hp -hm )/(2._rp*dx)*dx
      ay = -(hpy-hmy)/(2._rp*dy)*dy
      az = sgn*dz(0)
      s = abs(ax)+abs(ay)+abs(az)
      t = abs(az)/s
      if(t > tmax) then
        tmax = t
        mx = ax/s; my = ay/s; mz = az/s
      end if
    end if
  end subroutine normal_cc
  !
  pure subroutine hf_heights(nq,c,dq,eps,hgt,sgn,is_valid)
    !
    ! height functions from the volume fractions c(-1:1,-1:1,-nq:nq) in a 3x3 set of (2*nq+1)-cell columns along
    ! a direction, centred on the cell of interest; the last index runs along the columns, and dq(-nq:nq)
    ! are the cell sizes along them; the columns are valid if their end cells are pure, with phase 1 at
    ! the same end of all columns (sgn = +1 if phase 1 is at the lower end, -1 otherwise); the heights
    ! hgt(-1:1,-1:1) are the interface positions along the columns, relative to the centre of the central cell
    !
    ! references: Cummins et al., Comput. Struct. 83 (2005) 425--434; Popinet, J. Comput. Phys. 228 (2009) 5838--5866
    !
    implicit none
    integer , intent(in ) :: nq
    real(rp), intent(in ), dimension(-1:1,-1:1,-nq:nq) :: c
    real(rp), intent(in ), dimension(-nq:nq) :: dq
    real(rp), intent(in ) :: eps
    real(rp), intent(out), dimension(-1:1,-1:1) :: hgt
    real(rp), intent(out) :: sgn
    logical , intent(out) :: is_valid
    real(rp) :: lo,hi,vol
    integer :: a,b,q
    !$acc routine seq
    lo = 0.5_rp*dq(0) ! distance from the lower end of the columns to the central cell centre
    hi = 0.5_rp*dq(0) ! distance from the central cell centre to the upper end
    do q=1,nq
      lo = lo + dq(-q)
      hi = hi + dq( q)
    end do
    sgn = 0._rp
    if(     c(0,0,-nq) >= 1._rp-eps .and. c(0,0,nq) <= eps) then
      sgn =  1._rp
    else if(c(0,0,-nq) <= eps .and. c(0,0,nq) >= 1._rp-eps) then
      sgn = -1._rp
    end if
    is_valid = sgn /= 0._rp
    hgt(:,:) = 0._rp
    do b=-1,1
      do a=-1,1
        if(sgn > 0._rp) then
          is_valid = is_valid .and. c(a,b,-nq) >= 1._rp-eps .and. c(a,b,nq) <= eps
        else
          is_valid = is_valid .and. c(a,b,-nq) <= eps .and. c(a,b,nq) >= 1._rp-eps
        end if
        vol = 0._rp
        do q=-nq,nq
          vol = vol + c(a,b,q)*dq(q)
        end do
        hgt(a,b) = merge(vol-lo,hi-vol,sgn > 0._rp)
      end do
    end do
  end subroutine hf_heights
  !
  pure real(rp) function hf_curvature(hgt,sgn,dam,dap,dbm,dbp,is_act_a,is_act_b) result(kappa)
    !
    ! curvature from the heights hgt(-1:1,-1:1) of `hf_heights`, with the sign convention of Fizzy
    ! (kappa = -div(n), n pointing into phase 1: positive for a drop of phase 1); dam, dap (dbm, dbp) are
    ! the distances to the neighbouring columns along the first (second) transverse direction, and
    ! is_act_a, is_act_b are false along directions in which the problem is invariant
    !
    implicit none
    real(rp), intent(in), dimension(-1:1,-1:1) :: hgt
    real(rp), intent(in) :: sgn,dam,dap,dbm,dbp
    logical , intent(in) :: is_act_a,is_act_b
    real(rp) :: ha,hb,haa,hbb,hab,hbm,hbp
    !$acc routine seq
    ha  = 0._rp; haa = 0._rp
    hb  = 0._rp; hbb = 0._rp
    hab = 0._rp
    if(is_act_a) then
      ha  = dcentred(hgt(-1,0),hgt(0,0),hgt(1,0),dam,dap)
      haa = 2._rp*((hgt(1,0)-hgt(0,0))/dap-(hgt(0,0)-hgt(-1,0))/dam)/(dam+dap)
    end if
    if(is_act_b) then
      hb  = dcentred(hgt(0,-1),hgt(0,0),hgt(0,1),dbm,dbp)
      hbb = 2._rp*((hgt(0,1)-hgt(0,0))/dbp-(hgt(0,0)-hgt(0,-1))/dbm)/(dbm+dbp)
    end if
    if(is_act_a .and. is_act_b) then
      hbm = dcentred(hgt(-1,-1),hgt(-1,0),hgt(-1,1),dbm,dbp)
      hbp = dcentred(hgt( 1,-1),hgt( 1,0),hgt( 1,1),dbm,dbp)
      hab = dcentred(hbm,hb,hbp,dam,dap)
    end if
    kappa = -sgn*(haa*(1._rp+hb**2)+hbb*(1._rp+ha**2)-2._rp*hab*ha*hb)/(1._rp+ha**2+hb**2)**1.5_rp
  end function hf_curvature
  !
  pure subroutine hf_normal(hgt,sgn,dam,dap,dbm,dbp,is_act_a,is_act_b,na,nb,nq)
    !
    ! physical unit normal pointing out of phase 1, from the heights of `hf_heights`, with components
    ! along the first and second transverse directions (na, nb) and along the columns (nq)
    !
    implicit none
    real(rp), intent(in ), dimension(-1:1,-1:1) :: hgt
    real(rp), intent(in ) :: sgn,dam,dap,dbm,dbp
    logical , intent(in ) :: is_act_a,is_act_b
    real(rp), intent(out) :: na,nb,nq
    real(rp) :: ha,hb,s
    !$acc routine seq
    ha = 0._rp
    hb = 0._rp
    if(is_act_a) ha = dcentred(hgt(-1,0),hgt(0,0),hgt(1,0),dam,dap)
    if(is_act_b) hb = dcentred(hgt(0,-1),hgt(0,0),hgt(0,1),dbm,dbp)
    s  = sqrt(1._rp+ha**2+hb**2)
    na = -sgn*ha/s
    nb = -sgn*hb/s
    nq =  sgn/s
  end subroutine hf_normal
  !
  pure real(rp) function dcentred(qm,q0,qp,dzm,dzp) result(dq)
    !
    ! second-order centred first derivative on a (possibly) non-uniform three-point stencil;
    ! dzm and dzp are the distances from the central point to the lower and upper points
    !
    implicit none
    real(rp), intent(in) :: qm,q0,qp,dzm,dzp
    !$acc routine seq
    dq = (dzm**2*(qp-q0)+dzp**2*(q0-qm))/(dzm*dzp*(dzm+dzp))
  end function dcentred
  !
  pure subroutine sort3(x1,x2,x3,b1,b2,b3)
    implicit none
    real(rp), intent(in ) :: x1,x2,x3
    real(rp), intent(out) :: b1,b2,b3
    !$acc routine seq
    b1 = min(x1,x2,x3)
    b3 = max(x1,x2,x3)
    b2 = max(min(x1,x2),min(max(x1,x2),x3))
  end subroutine sort3
  !
  pure real(rp) function vol_half(b1,b2,b3,a) result(vol)
    !
    ! volume of the unit cube below the plane b.xi = a, with 0 <= b1 <= b2 <= b3, b1+b2+b3 = 1,
    ! and 0 <= a <= 1/2 (Scardovelli & Zaleski 2000); the expressions are arranged such that no
    ! cancellation errors arise for vanishing b1 (and b2), and no division by zero takes place
    !
    implicit none
    real(rp), intent(in) :: b1,b2,b3,a
    real(rp) :: b12,d2,d3
    !$acc routine seq
    b12 = b1+b2
    if(     a < b1) then
      vol = a*(a/b1)*(a/b2)/(6._rp*b3)
    else if(a < b2) then
      vol = (3._rp*a*(a-b1)+b1**2)/(6._rp*b2*b3)
    else if(a < min(b12,b3)) then
      d2 = a-b2
      vol = (3._rp*a*(a-b1)+b1**2-d2**2*(d2/b1))/(6._rp*b2*b3)
    else if(b12 <= b3) then
      vol = (a-0.5_rp*b12)/b3
    else
      d2 = a-b2
      d3 = a-b3
      vol = (3._rp*a*(a-b1)+b1**2-d2**2*(d2/b1)-d3**2*(d3/b1))/(6._rp*b2*b3)
    end if
  end function vol_half
  !
  pure real(rp) function dvol_half(b1,b2,b3,a) result(dvol)
    !
    ! derivative of `vol_half` with respect to a
    !
    implicit none
    real(rp), intent(in) :: b1,b2,b3,a
    real(rp) :: b12,d2,d3
    !$acc routine seq
    b12 = b1+b2
    if(     a < b1) then
      dvol = (a/b1)*(a/b2)/(2._rp*b3)
    else if(a < b2) then
      dvol = (2._rp*a-b1)/(2._rp*b2*b3)
    else if(a < min(b12,b3)) then
      d2 = a-b2
      dvol = (2._rp*a-b1-d2*(d2/b1))/(2._rp*b2*b3)
    else if(b12 <= b3) then
      dvol = 1._rp/b3
    else
      d2 = a-b2
      d3 = a-b3
      dvol = (2._rp*a-b1-d2*(d2/b1)-d3*(d3/b1))/(2._rp*b2*b3)
    end if
  end function dvol_half
  !
  pure real(rp) function alpha_half(b1,b2,b3,vol) result(a)
    !
    ! inverse of `vol_half`: plane constant for a volume 0 <= vol <= 1/2;
    ! closed form for the first two and the linear regimes, and bracketed Newton
    ! iterations for the cubic ones
    !
    implicit none
    integer , parameter :: itermax = 50
    real(rp), intent(in) :: b1,b2,b3,vol
    real(rp) :: b12,a3,v1,v2,v3,alo,ahi,vlo,vhi,f,df,anew
    integer :: iter
    !$acc routine seq
    b12 = b1+b2
    a3  = min(b12,b3)
    v1 = 0._rp
    v2 = 0._rp
    if(b1 > 0._rp) v1 = b1**2/(6._rp*b2*b3)
    if(b2 > 0._rp) v2 = (b2-b1)/(2._rp*b3)+b1**2/(6._rp*b2*b3)
    v3 = vol_half(b1,b2,b3,a3)
    if(     vol < v1) then
      a = (6._rp*b1*b2*b3*vol)**(1._rp/3._rp)
    else if(vol < v2) then
      a = 0.5_rp*b1+sqrt(max(2._rp*b2*b3*vol-b1**2/12._rp,0._rp))
    else if(vol >= v3 .and. b12 <= b3) then
      a = b3*vol+0.5_rp*b12
    else
      if(vol < v3) then
        alo = b2; vlo = v2
        ahi = a3; vhi = v3
      else
        alo = b3; vlo = v3
        ahi = 0.5_rp; vhi = 0.5_rp
      end if
      a = alo
      if(vhi > vlo) a = alo+(ahi-alo)*(vol-vlo)/(vhi-vlo)
      do iter=1,itermax
        f = vol_half(b1,b2,b3,a)-vol
        if(f > 0._rp) then
          ahi = a
        else
          alo = a
        end if
        df = dvol_half(b1,b2,b3,a)
        anew = 0.5_rp*(alo+ahi)
        if(df > 0._rp) anew = a-f/df
        if(anew <= alo .or. anew >= ahi) anew = 0.5_rp*(alo+ahi)
        if(abs(anew-a) <= epsilon(1._rp)) then
          a = anew
          exit
        end if
        a = anew
      end do
    end if
  end function alpha_half
end module mod_plic
