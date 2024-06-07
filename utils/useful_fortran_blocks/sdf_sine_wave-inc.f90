function sdist_to_wave(p,a,b,c,d) result(r)
  !
  ! adapted from: https://www.shadertoy.com/view/3t23WG
  !
  ! example usage: sdist = sdist_to_wave([x,z],0._rp,zfilm_max,2.*pi,0._rp)
  !
  use mod_param, only: pi
  implicit none
  real(rp), intent(in) :: p(2)     ! input point (x, y)
  real(rp), intent(in) :: a,b,c,d  ! parameters of the cosine wave
  real(rp) :: r                    ! resulting signed distance
  real(rp) :: p_transformed(2)
  real(rp) :: w,xa,xb,x,y,x_closest,y_closest,distance
  real(rp) :: qa(2),qb(2),pa(2),ba(2),h
  integer  :: i,isgn
  !
  ! transform to primitive cosine space
  !
  p_transformed = c*(p - [d,a])
  w = c*b
  !
  ! reduce to principal half cycle
  !
  p_transformed(1) = mod(p_transformed(1),2.*pi)
  if (p_transformed(1) > pi) p_transformed(1) = 2.*pi - p_transformed(1)
  !
  ! find zero of derivative using bisection
  !
  xa = 0.
  xb = pi
  do i=1,24
    x = 0.5*(xa + xb)
    y = x - p_transformed(1) + w*sin(x)*(p_transformed(2) - w*cos(x))
    if (y < 0.) then
      xa = x
    else
      xb = x
    end if
  end do
  !
  ! compute distance
  !
  x_closest = 0.5*(xa+xb)
  y_closest = w*cos(x_closest)
  qa = [x_closest,y_closest]
  qb = [xb,w*cos(xb)]
  pa = p_transformed - qa
  ba = qb - qa
  h = dot_product(pa,ba)/dot_product(ba,ba)
  h = max(0._rp,min(1._rp,h))
  distance = sqrt(dot_product(pa-ba*h,pa-ba*h))
  !
  ! determine sign (above or below the wave)
  !
  isgn = sign(1._rp,p_transformed(2)-y_closest)
  !
  ! convert back to non-primitive cosine space
  !
  r = isgn*distance/c
end function sdist_to_wave
