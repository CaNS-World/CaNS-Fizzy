block
  !
  ! total area of the interface, and profiles along z of the statistics of each phase;
  ! the two routines run on the host: include this block where the fields are up to date there,
  ! e.g., in `out1d.h90` (the profiles take a scalar field, here the first one)
  !
  use mod_output, only: cmpt_total_area,velstats
  real(rp) :: area
  call cmpt_total_area(n,nh_psi,nh_g,dli,dzci,dzfi,psi,area)
  var(1) = 1.*istep
  var(2) = time
  var(3) = area
  call out0d(trim(datadir)//'log_area.out',3,var)
  call velstats(trim(datadir)//'velstats_two_fluid_fld_'//fldnum//'.out',ng,nh_v,nh_p,nh_s,nh_psi,nh_g,lo,hi,l,dl,zc_g, &
                u,v,w,p,scalars(1)%val,psi)
end block
