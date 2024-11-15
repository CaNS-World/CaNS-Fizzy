block
  real(rp), dimension(2) :: en12,mass12
  call cmpt_total_energy(n,dl,dzf,rho12,psi,u,v,w,en12)
  call cmpt_total_mass(n,dl,dzf,rho12,psi,mass12)
  var(1) = 1.*istep
  var(2) = time
  var(3) = en12(1)
  var(4) = en12(2)
  call out0d(trim(datadir)//'log_energy.out',4,var)
  var(3) = mass12(1)
  var(4) = mass12(2)
  call out0d(trim(datadir)//'log_mass.out',4,var)
end block