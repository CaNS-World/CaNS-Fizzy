block
  use mod_output, only: cmpt_mean_mass,cmpt_mean_energy
  real(rp), dimension(2) :: en12,mass12
  call cmpt_mean_mass(n,l,dl,dzf,rho12,psi,mass12)
  call cmpt_mean_energy(n,l,dl,dzf,rho12,psi,u,v,w,en12)
  var(1) = 1.*istep
  var(2) = time
  var(3) = en12(1)
  var(4) = en12(2)
  call out0d(trim(datadir)//'log_energy.out',4,var)
  var(3) = mass12(1)
  var(4) = mass12(2)
  call out0d(trim(datadir)//'log_mass.out',4,var)
end block
