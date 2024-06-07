block
  real(rp) :: phic,phip,h
  integer  :: i,j,k
  !$acc wait
  !$acc update self(phi)
  h = 0.
  j = 1
  if(lo(1) == 1    ) then
    do k=1,n(3)
      phic = 0.5*(phi(0,j,k  )+phi(1,j,k  ))
      phip = 0.5*(phi(0,j,k+1)+phi(1,j,k+1))
      if(phic*phip < 0.) then
        h = h + zc(k)-phic*dzc(k)/(phip-phic)
      end if
    end do
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE,h,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  h = h - l(3)/2
  var(1) = time
  var(2) = h
  call out0d(trim(datadir)//'amplitude-cap-wav-1d.out',2,var)
end block
