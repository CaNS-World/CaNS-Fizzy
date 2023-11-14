program read_bubbles
  implicit none
  
  ! Declare the derived type
  type bubble
    real :: x, y, z, r
  end type bubble
  
  ! Declare variables
  character(len=100) :: filename
  integer :: num_bubbles, i
  type(bubble), allocatable :: bubbles(:)
  integer :: ierr
  character(len=100) :: line
  real :: x, y, z, r
  
  ! Set the filename
  filename = "bubbles.in"
  
  ! Open the file
  open(unit=10, file=filename, status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    write(*,*) "Error opening file:", trim(filename)
    stop
  endif
  
  ! Count the number of bubbles
  num_bubbles = 0
  ierr = 0
  do while(ierr == 0)
    read(10, '(A)', iostat=ierr) line
    if (ierr /= 0) exit
    num_bubbles = num_bubbles + 1
  end do
  
  ! Allocate the array
  allocate(bubbles(num_bubbles))
  
  ! Rewind the file
  rewind(10)
  
  ! Read the bubble data
  do i = 1, num_bubbles
    read(10, *) x, y, z, r
    bubbles(i)%x = x
    bubbles(i)%y = y
    bubbles(i)%z = z
    bubbles(i)%r = r
  end do
  
  ! Close the file
  close(10)
  
  ! Print the bubble data
  do i = 1, num_bubbles
    write(*, '(4F8.3)') bubbles(i)%x, bubbles(i)%y, bubbles(i)%z, bubbles(i)%r
  end do
  
  ! Deallocate the array
  deallocate(bubbles)
  
end program read_bubbles
