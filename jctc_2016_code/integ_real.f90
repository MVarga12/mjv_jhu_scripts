  integer :: N
  real, allocatable :: X(:), Y(:), Yfit(:)
  real, dimension(-3:3) :: W = (/ 1,3,6,7,6,3,1 /), WW

! input data must be in file "in" in a single column
  open (1,file='productregion.dat')
! find length of data in file "in" and read it in X
  N=0
  do
    read (1,*,iostat=ios) yyy; if (ios /= 0) exit; N=N+1
  enddo
  rewind (1)
  allocate ( X(N), Y(N), Yfit(N) )
  read (1,*) ( Y(i), i=1,N)

! do a weighted moving average of the input data
Yfit = 0
do i=1,N
   WW = W
do m=-3,3
   if (i+m>N .or. i+M<1) then  ! not enough points near edges
      WW(m) = 0.
      cycle
   endif
   Yfit(i) = Yfit(i) + W(m)*Y(i+m)/ sum(WW)
enddo
enddo

! write fitted values in file "fit"
open (2,file='fit')
do i=1,N
   write (2,*) i-1, Yfit(i)
enddo

! write the residuals: real values - fitted values
  open (3,file='resid')
  do i=1,N
     write (3,*) i-1, Y(i)-Yfit(i)
  enddo

! calculate the integral over fitted values using simpson's rule
  simps = 0.
  simps = simps + 1./3.*Yfit(1) + 1./3.*Yfit(N)
  do i=2,N-1
     if ( mod(i,2)==0 ) simps = simps + 4./3.*Yfit(i)
     if ( mod(i,2)==1 ) simps = simps + 2./3.*Yfit(i)
  enddo

  print *, simps

  end
