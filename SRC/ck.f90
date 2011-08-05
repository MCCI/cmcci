function ck(i_sx2,n_2p,kck)
  use precision
  integer       :: i_sx2, n_2p, kck
  integer       :: m1, m2, m3
  real(kind=pr) :: x, y, xm1, xm2, xmult, binomial, ck

  !Sanibel coefficients: principle case

  m1 = 2*i_sx2+2
  m1 = (-1)**kck*m1 
  m2 = n_2p+i_sx2+2

  xm1 = float(m1)
  xm2 = float(m2)
  xmult =  xm1 / xm2

  y = float(n_2p)
  y = y + float(i_sx2)
  m3 =nint( y/2.0_pr )

  x=nint(binomial(m3,kck))

  ck = xmult/x
  !ck = 1.0_pr  !for testing only!!


end function ck

function binomial(n,k)
  use precision
  integer       :: k,n
  real(kind=pr) :: binomial, factorial
  binomial = exp(factorial(n)-factorial(k)-factorial(n-k))
  return
end function binomial

function factorial(n)
  use precision
  integer       :: n
  real(kind=pr) :: a(100), gamma, factorial, xx
  save a
  data a/100*-1.0/
  if(n.lt.0) STOP'error: factorial of negative integer'
  xx =float(n+1)
  if(n.le.99) then
     if(a(n+1).lt.0.0) a(n+1) = gamma(xx)
     factorial=a(n+1)
  else
     factorial=gamma(xx)
  endif
  return
end function factorial

function gamma(xx)
  use precision
  real(kind=pr) :: gamma,xx
  real(kind=pr) :: ser, stp, tmp, x, y, cof(6)
  integer       :: j
  save    cof,stp
  data    cof,stp /76.18009172947146d0,&
       -86.50532032941677d0,&
       24.01409824083091d0,&
       -1.231739572450155d0,&
       .1208650973866179d-2,&
       -.5395239384953d-5,&
       2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp 
  ser=1.000000000190015d0
  do j=1,6
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gamma = tmp+log(stp*ser/x)
  return 
end function gamma
