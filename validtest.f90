program main
implicit none

real(kind=8)  :: a,b


b = atan(a/0.0d0)


print *,b


end program
