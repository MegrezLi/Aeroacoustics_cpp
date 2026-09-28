! Independent test oracle for published Howe-Chase Eqs. 17-18.
! Complex-step derivative of f, instead of the C++ expanded analytic derivative.
! This is new validation code, not an upstream OpenFAST routine.
module howe_reference
use iso_c_binding
implicit none
contains
function shape(q,h,lambda,delta) result(value) bind(C,name='howe_reference_shape')
real(c_double),value :: q,h,lambda,delta
real(c_double) :: value,eps
complex(c_double_complex) :: chi,r,b,t,f,e,ratio
eps=1.e-25_c_double
chi=cmplx(1.33_c_double,eps,kind=c_double)
r=sqrt(q*q+chi*chi)
b=q*q*(1+(4*h/lambda)**2)+chi*chi
f=1/b
if(h>0)then
  t=lambda/(2*delta)*r
  e=exp(-t)
  ratio=(1+e*e-2*e*cos(2*q*h/delta))/(1-e*e)
  f=f*(1+64*(h/lambda)**3*(delta/h)*q*q*ratio/(r*b))
endif
value=real(f,c_double)+0.5_c_double*1.33_c_double*aimag(f)/eps
end function
end module
