module cagd

implicit none
contains

subroutine bezier(t,n,d,P,C)
real,intent(in):: t
integer, intent(in):: n,d
real, dimension(0:n,1:d),intent(in):: P
real,dimension(1:d),intent(out):: C
integer:: i,j
real, dimension(0:n,1:d):: Q

Q(0:n,1:d)=P(0:n,1:d)

do j=1,n
	do i=0,n-j
		Q(i,1:d)=(1-t)*Q(i,1:d)+t*Q(i+1,1:d)
	end do
end do

C(1:d)=Q(0,1:d)

end subroutine bezier

subroutine deboor(t,n,k,d,tau,P,C)
real,intent(in) :: t
integer, intent(in):: n,k,d
real, dimension(0:n+k+1),intent(in)::tau
real, dimension(0:n,1:d),intent(in):: P
real,dimension(1:d),intent(out):: C
real, dimension(0:n,1:d):: Q
integer :: r,i,j
real :: ratio
r=n
do i=k,n-1
	if(t<tau(i+1)) then
		r=i
		exit
	end if
end do

Q(0:n,1:d)=P(0:n,1:d)

do j=1,k
	do i=r-k,r-j
		if (tau(i+k+1)==tau(i+j)) then
			Q(i,1:d)=0.0
		else
			ratio=(t-tau(i+j))/(tau(i+k+1)-tau(i+j))
			Q(i,1:d)=(1.0-ratio)*Q(i,1:d)+ratio*Q(i+1,1:d)		
		end if
	end do
end do
C(1:d)=Q(r-k,1:d)

end subroutine deboor

end module cagd


program curves

use cagd
implicit none
integer,parameter:: dm=160,d=2,np=512
real, dimension(0:dm,1:d)::P,Q
real, dimension(0:2*dm+1) :: tau
real, dimension(1:d)::C
integer:: n,k=3,l,i
real::t,h,a=0.0,b=1.0
character(len=1) :: ch

open(unit=2,file='p.txt')
read(2,*) ch,n
do i=0,n
	read(2,*) P(i,1:d)
end do


close(unit=2)

do i=0,k
	tau(i) =a
end do
h=(b-a)/real(n+1-k)
do i=k+1,n
	tau(i)= a+(i-k) *h
end do
do i=n+1,n+k+1
	tau(i)=b
end do

open(unit=1,file='c.txt')

h=(b-a)/np

do l=0,np
	t=a+l*h
	call deboor(t,n,k,d,tau(0:n+k+1),P(0:n,1:d),C(1:d))
	write(1,*) C(1:d)
end do

close(unit=1)

end program curves
