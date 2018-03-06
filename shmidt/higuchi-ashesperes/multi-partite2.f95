program payu
implicit none
integer,parameter::k=2
double precision,dimension(0:1,0:k-1)::A1,B1,A2
double precision,dimension(0:2**k-1,0:0)::B2,psi,psi2
double precision,dimension(0:0,0:0)::r1
double precision,dimension(0:1,0:0)::W1,W2
double precision::r2,r11
integer::i,j,l,m,n,seed,count
count=999
r2=0
psi(0:3,0)=(/1.0/sqrt(3.0),1.0/sqrt(3.0),0.0,1.0/sqrt(3.0)/)
call random_seed(seed)
do i=0,count
do j=0,k-1,1
do l=0,1
call random_number(A1(l,j))
A2(l,j)=sqrt(-2*log(A1(l,j)))*cos(8*atan(1.0)*A1(l,j))
enddo
enddo
do m=0,k-1
B1(0,m)=A2(0,m)/(sqrt(A2(0,m)**2+A2(1,m)**2))
B1(1,m)=A2(1,m)/(sqrt(A2(0,m)**2+A2(1,m)**2))
enddo
!print*,B1(1,0),B1(1,1)
call tensor_pro1(B1(0:1,0:0),B1(0:1,1:1),B2,2,1,2,1)
r1=matmul(transpose(B2),psi)
r11=r1(0,0)
if (r11**2>r2) then
r2=r11
W1=B1(0:1,0:0)
W2=B1(0:1,1:1)
else
continue
endif

enddo
print*,W1,W2,r2
do i=0,2**k-1,1
psi2(i,0)=psi(i,0)-r2
enddo

end program payu

subroutine tensor_pro1(T1,T2,T3,m1,n1,m2,n2)
implicit none
integer::m1,n1,m2,n2
double precision,dimension(0:m1-1,0:n1-1)::T1
double precision,dimension(0:m2-1,0:n2-1)::T2
double precision,dimension(0:m1*m2-1,0:n1*n2-1)::T3
integer::i,j,k,l

do i=0,m1-1,1
do j=0,n1-1,1
do k=0,m2-1,1
do l=0,n2-1,1
T3(i*m2+k,j*n2+l)=T1(i,j)*T2(k,l)
enddo
enddo
enddo    
enddo
    
end subroutine tensor_pro1
