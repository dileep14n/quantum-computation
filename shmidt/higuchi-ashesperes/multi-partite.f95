program multi_schmidt_decom
implicit none

integer,parameter::k=2                           
integer::i,j,i11,i12,m1,n1,m4,n4,l,INFO,s,seed,count    
real::s1,s2,s3,s4,res,res1,res2,res3 ,s11,s33,q1,qq1,q11,qq11,q2,qq2
real::r11,r22               
real,dimension(0:1,0:1)::p1,P,P2                              
real,dimension(0:2**k-1,0:0)::psi,psi2
real,dimension(0:2**k-1,0:2**k-1)::rhoABCDE
real,dimension(0:2**k-1,0:1)::F,F1
real,dimension(0:1,0:0)::B1,B2,B4,B5,B6,B7,B8,B9,A1,A2,A3,A4,A5
real,dimension(0:3,0:0)::B3,R1,R2
real,dimension(2)::W
real,dimension(5)::WORK
count=999
!psi(0:3,0)=(/1.0/sqrt(2.0),0.0,0.0,1.0/sqrt(2.0)/)
psi(0:3,0)=(/1.0/sqrt(2.0),0.0,0.0,1.0/sqrt(2.0)/)
res1=0
rhoABCDE=matmul(psi,transpose(psi))

call random_seed(seed)
do i11=0,count,1
call random_number(s1)
s11=sqrt(-2*log(s1))*cos(8*atan(1.0)*s1)
s2=sqrt(1-s11**2)
B1(0,0)=s11
B1(1,0)=s2
do i12=0,count,1
call random_number(s3)
s33=sqrt(-2*log(s3))*cos(8*atan(1.0)*s3)
s4=sqrt(1-s33**2)
B2(0,0)=s33
B2(1,0)=s4
call tensor_pro1(B1,B2,B3,2,1,2,1)
res=0
do j=0,3
res=res+B3(j,0)*psi(j,0)
enddo
if (res*res>=res1) then
res1=res*res
B4=B1
B5=B2
else
continue
endif
enddo
enddo
do i=0,1
B6(i,0)=B4(i,0)/sqrt(B4(1,0)**2+B4(0,0)**2)
B7(i,0)=B5(i,0)/sqrt(B5(1,0)**2+B5(0,0)**2)
enddo

print*,res1,B6,B7

call random_seed(seed)
call random_number(q1)
call random_number(qq1)
q11=sqrt(-2*log(q1))*cos(8*atan(1.0)*q1)
qq11=sqrt(-2*log(qq1))*cos(8*atan(1.0)*qq1)
q2=sqrt(1-q11**2)
qq2=sqrt(1-qq11**2)
A1(0,0)=q11
A1(1,0)=q2
A2(0,0)=qq11
A2(1,0)=qq2
call projector(B4,A1,A3)
call projector(B5,A2,A4)
!print*,A4,A3


call tensor_pro1(A3,A4,A5,2,1,2,1)
r11=0
do i=0,1
r11=r11+A5(i,0)*psi(i,0)
enddo
print*,r11

end program multi_schmidt_decom
subroutine tensor_pro1(T1,T2,T3,m1,n1,m2,n2)
implicit none
integer::m1,n1,m2,n2
real,dimension(0:m1-1,0:n1-1)::T1
real,dimension(0:m2-1,0:n2-1)::T2
real,dimension(0:m1*m2-1,0:n1*n2-1)::T3
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

subroutine projector(A1,A11,A3)
    implicit none
    real,dimension(0:1,0:0)::A1,A2,A3,A4,A11
    real,dimension(0:1,0:1)::B
    integer::i
    real::k
    B=matmul(A1,transpose(A1))
    A2=matmul(B,A11)
    A3=A11-A2
    k=0
    do i=0,1
    k=k+A3(i,0)**2
    enddo
    !print*,k
    !do i=0,1
    !A4(i,0)=A3(i,0)/k
    
    !enddo
    
end subroutine projector
