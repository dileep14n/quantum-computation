program gdv
implicit none
integer,parameter::k=3
double precision,dimension(0:1,0:k-1)::A1,B1,A2,W1,W3,W5,G1,G0,G2
double precision,dimension(0:2**k-1,0:0)::B2,psi,psi2,B3,B4
double precision,dimension(0:0,0:0)::r1,s1,k2
double precision,dimension(0:1,0:0)::W4,W6
double precision::r2,r11,r3,k1,k3,s11
integer::i,j,l,m,n,seed,count
psi(0:7,0)=(/1.0/sqrt(2.0),1.0/sqrt(2.0),0.0,0.0,0.0,0.0,0.0,0.0/)
seed=32
call random_seed(seed)
r3=0
do  
do i=0,1,1
do j=0,k-1
call random_number(A1(i,j))
A2(i,j)=sqrt(-2*log(A1(i,j)))*cos(8*atan(1.0)*A1(i,j))
enddo
enddo


do l=0,k-1
B1(0,l)=A2(0,l)/sqrt(A2(0,l)**2+A2(1,l)**2)
B1(1,l)=A2(1,l)/sqrt(A2(0,l)**2+A2(1,l)**2)
enddo

call tensor_pro(B1(0:1,0:0),B1(0:1,1:1),B1(0:1,2:2),B2,2,1,2,1,2,1)
r1=matmul(transpose(psi),B2)
r2=r1(0,0)
if (abs(r2)>0.9999) then
r3=r2
B3=B2
W1(0:1,0:0)=B1(0:1,0:0)
W1(0:1,1:1)=B1(0:1,1:1)
W1(0:1,2:2)=B1(0:1,2:2)
!if (r3>0.9999) then
!exit
!endif
exit
else
continue
endif
enddo
call gramschmidt(W1,G2)
call tensor_pro(G2(0:1,0:0),G2(0:1,1:1),G2(0:1,2:2),B4,2,1,2,1,2,1)
s1=matmul(transpose(B4),psi)
s11=s1(0,0)
do i=0,2**k-1
psi2(i,0)=r3*B3(i,0)+s11*B4(i,0)
enddo
print*,psi2
end program gdv

subroutine tensor_pro(T1,T2,T4,T5,m1,n1,m2,n2,m4,n4)
implicit none
integer::m1,n1,m2,n2,m3,n3,m4,n4
double precision,dimension(0:m1-1,0:n1-1)::T1
double precision,dimension(0:m2-1,0:n2-1)::T2
double precision,dimension(0:m1*m2-1,0:n1*n2-1)::T3
double precision,dimension(0:m4-1,0:n4-1)::T4
double precision,dimension(0:m1*m2*m4-1,0:n1*n2*n4-1)::T5
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

m3=m1*m2
n3=n1*n2

do i=0,m3-1,1
do j=0,n3-1,1
do k=0,m4-1,1
do l=0,n4-1,1
T5(i*m4+k,j*n4+l)=T3(i,j)*T4(k,l)
enddo
enddo
enddo    
enddo

    
end subroutine tensor_pro

subroutine gramschmidt(W0,W3)
implicit none
integer,parameter::k=3
double precision,dimension(0:1,0:k-1)::W0,W1,W2,v1,v2,y2,W3
integer::seed,i,j,l
double precision::r1
double precision,dimension(0:k-1,0:0)::k1
do i=0,k-1
    W1(0,i)=W0(0,i)/sqrt(W0(0,i)**2+W0(1,i)**2)
    W1(1,i)=W0(1,i)/sqrt(W0(0,i)**2+W0(1,i)**2)
enddo    
call random_seed(seed)
do i=0,k-1
   do j=0,1
      call random_number(v1(j,i))
   enddo
enddo

do i=0,k-1
   do j=0,1
     v2(j,i)=sqrt(-2*log(v1(j,i)))*cos(8*atan(1.0)*v1(j,i))
   enddo
enddo

do i=0,k-1
k1(i:i,0:0)=matmul(transpose(v2(0:1,i:i)),W1(0:1,i:i))
r1=k1(i,0)
!print*,r1
W2(0:1,i:i)=v2(0:1,i:i)-k1(i,0)*W1(0:1,i:i)
W3(0:0,i:i)=W2(0:0,i:i)/sqrt(W2(0:0,i:i)**2+W2(1:1,i:i)**2)
W3(1:1,i:i)=W2(1:1,i:i)/sqrt(W2(0:0,i:i)**2+W2(1:1,i:i)**2)
enddo


end subroutine gramschmidt