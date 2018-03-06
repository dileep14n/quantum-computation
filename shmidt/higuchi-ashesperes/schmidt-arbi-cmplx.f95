program gdv
implicit none
integer,parameter::k=3
complex*16,dimension(0:1,0:k-1)::B1,A2,W1,W5,G1,G0,G2,U1,U2,W,A12
double precision,dimension(0:1,0:k-1)::A1,A11
complex*16,dimension(0:2**k-1,0:0)::B2,psi,psi2,B3,B4,B5,A5,B6,A6,B7
complex*16,dimension(0:0,0:0)::r1,s1,k2,v1,v2,E,E1
complex*16,dimension(0:1,0:1)::sig_x,Id,Hd
complex*16,dimension(0:3,0:3)::o1
complex*16,dimension(0:7,0:7)::o2
complex*16,dimension(0:1,0:1,0:k-1)::Aq
complex*16,dimension(0:2**k-1,0:2**k-1)::rho,rho1
integer,dimension(0:k-1)::x
complex*16::r2,r11,r3,k1,k3,s11,n11
integer::i,j,l,n,seed,count,m1
integer*16::m,f
!!!!!!!!!!!!!!!!!!!!!initializing the k-qubit state vector

psi(0,0)=cmplx(0,0)
psi(1,0)=cmplx(1.0/sqrt(3.0),0)
psi(2,0)=cmplx(1.0/sqrt(3.0),0)
psi(3,0)=cmplx(0,0)
psi(4,0)=cmplx(1.0/sqrt(3.0),0)
psi(5,0)=cmplx(0,0)
psi(6,0)=cmplx(0,0)
psi(7,0)=cmplx(0,0)

seed=32
call random_seed(seed)
r3=cmplx(0,0)

do i=0,k-1
x(i)=0  
enddo

do m=0,99999                                               !loop for generating random numbers starts here

do i=0,1,1
do j=0,k-1
call random_number(A1(i,j))                                 !A1(0:2,0:k-1) stores random vector columns
call random_number(A11(i,j))
A12(i,j)=cmplx(A1(i,j),A11(i,j))
!A2(i,j)=sqrt(-2*log(A12(i,j)))*cos(8*atan(1.0)*A12(i,j))      ! box muller transformation of A1 elements are in A2
enddo
enddo

do i=0,k-1
B1(0,i)=A2(0,i)/sqrt((real(A2(0,i))**2+real(A2(1,i)**2))+aimag(A2(0,i)**2)+aimag(A2(1,i)**2))                 !columns of A2 are normalized and are in B1
B1(1,i)=A2(1,i)/sqrt((real(A2(0,i))**2+real(A2(1,i)**2))+aimag(A2(0,i)**2)+aimag(A2(1,i)**2))
enddo
call tensor_pro(B1,B2)                        ! tensor product of all unit columns are contained in B2(0:2**k-1,0:0)
r1=matmul(conjg(transpose(B2)),psi)                              !r1 is inner product 
r2=r1(0,0)
if(real(r2)**2>real(r3)**2.and.aimag(r2)==0.and.aimag(r3)==0)then                               ! condition to find unit vectors corresponding to maximized innerproduct 
r3=r2
print*,r2

!do i=0,k-1
!W1(0:1,i:i)=B1(0:1,i:i)                                    ! if condition is true, store all columns of B1 in W1
!enddo
else
continue
endif
enddo                                                     ! loop for generating random_numbers ends here
!call gramschmidt(W1,G2)                   ! columns of G2 contain orthonormal unit vectors of columns corresponding to W1    
!do i=0,2**k-1
!A5(i,0)=0
!enddo

!do i=0,k-1
!Aq(0:1,0,i)=W1(0:1,i)         ! logic to find all combinations of the unit vectors to find innerproducts and schmidt numbers
!Aq(0:1,1,i)=G2(0:1,i)         !Aq is a 3rd rank matrix contrains all unit vectors of the q qubit hilbert space.          
!enddo

!do i=0,2**k-1
!A5(i,0)=0
!enddo

!do j=0,2**k-1,1

!do i=0,2**k-1,1
!B5(i,0)=0   
!enddo
!call dm(j,x)
!do i=0,k-1,1
!W(:,i)=Aq(:,x(i),i)
!enddo

!call tensor_pro(W,B5)
!rho=matmul(B5,transpose(B5))
!do i=0,2**k-1
!do l=0,2**k-1
!rho1(i,l)=rho1(i,l)+rho(i,l)
!enddo
!enddo
!v1=matmul(transpose(psi),B5)
!do i=0,2**k-1
!A5(i,0)=A5(i,0)+v1(0,0)*B5(i,0)   
!enddo

!enddo
print*,r3

do i=0,2**k-1
write(123,*),rho1(:,i)
enddo
   
end program gdv

subroutine tensor_pro(T1,T4)
implicit none
integer,parameter::k=3
complex*16,dimension(0:1,0:k-1)::T1
complex*16,dimension(0:2**k-1,0:0)::T2,T3,T5,T4
complex*16,dimension(0:1,0:0)::A1,A2
integer::i,j,l,s

A1=T1(0:1,k-1:k-1)
A2=T1(0:1,k-2:k-2)
call tensor_pro1(A2,A1,T3,2,1,2,1)
T4=T3
s=4

do i=k-3,0,-1
call tensor_pro1(T1(0:1,i:i),T4,T5,2,1,s,1)
do j=0,2**k-1
T4(j,0)=0
enddo
T4=T5
s=s*2
do j=0,2**k-1
T5(j,0)=0
enddo
enddo
end subroutine tensor_pro
subroutine tensor_pro1(T1,T2,T3,m1,n1,m2,n2)
implicit none
integer::m1,n1,m2,n2
complex*16,dimension(0:m1-1,0:n1-1)::T1
complex*16,dimension(0:m2-1,0:n2-1)::T2
complex*16,dimension(0:m1*m2-1,0:n1*n2-1)::T3
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
subroutine gramschmidt(W0,W3)
implicit none                                
integer,parameter::k=3
complex*16,dimension(0:1,0:k-1)::W0,W1,W2,v2,y2,W3,v1
double precision,dimension(0:1,0:k-1)::v11,v22
integer::seed,i,j,l
complex*16::r1
complex*16,dimension(0:k-1,0:0)::k1
do i=0,k-1
W1(0,i)=W0(0,i)/sqrt(W0(0,i)**2+W0(1,i)**2)
W1(1,i)=W0(1,i)/sqrt(W0(0,i)**2+W0(1,i)**2)
enddo    
call random_seed(seed)
do i=0,k-1
do j=0,1
call random_number(v22(j,i))
call random_number(v11(j,i))
v1(i,j)=cmplx(v22(i,j),v11(i,j))
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
W2(0:1,i:i)=v2(0:1,i:i)-k1(i,0)*W1(0:1,i:i)
W3(0:0,i:i)=W2(0:0,i:i)/sqrt(W2(0:0,i:i)**2+W2(1:1,i:i)**2)
W3(1:1,i:i)=W2(1:1,i:i)/sqrt(W2(0:0,i:i)**2+W2(1:1,i:i)**2)
enddo
end subroutine gramschmidt
subroutine dm(i,x)
implicit none
integer,parameter::s=3
integer::i,a1,k,p
integer,dimension(0:s-1)::x

a1=i
p=0
do k = 0,s-1,1
if (mod(a1,2)==0) then
x(s-k-1) = 0
else
x(s-k-1) = 1
end if
a1 = a1/2
p=p+1
if (a1== 0) then
exit  
end if               
end do 
end subroutine
