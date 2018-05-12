program gdv
implicit none
integer,parameter::k=2
double precision,dimension(0:1,0:k-1)::A1,B1,A2,W1,W5,G1,G0,G2,U1,U2,W
double precision,dimension(0:2**k-1,0:0)::B2,psi,psi2,B3,B4,B5,A5,B6,A6,B7
double precision,dimension(0:0,0:0)::r1,s1,k2,v1,v2,E,E1
double precision,dimension(0:1,0:1)::sig_x,Id,Hd
double precision,dimension(0:3,0:3)::o1
double precision,dimension(0:7,0:7)::o2
double precision,dimension(0:1,0:1,0:k-1)::Aq
integer,dimension(0:k-1)::x
double precision::r2,r11,r3,k1,k3,s11,n11
integer::i,j,l,n,seed,count,m1
integer*16::m,f
!psi(0:7,0)=(/0.0/sqrt(2.0),1.0/sqrt(3.0),1.0/sqrt(3.0),0.0/sqrt(2.0),1.0/sqrt(3.0),0.0/sqrt(2.0),0.0/sqrt(2.0),0.0/sqrt(2.0)/)
!psi2(0:7,0)=(/0.0,0.0,0.0,0.40824829046386307,0.0,0.40824829046386307,0.81649658092772615,0.0/)
psi(0:3,0)=(/1.0/sqrt(7.0),2.0/sqrt(7.0),0.0,sqrt(2.0)/sqrt(7.0)/)
seed=32
call random_seed(seed)
r3=0
!do f=0,100,1
!n11=f/100.0d0
!do i=0,3
!do j=0,3
!o1(i,j)=0
!enddo
!enddo
!do i=0,7
!do j=1,7
!o2(i,j)=0
!enddo
!enddo
!do i=0,2**k-1
!psi(i,0)=0
!B2(i,0)=0
!B3(i,0)=0
!enddo
!do i=0,k-1
!do j=0,1
!A1(j,i)=0
!A2(j,i)=0
!B1(j,i)=0
!W1(j,i)=0
!G2(j,i)=0
!W(j,i)=0
!enddo
!enddo
!do i=0,1
!do j=0,1
!do l=0,k-1
!Aq(i,j,l)=0
!enddo
!enddo
!enddo
!r1(0,0)=0
!r3=0
!sig_x(0:1,0)=(/cos(n11*4*atan(1.0)),-sin(n11*4*atan(1.0))/)
!sig_x(0:1,1)=(/sin(n11*4*atan(1.0)),cos(n11*4*atan(1.0))/)
!!Hd(0:1,0)=(/1.0/2.0,-sqrt(3.0)/2.0/)
!Hd(0:1,1)=(/sqrt(3.0)/2.0,1.0/2.0/)
!Id(0:1,0)=(/1.0,0.0/)
!Id(0:1,1)=(/0.0,1.0/)
!call tensor_pro1(Id,sig_x,o1,2,2,2,2)
!call tensor_pro1(o1,Id,o2,4,4,2,2)
!psi=matmul(o2,psi2)

do i=0,k-1
x(i)=0  
enddo

do m=0,99999999     

do i=0,1,1
do j=0,k-1
call random_number(A1(i,j))
A2(i,j)=sqrt(-2*log(A1(i,j)))*cos(8*atan(1.0)*A1(i,j))         
enddo
enddo

do i=0,k-1
B1(0,i)=A2(0,i)/sqrt(A2(0,i)**2+A2(1,i)**2)   
B1(1,i)=A2(1,i)/sqrt(A2(0,i)**2+A2(1,i)**2)   
enddo
call tensor_pro(B1,B2)
r1=matmul(transpose(B2),psi)
r2=r1(0,0)
if (r2>r3) then
r3=r2
B3=B2

do i=0,k-1
W1(0:1,i:i)=B1(0:1,i:i)
enddo
else
continue
endif
enddo
call gramschmidt(W1,G2)               
do i=0,2**k-1
A5(i,0)=0
enddo

do i=0,k-1
Aq(0:1,0,i)=W1(0:1,i)
Aq(0:1,1,i)=G2(0:1,i)
enddo

do i=0,2**k-1
A5(i,0)=0
enddo

do j=0,2**k-1,1

do i=0,2**k-1,1
B5(i,0)=0   
enddo
call dm(j,x)
do i=0,k-1,1
W(:,i)=Aq(:,x(i),i)
enddo
call tensor_pro(W,B5)
v1=matmul(transpose(psi),B5)
do i=0,2**k-1
A5(i,0)=A5(i,0)+v1(0,0)*B5(i,0)   
enddo

!do i=0,1
if(v1(0,0)<0)then
write(12,*),-W(i,:),-v1
write(11,*),abs(v1)
else
write(12,*),W(i,:),v1
write(11,*),v1
endif     
!enddo

enddo

!enddo

!print*,A5
   
end program gdv

subroutine tensor_pro(T1,T4)
implicit none
integer,parameter::k=2
double precision,dimension(0:1,0:k-1)::T1
double precision,dimension(0:2**k-1,0:0)::T2,T3,T5,T4
double precision,dimension(0:1,0:0)::A1,A2
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
subroutine gramschmidt(W0,W3)
implicit none                                
integer,parameter::k=2
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
W2(0:1,i:i)=v2(0:1,i:i)-k1(i,0)*W1(0:1,i:i)
W3(0:0,i:i)=W2(0:0,i:i)/sqrt(W2(0:0,i:i)**2+W2(1:1,i:i)**2)
W3(1:1,i:i)=W2(1:1,i:i)/sqrt(W2(0:0,i:i)**2+W2(1:1,i:i)**2)
enddo
end subroutine gramschmidt
subroutine dm(i,x)
implicit none
integer,parameter::s=2
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
