!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Auther: Varikuti Naga Dileep!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!program for finding reduced density matrices!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program multi_schmidt_decom
implicit none

integer,parameter::k=3,r1=1             !k is no of qubits
integer::i,j,m1,n1,m4,n4,l,INFO         ! if r1 =1 rhoA is caluculated,same for all
real,dimension(0:1,0:1)::p1,P,P2
real,dimension(0:2**k-1,0:0)::psi
real,dimension(0:2**k-1,0:2**k-1)::rhoABCDE
real,dimension(0:2**k-1,0:1)::F,F1
real,dimension(0:2**(r1-1)-1,0:2**(r1-1)-1)::B1
real,dimension(0:2**(k-r1)-1,0:2**(k-r1)-1)::B2
real,dimension(0:1,0:1)::id
real,dimension(2)::W
real,dimension(5)::WORK
id(0,0)=1
id(0,1)=0
id(1,0)=0
id(1,1)=1
psi(0:7,0)=(/0.0,1.0/sqrt(3.0),1.0/sqrt(3.0),0.0,1.0/sqrt(3.0),0.0,0.0,0.0/)
! psi is given state vector
rhoABCDE=matmul(psi,transpose(psi))!!!rh0ABCDE is density matrix
!write(111,*),rhoABCDE
do i=0,2**(r1-1)-1,1
do j=0,2**(r1-1)-1,1
if (i==j) then
B1(i,j)=1
else
B1(i,j)=0
endif
enddo
enddo

do i=0,2**(k-r1)-1,1
do j=0,2**(k-r1)-1,1
if (i==j) then
B2(i,j)=1
else
B2(i,j)=0
endif
enddo

enddo

do i=0,1,1
do j=0,1,1
P(i,j)=0
enddo

enddo

if (r1==1) then
do i=0,2**(k-r1)-1
call tensor_pro1(id,B2(:,i),F,2,2,2**(k-r1),1)
F1=matmul(rhoABCDE,F)
p1=matmul(transpose(F),F1)
P=P+p1
enddo
else
if (r1==k) then
do i=0,2**(r1-1)-1
call tensor_pro1(B1(:,i),id,F,2**(r1-1),1,2,2)
F1=matmul(rhoABCDE,F)
p1=matmul(transpose(F),F1)
P=P+p1
enddo
write(111,*),P

else
do i=0,2**(r1-1)-1
do j=0,2**(k-r1)-1
call tensor_pro(B1(:,i),id,B2(:,j),F,2**(r1-1),1,2,2,2**(k-r1),1)
F1=matmul(rhoABCDE,F)
p1=matmul(transpose(F),F1)
P=P+p1
enddo
enddo
endif
endif
P2=matmul(P,P)    !!! P2 is square of density matrix
write(33,*),P     !!!P is reduced density matrix

end program multi_schmidt_decom
subroutine tensor_pro(T1,T2,T4,T5,m1,n1,m2,n2,m4,n4)
implicit none
integer::m1,n1,m2,n2,m3,n3,m4,n4
real,dimension(0:m1-1,0:n1-1)::T1
real,dimension(0:m2-1,0:n2-1)::T2
real,dimension(0:m1*m2-1,0:n1*n2-1)::T3
real,dimension(0:m4-1,0:n4-1)::T4
real,dimension(0:m1*m2*m4-1,0:n1*n2*n4-1)::T5
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

   !!!subroutines for calculating tensor product 
end subroutine tensor_pro
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
!!! end of program!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
