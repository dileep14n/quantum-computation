program multi_schmidt_decom
implicit none

integer,parameter::k=3,r1=2                                !k is no of qubits
integer::i,j,m1,n1,m4,n4,l,INFO,s                            ! if r1 =1 rhoA is caluculated and is stored in P, same for rhoB(r1=2), rhoc(r1=3) so on..
real,dimension(0:1,0:1)::p1,P                              !psi is any arbitrary column vector
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
!psi(0:11,0)=(/1.0/sqrt(8.0),1.0/sqrt(8.0),0.0,0.0,1.0/sqrt(8.0),1.0/sqrt(8.0),0.0,0.0,0.0,0.0,0.0,0.0/)
!psi(12:25,0)=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
!psi(26:31,0)=(/1.0/sqrt(8.0),-1.0/sqrt(8.0),0.0,0.0,-1.0/sqrt(8.0),1.0/sqrt(8.0)/)
psi(0:7,0)=(/0.0,1.0/sqrt(3.0),1.0/sqrt(3.0),0.0/sqrt(3.0),1.0/sqrt(3.0),0.0/sqrt(3.0),0.0/sqrt(3.0),0.0/)
rhoABCDE=matmul(psi,transpose(psi))

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

!do s=1,r1,1
    do i=0,1,1      
        do j=0,1,1
            P(i,j)=0
            p1(i,j)=0
        enddo
    enddo

    if (s==1) then
        do i=0,2**(k-s)-1
            call tensor_pro1(id,B2(:,i),F,2,2,2**(k-s),1)
            F1=matmul(rhoABCDE,F)
            p1=matmul(transpose(F),F1)
            P=P+p1
            write(111,*),p1
        enddo
    else
        if (s==k) then
            do i=0,2**(s-1)-1
                call tensor_pro1(B1(:,i),id,F,2**(s-1),1,2,2)
                F1=matmul(rhoABCDE,F)
                p1=matmul(transpose(F),F1)
                P=P+p1
                write(111,*),p1
            enddo
        else
            do i=0,2**(s-1)-1
                do j=0,2**(s-r1)-1
                    call tensor_pro(B1(:,i),id,B2(:,j),F,2**(s-1),1,2,2,2**(k-s),1)
                    F1=matmul(rhoABCDE,F)
                    p1=matmul(transpose(F),F1)
                    P=P+p1
                enddo
            enddo
            print*,p1
        endif
    endif
    call    ssyev ("N", "U", 2,P,2, W, WORK,5, INFO)
    if (INFO==0) then   
        write(33,*),W
    endif


!enddo
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
