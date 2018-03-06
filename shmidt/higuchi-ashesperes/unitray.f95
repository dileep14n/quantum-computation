program ppp
implicit none
integer,parameter::k=3
double precision,dimension(0:1,0:1)::sig_x,Id
double precision,dimension(0:3,0:3)::o1
double precision,dimension(0:7,0:7)::o2
double precision,dimension(0:2**k-1,0:0)::psi,psi2

!psi2(0:7,0)=(/0.0,1.0/sqrt(3.0),1.0/sqrt(3.0),0.0,1.0/sqrt(3.0),0.0,0.0,0.0/sqrt(8.0)/)
psi2(0:7,0)=(/1.0/sqrt(2.0),0.0/sqrt(4.0),0.0,0.0,0.0,0.0,0.0/sqrt(4.0),1.0/sqrt(2.0)/)

sig_x(0:1,0)=(/0.0,1.0/)
sig_x(0:1,1)=(/1.0,0.0/)
Id(0:1,0)=(/1.0,0.0/)
Id(0:1,1)=(/0.0,1.0/)
call tensor_pro1(Id,Id,o1,2,2,2,2)
call tensor_pro1(o1,sig_x,o2,4,4,2,2)
psi=matmul(o2,psi2)

print*,psi

    
end program ppp
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
