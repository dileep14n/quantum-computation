program tensor_pro
implicit none
integer,parameter::k=4
double precision,dimension(0:1,0:k-1)::T1
double precision,dimension(0:2**k-1,0:0)::T2,T3,T5,T4
double precision,dimension(0:1,0:0)::A1,A2
integer::i,j,l,s

T1(0:1,0)=(/1,1/)
T1(0:1,1)=(/1,0/)
T1(0:1,2)=(/1,0/)
T1(0:1,3)=(/1,1/)

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

print*,T4

end program
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

