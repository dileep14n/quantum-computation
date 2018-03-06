program gdv
implicit none
integer,parameter::k=3
double precision,dimension(0:1,0:1,0:k-1)::A
double precision,dimension(0:1,0:k-1)::A1,B1,A2,W1,W5,G1,G0,G2,U1,U2,W
double precision,dimension(0:2**k-1,0:0)::B2,psi,psi2,B3,B4,B5,A5,B6,A6,B7
double precision,dimension(0:0,0:0)::r1,s1,k2,v1,v2,E,E1
integer,dimension(0:k-1)::x1
double precision::r2,r11,r3,k1,k3,s11
integer*4::i,j,l,n,seed,count,f,m1
integer*16::m
-g-Wall-Wextra-Warray-temporaries-Wconversion-fimplicit-none-fbacktrace-ffree-line-length-0-fcheck=all-ffpe-trap=zero,overflow,underflow-finit-real=nan

!psi(6:15,0)=(/1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0/)
!psi(0:5,0)=(/1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/4.0/)
psi(0:7,0)=(/0.0,1.0/sqrt(3.0),1.0/sqrt(3.0),0.0,1.0/sqrt(3.0),0.0,0.0,0.0/sqrt(8.0)/)
!psi(0:7,0)=(/1.0/sqrt(3.0),1.0/sqrt(3.0),1.0/sqrt(3.0),0.0,0.0/sqrt(3.0),0.0,0.0,0.0/sqrt(3.0)/)
!psi2(0:7,0)=(/1.0/sqrt(2.0),0.0/sqrt(4.0),0.0,0.0,0.0,0.0,0.0/sqrt(4.0),1.0/sqrt(2.0)/)
seed=32


call random_seed(seed)
r3=0

do i=0,1
  do j=0,1
    do l=0,k-1
       A(i,j,l)=0      
    enddo
  enddo
enddo

do m=0,999999     

        do i=0,1,1
        do j=0,k-1,1
           call random_number(A1(i,j))
           A2(i,j)=sqrt(-2*log(A1(i,j)))*cos(8*atan(1.0)*A1(i,j))         
        enddo
        enddo


        do i=0,k-1,1
           B1(0,i)=A2(0,i)/sqrt(A2(0,i)**2+A2(1,i)**2)   
           B1(1,i)=A2(1,i)/sqrt(A2(0,i)**2+A2(1,i)**2)   
        enddo

    call tensor_pro(B1,B2)
    r1=matmul(transpose(B2),psi)
    r2=r1(0,0)
    if (r2>r3) then
       r3=r2
       write(36,*),B2
       B3=B2

       do i=0,k-1,1
          W1(0:1,i:i)=B1(0:1,i:i)
       enddo
    
    else
       continue
    endif

enddo

 call gramschmidt(W1,G2)                 

print*,W1,G2,r3

do i=0,1000
E(0,0)=i/10000d0
do j=0,k-1
do l=0,1
U1(l,j)=W1(l,j)+E(0,0)*G2(l,j)      
enddo
enddo
call tensor_pro(U1,B7)
v2=matmul(transpose(psi),B7)
write(335,*),v2,E
enddo    


 
    
    
end program gdv


! subroutine for calculating tensor product for n column vectors

subroutine tensor_pro(T1,T4)
implicit none
integer,parameter::k=3
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


subroutine gramschmidt(W0,W3)                 ! W0 is a matrix contains unit column vectors
implicit none                                 ! W3 contains unit column vectors that are orthonormal to those in W0
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
