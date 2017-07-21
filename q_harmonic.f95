program Q_SHO
   implicit none
   integer::m,flag
   parameter(m=100)
   real*8,dimension(0:m-1,0:m-1)::a1,a_drag,X,P

   !a1 is lowering operator and a_drag is raising operator
   !X is position operator
   !P is momentum operator

   integer::i,j,INFO                                  
   !lapack variable declaration
   double precision::A(0:m-1,0:m-1),W(0:m-1),WORK(3*m-1),B(0:m-1,0:m-1),C(0:m-1,0:m-1),D(0:m-1,0:m-1),M1(0:m-1,0:m-1)
   double precision::D1(0:m-1,0:m-1)
   flag=0
   do i=0,m-1,1
    do j=0,m-1,1
       call elements(i,j,a1)
    
    enddo
   
   enddo
   !print*,a
   a_drag=transpose(a1)
   do i=0,m-1,1
     do j=0,m-1,1
        X(i,j)=a1(i,j)+a_drag(i,j)                   !h cross, mass and angular frequency are made to unity  
     
     enddo
   
   enddo
   do i=0,m-1,1
     do j=0,m-1,1
        A(i,j)=X(i,j)
     
     enddo
   
   enddo
   B=transpose(A)
   C=matmul(A,B)
   D=matmul(B,A)

   !checking for matrix A if it is normal or not

   do i=0,m-1,1
     do j=0,m-1,1
        if (C(i,j)/=D(i,j)) then
        	flag=1
        endif
     	
     enddo
   
   enddo
   if (flag==1) then
     print*,"matrix cant be diagonalized"
     stop
   else
   	print*,"matrix  can be diagonalized"
   	call dsyev ('V','U',m,A,m,W,WORK,3*m-1,INFO)
   if (INFO==0) then
   	open(UNIT=10,FILE='Q_SHO.txt',STATUS='REPLACE',ACTION='WRITE')
   	WRITE(10,*),"no error occured"
   	!WRITE(10,*)A
   else 
   	WRITE(10,*),"error occured"
   	
   endif
   endif
   M1=matmul(transpose(A),X)          
   D1=matmul(M1,A)                                     !D1 is the diagonalized matrix of position operator
   WRITE(10,*),'the diagonalized matrix is'
   do i=0,m-1,1
   	WRITE(10,30),D1(:,i)
   	30 format(100f14.8)
   	
   enddo
   WRITE(10,*)
   close(10)
 
 end program Q_SHO

!subroutine for calculating the elements of lowering operator a1

subroutine elements(nd,n,matrix)  
    implicit none
    integer::nd,n
    integer, parameter :: m=100
    real*8,dimension(0:m-1,0:m-1)::matrix
    if (nd==n-1) then
       matrix(nd,n)=n**0.5
       !print*,matrix(nd,n)
   else
       matrix(nd,n)=0
    endif
    
end subroutine elements