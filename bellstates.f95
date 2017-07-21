program bell_states
   implicit none
   real,dimension(2,1)::zero,one,zero1,one1
   real,dimension(2,2)::hadamard,h
   real,dimension(4,4)::d1,d2,d3,d4,d,de
   integer,dimension(4,4)::cnot

   !zTz=tensor_product of zero and zero
   !zTo=tensor_product of zero and one
   !oTz=tensor_product of one and zero
   !oTo=tensor_product of one and one
   real,dimension(4,1)::zTz,zTo,oTz,oTo,b1,b2,b3,b4
   integer::i,j,m
   m=0

   !input data....
   zero=RESHAPE((/1,0/),(/2,1/))
   one=RESHAPE((/0,1/),(/2,1/))
   h=RESHAPE((/1.0,1.0,1.0,-1.0/),(/2,2/))
   do i=1,2,1
   	do j=1,2,1
   		hadamard(i,j)=(1.0/sqrt(2.0))*h(i,j)
   		
   	enddo
   	
   enddo

   cnot=RESHAPE((/1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0/),(/4,4/))

   !hadamard transform of 1st qubit
   zero1=matmul(hadamard,zero)
   !print*,zero1
   one1=matmul(hadamard,one)
   !print*,one1
   call tensor_product(zero1,zero,zTz)
   call tensor_product(zero1,one,zTo)
   call tensor_product(one1,zero,oTz)
   call tensor_product(one1,one,oTo)
   
   !cnot operation 
   b1=matmul(cnot,zTz)
   b2=matmul(cnot,zTo)
   b3=matmul(cnot,oTz)
   b4=matmul(cnot,oTo)


   !density matrix corresponding to pure states
   d1=matmul(b1,transpose(b1))
   d2=matmul(b2,transpose(b2))
   d3=matmul(b3,transpose(b3))
   d4=matmul(b4,transpose(b4))

   !density matrix corresponding to mixed states 

    d=d1+d2+d3+d4
    do i=1,4
      do j=1,4
        de(i,j)=(0.25)*d(i,j)
        
      enddo
      
    enddo

    open(unit=15,file='bells.txt',status='replace',action='write')
    write(15,*),'bell states are'
    write(15,*),'bell state 1 is'
    do i=1,4,1
    write(15,20),b1(i,:)
    20 format(10f8.4)
  enddo
  write(15,*),'bell state 2 is'
  do i=1,4,1
    write(15,21),b2(i,:)
    21 format(10f8.4)
  enddo
  write(15,*),'bell state 3 is'
do i=1,4,1
    write(15,22),b3(i,:)
    22 format(10f8.4)
  enddo
  write(15,*),'bell state 4 is'
do i=1,4,1
    write(15,23),b4(i,:)
    23 format(10f8.4)
  enddo
  write(15,*)'density matrix 1 is'
  do i=1,4,1
    write(15,24),d1(:,i)
    24 format(20f10.4)
    
  enddo
  write(15,*)'density matrix 2 is'
  do i=1,4,1
    write(15,25),d2(:,i)
    25 format(20f10.4)
    
  enddo
  write(15,*)'density matrix 3 is'
  do i=1,4,1
    write(15,26),d3(:,i)
    26 format(20f10.4)
    
  enddo
  write(15,*)'density matrix 4 is'
  do i=1,4,1
    write(15,24),d4(:,i)
    27 format(20f10.4)
    
  enddo
  write(15,*)'density matrix for mixed states is'
  do i=1,4,1
    write(15,28),de(:,i)
    28 format(20f10.6)
    
  enddo


    
end program bell_states
subroutine tensor_product(m1,m2,m3)
    implicit none
    integer::k,l
    integer, parameter :: n=2
    real,dimension(2,1)::m1,m2
    real,dimension(4,1)::m3
    do k=0,n-1,1
    	do l=1,n,1
    		m3(k*n+l,1)=m1(k+n-1,1)*m2(l,1)
    		
    	enddo
    	
    enddo
    
end subroutine tensor_product