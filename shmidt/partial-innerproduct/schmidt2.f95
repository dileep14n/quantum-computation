program schmidt
implicit none
integer,parameter::m=2,n=2,k=4,k1=5
real,dimension(0:m-1,0:n-1)::basis1,basis2,basis3,basis4,basis5
real,dimension(0:m-1,0:n*k+1)::BB
real,dimension(0:3,0:0)::t
real,dimension(0:2**k1-1,0:0)::state,state1,state2,statea,state1a,state2a
real,dimension(0:m-1,0:0)::state3a
integer::i,p,q,j,pa,qa,INFO
real,dimension(0:2**k1-1,0:0)::vec
real,dimension(0:2**k1-1,0:2**k1-1)::rhoABCDE
real,dimension(0:m-1,0:m-1)::rhoA,RhoB,rhoC,rhoD,rhoE
real,dimension(0:2**k1-1)::W
real,dimension(0:3*2**k1)::WORK
real,dimension(0:m-1)::eigen 
basis1=reshape((/1.0/sqrt(2.0),0.0,0.0,1.0/sqrt(2.0)/),(/2,2/))
basis2=reshape((/1.0,0.0,0.0,1.0/),(/2,2/))
basis3=reshape((/1.0/sqrt(2.0),1.0/sqrt(2.0),1.0/sqrt(2.0),-1.0/sqrt(2.0)/),(/2,2/))
basis4=reshape((/1.0/sqrt(2.0),1.0/sqrt(2.0),1.0/sqrt(2.0),-1.0/sqrt(2.0)/),(/2,2/))
basis5=reshape((/1.0,0.0,0.0,1.0/),(/2,2/))
BB(0:1,0:1)=basis1
BB(0:1,2:3)=basis2
BB(0:1,4:5)=basis3
BB(0:1,6:7)=basis4
BB(0:1,8:9)=basis5
do i=0,1
do j=0,9
!print*,BB(i,j)

enddo

enddo

do i=0,2**k1-1,1
state(i,0:0)=0
state1(i,0:0)=0
state2(i,0:0)=0

enddo

call tensor(BB(0:1,0:0),BB(0:1,2:2),state,2,2)
p=4
q=2
state1=state


do i=4,n*k,2
call tensor(state1,BB(0:1,i:i),state2,q,p)
state1=state2
p=p*2
enddo

!print*,state2



do i=0,2**k-1,1
statea(i,0:0)=0
state1a(i,0:0)=0
state2a(i,0:0)=0

enddo

call tensor(BB(0:1,1:1),BB(0:1,3:3),statea,2,2)
pa=4
qa=2
state1a=statea


do i=5,n*k+1,2
call tensor(state1a,BB(0:1,i:i),state2a,qa,pa)
state1a=state2a
pa=pa*2
!print*,state1
enddo

!print*,state2a

vec=state2a+state2

rhoABCDE=matmul(vec,transpose(vec))
!do i=0,2**k1-1,1
!write(34,*),rhoABCDE(i,:)
!enddo

call 	ssyev ("N", "U", 2**k1, rhoABCDE,2**k1, W, WORK,3*2**k1, INFO)
if (INFO==0) then
write(34,*),W
endif

state3a=BB(0:1,0:0)+BB(0:1,1:1)
call partial_trace(state3a,rhoA,2,eigen)
print*,eigen



end program schmidt
subroutine tensor(m1,m2,m3,m,n)
    implicit none
    integer::i,j,k,l
    integer :: m,n
    real,dimension(0:n-1,0:0)::m1
    real,dimension(0:m-1,0:0)::m2
    real,dimension(0:n*m-1,0:0)::m3
    real,dimension(0:2**4-1,0:2**4-1)::M6,m4



	do i=0,n-1,1
		do j=0,m-1,1
            !do k=0,n-1,1
            !    do l=0,m-1,1
            m3(i*m+j,0:0)=m1(i,0:0)*m2(j,0:0)
            !print*,m3(i*3+j,1)
        !enddo
        !enddo
		enddo
	enddo
    m4=transpose(m3)

    
end subroutine tensor

subroutine partial_trace(matrix1,matrix2,m,W)
    implicit none
    integer::m
    real,dimension(0:m-1,0:0)::matrix1
    real,dimension(0:m-1,0:m-1)::matrix2
    real,dimension(0:m-1)::W
    real,dimension(0:3*m-1)::WORK
    integer::INFO
    matrix2=matmul(matrix1,transpose(matrix1))
    call 	ssyev ("N", "U", m,matrix2,m, W, WORK,3*m, INFO)
    
end subroutine partial_trace