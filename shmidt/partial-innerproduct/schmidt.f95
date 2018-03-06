program multi_partite
implicit none
integer, parameter :: ma=2,na=1,k=4
real,dimension(0:ma-1,0:na-1)::KET0,KET1
real,dimension(0:ma**k-1,0:na-1)::psi_i
real,dimension(0:ma-1,0:ma-1)::basis1,basis2,basis3,basis4,basis5
real,dimension(0:ma-1,0:k*ma+1)::BB
real,dimension(0:2**4-1,0:2**4-1)::state
integer::i,j,m,n
KET0(0,0)=1
KET0(1,0)=0
KET1(0,0)=0
KET1(1,0)=1
basis1=reshape((/1,0,0,1/),(/2,2/))
basis2=reshape((/1,0,0,1/),(/2,2/))
basis3=reshape((/1,1,1,-1/),(/2,2/))
basis4=reshape((/1,1,1,-1/),(/2,2/))
basis5=reshape((/1,0,0,1/),(/2,2/))
BB(0:1,0:1)=basis1
BB(0:1,2:3)=basis2
BB(0:1,4:5)=basis3
BB(0:1,6:7)=basis4
BB(0:1,8:9)=basis5
do i=0,ma**k-1,1
do j=0,na**k-1,1
state(i,j)=0
enddo
enddo
m=2
n=1
do i=0,ma-1
do j=0,0
state(i,j)=BB(i,j)
enddo
enddo
do i=1,k-1,2
call tensor(BB(0:1,i:i),BB(0:1,i+1:i+1),state,m,n)
print*,state
m=m*2
enddo

    
end program multi_partite
subroutine tensor(m1,m2,m4,m,n)
    implicit none
    integer::i,j,k,l
    integer :: m,n
    real,dimension(0:n-1,0:n-1)::m1
    real,dimension(0:m-1,0:m-1)::m2
    real,dimension(0:n*m-1,0:n*m-1)::m3
    real,dimension(0:1,0:3)::m5
    real,dimension(0:2**4-1,0:2**4-1)::M6,m4
    m2=reshape((/1,0,0,-2/),(/2,2/))
    m1=reshape((/1,0,3,2,4,0,1,0,1/),(/3,3/))
    m5(0:1,0:1)=transpose(m2)
    m5(0:1,2:3)=transpose(m2)


	do i=0,n-1,1
		do j=0,m-1,1
            do k=0,n-1,1
                do l=0,m-1,1
            M6(i*m+j,k*m+l)=m1(i,k)*m2(l,j)
            !print*,m3(i*3+j,1)
        enddo
        enddo
		enddo
	enddo
    m4=transpose(m3)

    
end subroutine tensor