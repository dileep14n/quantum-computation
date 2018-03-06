program pp
implicit none
integer,parameter::k=3
integer::f,i,j
double precision,dimension(0:1,0:k-1)::G0,G1,W1,G2
W1(0:1,0)=(/1.0,0.0/)
w1(0:1,1)=(/1.0,0.0/)
W1(0:1,2)=(/1.0/sqrt(2.0),1.0/sqrt(2.0)/)
G2(0:1,0)=(/0.0,1.0/)
G2(0:1,1)=(/0.0,1.0/)
G2(0:1,2)=(/1.0/sqrt(2.0),-1.0/sqrt(2.0)/)
do f=k-1,0,-1
          do i=0,f,1
              G0(0:1,i:i)=W1(0:1,i:i)
              G1(0:1,i:i)=G2(0:1,i:i)
          enddo
          
          do j=f+1,k-1,1
              G0(0:1,j:j)=G2(0:1,j:j)
              G1(0:1,j:j)=W1(0:1,j:j)
          enddo
          print*,G0
          print*,G1
enddo
    
end program pp