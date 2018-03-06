program comb
implicit none
integer,parameter::k=3
integer,dimension(0:k-1)::x
double precision,dimension(0:1,0:1,0:k-1)::A
double precision,dimension(0:1,0:k-1)::u,v,W
integer::i,j,l,m

u(0:1,0)=(/1,0/)
u(0:1,1)=(/1,0/)
u(0:1,2)=(/1,1/)
v(0:1,0)=(/0,1/)
v(0:1,1)=(/0,1/)
v(0:1,2)=(/1,-1/)


do i=0,k-1
  A(0:1,0,i)=u(0:1,i)
  A(0:1,1,i)=v(0:1,i)
enddo

do j=0,2**k-1
do i=0,k-1
  call dm(j,x)
  W(0:1,i)=A(0:1,x(i),i)
enddo

print*,W
do l=0,k-1
  do m=0,1
    W(m,l)=0

  enddo
enddo

enddo

end program

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
