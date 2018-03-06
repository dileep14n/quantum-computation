program dm
implicit none
integer,parameter::s=3
integer::i,a1,k,p
integer,dimension(0:s-1)::x
      i=3

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
print*,x(0)
end program
