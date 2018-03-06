program dm
implicit none
integer,parameter::s=6
integer::i,j,k,l,m,w1,w2,a1,b,h,q,p,u,e1,e2,d1,d2,g1,g2,e3,e4,d3,d4,g3,g4
integer::xx,yy,zz,c1=1,c2=1,c3=1,c4=1,j1,j2,j3
real*4::j11,j22,j33
real*4,dimension(0:2**s-1,0:2**s-1)::HDM,HA1,HA2,HA4,HA3,H1,H2
real*4,dimension(0:s-1)::x,y,z,t
complex,dimension(0:2**s-1,0:2**s-1)::Hadm,HHDM,HAM
CHARACTER*1 ::JOBVL, JOBVR
INTEGER::INFO, LDA, LDVL, LDVR, LWORK, N
parameter(N=2**s,LDA=N,LDVL=N,LDVR=N,LWORK=2*N)
COMPLEX*16::A(0:LDA-1,0:N-1), VL(0:LDVL-1,0:N-1), VR(0:LDVR-1,0:N-1),W(0:N-1), WORK(0:0,0:LWORK-1),c(0:1,0:1)
DOUBLE PRECISION::RWORK(0:2*N-1)
JOBVL='V'
JOBVR='V'
! starting of loop of hamiltonian
do i=0,2**s-1,1
do j=0,2**s-1,1
HDM(i,j)=0
H1(i,j)=0
H2(i,j)=0
end do
end do


                do i=0,2**s-1,1
	        do j=0,2**s-1,1
                 !i=4
                 !j=2
                   do u=0,s-1
                   x(u)=0  
                   y(u)=0   !intialize the matrices which contains binary equivalent
                   end do
                         
    	                     !convert i to binary
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
                                 ! convert j to binary
                                 b=j
				do l = 0,s-1,1
                                 q=0
				if (mod(b,2)==0) then
				   y(s-l-1) = 0
				else
				   y(s-l-1) = 1
				end if
				b = b/2
				 q=q+1

				if (b == 0) then
				   exit  
				end if
		                end do

t=y 
e1=0
d1=0
e2=0
d2=0 
g1=0
g2=0
e3=0
d3=0
e4=0
d4=0 
g3=0
g4=0

do m=0,s-1,1

if(m.le.s-2) then
  if( (t(m)==1) .or. (t(m+1)==0)) then
       e1=0
  else
      t(m)=1
      t(m+1)=0
   call btod(t,w1)
      if(w1==i)  then
         d1=1
      else
         d1=0
      end if
      g1=g1+d1
  end if
else
   if( (t(s-1)==1) .or. (t(0)==0)) then
       e2=0
  else
      t(s-1)=1
      t(0)=0
   call btod(t,w1)
      if(w1==i)  then
         d2=1
      else
         d2=0
      end if
      g2=g2+d2
  end if     
end if
t=y
end do


H1(i,j)=e1+e2+g1+g2

do m=0,s-1,1

if(m.le.s-2) then
  if( (t(m)==0) .or. (t(m+1)==1)) then
       e3=0
  else
      t(m)=0
      t(m+1)=1
   call btod(t,w2)
      if(w2==i)  then
         d3=1
      else
         d3=0
      end if
      g3=g3+d3
  end if
else
   if( (t(s-1)==0) .or. (t(0)==1)) then
       e4=0
  else
      t(s-1)=0
      t(0)=1
   call btod(t,w2)
      if(w2==i)  then
         d4=1
      else
         d4=0
      end if
      g4=g4+d4
  end if     
end if
t=y
end do


H2(i,j)=e3+e4+g3+g4

end do
end do
HDM=H2-H1
do i=0,2**s-1,1
write(684,*),HDM(i,:)
end do
do i=0,2**s-1,1     
do j=0,2**s-1,1
Hadm(i,j)=cmplx(0,HDM(i,j)/2)
end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! nearest neighbour interaction!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

xx=1
yy=2
zz=3
call interaction(xx,HA1)
call interaction(yy,HA2)
call interaction(zz,HA3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!!!!! over!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !do j1=0,10,1
  j11=j1/10.0d0
 ! do j2=0,10,1
  j22=j2/10.0d0
  do j3=0,10,1
  j33=j3/10.0d0
  HAM=c1*HA1+0.5*c2*HA2+j33*c3*Hadm
A=HAM
call ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
             if (INFO==0) then
    	        print*,"no error occured"
                
                  do i=0,2**s-1,1    
                    write(685,*),j33,real(W(i))
                  end do
             endif
             if (INFO/=0) then
    	        print*,"error occured"
                stop
             endif
!end do
!end do
end do


end program


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! subroutine for nearest neighbour!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine interaction(xx,H)
integer,parameter::s=6
real*4,dimension(0:2**s-1,0:2**s-1)::H,HH
real*4,dimension(0:s-1)::x,y,z,t
integer::i,j,k,l,temp,ww,o,temp2,q,a1,b,p,u,j1,xx,yy,zz
real*4::h1,h2,h3,h4,h5,h6,j11,h7=0
 do i=0,2**s-1,1
	        do j=0,2**s-1,1
                    H(i,j)=0
                end do
                end do
! starting of loop of hamiltonian
                do i=0,2**s-1,1
	        do j=0,2**s-1,1
                  !i=1
                  !j=2
                   do u=0,s-1
                   x(u)=0  
                   y(u)=0   !intialize the matrices which contains binary equivalent
                   end do
    	                     !convert i to binary
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
                                 ! convert j to binary
                                 b=j
				do l = 0,s-1,1
                                 q=0
				if (mod(b,2)==0) then
				   y(s-l-1) = 0
				else
				   y(s-l-1) = 1
				end if
				b = b/2
				 q=q+1

				if (b == 0) then
				   exit  
				end if
				end do 

h4=0
h2=0
h3=0
h7=0
if(i==j)then
h5=1
else 
h5=0
end if
 !7.92
do o=0,s-1,1
  
    z=y
if(o.le.s-xx-1) then  
        temp=z(o)
        z(o)=z(o+xx)    
        z(o+xx)=temp
  call btod(z,ww)
      if(ww==i) then
         h1=1
      else
         h1=0
      end if
  h2=h2+h1

else
         z=y
         zz=0
        !print*,o,o+xx-s  
       temp=z(o)
       z(o)=z(o+xx-s)
       z(o+xx-s)=temp
  call btod(z,q)
        !print*,y
      if(q==i) then
         h7=1
      else
         h7=0
      end if
       
      h3=h3+h7     
end if

end do 
h4=h2+h3

H(i,j)=+h4/2-h5*s/4

end do
end do
end subroutine
! subroutine to convert binary to decimal
subroutine btod(m,n)              
integer*4,parameter::s=6
real*4,dimension(0:s-1)::m
integer::n,i,k
n=0
do i=0,s-1
  k=(2**(i))*(m(s-1-i))
n=n+k
end do
end subroutine btod










