subroutine Excer(nrad,rrad,rho,eta,Exc,vca)
	
	implicit none
	integer, intent (in) :: nrad
	real (8), intent(in):: rho(0:nrad),eta,rrad(0:nrad)
	real(8), intent(out):: Exc,vca(0:nrad)
	real(8):: ec(0:nrad),vcb,PI=4.d0*atan(1.d0)
	integer:: i
	
	exc=0
	
	do i=0, nrad
		if (rho(i)>1.d-20) then
			call LSD_PADE(rho(i),eta,ec(i),vca(i),vcb)
		else
			vca(i)=0
			ec(i)=0
			vcb=0
		end if
	end do
	
	do i=0,nrad-1
		exc=exc+ ((rrad(i+1)-rrad(i))/2) * ( rrad(i)**2*ec(i)*rho(i) +  rrad(i+1)**2*ec(i+1)*rho(i+1) )
	end do
		exc=exc*4.d0*pi
end subroutine Excer

subroutine Hartree (nrad,rrad,urad, rho, Vh, Eh)
	implicit none
	integer, intent (in) :: nrad
	real (8), intent (in) :: rho (0:nrad), rrad(0:nrad),urad(0:nrad)
	real (8), intent (out) :: Vh(0:nrad),Eh
	integer :: i,j
	real (8) :: PI=4.d0*atan(1.d0)
	
	!Hartree Potential *************************
	vh=0
	vh(0)=0

	do j=1,nrad-1
		vh(0)=vh(0)+((rrad(j)*rho(j)+rrad(j+1)*rho(j+1))) * (rrad(j+1)-rrad(j))/2
	end do
	vh(0)=(4d0*pi)*vh(0)
	
	do i=1,nrad
		
		do j=0,i-1
			vh(i)=vh(i)+ (1/rrad(i)) * (rho(j)*rrad(j)**2+(rho(j+1)*rrad(j+1)**2))*(rrad(j+1)-rrad(j))/2
		end do
		do j=i,nrad-1
				vh(i)=vh(i)+((rrad(j)*rho(j)+Rrad(j+1)*rho(j+1)))* (rrad(j+1)-rrad(j))/2
		end do
		
		vh (i)= vh(i)*(4d0*pi)
		!write (*,*) vh(i)
	end do
	
	
	!Hartree Energy *************************
	Eh=0
	do i=0,nrad-1
		Eh=Eh + (4*pi*(rrad(i)**2*rho(i)*vh(i)+ rrad(i+1)**2*rho(i+1)*vh(i+1))/2)*(rrad(i+1)-rrad(i))/2!(4*pi* rrad(i)**2)*(rrad(i+1)-rrad(i-1))/2
	end do
	

end subroutine Hartree


subroutine ro(nrad, urad,rho)
	implicit none
	integer, intent (in) :: nrad
	real (8), intent (in) :: urad(0:nrad)
	real (8), intent (out) :: rho(0:nrad)
	integer :: i
	real (8) :: PI=4d0*atan(1.d0)
	
	do i=0,nrad
		rho(i)=1*(urad(i)*urad(i))/(4d0*pi)
	end do
	
end subroutine ro

subroutine complexer (nrad,grad,fin)  
	
	implicit none 
	integer, intent(in) :: nrad
	real (8), intent(in) :: grad(0:nrad)
	complex (8), intent (out) :: fin (0:nrad)
	integer :: i

	do i=0,nrad
		fin(i)=cmplx(grad(i),0.d0)
	end do
	
end subroutine complexer



subroutine Sci (nrad,rrad,urad)  
	
	implicit none 
	integer, intent (in) :: nrad

	real (8), intent(out):: urad(0:nrad)
	real(8), intent (in) :: rrad(0:nrad)
	integer ::	 i
	!allocate (urad(0:nrad))
	!allocate (rrad(0:nrad))

	do i=0,nrad
		urad(i)=2.d0*exp(-1.2d0*rrad(i))
	end do
end subroutine Sci




subroutine Normalize (n,vec,coe) 
	implicit none
	integer, intent (in) :: n
	real (8), intent (inout) :: vec(0:n)
	real (8), intent (in) :: coe
	integer :: i	
	do i=0,n
		vec(i)=vec(i)/coe
	end do
	
end subroutine Normalize



subroutine pGradiant (nrad,rrad,urad, ener, SU, HU,vh,vca, grad)
	implicit none
	integer, intent (in)::nrad
	real (8), intent (in) :: SU(0:nrad), HU(0:nrad),vca(0:nrad),vh(0:nrad),rrad(0:nrad),urad(0:nrad),ener
	real(8), intent (out):: grad(0:nrad)
	real (8):: d(0:nrad),e
	integer :: i
	
	i=nrad
	d(i)=( (vh(i)+vca(i))*rrad(i)**2* ((rrad(i)-rrad(i-1))/2) * urad(i) )
	i=0
	d(i)=( (vh(i)+vca(i))*rrad(i)**2* ((rrad(i+1)-rrad(i))/2) * urad(i) )

	do i=1,nrad-1
		d(i)=( (vh(i)+vca(i))*rrad(i)**2* ((rrad(i+1)-rrad(i-1))/2) * urad(i) )
	end do
	
	e=ener
	do i=0,nrad
		e=e+d(i)*urad(i)
	end do
	
	grad=0;
	do i=0,nrad
		grad(i)= HU(i) + d(i) - e*SU(i)
	end do
		
end subroutine pGradiant

subroutine Gradiant (nrad, e, SU, HU, grad)
	implicit none
	integer, intent (in)::nrad

	real (8), intent (in) :: SU(0:nrad), HU(0:nrad)
	real(8), intent(in) :: e
	real(8), intent (out):: grad(0:nrad)
	integer :: i

	do i=0,nrad
		grad(i)=HU(i)-e*SU(i)
	end do

		
end subroutine Gradiant




subroutine Norm (n, vec, ans)
	implicit none
	integer, intent(in) :: n
	real (8), intent (in) :: vec(0:n)
	real(8), intent (out):: ans
	integer :: i
	
	ans=0
	
	do i=0,n
		ans=ans+vec(i)*vec(i)
	end do
	
	ans=sqrt(ans)

end subroutine norm
