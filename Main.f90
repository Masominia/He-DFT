program firstp
	implicit none
	real (8), allocatable :: rrad(:),urad(:),SU(:),HU(:),grad(:),ssp(:,:),hhp(:,:),rho(:),vh(:),vca(:)
	real (8) :: znuc,ener,sumn,gnorm,alpha,shti,eh,Etot,eta,Exc
	integer :: nrad,lang,iter,i
	complex (8), allocatable:: compg(:)
	
	alpha=1
	lang=0
	znuc=2.d0
	shti=0.1d0
	eta=1.d0
	
	read(*,*) nrad
	allocate (rrad(0:nrad))
	allocate (urad(0:nrad))
	allocate (SU(0:nrad))
	allocate (HU(0:nrad))
	allocate (grad(0:nrad))
	allocate (ssp(2,0:nrad))
	allocate (hhp(2,0:nrad))
	allocate (compg(0:nrad))
	allocate (rho(0:nrad))
	allocate (vh(0:nrad))
	allocate (vca(0:nrad))
	
	
	call radgrid(nrad,rrad)
	call Sci(nrad,rrad,urad)
	
	do iter=0,100000000
	
		call overlap (nrad,rrad,urad,sumn,SU)
		!write (*,*) sumn
		call normalize (nrad, urad, sqrt(sumn))
		call energr(nrad,lang,znuc,rrad,urad,ener,HU)
		call overlap (nrad,rrad,urad,sumn,SU)
		
		call ro(nrad,urad,rho)
		call Hartree(nrad,rrad,urad,rho,vh,Eh)
		call Excer(nrad,rrad,rho,eta,Exc,vca)
		Etot=ener+Eh+Exc
		
		
		call pGradiant (nrad,rrad,urad,ener,SU,HU,vh,vca,grad)
		!call Gradiant (nrad, ener, SU, HU, grad)
		!call Norm (nrad,grad,gnorm)
		call overlap(nrad,rrad,grad,gnorm,SU)
		gnorm=sqrt(gnorm)
		write (*,*) ener,Eh,exc,etot
		
		call crtssp(nrad, rrad, ssp)
		call crthhp (nrad, lang, znuc,rrad,hhp)
		
		call complexer (nrad,grad,compg)
		call ctridag (nrad+1, hhp,ssp,ener,shti,compg)
		
		if (gnorm<1.d-3) then
			exit
		end if
		
		
		do i=0,nrad
				urad(i)=urad(i)-alpha*real(compg(i))
		end do
		!write (*,*) iter,gnorm,ener
	end do
	
	write (*,*) iter,gnorm,ener
	
	call ro(nrad,urad,rho)
	call Hartree(nrad,rrad,urad,rho,vh,Eh)
	call Excer(nrad,rrad,rho,eta,Exc,vca)
	Etot=ener+Eh+Exc
	write (*,*) vh(0), Eh
	write (*,*) exc,Etot
	
end program
