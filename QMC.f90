program QMC  		
  implicit none 
  integer, parameter :: N=400 ! number of walkers
  integer, parameter :: MC=30000! number of MC steps 
  integer, parameter :: Itime=4000 ! Initialization steps 
  integer, parameter :: intbeta=5 ! number of different alpha values 
  real(8) :: walkers(6),testwalker(6),alpha,a, beta, Etest, Etotal     !do one experiment at a time, just walkers, testwalker, Etest(for current Elocal) and Etotal 
  real(8) :: R_L(3), R_R(3)         ! Position of nuclei  
  integer :: i,k,j 
  real(8) :: Emean(intbeta,N)
  real(8) :: s 

  open(UNIT=12, file="Emean(k-i).dat", status="replace")

  beta=0.3         !set first beta value 
  s=3d0            !Set distance between nuclei (test, incorporate later) +
  R_L=(/-0.5*s,0d0,0d0/)
  R_R=(/+0.5*s,0d0,0d0/) ! Place nuclei symmetrically around zero, and along x-axis (only distance matters)
  alpha=2            ! From cusp condition
  call NewtonRaphson(a,s)  ! Determine parameter with cusp condition (only need to do this once)
  Emean=0
!-------------------------------------------------------------------------------------------------------
!-------------------- MAIN PROGRAM AND SIMULATION-------------------------------------------------------

  do k=1,intbeta				!Parameter loop: \beta from 0.3 to 0.3+0.1*intbeta	
    call init_random_seed()
    call random_number(walkers)                 !First positions you choose "randomly" --> Get close to minimum
    walkers=walkers-0.5					
    write(12,*) beta				!to tell where are we, or we will lost ourselves in the sea of data...:D
    !for N independent Experiments
    do i=1, N
      !Initialization Loop
      do j=1,Itime	
        call Randwalk(walkers,testwalker)
        call Metropolis_ini(walkers,testwalker,alpha,beta,a,R_L,R_R)   !here is why I think we should use "Module"
      end do 
      Etest=localE(walkers,alpha,beta,a,R_L,R_R)		!get the first Elotal in case metropolis test fails for the first couple steps
      Etotal=0
      !MC Simulation
      do j=1,MC	
	!make loop for walkers  	  			!Start Monte Carlo simulation loop
	call Randwalk(walkers,testwalker)		  	! Let every particle take a small random 3D step 
	call Metropolis(walkers,testwalker,alpha,beta,a,R_L,R_R,Etest)    	
	Etotal=Etotal+Etest 					!accumulate Elocal for every step and get the average later
      end do 	

      !call testE(Elocal,MC,N,k,Etest)      !Etest(N,intbeta) will be mean Elocal for every walker for every beta
      Emean(k,i)=Etotal/MC                  !Use average over N independent walkers 
      write(12,*) Emean(k,i)		    !write the result of every experiment
      !print *, beta, Emean(k)
    end do

    beta=beta+0.1	             ! Choose different parameter 
  end do 

  close(12)
!--------------------------------------------------------------------------------------------------------
!-------------------POST PROCESSING, PRINT STATEMENTS ---------------------------------------------------
!do i=1,N
!print *, Etest(i,1), Elocal(i) , MC
!end do
 !  do i=1,2
  !  print *, alpha, Etest(:,i)
  ! end do
  
contains 
! -------------------------------------------------------------------------------------------------------
!--------------------------------SUBROUTINES AND FUNCTIONS----------------------------------------------- 
!--------------------------------------------------------------------------------------------------------

!Take a small step forward or backwards, for every particle. select x,y and z positions 
subroutine Randwalk(walkers,testwalker)
  real(8), intent(out):: testwalker(:)
  real(8), intent(inout) :: walkers(:)        
  real(8) :: A1(6) 
  call random_number(A1) 
  testwalker=walkers+1.0*(A1-0.5)
end subroutine Randwalk

!-------------------------------------------------------------------------------------------------------------
!Final Return :: Square of wavefunction |Psi|^2, a N-dimensional vector (for every walker)
real(8) function orbit(r,R_L,R_R,a)
  real(8), intent(in) :: r(:)
  real(8), intent(in) :: R_L(:), R_R(:)
  real(8), intent(in) :: a
  real(8) :: X_L(3), X_R(3)
    X_L=r-R_L       ! Matrix (3,N) with "distance" (difference only) to Left nucleus for every walker
    X_R=r-R_R    
  orbit=exp(-sqrt(sum(X_L**2))/a) + exp(-sqrt(sum(X_R**2))/a)
  return
end function orbit

real(8) function Jastrow(r1,r2,alpha,beta)
  real(8), intent(in) :: r1(:), r2(:)
  real(8), intent(in) :: alpha, beta 
  real(8) :: dr(3)
  dr=r1-r2
  Jastrow=exp(sqrt(sum(dr**2))/(alpha*(1+beta*sqrt(sum((dr**2))))))
  return
end function Jastrow

real(8) function Wavefunction(walkers,alpha,beta,R_L,R_R,a)
  real(8), intent(in) :: alpha, beta, a, R_L(:), R_R(:)
  real(8), intent(in) :: walkers(:)
  real(8) :: r1(3), r2(3)
  r1=walkers(1:3)
  r2=walkers(4:6)
  Wavefunction=orbit(r1,R_L,R_R,a)*orbit(r2,R_L,R_R,a)*Jastrow(r1,r2,alpha,beta)
  Wavefunction=Wavefunction**2
  return
end function Wavefunction

!---------------------------------------------------------------------------------------------------------------
!Metropolis with |psi|^2 probability   (NOTE Walkers inside the subroutine is a vector! (We input walkers==walkers(:,j)
subroutine Metropolis(walkers,testwalker,alpha,beta,a,R_L,R_R,Etest) 
  real(8), intent(in) :: testwalker(:), alpha, beta, a, R_L(:),R_R(:)
  real(8), intent(inout) :: walkers(:)
  real(8), intent(out) :: Etest 
  real(8) :: Psi1,Psi2,mu,B
  Psi1=Wavefunction(walkers,alpha,beta,R_L,R_R,a)
  Psi2=Wavefunction(testwalker,alpha,beta,R_L,R_R,a)  	!Wavefunction=Psi^2
  mu=psi2/psi1  
  call random_number(B)
  if (mu>B) then
    walkers=testwalker
    Etest=localE(walkers,alpha,beta,a,R_L,R_R)
  end if 	
end subroutine Metropolis

!Metropolis for the Initialization(without localE)
subroutine Metropolis_ini(walkers,testwalker,alpha,beta,a,R_L,R_R) 
  real(8), intent(in) :: testwalker(:),alpha,beta,a,R_L(:),R_R(:)
  real(8), intent(inout) :: walkers(:) 
  real(8) :: Psi1,Psi2,mu,B
  Psi1=Wavefunction(walkers,alpha,beta,R_L,R_R,a)
  Psi2=Wavefunction(testwalker,alpha,beta,R_L,R_R,a)  	!Wavefunction=Psi^2
  mu=psi2/psi1
  call random_number(B)  
  if (mu>B) then
    walkers=testwalker     				!accept new position only if mu>B
  end if 
end subroutine Metropolis_ini

!------------------------------------------------------------------------------------------------------------
!Local energy
real(8) function localE(walkers,alpha,beta,a,R_L,R_R)      
  real(8), intent(in) :: walkers(:), alpha, beta, a,R_L(:),R_R(:)
  real(8) :: r1(3), r2(3), r12(3), gama2
    r1=walkers(1:3)
    r2=walkers(4:6)
    r12=(r1-r2)/(sqrt(sum(r1-r2)**2))
    gama2=(1+beta*(sqrt(sum(r12**2))))**2
    localE=E1(r1,r12,alpha,gama2,a,R_L,R_R)+E2(r2,r12,alpha,gama2,a,R_L,R_R)+E12(r1,r2,alpha,gama2)-real(1)/(a**2) 	
  return
end function localE

real(8) function E1(r1,r12,alpha,gama2,a,R_L,R_R)
  real(8), intent(in) :: r1(:),r12(:),alpha,gama2,a,R_L(:),R_R(:)
  real(8) :: E1L,E1R,wf1l,wf1r,r1l(3),r1r(3)
    r1l=r1-R_L
    r1r=r1-R_R
    wf1l=exp(-sqrt(sum(r1l**2))/a)
    wf1r=exp(-sqrt(sum(r1r**2))/a)
    E1L=((1d0+(sum(r1l*r12))/(alpha*gama2))*wf1l/(a*(wf1l+wf1r))-1)/(sqrt(sum(r1l**2)))
    E1R=((1d0+(sum(r1r*r12))/(alpha*gama2))*wf1r/(a*(wf1l+wf1r))-1)/(sqrt(sum(r1r**2)))
    E1=E1L+E1R
  return
end function E1

real(8) function E2(r2,r12,alpha,gama2,a,R_L,R_R)
  real(8), intent(in) :: r2(:),r12(:),alpha,gama2,a,R_L(:),R_R(:)
  real(8) :: E2L,E2R,wf2l,wf2r,r2l(3),r2r(3)
    r2l=r2-R_L
    r2r=r2-R_R
    wf2l=exp(-sqrt(sum(r2l**2))/a)
    wf2r=exp(-sqrt(sum(r2r**2))/a)
    E2L=((1d0-(sum(r2l*r12))/(alpha*gama2))*wf2l/(a*(wf2l+wf2r))-1)/(sqrt(sum(r2l**2)))
    E2R=((1d0-(sum(r2r*r12))/(alpha*gama2))*wf2r/(a*(wf2l+wf2r))-1)/(sqrt(sum(r2r**2)))
    E2=E2L+E2R
  return
end function E2

real(8) function E12(r1,r2,alpha,gama2)
  real(8), intent(in) :: r1(:),r2(:),alpha,gama2
  real(8) :: r12(3),d12
    r12=r1-r2
    d12=sqrt(sum(r12**2))
    E12=(1-(2*alpha*sqrt(gama2)+d12)/((alpha*gama2)**2))/d12 
    return
end function E12

!-------------------------------------------------------------------------------------------------------------
!The test energy for a single alpha value, index k 
subroutine testE(Elocal,MC,N,k,Etest)
  integer, intent(in) :: MC,N,k
  real(8), intent(in) :: Elocal(:)
  real(8), intent(inout) :: Etest(:,:)
  integer :: i 
  do i=1,N
    Etest(i,k)=Elocal(i)/MC
  end do 
  !do i=1,N
   !print *, Etest(i,:)
  !end do
end subroutine testE

!---------------------------------------------------------------------------------------------------------------
! Newton-Raphson method, for cusp condition
real(8) function root(a,s)
  real(8), intent(in) :: a ,s
  root=(1+exp(-s/a))**(-1) - a
  return
end function root

real(8) function rootprime(a,s)
  real(8), intent(in) :: a ,s 
  rootprime=-s*exp(-s/a)*(a*(1+exp(-s/a)))**(-2)-1
  return
end function rootprime

subroutine NewtonRaphson(a,s)
  real(8), intent(in) :: s
  real(8), intent(out) :: a
  real(8)  :: leaveloop	!Stopping criterium 
  real(8)  :: f , f_prime 
  leaveloop=1d-10	
  a=-77			   ! Initial geuss for a 
  f=root(a,s)       ! The function at interest 
  do while (abs(f)>leaveloop) 
    f_prime=rootprime(a,s)       !Derivative of function 
    a=a-(f/f_prime)				 !NR-method 
    f=root(a,s)
  end do
end subroutine NewtonRaphson
 
!---------------------------------------------------------------------------------------------------------------
!This subroutine selects a random seed (starting point) for the speudo random number generator 
subroutine init_random_seed()
  integer :: i
  integer :: n
  integer :: clock
  integer, dimension(:), allocatable :: seed
  call random_seed(size = n)
  allocate(seed(n))
  call system_clock(count=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
  deallocate(seed)
end subroutine init_random_seed

end program QMC 
