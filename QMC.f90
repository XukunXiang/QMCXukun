program QMC  		
  implicit none 
  integer, parameter :: N=400 ! number of walkers
  integer, parameter :: MC=30000! number of MC steps 
  integer, parameter :: Itime=4000 ! Initialization steps 
  integer, parameter :: intbeta=5 ! number of different alpha values 
  real(8):: walkers(6,N),testwalker(6),alpha,a, beta, Etest(N,intbeta), Elocal(N)     ! r(3,N) is a 3x1 for every walker. Don't save history. Only save last MC step
  real(8) :: R_L(3), R_R(3)         ! Position of nuclei  
  integer :: i,k,j 
  real(8) :: Emean(intbeta)
  real(8) :: s 

  beta=0.3        !set first beta value 
  s=3d0            !Set distance between nuclei (test, incorporate later) +
  R_L=(/-0.5*s,0d0,0d0/)
  R_R=(/+0.5*s,0d0,0d0/) ! Place nuclei symmetrically around zero, and along x-axis (only distance matters)
  alpha=2            ! From cusp condition
  call NewtonRaphson(a,s)  ! Determine parameter with cusp condition (only need to do this once)

!-------------------------------------------------------------------------------------------------------
!-------------------- MAIN PROGRAM AND SIMULATION-------------------------------------------------------

  do k=1,intbeta				!Parameter loop: \beta from 0.3 to 0.3+0.1*intbeta	
    Elocal(:)=0
    Etest(:,k)=0
    call init_random_seed()
    call random_number(walkers)                 !First positions you choose "randomly" --> Get close to minimum
    walkers=walkers-0.5					

    !Initialization Loop
    do i=1,Itime	
      do j=1,N							 
	call Randwalk(walkers(:,j),testwalker)
	call Metropolis_ini(walkers(:,j),testwalker,alpha,beta,N)  
	walkers(:,j+1)=testwalkers
      end do     end do 

    !Simulation
    do i=1,MC	
      do j=1,N
	!make loop for walkers  	  	!Start Monte Carlo simulation loop
	call Randwalk(walkers(:,j),testwalker)		  ! Let every particle take a small random 3D step 
	call Metropolis(walkers(:,j),testwalker,alpha,beta,N)      ! updates the r to the "newly accepted r's" // Check acceptance
	walkers(:,j+1)=testwalker 		!after Metropolis, testwalker is the new position, give it to walker(j+1) for the next loop
      end do 	
    end do
		
    call testE(Elocal,MC,N,k,Etest)      ! Etest(N,intbeta) will be mean Elocal for every walker for every beta
    Emean(k)=sum(Etest,k)/N                 !Use average over N independent walkers 
    !print *, beta, Emean(k)		
    beta=beta+0.1	             ! Choose different parameter 
  end do 

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
subroutine Randwalk(testwalker)
  real(8), intent(out):: testwalker(:)        
  real(8) :: A1(6) 
  call random_number(A1) 
  testwalker=walkers+1.0*(A1-0.5)
end subroutine Randwalk

!-------------------------------------------------------------------------------------------------------------
!Final Return :: Square of wavefunction |Psi|^2, a N-dimensional vector (for every walker)
function orbit(r,R_L,R_R,a)
  real(8), intent(in) :: r(:)
  real(8), intent(in) :: R_L(:), R_R(:)
  real(8), intent(in) :: a
  real(8) :: X_L(3), X_R(3)
  real(8) :: orbit
    X_L=r-R_L       ! Matrix (3,N) with "distance" (difference only) to Left nucleus for every walker
    X_R=r-R_R    
  orbit=exp(-sqrt(sum(X_L**2))/a) + exp(-sqrt(sum(X_R**2))/a)
  return
end function orbit

function Jastrow(r1,r2,alpha,beta)
  real(8), intent(in) :: r1(:), r2(:)
  real(8), intent(in) :: alpha, beta 
  real(8) :: Jastrow
  real(8) :: dr(3)
  dr=r1-r2
  Jastrow=exp(sqrt(sum(dr**2))/(alpha*(1+beta*sqrt(sum((dr**2)))))
  return
end function Jastrow

function Wavefunction(walkers,alpha,beta,a,N)
  integer, intent(in) :: N
  real(8), intent(in) :: alpha, beta, a 
  real(8), intent(in) :: walkers(:)
  real(8) :: r1(3), r2(3)
  real(8) :: Wavefunction
  r1=walkers(1:3)
  r2=walkers(4:6)
  Wavefunction=orbit(r1,R_L,R_R,a)*orbit(r2,R_L,R_R,a)*Jastrow(r1,r2,alpha,beta)
  Wavefunction=Wavefunction**2
  return
end function Wavefunction

!---------------------------------------------------------------------------------------------------------------
!Metropolis with |psi|^2 probability   (NOTE Walkers inside the subroutine is a vector! (We input walkers==walkers(:,j)
subroutine Metropolis(walkers,testwalker,alpha,beta,N) 
  real(8), intent(in) :: alpha,beta,walkers(:) 
  integer, intent(in) :: N
  real(8), intent(inout) :: testwalker(:)
  real(8) :: Psi1,Psi2,mu,B
  Psi1=Wavefunction(walkers,alpha,beta,a,N)
  Psi2=Wavefunction(testwalker,alpha,beta,a,N)  	!Wavefunction=Psi^2
  !Psi1=exp(-2*alpha*(sum(r**2,1)))      !Compare |psi|^2 according to defined Trial function 
  !Psi2=exp(-2*alpha*(sum(rtest**2,1)))
  mu=psi2/psi1  
!  if (mu>1) then
!    walkers=testwalker 
    !CALL LOCALE HERE TO SAVE COMPUTATIONAL TIME 
!  else			! Check acceptance of the new step 
    call random_number(B)
!    if (B<mu) then
!      walkers=testwalker
      !CALL LocalE HERE 
!      else
!	walkers=walkers  
	!don't call LocalE here, that is what will save us time ;)
!   end if
!  end if
  if (mu<B) then
    testwalker=walkers     !Testwalker will be the new position to give to walkers(j+1), so only if mu<B then testwalker=walkers
  end if 		
end subroutine Metropolis

!Metropolis for the Initialization(without localE)
subroutine Metropolis_ini(walkers,testwalker,alpha,beta,N) 
  real(8), intent(in) :: alpha,beta,testwalker(:)
  integer, intent(in) :: N
  real(8), intent(inout) :: walkers(:) 
  real(8) :: Psi1,Psi2,mu,B
  Psi1=Wavefunction(walkers,alpha,beta,a,N)
  Psi2=Wavefunction(testwalker,alpha,beta,a,N)  	!Wavefunction=Psi^2
  mu=psi2/psi1
  call random_number(B)  
  if (mu<B) then
    testwalker=walkers     !Testwalker will be the new position to give to walkers(j+1), so only if mu<B then testwalker=walkers
  end if 
end subroutine Metropolis_ini

!------------------------------------------------------------------------------------------------------------
!Local energy
subroutine localE(walkers,alpha,beta,a,Elocal,k)      
  real(8), intent(in) :: walkers(:), alpha, beta, a
  integer, intent(in) :: k
  real(8), intent(out) :: Elocal(:)
    r1=walkers(1:3)
    r2=walkers(4:6)
    r12=(r1-r2)/(sqrt(sum(r1-r2)**2))
    gama2=(1+beta*(sqrt(sum(r12**2))))**2
    Elocal(k)=E1(r1,r12,alpha,gama2,a,R_L,R_R)+E2(r2,r12,alpha,gama2,a,R_L,R_R)+E12(r12,alpha,gama2)-real(1)/(a**2) 	!can be improved, maybe we can use one function to get both E1 and E2
end subroutine localE
!Elocal=3*alpha+(0.5-2*(alpha**2))*sum(r**2,1)+Elocal ! NOW Calculate Cum sum for something like  Elocal(j) 
!bla=alpha+(0.5-2*(alpha**2))*sum(r**2,1)
!do i=1,N
  !print *, bla(i),Elocal(i)
!end do

function E1(r1,r12,alpha,gama2,a,R_L,R_R)
  real(8), intent(in) :: r1(:),r12(:),alpha,gama2,a,R_L(:),R_R(:)
  real(8), intent(out):: E1
  real(8) :: E1L,E1R,wf1l,wf1r,r1l(3),r1r(3)
    r1l=r1-R_L
    r1r=r1-R_R
    wf1l=exp(-sqrt(sum(r1l**2))/a)
    wf1r=exp(-sqrt(sum(r1r**2))/a)
    E1L=((1+(r1l*r12)/(alpha*gama2))*wf1l/(a*(wf1l+wf1r))-1)/(sqrt(sum(r1l**2)))
    E1R=((1+(r1r*r12)/(alpha*gama2))*wf1r/(a*(wf1l+wf1r))-1)/(sqrt(sum(r1r**2)))
    E1=E1L+E1R
  return
end function E1

function E2(r2,r12,alpha,gama2,a,R_L,R_R)
  real(8), intent(in) :: r2(:),r12(:),alpha,gama2,a,R_L(:),R_R(:)
  real(8), intent(out):: E2
  real(8) :: E2L,E2R,wf2l,wf2r,r2l(3),r2r(3),gama2
    r2l=r2-R_L
    r2r=r2-R_R
    wf2l=exp(-sqrt(sum(r2l**2))/a)
    wf2r=exp(-sqrt(sum(r2r**2))/a)
    E2L=((1-(r2l*r12)/(alpha*gama2))*wf1l/(a*(wf1l+wf1r))-1)/(sqrt(sum(r2l**2)))
    E2R=((1-(r2r*r12)/(alpha*gama2))*wf1r/(a*(wf1l+wf1r))-1)/(sqrt(sum(r2r**2)))
    E2=E2L+E2R
  return
end function E2

function E12(r12,alpha,gama2)
  real(8), intent(in) :: r12(:),alpha,gama2
  real(8), intent(out):: E12
  real(8) :: d12
    d12=sqrt(sum(r12**2))
    E12=(1-(2*alpha*sqrt(gama2)+d12))/((alpha*gama2)**2))/d12 
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
function root(a,s)
  real(8), intent(in) :: a ,s
  real(8) :: root 
  root=(1+exp(-s/a))**(-1) - a
  return
end function root

function rootprime(a,s)
  real(8), intent(in) :: a ,s 
  real(8) :: rootprime 
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
end subroutine

end program QMC 
