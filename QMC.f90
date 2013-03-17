program QMC  		
 implicit none 
 integer, parameter :: N=400 ! number of walkers
 integer, parameter :: MC=30000 ! number of MC steps 
 integer, parameter :: Itime=4000 ! Initialization steps 
 integer, parameter :: intalpha=5 ! number of different alpha values 
 real(8):: r(3,N), rtest(3,N), alpha ,Etest(N,intalpha), Elocal(N)     ! r(3,N) is a 3x1 for every walker. Don't save history. Only save last MC step 
 integer :: i,k,j 
 real(8) :: Emean(intalpha)

  alpha=0.3        !set first alhpa value 
!-------------------------------------------------------------------------------------------------------------------------
!-------------------- MAIN PROGRAM AND SIMULATION--------------------------------------------------------------------------
do k=1,intalpha								  !Parameter loop	
Elocal(:)=0
Etest(:,k)=0
  call init_random_seed()
  call random_number(r)                       ! First positions you choose "randomly" --> Get close to minimum
   r=r-0.5 
   
		do j=1,Itime                          !Initialisation loop
	        call Randwalk(r,rtest)             
			call Metropolis(rtest,r,alpha,N)   
		end do 
		
		do i=1,MC                             ! Start Monte Carlo simulation loop
			call Randwalk(r,rtest)            ! Let every particle take a small random 3D step 
			call Metropolis(rtest,r,alpha,N)  ! updates the r to the "newly accepted r's" // Check acceptance
			call localE(r,alpha,Elocal,N)     !Cumulative sum of Elocal 
		end do
		
  call testE(Elocal,MC,N,k,Etest)      ! Etest(N,intalpha) will be mean Elocal for every walker for every alpha
  Emean=sum(Etest,1)/N                 !Use average over N independent walkers 
  print *, alpha, Emean(k)		
        alpha=alpha+0.1	             ! Choose different parameter 
end do
!-------------------------------------------------------------------------------------------------------------------------
!-------------------POST PROCESSING, PRINT STATEMENTS --------------------------------------------------------------------
!do i=1,N
!print *, Etest(i,1), Elocal(i) , MC
!end do


 !  do i=1,2
  !  print *, alpha, Etest(:,i)
  ! end do
  
contains 
! --------------------------------------------------------------------------------------------------------
!--------------------------------SUBROUTINES AND FUNCTIONS-------------------------------------------------------------------- 
!-------------------------------------------------------------------------------------------------------------------------

! Take a small step forward or backwards, for every particle. select x,y and z positions 
subroutine Randwalk(r,rtest)
real(8), intent(in) :: r(:,:)
real(8), intent(out):: rtest(:,:)
real(8) :: A1,A2,A3
integer :: i 
	do i=1,N
		call random_number(A1) 
		call random_number(A2)
		call random_number(A3)
		rtest(1,i)=r(1,i)+1.0*(A1-0.5)
		rtest(2,i)=r(2,i)+1.0*(A2-0.5)
		rtest(3,i)=r(3,i)+1.0*(A3-0.5)
	end do 
end subroutine Randwalk

!-------------------------------------------------------------------------------------------------------------
! Return:: Square of wavefunction |Psi|^2
function Wavefunction(r,N)
integer, intent(in) :: N
real(8), intent(in) :: r(:,:)
real(8) :: Wavefunction(N)
Wavefunction=exp(-2*alpha*(sum(r**2,1)))
return
end function
!---------------------------------------------------------------------------------------------------------------
!Metropolis with |psi|^2 probability 
subroutine Metropolis(rtest,r,alpha,N)
real(8), intent(in) :: alpha,rtest(:,:)
integer, intent(in) :: N
real(8), intent(inout) :: r(:,:)
integer :: i 
real(8) :: Psi1(N),Psi2(N),mu,B
 Psi1=Wavefunction(r,N)
 Psi2=Wavefunction(rtest,N)
 !Psi1=exp(-2*alpha*(sum(r**2,1)))      !Compare |psi|^2 according to defined Trial function 
 !Psi2=exp(-2*alpha*(sum(rtest**2,1)))
do i=1,N                     ! For every walker
  mu=psi2(i)/psi1(i)
	if (mu>1) then 
	r(:,i)=rtest(:,i)
  	else                           ! Check acceptance of the new step 
		call random_number(B)
			if (B<mu) then
			r(:,i)=rtest(:,i)
			else
			r(:,i)=r(:,i) 
			end if 
	end if 
end do		
end subroutine Metropolis
!------------------------------------------------------------------------------------------------------------

!Local energy
subroutine localE(r,alpha, Elocal,N)
real(8),intent(in) :: alpha, r(:,:)
integer, intent(in) :: N 
real(8), intent(inout) :: Elocal(:)	
integer :: i
real(8) :: bla(N)
	  Elocal=3*alpha+(0.5-2*(alpha**2))*sum(r**2,1)+Elocal ! Local energy for all the walkers , cummulative sum 
	  !bla=alpha+(0.5-2*(alpha**2))*sum(r**2,1)
	  !do i=1,N
	   !print *, bla(i),Elocal(i)
	  !end do
end subroutine localE
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

subroutine pushtest()
 
end subroutine
  
end program QMC 