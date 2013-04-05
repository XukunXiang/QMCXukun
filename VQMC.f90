module VQMC

implicit none
private

integer :: num_walkers
real(8), allocatable :: walkers(:, :), psi(:), E_L(:)
real(8) :: a, beta
real(8) :: R_L(3), R_R(3)

public init
public random_walk
public energy

contains
!---------------------------------------------------------------------------------------------------------------------------
subroutine init(n, s, b)
    integer, intent(in) :: n
    real(8), intent(in) :: s, b

    integer :: i

    num_walkers = n
    a = newton_raphson(s)
    beta = b

    R_L = [-s / 2, 0.0_8, 0.0_8]
    R_R = [ s / 2, 0.0_8, 0.0_8]

    if (allocated(walkers)) deallocate(walkers)
    allocate(walkers(6, num_walkers))

    if (allocated(psi)) deallocate(psi)
    allocate(psi(num_walkers))

    if (allocated(E_L)) deallocate(E_L)
    allocate(E_L(num_walkers))

    call random_number(walkers)
    walkers = 2 * walkers - 1
    do i = 1, num_walkers
        psi(i) = calc_psi(walkers(:, i))
        E_L(i) = calc_local_energy(walkers(:, i)) 
    end do
end subroutine
!----------------------------------------------------------------------------------------------------------------------
subroutine random_walk(dr,j)
    real(8), intent(inout) :: dr
    integer, intent(in) :: j
    integer :: i, accept
    real(8) :: testwalker(6) , u , psi_new, psi(num_walkers)

    accept = 0
    do i = 1, num_walkers
        call random_number(testwalker)
        testwalker = walkers(:, i) + (2 * testwalker - 1) * dr
        psi_new = calc_psi(testwalker)
        call random_number(u)
        if (u < psi_new**2 / psi(i)**2) then
            walkers(:, i) = testwalker
            psi(i) = psi_new
            E_L(i) = calc_local_energy(testwalker)
            accept = accept + 1
        end if
    end do
    dr = dr * accept / (num_walkers * 0.5_8)

    if (modulo(j,30000)==0) then
      do i = 1, num_walkers
	write (12,*) walkers(1,i),walkers(2,i),walkers(3,i)
        write (12,*) walkers(4,i),walkers(5,i),walkers(6,i)
      end do
    end if
end subroutine
!-----------------------------------------------------------------------------------------------------------------------------
!calculate up until the the Mu-th order moment of the energy 
subroutine energy(mu)
    real(8), intent(out) :: mu(:)
    integer :: i

    do i = 1, size(mu)
        mu(i) = sum(E_L**i) / size(E_L)
    end do
end subroutine
!------------------------------------------------------------------------------------------------------------------------------------
function calc_psi(position)
    real(8), intent(in) :: position(:)
    real(8) :: calc_psi
    real(8) :: r1(3), r2(3) , r1l(3), r2l(3), r1r(3), r2r(3), dr(3)

    r1 = position(1:3)
    r2 = position(4:6)
    
    r1l = r1-R_L
    r1r = r1-R_R
    r2l = r2-R_L
    r2r = r2-R_R
    dr  = r1-r2
    
 
    calc_psi = (exp(-sqrt(sum(r1l**2))/a) + exp(-sqrt(sum(r1r**2))/a))  &      ! Orbit1
              * (exp(-sqrt(sum(r2l**2))/a) + exp(-sqrt(sum(r2r**2))/a)) &      ! Orbit2
              * exp(sqrt(sum(dr**2))/(2*(1+beta*sqrt(sum((dr**2))))))          !Jastrow
    return 
end function

function calc_local_energy(position)
    real(8), intent(in) :: position(:)
    real(8) :: calc_local_energy
    
    real(8) ::  r1(3), r2(3), r12(3), gamma, d12, r1l(3),r1r(3),r2l(3),r2r(3),wf1l,wf1r,wf2r,wf2l 
    
    
    r1 = position(1:3)
    r2 = position(4:6)
    d12 = sqrt(sum((r1-r2)**2))
    r12 = (r1-r2)/d12
	    
    r1l = r1-R_L
    r1r = r1-R_R
    r2l = r2-R_L
    r2r = r2-R_R
    
    wf1l = exp(-sqrt(sum(r1l**2))/a)
    wf1r = exp(-sqrt(sum(r1r**2))/a)
    wf2l = exp(-sqrt(sum(r2l**2))/a)
    wf2r = exp(-sqrt(sum(r2r**2))/a)
    
    gamma = (1+beta*d12)
    
    calc_local_energy = ((1d0+(sum(r1l*r12))/(2*(gamma**2)))*wf1l/(a*(wf1l+wf1r))-1)/(sqrt(sum(r1l**2)))&
                       +((1d0+(sum(r1r*r12))/(2*(gamma**2)))*wf1r/(a*(wf1l+wf1r))-1)/(sqrt(sum(r1r**2)))&
	               +((1d0-(sum(r2l*r12))/(2*(gamma**2)))*wf2l/(a*(wf2l+wf2r))-1)/(sqrt(sum(r2l**2)))&
                       +((1d0-(sum(r2r*r12))/(2*(gamma**2)))*wf2r/(a*(wf2l+wf2r))-1)/(sqrt(sum(r2r**2)))& 
	               + (1-(2*2*gamma+d12)/(4*(gamma)**4))/d12                                         &
                      -real(1)/(a**2) 
    return
end function
!----------------------------------------------------------------------------------------------------------------------------------
function newton_raphson(s) result(a)
    real(8), intent(in) :: s
    real(8) :: a ,f ,f_prime , leaveloop
  leaveloop=1d-10	
  a=-77			   ! Initial geuss for a 
  f=(1+exp(-s/a))**(-1) - a      ! The function at interest 
  do while (abs(f)>leaveloop) 
    f_prime=-s*exp(-s/a)*(a*(1+exp(-s/a)))**(-2)-1       !Derivative of function 
    a=a-(f/f_prime)				 !NR-method 
    f=(1+exp(-s/a))**(-1) - a
  end do
end function

end module
