!======================================================================
! 1D Linear / Nonlinear Wave Propagation
! Methods: Lax, Upwind, Lax-Wendroff, Beam-Warming (Implicit)
! Author: Mohammad E. Heravifard
!======================================================================

program Project3
implicit none

integer :: i, N, Iter, Iter_max, k
real(8) :: c, CFL, Pi, t, tf, error, x_max, dx, dt
real(8) :: Aim1, Aip1, Eim1, Eip1
real(8), allocatable :: u_old(:), u_new(:), E(:), x(:)
real(8), allocatable :: low(:), ondiag(:), up(:), rhs(:), du(:)

!------------------------------------------------------------
! Parameters
!------------------------------------------------------------
N     = 51
c     = 1.0d0
x_max = 5.0d0
CFL   = 2.0d0
tf    = 2.0d0

Pi = acos(-1.0d0)

allocate(u_old(N), u_new(N), x(N), E(N), low(N), ondiag(N), up(N), rhs(N), du(N))

dx = x_max / (N - 1)
dt = dx*CFL / c
Iter_max = nint(tf / dt - 1e-6)

! Spatial grid
do i = 1, N
    x(i) = (i - 1) * dx
end do

!------------------------------------------------------------
! Initial & boundary conditions
!------------------------------------------------------------
u_old(1) = 0.0d0
u_new(1) = 0.0d0

! Discontinuous IC example (use Smooth IC by editing these)
do i = 2, 11
    u_old(i) = 1.0d0
end do

do i = 11, 22
    u_old(i) = sin(Pi*(x(i)-1.0d0))
end do

do i = 22, N
    u_old(i) = 1.0d0
end do

Iter = 0
t    = 0.0d0

!------------------------------------------------------------
! Output initial data
!------------------------------------------------------------
open (2020, file="Solution.plt")
write(2020,*) 'variables="x", "u"'
write(2020,*) 'zone T="t=', t, '"'
do i = 1, N
    write(2020,*) x(i), u_old(i)
end do

!------------------------------------------------------------
! Time marching loop
!------------------------------------------------------------
do while (Iter < Iter_max)

    t = t + dt
    Iter = Iter + 1
    write(*,*) 'Iteration = ', Iter

    do i = 1, N

        !----------------------------------------------------
        ! Uncomment ONE scheme only (linear or nonlinear)
        !----------------------------------------------------

        !-----------------------------
        ! Linear Lax method
        !-----------------------------
        !u_new(i) = 0.5d0*(u_old(i+1) + u_old(i-1)) - 0.5d0*CFL*(u_old(i+1) - u_old(i-1))

        !-----------------------------
        ! Nonlinear Lax method
        !-----------------------------
        !u_new(i) = 0.5d0*(u_old(i+1) + u_old(i-1)) - &
        !           0.25d0*CFL*((u_old(i+1)**2) - (u_old(i-1)**2))

        !-----------------------------
        ! Linear Upwind
        !-----------------------------
        !u_new(i) = u_old(i) - CFL*(u_old(i) - u_old(i-1))

        !-----------------------------
        ! Nonlinear Upwind
        !-----------------------------
        !u_new(i) = u_old(i) - 0.5*CFL*(u_old(i)**2 - u_old(i-1)**2)

        !-----------------------------
        ! Linear Lax-Wendroff
        !-----------------------------
        !u_new(i) = u_old(i) - 0.5*CFL*(u_old(i+1) - u_old(i-1)) + &
        !           0.5*CFL*CFL*(u_old(i+1)-2*u_old(i)+u_old(i-1))

        !-----------------------------
        ! Nonlinear Lax-Wendroff
        !-----------------------------
        !E(i) = 0.5*u_old(i)*u_old(i)
        !u_new(i) = u_old(i) - 0.5*CFL*(E(i+1)-E(i-1)) + &
        ! 0.25*CFL*CFL*((u_old(i+1)+u_old(i))*(E(i+1)-E(i)) - &
        ! (u_old(i)+u_old(i-1))*(E(i)-E(i-1)))

        !====================================================
        ! Beam-Warming (Implicit)
        !====================================================
        if(i == 1 .or. i == N) then
            ondiag(i) = 1.0d0
            up(i)     = 0.0d0
            low(i)    = 0.0d0
            rhs(i)    = 0.0d0
        else
            !--------- Linear BW ---------
            Aim1 = c
            Aip1 = c
            Eim1 = c*u_old(i-1)
            Eip1 = c*u_old(i+1)

            !--------- Nonlinear BW (uncomment below) --------
            !Aim1 = u_old(i-1)
            !Aip1 = u_old(i+1)
            !Eim1 = 0.5d0*u_old(i-1)**2
            !Eip1 = 0.5d0*u_old(i+1)**2

            ondiag(i) = 1.0d0
            up(i)     =  dt*Aip1/(4*dx)
            low(i)    = -dt*Aim1/(4*dx)
            rhs(i)    =  u_old(i) - dt*(Eip1-Eim1)/(2*dx) + &
                         dt*Aip1*u_old(i+1)/(4*dx) - dt*Aim1*u_old(i-1)/(4*dx)
        end if

    end do

    ! Solve tridiagonal system
    call TRIDAG(1, N, low, ondiag, up, rhs)

    du = rhs
    du(1) = 0.0d0
    du(N) = 0.0d0

    u_new = du
    u_old = u_new

end do

!------------------------------------------------------------
! Final output
!------------------------------------------------------------
write(2020,*) 'zone T="t=', t, '"'
do i = 1, N
    write(2020,*) x(i), u_old(i)
end do

pause
end program Project3
