!======================================================================
! TRIDAG – Solves tridiagonal linear systems using Thomas algorithm
! A * x = d
! Inputs:
!   IL : first equation index
!   IU : last  equation index
!   AA : subdiagonal (below main diag)
!   BB : diagonal
!   CC : superdiagonal
!   DD : RHS (solution overwrites DD)
!======================================================================

subroutine TRIDAG(IL, IU, AA, BB, CC, DD)
implicit real*8 (A-H, O-Z)

dimension AA(IU), BB(IU), CC(IU), DD(IU)

integer :: LP, I, J
real(8) :: R

! Forward elimination
LP = IL + 1
do I = LP, IU
    R     = AA(I) / BB(I-1)
    BB(I) = BB(I) - R * CC(I-1)
    DD(I) = DD(I) - R * DD(I-1)
end do

! Back-substitution
DD(IU) = DD(IU) / BB(IU)

do I = LP, IU
    J = IU - I + IL
    DD(J) = (DD(J) - CC(J)*DD(J+1)) / BB(J)
end do

return
end
