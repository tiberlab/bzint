!
! This file is part of bzint.
!
! bzint is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as
! published by the Free Software Foundation, either version 3 of
! the License, or (at your option) any later version.
!
! bzint is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General
! Public License along with bzint. If not, see
! <https://www.gnu.org/licenses/>.
!
!==============================================================
!
!           Module "realroots"
!
!==============================================================
!
! Author: Matthias Auf der Maur
!
! functions to calculate the real roots of different equations
!
!==============================================================


module realroots

  implicit none

contains

  !==============================================================
  !
  ! Subroutine "qroots":
  !
  ! Calculate the two real roots of a quadratic equation
  !    a*x^2 + b*x + c = 0
  !
  ! NOTE: The coefficients a or b are allowed to be zero.
  !       Non-existing (or complex) roots will be set to HUGE.
  !
  !==============================================================
  !
  ! INPUT
  !    coefficients a, b, c of quadratic equation
  !
  ! OUTPUT
  !     roots : the two roots
  !     n     : the number of roots (= 0,1,2)
  !
  !==============================================================
  subroutine qroots(a, b, c, roots, n)
    double precision, intent(in)                :: a, b, c
    double precision, dimension(2), intent(out) :: roots
    integer, intent(out), optional              :: n

    ! local variables for calculation
    double precision :: f, r, s, t, d
    ! temporary n
    integer :: tmpn

    ! a small number
    double precision, parameter :: EPS = 1.0D-12

    ! calculate scaled coefficients
    f = max(abs(a), abs(b), abs(c))
    r = a/f
    s = b/f
    t = c/f
    
    if (abs(r) > EPS) then
       !! We have a true quadratic equation
       !! ---------------------------------
       ! the determinant
       d = s*s - 4.0D0*r*t
       if (d >= 0.0D0) then
          d = sqrt(d)
          ! solve for first root of modified equation
          if (s > 0.0D0) then
             roots(1) = (-s-d)/(r+r)
          else
             roots(1) = (-s+d)/(r+r)
          end if
          ! second root
          roots(2) = (t/r)/roots(1)
          tmpn = 2
          if (abs(roots(2) - roots(1)) < EPS*abs(roots(1))) then
             ! a rather hypothetical case ...
             roots(2) = huge(1.0D0)
             tmpn = 1
          end if
       else
          ! the roots would be imaginary, which is considered in
          ! this context as non existent
          roots(1) = huge(1.0D0)
          roots(2) = huge(1.0D0)
          tmpn = 0
       end if
    else
       !! Hmm, seems to be linear
       !! -----------------------
       if (abs(s) >= EPS) then
          !! In fact, its linear
          roots(1) = -c/b
          roots(2) = huge(1.0D0)
          tmpn = 1
       else
          !! Well, no roots at all (we assume, that at least c is
          !! non-zero)
          !! ----------------------------------------------------
          roots(1) = huge(1.0D0)
          roots(2) = huge(1.0D0)
          tmpn = 0
       end if
    end if

    if (present(n)) then
       n = tmpn
    end if

  end subroutine qroots



  !==============================================================
  !
  ! Subroutine "trroots":
  !
  ! Calculate the two real roots of a trigonometric equation of
  ! the form a*cos(x) + b*sin(x) + c = 0 inside the interval
  ! -pi <= x < pi
  !
  ! NOTE: The coefficients a or b are allowed to be zero.
  !       Non-existing (or complex) roots will be set to HUGE.
  !
  !==============================================================
  !
  ! INPUT
  !    coefficients a, b, c of the trigonometric equation
  !
  ! OUTPUT
  !     roots : the two roots
  !     n     : (optional) the number of roots (= 0,1,2)
  !
  !==============================================================
  subroutine trroots(a, b, c, roots, n)
    double precision, intent(in)                :: a, b, c
    double precision, dimension(2), intent(out) :: roots
    integer, intent(out), optional              :: n

    ! the roots of the transformed quadratic eq
    double precision, dimension(2) :: tmproots
    ! the number of roots
    integer :: tmpn

    ! a small number
    double precision, parameter :: EPS = 1.0D-12


    ! Make a transformation xi = tan(x/2)
    ! With this sin(x) = 2*xi/(1+xi^2), cos(x) = (1-xi^2)/(1+xi^2)
    ! and (c-a)*xi^2 + 2*b*xi + a + c = 0
    call qroots(c - a, 2.0D0*b, a + c, tmproots, tmpn)

    select case (tmpn)
    case (0)
       ! there is no root
       roots = huge(1.0D0)
       tmpn = 0
    case (1)
       ! there's one root (or two roots, where the 2nd lies in 
       ! x = -pi. This is checked seperately)
       roots(1) = 2.0D0*atan(tmproots(1))
       ! check x = -pi: -a + c = 0
       if (abs(a - c) < EPS) then
          ! -pi is a root
          roots(2) = -acos(-1.0D0)
          tmpn = 2
       else
          roots(2) = huge(1.0D0)
          tmpn = 1
       end if
    case (2)
       ! there are two roots
       roots(1) = 2.0D0*atan(tmproots(1))
       roots(2) = 2.0D0*atan(tmproots(2))
       tmpn = 2
    end select

    if (present(n)) then
       n = tmpn
    end if

  end subroutine trroots
 
end module realroots
