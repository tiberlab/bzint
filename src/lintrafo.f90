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
!           Module "linearTrafo"
!
!==============================================================
!
! Author: Matthias Auf der Maur
!
! TODO: describe functionality
!
!
!==============================================================

module linearTrafo


  implicit none

  ! ltrafo defines the components of a linear
  ! transformation x = Ax' + b
  ! NOTE: A and b are initially undefined, so they should
  !       be explicitly defined by a call to init_ltrafo
  type :: ltrafo
     double precision, dimension(2,2) :: A
     double precision, dimension(2)   :: b
  end type ltrafo


  !################################################################
  !
  ! Member functions/subroutines
  !
  !################################################################
contains


  !==============================================================
  !
  ! function "det":
  !
  ! Returns the determinant of a given linear transformation
  !
  !==============================================================
  !
  ! INPUT
  !    tr : type(ltrafo) a linear transformation
  !
  ! OUTPUT
  !    the determinant
  !
  !==============================================================
  function det(t)
    type(ltrafo), intent(in) :: t
    double precision         :: det

    det = t%A(1,1)*t%A(2,2) - t%A(1,2)*t%A(2,1)

  end function det
    

  !==============================================================
  !
  ! function "invert":
  !
  ! Invert a given linear transformation
  !
  !==============================================================
  !
  ! INPUT
  !    tr : type(ltrafo) trafo to invert
  !
  ! OUTPUT
  !    type(ltrafo) inverted trafo
  !
  !==============================================================
  function invert(t)
    type(ltrafo), intent(in) :: t
    type(ltrafo)             :: invert

    double precision d

    d = det(t)
    invert%A(1,1) = t%A(2,2)/d
    invert%A(1,2) = -t%A(1,2)/d
    invert%A(2,1) = -t%A(2,1)/d
    invert%A(2,2) = t%A(1,1)/d

    invert%b(1) = -(invert%A(1,1)*t%b(1) + invert%A(2,1)*t%b(2))
    invert%b(2) = -(invert%A(1,2)*t%b(1) + invert%A(2,2)*t%b(2))

  end function invert
 

  !==============================================================
  !
  ! subroutine "init_trafo":
  !
  ! Initialize an instance of type(ltrafo) with unity trafo
  !
  !==============================================================
  !
  ! INPUT/OUTPUT
  !    tr : type(ltrafo) trafo to initialize
  !
  !==============================================================
  subroutine init_ltrafo(tr)
    type(ltrafo), intent(inout) :: tr

    tr%A = reshape( (/ 1, 0, &
                       0, 1 /), (/ 2, 2 /) )
    tr%b = (/ 0, 0 /)

  end subroutine init_ltrafo


end module linearTrafo

