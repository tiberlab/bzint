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
!           Module "qSurface"
!
!==============================================================
!
! Author: Matthias Auf der Maur
!
! TODO: describe functionality
!
!
!==============================================================

module qSurface

  ! implements routines for finding real roots of equations
  use realroots
  ! implements geometric types
  use geometry
  ! defines linear transformation
  use linearTrafo

  implicit none



  ! qcurve defines a quadratic curve in reduced form.
  ! The reduced form is obtained by applying an affine
  ! transformation x = Ax' + b, |A| = 1
  !
  !   coeff       : coefficients q1..q6 of the reduced
  !                 quadratic form
  !   trsf        : the affine transformation
  !                 that led to the reduced form
  !                 transformed simplex
  !   curvetype   : the type of the quadratic curve as defined
  !                 above
  type :: qcurve
     double precision, dimension(6) :: coeff
     type(ltrafo)                   :: trsf
     integer                        :: curvetype
  end type qcurve


  ! Some constants to identify the different types of quadratic
  ! curves.
  integer, parameter :: qELLIPSE   = 1, &
                        qHYPERBOLA = 2, &
                        qPARABOLA  = 3, &
                        qLINE      = 4, &
                        qDEG_LINES = 5, &
                        qDEG_DELTA = 6, &
                        qUNKNOWN   =-1


  !################################################################
  !
  ! Interfaces
  !
  !################################################################


  !################################################################
  !
  ! Member functions/subroutines
  !
  !################################################################
contains

  !==============================================================
  !
  ! Function "determine_curve_type":
  !
  ! Determines the type of the quadratic curve parametrized by
  ! the 6 coefficients exp_coeff:
  !    const = q1 + q2*x + q3*y + q4*x^2 + q5*x*y + q6*y^2
  ! The type can be determined by finding a affine transformation
  ! that transforms the quadratic form to an appropriate standard
  ! representation.
  !
  ! The type can be one of:
  !   qELLIPSE    -> q2, q3, q5 = 0, sgn(q4) = sgn(q6)
  !   qHYPERBOLA  -> q2, q3, q4, q6 = 0
  !   qPARABOLA   -> q2, q4, q6 = 0
  !   qLINE       -> q3, q4, q5, q6 = 0
  !   qDEG_LINES  -> q2, q3, q5, q6 = 0 (degenerate form, two
  !                                      parallel lines)
  !   qDEG_DELTA  -> q2, q3, q4, q5, q6 = 0 (degenerate form)
  !
  !==============================================================
  !
  ! INPUT
  !    exp_coeff : array containing the 6 expansion coefficients
  !
  ! OUTPUT
  !    type of curve including new coeff. and transformation
  !    (cf. definition of type qcurve)
  !
  !==============================================================
  function determine_curve_type(exp_coeff)
    double precision, dimension(6), intent(in) :: exp_coeff
    type(qcurve)                               :: determine_curve_type

    ! will store the transformation that transforms the general
    ! quadratic form of the intput to the reduced form.
    type(ltrafo) :: trafo

    ! for convenience
    double precision :: q1, q2, q3, q4, q5, q6

    ! the biggest of the coefficients
    double precision :: maxcoeff
    
    ! local temp variables for calculations
    double precision :: tmp1, tmp2, tmp3
    double precision, dimension(2)   :: eigs, b
    double precision, dimension(2,2) :: matrix

    ! a small number compared to 1
    double precision, parameter :: EPS = 1D-9

    ! We scale the coefficients by maxcoeff so they will
    ! always be between -1 <= qi <= 1
    maxcoeff =  max(abs(exp_coeff(1)),abs(exp_coeff(2)), &
                    abs(exp_coeff(3)),abs(exp_coeff(4)), &
                    abs(exp_coeff(5)),abs(exp_coeff(6)))


    ! the scaled coefficients
    q1 = exp_coeff(1)/maxcoeff
    q2 = exp_coeff(2)/maxcoeff
    q3 = exp_coeff(3)/maxcoeff
    q4 = exp_coeff(4)/maxcoeff
    q5 = exp_coeff(5)/maxcoeff
    q6 = exp_coeff(6)/maxcoeff

    ! initialize trafo to identity
    call init_ltrafo(trafo)

    !---------------------------------------------------
    ! Search affine transformation of quadratic form to
    ! one of the reduced forms
    !---------------------------------------------------
    if (abs(q5) > EPS) then
       ! there exists a trafo x = Ax', |A| =  1 such that
       ! q5 = 0. Write the expansion as
       !    e = q1 + k'*t + k'*S*k,
       ! where
       !    k' = transposed of k
       !    t  = [q2; q3]
       !    S  = 0.5*[2*q4 q5; q5 2*q6]
       ! The transformed equation is then
       !    e = q1 + q'*(A'*t) + q'*(A'*S*A)*q
       ! The searched A is therefore the Matrix, that
       ! diagonalizes S (and fullfills A' = A^-1)

       ! calculate eigenvalues of Matrix S
       call qroots(1.0D0, -q4-q6, q4*q6 - 0.25D0*q5*q5, &
                   eigs)
       ! NOTE: there are always two different eigenvalues,
       ! as q5 /= 0 and real

       ! trsf matrix A consists of orthonormal eigenvectors
       ! of S
       tmp1 = 2.0D0*(eigs(1) - q4)/q5
       tmp3 = sqrt(1.0D0 + tmp1*tmp1)
       trafo%A(1,1) = 1.0D0/tmp3
       trafo%A(2,1) = tmp1/tmp3

       tmp2 = 2.0D0*(eigs(2) - q4)/q5
       tmp3 = sqrt(1.0D0 + tmp2*tmp2)
       trafo%A(1,2) = 1.0D0/tmp3
       trafo%A(2,2) = tmp2/tmp3

       ! calculate the new qi's
       tmp1 = trafo%A(1,1)*q2 + trafo%A(2,1)*q3
       tmp2 = trafo%A(1,2)*q2 + trafo%A(2,2)*q3
       ! q1 remains
       q2 = tmp1
       q3 = tmp2
       q4 = eigs(1)
       q5 = 0.0D0
       q6 = eigs(2)
    end if

    if ((abs(q4) < EPS) .and. (abs(q6) > EPS)) then 
       ! x = Ax' with  x = y', y = x'
       ! update trafo (x = A1x', x' = A2x'' => x = A1*A2x'')
       matrix = reshape( (/ 0, 1, &
                            1, 0 /), (/ 2, 2 /) ) ! A2

       trafo%A = matmul(trafo%A, matrix)

       ! new qi's
       tmp1 = q2
       ! q1 remains
       q2 = q3
       q3 = tmp1
       q4 = q6
       ! q5 is zero
       q6 = 0.0D0
    end if

    if ((abs(q4) > EPS) .and. (abs(q2) > EPS)) then
       ! x = x' - q2/(2*q4), y = y'
       ! update trafo (x = Ax', x' = x'' + b => x = Ax'' + Ab)
       b(1) = -q2/(2.0D0*q4)
       b(2) = 0
       trafo%b = matmul(trafo%A, b)

       ! new qi's
       q1 = q1 - 0.25D0*q2*q2/q4
       q2 = 0.0D0
       ! q3, q4, q6 remain
    end if

    if ((abs(q3) > EPS) .and. (abs(q6) > EPS)) then
       ! x = x', y = y' - q3/(2*q6)
       ! update trafo (x = Ax' + b, x' = x'' + c => x = Ax'' + Ac + b)
       b(1) = 0.0D0
       b(2) = -q3/(2.0D0*q6)
       trafo%b = trafo%b + matmul(trafo%A, b)

       ! new qi's
       q1 = q1 - 0.25D0*q3*q3/q6
       q3 = 0.0D0
       ! q2, q4, q6 remain
    end if
    !---------------------------------------------------------
    ! now we have one of these forms:
    !  e = q1 + q4*x^2 + q6*y^2 , q4 /= 0 & q6 /= 0
    !  e = q1 + q3*y +q4*x^2    , q4 /= 0
    !  e = q1 + q2*x + q3*y
    !---------------------------------------------------------

    !---------------------------------------------------------
    ! Now examine the reduced quadratic form to get curve type
    ! of the iso-energy curve
    !---------------------------------------------------------
    if (abs(q6) > EPS) then
       if (sign(1.0D0, q4/q6)  > 0.0D0) then ! sgn(q4) = sgn(q6)
          !! => We have an ellipse
          !!
          determine_curve_type%curvetype = qELLIPSE
       else
          !! => We have a hyperbola
          !!
          determine_curve_type%curvetype = qHYPERBOLA
          ! in this case, we reduce further by means of
          ! x = Ax' such that q4 = q6 = 0, q5 /= 0
          ! x = x'/sqrt(2) + sqrt(|q6|/|2*q4|)*y'
          ! y =  sqrt(|q4|/|2*q6|)*x' - y'/sqrt(2)
          matrix(1,1) = sqrt(0.5D0)
          matrix(1,2) = sqrt(0.5D0*abs(q6)/abs(q4))
          matrix(2,1) = sqrt(0.5D0*abs(q4)/abs(q6))
          matrix(2,2) = -sqrt(0.5D0)
          ! update trafo structure
          trafo%A = matmul(trafo%A, matrix)

          ! new qi's
          ! q1 remains
          ! q2, q3 are already 0
          q5 = q4*sqrt(abs(q6)/abs(q4)) - q6*sqrt(abs(q4)/abs(q6))
          q4 = 0.0D0
          q6 = 0.0D0
          ! now we have e = q1 + q5*x*y
       end if
    else if (abs(q4) > EPS) then
       if (abs(q3) > EPS) then
          !! => We have a parabola
          !!
          determine_curve_type%curvetype = qPARABOLA
       else
          !! => We have two lines (degenerate case)
          !!
          determine_curve_type%curvetype = qDEG_LINES
       end if
    else if ((abs(q2) > EPS) .or. (abs(q3) > EPS)) then
       !! => We have a straight line
       !!
       ! NOTE: if q2 .or. q3 are very small, we treat this case
       !       as degenerated
       determine_curve_type%curvetype = qLINE
       ! in this case, we transform to always have q3 = 0
       ! x = x'*cos(t) + y'*sin(t), y = x'*sin(t) - y'*cos(t)
       tmp1 = sqrt(q2*q2 + q3*q3)
       tmp2 = q2/tmp1 ! cos(t)
       tmp3 = q3/tmp1 ! sin(t)
       ! trafo matrix
       matrix = reshape( (/ tmp2, tmp3, &
                            tmp3, -tmp2 /), (/ 2, 2 /) )
       ! update trafo structure
       trafo%A = matmul(trafo%A, matrix)

       ! new qi's
       ! q1 remains
       q2 = q2*tmp2 + q3*tmp3
       q3 = 0.0D0
       ! q4, q5, q6 are already 0
    else  
       !! => We have a degenerate case
       !!
       determine_curve_type%curvetype = qDEG_DELTA
    end if

    determine_curve_type%coeff = (/ q1, q2, q3, q4, q5, q6 /)*maxcoeff
    determine_curve_type%trsf = trafo

  end function determine_curve_type



  !==============================================================
  !
  ! subroutine "minmax":
  !
  ! Calculate the global minimum and maximum of a quadratic
  ! surface in the unit triangle
  !
  !==============================================================
  !
  ! INPUT
  !    coeff :  coefficients q1, ..., q6 of quadratic function
  !
  ! OUTPUT
  !    the minimum and maximum
  !
  !==============================================================
  subroutine minmax(coeff, minval, maxval)
    double precision, dimension(6), intent(in) :: coeff
    double precision, intent(out)              :: minval, maxval

    ! a determinant
    double precision :: D

    ! temporary variables
    double precision :: tmp, x, y


    ! look at the energy values at the corners
    minval = coeff(1) + min(0.0D0, coeff(2)+coeff(4), &
                            coeff(3)+coeff(6))
    maxval = coeff(1) + max(0.0D0, coeff(2)+coeff(4), &
                            coeff(3)+coeff(6))

    ! look for a local minimum/maximum
    D = 4.0D0*coeff(4)*coeff(6) - coeff(5)**2
    if (D > 0.0D0) then
      ! there is a local minimum or maximum
      x = (coeff(3)*coeff(5) - 2.0D0*coeff(2)*coeff(6))/D
      y = (coeff(2)*coeff(5) - 2.0D0*coeff(3)*coeff(4))/D

      tmp = coeff(1) + coeff(2)*x + coeff(3)*y + &
            coeff(4)*x*x + coeff(5)*x*y + coeff(6)*y*y
      if (coeff(4) > 0.0D0) then
        ! we have a minimum
        minval = min(minval, tmp)
      else
        ! we have a maximum
        maxval = max(maxval, tmp)
      end if
    end if

    ! look for minima/maxima with boundary condition
    ! look on 1. side: y = 0
    if (abs(coeff(4)) > 0.0D0) then
      x = - 0.5D0*coeff(2)/coeff(4)
      tmp = coeff(1) + coeff(2)*x + coeff(4)*x*x
      minval = min(minval, tmp)
      maxval = max(maxval, tmp)
    end if

    ! look on 2. side: x + y = 1
    tmp = coeff(4) - coeff(5) + coeff(6)
    if (abs(tmp) > 0.0D0) then
      x = 0.5D0*(coeff(3)-coeff(2)+2.0D0*coeff(6)-coeff(5))/tmp
      tmp = coeff(1) + coeff(3) + coeff(6) + &
            (coeff(2) - coeff(3) + coeff(5) - 2.0D0*coeff(6))*x + &
            tmp*x*x
      minval = min(minval, tmp)
      maxval = max(maxval, tmp)
    end if

    ! look on 3. side: -x = 0
    if (abs(coeff(6)) > 0.0D0) then
      y = - 0.5D0*coeff(3)/coeff(6)
      tmp = coeff(1) + coeff(3)*y + coeff(6)*y*y
      minval = min(minval, tmp)
      maxval = max(maxval, tmp)
    end if

  end subroutine minmax


  function is_on_the_left(p, coeff, curvetype)
    type(point), intent(in)                    :: p
    double precision, dimension(6), intent(in) :: coeff
    integer, intent(in)                        :: curvetype
    integer                                    :: is_on_the_left

    double precision :: tmp1, tmp2

    is_on_the_left = -1
    
    select case (curvetype)
    case (qELLIPSE)
       tmp1 = (-coeff(1) - coeff(4)*(p%x**2))/coeff(6)
       tmp2 = p%y**2
       if (tmp2 < tmp1) then
          is_on_the_left = 1
       elseif (tmp2 == tmp1) then
          is_on_the_left = 0
       end if

    case (qHYPERBOLA)
       ! do something
    end select
             
  end function is_on_the_left


end module qSurface

