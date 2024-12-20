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
!           Module "simplex2D"
!
!==============================================================
!
! Author: Matthias Auf der Maur
!
! This module implements the integration in a 2D simplex.
!
!
! General notes:
! --------------
!
! Definition of simplex:
!
!                   3
!                   x
!                  / \
!                 /   \
!              6 x     x 5
!               /       \
!              /         \
!             x-----x-----x
!            1      4      2
!
!     The points 4,5,6 are the midpoints of the triangle sides
!
!     The numbering corresponds with the indices to access the
!     function/energy/coordinate values at the respective points
!
! The algorithm follows the description in:
!   G Wiesenekker, G te Velde and E J Baerends: "Analytic
!     quadratic integration over the two dimensional Brillouin
!     zone". J. Phys. C: Solid State Phys. 21 (1988) 4263-4283
!
!==============================================================

module simplex2D

  ! implements routines for finding real roots of equations
  use realroots
  ! defines point, triangle
  use geometry
  ! defines quadratic surfaces/curves
  use qSurface
  ! defines linear transformations
  use linearTrafo

  implicit none

  ! as default, everything is hidden
  private
  ! these are the members that should be accessible
  public :: simplex2Dqint, THETA, DELTA, sort




  !################################################################
  !
  ! Public variables
  !
  !################################################################
  integer, parameter :: DELTA = 0, &
                        THETA = 1

  !################################################################
  !
  ! Private types
  !
  !################################################################

  !################################################################
  !
  ! Private variables
  !
  !################################################################

  ! exp_matrix is the inverse of the matrix containing the
  ! coefficients of the quadratic expansion of function f(x,y)
  ! in the unit triangle:
  ! TODO: comment further
  real, parameter, dimension(6,6) :: exp_matrix = &
       reshape( (/ 1.0,-3.0,-3.0, 2.0, 4.0, 2.0, &
                   0.0,-1.0, 0.0, 2.0, 0.0, 0.0, &
                   0.0, 0.0,-1.0, 0.0, 0.0, 2.0, &
                   0.0, 4.0, 0.0,-4.0,-4.0, 0.0, &
                   0.0, 0.0, 0.0, 0.0, 4.0, 0.0, &
                   0.0, 0.0, 4.0, 0.0,-4.0,-4.0 /), (/ 6, 6 /) )


  !################################################################
  !
  ! Interfaces
  !
  !################################################################
  interface simplex2Dqint
     ! NOTE: the procedures with indices are helper functions for
     !       different input parameter configurations
     ! general quadratic integration
     module procedure gen_qint
     module procedure gen_qint_0
     module procedure gen_qint_1
     module procedure gen_qint_2
     ! quadratic integration with constant f -> DOS
     module procedure dos_sqr
     module procedure dos_sqr_0
     module procedure dos_sqr_1
     module procedure dos_sqr_2
  end interface


  !################################################################
  !
  ! Member functions/subroutines
  !
  !################################################################
contains

  !==============================================================
  !
  ! Function "quadratic_expansion":
  !
  ! Calculate the 6 coefficients of a quadratic expansion given
  ! the function values at the vertices and the midpoints of a
  ! unit triangle:
  !                q(i) = exp_matrix(i,j)*fvals(j)
  !
  ! The expansion is defined as
  !     f_q = q1 + q2*x + q3*y + q4*x^2 + q5*x*y + q6*y^2
  !
  !==============================================================
  !
  ! INPUT
  !    fvals : array containing the 6 known function values
  !
  ! OUTPUT
  !    array containing the 6 expansion coefficients
  !
  !==============================================================
  function quadratic_expansion(fvals)
    double precision, dimension(6), intent(in) :: fvals
    double precision, dimension(6)             :: quadratic_expansion

    ! It's just a multiplication with a constant matrix
    quadratic_expansion = matmul(exp_matrix, fvals)
  end function quadratic_expansion



  !==============================================================
  !
  ! Function "integrate":
  !
  ! Evaluate the six integrals V_1, ..., V_6
  !
  !
  !          /
  !         |
  ! V_i =   | m_i*delta(E - e(x,y)) dxdy
  !         |
  !        /
  !        simplex
  !
  ! or
  !
  !          /
  !         |
  ! U_i =   | m_i*theta(E - e(x,y)) dxdy
  !         |
  !        /
  !        simplex
  !
  ! with m_i = 1, x, y, x^2, xy, y^2
  !
  !==============================================================
  !
  ! INPUT
  !    erg      : (vector) the energy values E for which the
  !               integral should be evaluated
  !    coeff    : coefficients of the quadratic expansion
  !               of the band dispersion (of a single band)
  !    inttype  : the type of the integrand (DELTA, THETA)
  !    flag     : if .TRUE. all integrals V_1 - V_6 are evaluated
  !               if .FALSE. only V_1 is calculated
  !    
  ! OUTPUT
  !     the 6 integrals in an array, where
  !       size(V,2) = size(erg,1) and
  !       size(V,1) = 6
  !
  !==============================================================
  function integrate(erg, coeff, inttype, flag)
    double precision, dimension(:), intent(in) :: erg
    double precision, dimension(6), intent(in) :: coeff
    integer, intent(in)                        :: inttype
    logical, intent(in), optional              :: flag
    double precision, dimension(6,size(erg))   :: integrate

    ! will represent the type of the (quadratic) energy curve
    type(qcurve) :: ergcurve
    ! the inverse of the linear transformation that leads from
    ! the general quadratic curve to one of the reduced forms.
    !It's needed to calculate the vertices of the new triangle
    type(ltrafo) :: invtrsf

    ! the jacobian of the parametrisation of the iso-energy curve
    double precision :: jacobian

    ! the integrals in the transformed triangle
    double precision, allocatable, dimension(:) :: v1, v2, v3, v4, v5, v6

    ! the roots of the equations fot calc of integration limits
    ! we initialize them as HUGE for convenience
    double precision, dimension(6) :: roots
    double precision, dimension(6) :: rootsides ! the corresponding sides
    double precision, dimension(2) :: tmproots

    ! midpoint between to roots
    type(point) :: midpoint
    ! position of midpoint
    !  +1 : inside
    !   0 : on boundary
    !  -1 : outside
    integer     :: midp_pos

    ! the global minimum and maximum of the energy inside the
    ! (transformed) triangle
    double precision :: emin, emax

    ! number of integration domains
    integer :: limitpairs
    ! the limits for the integration as pairs. There should be
    ! maximal 3 pairs, the number is not the same for each
    ! curve type, the first dimension is the pair index
    double precision, dimension(4,2) :: limits
    ! the corresponding triangle sides on which the limits ar
    integer, dimension(4,2) :: limitsides

    ! number of linesegments
    integer :: numofsegments
    ! for the theta function case we need also straight line
    ! segments
    type(point), dimension(4,2) :: linesegments

    ! the transformed triangle (ntr = newtriangle)
    type(triangle) :: ntr
    ! parametrisation of triangle sides as type(line)
    type(line), dimension(3) :: l

    ! for convenience, the transformation matrix and vector
    double precision :: a11, a12, a21, a22, b1, b2

    ! loop variables
    integer :: i,j

    ! vertex index, temp variable
    integer :: ip

    ! some temporary variables
    double precision               :: tmp1, tmp2, tmp3, &
                                      tmp4, tmp5, tmp6
    double precision, dimension(2) :: tmpvec
    integer                        :: tmpn, tmpm
    integer, dimension(3)          :: tmpints
    type(point)                    :: tmpp

    ! flag; if .TRUE. all integrals are evaluated
    logical :: calcall

    ! a small number
    double precision, parameter :: EPS = 1D-9


    !///////////!
    !/  BEGIN  /!
    !///////////!

    ! get value of flag (if present)
    if (present(flag)) then
       calcall = flag
    else
       calcall = .true.
    end if

    ! initialize
    roots = huge(1.0D0)

    ! allocate and initialize all needed integrals to zero
    allocate(v1(size(erg,1)))
    v1 = 0.0D0
    if (calcall) then
       allocate(v2(size(erg,1)))
       allocate(v3(size(erg,1)))
       allocate(v4(size(erg,1)))
       allocate(v5(size(erg,1)))
       allocate(v6(size(erg,1)))
       v2 = 0.0D0
       v3 = 0.0D0
       v4 = 0.0D0
       v5 = 0.0D0
       v6 = 0.0D0
    end if


    ! get global energy minimum and maximum inside triangle
    ! (including boundary!). If the desired energy range doesn't
    ! cover the energy range of the triangle we can return
    ! immediately
    call minmax(coeff, emin, emax)
    if (minval(erg) >= emax) then
       ! the triangle lies entirely in a region with E > e(x,y)
       ! so the delta integrals results to zero, the theta integrals
       ! to the integrals over the whole triangle
       select case (inttype)
       case (THETA)
          integrate(1,:) = 0.5D0
          if (calcall) then
 ! TODO
             integrate(2,:) = v2
             integrate(3,:) = v3
             integrate(4,:) = v4
             integrate(5,:) = v5
             integrate(6,:) = v6
          end if
          return
       case default
          integrate(1,:) = v1
          if (calcall) then
             integrate(2,:) = v2
             integrate(3,:) = v3
             integrate(4,:) = v4
             integrate(5,:) = v5
             integrate(6,:) = v6
          end if
          return
       end select
    end if
    ! the triangle lies entirely in a region with E < e(x,y)
    ! so all integrals evaluate to 0
    if (maxval(erg) <= emin) then
       integrate(1,:) = v1
       if (calcall) then
          integrate(2,:) = v2
          integrate(3,:) = v3
          integrate(4,:) = v4
          integrate(5,:) = v5
          integrate(6,:) = v6
       end if
       return
    end if

    ! determine the type of the quadratic curve represented by
    ! the coefficients coeff
    ergcurve = determine_curve_type(coeff)

    ! transform the triangle according to the transformation
    ! in ergcurve
    invtrsf = invert(ergcurve%trsf)
    tmpvec = matmul(invtrsf%A, (/ 0, 0/) - ergcurve%trsf%b)
    ntr%p(1)%x = tmpvec(1)
    ntr%p(1)%y = tmpvec(2)
    tmpvec = matmul(invtrsf%A, (/ 1, 0/) - ergcurve%trsf%b)
    ntr%p(2)%x = tmpvec(1)
    ntr%p(2)%y = tmpvec(2)
    tmpvec = matmul(invtrsf%A, (/ 0, 1/) - ergcurve%trsf%b)
    ntr%p(3)%x = tmpvec(1)
    ntr%p(3)%y = tmpvec(2)

    ! rearrange the vertices such that stepping through them in
    ! order 1-2-3 one passes the triangle counter clockwise
    call sort_vertices(ntr)

    ! find the triangle sides using a representation of the form
    !    nx*x + ny*y + c =0
    ! 1. side: p1 - p2
    l(1) = get_line(ntr%p(1), ntr%p(2))
    ! 2. side: p2 - p3
    l(2) = get_line(ntr%p(2), ntr%p(3))
    ! 3. side: p3 - p1
    l(3) = get_line(ntr%p(3), ntr%p(1))

    ! NOTE: We have to iterate through all energy values
    !       in erg. For each of these values we check if
    !       it is inside the energy range that is covered
    !       by the triangle and the current band.
    !       If it is inside the range we are assured that
    !       there must be some roots


    !----------------------------------------------------------------------!
    !                                                                      !
    ! The further procedure depends on curve type!                         !
    !                                                                      !
    !----------------------------------------------------------------------!
    select case (ergcurve%curvetype)
    case (qELLIPSE)
       ! the jacobian of the parametrisation
       jacobian = 0.5D0/sqrt((ergcurve%coeff(4)*ergcurve%coeff(6)))

       ! iteration through all energy values
       do i=1,size(erg,1)
          if ((erg(i) > emin) .and. (erg(i) < emax)) then
             ! a trigonometric equ needs to be solved for each side
             ! nx*a*cos(u) + ny*b*sin(u) = c with
             ! a = sqrt((E - q1)/q4), b = sqrt((E - q1)/q6)
             ! We consider the interval -pi <= u < pi
             tmp1 = sqrt((erg(i) - ergcurve%coeff(1))/ergcurve%coeff(4))
             ! tmp1 = a
             tmp2 = sqrt((erg(i) - ergcurve%coeff(1))/ergcurve%coeff(6))
             ! tmp2 = b
             ! 1. side
             call trroots(tmp1*l(1)%nx, tmp2*l(1)%ny, l(1)%c, tmproots, tmpn)
             ! only look at existent roots
             do j=1,tmpn
                roots(j) = tmproots(j)
                rootsides(j) = 1  ! this was found on triangle side 1
             end do
             tmpm = tmpn  ! total number of roots up to now
             ! 2. side
             call trroots(tmp1*l(2)%nx, tmp2*l(2)%ny, l(2)%c, tmproots, tmpn)
             ! only look at existent roots
             do j=1,tmpn
                roots(j+tmpm) = tmproots(j)
                rootsides(j+tmpm) = 2  ! this was found on triangle side 2
             end do
             tmpm = tmpm + tmpn  ! total number of roots up to now
             ! 3. side
             call trroots(tmp1*l(3)%nx, tmp2*l(3)%ny, l(3)%c, tmproots, tmpn)
             ! only look at existent roots
             do j=1,tmpn
                roots(j+tmpm) = tmproots(j)
                rootsides(j+tmpm) = 3  ! this was found on triangle side 3
             end do
             tmpm = tmpm + tmpn  ! total number of roots up to now

             ! sort the roots
             call sort(roots(1:tmpm), rootsides(1:tmpm))

             ! now build pairs, but consider only those for which the
             ! midpoint lies inside the triangle
             ! NOTE 1: We have maximally 6 roots and three pairs
             ! NOTE 2: As an ellipse is closed, N-1 could also be a pair
             ! NOTE 3: Special tests are nedded as the ellipse could
             !         lie entirely inside the triangle

             limitpairs = 0

             if (tmpm >= 2) then
                ! we have at least one integration domain
                ! loop through all possible pairs
                do j=1,tmpm-1
                   midpoint%x = tmp1*cos(0.5D0*(roots(j)+roots(j+1)))
                   midpoint%y = tmp2*sin(0.5D0*(roots(j)+roots(j+1)))
                   ip = is_inside(midpoint, ntr)
                   if (ip > 0) then
                      ! midpoint is inside the triangle, so it's a regular
                      ! integration domain
                      limitpairs = limitpairs + 1
                      limits(limitpairs,:) = (/ roots(j), roots(j+1) /)
                      limitsides(limitpairs,:) = &
                                        (/ rootsides(j), rootsides(j+1) /)
                   ! This is a very special case, and the code fails if a
                   ! triangle vertex lies on the ellipse
                   !elseif (ip == 0) then
                      ! it's only a regular domain if the limits don't lie
                      ! on the same triangle side
                   !   if (rootsides(j) /= rootsides(j+1)) then
                   !      limitpairs = limitpairs + 1
                   !      limits(limitpairs,:) = (/ roots(j), roots(j+1) /)
                   !      limitsides(limitpairs,:) = &
                   !                        (/ rootsides(j), rootsides(j+1) /)
                   !   end if
                   end if
                end do
                ! now check also the segment 1-N
                ! NOTE: we have to add pi to the angle to end up in the
                ! right segment
                midpoint%x = tmp1*cos(0.5D0*(roots(1)+roots(tmpm)) + PI)
                midpoint%y = tmp2*sin(0.5D0*(roots(1)+roots(tmpm)) + PI)
                ip = is_inside(midpoint, ntr)
                if (ip > 0) then
                   ! midpoint is inside the triangle, so it's a regular
                   ! integration domain
                   limitpairs = limitpairs + 1
                   limits(limitpairs,:) = &
                                (/ roots(tmpm), roots(1) + 2.0D0*PI /)
                   limitsides(limitpairs,:) = &
                                (/ rootsides(tmpm), rootsides(1) /)
                 ! see above
                 !elseif (ip == 0) then
                    ! it's only a regular domain if the limits don't lie
                    ! on the same triangle side
                 !   if (rootsides(tmpm) /= rootsides(1)) then
                 !      limitpairs = limitpairs + 1
                 !      limits(limitpairs,:) = &
                 !               (/ roots(tmpm), roots(1) + 2.0D0*PI /)
                 !      limitsides(limitpairs,:) = &
                 !               (/ rootsides(tmpm), rootsides(1) /)
                 !   end if
                end if
             else
                ! in both cases (tmpn = 0,1) the ellipse lies either
                ! entirely inside or outside the triangle
                ! check at two points (to exclude the hypothetic case
                ! where the check point coincides with a tangent point)
                ! check at u=0 and u=pi
                midpoint%x = tmp1
                midpoint%y = 0.0D0
                if (is_inside(midpoint, ntr) > 0) then
                   ! check point is truly inside the triangle, the
                   ! whole ellipse is an integration domain
                   limitpairs = 1
                   limits(limitpairs,:) = (/ 0.0D0, 2*PI /)
                elseif (is_inside(midpoint, ntr) < 0) then
                   ! the ellipse lies outside of our triangle
                   limitpairs = 0
                else
                   ! check a second point
                   midpoint%x = 0.0D0
                   midpoint%y = tmp2
                   if (is_inside(midpoint, ntr) > 0) then
                      ! check point is inside the triangle, so the
                      ! whole ellipse is an integration domain
                      limitpairs = 1
                      limits(limitpairs,:) = (/ 0.0D0, 2*PI /)
                   end if
                end if
             end if

             !-------------------------------------------------------!
             ! Now integrate                                         !
             !-------------------------------------------------------!
             ! which type of integral is asked?
             select case (inttype)
             case (THETA)
                !----------------------------------------------------!
                ! Solve for theta function                           !
                !----------------------------------------------------!
                ! Up to now we only know the ellipse segments.
                ! We have to find the straight line segments.
                ! NOTE: Because of 'limitsides' we know for 
                !       each pair of limits the corresponding
                !       triangle side.
                ! We will integrate over the surface inside the
                ! triangle and then check, if inside or outside
                ! should be returned (outside = triangle - inside).
                numofsegments = 0
                select case (limitpairs)
                case (0)
                   ! the ellipse is outside the triangle
                case (1)
                   ! 1 ellipse segment and 0/1/2/3 line segments
                   if ((limits(1,2) - limits(1,1)) /= 2*PI) then
                      ! We have at least 1 line segment
                      if (limitsides(1,1) == limitsides(1,2)) then
                         ! The ellipse is cut in half and we have 
                         ! exactly one line segment
                         numofsegments = 1
                         tmpp = point(tmp1*cos(limits(1,2)), &
                                      tmp2*sin(limits(1,2)))
                         linesegments(1,1) = tmpp
                         tmpp = point(tmp1*cos(limits(1,1)), &
                                      tmp2*sin(limits(1,1)))
                         linesegments(1,2) = tmpp
                      else
                      ip = sides_to_vertex(limitsides(1,1),limitsides(1,2))
                      if (is_on_the_left(ntr%p(ip), &
                                ergcurve%coeff - (/ erg(i),(0.0D0,j=1,5) /), &
                                qELLIPSE) == 1) then
                         ! This vertex lies inside the ellipse, therefore we get
                         ! two line segments
                         numofsegments = 2
                         tmpp = point(tmp1*cos(limits(1,1)), &
                                      tmp2*sin(limits(1,1)))
                         linesegments(1,:) = (/ ntr%p(ip), tmpp /)
                         tmpp = point(tmp1*cos(limits(1,2)), &
                                      tmp2*sin(limits(1,2)))
                         linesegments(2,:) = (/ tmpp, ntr%p(ip) /)
                      else
                         ! The two other vertices are inside, so we have three
                         ! line segments
                         numofsegments = 3
                         select case (ip)
                         case (1)
                            ! points 2 and 3 build segment
                            linesegments(1,:) = (/ ntr%p(2), ntr%p(3) /)
                         case (2)
                            ! points 3 and 1 build segment
                            linesegments(1,:) = (/ ntr%p(3), ntr%p(1) /)
                         case (3)
                            ! points 1 and 2 build segment
                            linesegments(1,:) = (/ ntr%p(1), ntr%p(2) /)
                         end select
                         ! now the other two segments
                         tmpp = point(tmp1*cos(limits(1,1)), &
                                      tmp2*sin(limits(1,1)))
                         linesegments(2,:) = (/ linesegments(1,2), tmpp /)
                         tmpp = point(tmp1*cos(limits(1,2)), &
                                      tmp2*sin(limits(1,2)))
                         linesegments(3,:) = (/ tmpp, linesegments(1,1) /)
                      end if
                   end if
                   end if
                case (2)
                   ! Two ellipse segments and three line segments.
                   ! The vertices that lie on the triangle side
                   ! with 2 intersections are outside the ellipse
                   tmpints = 0
                   do j=1,2
                      tmpints(limitsides(j,1)) = tmpints(limitsides(j,1)) + 1 
                      tmpints(limitsides(j,2)) = tmpints(limitsides(j,2)) + 1 
                   end do
                   ! on which triangle side do we have two intersections?
                   do j=1,3
                      if (tmpints(j) == 2) then
                         ! Get the point id inside
                         select case (j)
                         case (1)
                            ip = 3
                         case (2)
                            ip = 1
                         case (3)
                            ip = 2
                         end select
                         !continue
                      end if
                   end do

                   ! Now we can calmly set up the linesegments
                   numofsegments = 3
                   ! The first segment
                   tmpp = point(tmp1*cos(limits(1,1)), &
                                tmp2*sin(limits(1,1)))
                   linesegments(1,:) = (/ ntr%p(ip), tmpp /)
                   ! The second segment
                   tmpp = point(tmp1*cos(limits(1,2)), &
                                tmp2*sin(limits(1,2)))
                   linesegments(2,1) = tmpp
                   tmpp = point(tmp1*cos(limits(2,1)), &
                                tmp2*sin(limits(2,1)))
                   linesegments(2,2) = tmpp
                   ! The third segment
                   tmpp = point(tmp1*cos(limits(2,2)), &
                                tmp2*sin(limits(2,2)))
                   linesegments(3,:) = (/ tmpp, ntr%p(ip) /)

                case (3)
                   ! three ellipse segments and three line segments
                   ! All triangle points lie outside
                   numofsegments = 3
                   ! The first segment
                   tmpp = point(tmp1*cos(limits(1,2)), &
                                tmp2*sin(limits(1,2)))
                   linesegments(1,1) = tmpp
                   tmpp = point(tmp1*cos(limits(2,1)), &
                                tmp2*sin(limits(2,1)))
                   linesegments(1,2) = tmpp
                   ! The second segment
                   tmpp = point(tmp1*cos(limits(2,2)), &
                                tmp2*sin(limits(2,2)))
                   linesegments(2,1) = tmpp
                   tmpp = point(tmp1*cos(limits(3,1)), &
                                tmp2*sin(limits(3,1)))
                   linesegments(2,2) = tmpp
                   ! The third segment
                   tmpp = point(tmp1*cos(limits(3,2)), &
                                tmp2*sin(limits(3,2)))
                   linesegments(3,1) = tmpp
                   tmpp = point(tmp1*cos(limits(1,1)), &
                                tmp2*sin(limits(1,1)))
                   linesegments(3,2) = tmpp
                   
                end select
                ! Now we can integrate
                ! Do the ellipse segments. We do them first because we will
                ! multiply the result by a jacobian and other stuff which
                ! should not be done for the line elements.
                do j=1,limitpairs
                   v1(i) = v1(i) + (limits(j,2) - limits(j,1))
                end do
                v1(i) = v1(i)*jacobian*(erg(i) - ergcurve%coeff(1))
                ! Now the line segments
                do j=1,numofsegments
                   v1(i) = v1(i) + line_int(linesegments(j,1),linesegments(j,2))
                end do
   if (v1(i) < 0.0D0) then
      print *,""
      print *," **** v1(i) =",v1(i),"< 0 ****"
      print *," E =",erg(i)
      print *," q_i =",ergcurve%coeff(1),ergcurve%coeff(4),ergcurve%coeff(6)
      print *," The new triangle:"
      print *,"     p1:",ntr%p(1)
      print *,"     p2:",ntr%p(2)
      print *,"     p3:",ntr%p(3)
      print *," (original, reconstructed:"
      tmpvec = (/ ntr%p(1)%x, ntr%p(1)%y /)
      print *,"     p1:",matmul(ergcurve%trsf%A,tmpvec)+ergcurve%trsf%b
      tmpvec = (/ ntr%p(2)%x, ntr%p(2)%y /)
      print *,"     p2:",matmul(ergcurve%trsf%A,tmpvec)+ergcurve%trsf%b
      tmpvec = (/ ntr%p(3)%x, ntr%p(3)%y /)
      print *,"     p3:",matmul(ergcurve%trsf%A,tmpvec)+ergcurve%trsf%b,")"
      print *," Linesegments:"
      do j=1,numofsegments
        print *,"     p1",linesegments(j,1)
        print *,"     p2",linesegments(j,2)
        print *,"     int =",line_int(linesegments(j,1),linesegments(j,2))
      end do
      print *," Ellipse segments:"
      do j=1,limitpairs
        print *,"     u1,u2",limits(j,1),limits(j,2)
        print *,"     int =",(limits(j,2) - limits(j,1))*jacobian*(erg(i) - ergcurve%coeff(1))
      end do
   end if

  !HERE
                if (calcall) then
                   ! First the ellipse segments (c.f. above)
                   do j=1,limitpairs
  !TODO
                      v2(i) = v2(i) + (limits(j,2) - limits(j,1))
                   end do
                   v2(i) = v2(i)*jacobian*(erg(i) - ergcurve%coeff(1))
                   ! Now the line segments
                   do j=1,numofsegments
                      v2(i) = v2(i) + line_int_x(linesegments(j,1),linesegments(j,2))
                      v3(i) = v3(i) + line_int_y(linesegments(j,1),linesegments(j,2))
                      v4(i) = v4(i) + line_int_xx(linesegments(j,1),linesegments(j,2))
                      v5(i) = v5(i) + line_int_xy(linesegments(j,1),linesegments(j,2))
                      v6(i) = v6(i) + line_int_yy(linesegments(j,1),linesegments(j,2))
                   end do
                end if
                !-----------------!
                ! case THETA done !
                !-----------------!

             case default
                !----------------------------------------------------!
                ! solve for delta function                           !
                !----------------------------------------------------!
                do j=1,limitpairs
                   v1(i) = v1(i) + (limits(j,2) - limits(j,1))
                end do
                v1(i) = v1(i)*jacobian

                ! the others are only calculated when requested
                if (calcall) then
                   do j=1,limitpairs
                      tmp3 = sin(limits(j,1))
                      tmp4 = sin(limits(j,2))
                      tmp5 = cos(limits(j,1))
                      tmp6 = cos(limits(j,2))
                      v2(i) = v2(i) + tmp4 - tmp3
                      v3(i) = v3(i) + tmp6 - tmp5
                      v4(i) = v4(i) + tmp4*tmp6 - tmp3*tmp5
                      v5(i) = v5(i) + tmp6*tmp6 - tmp5*tmp5
                   end do
                   v6(i) = v4(i)
                   v2(i) = v2(i)*tmp1*jacobian
                   v3(i) = -v3(i)*tmp2*jacobian
                   v4(i) = 0.5D0*(v1(i) + v4(i)*jacobian)*(tmp1**2)
                   v5(i) = -0.5D0*v5(i)*tmp1*tmp2*jacobian
                   v6(i) = 0.5D0*(v1(i) - v6(i)*jacobian)*(tmp2**2)
                end if
                !-----------------!
                ! case DELTA done !
                !-----------------!
             end select

          end if
          ! if we calculated for THETA, then we must check if
          ! the outer or inner region should be returned
          ! Also we have to check, if erg(i) > max(eps)
          if (inttype == THETA) then
             if ((ergcurve%coeff(4) < 0.0D0) .or. (erg(i) >= emax)) then
                ! It's a maximum, so E - eps(x,y) is positiv outside
                ! We have to calculate the result for a Triangle
                ! The first integral over a triangle is the triangle
                ! area. Because we started with a unit triangle and
                ! applied just trafos with det = 1, the area is still 0.5
 ! TODO: check if erg - q1 is always negativ for this case!
                v1(i) = 0.5D0 + v1(i)
                if (calcall) then
 ! TODO
                end if
             end if
          end if
       end do
       !-----------------------------------------------------------
       ! End of 'qELLIPSE'
       !-----------------------------------------------------------

    case (qHYPERBOLA)
       ! the jacobian of the parametrisation
       jacobian = abs(1.0D0/ergcurve%coeff(5))

       ! iteration
       do i=1,size(erg,1)
          if ((erg(i) > emin) .and. (erg(i) < emax)) then
             ! a quadratic equ needs to be solved for each side
             ! nx*u^2 - c*u + ny*(E - q1)/q5 = 0
             tmp1 = (erg(i) - ergcurve%coeff(1))/ergcurve%coeff(5)
             ! 1. side
             call qroots(l(1)%nx, l(1)%c, l(1)%ny*tmp1, tmproots, tmpn)
             roots(1:2) = tmproots
             ! 2. side
             call qroots(l(2)%nx, l(2)%c, l(2)%ny*tmp1, tmproots)
             roots(3:4) = tmproots
             ! 3. side
             call qroots(l(3)%nx, l(3)%c, l(3)%ny*tmp1, tmproots)
             roots(5:6) = tmproots

             ! u = 0 is not a real root, but can be introduced
             ! artificially because of a multiplication by u
             do j=1,6
                if (abs(roots(j)) < EPS) then
                   roots(j) = huge(1.0D0)
                end if
             end do
             ! sort the roots
             call sort(roots(1:6))

             ! how many roots are there? (remember, HUGE is no valid
             ! root, and all HUGEs are at the high end of roots(:))
             do j=1,6
                if (roots(j) == huge(1.0D0)) then
                   tmpm = j-1
                   exit
                end if
             end do

             ! now build pairs, but consider only those for which the
             ! midpoint lies inside the triangle
             ! NOTE: We have maximally 4 roots, that means two pairs
             !       (there are two distinct branches u < 0, u > 0)

             limitpairs = 0

             do j=1,tmpm-1
                ! the two roots must have the same sign to build a pair
                ! as different signs belong to different branches of the
                ! hyperbola
                if (sign(1.0D0, roots(j)) == sign(1.0D0, roots(j+1))) then
                   midpoint%x = 0.5D0*(roots(j) + roots(j+1))
                   midpoint%y = tmp1/midpoint%x
                   if (is_inside(midpoint, ntr) >= 0) then
                      ! midpoint is inside the triangle, so it's a regular
                      ! integration domain
                      limitpairs = limitpairs + 1
                      if (roots(j) < 0.0D0) then
                         ! we change integration direction
                         limits(limitpairs,:) = (/ -roots(j+1), -roots(j) /)
                      else
                         limits(limitpairs,:) = (/ roots(j), roots(j+1) /)
                      end if
                   end if
                end if
             end do

             !---------------!
             ! Now integrate !
             !---------------!
             ! NOTE: The sign of the limits enters in the calculation,
             !       but this was already considered above (swapping the
             !       limits for negativ case and taking the absolute value)
             do j=1,limitpairs
                v1(i) = v1(i) + (log(limits(j,2)) - &
                     log(limits(j,1)))
             end do
             v1(i) = v1(i)*jacobian

             ! the others are only calculated when requested
             if (calcall) then
                do j=1,limitpairs
                   v2(i) = v2(i) + (limits(j,2) - limits(j,1))
                   v3(i) = v3(i) - &
                        tmp1*(1/limits(j,2) - 1/limits(j,1))
                   v4(i) = v4(i) + &
                        0.5D0*(limits(j,2)**2 - limits(j,1)**2)
                end do
                v2(i) = v2(i)*jacobian
                v3(i) = v3(i)*jacobian
                v4(i) = v4(i)*jacobian
                v5(i) = v1(i)*tmp1
                v6(i) = -v4(i)*tmp1**2
             end if
          end if
       end do
       !-----------------------------------------------------------
       ! End of 'qHYPERBOLA'
       !-----------------------------------------------------------

    case (qPARABOLA)
       ! the jacobian of the parametrisation
       jacobian = abs(1.0D0/ergcurve%coeff(3))

       ! iteration
       do i=1,size(erg,1)
          if ((erg(i) > emin) .and. (erg(i) < emax)) then
             ! a quadratic equ needs to be solved for each side
             ! -ny*(q4/q3)*u^2 + nx*u + ny*(E - q1)/q3 - c = 0
             tmp1 = (erg(i) - ergcurve%coeff(1))/ergcurve%coeff(3)
             tmp2 = ergcurve%coeff(4)/ergcurve%coeff(3)
             ! 1. side
             call qroots(-l(1)%ny*tmp2, l(1)%nx, l(1)%ny*tmp1 + l(1)%c, tmproots)
             roots(1:2) = tmproots
             ! 2. side
             call qroots(-l(2)%ny*tmp2, l(2)%nx, l(2)%ny*tmp1 + l(2)%c, tmproots)
             roots(3:4) = tmproots
             ! 3. side
             call qroots(-l(3)%ny*tmp2, l(3)%nx, l(3)%ny*tmp1 + l(3)%c, tmproots)
             roots(5:6) = tmproots

             ! sort the roots
             call sort(roots(1:6))

             ! now build pairs, but consider only those for which the
             ! midpoint lies inside the triangle
             limitpairs = 0

             ! we loop through all possible pairs
             do j=1,5
                midpoint%x = 0.5D0*(roots(j) + roots(j+1))
                midpoint%y = tmp1 - tmp2*(midpoint%x)**2
                if (is_inside(midpoint, ntr) >= 0) then
                   ! midpoint is inside the triangle, so it's a regular
                   ! integration domain
                   limitpairs = limitpairs + 1
                   limits(limitpairs,:) = (/ roots(j), roots(j+1) /)
                end if
             end do

             !---------------!
             ! Now integrate !
             !---------------!
             do j=1,limitpairs
                v1(i) = v1(i) + (limits(j,2) - limits(j,1))
             end do
             v1(i) = v1(i)*jacobian

             ! the others are only calculated when requested
             if (calcall) then
                do j=1,limitpairs
                   v2(i) = v2(i) + 0.50D0*(limits(j,2)**2 - &
                        limits(j,1)**2)
                   v4(i) = v4(i) + (limits(j,2)**3 - &
                        limits(j,1)**3)/3.0D0
                   v5(i) = v5(i) + 0.25D0*(limits(j,2)**4 - &
                        limits(j,1)**4)
                   v6(i) = v6(i) + 0.20D0*(limits(j,2)**5 - &
                        limits(j,1)**5)
                end do
                v2(i) = v2(i)*jacobian
                v4(i) = v4(i)*jacobian
                v3(i) = tmp1*v1(i) - tmp2*v4(i)
                v5(i) = tmp1*v2(i) - tmp2*v5(i)*jacobian
                v6(i) = (tmp1**2)*v1(i) - 2.0D0*tmp1*tmp2*v4(i) + &
                     (tmp2**2)*v6(i)*jacobian
             end if

          end if
       end do
       !-----------------------------------------------------------
       ! End of 'qPARABOLA'
       !-----------------------------------------------------------

    case (qLINE)
       ! the jacobian of the parametrisation
       jacobian = abs(1.0D0/ergcurve%coeff(2))

       ! iteration
       do i=1,size(erg,1)
          if ((erg(i) >= emin) .and. (erg(i) <= emax)) then
             ! a linear equ needs to be solved for each side
             ! ny*u + nx*(E - q1)/q2 - c = 0
             ! NOTE: as ny may be 0
             ! NOTE: tmp1 will be used also a few lines below
             tmp1 = (erg(i) - ergcurve%coeff(1))/ergcurve%coeff(2)

             ! 1. side
             call qroots(0.0D0, l(1)%ny, l(1)%nx*tmp1 + l(1)%c, tmproots)
             roots(1) = tmproots(1) ! there can be only one root
             ! 2. side
             call qroots(0.0D0, l(2)%ny, l(2)%nx*tmp1 + l(2)%c, tmproots)
             roots(2) = tmproots(1) ! there can be only one root
             ! 3. side
             call qroots(0.0D0, l(3)%ny, l(3)%nx*tmp1 + l(3)%c, tmproots)
             roots(3) = tmproots(1) ! there can be only one root

             ! sort the roots in ascending order
             call sort(roots(1:3))

             ! now build pairs, but consider only those for which the
             ! midpoint lies inside the triangle
             ! NOTE: There will be maximally 1 pair
             ! NOTE: If midpoint lies on the boundary, the integral
             !       will be divided by 2, as the adjacent triangle
             !       will contribute the same amount
             ! TODO: Is this true??? I think it should..

             limitpairs = 0
             ! check for (quasi-)coincidence
             if (roots(2) - roots(1) > EPS) then
                midpoint%x = tmp1
                midpoint%y = 0.5D0*(roots(1) + roots(2))
                midp_pos = is_inside(midpoint, ntr)
                if (midp_pos >= 0) then
                   ! midpoint is inside the triangle, so it's a regular
                   ! integration domain (and also the only one)
                   limits(1,:) = (/ roots(1), roots(2) /)
                   limitpairs = 1
                else
                   midpoint%y = 0.5D0*(roots(2) + roots(3))
                   midp_pos = is_inside(midpoint, ntr)
                   if (midp_pos >= 0) then
                      ! midpoint is inside the triangle, so it's a regular
                      ! integration domain (and also the only one)
                      limits(1,:) = (/ roots(2), roots(3) /)
                      limitpairs = 1
                   end if
                end if
             else
                ! roots 1 and 2 are (quasi-)coincident and we don't
                ! consider them
                midpoint%x = tmp1
                midpoint%y = 0.5D0*(roots(2) + roots(3))
                midp_pos = is_inside(midpoint, ntr)
                if (midp_pos >= 0) then
                   ! midpoint is inside the triangle, so it's a regular
                   ! integration domain (and also the only one)
                   limits(1,:) = (/ roots(2), roots(3) /)
                   limitpairs = 1
                end if
             end if

             !---------------!
             ! Now integrate !
             !---------------!
             if (limitpairs == 1) then
                v1(i) = (limits(1,2) - limits(1,1))*jacobian
                ! if integration domain lies on a triangle side,&
                ! divide by 2
                if (midp_pos == 0) then
                   v1(i) = 0.5D0*v1(i)
                end if

                ! the others are only calculated when requested
                if (calcall) then
                   v2(i) = v1(i)*tmp1
                   v3(i) = 0.5D0*(limits(1,2)**2 - &
                        limits(1,1)**2)*jacobian
                   v4(i) = v2(i)*tmp1
                   v5(i) = v3(i)*tmp1
                   v6(i) = (limits(1,2)**3 - &
                        limits(1,1)**3)*jacobian/3.0D0
                   ! if integration domain lies on a triangle side,&
                   ! divide by 2
                   if (midp_pos == 0) then
                      v2(i) = 0.5D0*v2(i)
                      v3(i) = 0.5D0*v3(i)
                      v4(i) = 0.5D0*v4(i)
                      v5(i) = 0.5D0*v5(i)
                      v6(i) = 0.5D0*v6(i)
                   end if
                end if
             end if

          end if
       end do
       !-----------------------------------------------------------
       ! End of 'qLINE'
       !-----------------------------------------------------------

    case (qDEG_LINES)
       ! iteration
       do i=1,size(erg,1)
          if ((erg(i) > emin) .and. (erg(i) < emax)) then
             ! the jacobian for our parametrisation
             jacobian = 0.5D0/sqrt(ergcurve%coeff(4)*(erg(i) - &
                  ergcurve%coeff(1)))

             ! a linear equ needs to be solved for each side
             ! ny*u +/- nx*sqrt((E - q1)/q4) - c = 0
             !!
             !! NOTE: for qDEG_LINES we have two lines, so 
             !!       everything has to be done two times
             tmp1 = sqrt((erg(i) - ergcurve%coeff(1))/ergcurve%coeff(4))

             !-------------------------------------!
             ! loop through the two parallel lines !
             !-------------------------------------!
             do j=1,2
                tmp1 = -1.0D0*tmp1
                ! 1. side
                call qroots(0.0D0, l(1)%ny, l(1)%nx*tmp1 + l(1)%c, tmproots)
                roots(1) = tmproots(1) ! there can be only one root
                ! 2. side
                call qroots(0.0D0, l(2)%ny, l(2)%nx*tmp1 + l(2)%c, tmproots)
                roots(2) = tmproots(1) ! there can be only one root
                ! 3. side
                call qroots(0.0D0, l(3)%ny, l(3)%nx*tmp1 + l(3)%c, tmproots)
                roots(3) = tmproots(1) ! there can be only one root

                ! sort the roots in ascending order
                call sort(roots(1:3))

                ! now build pairs, but consider only those for which the
                ! midpoint lies inside the triangle
                ! NOTE: There will be maximally 1 pair
                ! NOTE: If midpoint lies on the boundary, the integral
                !       will be divided by 2, as the adjacent triangle
                !       will contribute the same amount

                limitpairs = 0
                ! check for (quasi-)coincidence
                if (roots(2) - roots(1) > EPS) then
                   midpoint%x = tmp1
                   midpoint%y = 0.5D0*(roots(1) + roots(2))
                   midp_pos = is_inside(midpoint, ntr)
                   if (midp_pos >= 0) then
                      ! midpoint is inside the triangle, so it's a regular
                      ! integration domain (and also the only one)
                      limits(1,:) = (/ roots(1), roots(2) /)
                      limitpairs = 1
                   else
                      midpoint%y = 0.5D0*(roots(2) + roots(3))
                      midp_pos = is_inside(midpoint, ntr)
                      if (midp_pos >= 0) then
                         ! midpoint is inside the triangle, so it's a regular
                         ! integration domain (and also the only one)
                         limits(1,:) = (/ roots(2), roots(3) /)
                         limitpairs = 1
                      end if
                   end if
                else
                   ! roots 1 and 2 are (quasi-)coincident and we don't
                   ! consider them
                   midpoint%x = tmp1
                   midpoint%y = 0.5D0*(roots(2) + roots(3))
                   midp_pos = is_inside(midpoint, ntr)
                   if (midp_pos >= 0) then
                      ! midpoint is inside the triangle, so it's a regular
                      ! integration domain (and also the only one)
                      limits(1,:) = (/ roots(2), roots(3) /)
                      limitpairs = 1
                   end if
                end if

                ! if integration domain lies on a triangle side, divide by 2
                if (midp_pos == 0) then
                   tmp2 = 0.5D0
                else
                   tmp2 = 1.0D0
                end if

                !---------------!
                ! Now integrate !
                !---------------!
                if (limitpairs == 1) then
                   v1(i) = v1(i) + (limits(1,2) - limits(1,1))*tmp2

                   ! the others are only calculated when requested
                   if (calcall) then
                      v2(i) = v2(i) + (limits(1,2) - limits(1,1))*tmp1*tmp2
                      v3(i) = v3(i) + 0.5D0*(limits(1,2)**2 - &
                           limits(1,1)**2)*tmp2
                      v4(i) = v4(i) + (limits(1,2) - &
                           limits(1,1))*(tmp1**2)*tmp2
                      v5(i) = v5(i) + 0.5D0*(limits(1,2)**2 - &
                           limits(1,1)**2)*tmp1*tmp2
                      v6(i) = v6(i) + (limits(1,2)**3 - &
                           limits(1,1)**3)/3.0D0*tmp2
                   end if
                end if
             end do

             v1(i) = v1(i)*jacobian
             if (calcall) then
                v2(i) = v2(i)*jacobian
                v3(i) = v3(i)*jacobian
                v4(i) = v4(i)*jacobian
                v5(i) = v5(i)*jacobian
                v6(i) = v6(i)*jacobian
             end if

          end if
       end do
       !-----------------------------------------------------------
       ! End of 'qDEG_LINES'
       !-----------------------------------------------------------


    case (qDEG_DELTA)
       ! curvetype == qDEG_DELTA
       ! We will probably never hit exactly a flat band, so we
       ! don't do anything

       !-----------------------------------------------------------
       ! End of 'qDEG_DELTA'
       !-----------------------------------------------------------
    end select
    ! we have now the V_i integrals in the transformed triangle

    ! transform V_i to original (right unit) triangle
    a11 = ergcurve%trsf%A(1,1)
    a12 = ergcurve%trsf%A(1,2)
    a21 = ergcurve%trsf%A(2,1)
    a22 = ergcurve%trsf%A(2,2)
    b1  = ergcurve%trsf%b(1)
    b2  = ergcurve%trsf%b(2)
    integrate(1,:) = v1
    if (calcall) then
       integrate(2,:) = a11*v2 + a12*v3 + b1*v1
       integrate(3,:) = a21*v2 + a22*v3 + b2*v1
       integrate(4,:) = a11*a11*v4 + a12*a12*v6 + b1*b1*v1 + &
            2*(a11*a12*v5 + a11*b1*v2 + a12*b1*v3)
       integrate(5,:) = a11*a21*v4 + a12*a22*v6 + b1*b2*v1 + &
            (a11*a22 + a21*a12)*v5 + (a11*b2 + a21*b1)*v2 + &
            (a12*b2 + a22*b1)*v3
       integrate(6,:) = a21*a21*v4 + a22*a22*v6 + b2*b2*v1 + &
            2*(a21*a22*v5 + a21*b2*v2 + a22*b2*v3)
       ! deallocate allocated arrays
       deallocate(v2)
       deallocate(v3)
       deallocate(v4)
       deallocate(v5)
       deallocate(v6)
    end if
    deallocate(v1)

  end function integrate



  !==============================================================
  !
  ! Function "gen_qint":
  !
  ! Generalized quadratic integration over the simplex.
  ! 'generalized' means that a (quadratic) weighting function
  ! is considered:
  !
  !   --    /
  !   \    |
  !   /    | f_n(x) delta(E - e_n(x)) dxdy
  !   --   |
  !   n   /
  !       simplex
  !
  !
  !    --    /
  !    \    |
  !    /    | f_n(x) theta(E - e(x,y)) dxdy
  !    --   |
  !    n   /
  !        simplex
  !
  !==============================================================
  !
  ! INPUT
  !    erg      : ( vector ) the energy values E for which the
  !               integral should be evaluated
  !    simplex  : ( type(triangle) ) simplex (= triangle),
  !    e_disp   : the band energy values at the vertices and at
  !               the midpoints of the sides of the simplex
  !               e_disp(i,k) = e_k at vertex i (k = bandindex)
  !    func     : the function values at the vertices and at the
  !               midpoints of the sides of the simplex
  !               func(i,k) = f_k at vertex i (k = bandindex)
  !
  ! OUTPUT
  !     the value of the integral for each energy value erg(i)
  !     size(gen_qint,1) = size(erg,1)
  !     size(gen_qint,2) = 1 or 2 depending on integraltype
  !
  !==============================================================
  function gen_qint(erg, simplex, e_disp, func, integraltype)
    double precision, dimension(:),   intent(in) :: erg
    type(triangle),                   intent(in) :: simplex
    double precision, dimension(:,:), intent(in) :: e_disp
    double precision, dimension(:,:), intent(in) :: func
    integer, intent(in), optional                :: integraltype
    double precision, dimension(size(erg,1))     :: gen_qint

    ! the determinant of the transformation from the
    ! original triangle to a unit triangle
    double precision :: det

    ! the coefficients of the quadratic expansion of the
    ! energy dispersion inside the simplex
    double precision, dimension(6) :: qexp_coeff

    ! the integrals V_i
    double precision, dimension(6,size(erg,1)) :: V

    ! the weights w_j
    double precision, dimension(6) :: w

    ! the final 2D integral
    double precision, dimension(size(erg,1)) :: I

    ! the type of the integrand (DELTA, THETA)
    integer :: inttype

    ! loop variables
    integer :: band_id, erg_id

    if (present(integraltype)) then
      inttype = integraltype
    else
      inttype = DELTA
    end if

    !------------------------------------------------!
    ! Initialisations                                !
    !------------------------------------------------!
    ! this is important, as we use
    ! expressions I = I + ...
    I = 0.0D0

    !------------------------------------------------!
    ! Transform the original triangle into the unit  !
    ! triangle. We just need the determinant of the  !
    ! trafo, which equals 2 times the triangle area. !
    !------------------------------------------------!
    ! get the area of the triangle 
    det = 2.0D0*area(simplex)

    !! each energy band needs to be considered seperately,
    !! so loop through all bands
    do band_id=1,size(e_disp, 2)
    
       !---------------------------------------------!
       ! Inside the unit triangle, expand the        !
       ! energy dispersion quadratically.            !
       !---------------------------------------------!
       ! calculate expansion coefficients
       qexp_coeff = quadratic_expansion(e_disp(:,band_id))

       !---------------------------------------------!
       ! Now compute the integrals V_i that are      !
       ! needed for the 2D BZ integral.              !
       !---------------------------------------------!
       V = integrate(erg, qexp_coeff, inttype, .true.)

       !---------------------------------------------!
       ! calculate the weights w_j which are needed  !
       ! to evaluate the 2D BZ integral for the      !
       ! current band:                               !
       !                                             !
       !         I = sum(1..6) w_j*func_j            !
       !                                             !
       !     w_j = sum(1..6) exp_matrix(i,j)*V(i)    !
       !                                             !
       ! NOTE: the transposed of exp_matrix has to   !
       !       be used here.                          !
       !                                             !
       ! NOTE: for each energy erg(i), the w_j are   !
       !       different.                            !
       !---------------------------------------------!
       do erg_id=1,size(erg,1)
          w = matmul(transpose(exp_matrix), V(:,erg_id))
          I(erg_id) = I(erg_id) + dot_product(w,func(:,band_id))
       end do

    end do

    !------------------------------------------------!
    ! Transform back to the original triangle as     !
    ! given to the function: multiply by det         !
    !------------------------------------------------!
    gen_qint = det*I

  end function gen_qint
  !==============================================================
  ! Helper function, for the case where e_disp contains only one
  ! band
  !==============================================================
  function gen_qint_0(erg, simplex, e_disp, func, integraltype)
    double precision, dimension(:),   intent(in) :: erg
    type(triangle),                   intent(in) :: simplex
    double precision, dimension(6),   intent(in) :: e_disp
    double precision, dimension(:,:), intent(in) :: func
    integer, intent(in), optional                :: integraltype
    double precision, dimension(size(erg,1))     :: gen_qint_0

    double precision, dimension(6,1) :: disp

    disp(:,1) = e_disp
    ! call the function that does all the work
    if (present(integraltype)) then
       gen_qint_0 = gen_qint(erg, simplex, disp, func, integraltype)
    else
       gen_qint_0 = gen_qint(erg, simplex, disp, func)
    end if

  end function gen_qint_0
  !==============================================================
  ! Helper function, for the case where e_disp contains only one
  ! band, and vertices is a matrix
  ! vertices must be organized as
  !                       / x1 y1 \
  !                      |  x2 y2  |
  !                       \ x3 y3 /
  !==============================================================
  function gen_qint_1(erg, vertices, e_disp, func, integraltype)
    double precision, dimension(:),   intent(in) :: erg
    double precision, dimension(3,2), intent(in) :: vertices
    double precision, dimension(6),   intent(in) :: e_disp
    double precision, dimension(:,:), intent(in) :: func
    integer, intent(in), optional                :: integraltype
    double precision, dimension(size(erg,1))     :: gen_qint_1

    ! the simplex as type(triangle)
    type(triangle) :: tr

    double precision, dimension(6,1) :: disp

    ! convert vertices to a type(triangle)
    tr = triangle((/ point(vertices(1,1), vertices(1,2)), &
         point(vertices(2,1), vertices(2,2)), &
         point(vertices(3,1), vertices(3,2)) /))

    disp(:,1) = e_disp
    ! call the function that does all the work
    if (present(integraltype)) then
       gen_qint_1 = gen_qint(erg, tr, disp, func, integraltype)
    else
       gen_qint_1 = gen_qint(erg, tr, disp, func)
    end if

  end function gen_qint_1
  !==============================================================
  ! Helper function, for the case where vertices is a matrix
  ! vertices must be organized as
  !                       / x1 y1 \
  !                      |  x2 y2  |
  !                       \ x3 y3 /
  !==============================================================
  function gen_qint_2(erg, vertices, e_disp, func, integraltype)
    double precision, dimension(:),   intent(in) :: erg
    double precision, dimension(3,2), intent(in) :: vertices
    double precision, dimension(:,:), intent(in) :: e_disp
    double precision, dimension(:,:), intent(in) :: func
    integer, intent(in), optional                :: integraltype
    double precision, dimension(size(erg,1))     :: gen_qint_2

    ! the simplex as type(triangle)
    type(triangle) :: tr

    ! convert vertices to a type(triangle)
    tr = triangle((/ point(vertices(1,1), vertices(1,2)), &
         point(vertices(2,1), vertices(2,2)), &
         point(vertices(3,1), vertices(3,2)) /))

    ! call the function that does all the work
    if (present(integraltype)) then
       gen_qint_2 = gen_qint(erg, tr, e_disp, func, integraltype)
    else
       gen_qint_2 = gen_qint(erg, tr, e_disp, func)
    end if

  end function gen_qint_2


  !==============================================================
  !
  ! Function "dos_sqr":
  !
  ! Quadratic integration over the simplex
  ! The result is the DOS contribution of the simplex
  ! in the case of the delta-function
  !
  !   --    /
  !   \    |
  !   /    | delta(E - e_n(x)) dx
  !   --   |
  !   n   /
  !       simplex
  !
  !
  !    --    /
  !    \    |
  !    /    | theta(E - e(x,y)) dxdy
  !    --   |
  !    n   /
  !        simplex
  !
  !==============================================================
  !
  ! INPUT
  !    erg      : (vector) the energy values E for which the
  !               integral should be evaluated
  !    simplex  : ( type(triangle) ) simplex (= triangle),
  !    e_disp   : the band energy values at the vertices and
  !               at the midpoints of the sides
  !               e_disp(i,k) = e_k at vertex i
  !
  ! OUTPUT
  !     the value of the integral(s)
  !     size(dos_sqr,1) = size(erg,1)
  !     size(dos_sqr,2) = 1
  !
  !==============================================================
  function dos_sqr(erg, simplex, e_disp, integraltype)
    double precision, dimension(:),   intent(in) :: erg
    type(triangle),                   intent(in) :: simplex
    double precision, dimension(:,:), intent(in) :: e_disp
    integer, intent(in), optional                :: integraltype
    double precision, dimension(size(erg,1))     :: dos_sqr

    ! the determinant of the transformation from the
    ! original triangle to a unit triangle
    double precision :: det

    ! the coefficients of the quadratic expansion of the
    ! energy dispersion inside the simplex
    double precision, dimension(6) :: qexp_coeff

    ! the integral V_1 
    double precision, dimension(6,size(erg,1)) :: V1

    ! the type of the integrand (DELTA, THETA)
    integer :: inttype

    ! loop variables
    integer :: band_id

    if (present(integraltype)) then
      inttype = integraltype
    else
      inttype = DELTA
    end if


    !------------------------------------------------!
    ! Initialisations                                !
    !------------------------------------------------!
    ! this is important, as we use
    ! expressions V_1 = V_1 + ...
    V1 = 0.0D0

    !------------------------------------------------!
    ! Transform the original triangle into the unit  !
    ! triangle. We just need the determinant of the  !
    ! trafo, which equals 2 times the triangle area. !
    !------------------------------------------------!
    ! get the area of the triangle 
    det = 2.0D0*area(simplex)

    !! each energy band needs to be considered seperately,
    !! so loop through all bands
    do band_id=1,size(e_disp, 2)
    
       !---------------------------------------------!
       ! Inside the unit triangle, expand the        !
       ! energy dispersion quadratically.            !
       !---------------------------------------------!
       ! calculate expansion coefficients
       qexp_coeff = quadratic_expansion(e_disp(:,band_id))

       !---------------------------------------------!
       ! Now compute the integrals V_i that are      !
       ! needed for the 2D BZ integral.              !
       ! In the actual case, where we don't have     !
       ! weighting function, we just need V_1        !
       !---------------------------------------------!
       V1 = V1 + integrate(erg, qexp_coeff, inttype, .false.)
    end do

    !------------------------------------------------!
    ! In this case, where we have no wighting        !
    ! function, we are already done, because         !
    ! V_1 = DOS                                      !
    !------------------------------------------------!

    !------------------------------------------------!
    ! Transform back to the original triangle as     !
    ! given to the function: multiply by det         !
    !------------------------------------------------!
    dos_sqr = det*V1(1,:)

  end function dos_sqr
  !==============================================================
  ! Helper function, for the case where e_disp contains only one
  ! band
  !==============================================================
  function dos_sqr_0(erg, simplex, e_disp, integraltype)
    double precision, dimension(:),   intent(in) :: erg
    type(triangle),                   intent(in) :: simplex
    double precision, dimension(6),   intent(in) :: e_disp
    integer, intent(in), optional                :: integraltype
    double precision, dimension(size(erg,1))     :: dos_sqr_0

    double precision, dimension(6,1) :: disp

    disp(:,1) = e_disp
    ! call the function that does all the work
    if (present(integraltype)) then
       dos_sqr_0 = dos_sqr(erg, simplex, disp, integraltype)
    else
       dos_sqr_0 = dos_sqr(erg, simplex, disp)
    end if

  end function dos_sqr_0
  !==============================================================
  ! Helper function, for the case where e_disp contains only one
  ! band, and vertices is a matrix
  ! vertices must be organized as
  !                       / x1 y1 \
  !                      |  x2 y2  |
  !                       \ x3 y3 /
  !==============================================================
  function dos_sqr_1(erg, vertices, e_disp, integraltype)
    double precision, dimension(:),   intent(in) :: erg
    double precision, dimension(3,2), intent(in) :: vertices
    double precision, dimension(6),   intent(in) :: e_disp
    integer, intent(in), optional                :: integraltype
    double precision, dimension(size(erg,1))     :: dos_sqr_1

    ! the simplex as type(triangle)
    type(triangle) :: tr

    double precision, dimension(6,1) :: disp

    ! convert vertices to a type(triangle)
    tr = triangle((/ point(vertices(1,1), vertices(1,2)), &
         point(vertices(2,1), vertices(2,2)), &
         point(vertices(3,1), vertices(3,2)) /))

    disp(:,1) = e_disp
    ! call the function that does all the work
    if (present(integraltype)) then
       dos_sqr_1 = dos_sqr(erg, tr, disp, integraltype)
    else
       dos_sqr_1 = dos_sqr(erg, tr, disp)
    end if

  end function dos_sqr_1
  !==============================================================
  ! Helper function, for the case where vertices is a matrix
  ! vertices must be organized as
  !                       / x1 y1 \
  !                      |  x2 y2  |
  !                       \ x3 y3 /
  !==============================================================
  function dos_sqr_2(erg, vertices, e_disp, integraltype)
    double precision, dimension(:),   intent(in) :: erg
    double precision, dimension(3,2), intent(in) :: vertices
    double precision, dimension(:,:), intent(in) :: e_disp
    integer, intent(in), optional                :: integraltype
    double precision, dimension(size(erg,1))     :: dos_sqr_2

    ! the simplex as type(triangle)
    type(triangle) :: tr

    ! convert vertices to a type(triangle)
    tr = triangle((/ point(vertices(1,1), vertices(1,2)), &
         point(vertices(2,1), vertices(2,2)), &
         point(vertices(3,1), vertices(3,2)) /))

    ! call the function that does all the work
    if (present(integraltype)) then
       dos_sqr_2 = dos_sqr(erg, tr, e_disp, integraltype)
    else
       dos_sqr_2 = dos_sqr(erg, tr, e_disp)
    end if

  end function dos_sqr_2


  !==============================================================
  !
  ! subroutine "sort":
  !
  ! Sort a vector in ascending order.
  !
  ! NOTE: we use the bubble sort algorithm because our vectors
  !       to be sorted contain only a few items
  !
  !==============================================================
  !
  ! INPUT
  !    vec  : unsorted vector
  !    vec2 : second vector that will be reordered like vec
  !
  ! OUTPUT
  !    vec  : sorted vector
  !    vec2 : ordered like vec
  !
  !==============================================================
  subroutine sort(vec, vec2)
    ! at input:  unsorted
    ! at output: sorted
    double precision, dimension(:), intent(inout)           :: vec
    double precision, dimension(:), intent(inout), optional :: vec2

    ! loop variables
    integer :: i, j

    ! a temporary int variable
    integer :: id

    ! temp variable
    double precision :: temp

    ! a second vector
    double precision, dimension(size(vec)) :: v

    if (present(vec2)) then
       v = vec2
    else
       v = (/ (i,i=1,size(vec)) /)
    end if

    do i=size(vec),1,-1
       do j=1,i-1
          if (vec(j) > vec(i)) then
             temp = vec(j)
             vec(j) = vec(i)
             vec(i) = temp
             temp = v(j)
             v(j) = v(i)
             v(i) = temp
          end if
       end do
    end do

    if (present(vec2)) then
       id = min(size(vec2),size(v))
       vec2(1:id) = v(1:id)
    end if

  end subroutine sort


  function line_int(a, b)
      type(point), intent(in) :: a
      type(point), intent(in) :: b
      double precision        :: line_int

      line_int = 0.5D0*(a%x*b%y - a%y*b%x)
  end function line_int

!TODO
  function line_int_x(a, b)
      type(point), intent(in) :: a
      type(point), intent(in) :: b
      double precision        :: line_int_x

      line_int_x = 0.5D0*(a%x*b%y - a%y*b%x)
  end function line_int_x

  function line_int_y(a, b)
      type(point), intent(in) :: a
      type(point), intent(in) :: b
      double precision        :: line_int_y

      line_int_y = 0.5D0*(a%x*b%y - a%y*b%x)
  end function line_int_y

  function line_int_xx(a, b)
      type(point), intent(in) :: a
      type(point), intent(in) :: b
      double precision        :: line_int_xx

      line_int_xx = 0.5D0*(a%x*b%y - a%y*b%x)
  end function line_int_xx

  function line_int_xy(a, b)
      type(point), intent(in) :: a
      type(point), intent(in) :: b
      double precision        :: line_int_xy

      line_int_xy = 0.5D0*(a%x*b%y - a%y*b%x)
  end function line_int_xy

  function line_int_yy(a, b)
      type(point), intent(in) :: a
      type(point), intent(in) :: b
      double precision        :: line_int_yy

      line_int_yy = 0.5D0*(a%x*b%y - a%y*b%x)
  end function line_int_yy

end module simplex2D

