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
!           Module "geometry"
!
!==============================================================
!
! Author: Matthias Auf der Maur
!
! Defines useful types and functions for geometric objects.
! Currently only point and triangle are defined.
!
!==============================================================


module geometry

  implicit none

  ! Pi
  double precision, parameter :: PI = 3.14159265358979
  
  !==============================================================
  ! a point in the x-y-plane
  !==============================================================
  type :: point
     double precision :: x
     double precision :: y
  end type point

  !==============================================================
  ! a line in the x-y-plane
  ! the line representation is as follows:
  !  line : nx*x + ny*y + c = 0
  !==============================================================
  type :: line
     double precision :: nx
     double precision :: ny
     double precision :: c
  end type line


  !==============================================================
  ! triangle defines a triangle ...
  !==============================================================
  type :: triangle
     type(point), dimension(3) :: p ! the three vertices
  end type triangle


contains

  !==============================================================
  !
  ! function "get_line"
  !
  ! Returns the line definition as type(line) that is determined
  ! by two given points.
  !
  !==============================================================
  !
  ! INPUT
  !    p1 : the first point
  !    p2 : the second point
  !
  ! OUTPUT
  !    type(line) line representation
  !
  !==============================================================
  function get_line(p1, p2)
    type(point), intent(in) :: p1
    type(point), intent(in) :: p2
    type(line)              :: get_line

    ! for the two points p1 and p2, find the corresponding
    ! line representation nx*x + ny*y + c = 0
    !
    !    nx = y2 - y1
    !    ny = -(x2 - x1)
    !    c  = x2*y1 - x1*y2

    get_line%c  = p2%x*p1%y - p1%x*p2%y
    get_line%nx = p2%y - p1%y
    get_line%ny = p1%x - p2%x

  end function get_line
  

  !==============================================================
  !
  ! Function "area":
  !
  ! Calculate the area of a triangle given its vertices
  !
  !==============================================================
  !
  ! INPUT
  !    vertices : a triangle as type(triangle)
  !
  ! OUTPUT
  !     area of triangle
  !
  !
  !==============================================================
  function area(t)
    type(triangle), intent(in) :: t
    double precision           :: area

    
    ! if a,b are two triangle side vectors, then
    ! A = 0.5*|axb|
    area = 0.5D0*abs((t%p(2)%x - t%p(1)%x)*(t%p(3)%y - t%p(1)%y) - &
                     (t%p(3)%x - t%p(1)%x)*(t%p(2)%y - t%p(1)%y))
  end function area



  !==============================================================
  !
  !
  ! function "is_inside":
  !
  ! Check if a point p is inside a given triangle t
  !
  !==============================================================
  !
  ! INPUT
  !    p : the point to check
  !    t : the constraining triangle
  !
  ! OUTPUT
  !    1 : p is inside of t
  !    0 : p is on boundary of t
  !   -1 : p is outside of t
  !
  !==============================================================
  function is_inside(p, t)
    type(point), intent(in)    :: p
    type(triangle), intent(in) :: t
    integer                    :: is_inside

    ! parametrisation of triangle sides: nx*x + ny*y + c = 0
    ! where n = (nx,ny) is a normal vector of a straight line
    double precision :: nx, ny, c

    ! perpendicular in terms of normal vector
    double precision :: d, dref
    
    ! a small number compared to 1
    double precision, parameter :: EPS = 1.0D-9

    ! default is 'inside'
    is_inside = 1

    !! parametrize 1. side (p1-p2)
    ! c = xk*yi - xi*yk, nx = yk - yi, ny = -(xk - xi)
    c  = t%p(2)%x*t%p(1)%y - t%p(1)%x*t%p(2)%y
    nx = t%p(2)%y - t%p(1)%y
    ny = -t%p(2)%x + t%p(1)%x


    ! calculate the perpendicular such that p - d*n is element
    ! of side p1-p2. As we are just interested in the sign we
    ! don't divide by |n|^2
    ! Technical note: d is found by putting p - d*n into the eq
    !                 nx*x + ny*y + c = 0
    d = - c - nx*p%x - ny*p%y
    ! as we don't know in which direction n points we compare
    ! with the value of d for the third triangle vertex
    dref = - c - nx*t%p(3)%x - ny*t%p(3)%y

    ! now check if p lies on the 'right' (opposite side w.r. to
    ! p3) or on the 'left' (same side as p3) of side p1-p2
    ! NOTE: If d == 0, p could be on the boundary
    !       we check this first, as a comparison with sign(0.0)
    !       wouldn't make sense
    if (abs(d) < EPS) then
       ! p lies on the boundary or outside
       is_inside = 0
       d = dref
    end if
    if (sign(1.0D0, d) == sign(1.0D0, dref)) then
       ! they are on the same side => p is potentially inside

       !! parametrize 2. side (p2-p3)
       c  = t%p(3)%x*t%p(2)%y - t%p(2)%x*t%p(3)%y
       nx = t%p(3)%y - t%p(2)%y
       ny = -t%p(3)%x + t%p(2)%x

       ! calculate the perpendicular such that p - d*n is element
       ! of side p2-p3.
       d = - c - nx*p%x - ny*p%y
       ! dref for this side
       dref = - c - nx*t%p(1)%x - ny*t%p(1)%y

       ! check if p lies on the 'right' or on the 'left' of p2-p3
       if (abs(d) < EPS) then
          ! p lies on the boundary or outside
          is_inside = 0
          d = dref
       end if
       if (sign(1.0D0, d) == sign(1.0D0, dref)) then
          ! they are on the same side => p is potentially inside

          !! parametrize 3. side (p3-p1)
          c  = t%p(1)%x*t%p(3)%y - t%p(3)%x*t%p(1)%y
          nx = t%p(1)%y - t%p(3)%y
          ny = -t%p(1)%x + t%p(3)%x

          ! calculate the perpendicular such that p - d*n is element
          ! of side p3-p1.
          d = - c - nx*p%x - ny*p%y
          ! dref for this side
          dref = - c - nx*t%p(2)%x - ny*t%p(2)%y

          ! check if p lies on the 'right' or on the 'left' of p3-p1
          if (abs(d) < EPS) then
             ! p lies definitely on the boundary
             is_inside = 0
          else
             if (sign(1.0D0, d) == sign(1.0D0, dref)) then
                ! they are on the same side => p is definitely inside
             else
                is_inside = -1
             end if
          end if
       else
          is_inside = -1
       end if
    else
       is_inside = -1
    end if

  end function is_inside


  function sides_to_vertex(s1, s2)
    integer, intent(in) :: s1
    integer, intent(in) :: s2
    integer             :: sides_to_vertex

    select case (s1+s2)
    case (3)
        sides_to_vertex = 2
    case (4)
        sides_to_vertex = 1
    case (5)
        sides_to_vertex = 3
    end select

    end function sides_to_vertex


    subroutine sort_vertices(t)
      type(triangle), intent(inout) :: t

      type(point) :: tmp

      if (is_left_of_line(t%p(1), t%p(3), t%p(2)) == 1) then
          tmp = t%p(2)
          t%p(2) = t%p(3)
          t%p(3) = tmp
      end if
      
    end subroutine sort_vertices
      

    function is_left_of_line(p1, p2, p)
      type(point), intent(in) :: p1
      type(point), intent(in) :: p2
      type(point), intent(in) :: p
      integer                 :: is_left_of_line

      double precision :: d
      double precision :: nx, ny, c
      
      ! a small number compared to 1
      double precision, parameter :: EPS = 1.0D-9

      is_left_of_line = -1
      
      nx = -p2%y + p1%y
      ny =  p2%x - p1%x
      c  = -p2%x*p1%y + p1%x*p2%y

      d = c + nx*p%x + ny*p%y
      
      if (d > 0.0D0) then 
         is_left_of_line = 1
      elseif (abs(d) < EPS) then
         is_left_of_line = 0
      end if

      end function is_left_of_line


      function center_of_mass(t)
        type(triangle), intent(in) :: t
        type(point)                :: center_of_mass

        ! the midpoint of side 2
        type(point) :: mp2

        mp2%x = (t%p(2)%x + t%p(3)%x)/2.0D0
        mp2%y = (t%p(2)%y + t%p(3)%y)/2.0D0

        ! the center of mass can be found at 2/3 of the
        ! median to vertex 1
        center_of_mass%x = (t%p(1)%x + 2.0D0*mp2%x)/3.0D0
        center_of_mass%y = (t%p(1)%y + 2.0D0*mp2%y)/3.0D0

      end function center_of_mass


      function moment(t, m)
        type(triangle), intent(in) :: t
        character*2, intent(in)    :: m
        double precision           :: moment

        select case (m)
        case ("xx")
        case ("xy")
        case ("yy")
        end select

        moment = 0.0D0

      end function moment

end module geometry
