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
PROGRAM test
  USE simplex2D
  USE geometry
  
  IMPLICIT none

  ! test 2D BZ integration on a square domain

  ! the square:
  type(point) :: p1, p2, p3, p4

  ! the triangulation
  type(triangle), dimension(9) :: t

  ! a midpoint
  type(point) :: mp

  ! the number of energy values
  integer, parameter :: N = 4

  ! the number of bands
  integer, parameter :: NB = 1

  ! the energies we are interested in
  double precision, dimension(N) :: erg

  ! the energy dispersion values on a triangle
  double precision, dimension(6,NB) :: edisp

  ! the integral
  double precision, dimension(N) :: res, tmpres
  
  integer :: i

  interface
    function dispersion(p, m)
      use simplex2D
      use geometry
      implicit none
      type(point), intent(in) :: p
      integer, intent(in)     :: m
      double precision, dimension(1) :: dispersion
    end function dispersion
  end interface


  print *,''
  print *,'Test of 2D integration for LINE'
  print *,'==============================='
  print *,''

  ! define the square
  p1 = point(-1.0, -1.0)
  p2 = point(1.0, -1.0)
  p3 = point(1.0, 1.0)
  p4 = point(-1.0, 1.0)

  ! the energy values
  erg = (/ -0.5, 0.0, 0.7, -1.2 /)


  print *,'1. situation: e = x + y, -1 < x,y < 1'
  print *,'-----------------------------------------'
  print *,'Triangulation (A)'
  ! triangulate
  ! 1. triangle
  t(1)%p1 = point(-1.0,-1.0)
  t(1)%p2 = point(0.0,-1.0)
  t(1)%p3 = point(0.2,0.3)
  ! 2. triangle
  t(2)%p1 = point(0.0,-1.0)
  t(2)%p2 = point(1.0,-1.0)
  t(2)%p3 = point(0.2,0.3)
  ! 3. triangle
  t(3)%p1 = point(1.0,-1.0)
  t(3)%p2 = point(1.0,0.0)
  t(3)%p3 = point(0.2,0.3)
  ! 4. triangle
  t(4)%p1 = point(1.0,0.0)
  t(4)%p2 = point(1.0,1.0)
  t(4)%p3 = point(0.2,0.3)
  ! 5. triangle
  t(5)%p1 = point(1.0,1.0)
  t(5)%p2 = point(0.0,1.0)
  t(5)%p3 = point(0.2,0.3)
  ! 6. triangle
  t(6)%p1 = point(0.0,1.0)
  t(6)%p2 = point(-1.0,1.0)
  t(6)%p3 = point(0.2,0.3)
  ! 7. triangle
  t(7)%p1 = point(-1.0,1.0)
  t(7)%p2 = point(-1.0,0.0)
  t(7)%p3 = point(0.2,0.3)
  ! 8. triangle
  t(8)%p1 = point(-1.0,0.0)
  t(8)%p2 = point(-1.0,-1.0)
  t(8)%p3 = point(0.2,0.3)

  res = 0.0D0
  ! this would be very bad in a real case
  do i=1,8
     edisp(1,:) = dispersion(t(i)%p1, 0)
     edisp(2,:) = dispersion(t(i)%p2, 0)
     edisp(3,:) = dispersion(t(i)%p3, 0)

     mp%x = (t(i)%p1%x + t(i)%p2%x)/2.0
     mp%y = (t(i)%p1%y + t(i)%p2%y)/2.0
     edisp(4,:) = dispersion(mp, 0)
     mp%x = (t(i)%p2%x + t(i)%p3%x)/2.0
     mp%y = (t(i)%p2%y + t(i)%p3%y)/2.0
     edisp(5,:) = dispersion(mp, 0)
     mp%x = (t(i)%p3%x + t(i)%p1%x)/2.0
     mp%y = (t(i)%p3%y + t(i)%p1%y)/2.0
     edisp(6,:) = dispersion(mp, 0)
     
     print *,''
     print *,'TRIANGLE',i
     tmpres = simplex2Dqint(erg, t(i), edisp)
     print *,'Int on triangle',i,':',tmpres
     res = res + tmpres
  end do
  print *, ''
  print *, 'Result:', res
  print *, '-------'
  print *, ''
 
if (.TRUE.) then
  print *,''
  print *,'Triangulation (B)'
  ! triangulate
  ! 1. triangle
  t(1)%p1 = point(-1.0,-1.0)
  t(1)%p2 = point(0.0,-1.0)
  t(1)%p3 = point(-1.0,0.0)
! 9. triangle
  t(9)%p1 = point(-0.25,-0.75)
  t(9)%p2 = point(0.0,0.0)
  t(9)%p3 = point(-1.0,0.0)
  ! 2. triangle
  t(2)%p1 = point(0.0,-1.0)
  t(2)%p2 = point(0.0,0.0)
  t(2)%p3 = point(-0.25,-0.75)
  ! 3. triangle
  t(3)%p1 = point(0.0,-1.0)
  t(3)%p2 = point(1.0,-1.0)
  t(3)%p3 = point(1.0,0.0)
  ! 4. triangle
  t(4)%p1 = point(0.0,-1.0)
  t(4)%p2 = point(1.0,0.0)
  t(4)%p3 = point(0.0,0.0)
  ! 5. triangle
  t(5)%p1 = point(1.0,0.0)
  t(5)%p2 = point(1.0,1.0)
  t(5)%p3 = point(0.0,1.0)
  ! 6. triangle
  t(6)%p1 = point(0.0,0.0)
  t(6)%p2 = point(1.0,0.0)
  t(6)%p3 = point(0.0,1.0)
  ! 7. triangle
  t(7)%p1 = point(0.0,1.0)
  t(7)%p2 = point(-1.0,1.0)
  t(7)%p3 = point(-1.0,0.0)
  ! 8. triangle
  t(8)%p1 = point(-1.0,0.0)
  t(8)%p2 = point(0.0,0.0)
  t(8)%p3 = point(0.0,1.0)

  res = 0.0D0
  ! this would be very bad in a real case
  do i=1,9
     edisp(1,:) = dispersion(t(i)%p1, 0)
     edisp(2,:) = dispersion(t(i)%p2, 0)
     edisp(3,:) = dispersion(t(i)%p3, 0)

     mp%x = (t(i)%p1%x + t(i)%p2%x)/2.0
     mp%y = (t(i)%p1%y + t(i)%p2%y)/2.0
     edisp(4,:) = dispersion(mp, 0)
     mp%x = (t(i)%p2%x + t(i)%p3%x)/2.0
     mp%y = (t(i)%p2%y + t(i)%p3%y)/2.0
     edisp(5,:) = dispersion(mp, 0)
     mp%x = (t(i)%p3%x + t(i)%p1%x)/2.0
     mp%y = (t(i)%p3%y + t(i)%p1%y)/2.0
     edisp(6,:) = dispersion(mp, 0)

     print *,''
     print *,'TRIANGLE',i
     tmpres = simplex2Dqint(erg, t(i), edisp)
     print *,'Int on triangle',i,':',tmpres
     res = res + tmpres
  end do
  print *, ''
  print *, 'Result:', res
  print *, '-------'
  print *, ''
end if
 
 

END PROGRAM test



! defines our energy dispersion
function dispersion(p,m)
  use simplex2D
  use geometry

  implicit none

  type(point), intent(in) :: p
  integer, intent(in) :: m
  double precision, dimension(1) :: dispersion

  ! the expansion coefficients
  double precision :: q1, q2, q3, q4, q5, q6

  select case(m)
  case (0)
    q1 = 0.33
    q2 = 1.0
    q3 = 1.0
    q4 = 0.0
    q5 = 0.0
    q6 = 0.0
  case (1)
    q1 = 0.0
    q2 = 1.0
    q3 = 0.0
    q4 = 0.0
    q5 = 0.0
    q6 = 0.0
  end select

  dispersion(1) = q1 + q2*p%x + q3*p%y + &
               q4*p%x**2 + q5*p%x*p%y + q6*p%y**2
!  dispersion(2) = dispersion(1) + 0.5

end function dispersion



