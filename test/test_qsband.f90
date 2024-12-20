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
  type(point), dimension(4) :: p1, p2, p3, p4

  ! the triangulation
  type(triangle), dimension(8192) :: t

  ! a midpoint
  type(point) :: mp

  ! the number of energy values
  integer, parameter :: N = 501

  ! the energies we are interested in
  double precision, dimension(N) :: erg

  ! the energy dispersion values on a triangle
  double precision, dimension(6) :: edisp

  ! the integral
  double precision, dimension(N) :: res, tmpres, linres
  double precision :: D = 0.125
  double precision :: x, y, det
  
  integer :: i,j,id

  interface
    function dispersion(p)
      use simplex2D
      use geometry
      implicit none
      type(point), intent(in) :: p
      double precision :: dispersion
    end function dispersion
  end interface


  print *,''
  print *,'% Test of 2D integration for cubic s-band'
  print *,'% ======================================='
  print *,''

  print *, 'points=[25,81,289,1089];'
  
  ! define the square
  p1 = point(0.0, 0.0)
  p2 = point(1.0, 0.0)
  p3 = point(1.0, 1.0)
  p4 = point(0.0, 1.0)

  ! the energy values
  erg = (/ (i,i=0,500) /)
  erg = 0.002D0*erg
  
  ! the linear interpolation 0
  D = 0.25

  ! triangulate
  id=1
  y=0.0
  do i=1,4
    x=0.0
    do j=1,4
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% linear interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     ! linear approx
     call sort(edisp(1:3))
     det = 2*area(t(i))
     do j=1,N
       if (edisp(1) < erg(j) .and. erg(j) < edisp(2)) then
         linres(j) = linres(j) + (erg(j)-edisp(1))/((edisp(2)-edisp(1))* &
         (edisp(3)-edisp(1)))*det
       else if (edisp(2) < erg(j) .and. erg(j) < edisp(3)) then
         linres(j) = linres(j) + (edisp(3)-erg(j))/((edisp(3)-edisp(2))* &
         (edisp(3)-edisp(1)))*det
       end if
     end do
  end do
  print *, ''
  print *, 'Dlin0=['
  do i=1,N
    print *, 4*linres(i)
  end do
  print *, '];'
 
  ! the linear interpolation 1
  D = 0.125

  ! triangulate
  id=1
  y=0.0
  do i=1,8
    x=0.0
    do j=1,8
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% linear interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     ! linear approx
     call sort(edisp(1:3))
     det = 2*area(t(i))
     do j=1,N
       if (edisp(1) < erg(j) .and. erg(j) < edisp(2)) then
         linres(j) = linres(j) + (erg(j)-edisp(1))/((edisp(2)-edisp(1))* &
         (edisp(3)-edisp(1)))*det
       else if (edisp(2) < erg(j) .and. erg(j) < edisp(3)) then
         linres(j) = linres(j) + (edisp(3)-erg(j))/((edisp(3)-edisp(2))* &
         (edisp(3)-edisp(1)))*det
       end if
     end do
  end do
  print *, ''
  print *, 'Dlin1=['
  do i=1,N
    print *, 4*linres(i)
  end do
  print *, '];'
 
  ! the linear interpolation 2
  D = 0.0625

  ! triangulate
  id=1
  y=0.0
  do i=1,16
    x=0.0
    do j=1,16
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% linear interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     ! linear approx
     call sort(edisp(1:3))
     det = 2*area(t(i))
     do j=1,N
       if (edisp(1) < erg(j) .and. erg(j) < edisp(2)) then
         linres(j) = linres(j) + (erg(j)-edisp(1))/((edisp(2)-edisp(1))* &
         (edisp(3)-edisp(1)))*det
       else if (edisp(2) < erg(j) .and. erg(j) < edisp(3)) then
         linres(j) = linres(j) + (edisp(3)-erg(j))/((edisp(3)-edisp(2))* &
         (edisp(3)-edisp(1)))*det
       end if
     end do
  end do
  print *, ''
  print *, 'Dlin2=['
  do i=1,N
    print *, 4*linres(i)
  end do
  print *, '];'
  
  ! the linear interpolation 3
  D = 0.03125

  ! triangulate
  id=1
  y=0.0
  do i=1,32
    x=0.0
    do j=1,32
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% linear interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     ! linear approx
     call sort(edisp(1:3))
     det = 2*area(t(i))
     do j=1,N
       if (edisp(1) < erg(j) .and. erg(j) < edisp(2)) then
         linres(j) = linres(j) + (erg(j)-edisp(1))/((edisp(2)-edisp(1))* &
         (edisp(3)-edisp(1)))*det
       else if (edisp(2) < erg(j) .and. erg(j) < edisp(3)) then
         linres(j) = linres(j) + (edisp(3)-erg(j))/((edisp(3)-edisp(2))* &
         (edisp(3)-edisp(1)))*det
       end if
     end do
  end do
  print *, ''
  print *, 'Dlin3=['
  do i=1,N
    print *, 4*linres(i)
  end do
  print *, '];'
 
  ! quadratic interpolation 0
  D = 0.5

  ! triangulate
  id=1
  y=0.0
  do i=1,2
    x=0.0
    do j=1,2
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% quadratic interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     mp%x = (t(i)%p(1)%x + t(i)%p(2)%x)/2.0
     mp%y = (t(i)%p(1)%y + t(i)%p(2)%y)/2.0
     edisp(4) = dispersion(mp)
     mp%x = (t(i)%p(2)%x + t(i)%p(3)%x)/2.0
     mp%y = (t(i)%p(2)%y + t(i)%p(3)%y)/2.0
     edisp(5) = dispersion(mp)
     mp%x = (t(i)%p(3)%x + t(i)%p(1)%x)/2.0
     mp%y = (t(i)%p(3)%y + t(i)%p(1)%y)/2.0
     edisp(6) = dispersion(mp)
     
     tmpres = simplex2Dqint(erg, t(i), edisp)
     res = res + tmpres
  end do
  print *, ''
  print *, 'Dquad0=['
  do i=1,N
    print *, 4.0*res(i)
  end do
  print *, '];'
 

  ! quadratic interpolation 1
  D = 0.25

  ! triangulate
  id=1
  y=0.0
  do i=1,4
    x=0.0
    do j=1,4
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% quadratic interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     mp%x = (t(i)%p(1)%x + t(i)%p(2)%x)/2.0
     mp%y = (t(i)%p(1)%y + t(i)%p(2)%y)/2.0
     edisp(4) = dispersion(mp)
     mp%x = (t(i)%p(2)%x + t(i)%p(3)%x)/2.0
     mp%y = (t(i)%p(2)%y + t(i)%p(3)%y)/2.0
     edisp(5) = dispersion(mp)
     mp%x = (t(i)%p(3)%x + t(i)%p(1)%x)/2.0
     mp%y = (t(i)%p(3)%y + t(i)%p(1)%y)/2.0
     edisp(6) = dispersion(mp)
     
     tmpres = simplex2Dqint(erg, t(i), edisp)
     res = res + tmpres
  end do
  print *, ''
  print *, 'Dquad1=['
  do i=1,N
    print *, 4.0*res(i)
  end do
  print *, '];'
 
  ! quadratic interpolation 2
  D = 0.125

  ! triangulate
  id=1
  y=0.0
  do i=1,8
    x=0.0
    do j=1,8
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% quadratic interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     mp%x = (t(i)%p(1)%x + t(i)%p(2)%x)/2.0
     mp%y = (t(i)%p(1)%y + t(i)%p(2)%y)/2.0
     edisp(4) = dispersion(mp)
     mp%x = (t(i)%p(2)%x + t(i)%p(3)%x)/2.0
     mp%y = (t(i)%p(2)%y + t(i)%p(3)%y)/2.0
     edisp(5) = dispersion(mp)
     mp%x = (t(i)%p(3)%x + t(i)%p(1)%x)/2.0
     mp%y = (t(i)%p(3)%y + t(i)%p(1)%y)/2.0
     edisp(6) = dispersion(mp)
     
     tmpres = simplex2Dqint(erg, t(i), edisp)
     res = res + tmpres
  end do
  print *, ''
  print *, 'Dquad2=['
  do i=1,N
    print *, 4.0*res(i)
  end do
  print *, '];'
  
  ! quadratic interpolation 3
  D = 0.0625

  ! triangulate
  id=1
  y=0.0
  do i=1,16
    x=0.0
    do j=1,16
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% quadratic interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     mp%x = (t(i)%p(1)%x + t(i)%p(2)%x)/2.0
     mp%y = (t(i)%p(1)%y + t(i)%p(2)%y)/2.0
     edisp(4) = dispersion(mp)
     mp%x = (t(i)%p(2)%x + t(i)%p(3)%x)/2.0
     mp%y = (t(i)%p(2)%y + t(i)%p(3)%y)/2.0
     edisp(5) = dispersion(mp)
     mp%x = (t(i)%p(3)%x + t(i)%p(1)%x)/2.0
     mp%y = (t(i)%p(3)%y + t(i)%p(1)%y)/2.0
     edisp(6) = dispersion(mp)
     
     tmpres = simplex2Dqint(erg, t(i), edisp)
     res = res + tmpres
  end do
  print *, ''
  print *, 'Dquad3=['
  do i=1,N
    print *, 4.0*res(i)
  end do
  print *, '];'
  
  ! quadratic interpolation 4
  D = 0.03125D0

  ! triangulate
  id=1
  y=0.0
  do i=1,32
    x=0.0
    do j=1,32
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% quadratic interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     mp%x = (t(i)%p(1)%x + t(i)%p(2)%x)/2.0
     mp%y = (t(i)%p(1)%y + t(i)%p(2)%y)/2.0
     edisp(4) = dispersion(mp)
     mp%x = (t(i)%p(2)%x + t(i)%p(3)%x)/2.0
     mp%y = (t(i)%p(2)%y + t(i)%p(3)%y)/2.0
     edisp(5) = dispersion(mp)
     mp%x = (t(i)%p(3)%x + t(i)%p(1)%x)/2.0
     mp%y = (t(i)%p(3)%y + t(i)%p(1)%y)/2.0
     edisp(6) = dispersion(mp)
     
     tmpres = simplex2Dqint(erg, t(i), edisp)
     res = res + tmpres
  end do
  print *, ''
  print *, 'Dquad4=['
  do i=1,N
    print *, 4.0*res(i)
  end do
  print *, '];'
 

  ! quadratic interpolation (true)
  D = 0.015625D0

  ! triangulate
  id=1
  y=0.0
  do i=1,64
    x=0.0
    do j=1,64
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y)
      t(id)%p(3) = point(x+D,y+D)
      id=id+1
      
      t(id)%p(1) = point(x,y)
      t(id)%p(2) = point(x+D,y+D)
      t(id)%p(3) = point(x,y+D)
      id=id+1
      x=x+D
    end do
    y=y+D
  end do
  res = 0.0D0
  linres = 0.0D0
  print *, '% quadratic interpolation'
  print *, '% ',id-1,'triangles'
  do i=1,id-1
     edisp(1) = dispersion(t(i)%p(1))
     edisp(2) = dispersion(t(i)%p(2))
     edisp(3) = dispersion(t(i)%p(3))

     mp%x = (t(i)%p(1)%x + t(i)%p(2)%x)/2.0
     mp%y = (t(i)%p(1)%y + t(i)%p(2)%y)/2.0
     edisp(4) = dispersion(mp)
     mp%x = (t(i)%p(2)%x + t(i)%p(3)%x)/2.0
     mp%y = (t(i)%p(2)%y + t(i)%p(3)%y)/2.0
     edisp(5) = dispersion(mp)
     mp%x = (t(i)%p(3)%x + t(i)%p(1)%x)/2.0
     mp%y = (t(i)%p(3)%y + t(i)%p(1)%y)/2.0
     edisp(6) = dispersion(mp)
     
     tmpres = simplex2Dqint(erg, t(i), edisp)
     res = res + tmpres
  end do
  print *, ''
  print *, 'Dtrue=['
  do i=1,N
    print *, 4.0*res(i)
  end do
  print *, '];'
 

  print *, ''
  print *, 'E=['
  do i=1,N
    print *, erg(i)
  end do
  print *, '];'

 

END PROGRAM test



! defines our energy dispersion
function dispersion(p)
  use simplex2D
  use geometry

  implicit none

  type(point), intent(in) :: p
  double precision :: dispersion

  dispersion = -0.5D0*(cos(PI/2*p%x) + cos(PI/2*p%y))+1.0D0
  !dispersion = 0.5D0*(p%x**4 + p%y**4)
  !dispersion = sqrt(sqrt(p%x**2 + p%y**2))

end function dispersion



