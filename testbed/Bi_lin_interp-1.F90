!=======================================================================
subroutine interpolate(x1, x, Dx, y1, y, Dy, z11, z21, z12, z22, z)
!=======================================================================
! calculates the value of z=f(x,y) by linearly interpolating
! between the 4 closest data points on a uniform grid.  The subroutine
! requires a grid point (x1, y1), the grid spacing (Dx and Dy), and the 
! 4 closest data points (z11, z21, z12, and z22).

use kinds

! begin input variables
real x1  ! the x grid location of z11
real x   ! x-value at which you will interpolate z=f(x,y)
real Dx  ! grid spacing in the x direction
real y1  ! the y grid location of z11
!real y   ! y-value at which you will interpolate z=f(x,y)
real(kind=dbl_kind) :: y   ! y-value at which you will interpolate z=f(x,y)
real Dy  ! grid spacing in the y direction
real z11 ! f(x1, y1)
real z21 ! f(x1+Dx, y1)
real z12 ! f(x1, y1+Dy)
real z22 ! f(x1+Dx, y1+Dy)

! begin output variables
real z   ! f(x,y), the desired interpolated value

! begin internal variables
real zp  ! z'=first interpolated value at (x, y1)
real zpp ! z''=second interpolated value at (x, Y1+Dy)


   ! interpolate between z11 and z21 to calculate z' (zp) at (x, y1)
   zp=z11+(x-x1)*(z21-z11)/Dx

   ! interpolate between z12 and z22 to calculate z'' (zpp) at (x, Y1+Dy)
   zpp=z12+(x-x1)*(z22-z12)/Dx

   ! interpolate between zp and zpp to calculate z at (x,y)
   z=zp+(y-y1)*(zpp-zp)/Dy

   return
   
end subroutine interpolate
