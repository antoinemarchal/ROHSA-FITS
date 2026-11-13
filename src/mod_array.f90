!! This module contains tools (convolution/ravel) to manipulate 2D/3Darray 
module mod_array
  !! This module contains tools (convolution/ravel) to manipulate 2D/3Darray 
  use mod_constants
  implicit none

  private

  public :: convolution_2D_mirror, convolution_2D_mirror_adjoint, ravel_2D, ravel_3D, unravel_3D, mean, std, &
       std_2D, mean_2D, max_2D, ravel_3D_abs, unravel_3D_abs

  
contains

  subroutine convolution_2D_mirror(image, conv, dim_y, dim_x, kernel, dim_k)
    implicit none

    real(xp), intent(in), dimension(:,:), allocatable :: image
    real(xp), intent(inout), dimension(:,:), allocatable :: conv
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    integer, intent(in) :: dim_y, dim_x, dim_k

    integer :: i, j, ii, jj, m, n, mm, nn
    integer :: kCenterY, kCenterX
    real(xp), dimension(:,:), allocatable :: ext_conv, extended

    allocate(ext_conv(dim_y+4, dim_x+4))
    allocate(extended(dim_y+4, dim_x+4))

    ii = 0; jj = 0; mm = 0; nn = 0
    kCenterY = 0; kCenterX = 0
    ext_conv = 0._xp
    extended = 0._xp
    
    do j=1, dim_x
       do i=1, dim_y
          extended(2+i,2+j) = image(i,j)
       end do
    end do

    do j=1, 2
       do i=1, dim_y
          extended(2+i,j) = image(i,j);
       end do
    end do
    
    do i=1, 2
       do j=1, dim_x
          extended(i,2+j) = image(i,j);
       end do
    end do

    do j=dim_x+1, dim_x+2
       do i=1, dim_y
          extended(2+i,2+j) = image(i,j-2)
       end do
    end do

    do j=1, dim_x
       do i=dim_y+1, dim_y+2
          extended(2+i,2+j) = image(i-2,j)
       end do
    end do
    
    kCenterY = dim_k / 2 + 1
    kCenterX = kCenterY

    do j=1, dim_x+4
       do i=1, dim_y+4
          do m=1, dim_k
             mm = dim_k - m + 1
             
             do n=1, dim_k
                nn = dim_k - n + 1

                ii = i + (m - kCenterY)
                jj = j + (n - kCenterX)
                
                if( ii >= 1 .and. ii < dim_y+4 .and. jj >= 1 .and. jj < dim_x+4 ) then
                   ext_conv(i,j) = ext_conv(i,j) + extended(ii,jj) * kernel(mm,nn)
                end if
             end do
             
          end do
       end do
    end do

    do j=1, dim_x
       do i=1, dim_y
          conv(i,j) = ext_conv(2+i,2+j)
       end do
    end do    
    
  end subroutine convolution_2D_mirror


  !! out = D^T u, the exact adjoint of convolution_2D_mirror (same padding, same kernel):
  !! <D x, u> = <x, D^T u> for all x, u. The smoothing term is 0.5*lambda*|D x|^2, so its
  !! gradient is lambda * D^T D x. The mirror padding makes D non-symmetric in the two outer
  !! rows and columns, so D D x is wrong there, and on the coarse levels that is most pixels.
  subroutine convolution_2D_mirror_adjoint(u, out, dim_y, dim_x, kernel, dim_k)
    implicit none

    real(xp), intent(in), dimension(:,:), allocatable :: u
    real(xp), intent(inout), dimension(:,:), allocatable :: out
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    integer, intent(in) :: dim_y, dim_x, dim_k

    integer :: i, j, m, n, kc, p, q
    real(xp) :: acc
    real(xp), dimension(:,:), allocatable :: v

    kc = dim_k / 2 + 1
    ! v = C^T u on the padded grid, C being the convolution step of convolution_2D_mirror:
    ! output (i,j) reads extended(2+i+(m-kc), 2+j+(n-kc)) * kernel(dim_k-m+1, dim_k-n+1).
    allocate(v(dim_y+4, dim_x+4))
    do q=1, dim_x+4
       do p=1, dim_y+4
          acc = 0._xp
          do m=1, dim_k
             do n=1, dim_k
                i = p - 2 - (m-kc)
                j = q - 2 - (n-kc)
                if (i >= 1 .and. i <= dim_y .and. j >= 1 .and. j <= dim_x) then
                   acc = acc + u(i,j) * kernel(dim_k-m+1, dim_k-n+1)
                end if
             end do
          end do
          v(p,q) = acc
       end do
    end do

    ! P^T: add every padded cell back onto the image pixel it was copied from (same loops as the
    ! padding in convolution_2D_mirror; its corner cells are zero and carry no gradient).
    do j=1, dim_x
       do i=1, dim_y
          out(i,j) = v(2+i,2+j)
       end do
    end do
    do j=1, 2
       do i=1, dim_y
          out(i,j) = out(i,j) + v(2+i,j)
       end do
    end do
    do i=1, 2
       do j=1, dim_x
          out(i,j) = out(i,j) + v(i,2+j)
       end do
    end do
    do j=dim_x+1, dim_x+2
       do i=1, dim_y
          out(i,j-2) = out(i,j-2) + v(2+i,2+j)
       end do
    end do
    do j=1, dim_x
       do i=dim_y+1, dim_y+2
          out(i-2,j) = out(i-2,j) + v(2+i,2+j)
       end do
    end do
    deallocate(v)
  end subroutine convolution_2D_mirror_adjoint


  ! Return a contiguous flattened 1D array from a 2D array
  subroutine ravel_2D(map, vector, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_y, dim_x
    real(xp), intent(in), dimension(:,:), allocatable :: map
    real(xp), intent(inout), dimension(:), allocatable :: vector

    integer :: j, k, i__

    i__ = 1
    
    do k=1, dim_x
       do j=1, dim_y
             vector(i__) = map(j,k)
             i__ = i__ + 1
       end do
    end do
  end subroutine ravel_2D

  
  ! Return a contiguous flattened 1D array from a 3D array
  subroutine ravel_3D(cube, vector, dim_v, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(inout), dimension(:), allocatable :: vector

    integer :: i, j, k, i__

    i__ = 1
    
    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             vector(i__) = cube(i,j,k)
             i__ = i__ + 1
          end do
       end do
    end do
  end subroutine ravel_3D

  ! Return a 3D array from a contiguous flattened 1D array 
  subroutine unravel_3D(vector, cube, dim_v, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:), allocatable :: vector
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube

    integer :: i, j, k, i__

    i__ = 1
    
    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             cube(i,j,k) = vector(i__)
             i__ = i__ + 1
          end do
       end do
    end do
  end subroutine unravel_3D


  ! Return a contiguous flattened 1D array from a 3D array
  subroutine ravel_3D_abs(cube, cube_abs, vector, dim_v, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube, cube_abs
    real(xp), intent(inout), dimension(:), allocatable :: vector

    integer :: i, j, k, i__

    i__ = 1
    
    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             vector(i__) = cube(i,j,k)
             i__ = i__ + 1
          end do
       end do
    end do

    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             vector(i__) = cube_abs(i,j,k)
             i__ = i__ + 1
          end do
       end do
    end do

  end subroutine ravel_3D_abs


  ! Return a 3D array from a contiguous flattened 1D array 
  subroutine unravel_3D_abs(vector, cube, cube_abs, dim_v, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:), allocatable :: vector
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube, cube_abs

    integer :: i, j, k, i__

    i__ = 1
    
    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             cube(i,j,k) = vector(i__)
             i__ = i__ + 1
          end do
       end do
    end do

    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             cube_abs(i,j,k) = vector(i__)
             i__ = i__ + 1
          end do
       end do
    end do

  end subroutine unravel_3D_abs


  pure function mean(array)
    !! Compute the mean of a 1D array
    implicit none

    real(xp), intent(in), dimension(:) :: array !! 1D array
    real(xp) :: mean

    mean = sum(array)/(max(1,size(array)))
    
    return
  end function mean    


  pure function std(array)
    !! Compute the STD of a 1D array
    implicit none

    real(xp), intent(in), dimension(:) :: array !! 1D array
    integer :: i
    integer :: n
    real(xp) :: std !! standard deviation 
    real(xp) :: mean
    real(xp) :: var

    mean = 0._xp; var = 0._xp
    std = 0._xp

    n = size(array)
    mean = sum(array) / n

    do i=1, n
       var = var + (array(i) - mean)**2._xp
    end do
    
    var = var / (n - 1)
    std = sqrt(var)
    
    return
  end function std


  function std_2D(map, dim_y, dim_x)
    !! Compute the STD of a 2D map
    implicit none

    integer, intent(in) :: dim_y !! dimension along spatial axis y
    integer, intent(in) :: dim_x !! dimension along spatial axis x
    real(xp), intent(in), dimension(:,:), allocatable :: map !! 2D array
    real(xp), dimension(:), allocatable :: vector !! 1D array 
    real(xp) :: std_2D

    allocate(vector(dim_y*dim_x))

    call ravel_2D(map, vector, dim_y, dim_x)
    std_2D = std(vector)

    deallocate(vector)

  end function std_2D


  function max_2D(map, dim_y, dim_x)
    !! Compute the MAX of a 2D map
    implicit none

    integer, intent(in) :: dim_y !! dimension along spatial axis y
    integer, intent(in) :: dim_x !! dimension along spatial axis x
    real(xp), intent(in), dimension(:,:), allocatable :: map !! 2D array
    real(xp), dimension(:), allocatable :: vector !! 1D array 
    real(xp) :: max_2D

    allocate(vector(dim_y*dim_x))

    call ravel_2D(map, vector, dim_y, dim_x)
    max_2D = maxval(vector)

    deallocate(vector)

  end function max_2D


  function mean_2D(map, dim_y, dim_x)
    !! Compute the MEAN of a 2D map
    implicit none

    integer, intent(in) :: dim_y !! dimension along spatial axis y
    integer, intent(in) :: dim_x !! dimension along spatial axis x
    real(xp), intent(in), dimension(:,:), allocatable :: map !! 2D array
    real(xp), dimension(:), allocatable :: vector !! 1D array 
    real(xp) :: mean_2D

    allocate(vector(dim_y*dim_x))

    call ravel_2D(map, vector, dim_y, dim_x)
    mean_2D = mean(vector)

    deallocate(vector)

  end function mean_2D


end module mod_array
