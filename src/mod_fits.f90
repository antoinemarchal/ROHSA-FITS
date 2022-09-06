!!
module mod_fits
  !!
  use mod_constants
  use mod_convert
  
  implicit none
  
  private
  
  public :: read_fits2D, read_fits, readp_fits, writefits2D, writefits3D, roll_fits, unroll_fits, roll_fits2D

contains

  subroutine writefits2D(filename,fitmap,n1,n2)
    
    character(*),intent(in) :: filename
    integer, intent(in) :: n1,n2 
    real(KIND=4), dimension(n1,n2), intent(in) :: fitmap
    
    integer :: stat=0,uni=0,blocksize,group,fpixel,nelements,naxis,bitpix
    logical :: simple,extend
    integer, dimension(2) :: naxes
    
    stat=0
    call ftgiou(uni,stat)
    blocksize=1
    call ftinit(uni,filename,blocksize,stat)
    simple=.true.
    bitpix=-32
    naxis=2
    naxes(1)=n1
    naxes(2)=n2
    extend=.true.
    !
    call ftphpr(uni,simple,bitpix,naxis,naxes,0,1,extend,stat)
    group=1
    fpixel=1
    nelements=naxes(1)*naxes(2)
    call ftppre(uni,group,fpixel,nelements,fitmap,stat)
    call ftclos(uni, stat)
    call ftfiou(uni, stat)
    
  end subroutine writefits2D
  
  
  subroutine writefits3D(filename,fitmap,n1,n2,n3)
    
    character(*),intent(in) :: filename
    integer, intent(in) :: n1,n2,n3
    real(KIND=4), dimension(n1,n2,n3), intent(in) :: fitmap
    
    integer :: stat=0,uni=0,blocksize,group,fpixel,nelements,naxis,bitpix
    logical :: simple,extend
    integer, dimension(3) :: naxes
    
    call ftgiou(uni,stat)
    blocksize=1
    call ftinit(uni,filename,blocksize,stat)
    simple=.true.
    bitpix=-32
    naxis=3
    naxes(1)=n1
    naxes(2)=n2
    naxes(3)=n3
    extend=.true.
    !
    call ftphpr(uni,simple,bitpix,naxis,naxes,0,1,extend,stat)
    group=1
    fpixel=1
    nelements=naxes(1)*naxes(2)*naxes(3)
    call ftppre(uni,group,fpixel,nelements,fitmap,stat)
    call ftclos(uni, stat)
    call ftfiou(uni, stat)

  end subroutine writefits3D


  subroutine readp_fits(filename, array, fpix, lpix)
    implicit none
    
    character(len=512), intent(in) :: filename
    integer(kind=4), intent(in), dimension(:), allocatable :: fpix,lpix
    real(kind=4), intent(inout), dimension(:,:,:), allocatable :: array

    integer :: stat=0,uni=0,blocksize,naxes(3)
    character(len=80) comment
    logical                                     :: undef
    character errtext*30
    
    integer, dimension(:), allocatable :: incs

    allocate(incs(3))
    incs=1

    call ftgiou(uni,stat)
    if (stat .gt. 0) then
       call ftgerr(stat,errtext)
       print *,'FITSIO Error Status =',stat,': ',errtext
       stop
    end if
    call ftdkopn(uni,filename,0,blocksize,stat)      
    call ftgkyj(uni,'NAXIS1',naxes(1),comment,stat)
    call ftgkyj(uni,'NAXIS2',naxes(2),comment,stat)
    call ftgkyj(uni,'NAXIS3',naxes(3),comment,stat)
    ! write(*,*) 'Data is',naxes(1),'by',naxes(2),'by',naxes(3)
    call ftclos(uni,stat)
      
    print*,'read the sub matrix :'
    print*,'nx start',fpix(1),'ny start',fpix(2),'nz start',fpix(3)
    print*,'nx end  ',lpix(1),'ny end  ',lpix(2),'nz end  ',lpix(3)
    
    !Read array
    call ftdkopn(uni,filename,0,blocksize,stat)
    call FTGSVE(uni,1,3,naxes,fpix,lpix,incs,0, array(:,:,:),undef,stat)
    call ftclos(uni,stat)
    flush(6)
    
  end subroutine readp_fits

  subroutine read_fits(filename, array)
    implicit none
    
    character(len=512), intent(in) :: filename
    real(kind=4), intent(inout), dimension(:,:,:), allocatable :: array

    integer :: stat=0,uni=0,blocksize,naxes(3)
    character(len=80) comment
    logical                                     :: undef
    character errtext*30
    
    integer, dimension(:), allocatable :: incs

    allocate(incs(3))
    incs=1

    call ftgiou(uni,stat)
    if (stat .gt. 0) then
       call ftgerr(stat,errtext)
       print *,'FITSIO Error Status =',stat,': ',errtext
       stop
    end if
    call ftdkopn(uni,filename,0,blocksize,stat)      
    call ftgkyj(uni,'NAXIS1',naxes(1),comment,stat)
    call ftgkyj(uni,'NAXIS2',naxes(2),comment,stat)
    call ftgkyj(uni,'NAXIS3',naxes(3),comment,stat)

    allocate(array(naxes(1),naxes(2),naxes(3)))

    !Read array
    call ftgpve(uni,1,1,naxes(1)*naxes(2)*naxes(3),0,array,undef,stat)
    call ftclos(uni,stat)
    
  end subroutine read_fits


  subroutine read_fits2D(filename, array)
    implicit none
    
    character(len=512), intent(in) :: filename
    real(kind=4), intent(inout), dimension(:,:), allocatable :: array

    integer :: stat=0,uni=0,blocksize,naxes(2)
    character(len=80) comment
    logical                                     :: undef
    character errtext*30
    
    integer, dimension(:), allocatable :: incs

    allocate(incs(3))
    incs=1

    call ftgiou(uni,stat)
    if (stat .gt. 0) then
       call ftgerr(stat,errtext)
       print *,'FITSIO Error Status =',stat,': ',errtext
       stop
    end if
    call ftdkopn(uni,filename,0,blocksize,stat)      
    call ftgkyj(uni,'NAXIS1',naxes(1),comment,stat)
    call ftgkyj(uni,'NAXIS2',naxes(2),comment,stat)

    allocate(array(naxes(1),naxes(2)))

    !Read array
    call ftgpve(uni,1,1,naxes(1)*naxes(2),0,array,undef,stat)
    call ftclos(uni,stat)
    
  end subroutine read_fits2D


  subroutine roll_fits(array, rolled)
    implicit none

    real(kind=4), intent(in), dimension(:,:,:), allocatable :: array
    real(xp), intent(inout), dimension(:,:,:), allocatable :: rolled
    integer, dimension(3) :: dim_array
    integer :: i,j,k

    dim_array = shape(array)
    allocate(rolled(dim_array(3),dim_array(2),dim_array(1)))

    !Roll 3d array and make it real xp
    do k=1,dim_array(3)
       do j=1,dim_array(2)
          do i=1,dim_array(1)
             rolled(k,j,i) = real(array(i,j,k),xp)
          end do
       end do
    end do
    
  end subroutine roll_fits


  subroutine unroll_fits(array, unrolled)
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: array
    real(kind=4), intent(inout), dimension(:,:,:), allocatable :: unrolled
    integer, dimension(3) :: dim_array
    integer :: i,j,k

    dim_array = shape(array)

    !Roll 3d array and make it real xp
    do k=1,dim_array(1)
       do j=1,dim_array(2)
          do i=1,dim_array(3)
             unrolled(i,j,k) = real(array(k,j,i),kind(4))
          end do
       end do
    end do
    
  end subroutine unroll_fits


  subroutine roll_fits2D(array, rolled)
    implicit none

    real(kind=4), intent(in), dimension(:,:), allocatable :: array
    real(xp), intent(inout), dimension(:,:), allocatable :: rolled
    integer, dimension(2) :: dim_array
    integer :: i,j

    dim_array = shape(array)
    allocate(rolled(dim_array(2),dim_array(1)))

    !Roll 2d array and make it real xp
    do j=1,dim_array(2)
       do i=1,dim_array(1)
          rolled(j,i) = real(array(i,j),xp)
       end do
    end do
    
  end subroutine roll_fits2D

end module mod_fits
