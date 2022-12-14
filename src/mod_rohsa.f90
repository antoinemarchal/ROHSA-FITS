!! This module contains ROHSA subrtoutine
module mod_rohsa
  !! This module contains ROHSA subrtoutine
  use mod_constants
  use mod_convert
  use mod_array
  use mod_functions
  use mod_start
  use mod_optimize
  use mod_inout
  use mod_fits

  implicit none

  private
  
  public :: main_rohsa

contains

  subroutine main_rohsa(data, std_data, grid_params, fileout, timeout, n_gauss, lambda_amp, lambda_mu, lambda_sig, &
       lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, amp_fact_init, sig_init, lb_sig_init, &
       ub_sig_init, lb_sig, ub_sig, maxiter_init, maxiter, m, noise, regul, descent, lstd, ustd, init_option, &
       iprint, iprint_init, save_grid, lym, init_grid, fileinit, data_init, params_init, init_spec, norm_var)
    
    implicit none
    
    logical, intent(in) :: noise           !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
    logical, intent(in) :: regul           !! if true --> activate regulation
    logical, intent(in) :: descent         !! if true --> activate hierarchical descent to initiate the optimization
    logical, intent(in) :: save_grid       !! save grid of fitted parameters at each step of the multiresolution process
    logical, intent(in) :: lym             !! if true --> activate 2-Gaussian decomposition for Lyman alpha nebula emission
    logical, intent(in) :: init_grid       !! if true --> use fileinit to give the initialization of the last grid
    logical, intent(in) :: init_spec       !! if true --> use params mean spectrum with input
    logical, intent(in) :: norm_var        !! if true --> normalize the var sig energy term                                  

    integer, intent(in) :: m               !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: lstd            !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
    integer, intent(in) :: ustd            !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
    integer, intent(in) :: iprint          !! print option 
    integer, intent(in) :: iprint_init     !! print option init
    integer, intent(in) :: maxiter         !! max iteration for L-BFGS-B alogorithm
    integer, intent(in) :: maxiter_init    !! max iteration for L-BFGS-B alogorithm (init mean spectrum)

    real(xp), intent(in) :: lambda_amp     !! lambda for amplitude parameter
    real(xp), intent(in) :: lambda_mu      !! lamnda for mean position parameter
    real(xp), intent(in) :: lambda_sig     !! lambda for dispersion parameter

    real(xp), intent(in) :: lambda_var_amp !! lambda for amp dispersion parameter
    real(xp), intent(in) :: lambda_var_mu  !! lambda for mean position dispersion parameter
    real(xp), intent(in) :: lambda_var_sig !! lambda for variance dispersion parameter

    real(xp), intent(in) :: lambda_lym_sig !! lambda for difference dispersion parameter (2-gaussaian)

    real(xp), intent(in) :: amp_fact_init  !! times max amplitude of additional Gaussian
    real(xp), intent(in) :: sig_init       !! dispersion of additional Gaussian
    real(xp), intent(in) :: ub_sig_init    !! upper bound sigma init
    real(xp), intent(in) :: lb_sig_init    !! lower bound sigma init
    real(xp), intent(in) :: lb_sig         !! lower bound sigma
    real(xp), intent(in) :: ub_sig         !! upper bound sigma

    character(len=8), intent(in)   :: init_option !!Init ROHSA with the mean or the std spectrum    
    character(len=512), intent(in) :: fileout   !! name of the output result
    character(len=512), intent(in) :: timeout   !! name of the output result
    character(len=512), intent(in) :: fileinit  !!

    integer :: n_gauss      !! number of gaussian to fit
    integer :: nside        !! size of the reshaped data \(2^{nside}\)
    integer :: n            !! loop index
    integer :: power        !! loop index

    real(xp), intent(in), dimension(:,:,:), allocatable :: data        !! initial fits data
    real(xp), intent(in), dimension(:,:,:), allocatable :: data_init   !! initial fits data grid init full cube
    real(xp), intent(in), dimension(:,:), allocatable   :: std_data    !! standard deviation map fo the cube is given by the user 
    real(xp), intent(in), dimension(:), allocatable     :: params_init !!

    real(xp), intent(inout), dimension(:,:,:), allocatable :: grid_params !! parameters to optimize at final step (dim of initial cube)

    real(xp), dimension(:,:,:), allocatable :: cube            !! reshape data with nside --> cube
    real(xp), dimension(:,:,:), allocatable :: cube_mean       !! mean cube over spatial axis
    real(xp), dimension(:,:,:), allocatable :: fit_params      !! parameters to optimize with cube mean at each iteration
    real(xp), dimension(:,:), allocatable :: std_cube          !! standard deviation map fo the cube computed by ROHSA with lb and ub
    real(xp), dimension(:,:), allocatable :: std_map           !! standard deviation map fo the cube computed by ROHSA with lb and ub
    real(xp), dimension(:), allocatable :: b_params            !! unknow average sigma
    real(xp), dimension(:), allocatable :: mean_spect          !! mean spectrum of the observation
    real(xp), dimension(:), allocatable :: guess_spect         !! params obtain fi the optimization of the std spectrum of the observation

    real(xp) :: c_lym=1._xp !! minimized the variance of the ratio between dispersion 1 and dispersion of a 2-Gaussian model for Lym alpha nebula

    integer, dimension(3) :: dim_data !! dimension of original data
    integer, dimension(3) :: dim_cube !! dimension of reshape cube
    
    real(xp), dimension(:,:), allocatable :: kernel !! convolution kernel 
    real(xp), dimension(:,:), allocatable :: map

    character(len=512) :: fileout_nside, fileout_nside_m !! name of the output at level nside   
    real(kind=4), allocatable, dimension(:,:,:) :: grid_fits

    real :: lctime, uctime, start
    
    integer :: ios=0 !! ios integer
    integer :: i     !! loop index
    integer :: j     !! loop index
    integer :: k     !! loop index
    integer :: l     !! loop index
        
    print*, "fileout = '",trim(fileout),"'"
    print*, "timeout = '",trim(timeout),"'"
    ! print*, "fileinit = '",trim(fileinit),"'"
    
    print*, " "
    print*, "______Parameters_____"
    print*, "n_gauss = ", n_gauss

    print*, "lambda_amp = ", lambda_amp
    print*, "lambda_mu = ", lambda_mu
    print*, "lambda_sig = ", lambda_sig
    print*, "lambda_var_sig = ", lambda_var_sig
    print*, "lambda_lym_sig = ", lambda_lym_sig

    print*, "amp_fact_init = ", amp_fact_init
    print*, "sig_init = ", sig_init
    print*, "lb_sig_init = ", lb_sig_init
    print*, "ub_sig_init = ", ub_sig_init
    print*, "lb_sig = ", lb_sig
    print*, "ub_sig = ", ub_sig
    print*, "init_option = ", init_option
    print*, "maxiter_init = ", maxiter_init
    print*, "maxiter = ", maxiter
    print*, "lstd = ", lstd
    print*, "ustd = ", ustd
    print*, "noise = ", noise
    print*, "regul = ", regul
    print*, "descent = ", descent
    print*, "save_grid = ", save_grid

    print*, "lym = ", lym
    print*, "init_grid = ", init_grid
    print*, "init_spec = ", init_spec
    print*, "norm_var = ", norm_var

    print*, " "

    ! Check n_gauss = 2 for Lym akpha mode
    if (lym .eqv. .true.) then
       ! if (n_gauss .eq. 2) then
          print*, "Lym alpha mode activated"
       ! else 
       !    print*, "Lym alpha mode is based on a 2-Gaussian model / please select n_gauss = 2"
       !    stop
       ! end if
    end if

    print*, " "
    
    allocate(kernel(3, 3))
    
    kernel(1,1) = 0._xp
    kernel(1,2) = -0.25_xp
    kernel(1,3) = 0._xp
    kernel(2,1) = -0.25_xp
    kernel(2,2) = 1._xp
    kernel(2,3) = -0.25_xp
    kernel(3,1) = 0._xp
    kernel(3,2) = -0.25_xp
    kernel(3,3) = 0._xp
        
    dim_data = shape(data)
    
    write(*,*) "dim_v, dim_y, dim_x = ", dim_data
    write(*,*) ""
    write(*,*) "number of los = ", dim_data(2)*dim_data(3)
    
    nside = dim2nside(dim_data)
    
    write(*,*) "nside = ", nside

    allocate(b_params(n_gauss))

    !Use grid init
    if (init_grid .eqv. .true.) then
       print*, "Ignore logical keyword descent --> deactivated"
       !Init grid_params with data_init
       grid_params = 0._xp
       grid_params = data_init

       !Init b_params
       do i=1, n_gauss       
          allocate(map(dim_data(2), dim_data(3)))
          map = grid_params(3+(3*(i-1)),:,:)
          b_params(i) = mean_2D(map, dim_data(2), dim_data(3))
          deallocate(map)
       end do
    else
       !Use ROHSA algo       
       call dim_data2dim_cube(nside, dim_data, dim_cube)
    
       !Allocate memory for cube
       allocate(cube(dim_cube(1), dim_cube(2), dim_cube(3)))
       allocate(std_cube(dim_cube(2), dim_cube(3)))
    
       !Reshape the data (new cube of size nside)
       print*, " "
       write(*,*) "Reshape cube, new dimensions :"
       write(*,*) "dim_v, dim_y, dim_x = ", dim_cube
       print*, " "

       print*, "Compute mean and std spectrum"
       allocate(mean_spect(dim_data(1)))

       call mean_spectrum(data, mean_spect, dim_data(1), dim_data(2), dim_data(3))
       call reshape_up(data, cube, dim_data, dim_cube)

       !Allocate memory for parameters grids
       if (descent .eqv. .true.) then
          allocate(fit_params(3*n_gauss, 1, 1))
          !Init sigma = 1 to avoid Nan
          do i=1,n_gauss
             fit_params(1+(3*(i-1)),1,1) = 0._xp
             fit_params(2+(3*(i-1)),1,1) = 1._xp
             fit_params(3+(3*(i-1)),1,1) = 1._xp
          end do
       end if

       print*, "                    Start iteration"
       print*, " "

       if (descent .eqv. .true.) then
          print*, "Start hierarchical descent"

          if (save_grid .eqv. .true.) then
             !Open file time step
             open(unit=11, file=timeout, status='replace', access='append', iostat=ios)
             write(11,fmt=*) "# size grid, Time (s)"
             close(11)
             call cpu_time(start)
          end if

          !Start iteration
          do n=0,nside-1
             power = 2**n

             allocate(cube_mean(dim_cube(1), power, power))

             call mean_array(power, cube, cube_mean)

             if (n == 0) then
                if (init_spec .eqv. .true.) then
                   print*, "Use user init params"
                   fit_params(:,1,1) = params_init
                else
                   if (init_option .eq. "mean") then
                      print*, "Init mean spectrum"  
                      call init_spectrum(n_gauss, fit_params(:,1,1), dim_cube(1), mean_spect, amp_fact_init, sig_init, &
                           lb_sig_init, ub_sig_init, maxiter_init, m, iprint_init)
                      ! print*, fit_params(:,1,1)
                      ! stop
                   else 
                      print*, "init_option keyword should be 'mean' or 'std' or 'max' or 'maxnorm'"
                      stop
                   end if
                end if
                
                !Init b_params
                do i=1, n_gauss       
                   b_params(i) = fit_params(3+(3*(i-1)),1,1)
                end do
             end if
                
             if (regul .eqv. .false.) then
                call upgrade(cube_mean, fit_params, power, n_gauss, dim_cube(1), lb_sig, ub_sig, maxiter, m, iprint)
             end if

             if (regul .eqv. .true.) then
                if (n == 0) then                
                   print*,  "Update level", n
                   call upgrade(cube_mean, fit_params, power, n_gauss, dim_cube(1), lb_sig, ub_sig, maxiter, m, iprint)
                end if

                if (n > 0 .and. n < nside) then
                   allocate(std_map(power, power))

                   if (noise .eqv. .true.) then
                      call reshape_noise_up(std_data, std_cube, dim_data, dim_cube)
                      call sum_map_square(power, std_cube, std_map) 
                      std_map = sqrt(std_map) / real(2**(2*(nside-n)),xp)
                   else
                      call set_stdmap(std_map, cube_mean, lstd, ustd)
                   end if

                   ! Update parameters 
                   print*,  "Update level", n, ">", power
                   call cpu_time(lctime)
                   call update(cube_mean, fit_params, b_params, n_gauss, dim_cube(1), power, power, lambda_amp, lambda_mu, &
                        lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, lb_sig, ub_sig, maxiter, &
                        m, kernel, iprint, std_map, lym, c_lym, norm_var)
                   call cpu_time(uctime)
                   print*, "Time level = ", uctime-lctime, "seconds."
                   print*, "Total time since started = ", (uctime - start), "seconds."
                   ! print*, "Estimated total time", ((lctime - start) / 3600) + (((uctime-lctime) / 3600) * 4**(nside-n)), "hours."
                   print*, "Estimated time = ", ((uctime-lctime) / 3600) * 4**(nside-n), "hours."

                   !Write output time 
                   call cpu_time(uctime)
                   write(11,fmt=*) dim_cube(2), uctime-start
                   close(11)
                 
                   deallocate(std_map)
                end if
             end if

             deallocate(cube_mean)

             ! Propagate solution on new grid (higher resolution)
             call go_up_level(fit_params)
             write(*,*) " "
             write(*,*) "Interpolate parameters level ", n!, ">", power

             ! Save grid in file
             if (save_grid .eqv. .true.) then
                print*, "Save grid parameters"
                 ! call save_process(n, n_gauss, fit_params, power, fileout)
                !Save timestep
                if (n .ne. 0) then
                   open(unit=11, file=timeout, status='unknown', access='append', iostat=ios)
                   if (ios /= 0) stop "opening file error"
                   call cpu_time(uctime)
                   write(11,fmt=*) power, uctime-start
                   close(11)
                end if
             end if

          enddo

          print*, " "
          write(*,*) "Reshape cube, restore initial dimensions :"
          write(*,*) "dim_v, dim_y, dim_x = ", dim_data

          call reshape_down(fit_params, grid_params, (/3*n_gauss, dim_cube(2), dim_cube(3)/), &
               (/3*n_gauss, dim_data(2), dim_data(3)/))

          if (save_grid .eqv. .true.) then
             !Save previous last level in fits file
             allocate(grid_fits(dim_data(3), dim_data(2),3*n_gauss))
             fileout_nside = trim(fileout(:len_trim(fileout)-5)) // "_nside_" // trim(str(nside)) // ".fits"
             call unroll_fits(grid_params, grid_fits)
             do i=1, n_gauss
                grid_fits(:,:,2+(3*(i-1))) = grid_fits(:,:,2+(3*(i-1))) - 1._xp
             end do
             call writefits3D(fileout_nside, grid_fits, dim_data(3), dim_data(2), 3*n_gauss)
             deallocate(grid_fits)
             !Save m vector (b_params)
             fileout_nside_m = trim(fileout(2:len_trim(fileout)-5)) // "_nside_" // trim(str(nside)) // "_m_vector.dat"
             open(unit=12, file=fileout_nside_m, status='replace', access='append', iostat=ios)
             if (ios /= 0) stop "opening file error"
             do i=1, n_gauss
                write(12,fmt=*) b_params(i)
             end do
             close(12)
          end if

       else
          allocate(guess_spect(3*n_gauss))
          if (init_option .eq. "mean") then
             print*, "Use of the mean spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), mean_spect, amp_fact_init, sig_init, &
                  lb_sig_init, ub_sig_init, maxiter_init, m, iprint_init)
          else
             print*, "init_option keyword should be 'mean' or 'std' or 'max'"
             stop
          end if
          call init_grid_params(grid_params, guess_spect, dim_data(2), dim_data(3))

          deallocate(guess_spect)      
       end if

    end if
    
    !Update last level
    print*, " "
    print*, "Start updating last level."
    print*, " "
    
    allocate(std_map(dim_data(2), dim_data(3)))
    
    if (noise .eqv. .true.) then
       std_map = std_data
    else   
       call set_stdmap(std_map, data, lstd, ustd)
    end if
    
    if (regul .eqv. .true.) then
       call cpu_time(lctime)
       call update(data, grid_params, b_params, n_gauss, dim_data(1), dim_data(2), dim_data(3), lambda_amp, lambda_mu, &
            lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, lb_sig, ub_sig, maxiter, m, &
            kernel, iprint, std_map, lym, c_lym, norm_var)      
       call cpu_time(uctime)
       print*, "Time = ", uctime-lctime, "seconds."
       
       if (save_grid .eqv. .true.) then
          !Write last time and close file
          open(unit=11, file=timeout, status='unknown', access='append', iostat=ios)
          if (ios /= 0) stop "opening file error"
          call cpu_time(uctime)
          write(11,fmt=*) dim_data(2), uctime-start
          close(11)
       end if
    end if
            
  end subroutine main_rohsa
  
end module mod_rohsa
