module outputs
    use particle_trajectories_module
    implicit none

    contains

    subroutine export_RW_sf(res_k,res_x,phases,max_age,ages,alpha,tap,t_0,dt,max_t,name,path,s)

        !description    write Rossby wave streamfunction data to dat files at s intervals
        !! PARAMS[IN]
        !!! res_k:              spectral resolution                    [integer]
        !!! res_x:              spatial resolution                     [real array]
        !!! phases:             array of wave phases                   [real array]
        !!! max_age:            maximum value of the age               [real number]
        !!! ages:               array of wave ages                     [real array]
        !!! alpha:              plank-taper parameter                  [real number]
        !!! tap:                true when plank-taper is used          [logical]
        !!! t_0:                initial time                           [real number]    
        !!! dt:                 integration time step                  [real number]
        !!! max_t:              maximum value of time                  [real number]
        !!! name:               save file name                         [character, len=80]
        !!! path:               save destination path                  [character, len=100]
        !!! s:                  sampling time interval                 [real number]
    
        real(dp) :: t,max_age,alpha,t_0,dt,max_t,min,max,s
        integer(dp) :: res_k,res_x,num,ii,jj,index
        real(dp),dimension(res_k,res_k) :: phases,ages
        real(dp),dimension(res_x) :: x_array
        real(dp),dimension(res_x,res_x) :: sf
        real(dp),dimension(:),allocatable :: t_array
        logical :: tap
        character (len= 80) :: name
        character (len= 100):: filename
        character (len=200) :: path

        ! initialisation
        index = 0

        ! set up time interval
        num = size(range_real(t_0,max_t,dt))
        allocate(t_array(num))
        t_array = range_real(t_0,max_t,dt)

        ! set up x_array
        min = -5.0
        max = 5.0
        x_array = linspace(min,max,res_x)

        do ii=1,num
            t = t_array(ii)
            ! write at sampling frequency
            if (mod(t,s) == 0) then
                sf = RW_ensemble_streamfunction(x_array,t,res_k,res_x,phases,max_age,ages,alpha,tap)
                filename = index_filename(name,ii)
                open(unit=ii,file=file_path(path,filename),action='write',status='replace')
                do jj = 1,res_x
                    write(ii,"(*(G0,:,','))")sf(jj,:)
                end do
                close(ii)
                index = index + 1
            end if
            call update_array(ages,max_age,phases,res_k)
        end do
    end subroutine

    subroutine export_taper(alpha,max_age,res_k)

        !description    writes the plank-taper to dat file
        !! PARAMS[IN]
        !!! alpha:              plank-taper parameter                  [real number]
        !!! max_age:            maximum value of the age               [real number]
        !!! res_k:              spectral resolution                    [integer]


        integer(dp) :: res_k,ii,num
        real(dp) :: alpha,max_age,t_0,dt
        real(dp),dimension(res_k,res_k) :: taper,ages

        !!! Initialisation
        ages = 0.0
        t_0 = 0.
        dt = 0.01
        num = size(range_real(t_0,max_age,dt))

        open(1, file = 'taper_data.dat')

        !Write sdf data to CSV
        do ii = 1,num
            taper = plank_taper(ages,alpha,max_age,res_k)
            write(1,"(*(G0,:,','))")taper(10,20)
            ages = ages + dt
        end do

        close(1)

    end subroutine export_taper

    subroutine export_sdf(k_array,res_k)

        !description    writes the spectral density function to dat file
        !! PARAMS[IN]
        !!! k_array:            spectral array                         [real array]
        !!! res_k:              spectral resolution                    [integer]

        integer(dp) :: res_k,ii
        real(dp),dimension(res_k,res_k) :: sdf
        real(dp), dimension(res_k) :: k_array

        !!! Initialisation

        open(1, file = 'sdf_data.dat')

        !Write sdf data to CSV
    
        sdf = spectral_density(k_array)
        do ii = 1,res_k
            write(1,"(*(G0,:,','))")sdf(ii,:)
        end do

        close(1)

    end subroutine export_sdf

    subroutine export_vel(res_k,res_x,r_phases,d_phases,max_age,r_ages,d_ages,alpha,tap,t_0,gam,dt,max_t,name,path,s)
        
        !description    write velocity data to dat files at s intervals
        !! PARAMS[IN]
        !!! res_k:      spectral resolution                 [integer]
        !!! res_x:      spatial resolution                  [integer]
        !!! r_phases:   phases for rotational contribution  [real array]
        !!! d_phases:   phases for divergent contribution   [real array]
        !!! max_age:    maximum value of the age            [real number]
        !!! r_ages:     ages for rotational contribution    [real array]
        !!! d_ages:     ages for divergent contribution     [real array]
        !!! alpha:      plank-taper parameter               [real number]
        !!! tap:        true when plank-taper is used       [logical]
        !!! t_0:        initial time                        [real number]  
        !!! gam:        relative strength of divergent part [real number]  
        !!! dt:         integration time step                  [real number]
        !!! max_t:      maximum value of time                  [real number]
        !!! name:       save file name                         [character, len=80]
        !!! path:       save destination path                  [character, len=100]
        !!! s:          sampling time interval                 [real number]

        real(dp) :: max_age,alpha,t_0,dt,max_t,min,max,t,gam,s
        integer(dp) :: res_k,res_x,num,ii,jj,v_unit,index
        real(dp),dimension(res_k,res_k) :: r_phases,r_ages,d_phases,d_ages
        real(dp),dimension(:),allocatable :: t_array
        real(dp),dimension(res_x) :: x_array
        real(dp),dimension(2,res_x,res_x) :: vel
        real(dp),dimension(res_x,res_x) :: u,v,x_mat,y_mat
        logical :: tap
        character (len= 80) :: name_u,name_v,name
        character (len=200) :: path
        character (len= 100):: filename_u,filename_v

        name_u  = trim(name) // trim('_u')
        name_v = trim(name) // trim('_v')

        v_unit = 10000
        index = 0

        ! set up time interval
        num = size(range_real(t_0,max_t,dt))
        allocate(t_array(num))
        t_array = range_real(t_0,max_t,dt)

        ! set up x_array
        min = -5.0
        max = 5.0
        x_array = linspace(min,max,res_x)
        y_mat = vec_to_mat(x_array)
        x_mat = transpose(y_mat)

        do ii=1,num
            t = t_array(ii)

            ! write at sampling frequency
            if (mod(t,s) == 0) then
                vel = velocity_field(x_mat,y_mat,t,res_k,res_x,d_phases,r_phases,max_age,d_ages,r_ages,alpha,tap,gam)
                u = vel(1,:,:)
                v = vel(2,:,:)
                filename_u = index_filename(name_u,index)
                filename_v = index_filename(name_v,index)
                open(unit=ii,file=file_path(path,filename_u),action='write',status='replace')
                open(unit=v_unit-ii,file=file_path(path,filename_v),action='write',status='replace')
                do jj = 1,res_x
                    write(ii,"(*(G0,:,','))")u(jj,:)
                    write(v_unit-ii,"(*(G0,:,','))")v(jj,:)
                end do
                index = index + 1
                close(ii)
                close(v_unit-ii)
            end if
            call update_array(d_ages,max_age,d_phases,res_k)
            call update_array(r_ages,max_age,r_phases,res_k)
        end do

    end subroutine export_vel

    subroutine export_RK(n,res_k,res_x,r_phases,d_phases,max_age,r_ages,d_ages,alpha,tap,t_0,gam,dt,max_t,name,path,s)

        !description    write velocity data and particle position data to dat files at s intervals
        !! PARAMS[IN]
        !!! n:          sqrt number of particles            [integer]
        !!! res_k:      spectral resolution                 [integer]
        !!! res_x:      spatial resolution                  [integer]
        !!! r_phases:   phases for rotational contribution  [real array]
        !!! d_phases:   phases for divergent contribution   [real array]
        !!! max_age:    maximum value of the age            [real number]
        !!! r_ages:     ages for rotational contribution    [real array]
        !!! d_ages:     ages for divergent contribution     [real array]
        !!! alpha:      plank-taper parameter               [real number]
        !!! tap:        true when plank-taper is used       [logical]
        !!! t_0:        initial time                        [real number]  
        !!! gam:        relative strength of divergent part [real number]  
        !!! dt:         integration time step                  [real number]
        !!! max_t:      maximum value of time                  [real number]
        !!! name:       save file name                         [character, len=80]
        !!! path:       save destination path                  [character, len=100]
        !!! s:          sampling time interval                 [real number]
    
        real(dp) :: max_age,alpha,gam,dt,max_t,t_0,min,max,t,s
        integer(dp) :: res_k,res_x,y_unit,num,ii,jj,n,index,v_unit,kk
        real(dp),dimension(n) :: particle_array
        real(dp),dimension(res_x) :: x_array
        real(dp),dimension(res_k,res_k) :: r_phases,r_ages,d_phases,d_ages
        real(dp),dimension(res_x,res_x) :: x_mat,y_mat,u,v
        real(dp),dimension(n,n) :: pos_x,pos_y
        real(dp),dimension(2,n,n) :: pos
        real(dp),dimension(2,res_x,res_x) :: vel
        real(dp),dimension(:),allocatable :: t_array
        logical :: tap
        character (len= 80) :: name,name_x,name_y,name_u,name_v
        character (len=200) :: path
        character (len= 100):: filename_x,filename_y,filename_u,filename_v

        name_x  = trim(name) // trim('_x')
        name_y = trim(name) // trim('_y')
        name_u  = trim(name) // trim('_u')
        name_v = trim(name) // trim('_v')

        y_unit = 1000
        v_unit = 10000
        index = 0

        ! set up time interval
        num = size(range_real(t_0,max_t,dt))
        allocate(t_array(num))
        t_array = range_real(t_0,max_t,dt)

        ! basin edges
        min = -5.0
        max = 5.0

        ! basin grid points
        x_array = linspace(min,max,res_x) 
        y_mat = vec_to_mat(x_array)
        x_mat = transpose(y_mat)

        ! particle positions (uniform distribution)
        particle_array = linspace(min,max,n) 
        pos_y = vec_to_mat(particle_array)
        pos_x = transpose(vec_to_mat(particle_array))
        pos(2,:,:) = pos_y
        pos(1,:,:) = pos_x

        do ii = 1,num

            t = t_array(ii)
            
            ! to write or not to write
            if (mod(t,s) < 10e-5 .or. s - mod(t,s) < 10e-5) then

                ! write particles
                pos = RK_one_term(t,pos_x,pos_y,n,dt,res_k,r_phases,d_phases,tap,r_ages,d_ages,alpha,max_age,gam)
                call apply_dp(pos(1,:,:))
                call apply_dp(pos(2,:,:))
                pos_x = pos(1,:,:)
                pos_y = pos(2,:,:)

                filename_x = index_filename(name_x,index)
                filename_y = index_filename(name_y,index)
                open(unit=ii,file=file_path(path,filename_x),action='write',status='replace')
                open(unit=y_unit-ii,file=file_path(path,filename_y),action='write',status='replace')
                
                do jj = 1,n
                    write(ii,"(*(G0,:,','))")pos_x(jj,:)
                    write(y_unit-ii,"(*(G0,:,','))")pos_y(jj,:)
                end do
                close(ii)
                close(y_unit-ii)
                
                ! update velocity
                vel = velocity_field(x_mat,y_mat,t,res_k,res_x,d_phases,r_phases,max_age,d_ages,r_ages,alpha,tap,gam)
                u = vel(1,:,:)
                v = vel(2,:,:)

                ! write velocity
                filename_u = index_filename(name_u,index)
                filename_v = index_filename(name_v,index)
                open(unit= y_unit + 1+ ii,file=file_path(path,filename_u),action='write',status='replace')
                open(unit=v_unit-ii,file=file_path(path,filename_v),action='write',status='replace')
                do kk = 1,res_x
                    write(y_unit+1+ii,"(*(G0,:,','))")u(kk,:)
                    write(v_unit-ii,"(*(G0,:,','))")v(kk,:)
                end do
                index = index + 1
                close(y_unit+1+ii)
                close(v_unit-ii)

            end if

            call update_array(d_ages,max_age,d_phases,res_k)
            call update_array(r_ages,max_age,r_phases,res_k)

        end do

    end subroutine export_RK
end module outputs