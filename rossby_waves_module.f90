module rossby_waves_module 
    use admin_module
    implicit none
    !description    contains functions that calculate the values of attributes of (ensembles of) Rossby waves
    contains

    !!! ROSSBY WAVE ATTRIBUTES

    function dispersion(k_array)result(omega)

        !description    returns a square matrix of corresponding frequencies
        !! PARAM[IN]
        !!! k_array:    equally spaced wavenumbers [-3.2,3.2]   [real array]
        !! RETURNS
        !!! omega:      square matrix of frequencies            [real array]

        integer(dp) :: dim
        real(dp),dimension(:) :: k_array
        real(dp),dimension(:,:),allocatable :: k_mat,l_mat,omega 

        dim = size(k_array)
        ! allocate arrays
        allocate(k_mat(dim,dim))
        allocate(l_mat(dim,dim))
        allocate(omega(dim,dim))

        k_mat = vec_to_mat(k_array)
        l_mat = transpose(k_mat)
        omega = -beta * k_mat/(k_mat**2 + l_mat**2+Rd**(-2))

    end function dispersion
    
    function RW_velocity(x,y,t,phases,ages,res_k,max_age,tap,alpha)result(velocity)

        !description    returns the velocity of the Rossby wave ensemble at (x,y,t)
        !! PARAMS[IN]
        !!! x:          spatial x coordinate            [real number]
        !!! y:          spatial y coordinate            [real number]
        !!! t:          time                            [real number]
        !!! phases:     array of wave phases            [real array]
        !!! ages:       array of wave ages              [real array]
        !!! res_k:      spectral resolution             [integer]
        !!! max_age:    maximum value of the age        [real number]
        !!! tap:        true when plank-taper is used   [logical]
        !!! alpha:      plank-taper parameter           [real number]
        !! RETURNS
        !!! velocity:   velocity of ensemble            [real array]
        
        real(dp) :: x,y,alpha,max_age,t
        real(dp) :: min,max,psi_0
        integer(dp) :: res_k
        real(dp), dimension(res_k) :: spec_array,ones
        real(dp), dimension(res_k,res_k) :: phases,ages,l_mat,k_mat,sdf,pre_vel
        real(dp), dimension(2) :: velocity
        logical :: tap

        ! set up
        ones = 1.0

        !setting up spectral space domain
        min = -3.2
        max = 3.2
        spec_array = linspace(min,max,res_k)
        k_mat = vec_to_mat(spec_array)
        l_mat = transpose(k_mat)

        psi_0 = psi_0_val(alpha,res_k,tap)

        if (tap .eqv. .true.) then
            sdf = plank_taper(ages,alpha,max_age,res_k)*spectral_density(spec_array)
        else
            sdf = spectral_density(spec_array)
        end if

        pre_vel = 20.0*pi*sdf*sin(2*pi*k_mat*x+2*pi*l_mat*y-dispersion(spec_array)*t + phases)

        velocity(1) = psi_0*dot_product(matmul(ones,l_mat*pre_vel),ones)
        velocity(2) = psi_0*dot_product(matmul(ones,-k_mat*pre_vel),ones)

    end function RW_velocity

    function RW_streamfunction(x,y,t,phases,res_k,max_age,ages,alpha,tap)result(streamfunction)

        !description    returns the streamfunction of the Rossby wave ensemble at (x,y,t)
        !! PARAMS[IN]
        !!! x:                  spatial x coordinate            [real number]
        !!! y:                  spatial y coordinate            [real number]
        !!! t:                  time                            [real number]
        !!! phases:             array of wave phases            [real array]
        !!! res_k:              spectral resolution             [integer]
        !!! max_age:            maximum value of the age        [real number]
        !!! ages:               array of wave ages              [real array]
        !!! alpha:              plank-taper parameter           [real number]
        !!! tap:                true when plank-taper is used   [logical]
        !! RETURNS
        !!! streamfunction:     streamfunction of ensemble      [real number]
        
        real(dp) :: x,y,alpha,max_age,t
        real(dp) :: min,max,psi_0,streamfunction
        integer(dp) :: res_k
        real(dp), dimension(res_k) :: spec_array,ones
        real(dp), dimension(res_k,res_k) :: phases,ages,l_mat,k_mat,sdf,pre_sf
        logical :: tap

        ! set up
        ones = 1.0
        psi_0 = psi_0_val(alpha,res_k,tap)

        !setting up spectral space domain
        min = -3.2
        max = 3.2
        spec_array = linspace(min,max,res_k)
        k_mat = vec_to_mat(spec_array)
        l_mat = transpose(k_mat)

        ! choose appropriate sdf
        if (tap .eqv. .true.) then
            sdf = plank_taper(ages,alpha,max_age,res_k)*spectral_density(spec_array)
        else
            sdf = spectral_density(spec_array)
        end if

        pre_sf = 4e4*sdf*cos(2*pi*k_mat*x+2*pi*l_mat*y-dispersion(spec_array)*t + phases)
        streamfunction = psi_0*dot_product(matmul(ones,pre_sf),ones)

    end function RW_streamfunction

    !!! ROSSBY WAVE ENSEMBLE ATTRIBUTES

    function spectral_density(k_array) result(sdf)

        !description     returns square matrix of corresponding spectral densities
        !! PARAM[IN]
        !!! k_array:     equally spaced wavenumbers [-3.2,3.2]   [real array]
        !! RETURNS:
        !!! sdf:         square matrix of spectral densities     [real array]

        integer(dp) :: dim
        real(dp), dimension(:) :: k_array
        real(dp), dimension(:,:),allocatable :: sdf,k_sq

        dim = size(k_array)
        ! allocate arrays
        allocate(k_sq(dim,dim))
        allocate(sdf(dim,dim))

        k_sq = vec_to_mat(k_array**2)
        sdf = exp(-k_sq-transpose(k_sq))*(k_sq+transpose(k_sq))

    end function spectral_density

    function psi_0_val(alpha,res_k,tap) result(psi_0)

        !description    Returns the value of psi_0 so that the average speed is 1
        !! PARAMS[IN]   
        !!! alpha:      plank-taper parameter           [real number]
        !!! res_k:      spectral resolution             [integer]
        !!! tap:        true when plank-taper is used   [logical]
        !! RETURNS
        !!! psi_0       normalising constant            [real number]

        real(dp) :: alpha,psi_0
        integer(dp) :: res_k
        logical :: tap

        if (tap .eqv. .false.) then
            psi_0 = 1.0/427.27
        else
            if(res_k == 65 .and. alpha == 0.1) then
                psi_0 = 1.0/401.66
            else 
                psi_0 = 1.0
            end if
        end if   
        
    end function psi_0_val
    
    function RW_ensemble_streamfunction(x_array,t,res_k,res_x,phases,max_age,ages,alpha,tap)result(streamfunction)

        !description    returns square matrix with streamfunction of the Rossby wave ensemble at each coordinate
        !! PARAMS[IN]
        !!! x_array:            vector of x_coordinates                            [real array]
        !!! t:                  time                                               [real number]
        !!! res_k:              spectral resolution                                [integer]
        !!! res_x:              spatial resolution                                 [real array]
        !!! phases:             array of wave phases                               [real array]
        !!! max_age:            maximum value of the age                           [real number]
        !!! ages:               array of wave ages                                 [real array]
        !!! alpha:              plank-taper parameter                              [real number]
        !!! tap:                true when plank-taper is used                      [logical]
        !! RETURNS
        !!! streamfunction:     streamfunction of ensemble at each coordinate      [real array]

        real(dp) :: alpha,max_age,t
        integer(dp) :: ii,jj,res_x,res_k
        real(dp), dimension(res_x) :: x_array
        real(dp), dimension(res_x,res_x) :: x_mat,y_mat,streamfunction
        real(dp), dimension(res_k,res_k) :: phases,ages
        logical :: tap
      
        ! set up spatial array
        y_mat = vec_to_mat(x_array)
        x_mat = transpose(y_mat)

        do ii = 1,res_x
            do jj = 1,res_x
                streamfunction(ii,jj) = RW_streamfunction(x_mat(ii,jj),y_mat(ii,jj),t,phases,res_k,max_age,ages,alpha,tap)
            end do
        end do

      
    end function RW_ensemble_streamfunction

end module rossby_waves_module