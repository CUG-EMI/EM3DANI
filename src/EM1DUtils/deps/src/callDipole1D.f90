!
!  module InitDipole initialize parameters for dipole1D
! (c) PRH, 14 Aug, 2015
!
subroutine callDipole1D(txLoc, rxLoc, freqArray, sigma1D, depth1D, &
                                                nTx, nRx, nFreq, nZ, predEx1D, predEy1D, predEz1D)

    use dipole1d
    implicit none

    integer, intent(in) :: nTx, nRx, nFreq, nZ
    real(kind = 8), dimension(nTx, 5), intent(in) :: txLoc
    real(kind = 8), dimension(nRx, 3), intent(in) :: rxLoc
    real(kind = 8), dimension(nFreq), intent(in) :: freqArray
    real(kind = 8), dimension(nZ), intent(in) :: sigma1D, depth1D   ! conductivity structure, depth of top of layer
    complex(kind = 8), dimension(nRx,nTx,nFreq), intent(inout) :: predEx1D, predEy1D, predEz1D

    ! local variables
    integer :: i, j, iTx, iFreq
    real(8), dimension(:), allocatable :: sigsite

    ! number of transmitters, receivers, frequencies
    !nTx   = size(txLoc, 1)
    !nFreq = size(freqArray, 1)

    !
    ! Initialize the default values
    call init_defaults_Dipole1D

    ! Specify some parameters required by Dipole1D:
    HTmethod1D      = 'kk_ht_201'
    outputdomain1D  = 'spatial'
    phaseConvention = 'lead'
    lbcomp          = .false.
    sdm1D           = 1.0   ! (Am), dipole moment
    lUseSpline1D    = .false.
    linversion      = .false. ! Compute Derivatives with respect to sigma(layers)

    ! compute all electric field
    emFlag = 111

    !! allocate variables
    ! 1D model parameter
    nlay1D = size(depth1D, 1)
    allocate(sig1D(nlay1D), zlay1D(nlay1D))
    do i = 1, nlay1D
        sig1D(i)  = sigma1D(i)
        zlay1D(i) = depth1D(i)
    enddo

    !  receiver
    n1D  = size(rxLoc, 1)
    allocate(x1D(n1D), y1D(n1D), z1D(n1D), sigsite(n1D) )
    allocate ( ex1D(n1D), ey1D(n1D), jz1D(n1D), bx1D(n1D), by1D(n1D), bz1D(n1D) )
    do i = 1,n1D
        x1D(i) = rxLoc(i, 1)
        y1D(i) = rxLoc(i, 2)
        z1D(i) = rxLoc(i, 3)
    enddo

    !
    ! Now lets store the conductivity at the site location for later use.  Sites
    ! The assumption is that receivers on a boundary are actually in the top layer.
    ! So a sea-floor site would use the sea layer conductivity to convert to vertical
    ! current to electric field.
    do i = 1, n1D
        sigsite(i) = sig1D(i)
        do j = 2,nlay1D
            if (zlay1D(j) .lt. z1D(i)) then
                sigsite(i) = sig1D(j)
            endif
        enddo
    enddo

    ! do computation
    ! inner loop over frequencies
    do iFreq = 1, nFreq

        ftx1D = freqArray(iFreq)

        ! loop over each transmitter
        do iTx = 1, nTx

            ! assign Tx parameters:
            xTx1D = txLoc(iTx, 1)
            yTx1D = txLoc(iTx, 2)
            zTx1D = txLoc(iTx, 3)
            azimuthTx1D = txLoc(iTx, 4)
            dipTx1D     = txLoc(iTx, 5)

            ! compute CSEM fields
            call comp_dipole1D

            !
            ! Output the response for the current transmitter, note that  jz is converted to ez here
            ! The assumption is that receivers on a boundary are actually in the top layer.
            ! So a seafloor site would use the sea layer conductivity to convert to vertical
            ! current to electric field.
            do i = 1,n1D
                jz1D(i) = jz1D(i) / sigsite(i)
                predEx1D(i, iTx, iFreq) = ex1D(i)
                predEy1D(i, iTx, iFreq) = ey1D(i)
                predEz1D(i, iTx, iFreq) = jz1D(i)

            enddo

        enddo

    enddo

    !
    ! Deallocate arrays
    !
    deallocate ( x1D,y1D,z1D,sigsite )
    deallocate ( ex1D,ey1D,jz1D,bx1D,by1D,bz1D )
    deallocate ( sig1D, zlay1D )

    return

end subroutine callDipole1D

!
! callDipoleEB computes specific component of electromagnetic field of 1D CSEM
! compFlag indicates predefined component
!
subroutine callDipole1DEB(txLoc, rxLoc, freqArray, sigma1D, depth1D, &
                                                nTx, nRx, nFreq, nZ, compFlag, pred1D)

    use dipole1d
    implicit none

    integer, intent(in) :: nTx, nRx, nFreq, nZ, compFlag
    real(kind = 8), dimension(nTx, 5), intent(in) :: txLoc
    real(kind = 8), dimension(nRx, 3), intent(in) :: rxLoc
    real(kind = 8), dimension(nFreq), intent(in) :: freqArray
    real(kind = 8), dimension(nZ), intent(in) :: sigma1D, depth1D
    complex(kind = 8), dimension(nRx,nTx,nFreq), intent(inout) :: pred1D

    ! local variables
    integer :: i, j, iTx, iFreq
    real(8), dimension(:), allocatable :: sigsite

    ! Initialize the default values
    call init_defaults_Dipole1D

    ! Specify some parameters required by Dipole1D:
    HTmethod1D      = 'kk_ht_201'
    outputdomain1D  = 'spatial'
    phaseConvention = 'lead'
    lbcomp          = .false.
    sdm1D           = 1.0   ! (Am), dipole moment
    lUseSpline1D    = .true.
    linversion      = .false. ! Compute Derivatives with respect to sigma(layers)
    emFlag  = compFlag

    !! allocate variables
    ! 1D model parameter
    nlay1D = size(depth1D, 1)
    allocate(sig1D(nlay1D), zlay1D(nlay1D))
    do i = 1, nlay1D
        sig1D(i)  = sigma1D(i)
        zlay1D(i) = depth1D(i)
    enddo

    !  receiver
    n1D  = size(rxLoc, 1)
    allocate(x1D(n1D), y1D(n1D), z1D(n1D), sigsite(n1D) )
    !allocate ( ex1D(n1D), ey1D(n1D), jz1D(n1D), bx1D(n1D), by1D(n1D), bz1D(n1D) )
    select case (compFlag)
    case(1)
        allocate(ex1D(n1D))
    case(2)
        allocate(ey1D(n1D))
    case(3)
        allocate(jz1D(n1D))
    case(4)
        allocate(bx1D(n1D))
        lbcomp = .true.
    case(5)
        allocate(by1D(n1D))
        lbcomp = .true.
    case(6)
        allocate(bz1D(n1D))
        lbcomp = .true.
    end select

    do i = 1,n1D
        x1D(i) = rxLoc(i, 1)
        y1D(i) = rxLoc(i, 2)
        z1D(i) = rxLoc(i, 3)
    enddo

    !
    ! Now lets store the conductivity at the site location for later use.  Sites
    ! The assumption is that receivers on a boundary are actually in the top layer.
    ! So a sea-floor site would use the sea layer conductivity to convert to vertical
    ! current to electric field.
    do i = 1, n1D
        sigsite(i) = sig1D(i)
        do j = 2,nlay1D
            if (zlay1D(j) .lt. z1D(i)) then
                sigsite(i) = sig1D(j)
            endif
        enddo
    enddo

    ! do computation
    ! inner loop over frequencies
    do iFreq = 1, nFreq

        ftx1D = freqArray(iFreq)

        ! loop over each transmitter
        do iTx = 1, nTx

            ! assign Tx parameters:
            xTx1D = txLoc(iTx, 1)
            yTx1D = txLoc(iTx, 2)
            zTx1D = txLoc(iTx, 3)
            azimuthTx1D = txLoc(iTx, 4)
            dipTx1D     = txLoc(iTx, 5)

            ! compute CSEM fields
            call comp_dipole1D

            !
            ! Output the response for the current transmitter, note that  jz is converted to ez here
            ! The assumption is that receivers on a boundary are actually in the top layer.
            ! So a seafloor site would use the sea layer conductivity to convert to vertical
            ! current to electric field.
            do i = 1,n1D

                if (compFlag .eq. 1) then
                    pred1D(i, iTx, iFreq) = ex1D(i)
                elseif (compFlag .eq. 2) then
                    pred1D(i, iTx, iFreq) = ey1D(i)
                elseif (compFlag .eq. 3) then
                    jz1D(i) = jz1D(i) / sigsite(i)
                    pred1D(i, iTx, iFreq) = jz1D(i)
                elseif  (compFlag .eq. 4) then
                    pred1D(i, iTx, iFreq) = bx1D(i)
                elseif (compFlag .eq. 5) then
                    pred1D(i, iTx, iFreq) = by1D(i)
                elseif (compFlag .eq. 6) then
                    pred1D(i, iTx, iFreq) = bz1D(i)
                endif

            enddo

        enddo

    enddo

    !
    ! Deallocate arrays
    !
    deallocate ( x1D,y1D,z1D,sigsite )
    deallocate ( sig1D, zlay1D )
    select case(compFlag)
    case(1)
        deallocate ( ex1D)
    case(2)
        deallocate(ey1D)
    case(3)
        deallocate(jz1D)
    case(4)
        deallocate(bx1D)
    case(5)
        deallocate(by1D)
    case(6)
        deallocate(bz1D)
    end select

    return

end subroutine callDipole1DEB

subroutine callDipole1DFinite(txLoc, rxLoc, dipLen, nP, freqArray, sigma1D, depth1D, &
                                                nTx, nRx, nFreq, nZ, compFlag, pred1D)

    use dipole1d
    implicit none

    integer, intent(in) :: nP, nTx, nRx, nFreq, nZ, compFlag
    real(kind = 8), intent(in) :: dipLen
    real(kind = 8), dimension(nTx, 5), intent(in) :: txLoc
    real(kind = 8), dimension(nRx, 3), intent(in) :: rxLoc
    real(kind = 8), dimension(nFreq), intent(in) :: freqArray
    real(kind = 8), dimension(nZ), intent(in) :: sigma1D, depth1D
    complex(kind = 8), dimension(nRx,nTx,nFreq), intent(inout) :: pred1D

    ! local variables
    integer :: i, j, iTx, iFreq
    real(8), dimension(:), allocatable :: sigsite

    ! Initialize the default values
    call init_defaults_Dipole1D

    ! Specify some parameters required by Dipole1D:
    HTmethod1D      = 'kk_ht_201'
    outputdomain1D  = 'spatial'
    phaseConvention = 'lead'
    lbcomp          = .false.
    sdm1D           = 1.0   ! (Am), dipole moment
    lUseSpline1D    = .true.
    linversion      = .false. ! Compute Derivatives with respect to sigma(layers)

    ! initialize relevent parameters
    emFlag   = compFlag
    lenTx1D = dipLen
    numIntegPts = nP

    !! allocate variables
    ! 1D model parameter
    nlay1D = size(depth1D, 1)
    allocate(sig1D(nlay1D), zlay1D(nlay1D))
    do i = 1, nlay1D
        sig1D(i)  = sigma1D(i)
        zlay1D(i) = depth1D(i)
    enddo

    !  receiver
    n1D  = size(rxLoc, 1)
    allocate(x1D(n1D), y1D(n1D), z1D(n1D), sigsite(n1D) )
    !allocate ( ex1D(n1D), ey1D(n1D), jz1D(n1D), bx1D(n1D), by1D(n1D), bz1D(n1D) )
    select case (compFlag)
    case(1)
        allocate(ex1D(n1D))
    case(2)
        allocate(ey1D(n1D))
    case(3)
        allocate(jz1D(n1D))
    case(4)
        allocate(bx1D(n1D))
        lbcomp = .true.
    case(5)
        allocate(by1D(n1D))
        lbcomp = .true.
    case(6)
        allocate(bz1D(n1D))
        lbcomp = .true.
    end select

    do i = 1,n1D
        x1D(i) = rxLoc(i, 1)
        y1D(i) = rxLoc(i, 2)
        z1D(i) = rxLoc(i, 3)
    enddo

    !
    ! Now lets store the conductivity at the site location for later use.  Sites
    ! The assumption is that receivers on a boundary are actually in the top layer.
    ! So a sea-floor site would use the sea layer conductivity to convert to vertical
    ! current to electric field.
    do i = 1, n1D
        sigsite(i) = sig1D(i)
        do j = 2,nlay1D
            if (zlay1D(j) .lt. z1D(i)) then
                sigsite(i) = sig1D(j)
            endif
        enddo
    enddo

    ! do computation
    ! inner loop over frequencies
    do iFreq = 1, nFreq

        ftx1D = freqArray(iFreq)

        ! loop over each transmitter
        do iTx = 1, nTx

            ! assign Tx parameters:
            xTx1D = txLoc(iTx, 1)
            yTx1D = txLoc(iTx, 2)
            zTx1D = txLoc(iTx, 3)
            azimuthTx1D = txLoc(iTx, 4)
            dipTx1D     = txLoc(iTx, 5)

            ! compute CSEM fields
            call comp_dipole1D

            !
            ! Output the response for the current transmitter, note that  jz is converted to ez here
            ! The assumption is that receivers on a boundary are actually in the top layer.
            ! So a seafloor site would use the sea layer conductivity to convert to vertical
            ! current to electric field.
            do i = 1,n1D

                if (compFlag .eq. 1) then
                    pred1D(i, iTx, iFreq) = ex1D(i)
                elseif (compFlag .eq. 2) then
                    pred1D(i, iTx, iFreq) = ey1D(i)
                elseif (compFlag .eq. 3) then
                    jz1D(i) = jz1D(i) / sigsite(i)
                    pred1D(i, iTx, iFreq) = jz1D(i)
                elseif  (compFlag .eq. 4) then
                    pred1D(i, iTx, iFreq) = bx1D(i)
                elseif (compFlag .eq. 5) then
                    pred1D(i, iTx, iFreq) = by1D(i)
                elseif (compFlag .eq. 6) then
                    pred1D(i, iTx, iFreq) = bz1D(i)
                endif

            enddo

        enddo

    enddo

    !
    ! Deallocate arrays
    !
    deallocate ( x1D,y1D,z1D,sigsite )
    deallocate ( sig1D, zlay1D )
    select case(compFlag)
    case(1)
        deallocate ( ex1D)
    case(2)
        deallocate(ey1D)
    case(3)
        deallocate(jz1D)
    case(4)
        deallocate(bx1D)
    case(5)
        deallocate(by1D)
    case(6)
        deallocate(bz1D)
    end select

    return

end subroutine callDipole1DFinite






