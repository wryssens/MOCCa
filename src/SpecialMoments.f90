module SpecialMoments
!--------------------------------------------------------------------------------------
! Module to contain definitions of specially defined moments, that can be constrained.
!
! If you want to add weird constraints, you should add in this module
!
! a) A custom Iniroutine that takes care of the initialisation
! b) A custom Calculate routine that calculates your moment
! c) A custom PrintRoutine that prints info on your moment.
!
! Make sure that the SpherHarm variable is correctly initialised!
!
! Currently implemented:
!
! - Moment of the Schiff operator
!
!---------------------------------------------------------------------------------------


use Moments
use Densities
use Mesh

implicit none

    ! Contain the Schiff operator
    type(Moment), target :: Schiff

contains

    subroutine ReadSpecialMoments
        !----------------------------------------------------
        ! Subroutine that gets input on non-standard moments
        ! that might be desired.
        !
        !----------------------------------------------------

        integer      :: SchiffMoment = 0
        real(KIND=dp):: Constraint = 0.0_dp, Intensity=0.0_dp
        integer      :: ConstraintType=0

        Namelist /SpecialMoments/ SchiffMoment, Constraint, &
        &                         ConstraintType, Intensity

        read(unit=*, NML=SpecialMoments)

        if(SchiffMoment.eq.1) call IniSchiffOperator

        if(ConstraintType.ne.0) then
            allocate(Schiff%TrueConstraint(2) , Schiff%Constraint(2), Schiff%Intensity(2))
            Schiff%TrueConstraint(1) = 0
            Schiff%Constraint(1)     = 0
            ! Schiff should only constrain the proton value.

            Schiff%ConstraintType = ConstraintType

            Schiff%TrueConstraint(2) = Constraint
            Schiff%Constraint(2)     = Constraint

            Schiff%Isoswitch      = 2
        endif


    end subroutine ReadSpecialMoments


    subroutine IniSchiffOperator

        type(Moment), pointer :: current

        !------------------------------------------
        ! Note that -1 is taken as l to signify 1)
        ! that it is not a multipole moment
        ! and that 2) it is parity-odd.
        Schiff = NewMoment(-1, 0, 0)

        ! Notice that we don't initialize the spherical
        ! harmonic for the Schiff moment here.

        !----------------------------------------
        ! Find the last moment and associate the
        ! and put the Schiff moment behind it
        current => Root

        do while(associated(Current%Next))
            Current => Current%Next
        end do

        ! We now have the last moment
        Current%Next        => Schiff
        Schiff%Calculate   => CalculateSchiff
        Schiff%PrintMoment => PrintSchiff

    end subroutine IniSchiffOperator

    subroutine CalculateSchiff(SchiffMoment,SaveOld)
        !------------------------------------
        ! Calculate the Schiff moment.
        !
        !------------------------------------

        use Force, only : e2

        integer                      :: i,j,k
        type(Moment)                 :: Rms
        real(KIND=dp)                :: msr
        class(Moment), intent(inout) :: SchiffMoment
        integer                      :: SaveOld
        !----------------------------------------------
        ! This routine is somewhat weird, but this is
        ! necessary, since the Schiff moment is state
        ! dependent. This is why we recalculate the
        ! operator here every time.
        !----------------------------------------------

        Rms        = FindMoment(-2,0,.false.)
        msr        = 5.0_dp/3.0_dp * Rms%Value(2)

        if(SaveOld.ne.0) then
            do i=1,6
                SchiffMoment%OldValue(:,7-i+1) = SchiffMoment%OldValue(:,7-i)
            enddo
            SchiffMoment%OldValue(:,1)= SchiffMoment%Value
        endif

        SchiffMoment%Value   = 0.0_dp
        SchiffMoment%Squared = 0.0_dp

        do k=1,nz
            do j=1,ny
                do i=1,nx
                    Schiff%SpherHarm(i,j,k) = e2/10.0 *            &
                    &                        (( MeshX(i)**2        &
                    &                         + MeshY(j)**2        &
                    &                         + MeshZ(k)**2 )**2   &
                    &                         - msr ) * MeshZ(k)
                enddo
            enddo
        enddo

        ! Calculate the new value, but only for the protons.
        SchiffMoment%Value(2)      = sum(SchiffMoment%SpherHarm   *Density%Rho(:,:,:,2))*dv
        Schiffmoment%ValueCut(2)   = sum(SchiffMoment%SpherHarm   *Density%Rho(:,:,:,2) *Cutoff(:,:,:,2))*dv
        Schiffmoment%Squared(2)    = sum(SchiffMoment%SpherHarm**2*Density%Rho(:,:,:,2))*dv
        Schiffmoment%SquaredCut(2) = sum(SchiffMoment%SpherHarm**2*Density%Rho(:,:,:,2) *Cutoff(:,:,:,2))*dv

    end subroutine CalculateSchiff

    subroutine PrintSchiff(ToPrint)
        !-------------------------------------
        ! Printroutine for the Schiff moment.
        !
        !-------------------------------------
        class(Moment),       intent(in) :: ToPrint

        1 format (' Schiff Op.', 3(1x,f15.4))
        2 format ('Constrained',  2(1x,f15.4))
        3 format ('Parameter  ', 17x,  f15.4)

        print 1, ToPrint%Value, sum(ToPrint%Value)

        if(ToPrint%ConstraintType.ne.0) then
            print 2 , ToPrint%TrueConstraint
            print 3,  ToPrint%Intensity(2)
        endif

    end subroutine PrintSchiff


end module SpecialMoments
