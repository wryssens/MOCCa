module PotentialMixing

  use Energy
  use Meanfields
  use Coulomb, only: SolveCoulomb
  use Derivatives
  use Damping
  
  implicit none
  
  real(KIND=dp),allocatable  :: Upotold(:,:,:,:), Spotold(:,:,:,:,:)
  real(KIND=dp),allocatable  :: Wpotold(:,:,:,:,:,:)
  real(KIND=dp), allocatable :: du(:,:,:,:), ds(:,:,:,:,:), dw(:,:,:,:,:,:)
    
contains  

  subroutine ConstructPotentials
    !---------------------------------------------------------------------------
    ! Calculate the new mean-field potentials from the current densities. 
    ! 
    ! If mixing is required, mix the updates.
    !
    !---------------------------------------------------------------------------
  
    real(KIND=dp) :: alpha(nx,ny,nz,2), beta(nx,ny,nz,2) , update(nx,ny,nz,2)
    real(KIND=dp), save :: olddiff, aeff 
  
    integer :: i , P, S, TS
    
    !---------------------------------------------------------------------------
    !Allocate the potentials if this has not already happened
    if(.not.allocated(UPot)) then
      !-------------------------------------------------------------------------
      ! Time-even potentials
      allocate(BPot(nx,ny,nz,2), NablaBPot(nx,ny,nz,3,2), UPot(nx,ny,nz,2))
      allocate(UpotOld(nx,ny,nz,2)); allocate(du(nx,ny,nz,2))
      allocate(Wpot(nx,ny,nz,3,3,2))
      BPot = 0.0_dp ; NablaBPot=0.0_dp 
      Upot = 0.0_dp ; Upotold = 0.0_dp ; du = 0.0_dp
      Wpot = 0.0_dp 
      if ( PreConW .gt. 0.0_dp ) then 
        allocate(WpotOld(nx,ny,nz,3,3,2)) ; allocate(dw(nx,ny,nz,3,3,2)) 
        Wpotold = 0.0_dp ; dw = 0.0_dp
      endif
      if(.not. TRC) then
        !-----------------------------------------------------------------------
        ! Time-odd potentials
        allocate(SPot(nx,ny,nz,3,2),APot(nx,ny,nz,3,2))
        SPot     =0.0_dp ; Apot = 0.0_dp
        allocate(ds(nx,ny,nz,3,2)) ; ds = 0
        allocate(SpotOld(nx,ny,nz,3,2)) ; SpotOld = 0.0
      endif
     
      !-------------------------------------------------------------------------
      ! N2LO potentials, time-even 
      if(t1n2.ne.0.0_dp .or. t2n2.ne.0.0_dp) then
        allocate(Bmunu(nx,ny,nz,3,3,2))     ; allocate(DmuBmunu(nx,ny,nz,3,3,2))
        allocate(DN2LO(nx,ny,nz,2))         ; allocate(Xpot(nx,ny,nz,3,3,2))
        allocate(DXpot(nx,ny,nz,3,3,3,2))   ; allocate(ImTfield(nx,ny,nz,3,3,2))
        
        Bmunu = 0.0_dp    ; DmuBmunu  = 0.0_dp
        DN2LO = 0.0_dp    ; Xpot    = 0.0_dp
        DXpot = 0.0_dp    ; ImTField= 0.0_dp
        
        if(.not. TRC) then
            !-------------------------------------------------------------------
            ! Time-odd N2LO potentials
            allocate(SN2LOField(nx,ny,nz,3,2))   ; SN2LOField = 0.0_dp
            allocate(ReTfield(nx,ny,nz,3,3,3,2)) ; ReTfield   = 0.0_dp 
            allocate(ReDTfield(nx,ny,nz,3,3,2))  ; ReDTfield  = 0.0_dp
            allocate(PiField(nx,ny,nz,3,2))      ; PiField    = 0.0_dp
        endif
      endif
      
      !-------------------------------------------------------------------------
      ! Tensor potentials, all time-odd
      if(B15.ne.0.0_dp .or. B16.ne.0.0_dp .or. &
      &  B17.ne.0.0_dp .or. B18.ne.0.0_dp) then
        allocate(Cpot(nx,ny,nz,3,2),DPot(nx,ny,nz,3,2))
        allocate(DerCPot(nx,ny,nz,3,3,2),DerDPot(nx,ny,nz,3,3,2))
        allocate(DivDpot(nx,ny,nz,2))
        
        CPot     =0.0_dp ; DPot     =0.0_dp
        DerCPot  =0.0_dp ; DerDPot  =0.0_dp
        DivDpot  =0.0_dp
      endif  
      ! Define the preconditioning constant      
      aeff = -preconu
    endif 
    !---------------------------------------------------------------------------
    ! Decide on the way to calculate the action of the kinetic field B_mn.
    ! Done here for the ifort compiler
    if (BStack) then
      if(t1n2 .ne. 0.0_dp .or. t2n2 .ne. 0.0_dp) then
        ! Derivatives are stacked multiple times for the N2LO functional, 
        ! thus it is only safe to calculate it with Lagrange derivatives.
        call stp('Cannot stack derivatives for the B-field with N2LO.')
      else
        ActionOfB=> ActionOfBOld
      endif
    else
      if(t1n2 .ne. 0.0_dp .or. t2n2 .ne. 0.0_dp) then
        ActionOfB=> ActionOfBN2LO
      else
        ActionOfB=> ActionOfBNew
      endif
    endif
    !---------------------------------------------------------------------------    
    !Calculate the Coulomb Potential
    call SolveCoulomb()
    !---------------------------------------------------------------------------
    !Construct the mean-field potentials
    UpotOld = Upot
    Upot    = CalcUPot()!This includes the contribution of constraints & Coulomb
    if (allocated(WpotOld)) WpotOld = Wpot
    call CalcWPot()
    call CalcBPot()  
    
    if(.not.TRC) then 
        SpotOld = Spot
        Spot = CalcSPot()
        call CalcAPot()
        call CalcCPot()
        call CalcDPot()
    endif
    !---------------------------------------------------------------------------
    if(t1n2 .ne. 0.0_dp .or. t2n2 .ne. 0.0_dp) then
        ! All of the N2LO potentials
        call CalcBMunu()
        call CalcXpot()
        call calcDN2LO()
        call calcTfield()
        if(.not. TRC) then
            call calcPifield()
            call calcSN2LOField()
        endif
    endif
    !---------------------------------------------------------------------------
    ! Precondition the Upotential and Spotential
    !---------------------------------------------------------------------------
    if(MixingScheme .eq. 2 .and. (.not. all(UpotOld.eq.0.0))) then
    
      P = 0          ; S=0           ; TS=0
      if(PC) P = 1   ; if(SC) S = 1  ; if(TSC) TS = 1
    
      alpha =  aeff
      beta  =  1.0
      
      du    =  Upot-Upotold
      update=  PreconditionPotential(du,alpha,beta, P,S,TS,1)
      Upot  =  UpotOld + update
      
      if(.not. TRC) then
        ds = Spot - SpotOld
        ! First component
        update=  PreconditionPotential(ds(:,:,:,1,:),alpha,beta,P,-S,TS,1)
        Spot(:,:,:,1,:)  = Spotold(:,:,:,1,:) + update
        ! Second component
        update=  PreconditionPotential(ds(:,:,:,2,:),alpha,beta,P,-S,TS,2)
        Spot(:,:,:,2,:)  = Spotold(:,:,:,2,:) + update
        ! Third component
        update=  PreconditionPotential(ds(:,:,:,3,:),alpha,beta,P, S,TS,1)
        Spot(:,:,:,3,:)  = Spotold(:,:,:,3,:) + update

      endif

      if ( PreConW .gt. 0.0_dp ) then
        ! Define the preconditioning constant
      ! alpha = -PreConW
      ! beta  =  1.0

        dw = Wpot-Wpotold

        ! the following is organized by components having the same symmetry
        update=  PreconditionPotential(dw(:,:,:,1,1,:),alpha,beta,-P, S,TS,2)
        Wpot(:,:,:,1,1,:)  =  WpotOld(:,:,:,1,1,:) + update
        update=  PreconditionPotential(dw(:,:,:,2,2,:),alpha,beta,-P, S,TS,2)
        Wpot(:,:,:,2,2,:)  =  WpotOld(:,:,:,2,2,:) + update
        update=  PreconditionPotential(dw(:,:,:,3,3,:),alpha,beta,-P, S,TS,2)
        Wpot(:,:,:,3,3,:)  =  WpotOld(:,:,:,3,3,:) + update

        update=  PreconditionPotential(dw(:,:,:,1,2,:),alpha,beta,-P, S,TS,1)
        Wpot(:,:,:,1,2,:)  =  WpotOld(:,:,:,1,2,:) + update
        update=  PreconditionPotential(dw(:,:,:,2,1,:),alpha,beta,-P, S,TS,1)
        Wpot(:,:,:,2,1,:)  =  WpotOld(:,:,:,2,1,:) + update

        update=  PreconditionPotential(dw(:,:,:,1,3,:),alpha,beta,-P,-S,TS,2)
        Wpot(:,:,:,1,3,:)  =  WpotOld(:,:,:,1,3,:) + update
        update=  PreconditionPotential(dw(:,:,:,3,1,:),alpha,beta,-P,-S,TS,2)
        Wpot(:,:,:,3,1,:)  =  WpotOld(:,:,:,3,1,:) + update

        update=  PreconditionPotential(dw(:,:,:,2,3,:),alpha,beta,-P,-S,TS,1)
        Wpot(:,:,:,2,3,:)  =  WpotOld(:,:,:,2,3,:) + update
        update=  PreconditionPotential(dw(:,:,:,3,2,:),alpha,beta,-P,-S,TS,1)
        Wpot(:,:,:,3,2,:)  =  WpotOld(:,:,:,3,2,:) + update
      endif
    endif    
    !---------------------------------------------------------------------------
  end subroutine ConstructPotentials

end module PotentialMixing
