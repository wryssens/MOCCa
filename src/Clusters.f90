!-------------------------------------------------------------------------------
! A program that takes an input file with wavefunctions and copies them into
! different locations into the box.
!
!
!-------------------------------------------------------------------------------
! Note that:
! 1) The input file should have the correct symmetries and the correct nx/ny/nz.
!    This program can not change either of them. Use Ambrosya for symmetries.
!-------------------------------------------------------------------------------
!TODO
!
! - More input checking
!  -> Correct partners
!  -> correct symmetries
!
!-------------------------------------------------------------------------------
program Clusters

  use CompilationInfo
  use GenInfo
  use Wavefunctions
  use SpwfStorage
  use Moments
  use Densities  
  use InOutput       , only :  WFInput, InputFileName, output, OutputFileName
  use Transform      , only :  SymDisplace, Displace
  
  
  implicit none
  
  200 format (/,'___________________________________________________________', &
     &       /,'|                                                          |', &
     &       /,'|  Clusters - MOCCa Auxiliary                              |', &
     &       /,'|                                                          |', &
     &       /,'|  Copyright  P.-H. Heenen, M.Bender & W. Ryssens          |', &
     &       /,'|                 ; ,                                      |', &
     &       /,'|                ) ;( (                                    |', &
     &       /,'|               ( (  ) ;                                   |', &
     &       /,'|                ,-"""-.                                   |', &
     &       /,"|             ,-|`-...-'|                                  |", & 
     &       /,'|            ((_|       |                                  |', &   
     &       /,'|             `-\       /                                  |', &
     &       /,"|                `.___.'                                   |", &
     &       /,'|                                                          |', &
     &       /,'|__________________________________________________________|')
     
  1 format ( /, ' Cluster wavefunctions on the input file.') 
  2 format ( /, ' Wavefunctions on output.')
  3 format ( ' Cluster nr. ', i2)  
  4 format ( ' Location of the COM :', 3i3)  
  5 format ( ' Symmetrised: ', 2A3)
  
  integer,allocatable     :: Symmetrised(:,:), ClusterLocations(:,:)
  integer                 :: NumberOfClusters, Symx,SymY,SymZ, i,j, HFIndex,io
  integer                 :: newnwt,oldnwt, wave, X,Y,Z, l
  type(Spwf), allocatable :: Originals(:), Temp(:), Work(:)
  character(len=3)        :: Dir(3) = (/ ' X ', ' Y ', ' Z '/)
  
  NameList /ClusterInfo/ NumberOfClusters, InputFileName, nwt, OutputFileName
  NameList /Cluster    / X,Y,Z,SymX, SymY, SymZ
  
  print 200
  
  !-----------------------------------------------------------------------------
  ! Reading input
  call ReadGenInfo()

  read(NML=ClusterInfo, unit=*)

  allocate(ClusterLocations(3,NumberOfClusters)); ClusterLocations=0
  allocate(Symmetrised(3,NumberOfClusters))     ; Symmetrised=0
  do i=1,NumberOfClusters
    X = 0
    Y = 0
    Z = 0
    SymX = 0
    SymY = 0
    SymZ = 0
    read(NML=Cluster, unit=*,iostat=io)
    if(io.eq.0) then
      ClusterLocations(1,i) = X
      ClusterLocations(2,Y) = i
      ClusterLocations(3,i) = Z
      Symmetrised(1,i)      = SymX
      Symmetrised(2,i)      = SymY
      Symmetrised(3,i)      = SymZ
    else
      call stp('Read error for the Cluster namelists!')
    endif
  enddo
  !-----------------------------------------------------------------------------
  ! Reading the input wavefunctions.
  call WFInput()
  do i=1,nwt
    call HFBasis(i)%SymmetryOperators()
  enddo
  print 1
  !call PrintSpwf(1)
  !-----------------------------------------------------------------------------
  ! Print some info to the output.
  do i=1,NumberOfClusters
    print 3, i
    print 4, ClusterLocations(:,i)
    do j=1,3
      if(Symmetrised(j,i).ne.0) then
        print 5, Dir(j), 'Yes'
      else
        print 5, Dir(j), 'No '
      endif
    enddo
  enddo
  !-----------------------------------------------------------------------------
  ! Copying for further actions
  allocate(Originals(nwt))
  Originals = HFBasis 
  !-----------------------------------------------------------------------------
  !Making space for the new wavefunctions
  NewNwt = 0
  do i=1,NumberOfClusters
    !Space for the original cluster
    NewNwt  = NewNwt + nwt
    
    do j=1,3
      if(Symmetrised(j,i).eq.1) then
        if(ClusterLocations(j,i).ne.0) then
          !Space for one symmetrisation of the cluster
          NewNwt = NewNwt + nwt
        endif   
      endif
    enddo
  enddo  
  Oldnwt = nwt
  nwt    = NewNwt
  call ChangeNumberWaveFunctions(newnwt)
  !-----------------------------------------------------------------------------
  ! Big loop over clusters
  HFIndex=1
  do i=1,NumberOfClusters
    !---------------------------------------------------------------------------
    ! Loop over all wavefunctions in the cluster
    do wave=1,oldnwt
      allocate(Work(1)) ; Work(1) = Originals(wave)
      !-------------------------------------------------------------------------
      ! Loop over cartesian directions
      do j=1,3
        if(Symmetrised(j,i).eq.1 .and. ClusterLocations(j,i).ne.0) then
          !---------------------------------------------------------------------
          ! Displace symetrically in this direction.
          allocate(Temp(size(Work)))
          do l=1, size(Work)
            Temp(l) = Work(l)
          enddo
          deallocate(Work) ; allocate(Work(2*Size(Temp)))
          do l=1, size(Temp)
            call SymDisplace(j,Temp(l), Work(l), work(size(Temp) + l),         &
                                                          ClusterLocations(j,i))
          enddo
          deallocate(Temp)
        else
          !---------------------------------------------------------------------
          !Just displace, don't symmetrise in this direction.
          do l=1,size(Work)
            call Displace(j,Work(l), ClusterLocations(j,i))
          enddo
        endif
        if(allocated(Temp))deallocate(Temp)
        !-------------------------------------------------------------------------     
      enddo
      do l=1,Size(Work)
        HFBasis(HFIndex) = Work(l)
        HFIndex = HfIndex + 1
      enddo
      deallocate(Work)
      !-------------------------------------------------------------------------
    enddo
    !---------------------------------------------------------------------------
  enddo
  if(HFIndex.ne.nwt+1) call stp('Wrong HFINDEX.', 'HFINDEX', HFINDEX, 'nwt', nwt)
  !-----------------------------------------------------------------------------
  !
  !call GramSchmidt
  do wave=1,nwt
    call HFBasis(wave)%SymmetryOperators()
  enddo
  print 2
  
  call Output
  call stp('Correct ending for Clusters.f90.')
end program Clusters
