module Testing
!    !-------------------------
!    ! Module containing various Testing Routines
!    !---------------------
    use CompilationInfo
    use GenInfo
    use WaveFunctions
    use SpwfStorage
    use InOutput
    use Densities
    use Derivatives
    use OptimizedDerivatives
    use Moments
    use Mesh
    use InOutput
    use Energy
    use Force
    use Coulomb
    use Pairing 
    use MeanFields
    use Spinors
    use Interfaces

    implicit none

    public

    contains

    subroutine TestJ2

      integer :: i,j
      real(KIND=dp) :: AngMom(6,2)

      print *, '< Jx >'
      do i=1,nwt
        j = i
        !do j=1,nwt
          AngMom = Angularmomentum(HFBasis(i), HFBasis(j),.false.,.false.,.false.,.false.)
          print *, i,j,AngMom(1,1),InproductSpinorReal(HFBasis(i)%GetValue(), AngMomOperator(HFBasis(j),1))
        !enddo
      enddo
      print *

      print *, '< Jy >'
      do i=1,nwt
        j= i
        !do j=1,nwt
          AngMom = Angularmomentum(HFBasis(i), HFBasis(j),.false.,.false.,.false.,.false.)
          print *, i,j,AngMom(2,1),InproductSpinorReal(HFBasis(i)%GetValue(), AngMomOperator(HFBasis(j),2))
        !enddo
      enddo
      print *

      print *, '< Jz >'
      do i=1,nwt
        j = i
        !do j=1,nwt
          AngMom = Angularmomentum(HFBasis(i), HFBasis(j),.false.,.false.,.false.,.false.)
          print *, i,j,AngMom(3,1),InproductSpinorReal(HFBasis(i)%GetValue(), AngMomOperator(HFBasis(j),3))
        !enddo
      enddo
      print *

      print *, '< Jx | T >'
      do i=1,nwt
        j = i
        !do j=1,nwt
          AngMom = Angularmomentum(HFBasis(i), HFBasis(j),.false.,.true.,.true.,.true.)
          print *, i,j,AngMom(1,1),InproductSpinorReal(TimeReverse(HFBasis(i)%GetValue()), AngMomOperator(HFBasis(j),1))
        !enddo
      enddo
      print *

      print *, '< Jy | T >'
      do i=1,nwt
        j= i
        !do j=1,nwt
          AngMom = Angularmomentum(HFBasis(i), HFBasis(j),.false.,.true.,.true.,.true.)
          print *, i,j,AngMom(2,1),InproductSpinorReal(TimeReverse(HFBasis(i)%GetValue()), AngMomOperator(HFBasis(j),2))
        !enddo
      enddo
      print *

      print *, '< Jx**2 >'
      do i=1,nwt
        j = i
        !do j=1,nwt
          AngMom = Angularmomentum(HFBasis(i), HFBasis(j),.true.,.true.,.true.,.true.)
          print *, i,j,AngMom(4,1),InproductSpinorReal(AngMomOperator(HFBasis(i),1), AngMomOperator(HFBasis(j),1))
        !enddo
      enddo
      print *

      print *, '< Jy**2 >'
      do i=1,nwt
        j = i
        !do j=1,nwt
          AngMom = Angularmomentum(HFBasis(i), HFBasis(j),.true.,.true.,.true.,.true.)
          print *, i,j,AngMom(5,1),InproductSpinorReal(AngMomOperator(HFBasis(i),2), AngMomOperator(HFBasis(j),2))
        !enddo
      enddo
      print *

      print *, '< Jz**2 >'
      do i=1,nwt
        j = i
        !do j=1,nwt
          AngMom = Angularmomentum(HFBasis(i), HFBasis(j),.true.,.true.,.true.,.true.)
          print *, i,j,AngMom(6,1),InproductSpinorReal(AngMomOperator(HFBasis(i),3), AngMomOperator(HFBasis(j),3))
        !enddo
      enddo
      print *


      stop

    end subroutine TestJ2



    subroutine OddTpot

      integer :: i,j,k,l,it

      !Testing S potential

      print *,'--------------------------------------'
      print *,'U potential....'
      print *,'--------------------------------------'
      do it=1,2
        do k=1,nz
          do j=1,ny
            do i=1,nx
              if(abs(cr8U(i,j,k,it) - Upot(i,j,k,it))/abs(cr8U(i,j,k,it)) .gt.1d-10) then
                print *, i,j,k,it,cr8U(i,j,k,it),Upot(i,j,k,it), &
                 & abs(cr8U(i,j,k,it) - Upot(i,j,k,it))/abs(cr8U(i,j,k,it))
              endif
            enddo
          enddo
        enddo
      enddo

      print *,'--------------------------------------'
      print *,'A potential....'
      print *,'--------------------------------------'
      do it=1,2
        do l=1,3
          do k=1,nz
            do j=1,ny
              do i=1,nx
                if(abs(cr8A(i,j,k,l,it) - Apot(i,j,k,l,it))/abs(cr8A(i,j,k,l,it)) .gt.1d-10) then
                  print *, i,j,k,l,it,cr8A(i,j,k,l,it),Apot(i,j,k,l,it), &
                   & abs(cr8A(i,j,k,l,it) - Apot(i,j,k,l,it))/abs(cr8A(i,j,k,l,it))
                endif
              enddo
            enddo
          enddo
        enddo
      enddo

      print *,'--------------------------------------'
      print *,'S potential....'
      print *,'--------------------------------------'
      do it=1,2
        do l=1,3
          do k=1,nz
            do j=1,ny
              do i=1,nx
                if(abs(cr8S(i,j,k,l,it) - Spot(i,j,k,l,it))/abs(cr8S(i,j,k,l,it)) .gt.1d-10) then
                  print *, i,j,k,l,it,cr8S(i,j,k,l,it),Spot(i,j,k,l,it), &
                   & abs(cr8S(i,j,k,l,it) - Spot(i,j,k,l,it))/abs(cr8S(i,j,k,l,it))
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
      stop
    end subroutine OddTpot

!    subroutine TestDiag
!    
!      real(KIND=dp) :: N(2), Trace(2)
!      real(KIND=dp) ::  S
!      integer       :: i,j,k,it, at, ii, jj
!      complex(KIND=dp) :: OldRho(nwt,nwt)
!      
!      OldRho = RhoHFB
!      
!      call constructHFBHamiltonian(Fermi, Delta, LNLambda)
!      print *, 'Positive Parity, positive signature,  Protons'

!      call DeriveAll()
!      
!    
!      do i=1,nwt
!        call HFBasis(i)%CompAngMoment()
!        call HFBasis(i)%SymmetryOperators()
!      enddo


!      do i=1,2*HFBSize
!        ii = mod(i-1,nwt)+1        
!        if(HFBasis(ii)%GetIsospin().ne.1) cycle
!        if(HFBasis(ii)%GetParity().ne.1) cycle
!        if(HFBasis(ii)%GetSignature().ne.1) cycle        
!        do j=1,2*HFBSize
!          jj = mod(j-1,nwt)+1
!          if(HFBasis(ii)%GetIsospin().ne.HFBasis(jj)%GetIsospin()) cycle
!          if(HFBasis(ii)%GetParity().ne.HFBasis(jj)%GetParity()) cycle
!          if(j.gt.nwt) then ! .and. .not. i.gt.nwt) then
!            if (HFBasis(ii)%GetSignature().eq.HFBasis(jj)%GetSignature()) cycle
!          else
!            if (HFBasis(ii)%GetSignature().ne.HFBasis(jj)%GetSignature()) cycle
!          endif
!          
!          write(*,'(f7.3)',advance='no'), real(HFBHamil(i,j))
!        enddo
!        print *
!      enddo
!      
!      print *
!      print *
!      print *, 'Positive Parity, negative signature,  Protons'

!      do i=1,2*HFBSize
!        ii = mod(i-1,nwt)+1        
!        if(HFBasis(ii)%GetIsospin().ne.1) cycle
!        if(HFBasis(ii)%GetParity().ne.1) cycle
!        if(HFBasis(ii)%GetSignature().ne.-1) cycle        
!        do j=1,2*HFBSize
!          jj = mod(j-1,nwt)+1
!          if(HFBasis(ii)%GetIsospin().ne.HFBasis(jj)%GetIsospin()) cycle
!          if(HFBasis(ii)%GetParity().ne.HFBasis(jj)%GetParity()) cycle
!          if(j.gt.nwt) then !.and. .not. i.gt.nwt) then
!            if (HFBasis(ii)%GetSignature().eq.HFBasis(jj)%GetSignature()) cycle
!          else
!            if (HFBasis(ii)%GetSignature().ne.HFBasis(jj)%GetSignature()) cycle
!          endif
!          
!          write(*,'(f7.3)',advance='no'), real(HFBHamil(i,j))
!        enddo
!        print *
!      enddo          
!      
!      print *
!      print *, 'Negative Parity, positive signature, Protons'
!      
!      do i=1,2*HFBSize
!        ii = mod(i-1,nwt)+1
!        
!        if(HFBasis(ii)%GetIsospin().ne.1) cycle
!        if(HFBasis(ii)%GetParity().ne.-1) cycle
!        if(HFBasis(ii)%GetSignature().ne.1) cycle 
!        do j=1,2*HFBSize
!          jj = mod(j-1,nwt)+1
!          if(HFBasis(ii)%GetIsospin().ne.HFBasis(jj)%GetIsospin()) cycle
!          if(HFBasis(ii)%GetParity().ne.HFBasis(jj)%GetParity()) cycle
!          if(j.gt.nwt ) then ! .and..not. i.gt.nwt) then
!            if (HFBasis(ii)%GetSignature().eq.HFBasis(jj)%GetSignature()) cycle
!          else
!            if (HFBasis(ii)%GetSignature().ne.HFBasis(jj)%GetSignature()) cycle
!          endif
!          
!          write(*,'(f7.3)',advance='no'), real(HFBHamil(i,j))
!        enddo
!        print *
!      enddo     
!      
!      print *
!      print *, 'Negative Parity, negative signature, Protons'
!      
!      do i=1,2*HFBSize
!        ii = mod(i-1,nwt)+1
!        
!        if(HFBasis(ii)%GetIsospin().ne.1) cycle
!        if(HFBasis(ii)%GetParity().ne.-1) cycle
!        if(HFBasis(ii)%GetSignature().ne.-1) cycle 
!        do j=1,2*HFBSize
!          jj = mod(j-1,nwt)+1
!          if(HFBasis(ii)%GetIsospin().ne.HFBasis(jj)%GetIsospin()) cycle
!          if(HFBasis(ii)%GetParity().ne.HFBasis(jj)%GetParity()) cycle
!          if(j.gt.nwt ) then !.and..not. i.gt.nwt) then
!            if (HFBasis(ii)%GetSignature().eq.HFBasis(jj)%GetSignature()) cycle
!          else
!            if (HFBasis(ii)%GetSignature().ne.HFBasis(jj)%GetSignature()) cycle
!          endif
!          
!          write(*,'(f7.3)',advance='no'), real(HFBHamil(i,j))
!        enddo
!        print *
!      enddo     
!      
!      
!      print *
!      print *, 'Positive Parity, positive signature, Neutrons'
!      
!      do i=1,2*HFBSize
!        ii = mod(i-1,nwt)+1
!        
!        if(HFBasis(ii)%GetIsospin().ne.-1) cycle
!        if(HFBasis(ii)%GetParity().ne. 1) cycle
!        if(HFBasis(ii)%GetSignature().ne.1) cycle 
!        do j=1,2*HFBSize
!          jj = mod(j-1,nwt)+1
!          if(HFBasis(ii)%GetIsospin().ne.HFBasis(jj)%GetIsospin()) cycle
!          if(HFBasis(ii)%GetParity().ne.HFBasis(jj)%GetParity()) cycle
!          if(j.gt.nwt) then ! .and. .not. i.gt.nwt) then
!            if (HFBasis(ii)%GetSignature().eq.HFBasis(jj)%GetSignature()) cycle
!          else
!            if (HFBasis(ii)%GetSignature().ne.HFBasis(jj)%GetSignature()) cycle
!          endif
!          write(*,'(f7.3)',advance='no'), real(HFBHamil(i,j))
!        enddo
!        print *
!      enddo 
!      
!      print *
!      print *, 'Positive Parity, negative signature, Neutrons'
!      
!      do i=1,2*HFBSize
!        ii = mod(i-1,nwt)+1
!        
!        if(HFBasis(ii)%GetIsospin().ne.-1) cycle
!        if(HFBasis(ii)%GetParity().ne. 1) cycle
!        if(HFBasis(ii)%GetSignature().ne.-1) cycle 
!        do j=1,2*HFBSize
!          jj = mod(j-1,nwt)+1
!          if(HFBasis(ii)%GetIsospin().ne.HFBasis(jj)%GetIsospin()) cycle
!          if(HFBasis(ii)%GetParity().ne.HFBasis(jj)%GetParity()) cycle
!          if(j.gt.nwt) then ! .and. .not. i.gt.nwt) then
!            if (HFBasis(ii)%GetSignature().eq.HFBasis(jj)%GetSignature()) cycle
!          else
!            if (HFBasis(ii)%GetSignature().ne.HFBasis(jj)%GetSignature()) cycle
!          endif
!          write(*,'(f7.3)',advance='no'), real(HFBHamil(i,j))
!        enddo
!        print *
!      enddo    
!      
!      print *
!      print *
!         
!      print *
!      print *, 'Negative Parity, positive signature Neutrons'
!      
!      do i=1,2*HFBSize
!        ii = mod(i-1,nwt)+1
!        
!        if(HFBasis(ii)%GetIsospin().ne.-1) cycle
!        if(HFBasis(ii)%GetParity().ne.-1) cycle
!        if(HFBasis(ii)%GetSignature().ne.1) cycle 
!        do j=1,2*HFBSize
!          jj = mod(j-1,nwt)+1
!          if(HFBasis(ii)%GetIsospin().ne.HFBasis(jj)%GetIsospin()) cycle
!          if(HFBasis(ii)%GetParity().ne.HFBasis(jj)%GetParity()) cycle
!          if(j.gt.nwt ) then !.and. .not. i.gt.nwt ) then
!            if (HFBasis(ii)%GetSignature().eq.HFBasis(jj)%GetSignature()) cycle
!          else
!            if (HFBasis(ii)%GetSignature().ne.HFBasis(jj)%GetSignature()) cycle
!          endif
!          write(*,'(f7.3)',advance='no'), real(HFBHamil(i,j))
!        enddo
!        print *
!      enddo    
!      
!       print *
!      print *, 'Negative Parity, negative signature Neutrons'
!      
!      do i=1,2*HFBSize
!        ii = mod(i-1,nwt)+1
!        
!        if(HFBasis(ii)%GetIsospin().ne.-1) cycle
!        if(HFBasis(ii)%GetParity().ne.-1) cycle
!        if(HFBasis(ii)%GetSignature().ne.-1) cycle 
!        do j=1,2*HFBSize
!          jj = mod(j-1,nwt)+1
!          if(HFBasis(ii)%GetIsospin().ne.HFBasis(jj)%GetIsospin()) cycle
!          if(HFBasis(ii)%GetParity().ne.HFBasis(jj)%GetParity()) cycle
!          if(j.gt.nwt ) then !.and. .not. i.gt.nwt ) then
!            if (HFBasis(ii)%GetSignature().eq.HFBasis(jj)%GetSignature()) cycle
!          else
!            if (HFBasis(ii)%GetSignature().ne.HFBasis(jj)%GetSignature()) cycle
!          endif
!          write(*,'(f7.3)',advance='no'), real(HFBHamil(i,j))
!        enddo
!        print *
!      enddo    
!      
!      print *, 'Single line'
!      print *, real(HFBHamil(1,:))
!      print *      

!      
!      call ComputePairingCutoffs(Fermi)
!      
!      print *, 'Pairingmu', pairingmu(:)
!      print *, 'Pairingcut', pairingcut
!      print *, CosineCut(4.86_dp,0.0_dp,1)
!      
!      N = HFBNumberOfParticles(Fermi, Delta, LNLambda)
!      print *, 'N', N
!      
!      trace = 0
!      do i=1,HFBSize
!        it =( HFBasis(i)%GetIsospin()+3)/2
!        Trace(it) = Trace(it) + real(RhoHFB(i,i))
!      enddo
!      print *, 'trace after', trace
!      
!      do j=1,nwt
!        do i=1,nwt
!          if(abs(RhoHFB(i,j)-OldRho(i,j)).gt.1d-8 .and. abs(RhoHFB(i,j)).gt.1d-4) then
!            print *, i,j,RhoHFB(i,j), OldRho(i,j),(RhoHFB(i,j)-OldRho(i,j))/OldRho(i,j)          
!          endif
!        enddo
!      enddo
!      
!      do i=1,2*HFBSize
!        print *, QuasiEnergies(i)
!      enddo

!      print *, 'Diagonal'
!      
!      do i=1,HFBSize-12
!        print *, HFBHamil(i,i+HFBSize+12)
!      enddo
!     
!    end subroutine TestDiag
!    
   subroutine TestPairingFields
   
     integer :: i ,j , k , it, P
     real(KIND=dp) :: N(2)
   
     call ComputePairingCutoffs(Fermi)
     call CompDensityFactor

     call HFBPairingField(PairingField,PairingFieldLN, Delta)

     !call ConstructHFBHamiltonian(Fermi, Delta, LNLambda)
     !call DiagonaliseHFBHamiltonian()
     !call ConstructHFBstate
     N = 0.0
     do it=1,2
      do P=1,2
        do i=1,blocksizes(P,it)
          N(it) = N(it) + RhoHFB(i,i,P,it)
        enddo   
      enddo 
    enddo
    print *, N
     
      do it=1,2
        print *
        print *, ' Isospin it = ', it
        print *, '-------------------------------'
        print *, ' Real part '
        print *, '-------------------------------'
        do k=1,nz
          do j=1,ny
            do i=1,nx
              if(abs(real(PairingField(i,j,k,it))-dr(i,j,k,it))/abs(dr(i,j,k,it)).gt. 0.00001_dp) then
                write (*,'( 4i3, 3e15.5)'), i,j,k,it,real(PairingField(i,j,k,it)), dr(i,j,k,it), &
                &        abs(real(PairingField(i,j,k,it)) - dr(i,j,k,it))/abs(dr(i,j,k,it))
              endif
            enddo
          enddo
        enddo
      
        print *, '-------------------------------'      
        print *, ' Imaginary part '
        print *, '-------------------------------'
        
        do k=1,nz
          do j=1,ny
            do i=1,nx
              if(abs(AImag(PairingField(i,j,k,it))-di(i,j,k,it))/abs(di(i,j,k,it)).gt. 0.00001_dp) then
                write (*,'( 4i3, 3e15.5)'), i,j,k,it,AImag(PairingField(i,j,k,it)), di(i,j,k,it), &
                &        abs(AImag(PairingField(i,j,k,it)) - di(i,j,k,it))/abs(di(i,j,k,it))
              endif
            enddo
          enddo
        enddo
      enddo
      print *
      
      if(Lipkin) then
        print  *, 'Testing LN modified pairing fields'
        do it=1,2
          print *
          print *, ' Isospin it = ', it
          
          print *, '-------------------------------'
          print *, ' Real part '
          print *, '-------------------------------'
          
          do k=1,nz
            do j=1,ny
              do i=1,nx
                if(abs(real(PairingFieldLN(i,j,k,it))-dlnr(i,j,k,it))/abs(dlnr(i,j,k,it)).gt. 0.00001_dp) then
                write (*,'(4i3, 3f10.5)'), i,j,k,it,real(PairingFieldln(i,j,k,it)), dlnr(i,j,k,it), &
                &        abs(real(PairingFieldln(i,j,k,it)) - dlnr(i,j,k,it))/abs(dlnr(i,j,k,it))
              endif
              enddo
            enddo
          enddo
          
          print *, '-------------------------------'      
          print *, ' Imaginary part '
          print *, '-------------------------------'
          
          do k=1,nz
            do j=1,ny
              do i=1,nx
                if(abs(AImag(PairingFieldLN(i,j,k,it))-dlni(i,j,k,it))/abs(dlni(i,j,k,it)).gt. 0.00001_dp) then
                  write (*,'( 4i3, 3e15.5)'), i,j,k,it,AImag(PairingFieldLN(i,j,k,it)), dlni(i,j,k,it), &
                  &        abs(AIMag(PairingFieldLN(i,j,k,it)) - dlni(i,j,k,it))/abs(dlni(i,j,k,it))
                endif
              enddo
            enddo
          enddo
 
        enddo
      endif
   end subroutine TestPairingFields
   
   subroutine TestDelta
   
     integer :: i ,j , k , it, P
     real(KIND=dp), allocatable :: copydelta(:,:,:,:)
     real(KIND=dp) :: Test(2)
          
     call ComputePairingCutoffs(Fermi)
         
     allocate(copydelta(hfbsize,hfbsize,2,2))
     CopyDelta = Delta
     call HFBGaps(Delta,DeltaLN,PairingField,PairingFieldLN,Gaps,ConstantGap)
     
     print *
     print *, 'Testing calculation of the gaps'
     print *
     
     print *, 'Real part, ordinary delta'
     
     do it=1,2
      do P=1,2
        do i=1, blocksizes(P,it)
          do j=1, blocksizes(P,it)
            if(abs(Delta(i,j,P,it) - copydelta(i,j,P,it))/abs(copydelta(i,j,P,it)).gt.1d-8) then
              print *, i,j, abs(Delta(i,j,P,it)), abs(copydelta(i,j,P,it)),  &
               & abs(Delta(i,j,P,it) - copydelta(i,j,P,it))/abs(copydelta(i,j,P,it))
            endif
          enddo
        enddo
      enddo
    enddo

    print *, 'Energy', HFBEnergy(Delta)

    print *, HFBNumberOfParticles(Fermi,Delta,LnLambda)
   end subroutine TestDelta


!    subroutine TestPotentials
!    
!      integer :: i,j,k,it,m,Loca(4)
!      real(KIND=dp) :: Val
!      
!      call DeriveAll()
!      call Densit(0)
!      call ConstructPotentials()
!    
!      if(any(BPot .ne. BPotTest)) then
!        print *, maxval(abs(Bpot - BpotTest))
!        call stp('Bpot!')
!      else
!        print *, 'Bpot is ok!'
!      endif
!  
!      if(any(WPot .ne. WPotTest)) then 
!        print *, maxval(abs(Wpot - WpotTest))
!        call stp('Wpot!')
!      else
!        print *, 'Wpot is ok!'
!      endif

!      if(any(UPot .ne. UPotTest)) then
!        Loca= MaxLoc(abs(Upot - UpotTest))
!        Val= maxval(abs(Upot - UpotTest))
!        print *, Upot(Loca(1), Loca(2), Loca(3), Loca(4)),Val
!        !call stp('Upot!')
!      else
!        print *, 'Upot is ok!'
!      endif      
!    end subroutine TestPotentials

!    subroutine ActionOfTime()
!      type(Spwf) :: WF
!      type(Spinor):: AA, AS
!      integer    :: i,j,k,l

!      call DeriveAll()
!      call densit(0)
!      call ConstructPotentials()
!      call DerivePotentials()
!      
!      if (all(APot.eq.0.0_dp)) then
!        call stp('Apot is zero!')
!      endif
!      
!      WF = HFBasis(1)
!      
!      AA = ActionOFA(WF)
!      AS = ActionOfS(WF)
!      
!      do l=1,4
!        do k=1,nz
!          do j=1,ny 
!            do i=1,nx
!              if(abs(AA%Grid(i,j,k,l,1) - cr8ActionOfA%Grid(i,j,k,l,1)) .gt. 1.0d-14) then
!                print *, abs(AA%Grid(i,j,k,l,1) - cr8ActionOfA%Grid(i,j,k,l,1)), AA%Grid(i,j,k,l,1), cr8ActionOfA%Grid(i,j,k,l,1)
!              endif
!            enddo
!          enddo
!        enddo
!      enddo
!      
!      print *, 'Action Of S'
!      
!            do l=1,4
!        do k=1,nz
!          do j=1,ny 
!            do i=1,nx
!              if(abs(AS%Grid(i,j,k,l,1) - cr8ActionOfS%Grid(i,j,k,l,1)) .gt. 1.0d-14) then
!                print *, abs(AS%Grid(i,j,k,l,1) - cr8ActionOfS%Grid(i,j,k,l,1)), AS%Grid(i,j,k,l,1), cr8ActionOfS%Grid(i,j,k,l,1)
!              endif
!            enddo
!          enddo
!        enddo
!      enddo
!      
!    
!    end subroutine ActionOfTime
!       
!    subroutine TestIsospinSort()
!    	integer, allocatable :: Indices(:)
!    	
!    	Indices = OrderSpwfsISO(1)
!    	print*, Indices
!    
!    end subroutine TestIsospinSort

!    subroutine DrawAllWaveFunctions
!        integer :: wave, i
!        character(len=3)  :: strwave
!        character(len=20) :: Title
!        do wave=1,nwt
!            write(strwave, '(i3)') wave
!            Title='DensityX.wv=' // trim(adjustl(strwave))
!            call DrawWaveFunction(wave, 1, Title)
!            call sleep(3)
!            Title='DensityY.wv=' // trim(adjustl(strwave))
!            call DrawWaveFunction(wave, 2, Title)
!            call sleep(3)
!            Title='DensityZ.wv=' // trim(adjustl(strwave))
!            call DrawWaveFunction(wave, 3, Title)
!            call sleep(3)
!        enddo

!    end subroutine DrawAllWaveFunctions

!    subroutine TestActionOfB
!        integer :: i,j,k,l, it
!        type(Spinor) :: ActionOfB1, ActionOfB2, Lap, Der(3)
!        real(Kind=dp):: a, b, Correction

!        call DeriveAll()
!        call Densit(0)
!        call ConstructPotentials()
!        it = (HFBasis(1)%GetIsospin() + 3)/2
!        print *, it
!        ActionOfB1 = ActionOfBNew(HFBasis(1))
!        ActionOfB2 = ActionOfBOld(HFBasis(1))
!        Der(1) = HFBasis(1)%GetDer(1)
!        Der(2) = HFBasis(1)%GetDer(2)
!        Der(3) = HFBasis(1)%GetDer(3)
!        Lap = HFBasis(1)%GetLap()
!        print *, hbm(it) , xm(it), xm(3)
!        b = hbm(it)/2.0_dp * (1.0_dp - xm(it)/xm(3))
!        print*, b
!        do l=1,4
!            do k=1,nz
!                do j=1,ny
!                    do i=1,nx
!                        a = abs(ActionOfB1%Grid(i,j,k,l,1) - ActionOfB2%Grid(i,j,k,l,1))/abs(ActionOfB1%Grid(i,j,k,l,1))
!                        Correction = BPot(i,j,k,it) * Lap%Grid(i,j,k,l,1) 
!                        Correction = Correction + Der(1)%Grid(i,j,k,l,1)*NablaBPot(i,j,k,1,it)
!                        Correction = Correction + Der(2)%Grid(i,j,k,l,1)*NablaBPot(i,j,k,2,it)
!                        Correction = Correction + Der(3)%Grid(i,j,k,l,1)*NablaBPot(i,j,k,3,it)
!                        if (a .gt. 1d-14) then
!                            print *, i,j,k,l, a, &
!             & ActionOfB1%Grid(i,j,k,l,1) + Correction,&
!             & ActionOfB2%Grid(i,j,k,l,1) + Correction
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo

!    end subroutine TestActionOfB

!    subroutine TestSignatureCutoff

!        integer:: i,j,k, it
!        real(Kind=dp) :: a

!        call deriveall()
!        call Densit(0)
!        call CalculateAllMoments

!        do k=1,nz
!            do j=1,ny
!                do i=1,nx
!                    do it=1,2
!                        a = abs(Cutoff(i,j,k,it) - Cutoff(nx-i+1,j,k,it))/abs(Cutoff(i,j,k,it))
!                        if(a .gt. 1d-14) then
!                            print *, i,j,k,it, a
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo

!    end subroutine TestSignatureCutoff

!     subroutine CheckForDensityIslands

!         integer :: i,j,k
!         call DeriveAll()
!         !Calculating the density
!         call densit

!         do k=1,nz-1
!             do j=1,ny-1
!                 do i=nx/2,nx-1
!                     if(any(rho(i,j,k,:) .le. rho(i+1,j+1,k+1,:))) then
!                         print *, i,j,k, rho(i,j,k,:), rho(i+1,j+1,k+1,:)
!                     endif
!                 enddo
!             enddo
!         enddo

!     end subroutine CheckForDensityIslands

!     subroutine CheckForEqual

!         integer:: i,j,k
!         real(Kind=dp) :: NormMatrix(nwt,nwt,2)

!         do i=1,nwt
!             do j=1,nwt
!                 NormMatrix(i,j,:) = Inproduct(HFBasis(i), HFBasis(j))
!             enddo
!         enddo

!         print *,NormMatrix(:,:,1)
!         print *, 'Imag' 
!         print *,NormMatrix(:,:,2)

!     end subroutine CheckForEqual

    
!     subroutine EquivalentLaplacians
!     !
! 	! Testing if Coulomb & standard Laplacians give the same answer.
! 	!
    
!     	real(KIND=dp) :: Test(nx,ny,nz), CResp(nx,ny,nz), NResp(nx,ny,nz)
!     	integer       :: i,j,k
    
!     	call AssignFDCoefs(1,1,1)
    
!     	print *, FDCoulomb
!     	print *, FDLap
    
!     	Test=OneR
    
!     	CResp = CoulombLaplacian(Test, nx,ny,nz,ParityInt,SignatureInt,TimeSimplexInt,1)
!       NResp = Laplacian_Central(Test,ParityInt,SignatureInt,TimeSimplexInt,1)
    
!     	do k=1,nz
!     		do j=1,ny
!     			do i=1,nx
!     				if(abs(CResp(i,j,k) - NResp(i,j,k)) .ge. 1.d-8) then
!     					print *, 'Nonequivalent Laplacians', i,j,k, CResp(i,j,k), NResp(i,j,k)
! 				endif
! 			enddo
! 		enddo
! 	enddo
    
!     end subroutine EquivalentLaplacians
    
!     subroutine SignatureBreakingCoulomb()
!     !
!     ! Subroutine that checks if the solution of Poissions equation breaks signature.
!     !
!     !
!     implicit none
    
!     real(KIND=dp) :: Copy(nx+2,ny+1,nz+1), DensityCopy(nx,ny,nz,2), Test(34), Ext(1), TestShort(17), CPotLap(34,17,17)
!     integer       :: i,j,k, C
    
! 	    call deriveall()
! 	    call Densit(0)
! 	    call CalculateAllMoments()
! 	    call PrintAllMoments()
! 	    call SolveCoulomb()
	    
! 	    do k=1,nz
! 	    	do j=1,ny
! 	    		do i=1,nx
! 	    			DensityCopy(i,j,k,:) = Rho(nx-i+1,j,k,:)
! 	    			if(any(abs(DensityCopy(i,j,k,:) - Rho(i,j,k,:)).ge.1d-8)) then
! 	    				print *, 'Density not invariant', i,j,k, DensityCopy(i,j,k,:), Rho(i,j,k,:)
!     				endif
!     			enddo
! 		enddo
! 	    enddo
! 	    DensityCopy(:,:,:,1)=CoulombLaplacian(Rho(:,:,:,1),nx,ny,nz,1,0,1,1,.true.,2,1,1,nx-1,ny-1,nz-1)
! 	    if(all(DensityCopy.eq.ZeroR)) stop
! 	    do k=1,nz
! 		print *, k
! 		do j=1,ny
! 			do i=1,nx
! 				if(abs(DensityCopy(i,j,k,1) - DensityCopy(nx - i + 1,j,k,1)).ge.1d-8) then
! 					print *, 'CoulombLap problem', i,j,k
! 				endif
! 			enddo
! 		enddo
! 	    enddo
! 	    print *, 'sizes in Test', size(CoulombPotential,1),size(CoulombPotential,2),size(CoulombPotential,3)
! 		do k=1,nz+1
! 			do j=1,ny+1
! 				do i=1,nx+2
! 					Copy(i,j,k) = CPotComplete(nx - i + 3, j, k)
! 				enddo
! 			enddo
! 		enddo

! 	open(unit=25, file='Coulomb.dat')
! 	do k=2,2
! 		do j=3,3
! 			do i=1,nx+2
! 				write(fmt='(i2,2f12.3)', unit=25),i, CPotComplete(i,j,k), Copy(i,j,k)
! 			enddo
! 		enddo
! 	enddo
!     	close(25)
    
! !    	C=0
! !	do k=1,nz+1
! !		do j=1,ny+1
! !			do i=1,nx+2
! !				if(abs(Copy(i,j,k) - CPotComplete(i,j,k)).ge.1.d-5) then
! !					print *, 'Error', i,j,k,Copy(i,j,k), CPotComplete(i,j,k), BoundaryConditions(i,j,k)
! !					C = C +1
! !				endif
! !				if (C.ge.50) stop	
! !			enddo
! !		enddo
! !	enddo
	
	
! 	!Testing the Coulomb_1D procedures for signature inversion
! 	Ext=ZeroR
! 	Test=ZeroR
! 	Test(2:33) = 4 * pi *e2*Rho(:,1,1,1)
	
! 	do i=1,34
! 		if(abs(Test(i) - Test(34-i+1)).ge.1d-8) then
! 			print *, 'Problem before ', i, Test(i), Test(34-i+1)
! 		else
! 			print *, 'No Problem before',i, Test(i), Test(34-i+1)
! 		endif
! 	enddo
	
	
! 	Test = Coulomb_1D(Test,34,Ext,1,0,2,33)
	
	
! 	do i=1,34
! 		if(abs(Test(i) - Test(34-i+1)).ge.1d-8) then
! 			print *, 'Problem ', i, Test(i), Test(34-i+1)
! 		else
! 			print *, 'No Problem',i, Test(i), Test(34-i+1)
! 		endif
! 	enddo
	
! 	call assignFDCoefs(1,1,1)
	
	
! 	CPOtComplete=ZeroR
! 	CPotComplete(2:33,1:16,1:16) = 4 * pi * e2 * Rho(:,:,:,1)
!       CPotComplete(2:32,1:16,1:16)=  CoulombLaplacian(CPotComplete(2:33,1:16,1:16),32,16,16,1,0,1,OneI, &
!                 & .true., 1,1,1,32,16,16)              
                
                
		
		
! !		!Debugging
! !		do k=1,mz
! !			do j=1,my
! !				do i=1,mx
! !					if(abs(CPotComplete(i,j,k) - CPotComplete(mx-i+1,j,k)).ge.1d-12) then
! !						print *, 'Starting Residual problem', i,j,k
! !					endif
! !				enddo
! !			ENDDO
! !		enddo                
! 	CPotLap=CpotComplete
! 	CPotComplete=ZeroR
! 	CPotComplete(2:33,1:16,1:16) = 4 * pi * e2 * Rho(:,:,:,1)
! 	CPotComplete(2:33,1:16,1:16)=Laplacian_Central(CpotComplete(2:33,1:16,1:16),1, 0, 1,1)
! 	!Debugging
! 		do k=1,mz
! 			do j=1,my
! 				do i=1,mx
! 					if(abs(CPotComplete(i,j,k) - CPotComplete(mx-i+1,j,k)).ge.1d-12) then
! 						print *, 'Final', i,j,k
! 					endif
! 				enddo
! 			ENDDO
! 		enddo                

! 	!Debugging
! 		do k=1,mz
! 			do j=1,my
! 				do i=1,mx
! 					if(abs(CPotComplete(i,j,k) - CPotLap(i,j,k)).ge.1d-12) then
! 						print *, 'LapError', i,j,k, CPotComplete(i,j,k), CpotLap(i,j,k),CPotComplete(i,j,k)/CpotLap(i,j,k)
! 					endif
! 				enddo
! 			ENDDO
! 		enddo                

	
	
! !	
! !	TestShort=ZeroR
! !	TestShort(1:16)=4* pi *e2*Rho(1,:,1,1)
! !	Ext = Test(1)
! !	TestShort = Coulomb_1D(TestShort,17,Ext,1,1,1,16)
! !		
! !	do i=1,16
! !		if(abs(TestShort(i) - TestShort(34-i+1)).ge.1d-8) then
! !			print *, 'Problem ', i, TestShort(i), TestShort(16-i+1)
! !		else
! !			print *, 'No Problem',i, TestShort(i), TestShort(16-i+1)
! !		endif
! !	enddo
    
!     end subroutine SignatureBreakingCoulomb
    
    
!     subroutine PrintSomeDerivatives
    
!     	integer :: i, j,k,l
!     	type(Spinor) :: Lap, Value
    
    
    	
!     	call DeriveAll()
    	
!     	print *, 'Laplacian of wavefunction nr. 1'
!     	Lap = HFBasis(1)%GetLap()
!     	Value=HFBasis(1)%GetValue()
!     	do l=1,4
!     	do k=1,nz
!     		do j=1,ny
!     			do i=1,nx
!     				print *, i,j,k,l, Lap%Grid(i,j,k,l,1), Value%Grid(i,j,k,l,1)
! 			enddo
! 			print *,''
! 		enddo
! 		print *,''
! 	enddo
!       enddo
      
      
!     end subroutine PrintSomeDerivatives
    
    
!     subroutine TestAngularMomentum 
!     !
!     !
!     !
!     !
!     integer :: i,j,k,iwa,l,m, o
!     type(Spinor) :: Spin, Orbital, Value, Der(3), Temp2, Temp3
!     type(Spwf)   :: WF
!     real(KIND=dp):: orb,sp, orb2, orb3, OrbEigenFunction(nx,ny,nz,4,1), phi, theta, r
    
!     OrbEigenFunction=ZeroR
!     do k=1,nz
!     	do j=1,ny
!     		do i=1,nx
!     			r = MeshX(i)**2 + MeshY(j)**2 + MeshZ(k)**2
!     			theta = acos(MeshZ(k)/(MeshX(i)**2 + MeshY(j)**2 + MeshZ(k)**2))
!     			phi   = atan2(MeshY(j),MeshX(i))
! 			OrbEigenFunction(i,j,k,1,1) = exp(-r)!(2*MeshZ(k)**2 - MeshY(j)**2 - MeshX(i)**2)/r  *exp(-r)   !sin(theta)**2 * cos(2 * phi) * exp( - r )
! 			OrbEigenFunction(i,j,k,3,1) = exp(-r)!(2*MeshZ(k)**2 - MeshY(j)**2 - MeshX(i)**2)/r *exp(-r) !sin(theta)**2 * cos(2 * phi) * exp( - r )     			
!     			OrbEigenFunction(i,j,k,2,1) = exp(-r)!(2*MeshZ(k)**2 - MeshY(j)**2 - MeshX(i)**2)/r *exp(-r)!sin(theta)**2 * sin(2 * phi) * exp( - r )
!     			OrbEigenFunction(i,j,k,4,1) = exp(-r)!(2*MeshZ(k)**2 - MeshY(j)**2 - MeshX(i)**2)/r *exp(-r) !sin(theta)**2 * sin(2 * phi) * exp( - r )
!     		enddo
!     	enddo		
!     enddo

!     OrbEigenFunction = OrbEigenFunction/sqrt((sum(OrbEigenFunction**2)*dv))
!     print *, sum(OrbEigenFunction**2)*dv
    
!     call HfBasis(2)%SetGrid(OrbEigenFunction)
    
!     call DeriveAll() 
    
!     do iwa=1,nwt
! 	    WF = HFBasis(iwa)
! 	    Value = WF%GetValue()
! 	    Spin = Pauli(Value, 3)
! 	    Spin = 1/TwoR * Spin
	    
! 	    Orbital = AngMomOperator(WF,3)
! 	    Orbital = Orbital - Spin
! !	    
! 	    sp = InproductSpinorReal(Value,Spin)
! 	    orb = InproductSpinorReal(Value,Orbital)
! !	    
! 	    do i=1,3
! 	      Der(i) = WF%GetDer(i)
! 	    enddo	

!           Temp2 = WF%GetDer(2)          
          
!           Temp2 = Mesh3Dx * Temp2
!           Temp3 = WF%GetDer(1)
!           Temp3 = Mesh3DY*Temp3	

! 	    Temp2 = Temp3 - Temp2
! 	    Temp2 = MultiplyI(Temp2)
	    
! 	    orb2 =  InproductSpinorReal(Value,Orbital)
! 	    orb3 = ZeroR
! 	    do k=1,nz
! 	    	do j=1,ny
! 	    		do i=1,nx
! 			orb3 = orb3 + MeshX(i)*TwoR*(Value%Grid(i,j,k,1,1) * Der(2)%Grid(i,j,k,2,1) + &
! 			&                            Value%Grid(i,j,k,3,1) * Der(2)%Grid(i,j,k,4,1))  &
! 			&           - MeshY(j)*TwoR*(Value%Grid(i,j,k,1,1) * Der(1)%Grid(i,j,k,2,1) + &
! 			&                            Value%Grid(i,j,k,3,1) * Der(1)%Grid(i,j,k,4,1))   
! 	    		enddo
! 	       enddo
! 	    enddo
! 	    print *, iwa, sp, orb, orb2,orb3*dv
!     enddo    
    
!     do i=1,nx
!     print *, MeshX(i), MeshY(i), MeshZ(i), Mesh3DX(i,i,i), Mesh3DY(i,i,i), Mesh3DZ(i,i,i)
!     enddo
    
!     end subroutine TestAngularMomentum
    
!     subroutine TestConstraintEnergy
!     !
!     !
!     !
!     !
!   1 format(2e15.8)
!   2 format(78('_'))  
!     real*8  :: EV8ConstraintEnergy(nx,ny,nz,2), rms(3),r2
!     integer :: i,j,k,l,m
!     type(Moment), pointer :: Current
    
!     call CompRho(ZeroR)
!     call CalculateAllMoments()
!     !Calculating constraint energy
!     call CalcConstraintEnergy()
 
!     call PrintAllMoments()
    
!     open (unit=22, file='ConstraintEnergy')
!     do k=1,nz
!       do j=1,ny
!             do i=1,nx
!                  read (22, 1) EV8ConstraintEnergy(i,j,k,:) 
!             enddo
!       enddo
!     enddo
    
!     print 2

!    l=1
!    m=1
    
   
  
!     do i=1,nx
!       print *,'N', ConstraintEnergy(i,l,m,1), EV8ConstraintEnergy(i,l,m,1), &
!    &   abs( (ConstraintEnergy(i,l,m,1) -  EV8ConstraintEnergy(i,l,m,1))/ EV8ConstraintEnergy(i,l,m,1))
!       print *,'p', ConstraintEnergy(i,l,m,2), EV8ConstraintEnergy(i,l,m,2), &
!    &  abs( (ConstraintEnergy(i,l,m,2) -  EV8ConstraintEnergy(i,l,m,1))/ EV8ConstraintEnergy(i,l,m,2))
!     enddo
!     print 2 
    
!     if(all(ConstraintEnergy.eq.Zeror)) print *, 'ConstraintEnergy Zero!'
    
    
    
!     !Print some stuff to compare with ev8
!     rms =zeroR
    
!     do k=1,nz
!     	do j=1,ny
!     		do i=1,nx
!     			r2 = MeshX(i)**2 + MeshY(j)**2 + MeshZ(k)**2
!     			rms(1:2) = rms(1:2) + r2 * rho(i,j,k,1:2)
!     		enddo
! 	enddo
!     enddo
!     rms = rms * dv
!     rms(3) = rms(1) + rms(2)
    
!     rms(1) = sqrt(rms(1)/neutrons)
!     rms(2) = sqrt(rms(2)/Protons)
!     rms(3) = sqrt(rms(3)/(neutrons + protons))
!     print *, 'Rms, n, p, t',rms
    
!     Current => FindMoment(2,0,.false.)
!     print*, 'Re Q20', Current%Value*2*sqrt(5.0d0/(16.d0 * pi)), sum(Current%ValueCut)*2*sqrt(5.0d0/(16.d0 * pi))
    
    
!     print*, 'Re Q20', sum(Current%SpherHarm(:,:,:) * Rho(:,:,:,1)) *dv, sum(Current%SpherHarm(:,:,:) * Rho(:,:,:,2)) *dv 
!     print *, ''
!     print *, Current%SpherHarm(4,5,2), (2 * MeshZ(2)**2 - MeshX(4)**2 - MeshY(5)**2)/2.d0 
    
    
!     Current => FindMoment(4,0,.false.)
!     print*, 'Re Q40', Current%Value*sqrt(9.0d0/(4*pi)), sum(Current%ValueCut)*sqrt(9.0d0/(4*pi))
!     Current => FindMoment(4,2,.false.)
!     print*, 'Re Q42', Current%Value*sqrt(9.0d0/(4*pi)), sum(Current%ValueCut)*sqrt(9.0d0/(4*pi))
!     Current => FindMoment(4,4,.false.)
!     print*, 'Re Q44', Current%Value*sqrt(9.0d0/(4*pi)), sum(Current%ValueCut)*sqrt(9.0d0/(4*pi))
    
    
!     print *, ''
!     print *, ''
    
!     print *, 
    
    
!     end subroutine TestConstraintEnergy
    
    
    
    


!     subroutine TestMOCCaDer
        
!         integer       :: i,j            
!         type(Spinor)  :: Temp
!         real(Kind=dp) :: Derivatives(nx,ny,nz,4,3,nwt), Laplacians(nx,ny,nz,4,nwt)    

!         print *, "--------------------------------------------------------------------------"
!         print *, " Comparing derivatives from file with the calculated ones!"
!         print *, "--------------------------------------------------------------------------"
        


!         !Save the derivatives from the previous code
!         do i=1,nwt
!                 do j=1,3
!                     Temp = HFBasis(i)%GetDer(j)
!                     Derivatives(:,:,:,:,j,i) = Temp%Grid(:,:,:,:,1)
!                 enddo
!                 temp = HFBasis(i)%GetLap()
!                 Laplacians(:,:,:,:,i) = Temp%Grid(:,:,:,:,1)
!         enddo
    
!         call DeriveAll()

!         do i=1,nwt
!                 do j=1,3
!                     Temp = HFBasis(i)%GetDer(j)
!                     if(any(Derivatives(:,:,:,:,j,i).ne.Temp%Grid(:,:,:,:,1))) print*, 'Der', i,j
!                 enddo
!                 temp = HFBasis(i)%GetLap()
!                 if(any(Laplacians(:,:,:,:,i).ne.Temp%Grid(:,:,:,:,1))) print*, 'Lap', i,j 
!         enddo


!     end subroutine TestMOCCaDer


!     subroutine TestAdd
            
!             real(Kind=dp)             :: Time0, Time1, Time2, Time3
!             type(Spinor)              :: Test, Test1, Test2
!             real(KIND=dp),allocatable :: RATest(:,:,:,:,:), RATest1(:,:,:,:,:), RATest2(:,:,:,:,:)
!             real(Kind=dp)             :: RTest(nx,ny,nz,4,1), RTest1(nx,ny,nz,4,1), RTest2(nx,ny,nz,4,1)
!             integer                   :: i

!             call NewSpinor(Test); call NewSpinor(Test1); call NewSpinor(Test2) 
!             allocate(RATest(nx,ny,nz,4,1)); allocate(RATest1(nx,ny,nz,4,1)); allocate(RATest2(nx,ny,nz,4,1))
!             call CPU_Time(Time0)
!             do i=1,1000
!                 Test = Test1 + Test2
!             enddo
!             call CPU_Time(Time1)
!             do i=1,1000
!                 RATest = RATest1 + RATest2
!             enddo
!             call CPU_Time(Time2)
!             do i=1,1000
!                 RTest = RTest1 + RTest2
!             enddo
!             call CPU_Time(Time3)

!             print*, Time3-Time2,Time2-Time1, Time1-Time0

!     end subroutine TestAdd

!     subroutine CompareHFBPotentials
!         !--------------------------------------------------------------------------------------------------------------------------------------------
!         ! Subroutine testing if the HFB potentials are the same as the ones read from. Note that the ....TEST potentials need to have been allocated!
!         !--------------------------------------------------------------------------------------------------------------------------------------------

!         call DeriveAll()
!         call Densit(0)
!         call ConstructPotentials()
!         call CompEnergy()
!         call PrintEnergy()

!         print *, "--------------------------------------------------------------------------"
!         print *, " Comparing potentials from file with the calculated ones!"
!         print *, "--------------------------------------------------------------------------"
        
!        if(any(BPot.ne.BPotTest)) print*, "BPot is not correct!" 
!        if(any(NablaBPot.ne.NablaBPotTest)) print*, "NablaBPot is not correct!"
!        if(any(UPot(:,:,:,1).ne.UPotTest(:,:,:,1))) print*, "Upot is not correct for neutrons!"
!        if(any(UPot(:,:,:,2).ne.UPotTest(:,:,:,2))) print*, "Upot is not correct for protons!"
!        if(any(Apot.ne.ApotTest)) print*, "APot is not correct!"
!        if(any(Spot.ne.SpotTest)) print*, "SPot is not correct!"
!        if(any(CPot.ne.CPotTest)) print*, "CPot is not correct!"
!        if(any(WPot.ne.WPotTest)) print*, "WPot is not correct!"
!        if(any(DPot.ne.DpotTest)) print*, "Dpot is not correct!"
!        if(any(DerCPot.ne.DerCPotTest)) print*, "DerCPot is not correct!"
!        if(any(DivDPot.ne.DivDPotTest)) print*, "DivDPot is not correct!" 

!        print*, maxval(abs(abs(ApotTest)), maxval(abs(abs(Apot))
!        print*, maxval(abs(abs(SPotTest)), maxval(abs(abs(Spot))       
!        print*, maxval(abs(abs(DivDPot)), maxval(abs(abs(DivDPotTest))


!     end subroutine CompareHFBPotentials
! !     subroutine TestTimeDer
! !            
! !            type(Spinor):: Test, Der(3)
! !            real(KIND=DP), allocatable :: allocTest(:,:,:), allocDer(:,:,:,:)
! !            real(Kind=dp)              :: realTest(nx,ny,nz), realDer(nx,ny,nz,3), Time0, Time1, Time2, Time3, Time4
! !            real(KIND=dp)              :: fixed(16,16,16), FixedDer(16,16,16,3)
! !            integer                    :: i 

! !            allocate(allocTest(nx,ny,nz))
! !            allocate(allocDer(nx,ny,nz,3))
! !            call newSpinor(Test)
! !            allocTest=ZeroR
! !            allocDer=ZeroR
! !            realTest=ZeroR
! !            RealDer=ZeroR
! !            fixed=ZeroR
! !            FixedDer=ZeroR

! !            call CPU_Time(Time0)
! !            do i=1,100
! !                Der = DeriveSpinor(Test,1,1,1)
! !            enddo
! !            call CPU_Time(Time1)
! !            do i=1,100
! !                realDer = Central(RealTest,1,1,1,1)
! !            enddo
! !            call CPU_TIme(Time2)
! !            do i=1,100
! !                allocDer = Central(allocTest,1,1,1,1)
! !            enddo
! !            call CPU_Time(Time3)
! !            do i=1,100
! !                    FixedDer = Central(Fixed,1,1,1,1)
! !            enddo
! !            call CPU_Time(Time4)

! !            print*, (Time1-Time0)/FourR, Time2-Time1, Time3-Time2, Time4-Time3

! !     end subroutine TestTimeDer


! !    
! !    subroutine TestSort
! !        !----------------
! !        ! Small subroutine testing and timing the sorting of Spwf according to energy
! !        !----------------
! !        
! !        real(KIND=dp)    :: Time1, Time0
! !        type(Spwf)       :: HFBasisCopy(nwt), BubbleCopy(nwt)
! !        
! !        print*, "------------------------------------------------"
! !        print*, " Starting TestSort"
! !        print*, "------------------------------------------------"
! !        
! !        HFBasisCopy = HFBasis
! !        BubbleCopy = HFBasis
! !            
! !        call CPU_Time(Time0)
! !        call QuickSort(HFBasisCopy, nwt,1, nwt)
! !        call CPU_Time(Time1)
! !    
! !        print*, "QuickSort takes: ", Time1-Time0
! !        
! !        call CPU_Time(Time0)
! !        call Sort
! !        call CPU_Time(Time1)
! !    
! !        print*, "Normal Sort takes: ", Time1-Time0
! !        
! !        call CPU_Time(Time0)
! !        call BubbleSort(BubbleCopy,nwt)
! !        call CPU_Time(Time1)
! !    
! !        print*, "Bubble Sort takes: ", Time1-Time0
! !        
! !        
! !        call CPU_Time(Time0)
! !        call QuickSort(HFBasisCopy, nwt,1, nwt)    
! !        call CPU_Time(Time1)
! !        
! !        print*, "QuickSort on an already sorted list takes: ", Time1-Time0
! !        
! !        call CPU_Time(Time0)
! !        call Sort    
! !        call CPU_Time(Time1)
! !        
! !        print*, "Normal Sort on an already sorted list takes: ", Time1-Time0
! !        
! !        call CPU_Time(Time0)
! !        call BubbleSort(BubbleCopy,nwt)
! !        call CPU_Time(Time1)
! !    
! !        print*, "Bubble Sort on an already sorted list takes: ", Time1-Time0
! !        
! !    end subroutine TestSort
! !    
! !    
! !    subroutine TestHFFill
! !        !-----------------------------------------------------------------------------------
! !        ! Subroutine rapidly testing the filling of energy levels according to a HF scheme.
! !        !
! !        !-----------------------------------------------------------------------------------
! !    
! !        real(KIND=dp) :: Time1, Time0
! !    
! !        print*, "------------------------------------------------"
! !        print*, " Starting TestHFFill"
! !        print*, "------------------------------------------------"
! !    
! !    
! !        call CPU_Time(Time0)
! !        call Sort
! !        call CPU_Time(Time1)
! !        print*, "Sort takes: ", Time1-Time0
! !    
! !        print*, "Before HFFill"
! !        call PrintSpwf()
! !        
! !        call HFFill()
! !        
! !        print*, "========================================="
! !        print*, "After HFFill"
! !        call PrintSpwf
! !        
! !        print*, "------------------------------------------------"
! !        print*, " Ending TestHFFill"
! !        print*, "------------------------------------------------"
! !    
! !    end subroutine TestHFFill
! !    
!     subroutine TestOrtho
!     !-------------------
!     ! A subroutine designed for testing the orthogonalisation procedure Ortho.
!     !
!     !-------------------
    
!         real(KIND=dp) :: Time0, Time1, NormMatrix(nwt,nwt,2)
!         integer :: i,j
    
!         call DeriveAll()
!         call Densit(ZeroR)
!         !call SolveCoulomb()
!        ! call CompEnergy()
        
!         print*, "Total Energy before orthonormalisation: ", TotalEnergy
        
!         !call QuickSort(HFBasis, nwt,1,nwt)
        
!         print*,  "========================================================="
!         print*, "Places where the normmatrix is not sufficiently close to the identity matrix, on input"
!         do i=1,nwt
!                 do j=i,nwt
!                         NormMatrix(i,j,:) = Inproduct(HFBasis(i), HFBasis(j))
!                         if(i.eq.j .and. (abs(NormMatrix(i,j,1)-1).gt.0.00001)) then
!                                 print*, i,j, NormMatrix(i,j,1)
!                         elseif(i.ne.j .and. abs(NormMatrix(i,j,1)).gt.0.00001) then
!                                 print*, i,j, NormMatrix(i,j,1)
!                         endif                              
!                 enddo
!         enddo
!         print*,  "========================================================="
        
!         call PrintSpwf
!         call CPU_time(Time0)
!         call Ortho
!         call CPU_Time(Time1)
        
!         call HFFill
        
!         call DeriveAll()
!         call Densit(ZeroR) ! No Smoothing
!         !call SolveCoulomb()
!         !call CompEnergy()
        
!         call PrintSpwf
        
!         print*, "Time taken to orthonormalise ", Time1-Time0
        
!         print*, "Total Energy after orthonormalisation: ", TotalEnergy
        
!         print*,  "========================================================="
!         print*, "Places where the normmatrix is not sufficiently close to the identity matrix, on output"
!          do i=1,nwt
!                 do j=i,nwt
!                         NormMatrix(i,j,:) = Inproduct(HFBasis(i), HFBasis(j))
!                         if(i.eq.j .and. (abs(NormMatrix(i,j,1)-1).gt.0.00001)) then
!                                 print*, i,j, NormMatrix(i,j,1)
!                         elseif(i.ne.j .and. abs(NormMatrix(i,j,1)).gt.0.00001) then
!                                 print*, i,j, NormMatrix(i,j,1)
!                         endif                              
!                 enddo
!         enddo
!         print*,  "========================================================="
    
!     end subroutine TestOrtho
! !    
! !        
!     subroutine TestCutOff(it)
! !        ------------------------------
! !         This subroutine implements a naive test of the CompCutoff subroutine in the Moments module.
! !        ------------------------------
    
!         integer :: i,j,k, it,a
    
!         1 format ("Cut-off at (",i2,",",i2,",",i2,")", "=", f18.16) 
    
!         call CompRho(ZeroR)
!         call CompCutOff()
    
!         do i=1,nx
!                 print 1, i,1,1, Cutoff(i,1,1,it)
!         enddo
        
!         do a=1,2
!         do k=1,nz
!                 do j=1,ny
!                         do i=1,nx
!                                 if(CutOff(i,j,k,a).gt.OneR) then
!                                         print*, "Cutoff > 1 at (", i,j,k,") for it=",a
!                                         print*,  CutOff(i,j,k,a)
!                                 endif
!                         enddo
!                 enddo
!         enddo
!         enddo
    
!     end subroutine TestCutOff
! !    
! !    
! !    subroutine TestMesh()
! !        !--------------------------------------
! !        ! This subroutine prints out the values of the quantities in the module mesh. 
! !        ! Note that this uses the symmetries of the input.
! !        !--------------------------------------
! !    
! !        integer :: i,j,k
! !        
! !        1 format ("X", e15.8)
! !        2 format ("Y", e15.8)
! !        3 format ("Z", e15.8)
! !        
! !        print*, "------------------------------------------"
! !        print*, "Starting TestMesh"
! !        print*, "------------------------------------------"
! !        
! !        
! !        do i=1,nx
! !                print 1, MeshX(i)
! !        enddo
! !        
! !        print*
! !        
! !        do j=1,ny
! !                print 2, MeshY(j)
! !        enddo
! !        
! !        print*
! !        
! !        do k=1,nz
! !                print 3, MeshY(k)
! !        enddo
! !        
! !        return
! !    end subroutine TestMesh
! !      
! !    subroutine TestDerivingvsSumming()
! !        !-----------------------------
! !        ! This subroutine compares deriving densities as opposed to summing them.
! !        !  
! !        !      1) The derivative of rho is computed using two different methods:
! !        !               a) Calling Central with Rho as argument
! !        !               b) Calling GetDerRho
! !        !      2) Same analysis for the laplacian of Rho.
! !        !      
! !        !
! !        !-----------------------------
! !        real(KIND=dp)  :: DerRhoDer(nx,ny,nz,3,2),Time0,Time1,Time2,X,Y,Z, LapRhoDer(nx,ny,nz,2),A,B
! !        integer        :: i,k,l,j,m,it
! !        
! !        print*, "*******************************************"
! !        print*, "Starting Deriving vs Summing"
! !        print*, "*******************************************"
! !        
! !        print*, "-------------------------------------------"
! !        print*, "Step 1: Deriving and summing grad rho"
! !        print*, "-------------------------------------------"
! !        
! !        call PrintSpwf
! !        call DeriveAll()
! !        call cpu_time(Time0)
! !        call CompRho(ZeroR)
! !        call Cpu_time(Time1)
! !        do it=1,2
! !                DerRhoDer(:,:,:,:,it)=Central(Rho(:,:,:,it), nx,ny,nz,OneI,OneI,OneI,OneI) 
! !                LapRhoDer(:,:,:,it) = Laplacian_Central(Rho(:,:,:,it),nx,ny,nz,OneI,OneI,OneI,OneI)
! !        enddo
! !        call cpu_time(Time2)
! !        print*, "Timing"
! !        print*, "Sum", Time1-Time0
! !        print*, "Deriving", Time2-Time1
! !        
! !        
! !        print*, ""
! !        print*, "----------------------------------------------------------"
! !        print*, "1) \nabla \rho"
! !        print*, "----------------------------------------------------------"
! !        
! !        print*, "Neutrons"
! !        print*, "Sum", "                        ", "Derived", "                             ", "Difference"
! !        do i=1,nx
! !                do k=1,3
! !                        print*, DerRho(i,i,i,k,1), DerRhoDer(i,i,i,k,1), abs(DerRho(i,i,i,k,1) - DerRhoDer(i,i,i,k,1))
! !                enddo
! !                print*, ""
! !        enddo
! !        print*,""
! !        print*, "Protons" 
! !        print*, "Sum", "                         ", "Derived", "                             ", "Difference"
! !        do i=1,nx
! !                do k=1,3
! !                        print*, DerRho(i,i,i,k,2), DerRhoDer(i,i,i,k,2), abs(DerRho(i,i,i,k,2) - DerRhoDer(i,i,i,k,2))
! !                enddo
! !                print*,""
! !        enddo
! !        
! !        print*, "----------------------------------------------------------"
! !        print*, "2) \Delta \rho"
! !        print*, "----------------------------------------------------------"
! !        
! !        
! !        print*, "Neutrons"
! !        print*, "Sum", "                        ", "Derived", "                                 CR8        "&
! !        &       ,"      ", "Difference"
! !        do i=1,nx
! !               
! !                        print*, LapRho(i,i,i,1), LapRhoDer(i,i,i,1), LapRhoCR8(i,i,i,1), abs(LapRho(i,i,i,1) - LapRhoDer(i,i,i,1))
! !                print*, ""
! !        enddo
! !        print*,""
! !        print*, "Protons" 
! !        print*, "Sum", "                         ", "Derived", "                                 CR8        "&
! !        &       ,"                           ", "Difference"
! !        do i=1,nx
! !                
! !                        print*, LapRho(i,i,i,2), LapRhoDer(i,i,i,2),LapRhoCR8(i,i,i,2), abs(LapRho(i,i,i,2) - LapRhoDer(i,i,i,2))
! !                print*,""
! !        enddo
! !        
! !        A = sum(sum(Rho,4)*sum(LapRhoDer,4))*dv*B5
! !        B = sum(sum(Rho,4)*sum(LapRho,4))*dv*B5
! !       
! !        print*, "For example, the two different approaches result in an energy difference due to the B5."
! !        print*, "Deriving", A
! !        print*, "Summing", B
! !        print*, "Difference", (A-B)
! !       
! !    
! !        print*, "*******************************************"
! !        print*, "End of Deriving vs Summing"
! !        print*, "*******************************************"
! !        
! !                  
! !    end subroutine TestDerivingvsSumming


! !    subroutine TestSkyrmeTerms(EV8)
! !        !---------------------------------------------------------------------------------------
! !        ! This subroutine will test the calculation of Skyrme contributions to the total energy.
! !        !       EV8 determines the type of input
! !        !       1) EV8 input
! !        !       2) CR8 Input
! !        !---------------------------------------------------------------------------------------
! !        integer :: i,EV8
! !    
! !        print*, "*******************************************"
! !        print*, "Starting TestSkyrmeTerms"
! !        print*, "*******************************************"
! !        
! !        
! !        if(EV8.eq.OneI) then
! !                
! !        else
! !                call EnergyCR8Style
! !        endif
! !        
! ! 
! !        print*, "*******************************************"
! !        print*, "End Of TestSkyrmeTerms"
! !        print*, "*******************************************"

! !    end subroutine


! !!    subroutine TestLineExtension(wfnumber)
! !!        
! !!        real(KIND=dp), allocatable :: LineExtensionX(:),LineExtensionY(:),LineExtensionZ(:) 
! !!        integer :: wfnumber, Comp,i, parity, timesimplex, signature

! !!        print*, "*******************************************"
! !!        print*, "Starting TestLineExtension"
! !!        print*, "*******************************************"

! !!        call ReadCR8(ZeroI)
! !!        allocate(LineExtensionX(nx),LineExtensionY(ny),LineExtensionZ(nz) )

! !!        parity=HFBasis(wfnumber)%GetParity()
! !!        signature=HFBasis(wfnumber)%GetSignature()
! !!        timesimplex=HFBasis(wfnumber)%GetTimeSimplex()

! !!        print*, "QNumbers"
! !!        print*, "P=", Parity
! !!        print*, "TS=", TimeSimplex
! !!        print*, "S=", Signature

! !!        do Comp=1,4
! !!            print*, ""
! !!            print*, "Component", Comp
! !!            
! !!            LineExtensionX = PointToLineExtensionX(HFBasis(wfnumber)%GetGrid(Comp),nx,ny,nz,parity,timesimplex,signature,1,1)
! !!            LineExtensionY = PointToLineExtensionY(HFBasis(wfnumber)%GetGrid(Comp),nx,ny,nz,parity,timesimplex,signature,1,1)
! !!            LineExtensionZ = PointToLineExtensionZ(HFBasis(wfnumber)%GetGrid(Comp),nx,ny,nz,parity,timesimplex,signature,1,1)

! !!            do i=1,nx
! !!            print*, "i", i
! !!            print*, "X", LineExtensionX(i), HFBasis(wfnumber)%GetValue(i,1,1,Comp)        
! !!            print*, "Y", LineExtensionY(i), HFBasis(wfnumber)%GetValue(1,i,1,Comp)
! !!            print*, "Z", LineExtensionZ(i), HFBasis(wfnumber)%GetValue(1,1,i,Comp)
! !!            print*, ""
! !!            enddo
! !!        enddo

! !!        do Comp=1,4
! !!            print*, "Component", Comp
! !!            print*, "SignExt X,Y,Z", CompSignExtension(1,Parity,Signature,TimeSimplex,Comp) &
! !!                & ,CompSignExtension(2,Parity,Signature,TimeSimplex,Comp), &
! !!                &  CompSignExtension(3,Parity,Signature,TimeSimplex,Comp)

! !!        enddo


! !!        print*, "*******************************************"
! !!        print*, "End Of TestLineExtension"
! !!        print*, "*******************************************"

! !!    end subroutine TestLineExtension
! !!        
! !!    subroutine TestSingleTau()
! !!        integer :: i,j,k, wfnumber,l,m
! !!        real(KIND=dp)            :: CFDDer(nx,ny,nz,3,4), MFDDer(nx,ny,nz,3,4), LagDer(nx,ny,nz,3,4), random,GetTauResult(nx,ny,nz,3)
! !!        real(KIND=dp)            :: SumTest(nx,ny,nz)
! !!        print*, "*******************************************"
! !!        print*, "Starting TestSingleTau"
! !!        print*, "*******************************************"

! !!        
! !!!        call random_seed()
! !!!        call random_number(random)
! !!!        wfnumber=ceiling(random*nwt)
! !!        wfnumber=135
! !!        call ReadNil(ZeroI)

! !!        call HFBasis(wfnumber)%CompDer(ZeroI)
! !!        do l=1,4
! !!            CFDDer(:,:,:,1,l) = HFBasis(wfnumber)%GetDer(l,1)
! !!            CFDDer(:,:,:,2,l) = HFBasis(wfnumber)%GetDer(l,2)
! !!            CFDDer(:,:,:,3,l) = HFBasis(wfnumber)%GetDer(l,3)
! !!        enddo
! !!        GetTauResult(:,:,:,1) = HFBasis(wfnumber)%GetTau()

! !!        call HFBasis(wfnumber)%CompDer(OneI)
! !!        do l=1,4
! !!            MFDDer(:,:,:,1,l) = HFBasis(wfnumber)%GetDer(l,1)
! !!            MFDDer(:,:,:,2,l) = HFBasis(wfnumber)%GetDer(l,2)
! !!            MFDDer(:,:,:,3,l) = HFBasis(wfnumber)%GetDer(l,3)
! !!        enddo
! !!        GetTauResult(:,:,:,2) = HFBasis(wfnumber)%GetTau()
! !!        call HFBasis(wfnumber)%CompDer(TwoI)
! !!        do l=1,4
! !!            LagDer(:,:,:,1,l) = HFBasis(wfnumber)%GetDer(l,1)
! !!            LagDer(:,:,:,2,l) = HFBasis(wfnumber)%GetDer(l,2)
! !!            LagDer(:,:,:,3,l) = HFBasis(wfnumber)%GetDer(l,3)
! !!        enddo
! !!        GetTauResult(:,:,:,3) = HFBasis(wfnumber)%GetTau()

! !!        print*, "Testing wavefunction nr." , wfnumber
! !!        
! !!        
! !!        j=1
! !!        k=1

! !!        print*, ""
! !!        print*, "Testing Derivatives"


! !!        print*, "Component  Direction             CFD                   MFD                    Lag"
! !!        do i=1,10
! !!        do l=1,4
! !!            do m=1,3
! !!                print*, l,m, CFDDer(i,j,k,m,l)**2,MFDDer(i,j,k,m,l)**2,LagDer(i,j,k,m,l)**2
! !!                print*, ""
! !!            enddo
! !!        enddo
! !!        enddo
! !!        print*, ""
! !!        print*, "================================================="
! !!        print*, "Testing GetTau"
! !!        
! !!        do l=1,10
! !!            print*, GetTauResult(l,1,1,:)
! !!            print*, ""
! !!        enddo

! !!        

! !!    end subroutine TestSingleTau
! !!    
!     subroutine TestDensitiesCR8(NumberOfPoints)
!         !--------------------------------------------------------------------
!         ! This subroutine tests the computation of several densities.
!         ! A comparison is made with CR8 for a <NumberOfPoints> locations.
!         !
!         !    1) Normal density \rho
!         !    2) Kinetic energy density \tau
!         !    3) Nabla J
!         !--------------------------------------------------------------------
!    
!         use Spinors
!        
!         integer, intent(in)   :: NumberOfPoints
!         integer               :: i,j,k,l,m, iwa, NNN, mu,nu
!         real(KIND=dp)         :: MoccaDer(nx,ny,nz,4,3), Temp(nx,ny,nz,3), MoccaLag(nx,ny,nz,3,4), MoccaLap(nx,ny,nz,4)
!         real(KIND=dp)         :: TempLap(nx,ny,nz), CompareVecJ(nx,ny,nz,3,2)
!         type(Spinor)          :: TempSp, AA, AS, AU, AW
!    
!  			101     format("t0 =", es26.17 , " x0  =", es26.17 , /, &
!			&              't1 =', es26.17,  " x2  =", es26.17 , /, &
!			&              "t2 =", es26.17 , " x2  =", es26.17 , / ,&
!			&              "t3a=", es26.17 , " x3a =", es26.17 , " yt3a=", es26.17 ,/ ,&
!			&              "t3b=", es26.17 , " x3b =", es26.17 , ' yt3b=', es26.17 ,/, &
!			&              "te =", es26.17 , " to  =", es26.17 ,                  /, &
!			&              "wso=", es26.17 , " wsoq=", es26.17)     
!    
!    
!         print*, "*******************************************"
!         print*, "Starting TestDensitiesCR8"
!         print*, "*******************************************"
!        
!         print*, "Kind",kind(1._dp), kind(1.d0)
!         print*, 2._dp, 2.d0
!        
!         print*, "REMINDER!!!"
!         print*, "Machine Epsilon is ", epsilon(ZeroR)
!         print *, 'Force Used'
!         print 101 , t0, x0, t1, x1, t2, x2, t3a, x3a, yt3a, t3b, x3b, yt3b, te, to, wso, wsoq
!         
!         
!         
!         
        
 !        !Checking all the CR8 Derivatives
 !        print*, ""
 !        print*, "--------------------------------------------------------"
 !        print*, "Testing CR8 Derivatives"
 !        print*, "--------------------------------------------------------"
 !        do iwa=1,nwt

 !                
 !                if(HFBasis(iwa)%GetSignature().eq.-OneI) then
 !                        Temp = DerCR8(:,:,:,:,1,iwa)
 !                        DerCR8(:,:,:,:,1,iwa) = DerCR8(:,:,:,:,3,iwa)
 !                        DerCR8(:,:,:,:,3,iwa) = Temp
 !                        Temp = DerCR8(:,:,:,:,2, iwa)
 !                        DerCR8 (:,:,:,:,2,iwa) = DerCR8(:,:,:,:,4,iwa)
 !                        DerCR8 (:,:,:,:,4,iwa) = Temp
 !                endif
 !                
 !        
 !                do m=1,3
 !                
 !                TempSp = HFBasis(iwa)%GetDer(m)                    
 !                    do l=1,4
 !                        do k = 1,nz
 !                                do j = 1,ny
 !                                        do i=1,nx
 !!                                                print *, DerCR8(i,j,k,m,l,iwa)
 !                                                print *, TempSp%Grid(i,j,k,l,1)                                              
 !!                                        if( abs(TempSp%Grid(i,j,k,l,1) - DerCR8(i,j,k,m,l,iwa))   .ge. prec8  ) then
 !!                                                print*, "WFnumber", iwa, "Component", l, "Direction", m
 !!                                                print*, i,j,k
 !!                                                print*, abs(TempSp%Grid(i,j,k,l,1) - DerCR8(i,j,k,m,l,iwa))
 !!                                        endif
 !                                        enddo     
 !                                enddo
 !                        enddo
 !                    enddo
 !                enddo
 !        enddo  
 !        
!         print*, ""
!         print*, "--------------------------------------------------------"
!         print*, "Testing CR8 Laplacians"
!         print*, "--------------------------------------------------------"
             
         !Checking all the CR8 Laplacians
 !        do iwa=1,nwt
 !                
 !            TempSp =  HFBasis(iwa)%GetLap()                        
 !            MoccaLap =TempSp%Grid(:,:,:,:,1)                
 !            
 !!            if(HFBasis(iwa)%GetSignature().eq.-OneI) then
 !!                    TempLap = LapCR8(:,:,:,1,iwa)
 !!                    LapCR8(:,:,:,1,iwa) = LapCR8(:,:,:,3,iwa)
 !!                    LapCR8(:,:,:,3,iwa) = TempLap
 !!                    TempLap = LapCR8(:,:,:,2, iwa)
 !!                    LapCR8 (:,:,:,2,iwa) = LapCR8(:,:,:,4,iwa)
 !!                    LapCR8 (:,:,:,4,iwa) = TempLap
 !!            endif
 !                
 !            NNN= ZeroI
 !            do l=1,4             
 !            do k = 1,nz
 !                    do j = 1,ny
 !                            do i=1,nz
 !                                    
 !!                            if( (abs(MoccaLap(i,j,k,l) - LapCR8(i,j,k,l,iwa))   .ge. prec8).and.(i.eq.OneI)&
 !!                            & .and.(j.eq.OneI).and. (k.eq.OneI )) then
 !!                                    print*, "WFnumber", iwa, "Component", l 
 !!                                    print*, i,j,k
 !!                                    print*, abs(MoccaLap(i,j,k,l) - LapCR8(i,j,k,l,iwa))
 !!                                    NNN = NNN + 1
 !!                                    
 !!                                    if(NNN.ge.5*nx) then
 !!                                            call stp("Enough Errors")
 !!                                    endif
 !!                            endif
 !                            enddo     
 !                    enddo
 !            enddo
 !            enddo

 !        enddo       
 !       
!         call DeriveAll() 
!         call Densit(0)
!         call CalculateAllMoments()
!        
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of rho"
!         print*, "----------------------------------------"
!        
!        
!         do l=1, NumberOfPoints
!             i=1
!             j=1
!             k=l
!        
!             print*, "Isospin  MOCCa                CR8                MOCCa - CR8"
!             print*, "N", Rho(i,j,k,1), cr8rhon(i,j,k), abs(Rho(i,j,k,1) - cr8rhon(i,j,k)), abs(Rho(i,j,k,1)/cr8rhon(i,j,k))
!             print*, "P", Rho(i,j,k,2), cr8rhop(i,j,k), abs(Rho(i,j,k,2) - cr8rhop(i,j,k)), abs(Rho(i,j,k,2)/ cr8rhop(i,j,k))
!             print*, ""
!         enddo
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of lap(rho)"
!         print*, "----------------------------------------"
!        
!        
!         do l=1, NumberOfPoints
!             i=1
!             j=1
!             k=l
!        
!             print*, "Isospin  MOCCa                CR8                MOCCa - CR8"
!             print*, "N", LapRho(i,j,k,1), LapRhoCr8(i,j,k,1), abs(LapRho(i,j,k,1) - LapRhoCr8(i,j,k,1))& 
!             & , abs(LapRho(i,j,k,1)/LapRhoCr8(i,j,k,1))
!             print*, "P", LapRho(i,j,k,2), LapRhoCr8(i,j,k,2), abs(LapRho(i,j,k,2) - LapRhoCr8(i,j,k,2))&
!             & , abs(LapRho(i,j,k,2)/ LapRhoCr8(i,j,k,2))
!             print*, ""
!         enddo
!        
!        
!        
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of tau"
!         print*, "----------------------------------------"
!         do l=1, NumberOfPoints

!             i=1
!             j=l
!             k=1

!             print*, ""
!             Print*, "Checking at (i,j,k)=", i,j,k
!             print*, ""
!             print*, "Isospin       CFD                   CR8                     Difference" 
!             print*, "N", Tau(i,j,k,1), cr8vtau(i,j,k,1),abs(cr8vtau(i,j,k,1) - Tau(i,j,k,1))
!             print*, "P", Tau(i,j,k,2), cr8vtau(i,j,k,2),abs(cr8vtau(i,j,k,2) - Tau(i,j,k,2))
!             print*, ""
!         enddo
!        
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of nabla J"
!         print*, "----------------------------------------"


!    
!         do l=1, NumberOfPoints

!             i=1
!             j=1
!             k=l


!             print*, ""
!             Print*, "Checking at (i,j,k)=", i,j,k
!             print*, ""
!             print*, "Isospin       CFD                         CR8                   Difference" 
!             print*, "N", NablaJ(i,j,k,1), cr8vdiv(i,j,k,1), abs(cr8vdiv(i,j,k,1) - NablaJ(i,j,k,1))
!             print*, "P", NablaJ(i,j,k,2), cr8vdiv(i,j,k,2), abs(cr8vdiv(i,j,k,2) - NablaJ(i,j,k,2))
!             print*, ""                   
!         enddo
!        
!        
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of vec(j)"
!         print*, "----------------------------------------"
!        
!         do l=1,2
!             CompareVecJ(:,:,:,:,l) = B3 * sum( VecJ, 5) + B4 * Vecj(:,:,:,:,l)  
!         enddo

!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", CompareVecj(i,j,k,1,1), xr(i,j,k,1), abs(xr(i,j,k,1) - Comparevecj(i,j,k,1,1))
!                     print*, "P", Comparevecj(i,j,k,1,2), xr(i,j,k,2), abs(xr(i,j,k,2) - Comparevecj(i,j,k,1,2))
!                     print*, ""          

!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Comparevecj(i,j,k,2,1), yr(i,j,k,1), abs(yr(i,j,k,1) - Comparevecj(i,j,k,2,1))
!                     print*, "P", Comparevecj(i,j,k,2,2), yr(i,j,k,2), abs(yr(i,j,k,2) - Comparevecj(i,j,k,2,2))
!                     print*, ""         

!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Comparevecj(i,j,k,3,1), zr(i,j,k,1), abs(zr(i,j,k,1) - Comparevecj(i,j,k,3,1))
!                     print*, "P", Comparevecj(i,j,k,3,2), zr(i,j,k,2), abs(zr(i,j,k,2) - Comparevecj(i,j,k,3,2))
!                     print*, ""     

!                
!         enddo
!        
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of vec(s)"
!         print*, "----------------------------------------"
!        
!         print*, "X component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", vecs(i,j,k,1,1), rvx(i,j,k,1), abs(rvx(i,j,k,1) - vecs(i,j,k,1,1))
!                     print*, "P", vecs(i,j,k,1,2), rvx(i,j,k,2), abs(rvx(i,j,k,2) - vecs(i,j,k,1,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "Y component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", vecs(i,j,k,2,1), rvy(i,j,k,1), abs(rvy(i,j,k,1) - vecs(i,j,k,2,1))
!                     print*, "P", vecs(i,j,k,2,2), rvy(i,j,k,2), abs(rvy(i,j,k,2) - vecs(i,j,k,2,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "Z component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", vecs(i,j,k,3,1), rvz(i,j,k,1), abs(rvz(i,j,k,1) - vecs(i,j,k,3,1))
!                     print*, "P", vecs(i,j,k,3,2), rvz(i,j,k,2), abs(rvz(i,j,k,2) - vecs(i,j,k,3,2))
!                     print*, ""          
!                
!         enddo
!         
!        print*, "----------------------------------------"
!         print*, "Testing the calculation of rot(j)"
!         print*, "----------------------------------------"
!        
!         print*, "X component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", rotvecj(i,j,k,1,1), cr8rotj(i,j,k,1,1), abs(cr8rotj(i,j,k,1,1) - rotvecj(i,j,k,1,1))
!                     print*, "P", rotvecj(i,j,k,1,2), cr8rotj(i,j,k,1,2), abs(cr8rotj(i,j,k,1,2) - rotvecj(i,j,k,1,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "Y component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", rotvecj(i,j,k,2,1), cr8rotj(i,j,k,2,1), abs(cr8rotj(i,j,k,2,1) - rotvecj(i,j,k,2,1))
!                     print*, "P", rotvecj(i,j,k,2,2), cr8rotj(i,j,k,2,2), abs(cr8rotj(i,j,k,2,2) - rotvecj(i,j,k,2,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "Z component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", rotvecj(i,j,k,3,1), cr8rotj(i,j,k,3,1), abs(cr8rotj(i,j,k,3,1) - rotvecj(i,j,k,3,1))
!                     print*, "P", rotvecj(i,j,k,3,2), cr8rotj(i,j,k,3,2), abs(cr8rotj(i,j,k,3,2) - rotvecj(i,j,k,3,2))
!                     print*, ""          
!                
!         enddo
!         
!         
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of Jmunu"
!         print*, "----------------------------------------"
!        
!        
!         do l=1, NumberOfPoints
!             i=1
!             j=1
!             k=l
!             do mu=1,3
!              print *,' Mu = ' , mu
!                
!              do nu=1,3 
!               print *,' Nu = ' , nu
!               print*, "Isospin  MOCCa                CR8                MOCCa - CR8"
!               print*, "N", Jmunu(i,j,k,mu,nu,1) , cr8jmunu(i,j,k,mu,nu,1),              &
!               &            abs(Jmunu(i,j,k,mu,nu,1)- cr8jmunu(i,j,k,mu,nu,1))           &
!               &           ,Jmunu(i,j,k,mu,nu,1)/cr8jmunu(i,j,k,mu,nu,1)
!               print*, "P", Jmunu(i,j,k,mu,nu,2) , cr8jmunu(i,j,k,mu,nu,2),              &
!               &            abs(Jmunu(i,j,k,mu,nu,2)- cr8jmunu(i,j,k,mu,nu,2))           &
!               &           ,Jmunu(i,j,k,mu,nu,2)/ cr8jmunu(i,j,k,mu,nu,2)
!               print*, ""
!              enddo
!             ENDDO
!         enddo
!         
!         
!         
!         
!         
!         
!         call ConstructPotentials()
!         call DerivePotentials()
!         
!        print*, "----------------------------------------"
!         print*, "Testing the calculation of Upot"
!         print*, "----------------------------------------"
!        
!         
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Upot(i,j,k,1), wnn(i,j,k), abs(wnn(i,j,k) - Upot(i,j,k,1))
!                     print*, "P", Upot(i,j,k,2), wpp(i,j,k), abs(wpp(i,j,k) - Upot(i,j,k,2))
!                     print*, ""          
!                
!         enddo
!         
!         
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of Apot"
!         print*, "----------------------------------------"
!        
!         print*, "X component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Apot(i,j,k,1,1), pmx(i,j,k,1), abs(pmx(i,j,k,1) - Apot(i,j,k,1,1))
!                     print*, "P", Apot(i,j,k,1,2), pmx(i,j,k,2), abs(pmx(i,j,k,2) - Apot(i,j,k,1,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "Y component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Apot(i,j,k,2,1), pmy(i,j,k,1), abs(pmy(i,j,k,1) - Apot(i,j,k,2,1))
!                     print*, "P", Apot(i,j,k,2,2), pmy(i,j,k,2), abs(pmy(i,j,k,2) - Apot(i,j,k,2,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "Z component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Apot(i,j,k,3,1), pmz(i,j,k,1), abs(pmz(i,j,k,1) - Apot(i,j,k,3,1))
!                     print*, "P", Apot(i,j,k,3,2), pmz(i,j,k,2), abs(pmz(i,j,k,2) - Apot(i,j,k,3,2))
!                     print*, ""          
!                
!         enddo
!         
!         print *, 'byt3a', byt3a
!         
!         print*, "----------------------------------------"
!         print*, "Testing the calculation of Spot"
!         print*, "----------------------------------------"
!        
!         print*, "X component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Spot(i,j,k,1,1), vsx(i,j,k,1), abs(vsx(i,j,k,1) - Spot(i,j,k,1,1))
!                     print*, "P", Spot(i,j,k,1,2), vsx(i,j,k,2), abs(vsx(i,j,k,2) - Spot(i,j,k,1,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "Y component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Spot(i,j,k,2,1), vsy(i,j,k,1), abs(vsy(i,j,k,1) - Spot(i,j,k,2,1))
!                     print*, "P", Spot(i,j,k,2,2), vsy(i,j,k,2), abs(vsy(i,j,k,2) - Spot(i,j,k,2,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "Z component"
!         do l=1,NumberOfPoints
!                
!                 i=1
!                 j=1
!                 k=l
!                
!                 print*, ""
!                     Print*, "Checking at (i,j,k)=", i,j,k
!                     print*, ""
!                     print*, "Isospin       CFD                         CR8                   Difference" 
!                     print*, "N", Spot(i,j,k,3,1), vsz(i,j,k,1), abs(vsz(i,j,k,1) - Spot(i,j,k,3,1))
!                     print*, "P", Spot(i,j,k,3,2), vsz(i,j,k,2), abs(vsz(i,j,k,2) - Spot(i,j,k,3,2))
!                     print*, ""          
!                
!         enddo
!        
!         print*, "----------------------------------------"
!         print*, "Testing the action of Apot"
!         print*, "----------------------------------------"
!         
!         AA = ActionOfA(HFBasis(60))
!        
!         do l=1,NumberOfPoints
!          i=1;j=1;k=l
!          
!          print *, 'Checking at (i,j,k) = ', i,j,k
!          print *, ''
!          print*, "Comp       CFD                         CR8                   Difference" 
!          do m=1,4
!            print *, m, AA%Grid(i,j,k,m,1), CR8ActionOfA%Grid(i,j,k,m,1), abs(AA%Grid(i,j,k,m,1) - CR8ActionOfA%Grid(i,j,k,m,1))
!          enddo
!         
!         enddo
!         
!         print*, "----------------------------------------"
!         print*, "Testing the action of Spot"
!         print*, "----------------------------------------"
!         
!         AS = ActionOfS(HFBasis(60))
!        
!         do l=1,NumberOfPoints
!          i=1;j=1;k=l
!          
!          print *, 'Checking at (i,j,k) = ', i,j,k
!          print *, ''
!          print*, "Comp       CFD                         CR8                   Difference" 
!          do m=1,4
!            print *, m, AS%Grid(i,j,k,m,1), CR8ActionOfS%Grid(i,j,k,m,1), abs(AS%Grid(i,j,k,m,1) - CR8ActionOfS%Grid(i,j,k,m,1))
!          enddo
!         
!         enddo
!         
!         
!         print*, "----------------------------------------"
!         print*, "Testing the action of Upot"
!         print*, "----------------------------------------"
!         
!         AU = ActionOfU(HFBasis(60))
!        
!         do l=1,NumberOfPoints
!          i=1;j=1;k=l
!          
!          print *, 'Checking at (i,j,k) = ', i,j,k
!          print *,''
!          print*, "Comp       CFD                         CR8                   Difference" 
!          do m=1,4
!            print *, m, AU%Grid(i,j,k,m,1), CR8ActionOfU%Grid(i,j,k,m,1), abs(AU%Grid(i,j,k,m,1) - CR8ActionOfU%Grid(i,j,k,m,1))
!          enddo
!         
!         enddo
!         
!                  
!         print*, "----------------------------------------"
!         print*, "Testing the action of Wpot"
!         print*, "----------------------------------------"
!         
!         AW = ActionOfW(HFBasis(60))
!        
!         do l=1,NumberOfPoints
!          i=1;j=1;k=l
!          
!          print *, 'Checking at (i,j,k) = ', i,j,k
!          print *,''
!          print*, "Comp       CFD                         CR8                   Difference" 
!          do m=1,4
!            print *, m, AW%Grid(i,j,k,m,1), CR8ActionOfW%Grid(i,j,k,m,1), abs(AW%Grid(i,j,k,m,1) - CR8ActionOfW%Grid(i,j,k,m,1))
!          enddo
!         
!         enddo
!         
!         
!        
!         print*, "*******************************************"
!         print*, "End of TestDensitiesCR8"
!         print*, "*******************************************"
!        
!        
!     end subroutine TestDensitiesCR8

! !    subroutine TestDensitiesEV8(NumberOfPoints)
! !        !--------------------------------------------------------------------
! !        ! This subroutine tests the computation of several densities and their dependence on the derivation routines.
! !        ! A comparison is made with EV8, if possible, for a <NumberOfPoints> number of random locations.
! !        ! In addition, some timing experiments surrounding the derivative methods are executed.    
! !        !    
! !        !
! !        !    1) Normal density \rho
! !        !    2) Kinetic energy density \tau
! !        !    3) Nabla J
! !        !--------------------------------------------------------------------

! !        integer, intent(in)   :: NumberOfPoints
! !        real(KIND=dp)                :: CFDTau(nx,ny,nz,2), MFDTau(nx,ny,nz,2), LagTau(nx,ny,nz,2), random
! !        real(KIND=dp)                :: GetTauResult(nx,ny,nz,nwt,3)
! !        real(KIND=dp)                :: CFDNablaJ(nx,ny,nz,2), MFDNablaJ(nx,ny,nz,2), LagNablaJ(nx,ny,nz,2)
! !        real(KIND=dp)                :: time0,time1,time2,time3
! !        integer               :: i,j,k,l,m

! !        print*, "*******************************************"
! !        print*, "Starting TestDensitiesEV8"
! !        print*, "*******************************************"

! !!        call ReadNil(ZeroI)
! !        
! !        ev8rhon=2*ev8rhon
! !        ev8rhop=2*ev8rhop
! !        vtau=2*vtau
! !        vdiv=2*vdiv

! !        ! Computation of the densities using CFD Derivations.
! !        call cpu_time(time0)
! !        call iniden
! !        call DeriveAll(ZeroI)
! !        call densit(ZeroR)
! !        CFDTau(:,:,:,1) = Tau(:,:,:,1)
! !        CFDTau(:,:,:,2) = Tau(:,:,:,2) 
! !        CFDNablaJ(:,:,:,1) = NablaJ(:,:,:,1)
! !        CFDNablaJ(:,:,:,2) = NablaJ(:,:,:,2)
! !        call cpu_time(time1)

! !        ! Computation of the densities using MFD Derivations.
! !        call iniden
! !        call DeriveAll()
! !        call densit(ZeroR)
! !        MFDTau = Tau
! !        MFDNablaJ = NablaJ

! !        call cpu_time(time2)

! !        ! Computation of the densities using Lagrange Derivations.
! !        call iniden
! !        call DeriveAll(TwoI)
! !        call densit(ZeroR)
! !        LagTau = Tau
! !        LagNablaJ = NablaJ
! !        call cpu_time(time3)

! !        print*, "Timing of the methods"
! !        print*, "CFD", time1-time0
! !        print*, "MFD", time2-time1
! !        print*, "Lagrange", time3-time2

! !        !1) Checking \rho
! !        print*, "----------------------------------------"
! !        print*, "Testing the calculation of rho"
! !        print*, "----------------------------------------"

! !        do l=1, NumberOfPoints
! !!            call Random_Number(random)
! !!                    i=Ceiling(random*nx)
! !!            call Random_Number(random)
! !!            j=Ceiling(random*ny)
! !!            call Random_Number(random) 
! !!            k=Ceiling(random*nz)
! !!            

! !            i=1
! !            j=1
! !            k=l
! !        
! !            print*, "Isospin  MOCCa                EV8                MOCCa - EV8"
! !            print*, "N", Rho(i,j,k,1), ev8rhon(i,j,k), abs(Rho(i,j,k,1) - ev8rhon(i,j,k))
! !            print*, "P", Rho(i,j,k,2), ev8rhop(i,j,k), abs(Rho(i,j,k,2) - ev8rhop(i,j,k))
! !            print*, ""
! !        enddo



! !        !2) Checking \tau 
! !        print*, "----------------------------------------"
! !        print*, "Testing the calculation of tau"
! !        print*, "----------------------------------------"
! !        do l=1, NumberOfPoints
! !!            call Random_Number(random)
! !!                    i=Ceiling(random*nx)
! !!            call Random_Number(random)
! !!            j=Ceiling(random*ny)
! !!            call Random_Number(random) 
! !!            k=Ceiling(random*nz)
! !            
! !            i=1
! !            j=1
! !            k=l

! !            print*, ""
! !            Print*, "Checking at (i,j,k)=", i,j,k
! !            print*, ""
! !            print*, "Isospin       CFD                                  MFD                                     Lag" 
! !            print*, "N", CFDTau(i,j,k,1), MFDTau(i,j,k,1), LagTau(i,j,k,1)
! !            print*, "P", CFDTau(i,j,k,2), MFDTau(i,j,k,2), LagTau(i,j,k,2)
! !            print*, ""
! !            print*, "Isospin       EV8"
! !            print*, "N", vtau(i,j,k,1)
! !            print*, "P", vtau(i,j,k,2)
! !            print*, ""
! !            print*, "Isospin       CFD-Lag                              MFD-Lag"
! !            print*, "N", abs(CFDTau(i,j,k,1) - LagTau(i,j,k,1)), abs(MFDTau(i,j,k,1) - LagTau(i,j,k,1))
! !            print*, "P", abs(CFDTau(i,j,k,2) - LagTau(i,j,k,2)), abs(MFDTau(i,j,k,2) - LagTau(i,j,k,2))
! !            print*, "Isospin       CFD-EV8                              MFD-EV8"
! !            print*, "N", abs(vtau(i,j,k,1) - CFDTau(i,j,k,1)),abs(vtau(i,j,k,1) - MFDTau(i,j,k,1))
! !            print*, "P", abs(vtau(i,j,k,2) - CFDTau(i,j,k,2)),abs(vtau(i,j,k,2) - MFDTau(i,j,k,2))
! !            print*, ""

! !        enddo

! !        !3) Checking \nabla J
! !        print*, "----------------------------------------"
! !        print*, "Testing the calculation of nabla J"
! !        print*, "----------------------------------------"


! !    
! !        do l=1, NumberOfPoints
! !        !    call Random_Number(random)
! !        !      i=Ceiling(random*nx)
! !        !    call Random_Number(random)
! !        !    j=Ceiling(random*ny)
! !        !    call Random_Number(random) 
! !        !    k=Ceiling(random*nz)
! !            
! !            i=1
! !            j=1
! !            k=l


! !            print*, ""
! !            Print*, "Checking at (i,j,k)=", i,j,k
! !            print*, ""
! !            print*, "Isospin       CFD                                  MFD                                     Lag" 
! !            print*, "N", CFDNablaJ(i,j,k,1), MFDNablaJ(i,j,k,1), LagNablaJ(i,j,k,1)
! !            print*, "P", CFDNablaJ(i,j,k,2), MFDNablaJ(i,j,k,2), LagNablaJ(i,j,k,2)
! !            print*, ""
! !            print*, "Isospin       EV8"
! !            print*, "N", vdiv(i,j,k,1)
! !            print*, "P", vdiv(i,j,k,2)
! !            print*, ""
! !            print*, "Isospin       CFD-Lag                              MFD-Lag"
! !            print*, "N", abs(CFDNablaJ(i,j,k,1) - LagNablaJ(i,j,k,1)), abs(MFDTau(i,j,k,1) - LagTau(i,j,k,1))
! !            print*, "P", abs(CFDNablaJ(i,j,k,2) - LagNablaJ(i,j,k,2)), abs(MFDNablaJ(i,j,k,2) - LagNablaJ(i,j,k,2))
! !            print*, "Isospin       CFD-EV8                              MFD-EV8"
! !            print*, "N", abs(vdiv(i,j,k,1) - CFDNablaJ(i,j,k,1)),abs(vdiv(i,j,k,1) - MFDNablaJ(i,j,k,1))
! !            print*, "P", abs(vdiv(i,j,k,2) - CFDNablaJ(i,j,k,2)),abs(vdiv(i,j,k,2) - MFDNablaJ(i,j,k,2))

! !                        

! !        enddo
! !    
! !        print*, "*******************************************"
! !        print*, "End of TestDensitiesEV8"
! !        print*, "*******************************************"
! !    end subroutine TestDensitiesEV8

! !!    subroutine TestDerivatives()
! !!        !-------------------------------------------------------------------
! !!        ! This subroutine tests the validity of different derivatives.
! !!        !    1) Deriving all of the wavefunctions read by NIL8, random Direction and random component
! !!        !    2) Calculating all of the Laplacians.
! !!        !    
! !!        !-------------------------------------------------------------------

! !!        real(KIND=dp) :: random, CFD(nx,ny,nz), MFD(nx,ny,nz), Lag(nx,ny,nz)
! !!        real(KIND=dp) :: CFDLap(nx,ny,nz), MFDLap(nx,ny,nz), LagLap(nx,ny,nz)
! !!        integer:: wfnumber,WFsToCheck,Component,Direction,i,j

! !!        call ReadNil(ZeroI)
! !!        

! !!        Component=1
! !!        Direction=1
! !!        print*, "***********************************************"
! !!        print*, "Starting TestDerivatives"
! !!        print*, "***********************************************"

! !!        do j=1,nwt
! !!            wfnumber=j
! !!            print*, "-----------------------------------------------"
! !!            print*, "Wavefunction nr.", wfnumber, "Component", component,"Direction", Direction, " Will be Tested"
! !!            print*, "Quantum Numbers", HFBasis(wfnumber)%GetQNumbers()

! !!            ! Retrieving Central Finite Difference Values
! !!            call HFBasis(wfnumber)%CompDer(ZeroI)
! !!            CFD(:,:,:)=HFBasis(wfnumber)%GetDer(Direction,Component)
! !!            CFDLap(:,:,:)=HFBasis(wfnumber)%GetLap(Component)
! !!            !Retrieving Mixed Finite Difference Values        
! !!            call HFBasis(wfnumber)%CompDer(OneI)
! !!            MFD(:,:,:)=HFBasis(wfnumber)%GetDer(Direction,Component)
! !!            MFDLap(:,:,:)=HFBasis(wfnumber)%GetLap(Component)
! !!            !Retrieving Lagrange Mesh Derivatives
! !!            call HFBasis(wfnumber)%CompDer(TwoI)
! !!            Lag(:,:,:)=HFBasis(wfnumber)%GetDer(Direction,Component)
! !!            LagLap(:,:,:)=HFBasis(wfnumber)%GetLap(Component)
! !!            
! !!            !1) Printing all the derivatives
! !!            
! !!            print*, "First Order Derivatives"
! !!            Print*, ""
! !!            print*, "CFD                                MFD                                Lagrange Mesh"
! !!            print*, "*****X****"
! !!            do i=1,nx
! !!                print*, CFD(i,1,2), MFD(i,1,2), Lag(i,1,2)
! !!            enddo
! !!            print*, "*****Y****"
! !!            do i=1,ny
! !!                print*, CFD(1,i,2), MFD(1,i,2), Lag(1,i,2)

! !!            enddo
! !!            print*, "*****Z****"
! !!            do i=1,nz
! !!                print*, CFD(1,1,i), MFD(1,1,i), Lag(1,1,i)

! !!            enddo
! !!            
! !!            !2) Printing all the Laplacians


! !!            print*, "Laplacians Derivatives"
! !!            Print*, ""
! !!            print*, "CFD                                MFD                                Lagrange Mesh"
! !!            print*, "*****X****"
! !!            do i=1,nx
! !!                print*, CFDLap(i,1,2), MFDLap(i,1,2), LagLap(i,1,2)
! !!            enddo
! !!            print*, "*****Y****"
! !!            do i=1,ny
! !!                print*, CFDLap(1,i,2), MFDLap(1,i,2), LagLap(1,i,2)

! !!            enddo
! !!            print*, "*****Z****"
! !!            do i=1,nz
! !!                print*, CFDLap(1,1,i), MFDLap(1,1,i), LagLap(1,1,i)
! !!            enddo


! !!        enddo

! !!        print*, "**************************"
! !!        print*, "End of TestDerivatives"
! !!        print*, "**************************"
! !!        print*, ""
! !!        return
! !!    end subroutine TestDerivatives    

! !!    
! !!    subroutine TestGnuForInterface()
! !!        !--------------------------------------------------------------------
! !!        ! This subroutine tests if the PlotDensity commands work as expected.
! !!        !--------------------------------------------------------------------
! !!        call IniMesh
! !!        call ReadNil(ZeroI)
! !!        call iniden
! !!        call densit(ZeroR)
! !!        call PlotDensity(1,1,0)
! !!        
! !!    end subroutine
! !!    
! !    subroutine TestMomentGeneration(PlacesToCheck)
! !        !-------------------------------------------------------------------------------
! !        ! This Subroutine Tests the generation of the moments in the Moments module.
! !        ! It also tests if the correct values are assigned to their SpherHarm Variables.
! !        !-------------------------------------------------------------------------------
! !        
! !        real(KIND=dp) :: X,Y,Z,random,f
! !        integer       :: i,j,k,q,l,m
! !        integer       :: PlacesToCheck
! !        class(Moment),pointer :: Test
! !        logical       :: Imaginary

! !        PC = .false.
! !        TRC= .false.
! !        TSC= .false.
! !        SC = .false.

! !        print*, "Testing the Generation of moments"
! !        print*, "Symmetries Present:"
! !        print*, "TimeReversal", TRC, "Time Simplex", TSC
! !        print*, "Parity", PC, "Signature",SC

! !        call IniMesh
! !        call IniMoments
! !        
! !        call PrintAllMoments
! ! 
! !        call random_seed
! !        print*, ''
! !        print*, 'Checking the SpherHarm variables'
! !        do q=1,PlacesToCheck 
! !            call Random_Number(random)
! !            i=Ceiling(random*(nx))
! !            X=MeshX(i)
! !            call Random_Number(random)
! !            j=Ceiling(random*(ny))
! !            Y=MeshY(j)
! !            call Random_Number(random) 
! !            k=Ceiling(random*(nz))
! !            Z=MeshZ(k)
! !            Test => Root%Next
! !            Print*, "Checking at (i,j,k)=", i,j,k
! !                
! !            do while(associated(Test))
! !                
! !                l=Test%l
! !                m=Test%m
! !                Imaginary=Test%Impart
! !                print*, "l,m,Imaginary?", l,m,Imaginary
! !                
! !                !l=1 Spherical Harmonics
! !                if(l.eq.1) then
! !                    if(m.eq.1) then
! !                        if(Imaginary) then
! !                            f=y/(sqrt(TwoR))
! !                        else
! !                            f=x/(sqrt(TwoR))
! !                        endif
! !                    else
! !                        f=z
! !                    endif
! !                !l=2 Spherical Harmonics    
! !                elseif(l.eq.2) then
! !                    if(m.eq.2) then
! !                        if(Imaginary) then
! !                            f=TwoR*x*y*sqrt(ThreeR/8._dp)
! !                        else
! !                            f=(x**2-y**2)*sqrt(ThreeR/8._dp)
! !                        endif
! !                    
! !                    elseif(m.eq.1) then
! !                        if(Imaginary) then
! !                            f=z*y*sqrt(ThreeR/8._dp)*2
! !                        else
! !                            f=x*z*sqrt(ThreeR/8._dp)*2
! !                        endif
! !                        
! !                    else
! !                        f=(2*z**2-x**2-y**2)/2
! !                    endif
! !                    !l=3 Spherical Harmonics
! !                elseif(l.eq.3) then
! !                    select case (m)
! !                        case (3)
! !                            if(Imaginary)then
! !                                f=(3*x**2 - y**2)*y*sqrt(5._dp)/4._dp
! !                            else
! !                                f=(x**2 - 3*y**2)*x*sqrt(5._dp)/4._dp
! !                            endif
! !                        case (2)
! !                            if(Imaginary)then
! !                                f=TwoR*x*y*z/2*sqrt(15._dp/2)
! !                            else
! !                                f=(x**2-y**2)*z/2*sqrt(15._dp/2)
! !                            endif
! !                        case (1)
! !                            if(Imaginary)then
! !                                f=y*(4*z**2-x**2-y**2)*sqrt(ThreeR)/4._dp
! !                            else
! !                                f=x*(4*z**2-x**2-y**2)*sqrt(ThreeR)/4._dp
! !                            endif
! !                        case (0)
! !                            f = Z*(2*Z**2 - 3*X**2 - 3*Y**2)/TwoR      
! !                    end select
! !                
! !                else
! !                    call stp('This subroutine is not suited to test for l>3')
! !                                                    
! !                endif
! !                        
! !                print*, "          SpherHarm          ", " What It should be           ","        Difference"
! !                print*, Test%SpherHarm(i,j,k),f, abs(f-  Test%SpherHarm(i,j,k)) 
! !                Test => Test%Next
! !            enddo
! !        enddo

! !        call PrintAllMoments

! !        print*, "End of TestMomentGeneration"
! !    end subroutine TestMomentGeneration
! !!    
! !!    
! !    subroutine TestMomentCalculation
! !        !-----------------------
! !        ! A Subroutine that compares the calculation of moments with EV8.
! !        !-----------------------
! !        !Variables containing the calculated moments of EV8
! !        real(KIND=dp) :: EV8l2m0(3), EV8l2m2(3), EV8l2m0c(3), EV8l2m2c(3)
! !        real(KIND=dp) :: CR8l2m0(3), CR8l2m2(3), CR8l2m0c(3), CR8l2m2c(3)
! !        
! !        print*, "Testing the Calculation of moments"
! !        print*, "Symmetries Present:"
! !        print*, "TimeReversal", TRC, "Time Simplex", TSC
! !        print*, "Parity", PC, "Signature",SC
! !             
! !        call densit(ZeroR)
! !        
! !        call CalculateAllMoments
! !        
! !        call PrintAllMoments
! !        
! !        if(InputFileName(1:3).eq."ev8" ) then        
! !        
! !                EV8l2m0(1) = qznc
! !                EV8l2m0(2) = qzpc
! !                EV8l2m0(3) = (qznc+qzpc)
! !                
! !                EV8l2m2(1) = (qxnc - qync)/ThreeR*sqrt(ThreeR/8._dp)*TwoR
! !                EV8l2m2(2) = (qxpc - qypc)/ThreeR*sqrt(ThreeR/8._dp)*TwoR
! !                EV8l2m2(3) = EV8l2m2(1) + EV8l2m2(2)
! !                
! !                EV8l2m2 = EV8l2m2

! !                print*, ""
! !                print*, "EV8 Moments, with cutoff"
! !                print*, "L= 2, M= 0"
! !                print*, "N", EV8l2m0(1)
! !                print*, "P", EV8l2m0(2)
! !                print*, "T", EV8l2m0(3)
! !                print*, "L= 2, M= 2"
! !                print*, "N", EV8l2m2(1)
! !                print*, "P", EV8l2m2(2)
! !                print*, "T", EV8l2m2(3)
! !        
! !        elseif(InputFileName(1:3).eq."cr8") then
! !                CR8l2m0(1) = 2*z2n - x2n - y2n
! !                CR8l2m0(2) = 2*z2p - x2p - y2p
! !                CR8l2m0(3) = CR8l2m0(1) + CR8l2m0(2)
! !                
! !                CR8l2m2(1) = x2n - y2n
! !                CR8l2m2(2) = x2p - y2p
! !                CR8l2m2(3) = CR8l2m2(1) + CR8l2m2(2) 
! !        
! !                print*, ""
! !                print*, "CR8 Moments, without cutoff"
! !                print*, "L= 2, M= 0"
! !                print*, "N", CR8l2m0(1)/TwoR
! !                print*, "P", CR8l2m0(2)/TwoR
! !                print*, "T", CR8l2m0(3)/TwoR
! !                print*, "L= 2, M= 2"
! !                print*, "N", CR8l2m2(1)*sqrt(ThreeR/8._dp)
! !                print*, "P", CR8l2m2(2)*sqrt(ThreeR/8._dp)
! !                print*, "T", CR8l2m2(3)*sqrt(ThreeR/8._dp)
! !                
! !                CR8l2m0(1) = 2*z2nc - x2nc - y2nc
! !                CR8l2m0(2) = 2*z2pc - x2pc - y2pc
! !                CR8l2m0(3) = CR8l2m0(1) + CR8l2m0(2)
! !                
! !                CR8l2m2(1) = x2nc - y2nc
! !                CR8l2m2(2) = x2pc - y2pc
! !                CR8l2m2(3) = CR8l2m2(1) + CR8l2m2(2) 
! !                
! !                print*, ""
! !                print*, "CR8 Moments, with cutoff"
! !                print*, "L= 2, M= 0"
! !                print*, "N", CR8l2m0(1)/TwoR
! !                print*, "P", CR8l2m0(2)/TwoR
! !                print*, "T", CR8l2m0(3)/TwoR
! !                print*, "L= 2, M= 2"
! !                print*, "N", CR8l2m2(1)*sqrt(ThreeR/8._dp)
! !                print*, "P", CR8l2m2(2)*sqrt(ThreeR/8._dp)
! !                print*, "T", CR8l2m2(3)*sqrt(ThreeR/8._dp)
! !                
! !        endif
! !        
! !    end subroutine TestMomentCalculation
    
!     subroutine TestMomentConstraints
    
! !        real(KIND=dp) :: Q20Constraint(2), Q22Constraint(2), Q20Intensity, Q20Value(2)
! !        real(KIND=dp) :: TestCEnergy(nx,ny,nz,2)
! !        integer       :: it, i
! !        type(Moment), pointer :: Current=>null()
! !        
! !        call densit
! !        call CalculateAllMoments
! !        call CalcConstraintEnergy
! !    
! !        print*, "Testing Constraints on Multipole Moments"
! !        
! !        Current => FindMoment(2,0,.false.)
! !        Q20Constraint= Current%Constraint
! !        Q20Intensity = Current%Intensity
! !        Q20Value     = Current%Value
! !        
! !        do it=1,2
! !                TestCEnergy(:,:,:,it) = Q20Intensity*TwoR*Current%SpherHarm(:,:,:)*(Q20Value(it) - Q20Constraint(it) ) &
! !                & *Cutoff(:,:,:,it)
! !        enddo 
! !        
! !        do i=1,nx
! !                print*, i, ConstraintEnergy(i,1,1,1), TestCEnergy(i,1,1,1), ConstraintEnergy(i,1,1,1)/TestCEnergy(i,1,1,1)
! !                
! !        enddo
! !        
! !        !Current => FindMoment(2,2,.false.)
! !        !Q22Constraint = Current%Constraint
! !        
! !        
    
!     end subroutine TestMomentConstraints
    
    
!subroutine TestSignatureBreaking
!	implicit none
!	
!	integer :: iwa, s, i,j,k,l, C
!	type(Spinor) :: Value, Der(3), Lap, RLap, RDer(3), Test, U, B, W
!	type(Spwf)   :: Psi
!	
!	print *, '----------------------------------------------------------'
!	print *, ' Test what things break signature'
!	print *, '----------------------------------------------------------'
!	
!	print *, '----------------------------------------------------------'
!	print *, ' First Derivatives'
!	print *, '----------------------------------------------------------'

!	call GramSchmidt
!	
!		!Finding out if the wf's break signature after ortho
!		do i=1,nwt
!			print *, 'Testing iwa=', i
!			call SpinorBreakSignature(HFBasis(i)%GetValue())
!		enddo
!	
!	
!	call DeriveAll()
!	
!	!Find out the actual signature quantum number of a state
!	do iwa=1,nwt
!	
!	print *, ''
!	print *, 'Processing iwa ', iwa
!	
!	Value = HFBasis(iwa)%GetValue()
!	
!	!s = Value%Grid(1,ny,nz,2,1)/abs(Value%Grid(nx,ny,nz,1,1))	
!	
!	do i=1,3
!		Der(i) = HFBasis(iwa)%GetDer(i)
!		RDer(i) = Conj(Der(i))
!		RDer(i)%Grid(:,:,:,3:4,:) = - RDer(i)%Grid(:,:,:,3:4,:)
!	enddo
!	RDer(1)%Grid = - RDer(1)%Grid
!	RDer(2)%Grid = - RDer(2)%Grid
!	
!	
!	Lap = HFBasis(iwa)%GetLap()
!	
!	
!	!Testing if Laplacian is invariant under Signature
!	C = 0
!	
!	!RLap = MultiplyI(Lap)
!	RLap = Conj(Lap)
!	RLap%Grid(:,:,:,3:4,:) = - RLap%Grid(:,:,:,3:4,:)

!!	do i=1,nx
!!		print*, Lap%Grid(i,1,1,1,1), RLap%Grid(nx-i+1,1,1,1,1)
!!		print*, Lap%Grid(i,1,1,2,1), RLap%Grid(nx-i+1,1,1,2,1)
!!		print*, Lap%Grid(i,1,1,3,1), RLap%Grid(nx-i+1,1,1,3,1)
!!		print*, Lap%Grid(i,1,1,4,1), RLap%Grid(nx-i+1,1,1,4,1)
!!		print *
!!	
!!	enddo
!	
!	do k=1,nz
!		do j=1,ny
!			do i=1,nx
!				if (abs((abs(Lap%Grid(i,j,k,1,1)) - abs(RLap%Grid(nx-i+1,j,k,1,1)))).ge. 1d-16 ) &
!				& print *, 'Lap1', i,j,k ,Lap%Grid(i,j,k,1,1), RLap%Grid(nx-i+1,j,k,1,1)
!				if (abs((abs(Lap%Grid(i,j,k,2,1)) - abs(RLap%Grid(nx-i+1,j,k,2,1)))).ge. 1d-16) &
!				& print *, 'Lap2', i,j,k,Lap%Grid(i,j,k,2,1), RLap%Grid(nx-i+1,j,k,2,1)
!				if (abs((abs(Lap%Grid(i,j,k,3,1)) - abs(RLap%Grid(nx-i+1,j,k,3,1)))).ge. 1d-16 ) &
!				& print *, 'Lap3', i,j,k,Lap%Grid(i,j,k,3,1), RLap%Grid(nx-i+1,j,k,3,1)
!				if (abs((abs(Lap%Grid(i,j,k,4,1)) - abs(RLap%Grid(nx-i+1,j,k,4,1)))).ge. 1d-16 ) &
!				& print *, 'Lap4', i,j,k,Lap%Grid(i,j,k,4,1), RLap%Grid(nx-i+1,j,k,4,1)
!			enddo
!		enddo
!	enddo
!	do l=1,3
!	do k=1,nz
!		do j=1,ny
!			do i=1,nx
!				if (abs((abs(Der(l)%Grid(i,j,k,1,1)) - abs(RDer(l)%Grid(nx-i+1,j,k,1,1)))).ge. 1d-16 ) &
!				& print *, 'Der1', i,j,k,l ,Der(l)%Grid(i,j,k,1,1), RDer(l)%Grid(nx-i+1,j,k,1,1)
!				if (abs((abs(Der(l)%Grid(i,j,k,2,1)) - abs(RDer(l)%Grid(nx-i+1,j,k,2,1)))).ge. 1d-16 ) &
!				& print *, 'Der2', i,j,k,l,Der(l)%Grid(i,j,k,2,1), RDer(l)%Grid(nx-i+1,j,k,2,1)
!				if (abs((abs(Der(l)%Grid(i,j,k,3,1)) - abs(RDer(l)%Grid(nx-i+1,j,k,3,1)))).ge. 1d-16 ) &
!				& print *, 'Der3', i,j,k,l,Der(l)%Grid(i,j,k,3,1), RDer(l)%Grid(nx-i+1,j,k,3,1)
!				if (abs((abs(Der(l)%Grid(i,j,k,4,1)) - abs(RDer(l)%Grid(nx-i+1,j,k,4,1)))).ge. 1d-16 ) &
!				& print *, 'Der4', i,j,k,l,Der(l)%Grid(i,j,k,4,1), RDer(l)%Grid(nx-i+1,j,k,4,1)
!			enddo
!		enddo
!	enddo
!	enddo
!	enddo
!	
!	print *, '----------------------------------------------------------'
!	print *, ' Then parts of single particle hamiltonian'
!	print *, '----------------------------------------------------------'
!	
!	call Densit(0)
!	call CalculateAllMoments()
!	
!	call ConstructPotentials()
!	
!	Psi = HFbasis(20)
!	
!	B = ActionOfB(Psi)                         
!      U = ActionOfU(Psi) 
!      W = ActionOfW(Psi)
!      
!      print *, 'Testing B'
!      call SpinorBReakSignature(B)
!      print *, 'Testing U'
!      call SpinorBReakSignature(U)      
!      print *, 'Testing W'
!	call SpinorBreakSignature(W)
!      
!	
!end subroutine TestSignatureBreaking


!subroutine SpinorBreakSignature(Psi)
!	type(Spinor), intent(in) :: Psi
!	integer :: i,j,k,l
!	
!	!Checks if Psi breaks signature and, if so, where
!	do l=1,4
!	do k=1,nz
!		do j=1,ny
!			do i=1,nx
!				if(abs((abs(Psi%Grid(i,j,k,l,1)) - abs(Psi%Grid(nx-i+1,j,k,l,1)))/(abs(Psi%Grid(i,j,k,l,1)))).ge.1d-14) then
!					print *, 'Problem at', i,j,k,l
!					print *,Psi%Grid(i,j,k,l,1), Psi%Grid(nx-i+1,j,k,l,1), &
!                     & abs((abs(Psi%Grid(i,j,k,l,1)) - abs(Psi%Grid(nx-i+1,j,k,l,1)))/(abs(Psi%Grid(i,j,k,l,1))))
!				endif
!			enddo
!		enddo
!	enddo
!	enddo
!end subroutine SpinorBreakSignature


! subroutine TestTimeReversal
!     !----------------------------------------------------------------------------------------------------
!     ! This subroutine tests which parts of the single particle Hamiltonian break Time Reversal Symmetry
!     !   
!     !----------------------------------------------------------------------------------------------------    

!     use Derivatives
!     use GenInfo
!     use Densities
!     use MeanFields
!     use SpwfStorage

!     implicit none

!     integer      :: wave(2),i, iwa
!     type(Spwf)   :: wf1, wf2
!     type(Spinor) :: Spin(2), B(2), A(2),W(2), S(2), U(2), Der(2,3), Lap(2)
!     logical      :: Found=.false.
!    
!     print *, "----------------------------------------------------------"
!     print *, 'Testing Time Reversal Invariance'
!     print *, '----------------------------------------------------------'
!     

!     call DeriveAll()

!     call printSpwf()

!     !Calculating the densities
!     call Densit(0)

!     !Calculating the Potentials
!     call ConstructPotentials()

!     if (any (abs(Spot).gt.1.0d-14)) then
!      print *, 'SPot is not zero!'
!      print *, maxval(abs(Spot))
!     endif
!     
!     if (any (abs(Apot).gt.1.0d-14)) then
!      print *, 'APot is not zero!'
!      print *, maxval(abs(Apot))
!     endif
!     
!     if (any (abs(VecJ).gt.1.0d-14)) then
!      print *, 'vecj is not zero!'
!      print *, maxval(abs(VecJ))
!     endif
!     if (any (abs(RotVecJ).gt.1.0d-14)) then
!      print *, 'rotvecj is not zero!'
!      print *, maxval(abs(RotVecJ(:,:,:,1,:)))
!      print *, maxval(abs(RotVecJ(:,:,:,2,:)))
!      print *, maxval(abs(RotVecJ(:,:,:,3,:)))
!     endif
!     if (any (abs(Vecs).gt.1.0d-14)) then
!      print *, 'Vecs is not zero!'
!      print *, maxval(abs(vecs))
!     endif
!     if (any (abs(vect).gt.1.0d-14)) then
!       print *, 'VecT is not zero!'
!       print *, maxval(abs(vect))
!       print *, maxval(abs(vect(:,:,:,1,:)))
!       print *, maxval(abs(vect(:,:,:,2,:)))
!       print *, maxval(abs(vect(:,:,:,3,:)))

!     endif
!     if (any (abs(vecf).gt.1.0d-14)) then
!      print *, 'vecF is not zero!'
!      print *, maxval(abs(vecF))
!     endif
!     if (any (abs(laps).gt.1.0d-14)) then
!      print *, 'Laps is not zero!'
!      print *, maxval(abs(lapS))
!    endif
!     if (any (abs(graddivs).gt.1.0d-14)) then
!       print *, 'graddivs is not zero!'
!       print *, maxval(abs(graddivs))
!     endif
!     do iwa=1,nwt
!         print*, 'Testing WaveFunction ', iwa    
!         Found=.false.
!         wave(1)=iwa
!         wave(2)=-1


!         Spin(1) = HFBasis(wave(1))%GetValue()
!         do i=1,nwt
!             Spin(2) = HFBasis(i)%GetValue()
!             Spin(2) = TimeReverse(Spin(2))

!             if(any(abs(Spin(2)%Grid - Spin(1)%Grid).ge.prec13)) then
!                
!             else
!                  wave(2) = i
!                  print *, 'Found Time Reversed partner at ', i
!                  Found=.true.
!                  exit
!             endif
!         enddo
!         if(.not.Found) then
!                 print *, 'No Time-Reversed Partners'
!                 cycle
!         endif

!         !Checking if we are dealing with time reversed partners
!         do i=1,2
!            B(i) = ActionOfB(HFBasis(wave(i))) 
!            U(i) = ActionOfU(HFBasis(wave(i))) 
!            A(i) = ActionOfA(HFBasis(wave(i)))
!            S(i) = ActionOfS(HFBasis(wave(i))) 
!            W(i) = ActionOfW(HFBasis(wave(i)))
!         enddo

!         B(2) = TimeReverse(B(2))
!         U(2) = TimeReverse(U(2))
!         S(2) = TimeReverse(S(2))
!         A(2) = TimeReverse(A(2))
!         W(2) = TimeReverse(W(2))

!         !Testing the Actions
!         if(any(abs(B(1)%Grid - B(2)%Grid).ge.prec13 )) then
!                 print *, 'B is breaking time reversal!'
!         else
!                 print *, 'B does not break Time Reversal'
!         endif

!          if(any(abs(U(1)%Grid - U(2)%Grid).ge. prec13 )) then
!                 print *, 'U is breaking time reversal!'
!         else
!                 print *, 'U does not break time reversal'
!         endif

!          if(any(abs(S(1)%Grid - S(2)%Grid).ge. prec13 )) then
!                 print *, 'S is breaking time reversal!'
!         else
!                 print *, 'S does not break timereversal'
!         endif
!          if(any(abs(W(1)%Grid - W(2)%Grid).ge. prec13 )) then
!                 print *, 'W is breaking time reversal!'
!          else
!                 print *, 'W is not breaking time reversal'
!         endif

!         if(any(abs(A(1)%Grid - A(2)%Grid).ge. prec13 )) then
!                 print *, 'A is breaking time reversal!'
!         else
!                 print *, 'A is not breaking Time reversal.'
!         endif

!         !Testing the derivatives

!         Spin(2) =  HfBasis(wave(2))%GetValue()
!         do i=1,3
!             Der(1,i)  =  HfBasis(wave(1))%GetDer(i)
!             Der(2,i)  =  HfBasis(wave(2))%GetDer(i)
!             Der(2,i)  = TimeReverse(Der(2,i))
!             if(any(abs(Der(1,i)%Grid - Der(2,i)%Grid).ge. prec13 )) then
!                 print*, 'NonCorrect derivatives for direction ', i 
!             else
!                 print*, "Derivatives are time-reversal invariant in direction ", i
!             endif
!         enddo

!         Lap(2) =  HfBasis(wave(2))%GetValue()
!         Lap(1)  =  HfBasis(wave(1))%GetLap()
!         Lap(2)  =  HfBasis(wave(2))%GetLap()
!         Lap(2)  = TimeReverse(Lap(2))
!         if(any(abs(Lap(1)%Grid - Lap(2)%Grid).ge. prec13 )) then
!             print*, 'NonCorrect laplacians for direction ', i 
!         else
!             print*, "Laplacian is time-reversal invariant"
!         endif
!     enddo

!     print *, 'Evolving for one Iteration!'

!     call Evolve(1,1)

!     call printSpwf

!     if (any (abs(Spot).gt.1.0d-14)) then
!      print *, 'SPot is not zero!'
!      print *, maxval(abs(Spot))
!     endif
!     
!     if (any (abs(Apot).gt.1.0d-14)) then
!      print *, 'APot is not zero!'
!      print *, maxval(abs(Apot))
!     endif
!     
!     if (any (abs(VecJ).gt.1.0d-14)) then
!      print *, 'vecj is not zero!'
!      print *, maxval(abs(VecJ))
!     endif
!     if (any (abs(RotVecJ).gt.1.0d-14)) then
!      print *, 'rotvecj is not zero!'
!      print *, maxval(abs(RotVecJ(:,:,:,1,:)))
!      print *, maxval(abs(RotVecJ(:,:,:,2,:)))
!      print *, maxval(abs(RotVecJ(:,:,:,3,:)))
!     endif
!     if (any (abs(Vecs).gt.1.0d-14)) then
!      print *, 'Vecs is not zero!'
!      print *, maxval(abs(vecs))
!     endif
!     if (any (abs(vect).gt.1.0d-14)) then
!       print *, 'VecT is not zero!'
!       print *, maxval(abs(vect))
!       print *, maxval(abs(vect(:,:,:,1,:)))
!       print *, maxval(abs(vect(:,:,:,2,:)))
!       print *, maxval(abs(vect(:,:,:,3,:)))

!     endif
!     if (any (abs(vecf).gt.1.0d-14)) then
!      print *, 'vecF is not zero!'
!      print *, maxval(abs(vecF))
!     endif
!     if (any (abs(laps).gt.1.0d-14)) then
!      print *, 'Laps is not zero!'
!      print *, maxval(abs(lapS))
!    endif
!     if (any (abs(graddivs).gt.1.0d-14)) then
!       print *, 'graddivs is not zero!'
!       print *, maxval(abs(graddivs))
!     endif

!     do iwa=1,nwt
!         print*, 'Testing WaveFunction ', iwa    
!         Found=.false.
!         wave(1)=iwa
!         wave(2)=-1

!          Spin(1) = HFBasis(wave(1))%GetValue()

!         do i=1,nwt
!             Spin(2) = HFBasis(i)%GetValue()
!             Spin(2) = TimeReverse(Spin(2))

!             if(any(abs(Spin(2)%Grid - Spin(1)%Grid).ge.prec4)) then
!               
!             else
!                  wave(2) = i
!                  print *, 'Found Time Reversed partner at ', i
!                  Found=.true.
!                  exit
!             endif

!         enddo

!        
!             if(.not.Found) then
!                 print*, ('No Time-Reversed Partners')
!                 cycle
!             endif

!             !Checking if we are dealing with time reversed partners
!             do i=1,2
!                B(i) = ActionOfB(HFBasis(wave(i))) 
!                U(i) = ActionOfU(HFBasis(wave(i))) 
!                A(i) = ActionOfA(HFBasis(wave(i)))
!                S(i) = ActionOfS(HFBasis(wave(i))) 
!                W(i) = ActionOfW(HFBasis(wave(i)))
!             enddo

!             B(2) = TimeReverse(B(2))
!             U(2) = TimeReverse(U(2))
!             S(2) = TimeReverse(S(2))
!             A(2) = TimeReverse(A(2))
!             W(2) = TimeReverse(W(2))

!             !Testing the Actions
!             if(any(abs(B(1)%Grid - B(2)%Grid).ge.prec13 )) then
!                     print *, 'B is breaking time reversal!'
!             endif

!              if(any(abs(U(1)%Grid - U(2)%Grid).ge. prec13 )) then
!                     print *, 'U is breaking time reversal!'
!             endif

!              if(any(abs(S(1)%Grid - S(2)%Grid).ge. prec13 )) then
!                     print *, 'S is breaking time reversal!'
!             endif
!              if(any(abs(W(1)%Grid - W(2)%Grid).ge. prec13 )) then
!                     print *, 'W is breaking time reversal!'
!             endif

!             if(any(abs(A(1)%Grid - A(2)%Grid).ge. prec13 )) then
!                     print *, 'A is breaking time reversal!'
!             endif
!         enddo
!       

!    


!     print *, "----------------------------------------------------------"
!     print *, 'End of the test of Time Reversal Invariance'
!     print *, '----------------------------------------------------------'


!end subroutine TestTimeReversal


! subroutine TestParity
!     !-------------------------------   
!     ! This routine tests if MOCCa spontaneaously breaks parity when it should not.
!     !
!     !-------------------------------

!     use Derivatives
!     use GenInfo
!     use Densities
!     use MeanFields
!     use SpwfStorage

!     implicit none

!     integer      :: wave(2),i, iwa, j
!     type(Spwf)   :: wf1, wf2
!     type(Spinor) :: Spin, B, A,W, S, U, Der(3), Lap
!     logical      :: ParityBroken(3)=.false.
    
!     print *, "----------------------------------------------------------"
!     print *, 'Testing Parity Invariance'
!     print *, '----------------------------------------------------------'

!     !Deriving all wavefunctions
!     call DeriveAll()

!     do iwa=1,nwt
        
!         ParityBroken(:) = .false.

!         print *, 'Wavefunction nr. ', iwa

!         !Check if each wavefunction if it is still invariant

!         Spin = HfBasis(iwa)%GetValue()
        
!         do i=1,nz/2
!                 if(any(abs(abs(Spin%Grid(:,:,i,:,1)) - abs(Spin%Grid(:,:,nz-i+1,:,1))) .ge. Prec13 )) then
!                     ParityBroken(1)=.true.
!                 endif
!         enddo
            
!         if(ParityBroken(1)) then
!             print *, 'Parity broken on wf level!'
!             cycle !no need to compare derivatives
!         endif
        
!         !Check for each wavefunction if the derivatives are also parity invariant
!         do j=1,3
!             Der(j) = HFBasis(iwa)%GetDer(j)       
!              do i=1,nz/2
!                 if(any(abs(abs(Der(j)%Grid(:,:,i,:,1)) - abs(Der(j)%Grid(:,:,nz-i+1,:,1))) .ge. Prec13 )) then
!                     ParityBroken(j)=.true.
!                 endif
!              enddo
             
!             if(ParityBroken(j)) then
!                 print *, 'Parity broken for ', j , ' derivative.' 

!                 do i=1,nz/2

!                     print *, Der(j)%Grid(1,1,i,1,1), Der(j)%Grid(1,1,nz-i+1,1,1), Der(j)%Grid(1,1,nz-i+1,1,1)/Der(j)%Grid(1,1,i,1,1)
!                 enddo
!                 exit
!             endif
!         enddo

       




!     enddo


!     call PrintSpwf()


!     print *, "----------------------------------------------------------"
!     print *, 'End of Testing of Parity Invariance'
!     print *, '----------------------------------------------------------'


! end subroutine TestParity



    
 end module Testing
