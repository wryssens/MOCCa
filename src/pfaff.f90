module pfaff


  implicit none
  
contains 


!  11 May 2011
!
!+---------------------------------------------------------------------+
!|   Computes the Pfaffian of a skew-symmetric matrix SK using         |
!|   the ideas of Aitken's block diagonalization formula.              |
!|   Full pivoting is implemented to make the whole procedure stable   |
!+---------------------------------------------------------------------+
!|   Only the upper triangle of SK is used                             |
!|   In output the lower triangle contains the transformation matrix   | 
!+---------------------------------------------------------------------+
!|   SK .....Skew symmetric input matrix                               |

!|   LDS ....Leading dimension of matrix SK                            |

!|   N  .....Dimension of SK                                           |

!|   Ipiv ...Integer vector to hold the pivoted columns                |

!+---------------------------------------------------------------------+
!| This program may be freely used providing such use cites the source,|
!|      Numeric and symbolic evaluation of the pfaffian of general     |
!|      skew-symmetric matrices                                        |
!|                                                                     |
!|      C. Gonzalez-Ballestero, L.M.Robledo and G.F. Bertsch           |
!|      Computer Physics Communications (2011)                         |
!+---------------------------------------------------------------------+
!
       Subroutine PfaffianF(SK,LDS,N,Ipiv,Pf)
       Implicit none
       Double Precision SK(LDS,N)
       Double Precision Pf,big,SS,epsln,phas,SW
       Integer Ipiv(N,2),N,LDS,NB,IB,i,j,k,ip,NR,NC,I1,I2

       epsln = 1.0d-13   ! smallest number such that 1+epsln=1
       
       if(mod(N,2).eq.1) then ! N odd
       
       Pf = 0.0d+00 
       
       else
          
       NB= N / 2   ! Number of 2x2 blocks
              
       if (NB.eq.1) then 

! The pfaffian of a 2x2 matrix is Pf=Sk(1,2)
       
       Pf = SK(1,2)
       
       else   ! NB.gt.1
       
       Pf = 1.0d+00
       
       do IB=1,NB-1
          
          NR = IB*2 - 1 ! row numb of the 1,2 element of the 2x2 block
          NC = NR + 1

          big = 0.0d+00

          I1 = NR
          I2 = NC
                    
          do i=NR,N-1 ! all rows
             do j=i+1,N
          
                if(dabs(SK( i , j )).gt. big) then
                   big = dabs(SK( i , j ))
                   I1 = i
                   I2 = j 
                end if
                
             end do
          end do
          
          Ipiv( IB , 1) = I1 ! to initialize
          Ipiv( IB , 2) = I2 ! to initialize

          phas = 1.0d+00
          if(I1.eq.NR) phas = -phas
          if(I2.eq.NC) phas = -phas 
          
! Pivoting of element NR,NC with ip1,ip2

! 
!   Pivoting for a skew-symmetric matrix (Upper)
!
          if(I1.ne.NR) Call exch(SK,LDS,N,NR,I1)
          if(I2.ne.NC) Call exch(SK,LDS,N,NC,I2)
          
       ss = Sk( NR , NC )

       if(dabs(ss).gt.epsln) then
!
!      Updating the Schur complement  matrix
!

       do i = 2*IB+1 , N-1
       
          do j = i+1 , N

             Sk(i,j) = Sk(i,j) + &
     &                (Sk(NC,i)*Sk(NR,j)-Sk(NR,i)*Sk(NC,j))/SS
          end do
       end do
!
! Storing X and Y vectors in the lower part of the matrix
!       
       do i = 2*IB+1 , N
       
          Sk( i , NR ) =  -Sk( NC , i )/SS
          Sk( i , NC ) =   Sk( NR , i )/SS
          
       end do
        
       if (IB.ge.2) then  ! swap

!
!   Swapping
!       
       do j= 1 , 2*IB
       
          SW           = Sk( NR , j )
          Sk( NR , j ) = Sk( I1 , j ) 
          Sk( I1 , j ) = SW
       
          SW           = Sk( NC , j )
          Sk( NC , j ) = Sk( I2 , j ) 
          Sk( I2 , j ) = SW

       end do       
       
       end if
                      
       else
       
       Pf = 0.0d+00
       
       return
       
       end if ! dabs(ss).gt.epsln
       
       Pf = Pf * SS * phas
       
       end do ! IB
       
       Pf = Pf * Sk(N-1,N)
       
       end if !    NB.eq.1
          
       end if !    mod(N,2).eq.1
       
       return
       end
!+---------------------------------------------------------------------+
!  Exchange of rows i1 i2 and columns i1 i2
!  SK is assumed upper triangular
!  i1 < i2
!+---------------------------------------------------------------------+
       Subroutine Exch(SK,LDS,N,I1,I2)
       Implicit none
       Double Precision SK(LDS,N)
       Double Precision SS
       Integer LDS,N,I1,I2,I
              
       SK(i1,i2) = -SK(i1,i2)
                              
       if(I1.ne.1) then
          
          do i=1,I1-1
             SS           = SK( i , I1 )
             SK( i , I1 ) = SK( i , I2 )
             SK( i , I2 ) = SS
          end do
       end if

       if(I2.ne.N) then
          
          do i=I2+1,N
             SS           = SK( I1 , i )
             SK( I1 , i ) = SK( I2 , i )
             SK( I2 , i ) = SS
          end do
       end if

       if(I2.ge.I1+2) then
          
          do i=I1+1,I2-1
             SS           =  SK( I1 , i )
             SK( I1 , i ) = -SK( i , I2 )
             SK( i , I2 ) = -SS
          end do
       end if
       
       return
       end
!c+--------------------------------------------------------------------+
!c|   C o m p l e x    V e r s i o n                                   |
!c|                                                                    |
!c+--------------------------------------------------------------------+
       Subroutine ZPfaffianF(SK,LDS,N,Ipiv,Pf)
       Implicit none
       Double Complex SK(LDS,N)
       Double Complex Pf,SS,SW,one,zero
       Double Precision epsln,phas,big
       Integer Ipiv(N,2),N,LDS,NB,IB,i,j,k,ip,NR,NC,I1,I2

       one = dcmplx(1.0d+00,0.0d+00)
       zero= dcmplx(0.0d+00,0.0d+00)
       
       epsln = 1.0d-13   ! smallest number such that 1+epsln=1
       
       if(mod(N,2).eq.1) then ! N odd
       
       Pf = zero
       
       else
          
       NB= N / 2   ! Number of 2x2 blocks
              
       if (NB.eq.1) then 

! The pfaffian of a 2x2 matrix is Pf=Sk(1,2)
       
       Pf = SK(1,2)
       
       else   ! NB.gt.1
       
       Pf = one
       
       do IB=1,NB-1
          
          NR = IB*2 - 1 ! row numb of the 1,2 element of the 2x2 block
          NC = NR + 1

          big = 0.0d+00

          I1 = NR
          I2 = NC
                    
          do i=NR,N-1 ! all rows
             do j=i+1,N
          
                if(abs(SK( i , j )).gt. big) then
                   big = abs(SK( i , j ))
                   I1 = i
                   I2 = j 
                end if
                
             end do
          end do
          
          Ipiv( IB , 1) = I1 ! to initialize
          Ipiv( IB , 2) = I2 ! to initialize

          phas = 1.0d+00
          if(I1.eq.NR) phas = -phas
          if(I2.eq.NC) phas = -phas 
          
! Pivoting of element NR,NC with ip1,ip2

! 
!   Pivoting for a skew-symmetric matrix (Upper)
!
          if(I1.ne.NR) Call Zexch(SK,LDS,N,NR,I1)
          if(I2.ne.NC) Call Zexch(SK,LDS,N,NC,I2)
          
       ss = Sk( NR , NC )


       if(abs(ss).gt.epsln) then
!
!      Updating the Schur complement  matrix
!

       do i = 2*IB+1 , N-1
       
          do j = i+1 , N

             Sk(i,j) = Sk(i,j) + &
     &                (Sk(NC,i)*Sk(NR,j)-Sk(NR,i)*Sk(NC,j))/SS
          end do
       end do
!
! Storing X and Y vectors in the lower part of the matrix
!       
       do i = 2*IB+1 , N
       
          Sk( i , NR ) =  -Sk( NC , i )/SS
          Sk( i , NC ) =   Sk( NR , i )/SS
          
       end do
        
       if (IB.ge.2) then  ! swap

!
!   Swapping
!       
       do j= 1 , 2*IB
       
          SW           = Sk( NR , j )
          Sk( NR , j ) = Sk( I1 , j ) 
          Sk( I1 , j ) = SW
       
          SW           = Sk( NC , j )
          Sk( NC , j ) = Sk( I2 , j ) 
          Sk( I2 , j ) = SW

       end do       
       
       end if
                      
       else
       
       Pf = zero
       
       return
       
       end if ! dabs(ss).gt.epsln
       
       Pf = Pf * SS * phas
       
       end do ! IB
       
       Pf = Pf * Sk(N-1,N)
       
       end if !    NB.eq.1
          
       end if !    mod(N,2).eq.1
       
       return
       end
!
!  Exchange of rows i1 i2 and columns i1 i2
!  SK is assumed upper triangular
!  i1 < i2
!
       Subroutine ZExch(SK,LDS,N,I1,I2)
       Implicit none
       Double Complex SK(LDS,N)
       Double Complex SS
       Integer LDS,N,I1,I2,I
              
       SK(i1,i2) = -SK(i1,i2)
                              
       if(I1.ne.1) then
          
          do i=1,I1-1
             SS           = SK( i , I1 )
             SK( i , I1 ) = SK( i , I2 )
             SK( i , I2 ) = SS
          end do
       end if

       if(I2.ne.N) then
          
          do i=I2+1,N
             SS           = SK( I1 , i )
             SK( I1 , i ) = SK( I2 , i )
             SK( I2 , i ) = SS
          end do
       end if

       if(I2.ge.I1+2) then
          
          do i=I1+1,I2-1
             SS           =  SK( I1 , i )
             SK( I1 , i ) = -SK( i , I2 )
             SK( i , I2 ) = -SS
          end do
       end if
       
       return
       end
       
end module

