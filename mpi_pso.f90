!*******************************************************************************
! PROGRAM MPI_PSO
! PSO optimization of function using MPI
!
! Program reads user input for number of particles in swarm, then optimizes
! defined in code. Convergence criteria are hard coded.
! Func= (1-sqrt(x**2+y**2))+0.25*cos(x*3*pi)+0.25*cos(y*3*pi) has max of 1.5 at (0,0)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Author: S J Davie
!         stuart.davie@y7mail.com
!
! Written: 30JUN2016
!*******************************************************************************
program MPI_PSO

!-------------------Declare Variables-----------------------------
implicit none
include 'mpif.h'

!PSO Variables
integer							:: TPart,NPart,i,j,MaxIter,ConvIt,CurConvIt,CurTotIt
double precision, allocatable 	:: PartPos(:,:), PartVel(:,:),PartVal(:)
double precision				:: Func, Conv, GlobVal, GlobPos(2), ConvDiff
parameter (MaxIter=1000,Conv=0.001,ConvIt=50)

!MPI variables
integer 			:: ierr, MyId, NumProcs,Master
double precision		:: MaxVal_Loc(2), MaxVal_Glob(2)
parameter (Master = 0)
!-----------------------------------------------------------------

!----Initialize MPI Environment and obtain process ID's etc.------
call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world, MyId, ierr)
call mpi_comm_size(mpi_comm_world, NumProcs, ierr)
!-----------------------------------------------------------------

!--------------Obtain number of particles-------------------------
if(MyId.eq.Master) then
   print*,'No. of particles?'
   read*,TPart
endif
!-----------------------------------------------------------------

!-------Divide particles between processors-----------------------
call mpi_bcast(TPart,1,mpi_double_precision,Master,mpi_comm_world,ierr) 
if(MyId.le.mod(TPart,NumProcs)-1) then
   Npart = TPart/NumProcs+1
else
   Npart = TPart/NumProcs
end if
if(NPart.eq.0)Npart=1

print*,'Process',MyId,'has',Npart,'particles'
!-----------------------------------------------------------------

!-------------Allocate variables, seed rand-----------------------
allocate(PartPos(NPart,2))
allocate(PartVel(NPart,2))
allocate(PartVal(NPart))
CurConvIt=0;CurTotIt=0
call srand(time()*(10*MyId+1)*(10*MyId+1))
!-----------------------------------------------------------------

!-------------Initialize PSO variables----------------------------
Do i=1,Npart
   Do j=1,2
      PartPos(i,j)=2*rand()-1
      PartVel(i,j)=(2*rand()-1)/2
   enddo
   PartVal(i)=Func(PartPos(i,1),PartPos(i,2))
enddo

MaxVal_Loc(1)=MaxVal(PartVal,1) 									!Find initial, local max val
MaxVal_Loc(2)=MyId             		
GlobPos(:)=PartPos(MaxLoc(PartVal,1),:)									!Find intial, local max positions
call mpi_allreduce(MaxVal_Loc,MaxVal_Glob,1,mpi_2double_precision,mpi_maxloc,mpi_comm_world,ierr) 	!Find initial, global max val  
GlobVal= MaxVal_Glob(1)								
call mpi_bcast(GlobPos,2,mpi_double_precision,int(MaxVal_Glob(2)),mpi_comm_world,ierr)			!Share initial, global max pos.
!-----------------------------------------------------------------

!*******************Main PSO loop*********************************
Do while((CurTotIt.lt.MaxIter).and.(CurConvIt.lt.ConvIt))
   CurTotIt=CurTotIt+1											!Update max iterations counter

!-------------------update swarm positions and velocities---------
   Do i=1,NPart
      Do j=1,2
         PartVel(i,j)=0.95*PartVel(i,j)+0.1*(GlobPos(j)-PartPos(i,j))					!Update velocities
         PartPos(i,j)=PartPos(i,j)+PartVel(i,j)								!Update positions
      enddo
   PartVal(i)=Func(PartPos(i,1),PartPos(i,2))								!Update function val
   enddo
!-----------------------------------------------------------------

!-------------------update global maximum details-----------------
   MaxVal_Loc(1)=MaxVal(PartVal,1);MaxVal_Loc(2)=MyId
call mpi_allreduce(MaxVal_Loc,MaxVal_Glob,1,mpi_2double_precision,mpi_maxloc,mpi_comm_world,ierr)  	!Find max local val
   If(GlobVal .lt. MaxVal_Glob(1)) then							
      If(abs(GlobVal-MaxVal_Glob(1)).lt.Conv) then
         CurConvIt=CurConvIt+1										!
      else												!Update convergence counter
         CurConvIt=0											!
      endif
      GlobVal=MaxVal_Glob(1)										!Update global max val
      GlobPos(:)=PartPos(MaxLoc(PartVal,1),:)
call mpi_bcast(GlobPos,2,mpi_double_precision,int(MaxVal_Glob(2)),mpi_comm_world,ierr)			!Update global max pos
   else
      CurConvIt=CurConvIt+1
   endif
!-----------------------------------------------------------------

enddo
!*****************************************************************

!--------------------------Output---------------------------------
call sleep(1)	
if(MyId.eq.Master) then
print*,'Global maximum of',GlobVal,'at',GlobPos(:),'found after',CurTotIt,'iterations'
endif
!-----------------------------------------------------------------

call mpi_finalize(ierr)
end program


!---------------------Function------------------------------------
double precision Function Func(x,y)
double precision :: x,y,Val,pi
parameter (pi = 3.14159)
Func= (1-sqrt(x**2+y**2))+0.25*cos(x*3*pi)+0.25*cos(y*3*pi)
return
end function
!-----------------------------------------------------------------
    

