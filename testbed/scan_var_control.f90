!
!=======================================================================
      subroutine Plot2D3D(Action)
!=======================================================================      
! Calculates X and Y variables for plotting
!
! Modifications:
!  Created by Kevin Schaefer (7/1/01)
!  Kevin Schaefer switched to arrays to store scanning info
!----------------------------------------------------------------------
!
      use scan_Variables
!
      implicit none
!
! input variables
      character*8 Action ! flag indicating action to take when executing
!
! read in variable data
      if(Action=='init    ') then
        call AssignVariableValues
      endif
!
! calculate Delta X
      if(Action=='DeltaX  ') then
          LabX=Label(Xvar)
          Dx=Range(Xvar)/real(numpts-1)
          Xmin=PlotMin(Xvar)
          Xmax=PlotMax(Xvar)
      endif
!
! calculate Delta Y
      if(Action=='DeltaY  ') then
          LabY=Label(Yvar)
          DY=Range(Yvar)/real(numpts-1)
          Ymin=PlotMin(Yvar)
          Ymax=PlotMax(Yvar)
      endif
!
! initialize data values
      if(Action=='StartX  ') then
!       reset scanned variable to start value
        Value(Xvar)=Start(Xvar)
!
!       reset all variables
        call AssignVariableValues
      endif
!
      if(Action=='StartY  ') then
!       reset scanned variable to start value
        Value(Yvar)=Start(Yvar)
!
!       reset all variables
        call AssignVariableValues
      endif
!
! assign x values to X variable
      if(Action=='AssXval ') Xval=Value(Xvar)
!
! increment X variable
      if(Action=='IncrX   ') then
        Value(Xvar)=Value(Xvar)+Dx
!
!       reset all variables
        call AssignVariableValues
      endif
!
! increment Y variable
      if(Action=='IncrY   ') then
        Value(Yvar)=Value(Yvar)+Dy
!
!       reset all variables
        call AssignVariableValues
      endif      
!
      return                                                                    
      end
