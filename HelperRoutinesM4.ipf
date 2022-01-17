#pragma rtGlobals=1		// Use modern global access method.

// 12/10/13  - added  "Pick Input Movie" .   Deleted obsolete mask-defining routines

//----------- PickInputMovie() -------------
// 12/10/13 - opens movie with raw data to analyze, and determines how many frames it has.

Function PickInputMovie(ctrlName) : ButtonControl   //--------Set Working Folder
	string ctrlName
NVAR InputMovieID_G, NmovieFrames_G

//   playmovie  //open input movie OBSOLETE IN OS11 
// 9/8/21 DHR added for Igor9
   String/G MoviePath
   MoviePath = DoOpenFileDialog()
   PlayMovieAction open = MoviePath
   
   PlayMovieAction getID
   inputmovieID_G = V_Value  // get  input movie id
   
   // get number of frames in input movie (from Igor examples, uses time per frame)
   PlayMovieAction gotoEnd,getTime
   Variable tend= V_Value
   PlayMovieAction step=-1,getTime
   NmovieFrames_G = floor(tend/(tend-V_value))
   
   PlayMovieAction gotoBeginning   // reset to first frame

end

//---------------DoOpenFileDialog() ----------
// 9/7/21 DHR -  From Igor examples.  
Function/S DoOpenFileDialog()
	Variable refNum
	String message = "Select a file"
	String outputPath
	//String fileFilters = "Data Files (*.txt,*.dat,*.csv):.txt,.dat,.csv;"	
	String fileFilters = "All Files:.*;"
	//fileFilters += "All Files:.*;"

	Open /D /R /F=fileFilters /M=message refNum
	outputPath = S_fileName
	
	return outputPath		// Will be empty if user canceled
End

///--------------------------
Function SetWorkingFolder(ctrlName) : ButtonControl   //--------Set Working Folder
	string ctrlName
	//NewPath/O ToPractice
	NewPath/C/O MovieFolder
end

///---------------Helper function to Make it quiver
//returns 0 < theta < 2pi
function angle(x,y)
    variable x,y
     variable val,sgn
     if (x == 0)
         if (y >0)
             val = pi/2
         else
             val = -pi/2
         endif
     elseif (x >0)
         val = atan(y/x)
     else  // x < 0
         val = atan(y/x) + pi
     endif
return(val)
end
///--------------------------
///---------------------------------------------------------------------

Function CancelButton(ctrlName) : ButtonControl  //Done Button 
string ctrlName
NVAR WindowFlag, ClickVal
	if (WindowFlag==1)
		DoWindow/K WdWd
	elseif (WindowFlag==2)
     		DoWindow/K ITP
     		WindowFlag=1
      elseif (WindowFlag==3)
		DoWindow/K DefArr
		WindowFlag=1
	endif
	ClickVal=1
end
///---------------------------------------------------------------------

Function DoneButton(ctrlName) : ButtonControl  //Done Button 
string ctrlName
NVAR WindowFlag, ClickVal
	if (WindowFlag==1)
		DoWindow/K WdWd
	elseif (WindowFlag==2)
     		DoWindow/K ITP
     		WindowFlag=1
      elseif (WindowFlag==3)
		DoWindow/K DefArr
		WindowFlag=1
	endif
	ClickVal=0
end
///--------------------------------

function AcceptName(ctrlName):ButtonControl
string CtrlName
ControlUpdate RootStringValue
SVAR RootString,MaskString, PFS_String, LUTr_String, LUT_String,MaskString2
SVAR  FitBottom, FitErrBottom, FitRes_String, FitErr_String,  PFS_Arr
SVAR AttributeString, ShiftString, AttributeString2,bgshiftstring

PathInfo MovieFolder // Modified 06/08/06
NewPath/C/O ToPractice S_Path+RootString

//PFS_String="PFS_"+RootString
PFS_String="PSF_"+RootString
LUTr_String="LUTr_"+RootString
LUT_String="LUT_"+RootString
MaskString="Mask_"+RootString
MaskString2= MaskString +"2"
FitBottom="FB_"+RootString
FitErrBottom="FBErr_"+RootString
FitRes_String="FitRes_"+RootString
FitErr_String="FitErr_"+RootString
AttributeString="Att_"+RootString
AttributeString2="Att_"+RootString +"Mask2"
ShiftString="Shift"+RootString
//PFS_Arr="PFSArr_"+RootString
PFS_Arr="PSFArr_"+RootString
bgshiftstring="bgshift_"+RootString 
end


