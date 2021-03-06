#pragma rtGlobals=1		// Use modern global access method.
#pragma rtGlobals=1		// Use modern global access method.

#include <Keyword-Value>
#include <Graph Utility Procs>
#include <Image Common>
#include  <All IP Procedures>

// *************MODIFICATION HISTORY  - EFFECTIVE MAY 19, 2006

//   Rules:   List any procedure files that are changed  - noting the date of the change
//                Details on the changes should be given in the procedure file itself.
//                The changes should be listed at the top of the file (date stamped) and also noted at the point of change in the code

//12/7/21 modified to open/close movie

//12/5/21  DHR In Initialize and Pick ROI GUI:   reordered buttons to fit process flow. Removed Total number of frames variable
// as not changed ever.   Removed "change intensities" button and controls.  Set Low bound and High bound to min/max intensity values of image.
// In DefineMaskGUI:  got rid of "2 mask" buttons.   Got rid of force vector plotting. Rearranged buttons for adding cells.

//12/3/21  initial guess stuff removed
// in DefineArray2M() get rid of "choice" of grid method.
// 09/13/21 DHR: modified ROI Select - to make wave references work

// 12/10/13 - Modified to take input from movie file.  Movie previously "opened" by PickInputMovie()

//1/18/2012 Craig added fourth point in ROI select resulting in better initial guesses for post positions

// 12/11/08   GUI MASKMAKER
// Maskwave no longer made in ROI Select - now in GUI maskmaker 

// 12/5/08  Lower left point in define array
//     modern clicking code in 

// 8/28/08   First attempt to "modernized "Clicking Code"  

//static strconstant  ksRhodFileRoot = "RHOD"  // for 2-digit names
static strconstant  ksRhodFileRoot = "RHOD"  // for 3-digit names

///---------------------This is envelope function which pops up a headline menu-------------------
///------------------------------------------------------------------------------------------------------------------------
Function DefineArray2M()
	Make/O/N=1 WMPositionMeasurementX,WMPositionMeasurementY
	NVAR WindowFlag, ClickVal; WindowFlag=3

	Variable/G LPX, LPY, RPX, RPY
	Variable/G  LLPX, LLPY,LRPX,LRPY  // 12/5/08  lower left points ,1/18/12 lower right points
	String/G MaskString,  AttributeString, ShiftString

	WMInitCalibrations()
	DoWindow/K DefArr
	NewPanel/K=1/N=DefArr/W=(200,400,600,570) //as "Now Define The Array?"
	SetDrawLayer UserBack
	ControlInfo GridMethod   //This is where the choice is made on the grid control method  12/3/21 not used anywhere else

		DrawText 1,15,"Grid Method "
		Button ROIselect,pos={75,23},size={225,20},proc=ROIselectM,title="Initialize and Pick ROI"
//		Button InitShifts,pos={75,63},size={225,20},proc=InitCalibM,title="Initialize shifts "
		Button MaskGUI,pos={75,53},size={225,20},proc=MakeMaskGUIM,title="DefineMaskGUI"
		Button DoneButton,pos={75,83},size={225,20},proc=DoneButton,title="\\K(65280,0,0)\f01Done"
		Button CancelButton,pos={75,113},size={225,25},proc=CancelButton,title="\\K(0,65280,0)\f01Cancel"
	

End
///-------------------Below are the functions for each button------------------
///--------------------------------------------------------------------------------------------
///--------------------------------------------------------------------------------------------

Function IntensitiesButton(ctrlName): ButtonControl
string ctrlName
NVAR  IntenseLowVar, IntenseHighVar
ControlUpdate IntenseLow
ControlUpdate IntenseHighVar
ModifyImage Step1Image ctab= {IntenseLowVar, IntenseHighVar,Grays,0}
doupdate

end

///----------------ROIselectM() ------------
//
//12/7/21 modified to open/close movie
// 9/13/21 DHR - fixed wave references
// 12/10/13 - movie input modification
//  8/28/08  Modified for modern Igor clicking

Function ROIselectM(ctrlName) : ButtonControl              ///-------------Selecting ROI, saves a whole bunch of data
	string ctrlName								///--------------key program!!!----------------------

// new 12/10/13 for movie input
   NVAR InputMovieID_G, NmovieFrames_G   //12/10/13
   SVAR MoviePath          // added 12/7/21 

	Variable/G 	TNFrames=NmovieFrames_G  //12/10/13
	Variable/G PostsVert, PostsHoriz, Bangle=0, doneflag, magnification = 60, hexcp = 1
	Variable/G IntenseLowVar, IntenseHighVar


	SVAR RootString, PWF
	NVAR LPX, LPY, RPX, RPY, ClickVal
	NVAR  LLPX, LLPY,LRPX,LRPY  // 12/5/08  lower left points, 1/18/12 lower right points
	SVAR MaskString, ShiftString, AttributeString
	ShiftString="Shift"+RootString
	MaskString="Mask_"+RootString
	AttributeString="Att_"+RootString
	
// new 12/10/13 for movie input
//9/13/21   Fixed for Igor 9 - have to have wave assignments after call that creates them
   PlayMovieAction Open = MoviePath  //12/7/21

   PlayMovieAction  gotoBeginning                 // 9/13/21  remove stop Separate commands
   PlayMovieAction  extract                 // extract one frame
   PlayMovieAction kill    // 12/7/21 close movie
   WAVE M_MovieFrame					//Identify wave after it is created
   ImageTransform rgb2gray M_MovieFrame   // convert to grayscale in M_rgb2gray
   WAVE M_rgb2gray
   Duplicate/O M_rgb2gray Step1Image

//12/5/21  to replace ChangeIntensities
   wavestats/Q Step1Image
   IntenseLowVar = V_min
   IntenseHighVar = V_max
       
	NewImage/N=Step1 Step1Image
	String dfSav= GetDataFolder(1)
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S root:Packages:WMImProcess
	NewDataFolder/O/S root:Packages:WMImProcess:ImageRotate

	Variable/G angle,fill,inited
	variable/G xxx
	if( inited!=1 )
		angle=10
		fill=0
	endif

	SetDataFolder dfSav
	
	variable doneflag2
	
	prompt doneflag2  "Set flag = 0 for a new analysis or to reset values.  Set = 1 to use previously input values."
	doprompt "Define ROI and post grid:" doneflag2
	wave tosave
	
	if(doneflag2==0)
	  //if(doneflag==0)
		 
		Make/O/N=4 XT,YT  // 8/28/08  temp waves for leftpoint (0) and rightpoint (1) 12/5/08 lowerleftpoint(2)
		XT[0] = 10 // initial values
		YT[0] = 10
		XT[1] = 500
		YT[1] = 10
		XT[2] = 10
		YT[2] = 250
		XT[3] = 500
		YT[3] = 250
	
	 elseif(doneflag2==1)
	       Make/O/N=4 XT,YT
	       wave XT2,YT2
		XT[0] = XT2[0]
		YT[0] = YT2[0]
		XT[1] = XT2[1]
		YT[1] = YT2[1]
		XT[2] = XT2[2]
		YT[2] = YT2[2]
		XT[3] = XT2[3]
		YT[3] = YT2[3]
	// endif
	endif
   // overplotleft and right point positions
             Appendtograph/T YT vs XT
             ModifyGraph mode(YT)=3, rgb(YT)=(65535,0,0)   //red

	Variable/G ActiveTempPoint1 
	ActiveTempPoint1 = -1
		
	NVAR WindowFlag
	NVAR Bangle= root:Packages:WMImProcess:ImageRotate:angle
	
	Bangle=0

	WindowFlag=2
	DoWindow/K ITP
	NewPanel/K=1/N=ITP/W=(200,400,650,800) as "Input Parameters"
	///--------------
			
	if(WMIsCalibrationApplied()==0)
	WMInitCalibrations()
	endif
	DoWindow/F WMMarkPositionsPanel
	if(V_Flag==1)
		AutoPositionWindow /E/R=$WMTopImageGraph()
		return 0
	endif
		
	String topImage=WMTopImageGraph()
	if(strlen(topImage)>0)
		AutoPositionWindow /E/R=$topImage
	endif
	SetDrawEnv linefgc= (43520,43520,0),fillpat= 0;	DrawRect 10,20,430, 70
	SetDrawEnv linefgc= (43520,43520,0),fillpat= 0;	DrawRect 10,100,430,240
//	SetDrawEnv linefgc= (43520,43520,0),fillpat= 0;	DrawRect 10,250,430,260	
	SetDrawEnv linefgc= (43520,43520,0),fillpat= 0;	DrawRect 10,270,430,330
	DrawText 15,154,"Upper Left:"
	DrawText 15,179,"Upper Right:"
	DrawText 15,204,"Lower Left:"
	DrawText 15,229,"Lower Right:"
	
	SetDrawLayer UserBack

	SetVariable angle,pos={16,33},size={130,15},title="Angle (degrees)",limits={-360,360,1},value= root:Packages:WMImProcess:ImageRotate:angle
	Button rotate,pos={160,30},size={109,20},proc=WMRotateButtonProc,title="Rotate"
	Button buttonRestore,pos={290,30},size={80,20},proc=WMImageBPUndo,title="Revert"

   Button MoveLeftRightPoints,pos={120,110},size={180,20},proc=MoveLRPointsDraw,title="Draw Rectangular ROI"   // 9/14/21 DHR changed from "Move L and R points"
      SetVariable LXtemp, pos ={95,140}, size = {150,18}, limits={0,Inf,1},value=XT[0]
      SetVariable LYtemp, pos ={265,140}, size = {150,18}, limits={0,Inf,1},value=YT[0]
      SetVariable RXtemp, pos ={95,165}, size = {150,18}, limits={0,Inf,1},value=XT[1]
      SetVariable RYtemp, pos ={265,165}, size = {150,18}, limits={0,Inf,1},value=YT[1]
      SetVariable LLXtemp, pos ={95,190}, size = {150,18}, limits={0,Inf,1},value=XT[2]   //12/5/08
      SetVariable LLYtemp, pos ={265,190}, size = {150,18}, limits={0,Inf,1},value=YT[2]
      SetVariable LRXtemp, pos={95,215},size = {150,18}, limits={0,Inf,1},value=XT[3]
      SetVariable LRYtemp, pos={265,215},size = {150,18}, limits={0,Inf,1},value=YT[3]
     
	
// 12/5/21	SetVariable ITP1,format="%d",pos={80,280},size={225,18},title="Total Number of Frames", limits={1,Inf,1},value= TNFrames
//***********
	SetVariable ITP2,pos={226,290},size={194,18},title="# posts in vertical dimension", format="%d",limits={1,Inf,1},value=PostsVert
	SetVariable ITP3,pos={16,290},size={194,18},title="# posts in horizontal dimension", format="%d", limits={1,Inf,1},value=PostsHoriz

	//SetVariable ITP2,pos={226,48},size={194,18},title="Enter 1 for Hex CP Array", format="%d",limits={0,Inf,1},value=hexcp
	//SetVariable ITP3,pos={16,48},size={194,18},title="Which objective did you use?", format="%d", limits={1,Inf,1},value=magnification
	hexcp = 1   // Flag indicating hexagonal array of posts
//***********

	
	ModifyPanel fixedSize=1
	Button DoneButton,pos={120,360},size={180,20},proc=DoneButton,title="\\K(65280,0,0)Done"
	PauseForUser ITP, Step1   // 9/14/21 DHR added Step1 here to make both windows accessible with mouse
	DoWindow/K Step1
	WindowFlag=3
	if (ClickVal==0)
	      LPX = XT[0]   // 8/28/08
	      LPY = YT[0]
	      RPX = XT[1]
	      RPY = YT[1]
	      LLPX = XT[2]
	      LLPY = YT[2]
	      LRPX = XT[3]
	      LRPY = YT[3]
		Make/O/N=(TNFrames, 14)  ToWorkWith=0   //12/3/21  this will be the old ShiftWave - maybe can be deleted?
		ToWorkWith[*][7]=Bangle
//		variable latticespacingpix
//		if(hexcp == 1)
//		
//			if(magnification == 60)
//				latticespacingpix = 37.2
//			elseif(magnification == 40)
//				latticespacingpix = 24
//			elseif(magnification == 20)
//				latticespacingpix = 12
//			endif
//			PostsHoriz = ceil(abs(RPX-LPX)/latticespacingpix +1 )//testing so I don't have to count CMK 2/12/09 37 is horiz lattice spacing in pixels for 60 x only!
//			PostsVert = round(abs(LLPY-LPY)/(latticespacingpix*sqrt(3)/2)+1)
//		else
//			if(magnification == 60)
//				latticespacingpix = 81
//			elseif(magnification == 40)
//				latticespacingpix = 54
//			elseif(magnification == 20)
//				latticespacingpix = 27
//			endif
//		PostsHoriz = round (abs(RPX-LPX)/latticespacingpix +1) //testing so I don't have to count CMK 2/12/09 37 is horiz lattice spacing in pixels for 60 x only!
//		PostsVert = round(abs(LLPY-LPY)/latticespacingpix+1)
//		endif
		//*************
		variable avgxdist
		variable avgydist
		avgxdist = (abs(RPX-LPX)+abs(LRPX-LLPX))/2
		avgydist = (abs(LLPY-LPY)+abs(LRPY-RPY))/2
		ToWorkWith[*][8]= avgxdist/(PostsHoriz-1) //kia = 8
		ToWorkWith[*][9] = avgydist/(PostsVert-1) // kia=9   12/5/08  clicking to get lattice-B
//12/3/21 added when removing initialize shifts button
      ToWorkWith[*][10] = 0
      ToWorkWith[*][11] = 0
		Duplicate/O ToWorkWith $ShiftString

		Save/O/C/P=ToPractice $ShiftString as (ShiftString+".ibw")  
//		Save/O/C/P=ToPractice $MaskString as (MaskString+".ibw")   // 12/11/08
     		Make/O/N=14 ToSave
     		ToSave[0]=LPX;      ToSave[1]=LPY;  //Coordinates of the left top post in the RHODBASE
     		ToSave[2]= Bangle;  // Bangle has been set before
     		ToSave[3]= PostsHoriz;      ToSave[4]= PostsVert; // Number of posts in horizontal and vertical directions
     		ToSave[5]= abs(RPX-LPX)/(PostsHoriz-1);  //Lattice constant
     		ToSave[6]=IntenseLowVar
     		ToSave[7]=IntenseHighVar
     		ToSave[8]=RPX; ToSave[9]= RPY
     		ToSave[10] = LLPX; ToSave[11] = LLPY
     		ToSave[11] = LRPX; ToSave[13] = LRPY
     		duplicate/o xt xt2
     		duplicate/o yt yt2
     		Duplicate/O ToSave $AttributeString 
     		Save/O/C/P=ToPractice $AttributeString as (AttributeString+".ibw")
     		Save/G/O/P=ToPractice ToSave as (AttributeString+".txt")
     		    		
      	endif

doneflag = 1

end
///------------------------------
Function ActivateGraph(ctrlName):ButtonControl
	string ctrlName
	DoWindow Step1
	if (V_Flag==1)
		DoWindow/F Step1
	endif
	DoWindow GraphProca
	if (V_Flag==1)
		DoWindow/F GraphProca
	endif
end
///---------------------------
Function ActivatePanel(ctrlName):ButtonControl
	string ctrlName
	Print "Activate Panel"
	DoWindow/F WMMarkPositionsPanel
end



////------------------------------InitCalib()
//12/3/21 - removed from DefArr
// 12/5/08  Modified for modern clicking
//12/10/13 - for Movie input - gutted so that shifts are all set to Zero.

Function InitCalibM(ctrlName) : ButtonControl         ///To imput zero displacements
string ctrlName									///Right button doesn't do anything at this point

NVAR NmovieFrames_G   //12/10/13
Variable f1=0, f2=NmovieFrames_G   //12/10/13
Prompt f1, "Enter First File Index"
Prompt f2, "Enter Last File Index"
//wave WMPositionMeasurementX, WMPositionMeasurementY
Make/O/N=0 WMPositionMeasurementX,WMPositionMeasurementY
NVAR LPX, LPY, RPX, RPY
SVAR RootString, ShiftString
string AttTemp="Att_"+RootString
ShiftString="Shift"+RootString

LoadWave/Q/O/P=ToPractice (ShiftString+".ibw")
LoadWave/G/Q/N='Abts'/O/P=ToPractice (AttTemp+".txt")
wave Abts0
Variable RAngle=Abts0[2]

Duplicate/O $ShiftString SSInitCalib
variable ifile, shX, shY, X0=0, Y0=0
string s1

DoPrompt "Indexes", f1, f2
if (V_flag==1)	
	return 0
endif

 
for (ifile=f1; ifile <=f2; ifile +=1)   // loop over raw data files
     SSInitCalib[ifile][10]=0//ShX
     SSInitCalib[ifile][11]=0//ShY
endfor 

Duplicate/O SSInitCalib $ShiftString
Save/O/C/P=ToPractice  $ShiftString as (ShiftString+".ibw")  // The path is temporary
return 0
end 

///------------------------------
Function LeftButtonProc(ctrlName) : ButtonControl
string ctrlName
NVAR LPX, LPY
wave WMPositionMeasurementX, WMPositionMeasurementY
	CustomAddPositionButtonProc(ctrlName)
	variable WMPX=numpnts(WMPositionMeasurementX), WMPY=numpnts(WMPositionMeasurementY)
	LPX=WMPositionMeasurementX[WMPX]; //Print "LPX", LPX
	LPY=WMPositionMeasurementY[WMPY]
	SetAxis/A/R left
	SetAxis/A top
end
////------------------------------------------------
Function RightButtonProc(ctrlName) : ButtonControl
string ctrlName
NVAR RPX, RPY
wave WMPositionMeasurementX, WMPositionMeasurementY
	CustomAddPositionButtonProc(ctrlName)
	variable WMPX=numpnts(WMPositionMeasurementX), WMPY=numpnts(WMPositionMeasurementY)
	RPX=WMPositionMeasurementX[WMPX];// Print "RPX", RPX
	RPY=WMPositionMeasurementY[WMPY]
	SetAxis/A/R left
	SetAxis/A top
end
///----------------------------------------------------------
Function ErasePointsHere(ctrlName) : ButtonControl
string ctrlName
ErasePoints(ctrlName)
Make/O/N=1 WMPositionMeasurementX,WMPositionMeasurementY
end

////-------------------------------------------------------------------------------------------------------------------
///----------------------------The part below refers to the graphic readout of cursors---------------
///------------------------------------------------------------------------------------------------------------------------

function proca(ProcaImage, zmin, zmax) // proca is the envelope function, all it does is it calls up the graph and also PositionButtonSasha()
wave ProcaImage
variable zmin, zmax
Duplicate/O ProcaImage OrigImage1
	DoWindow/K GraphProca
	NewImage/N=GraphProca Origimage1
	ModifyImage OrigImage1 ctab= {zmin,zmax,Grays,0}
	Button activatePanel, pos={100,55}, size={225,20}, proc=ActivatePanel, title="Activate Panel"
	ShowInfo
	//KillWaves removed on March 03
	//KillWaves/Z WMPositionMeasurementX,WMPositionMeasurementY
	PositionButtonSasha()
	
end
///---------------------- This one has both left and right button
function procaLeftRight(ProcaImage, zmin, zmax)
wave ProcaImage
variable zmin, zmax
NVAR LPX, LPY, RPX, RPY
Duplicate/O ProcaImage OrigImage
	DoWindow/K GraphProca
	NewImage/N=GraphProca Origimage
	ModifyImage Origimage ctab= {zmin,zmax,Grays,0}
	ShowInfo
	NewPanel/K=1/N=ITP/W=(200,200,650,400) as "Input Parameters"
	
	///--------------
	 		
	if(WMIsCalibrationApplied()==0)
	WMInitCalibrations()
	endif
	DoWindow/F WMMarkPositionsPanel
	//DoWindow/F ITP
	if(V_Flag==1)
		AutoPositionWindow /E/R=$WMTopImageGraph()
		return 0
	endif
		
	String topImage=WMTopImageGraph()
	if(strlen(topImage)>0)
		AutoPositionWindow /E/R=$topImage
	endif
		
	SetDrawLayer UserBack
	DrawText 13,18,"Move the red cursor to the location that you want to mark"
		
	ValDisplay curXValdisp,pos={16,20},size={150,18},title="X:"
	ValDisplay curXValdisp,limits={0,0,0},barmisc={0,1000}
	ValDisplay curXValdisp,value= #"root:Packages:WMCalibrations:curX"
      Button LeftPositionButton,pos={16,50},size={150,20},proc=LeftButtonProc,title="Left Point"
      Button RightPositionButton,pos={256,50},size={150,20},proc=RightButtonProc,title="Right Point"
	CheckBox positionTagCheck,disable=1, pos={74,213},size={120,20},title="Add with a tag", value=1
	CheckBox showLengthCheck, disable=1, pos={73,235},size={160,20},title="Show position in tag", value=1
	ValDisplay curYValdisp,pos={233,20},size={150,18},title="Y:"
	ValDisplay curYValdisp,format="%g",limits={0,0,0},barmisc={0,1000}
	ValDisplay curYValdisp,value= #"root:Packages:WMCalibrations:curY"
	Button erasemeasurementTablebutton,pos={100,80},size={225,20},proc=ErasePointsHere,title="Erase measurements"
	Button erasemeasurementTablebutton,help={"Removes traces from the image and kills the associated waves."}
	Button markPositionDoneButton,pos={213,110},size={180,20},proc=CustomMarkPositionDoneButton,title="Accept This Frame"

	String cdf=GetDataFolder(1)
	String newDF=CustImageCursors(1,3)	
	SetDataFolder  newDF
	WM_DrawCursors(1)
	SetDataFolder cdf
	///--------
	PauseForUser GraphProca
	DoWindow/K ITP
end

Function PositionButtonSasha()//(ctrlName) : ButtonControl // Calls up a pannel and pops up a cursor in the left upper corner on the screen. 
//You may then move the cursor to where you want it and click "add to image". A red cross is added to the image. You may repeat the operation several times.
// Results are stored in WMPositionMeasurementX,WMPositionMeasurementY- these readouts are actually almost correct !!!
//The function annoyingly adds white fields to the image. 
	//String ctrlName
	
	if(WMIsCalibrationApplied()==0)
	WMInitCalibrations()
	endif
	DoWindow/F WMMarkPositionsPanel
	if(V_Flag==1)
		AutoPositionWindow /E/R=$WMTopImageGraph()
		return 0
	endif
	
	NewPanel /K=1/W=(196.8,374,613.8,546) as "Marking Positions"
	
	DoWindow/C WMMarkPositionsPanel
	DoWindow/F WMMarkPositionsPanel
	
	String topImage=WMTopImageGraph()

	if(strlen(topImage)>0)
		AutoPositionWindow /E/R=$topImage
	endif
	SetDrawLayer UserBack
	DrawText 13,18,"Move the red cursor to the location that you want to mark"
	///-------New as of Feb 28, 2006
	Button activateGraph,pos={2,5},size={225,20},proc=ActivateGraph,title="Activate Graph"
	///----------------------------
	ValDisplay curXValdisp,pos={16,46},size={150,18},title="X:",limits={0,0,0},barmisc={0,1000},value= #"root:Packages:WMCalibrations:curX"
	Button addPositionButton,pos={10,83},size={150,20},proc=CustomAddPositionButtonProc,title="Add to image"
	CheckBox positionTagCheck,pos={74,112},size={120,20},title="Add with a tag", value=1
	CheckBox showLengthCheck,pos={73,139},size={160,20},title="Show position in tag", value=1
	ValDisplay curYValdisp,pos={233,46},size={150,18},title="Y:",format="%g",limits={0,0,0},barmisc={0,1000},value= #"root:Packages:WMCalibrations:curY"
	Button erasemeasurementTablebutton,pos={177,83},size={225,20},proc=ErasePoints,title="Erase measurements"
	Button markPositionDoneButton,pos={213,123},size={180,20},proc=CustomMarkPositionDoneButton,title="Accept This Frame"
	String cdf=GetDataFolder(1)
	String newDF=CustImageCursors(1,3)	
	SetDataFolder  newDF
	WM_DrawCursors(1)
	WMmakeMeasurementTextWave()			// saves text info
	SetDataFolder cdf
End

//************************************************************************************************
// The waves WMPositionMeasurementX & Y just save positions; they do not have NaNs
// separating segments because a segment consists of one point in each wave.

Function CustomAddPositionButtonProc(ctrlName) : ButtonControl
	String ctrlName
	
	String imageWaveDF=WMGetImageDF()
	String cdf=GetDataFolder(1)
	SetDataFolder imageWaveDF
	// see if the proper waves exist
	if((WaveExists(WMPositionMeasurementX)==0)  %|  (WaveExists(WMPositionMeasurementY)==0))
		// if the waves do not exist create them
		Make/O/N=1	WMPositionMeasurementX
		Make/O/N=1	WMPositionMeasurementY
	else
		// the waves exist so add to each 3 points to accomodate the new data
		Variable numPoints=numpnts(WMPositionMeasurementX)
		Redimension/N=(numPoints+1) WMPositionMeasurementX,WMPositionMeasurementY
	endif

	String cursorDF="root:WinGlobals:"+WMTopImageGraph()
	NVAR curX=root:Packages:WMCalibrations:curX
	NVAR curY=root:Packages:WMCalibrations:curY
	NVAR 	YOrigin=root:Packages:WMCalibrations:YOrigin
	NVAR 	XOrigin=root:Packages:WMCalibrations:XOrigin
		
	WMPositionMeasurementX[numPoints]=curX+XOrigin
	WMPositionMeasurementY[numPoints]=curY+YOrigin
	
	// now that we added the points, it is worthwhile to check if the waves are displayed.
	CheckDisplayed /W=$WMTopImageGraph() WMPositionMeasurementY
	if(V_flag==0)
		String cmd="AppendToGraph "+WMGetAxesForCommand()+" WMPositionMeasurementY vs WMPositionMeasurementX"
		Execute cmd
		ModifyGraph marker(WMPositionMeasurementY)=0,mode(WMPositionMeasurementY)=3,msize(WMPositionMeasurementY)=5
	       ModifyGraph rgb(WMPositionMeasurementY)=(16384,16384,65280)
	endif

	SVAR S_curMeasurement=root:Packages:WMCalibrations:S_curMeasurement
	SVAR S_curNote=root:Packages:WMCalibrations:S_curNote

	// now check if we need to add a tag based on the checkbox in the panel
	String position="("+num2str(curX+XOrigin)+","+num2str(curY+YOrigin)+")"
	ControlInfo positionTagCheck
	if(V_Value==1)
		String tagString=S_curMeasurement
		if(strlen(tagString)<=0)
			tagString="\Z09"+position
		else
			tagString+="\r\Z09"+position
		endif
		Tag /F=0/X=2 WMPositionMeasurementY,Pnt2x(WMPositionMeasurementY,numPoints), tagString	
	endif

	RemoveFromGraph/Z marker1y,markerDummyy
	KillWaves/Z marker1y
	String newDF=CustImageCursors(1,3)	// move it back to the original position
	SetDataFolder  newDF
	WM_DrawCursors(1)
	SetDataFolder cdf
End
///----------------------------------------------------------------------------------------------------------------------------------------
Function/S CustImageCursors(numCursors,type)
	Variable numCursors,type
	
	NVAR 	YOrigin=root:Packages:WMCalibrations:YOrigin
	NVAR 	XOrigin=root:Packages:WMCalibrations:XOrigin
	XOrigin=0
	YOrigin=0
	NVAR 	XUnitsPerPixel=root:Packages:WMCalibrations:XUnitsPerPixel
	NVAR	YUnitsPerPixel=root:Packages:WMCalibrations:YUnitsPerPixel
	XUnitsPerPixel=1
	YUnitsPerPixel=1
	// make sure the axes do not interfere with measurements; Because we are using markers
	// there will be enough room between the axes and the image.
	//ModifyGraph standoff=1
	NewDataFolder/O root:WinGlobals
	String cursorDF="root:WinGlobals:"+WMTopImageGraph()
	String cdf=GetDataFolder(1)
	NewDataFolder/O/S	$cursorDF
	
	SetWindow $WMTopImageGraph(), hook=CustomCleanOnKill
	String/G S_TraceOffsetInfo=""
	Variable/G dependentVar

	// the following dummy waves are used to shift the axes away from the image.
	Variable xmin,xmax,ymin,ymax
	Wave theImage=$WMGetImageWave(WMTopImageGraph())

	xmin=XOrigin-0.05*DimSize( theImage,0)*XUnitsPerPixel
	ymin=YOrigin-0.05*DimSize(theImage,1)*YUnitsPerPixel

	xmax=XOrigin+1.05*DimSize( theImage,0)*XUnitsPerPixel
	ymax=YOrigin+1.05*DimSize(theImage,1)*YUnitsPerPixel
	Make/O /N=3 markerDummyy={ymin,NaN,ymax}
	Make/O /N=3 markerDummyx={xmin,NaN,xmax}

	// the following is the actual cursor cross.  It may be drawn too small to be see.
		
	if(numCursors==1)
		Make/O marker1y={15,25,NaN,20,20},marker1x={10,10,NaN,5,15}
		
		marker1y=marker1y*YUnitsPerPixel+YOrigin
		marker1x=marker1x*XUnitsPerPixel+XOrigin
		Make/O/N=1 markerLineX,markerLineY
			markerLineX[0]=marker1x[0]
			markerLineY[0]=marker1y[0]
	
		if(type==0)
			SetFormula dependentVar,"WMUpdateMarkerOrigin(S_TraceOffsetInfo)"
		else
			if(type==3)
				SetFormula dependentVar,"WMUpdateMarkerPosition(S_TraceOffsetInfo)"
			endif
		endif
	else
		Make/O marker1y={15,17,NaN,16,16},marker1x={10,10,NaN,9,11}
	endif
	
	SetDataFolder cdf
	
	return cursorDF
End
//************************************************************************************************
Function ErasePoints(ctrlName) : ButtonControl
	String ctrlName

	String imageWaveDF=WMGetImageDF()
	String cdf=GetDataFolder(1)
	SetDataFolder imageWaveDF
	// first remove from the graph
	CheckDisplayed /W=$WMTopImageGraph() WMLengthMeasurementY
	if(V_flag==1)
		RemoveFromGraph/Z WMLengthMeasurementY
	endif
	CheckDisplayed /W=$WMTopImageGraph() WMAngleMeasurementY
	if(V_flag==1)
		RemoveFromGraph/Z WMAngleMeasurementY
	endif
	CheckDisplayed/W=$WMTopImageGraph() WMPositionMeasurementY
	if(V_flag==1)
		RemoveFromGraph/Z WMPositionMeasurementY
	endif
	
	CheckDisplayed/W=$WMTopImageGraph() markerDummyy
	if(V_flag==1)
		RemoveFromGraph/Z markerDummyy
	endif
	
	// now we kill the actual waves in which the measurement data are stored
	KillWaves/Z markerDummyy,markerDummyx//,WMPositionMeasurementX,WMPositionMeasurementY
	Make/O/N=0 WMPositionMeasurementX,WMPositionMeasurementY    ///Changed as of March 03, 2006

	SetDataFolder cdf
End

//************************************************************************************************
Function CustomMarkPositionDoneButton(ctrlName) : ButtonControl
	String ctrlName
	String cursorDF="root:WinGlobals:"+WMTopImageGraph()
	String cdf=GetDataFolder(1)
	SetDataFolder cursorDF
	RemoveFromGraph/Z marker1y,markerDummyy
	DoWindow/K WMMarkPositionsPanel
	SetFormula dependentVar,""
	SetDataFolder cdf
	ModifyGraph standoff=0
	//DoWindow/K ProcaWindow
	DoWindow/K GraphProca	
End

///__________________
Function CustomCleanOnKill (infoStr)
	String infoStr
	String event= StringByKey("EVENT",infoStr)
	if(cmpstr(event,"kill")==0)
		WMSetOffsetAngleDone("")
		WMMeasureAngleDoneButtonProc("")
		WMMarkPositionDoneButtonProc("")
		WMMeasureLengthDoneButtonProc("")
		//WMSetOriginDone("")
		DoWindow/K WMSpatialMeasurementPanel
	endif
	return 0				 
End