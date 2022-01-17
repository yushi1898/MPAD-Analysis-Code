
gfitgridall()#pragma rtGlobals=1		// Use modern global access method.

#pragma rtGlobals=1		// Use modern global access method.#pragma rtGlobals=1		// Use modern global access method.
#include  <All IP Procedures>

// *************MODIFICATION HISTORY  - EFFECTIVE MAY 19, 2006



//   Rules:   List any procedure files that are changed  - noting the date of the change
//                Details on the changes should be given in the procedure file itself.
//                The changes should be listed at the top of the file (date stamped) and also noted at the point of change in the code

//12/28/21 Modified panel layout to move MSD plot after generating bifurcated mask
//12/28/21 Made updates in CalcMSD5M to include generating figure 2 in PNAS
//12/28/21 Chnaged the button name of "cont Fit All Posts"
//12/28/21 Changed Double Freq button to allow user to choose which video to perform DLI upon.
//12/28/21 Changed Digital Lock in to only plot heatmap of f = 0.1 Hz in cont double freq 
//12/28/21 Changed CalcPostCentersAndShifts to remove calculation of average x and y spacing
//12/12/21  started changes to magnetic rheology analysis
//12/10/21 Change hard-coding of microns per pixel.  Now defined by global MicronsPerPixel and referenced by NVARs throughout all routines. User should change as needed.
//12/8/21 Change GFitAll to allow suppression of images and diagnostic printout
//12/07/21  Change movie handling to kill movie after each routine uses it - must then be reopened.
//do by passing the variable MoviePath as a Global String
//12/6/21 saving of error files in GfitAll() and CalcPostCentersAndShifts()  deleted as no longer created by centroid fit
// //  removed: viewGfitErrors() (and associated .ipf file from includes), LoadGFit(), PostCenterLoad(), 

//12/5/21 YS moved MSD plotting functions from MSD to PlotALot; changed PlotALotG and (renamed) ComputeMSD
// 12/3/21 Change button names and re-order positions
// 12/3/21  Started removing Gfit routines from here and the FitAllPosts....ipf file
//09/16/21  Begin modifications by DHR to Work in Igor 9
//  
// 09/09/21   Modified by YS to simplify main panel

//9/7/21 DHR:  Fixed calccrosscorrelation2.ipf and CalcMSD1_9Ma.ipf to compile in Igor 9.

//latest version for large posts on passive data            F//11/29/18
//quick calulcation of MSD at 100s   11/3/2018


//5/8/2018     suit for 100 fps to 10 fps averaged.
//11/10/2016  changed MSD to enable both frequency dependence and single file analysis
//5/10/2016    changed makeitquiver function to calculate avg force 
//9/20/15    change original MSD to rescaled MSD, changes made in calcMSD function
//9/9/2015  added initial shift when calculating digital lock in 
//11/3/14    added allowing analysis videos continuously

//10/30/14  added button for DLI calculation of double frequency measurement and correaltion measurement

//7/14/14    enlarge panel size and add buttom for heatmap
//6/13/14    change from gaussian fitting to centroid fitting
//1/3/14      added in plot a lot to enable whether include bgshift

//12/12/13   added the bottom for calculating background shift and MSD             Yu Shi
//12/12/13   added wave for storing background shift wave                                 Yu Shi
//12/12/13   good for analysis white post


// 12/10/13  Modifications by DHR to allow working with movie files - immediate goal is Yu Shi's mPAD microrheology data
//        .ipf files with "M" at end of name...

// 12/9/13
// Line temporarily commented out in Function FindForceCutoff - ASK CRAIG!!


// 2/13/2012 Craig- added force filtering to filter out "noise" posts when calculating tugging forces

//1/27/2012 Craig - makeitquiver now treats guide posts as ignored posts.  This is to visually differentiate, in makemaskGUI, between posts being ignored outright, and posts being ignored to 
					//   separate a cell pair.



// 12/15/2011 Craig - MakeMaskGUI02,  allows for drawing of polygon around cell to then be filled with cell posts


// v 3.5 4/14/2010   Craig and Dan's initial modification to do "update Gfit" to all posts and frames automatically
//          adds "UpdateGFitAll".  Also added prompt and flagging for empty post vectors (for finding errors), and post height (spring constant)


// 4/23/10 plotalot now uses min,max of data to set limits for plots

//7-13-09 Writing "dIspArray" waves to file, for use in separate procedure to do batch analysis 
	//Saved "Topractice" as "disparray_string"
//7-10-09 Changed error in energies. Wasn't dividing by total number of cell posts!!!

// 04/04/07    Error Check put in FItAllPostsGErrChk
//                       NEW:   ViewFitErrors - for scanning for large fit errors
//                                  Creates grayscale image of errors post vs  , and a histogram of errorsize for finding outlyers

// 04/03/07  "RHOD000" or "RHOD00" naming convention set by static string constant
//                          ksRhodFileRoot = "RHOD0" for 3-digit names
//                                                    = "RHOD" for 2-digit names
// this constant must be set in 
//       GfitCompleteInterpAll.ipf, FitAllPostsG.ipf,  DefineArray.ipf, and ExamineFitG.ipf
//                   WARNING: WON'T WORK FOR RUNS WITH > 99 FRAMES

// 3/28/07  - Revising MakeItQuiver for new GFit paradigm
//                   - shifts calculated with Calc Centers and Shifts button only  - not in fitall

// 3/24/07  Modified to do grid interpolation for every frame

//   LIST OF CHANGES TO OTHER PROCEDURE FILES
// 7/10/06 Various changes made to Examine Bottom Fit and Examine Fit segments.
// Now files are stored in a separate folder, just like in a PSF routine

// 6/28/06  Modifications begun to bring Gfit up to date to parallel PSFfit 
 //                      error bar propagation etc.  
 
// 5/31/06 ExamineFitBaseG.ipf replaces ExamineFitBase.ipf
//                  ExamineFitG.ipf   replaces ExamineFitLUT.ipf
//                  FitGuidePostsG.ipf - doGfit() :  fit ROI clipped to remain in image

//    LIST OF CHANGES TO THIS FILE

//5/31/06  
//                  Routines in this file that are changed from their counterparts in PSFfitComplete.ipf have "G" suffix added

//5/26/06 MakeItQuiverPSF changed to MakeItQuiverG. Significant differences in assinglements of variables. Also using fitBaseWave
//  PlotALot changed to PlotALotG: ApplyShiftPSF changed to ApplyShift
// 5/26/06   GfitComplete created from 5/25/06 version of PSFfitComplete
//*******************************************

#include ":FitAllPostsGerrChkHexMv07"   ////12/12/21 movie suffix passed 12/8/21 updates DHR 12/3/21   DHR modified to remove old Gaussian fit (GfitAllPosts)
#include  ":HelperRoutinesM5"   ////12/7/21 updated 9/7/21 DHR updated for Igor9  12/3/21 modified for all files in movie prompt
#include  ":DefineArray09M"   //9/13/21, 12/3/21 DHR updated 12/5/21 DHR updated
//#include  ":ExamineFitGM1.12"   // 12/6/21 DHR removed - no longer used
#include  ":PlotAlot8M"    //12/9/21 YS modified; 12/7/21 DHR modified 12/5/21 YS modified
#include  ":MakeTime3"
// #include ":ViewGFitErrors"  // 12/6/21 DHR removed - no longer used
#include ":ModernClicking04"  //9/16/21 DHR updated
#include ":MakeMaskGUI06M"  //12/7/21 DHR 9/9/21 modified by YS  12/5/21  modified by DHR
#include ":Correlation"
#include ":CalcMSD5M"   //12/10/21 nm/um to pixels fixed 12/6/21 DHR started fixing wave refs for Igor 9 12/5/21  YS updated.  9/7/21  DHR updated this file
#include ":Digital lock in4"  //12/28/21
#include ":Histogram2"  //  12/10/21 nm to pixels fixed  12/6/21  YS NEEDS TO REVIEW AND UPDATE
#include ":activestress calculation2"  //12/10/21 nm to pixels fixed
#include ":calccrosscorrelation2"   // 9/7/21 DHR updated
#include ":Frequency Dependence3"  //12/16/21 
#include ":PlotMSDHeatMapV1"  //12/28/21      Yu updated to generate PNAS heatmap plots
//#include ":Levyanalysis_v1_5"  This is not included in the published version.  Contains code to further analysis
//  van Hoves,  avalance detection etc. 



//make all files writeable
newpath writepath, "Mollenhauer HD:Users:Corinne:GFitAnalysis:GfitAllPosts:GfitComplete"
setfilefolderinfo /p=writepath /ro=0 /r=1

// constant for datafile naming convention control

//static strconstant  ksRhodFileRoot = "RHOD"  // for 2-digit names
static strconstant  ksRhodFileRoot = "RHOD"  // for 3-digit names

// maskwave constants  added 6/28/06

static constant kEmpty = 0
static constant kCell = 1
static constant kGuide = 2
static constant kIgnore = 3
static constant kWirepost = 4
static constant kFitAll = -1
static constant kInterp = 5     // 5/24/06  new constant interp base fit

// indices for fitparams
static constant kNfitParamsG  = 9  // includes columns for chisq and FitErrorReporting

static constant kFbkgnd = 0
static constant kFAmpl = 1
static constant KFMSD = 1
static constant kFX0 = 2
static constant kFSigX= 3
static constant kFY0 = 4
static constant kFSigY = 5
static constant kFCorr = 6
static constant kFChisqG = 7
static constant kFerrorFlag = 8

// constant value for DLI heatmap plot range
static constant DLIHeatmap_low = 0.0001   //12/28/21 minimum value of DLI deflection magnitude shown in heatmap
static constant DLIHeatmap_high = 0.01    //12/28/21 maximum value for DLI deflection magnitude shown in heatmap


function GfitGridAllM()  ///------------This is an envelope function popping up menu with buttons
						///-------------The rest is going to be procedures for each button

// 12/10/13 - variables for input (data) movie handling
Variable/G InputMovieID_G
Variable/G NmovieFrames_G   // number of frames in movie
String/G MoviePath          // added 12/7/21 

Variable/G MicronsPerPixel = 0.125   //12/10/21 lengthscale conversion for 40x objective, Nikon T/E 2000/E and Allied Vision GX1050 camera at JHU



Variable/G WindowFlag=1, ClickVal, doneflag = 0
String/G RootString="Default"
String/G PFS_String="PSF_"+RootString
String/G LUTr_String="LUTr_"+RootString
String/G LUT_String="LUT_"+RootString
String/G MaskString="Mask_"+RootString
String/G MaskString2= MaskString +"2"
String/G FitBottom="FB_"+RootString
String/G FitErrBottom="FBErr_"+RootString
String/G FitRes_String="FitRes_"+RootString
String/G FitErr_String="FitErr_"+RootString
String/G AttributeString="Att_"+RootString
String/G AttributeString2="Att_"+RootString +"Mask2"
String/G ShiftString="Shift"+RootString
String/G PFS_Arr="PSFArr_"+RootString

String/G FitresGbaseStr = "FRGbase_" + RootString   //5/31/06 for full results of G fit for base (currently on Rhod000)
String/G FiterrGbaseStr = "FEGbase_" + RootString



//make global string for background shift
string/G  bgshiftstring="bgshift_"+RootString 



DoWindow/K GFit_GridEveryFrame
NewPanel/K=1/N=GFit_GridEveryFrame/W=(100,100,385, 700)	///--------------------------------------------------------------------------------------
SetDrawEnv linefgc= (0,0, 53520),fillpat= 0, linethick=2;	DrawRect 5,5,265,65
SetDrawEnv linefgc= (0,0, 53520),fillpat= 0, linethick=2;	DrawRect 5,120,265,275
SetDrawEnv linefgc= (0,53520,0),fillpat= 0, linethick=2;	DrawRect 5,305,245,520
SetDrawEnv linefgc= (64512,14848,14848),arrow= 1,linethick= 4.00,arrowfat= 1.00; DrawLine 275,10,275,370

	SetDrawLayer UserBack
      String/G PWF="Huller HD:Users:Corinne:GfitAnalysis:Data:Cell20" ///This can be changed (computer specific)

// 12/10/13 -anything with "M" suffix changed for movie input data
	Button SWF,pos={15,10},size={110,20},proc=SetWorkingFolder,title="Set Folder "  // 12/10/13 moved over
	Button InputMovie,pos={135,10},size={120,20},proc=PickInputMovie,title="Pick Input Movie "  //12/10/13 new for movies
	SetVariable RootStringValue,pos={10,40},size={180,20},title="Working Name", Value= RootString
	Button AcceptName, pos={200,40}, size={60,20}, title="Accept", proc=AcceptName
	Button ArrayDefine, pos={40, 90}, size={180, 20}, title="Define Array and Mask", proc=ArrayDefineM    // changed 12/10/13 Find  array lattice constants,  micronsperpixelImage, angles, initial shifts, initial guesses for 

	Button GFitAll, pos={10,130}, size={130,20}, title="Fit All Posts", proc=GFitAll 
	Button Calcshift, pos={150,130}, size={110,20}, title="Calc_bgshift",proc=calcbgshift

	Button CalcPostCenters, pos={10,160}, size={130,20}, title="Calc. Grid Center", proc=CalcPostCentersAndShifts // 	
	//adeded 12/12/13  Yu Shi
	Button ComputeMSDButton, pos={150,160}, size={110,20},title="MSD", proc=ComputeMSD

	Button MSDmap, pos={10, 220},size={130,20}, title="MSD Heat Map",proc=calcMSDheatmap
	Button MSDstatscalc, pos={150,220}, size={110,20}, title="MSD stats Plot", proc=calcMSDstats

	Button  CELLMASK, pos={10,190},size={130,20}, proc = gencellmask, title="Make Cell Mask MSD"
	Button  FORCEBIFURCATION, pos = {150,190}, size = {110,20}, proc = bifcell, title = "Force Bifurcation"
		
	Button PlotALot, pos={10,250}, size={130,20}, title="Plot-A-Lot", proc=PlotALotG
	Button MakeItQuiverPSF,pos={150,250},size={110,20},proc=MakeItQuiverG,title="Make it Quiver"  // Modified May 29, 2006
	
	
	Button DFDLIheatmap, pos={130, 340}, size = {110, 20}, title = "Double Freq", proc = calcDFDLIheatmap
	
	Button CONcalshift, pos={130, 310}, size = {110,20}, title = "cont calc bgshift", proc= COUNTCalcbgshift
	Button CONgfitall, pos = {10,310}, size = {110,20}, title = "con Fit All Posts", proc = COUNTGFitAll
	Button CONDLIDFheatmap, pos={10,340}, size = {110,20}, title = "cont double freq", proc=COUNTcalcDFDLIheatmap
	
	Button  READWAVE, pos={30,370},size={200,20}, proc = b_readwaves, title="Load posts DLI"
	Button  RELATIVESIG, pos = {30,400}, size = {200,20}, proc = b_relativeplot, title = "Correct time variance"
	Button  READHALLSEN pos={30,430},size={200,20}, proc = b_readhall, title="Load Hall sensor readout"
	Button  CALCDLIMAG, pos = {30,460}, size = {200,20}, proc = b_calcmagDLI, title = "Calc Hall sensor DLI"
	//Button  Vanhoffcurve, pos={10,80},size={250,20}, proc=genvanhoff, title = "Generate Vanhove"
	Button  CALCMODU, pos={30,490},size={200,20}, proc=b_calcmodulus, title = "Calc Modulus"		
	
	
	
	
	Button doneButton, pos={10, 560}, size={240, 25}, title="\\K(65280,0,0)\f01Done", proc=DoneButton  //Done
	
end
///---------------ArrayDefineM -------   
//12/10/13 - changed for movie input      
function ArrayDefineM(ctrlName): ButtonControl //1. Button ArrayDefine, pos={40, 40}, size={180, 20}, title="Define Array", proc=ArrayDefine  //Done  //Find  array lattice constants,  micronsperpixelImage, angles, initial shifts, initial guesses for 
string CtrlName //amplitudes, backgrounds, etc.       //OUTPUT :  constantswave to hold various useful stuff
	DefineArray2M()              
end           

// --------------GFitAll ------------

// does centroid fitting.
//12/8/21 Change  to allow suppression of images and diagnostic printout:  centroidgfit2()
//12/6/21  removed saving of fiterrors wave as no longer generated in centroid method
// 12/3/21   initial guess removed from centroidgfit()
function GFitAll(ctrlName): ButtonControl
    string CtrlName
    
    NVAR NmovieFrames_G  // total number of movie frames
    string s1,s2
    variable flast = NmovieFrames_G-1
    variable showimageflag = 0 
    variable printeveryNframes
    sprintf s1, "Set frames to fit (%d total)", NmovieFrames_G
    //try to set this for centroid fitting
    variable f1=0,f2=0
    Prompt f1, "Enter first frame to fit (0?):"     
    sprintf s2, "Enter last frame to fit (%d?):", flast
    Prompt f2, s2
    Prompt showimageflag, "Set = 1 to show image of each frame during fit (will be slower!)" 
    printeveryNframes = max(floor(NmovieFrames_G/10),1)
    Prompt printeveryNframes, " Print frame number every N frames to track progress (= -1 to suppress). N = "   
    DoPrompt s1, f1, f2, showimageflag, printeveryNframes
//    if (V_Flag==1)
//		return 0
//    endif          
    SVAR MaskString, FitRes_String, FitErr_String, ShiftString      
    SVAR  AttributeString
    LoadWave/O/P=ToPractice AttributeString+".ibw"
    LoadWave/O/P=ToPractice ShiftString+".ibw"
    Duplicate/O $AttributeString ASTr
 	

   centroidgfit2(f1,f2,ASTr[0],ASTr[1],ASTr[3],ASTr[4],ASTr[6],ASTr[7],$MaskString,$ShiftString,showimageflag,printeveryNframes) //for back compatible code
//12/5/21  removed saving of fiterrors wave as no longer generated in centroid method
   WAVE FitResults
	Duplicate/O FitResults $FitRes_String
	Save/O/P=ToPractice $FitRes_String  as FitRes_String+".ibw"	
	
end 




////---------------------------------------------------------
//** CalcPostCentersAndShifts() 
//12/6/21 saving of error files deleted as no longer created by centroid fit
//   3/24/07 - Calculates deflections by doing gridding for each frame
// For backward compatibility, uses "fitbottom"  ie FB_ and FBerr_ strings and filenames
//   to store data.  
// actually computes Undeflected positions
// 05/19/06  Autorefit capability added (refits each cell post with kfitrange = 0.3 instead of 0.4)

function CalcPostCentersAndShifts(ctrlName): ButtonControl 
string CtrlName   
    NVAR NmovieFrames_G  // total number of movie frames
    string s1,s2
    variable flast = NmovieFrames_G-1
    sprintf s1, "Set frames to analyze (%d total)", NmovieFrames_G
    sprintf s2, "Enter last frame (%d?):", flast
    //try to set this for centroid fitting
    variable f1=0,f2=0
    Prompt f1, "Enter first frame (0?):"     
    Prompt f2, s2     
    DoPrompt s1, f1, f2
 	 if (V_Flag==1)
		return 0
	 endif          


	SVAR MaskString, FitRes_String, FitErr_String, ShiftString      
	SVAR  AttributeString,FitBottom, FitErrBottom
	

        variable npost ,nframe
        nframe = f2 - f1
        npost = dimsize($fitres_string, 0)	
	

	LoadWave/O/P=ToPractice AttributeString+".ibw"
 	Duplicate/O $AttributeString ASTr
	LoadWave/O/P=ToPractice ShiftString+".ibw"
	LoadWave/O/P=ToPractice FitRes_String+".ibw"
//	LoadWave/O/P=ToPractice FitErr_String+".ibw"
	Duplicate/O $FitRes_String FitResults
//	Duplicate/O $FitErr_String FitErrors
//	Duplicate/O $FitRes_String ZeroPos
        make/o/n=(npost,2,nframe)  zeropos
        duplicate/o/R=(*)(*)(f2-f1) $fitres_string zeropos
//	Duplicate/O $FitErr_String ZeroPosErrors
	//Duplicate/O $ShiftString MyShiftWave

//changed for centroid method   7/11/14

       CalcInterpPostCenters(f1,f2,ASTr[0],ASTr[1],ASTr[3],ASTr[4],$MaskString,ZeroPos,0,1)
       
       	// calculate shifts relative to first frame for backward compatibility
       calcShiftsGrelF1rot(f1,f2,ASTr[3],ASTr[4],$MaskString,$ShiftString,FitResults)

       
       Duplicate/O ZeroPos $FitBottom
//       Duplicate/O ZeroPosErrors $FitErrBottom  //12/6/21   Deleted as ZeroPosErrors no longer created
	Save/O/P=ToPractice $FitBottom  as FitBottom+".ibw"	
//	Save/O/P=ToPractice $FitErrBottom as FitErrBottom+".ibw"
	Save/O/P=ToPractice $ShiftString as ShiftString+".ibw"  
	

     //**********************************************************************************      
	
end 

//------------------ PostCenterLoad ---------------------

// 3/25/07
function PostCenterLoad(ctrlName): ButtonControl 
string CtrlName  
	SVAR RootString, FitRes_String, FitErr_String, MaskString, ShiftString
	SVAR  AttributeString, FitBottom, FitErrBottom
	NVAR hexcp
	LoadWave/O/P=ToPractice AttributeString+".ibw"
	Duplicate/O $AttributeString ASTr
	
		Variable/G kfitrange=floor(ASTR[5]*0.4)
	
	
	LoadWave/O/P=ToPractice MaskString+".ibw"
	LoadWave/O/P=ToPractice ShiftString+".ibw"
	LoadWave/O/P=ToPractice FitRes_String+".ibw"
	LoadWave/O/P=ToPractice FitErr_String+".ibw"
	Duplicate/O $FitRes_String FitResults
	Duplicate/O $FitErr_String FitErrors
	LoadWave/O/P=ToPractice FitBottom+".ibw"
	LoadWave/O/P=ToPractice FitErrBottom+".ibw"
	Duplicate/O $FitBottom ZeroPos
	Duplicate/O $FitErrBottom ZeroPosErrors
end



///------------------------------------------------------
// PlotAlotG
//12/5/21  YS edited to include MSD plots
// 3/7/06 - modified to allow adding shifts
// 3/25/07  Updated for GRID EVERY FRAME

///------------------------------------------------------
Function PlotAlotG(ctrlName):ButtonControl
string CtrlName
SVAR FitRes_String, FitErr_String,ShiftString, AttributeString,MaskString,FitBottom,FitErrBottom,bgshiftstring
variable posti=0,postf=1 ,param1=3,param2=4,shiftflag = 1,printflag=0, notebookflag = 1, timeflag = 0,cellpostflag=0, bgshiftflag=1, spikeflag=0
//string p1string="X",p2string="Y"
variable framerate=10
string p1string,p2string,masktype_string,tracetypestring
variable f1,f2
wave linmaskval  // manually created mask wave
wave linmaskval_intercept_int  // mask wave after force bifurcation

LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

postf = AStr[3]*Astr[4]-1

String LastPost="Last post"

Prompt posti,"First post"
Prompt postf, LastPost
Prompt p1string,"Print in Left Column",  popup "X"
Prompt p2string,"Print in Right Column",  popup "Y;MSD"
Prompt masktype_string,"Mask Type", popup "Manual;Automatic"  // manual uses original manually created mask;  Automatic uses mask made in Force Bifurcation.


Prompt notebookflag, "type 1 to make notebook"
Prompt cellpostflag, "Type 1 to plot only cell posts"
Prompt tracetypestring, "Trajectory to plot", popup "Raw;Dedrifted;DeflectionFromGrid"
Prompt framerate, "type framerate"
prompt spikeflag, "type 1 to remove spikes"
DoPrompt "PlotALot Parameters", posti,postf,p1string,p2string,cellpostflag, notebookflag, tracetypestring, framerate,masktype_string,spikeflag   


if (V_flag==1)	
		return 0
endif
//Added July 27, 06
if (stringmatch(p1string, "X")==1)
	param1=kFX0
endif

if (stringmatch(p2string, "Y")==1)
	param2=kFY0
elseif (stringmatch(p2string, "MSD")==1)
   param2=kFMSD
endif

if (stringmatch(tracetypestring,"Raw")==1)
	bgshiftflag = 0
elseif (stringmatch(tracetypestring,"Dedrifted")==1)
	bgshiftflag = 1
elseif (stringmatch(tracetypestring,"DeflectionFromGrid")==1)
	bgshiftflag = 2
endif



duplicate/O $FitRes_String  fitres_plotAlot


f1 = dimoffset(fitres_plotAlot,2)
f2 = dimsize(fitres_plotAlot,2) + dimoffset(fitres_plotAlot,2) -1

if (stringmatch(masktype_string,"Manual")==1)
	PlotALot4EB(posti,postf,param1,param2,p1string,p2string,ASTr,linmaskval,fitres_plotAlot,$bgshiftstring,cellpostflag,notebookflag,bgshiftflag, framerate, spikeflag)  // //6/02/06 - eb added
elseif(stringmatch(masktype_string,"Automatic")==1)
	PlotALot4EB(posti,postf,param1,param2,p1string,p2string,ASTr,linmaskval_intercept_int,fitres_plotAlot,$bgshiftstring,cellpostflag,notebookflag, bgshiftflag, framerate, spikeflag)  // //6/02/06 - eb added
endif
end


// ---------------------- MakeItQuiverG() ------------

//4/26/10 (Craig) Added prompt and flagging for empty post vectors (for finding errors), and post height (spring constant)

// 3/28/07 Adapted for new method - GRID EACH FRAME

//7/7/06  - Error bars and weighted averaging done for avg. displacement calcs - see DHR Book 7, p. 73
// plot vs timewave option added
///----------------------Static Constants added on May 29,2006
// static constants moved to top of file 5/31/06
//4/23/09 added option to plot index number to identify misfitted posts
//1/24/14  changed to analyse movies    Yu Shi    

//12/7/21  WILL HAVE TO BE REWRITTEN FOR IGOR9 MOVIE HANDLING! 

Function MakeItQuiverG(ctrlName) : ButtonControl  ///-------------------Makes Quiver Movie of Relative Displacements
string ctrlName
//---------
Variable f1=0, f2=0, fps=1, shiftflag=0, MFact=1, timeflag = 0, Fscale=10, Afact=1,AX=50,AY=0,fzero=40,indexflag=1,bArrowflag=1,kspringflag=6,interplineflag=0, freq=100
Prompt f1, "Enter First File Index"
Prompt f2, "Enter Last File Index"
Prompt fps, "Frames per Second"
Prompt MFact, "Multiplication Factor"
Prompt shiftflag, "0 for GridEach,  1 for Shift"  //3/28/07
Prompt timeflag, "Type 1 to plot avgdisp vs timewave"
Prompt fzero, "Enter Frame for Field ON (t = 0)"
//Prompt Fscale, "Force Scale Bar [10-99nN]" //Temporarily removing prompt for errorbar size
Prompt Afact, "Arrowhead Scale [0-3]"
Prompt AX,"Scale Bar X offset"
Prompt AY,"Scale Bar Y offset"
Prompt indexflag,"Type 1 to plot indexnumbers"

DoPrompt "Indexes", f1, f2, fps,MFact, shiftflag, timeflag, Afact,AX,AY,indexflag

//Craig added prompt 4/26/10, and 1/30/12
Prompt bArrowflag,"Enter 1 to plot empty post arrows"
Prompt kspringflag,"Which size posts did you use? enter 4,5,6,7,or 10"
Prompt interplineflag, "Show interpolation lines and base positions?"
Prompt  freq, "enter analysis frequency"
DoPrompt "More Options", bArrowflag,kspringflag,interplineflag, freq
if (V_flag==1)	
		return 0
endif	

SVAR MaskString, FitRes_String, FitErr_String, RootString, ShiftString, FitBottom, FitErrBottom
wave ASTR
wave timewave   // must be created (by hand?) elsewhere
wave fitresultsbase  // craig added 1/30/12
NVAR InputMovieID_G, Nmovieframes_G


Duplicate/O $ShiftString ShiftWave
Duplicate/O $FitRes_String FRS, FRdisp  // 3/28/07
Duplicate/O $FitErr_String FRE, FREdisp // 3/28/07
Duplicate/O $FitBottom fitresultsbase   // these next two are now 3d waves!!! for compatiblity with shift use [][][f1]
Duplicate/O $FitErrBottom fiterrorsbase

 calcdispGridEach(f1,f2,ASTr[3],ASTr[4],FRdisp,FREdisp,$FitBottom,$FitErrBottom)
// FRdisp and FREdisp now hold grid-each-frame displacements and errors

	ApplyShiftsG(f1,f2,ASTr[3],ASTr[4],$ShiftString,FRS,FRE) //Added May 29, 2006
// note that this is different from the PSF version, which doesn
	Variable dx,dy,length // added 7/7/06
	Variable sumlen, sumlenerr, sumlenblue, sumlenerrblue
	Variable energy, energyerr, energyblue, energyerrblue // 8/28/06
	variable kspring   //added by CMK 9/9/08 //changed 
	variable lattice //lattice in um
	
	variable energy_c2, energyerr_c2, sumlen_c2, sumlenerr_c2
	

	 variable energy2, energy0

	//This depends on your substrate!!! 18.1 for M5, 15.7 for M6 and 22.3 for M4...craig added 4/26/10
	if(kspringflag==3)
		kspring = 33.03
	elseif(kspringflag==4)
		kspring = 22.3
	elseif(kspringflag==5)
		kspring = 18.1
	elseif(kspringflag==6)
		kspring = 15.7
	elseif(kspringflag==7)
		kspring = 11.5
	elseif(kspringflag==10)
		kspring = 5.4
	endif
	
	
	lattice =  4
	
	
	Variable i, j, k, HorizPosts, VertPosts
	Variable TotalNumPosts, QuiV, CellOccupado, BlueOccupado
	Variable AngleOfRotation=ASTR[2]
	variable numfiles = f2-f1+1
	
	
	HorizPosts=ASTr[3]
	VertPosts=ASTR[4]
	TotalNumPosts=HorizPosts*VertPosts
	
	//write to file dx and dy, deltax, deltay into DispArray_filestring
	make/o /n=(TotalNumPosts,4,numfiles) disparray=0

	make/O/N=(TotalNumPosts) makeshiftX, makeshiftY, blueshiftX, blueshiftY, makeshiftX_c2, makeshiftY_c2
	make/O/N=(TotalNumPosts,2) lengthandangle, lengthandangle_c2     //added 6/1/2016
	make/O/N=(TotalNumPosts,2) bluelengthandangle
	make/o/n=(TotalNumPosts) lengthchangewave, lengthchange0, lengthchangemid //added 2/17/09
	make/O/N=(TotalNumPosts) lengtherr,bluelengtherr   // added 7/7/06
	make/O/N=(f2-f1+1) AvDisp,  BlueDisp, AvErr,BlueErr, TotalE, TotalE_c2, TotalE2, TotalEerr, BlueTotalE, BlueTotalEerr, EnergyChangeWave// added 6/2/06, 8/28/06
        make/O/N=((f2-f1+1), 2) totalforce=0, totalforce_c2=0
        make/O/N=(f2-f1+1) TotalF                                           //added 4/16/18
        variable avgE, avgF


//*************************************************************************************************
// Craig 2/16/12 changing the below line of code so that makitquiver uses the mask from linmaskval, which is modified by force filtering
	//Duplicate/O $MaskString CellMatrix
	wave linmaskval
	Duplicate/O linmaskval CellMatrix
//*************************************************************************************************	
	//*****
//	wave maskwave
	//Duplicate/o maskwave cellmatrix
	//*******
	//wave tempmask
	//Duplicate/o tempmask cellmatrix
	//****************************************************************
	// Since cellmatrix is no longer from maskstring, i remove this line too
//	matrixtranspose CellMatrix   
	//****************************************************************
	string s1,s10
	//added 1/23/14    suitable for movie analysis   YuShi
      playmovieaction setfrontmovie= inputMovieID_G
      playmovieaction  stop, frame=f1
      playmovieaction  extract
      imagetransform rgb2gray M_movieframe
      duplicate/O M_rgb2gray image1 
	Duplicate/O image1 OrigImage1
	ImageRotate/A=(AngleOfRotation) OrigImage1
	Duplicate/O M_RotatedImage RotatedOrig1
	DoWindow/K RotatedImage
	NewImage/K=1/N=RotatedImage RotatedOrig1
//	ASTR[6] = 1000
//	ASTR[7] = 2200
	ModifyImage RotatedOrig1 ctab= {ASTR[6],ASTR[7],Grays,0}
	Variable/G ImScale = 0.5
	ModifyGraph width=ASTR[3]*ASTR[5]*ImScale,height=ASTR[4]*ASTR[5]*ImScale
	NewMovie/O/P=ToPractice/F=(fps)/L as RootString+".mov" // loop over raw data files
//	NewMovie/O/P=ToPractice/F=(fps)/L as RootString+"_DIC.mov" // loop over raw data files
	//changed 5/31/2014 by Yu
	duplicate/O frdisp corrwave //
	corrwave = 0
	

if (timeflag == 1)
     MakeTime(f1,f2)
     wave vTime
     variable timezero = vtime(fzero)
     vTime -= timezero
endif

	      disparray[][0][]= frdisp[p][0][r]
		disparray[][1][] = frdisp[p][1][r]
	
		disparray[][2][] = fredisp[p][0][r]
		disparray[][3][] = fredisp[p][1][r]
		
		String/G DispArray_String="DispArray_"+RootString
	duplicate/o disparray $DispArray_String
	Save/O/P=ToPractice $DispArray_String  as "DispArray_"+RootString+".ibw"

      playmovieaction setfrontmovie= inputMovieID_G
      playmovieaction  stop, frame=f1

	
	for (j=f1; j<=f2; j+=1)
	
	
		
	
       sprintf s1, "%s%03d.tif", ksRhodFileRoot,j   // 4/3/07
//		sprintf s1, "RHOD%02d.tif",j
//		sprintf s1, "DIC%02d.tif",j
//added to analyse movie     Yu Shi     1/23/14
             playmovieaction  extract
             if( j<f2) 
             playmovieaction step=1
             endif
             imagetransform rgb2gray M_movieframe
             duplicate/O M_rgb2gray image1 
		Duplicate/O image1 OrigImage1
		ImageRotate/A=(AngleOfRotation) OrigImage1
		Duplicate/O M_RotatedImage RotatedOrig1
		DoWindow/K RotatedImage
		NewImage/K=1/N=RotatedImage/S=1 RotatedOrig1
	
	//	ModifyImage RotatedOrig1 ctab= {ASTR[6],Astr[7],Grays,0} // 9/7/10 CC, below modifies this line for changing image color *********************************
		ModifyImage RotatedOrig1 ctab= {ASTR[6],Astr[7],Grays,0}
	
		ModifyGraph width=ASTR[3]*ASTR[5]*ImScale,height=ASTR[4]*ASTR[5]*ImScale
		 
		lengthandangle=0; lengthandangle_c2 = 0; bluelengthandangle=0; k=0; CellOccupado=0; BlueOccupado=0
		sumlen = 0; sumlenerr = 0; sumlenblue = 0; sumlenerrblue=0;
		energy = 0; energy2=0;energyerr = 0 ; energyblue = 0; energyerrblue = 0  ; energy_c2 = 0// 8/26/06
		makeshiftX=0; makeshiftY=0

		for (i=0; i<=(TotalNumPosts-1); i+=1)
					
			// 7/10/06 : fixed indexing of fitresultsbase (no 3rd index)
			if ((Cellmatrix[i]== kCell) || (Cellmatrix[i]== kEmpty)|| (Cellmatrix[i] == kwirepost)) //changed 1/7/09 CMK to include guide posts
	
			    if (shiftflag==0)  // grideachframe
			        dx = FRdisp[i][0][j]
			        dy = FRdisp[i][1][j]
			        
			        length =sqrt(dx^2+dy^2)
			   
			   		    corrwave[i][0][j] = length //for use in correlation calculations -- just a copy of the energy
			        
			        if(j==0) 
			        	lengthchange0[i] = length
			        elseif(j==(f2-f1+1)/2)
			        	lengthchangemid[i] = length
			        endif
			        if(j< (f2-f1+1)/2)
			        	
			        	lengthchangewave[i] =abs( length - lengthchange0[i]  )
			        elseif(j>= (f2-f1+1)/2)
			            
			        	lengthchangewave[i] = abs(length - lengthchange0[i])
				  endif
//					if(j==(f2-f1+1)/2)		//this is written assuming that field is turned on halfway through experiment CMK
//						energychangewave[j] +=length^2
//					endif
			        
                            if (length == 0)
					    printf, " Warning: Length = 0 ;  post = %d nimage = %d\r", i,j
					    lengtherr[i] = sqrt(FREdisp[i][0][j]^2 + FREdisp[i][1][j]^2 )
			        else
                                lengtherr[i] = sqrt(FREdisp[i][0][j]^2*dx^2 + FREdisp[i][1][j]^2 * dy^2)
				      lengtherr[i] /= length    
				 endif
		          elseif (shiftflag == 1)  // use old shift method  
				       dx = FRS[i][0][j]-fitresultsbase[i][0][f1]   // note f1 here
					dy = FRS[i][1][j]-fitresultsbase[i][1][f1]
					length =sqrt(dx^2+dy^2)
					if (length == 0)
					    printf, " Warning: Length = 0 ;  post = %d nimage = %d\r", i,j
					     lengtherr[i] =   sqrt(FRE[i][0][j]^2 + FRE[i][1][j]^2 )
					    // lengtherr[i] =   sqrt(FRE[i][kFX0][j]^2 + fiterrorsbase[i][kFX0][f1]^2 + FRE[i][kFY0][j]^2 +  fiterrorsbase[i][kFY0][f1]^2)
					else
			                  lengtherr[i] =   sqrt((FRE[i][0][j]^2 )*dx^2 + (FRE[i][1][j]^2)*dy^2)
			                  //lengtherr[i] =   sqrt((FRE[i][kFX0][j]^2 + fiterrorsbase[i][kFX0][f1]^2)*dx^2 + (FRE[i][kFY0][j]^2 + fiterrorsbase[i][kFY0][f1]^2)*dy^2)
				           lengtherr[i] /= length    
					endif
		          endif
			   if (Cellmatrix[i]==kCell)  // || Cellmatrix[i] == kwirepost)      commented out for making cell pair  6/1/2016

					//makeshiftX[i]=FRS[i][kFX0][j]; makeshiftY[i]=FRS[i][kFY0][j]
					makeshiftX[i]=FRS[i][0][j]+ShiftWave[j][12]; makeshiftY[i]=FRS[i][1][j]+ShiftWave[j][13]  //Modified Jul31, 2006
					sumlen += length  // changed 8/28/06
					sumlenerr += lengtherr[i]^2       // these ''errors" are actually error^2 at this point
					energy += length^2
					energy2 += (lengthchangewave[i])^2
				
//					if(j==0)
//						energychangewave[j] +=length^2
//					endif
//					if(j==(f2-f1+1)/2)		//this is written in assuming that field is turned on halfway through experiment CMK
//						energychangewave[j] +=length^2
//					endif
					//energyerr += 4*length^2 * lengtherr[i]^2 //CMK 7-10-09 length should not be multiplied
					energyerr += 2* lengtherr[i]*length
					
				
					
					lengthandangle[i][0]  = MFact * length   * kspring * lattice/Astr[5]// scale for quiver plot
					lengthandangle[i][1]=-angle(dx, dy)
					
					totalforce[j][0] += dx * kspring * lattice/Astr[5]
					totalforce[j][1] += dy * kspring * lattice/Astr[5]
					
					CellOccupado=CellOccupado+1    
			   //add for cell pairs computation   //  6/1/2016
			   elseif (Cellmatrix[i] == kwirepost)
			   
					makeshiftX_c2[i]=FRS[i][0][j]+ShiftWave[j][12]; makeshiftY_c2[i]=FRS[i][1][j]+ShiftWave[j][13]  //Modified Jul31, 2006
					sumlen_c2 += length  // changed 8/28/06
					sumlenerr_c2 += lengtherr[i]^2       // these ''errors" are actually error^2 at this point
					energy_c2 += length^2
					energy2 += (lengthchangewave[i])^2
				
//					if(j==0)
//						energychangewave[j] +=length^2
//					endif
//					if(j==(f2-f1+1)/2)		//this is written in assuming that field is turned on halfway through experiment CMK
//						energychangewave[j] +=length^2
//					endif
					//energyerr += 4*length^2 * lengtherr[i]^2 //CMK 7-10-09 length should not be multiplied
					energyerr_c2 += 2* lengtherr[i]*length
					
				
					
					lengthandangle_c2[i][0]  = MFact * length   * kspring * lattice/Astr[5]// scale for quiver plot
					lengthandangle_c2[i][1]=-angle(dx, dy)
					
					totalforce_c2[j][0] += dx * kspring * lattice/Astr[5]
					totalforce_c2[j][1] += dy * kspring * lattice/Astr[5]
					
					CellOccupado=CellOccupado+1    			            		
					
			   elseif (Cellmatrix[i]==kEmpty)
					
					blueshiftX[i]=FRS[i][0][j]+ShiftWave[j][12]; blueshiftY[i]=FRS[i][1][j]+ShiftWave[j][13]
					sumlenblue += length  // changed 8/28/06
					bluelengtherr[i] = lengtherr[i]
					lengtherr[i] = 0   // 2/16/12 bookkeeping to separate cell and non cell errors
					
	// 2/15/12 - bug fix				sumlenerrblue += bluelengtherr[i]^2       // these ''errors" are actually error^2 at this point
					sumlenerrblue += bluelengtherr[i]^2       // these ''errors" are actually error^2 at this point
					energyblue += length^2
// 2/15/12 bux fix - Craig and Dan					energyerrblue += 4*length^2 * bluelengtherr[i]^2
					energyerrblue += 2*length * bluelengtherr[i]

					bluelengthandangle[i][0]  = MFact * length  * kspring * lattice/Astr[5]  // scale for quiver plot
					bluelengthandangle[i][1]=-angle(dx, dy)
					BlueOccupado+=1	
				endif
			endif
		endfor
	//**********2/15/12*********
	duplicate/o lengtherr lengtherrF
	duplicate/o bluelengtherr bluelengtherrF
	lengtherrF *=MFact * kspring * lattice/Astr[5]
	BluelengtherrF *=MFact * kspring * lattice/Astr[5]
	//**********************
//8/28/06 TOTAL DISPLACEMENT NOT AVG

	
	//Print "Energy = ", energy
       AvDisp[j] = sumlen
       AvErr[j] = sqrt(sumlenerr)
       //TotalE[j] = 0.5 * kspring * energy * (lattice /ASTr[5])^2   // 8/28/06 "U = 0.5 K x^2"  30 = nN/um , 9 = lattice cte in um, ASTr[5] = lattice cte in pixel  :  gives energy in attoJoules "energy ==  length^2
       //showing avg force
        TotalE[j] = kspring * sqrt(energy) * (lattice /ASTr[5])
        TotalE_c2[j] = kspring * sqrt(energy_c2) * (lattice /ASTr[5])         //added 6/1/2016
       TotalEerr[j] =  0.5 * kspring  *(lattice /ASTr[5])* (energyerr)
       TotalF[j] = 0.5 * kspring * energy * (lattice /ASTr[5])^2                   //added 4/16/18    
       TotalE2[j] = 0.5 * kspring * energy2 * (lattice /ASTr[5])
       
             
       
//       if(j<(f2-f1+1)/2)
//       	TotalE2[j] = 0.5 * kspring * (energy-energychangewave[0]) * (lattice /ASTr[5])^2 
//       elseif(j>=(f2-f1+1)/2)
//       	TotalE2[j] = 0.5 * kspring * (energy-energychangewave[(f2-f1+1)/2]) * (lattice /ASTr[5])^2 
//       endif
       //  Add scale bar
		makeshiftX[0]=ASTR[0]+AX+ShiftWave[j][12]; makeshiftY[0]=ASTR[1]+AY+ShiftWave[j][13]
		
		lengthandangle[0][0] = MFACT*Fscale*kspring*lattice/Astr[5] //CMK what is this number?-- not important that's just the scale bar...
		lengthandangle[0][1] = 0 // 
       
	AppendToGraph/T  makeshiftY vs makeshiftX

if (timeflag == 1)
     TextBox/C/N=text0/A=RT/X=1.00/Y=1.00 num2str(vtime(j))	
endif
	
	//TextBox/C/N=text1/A=LT/B=1/G=(65200,0,0)/X=1.20/Y=1.50/T=20 num2str(Fscale)+" nN"//9/7/10 CC, below modifies this line for changing scale bar *************************
	TextBox/C/N=text1/A=LT/B=1/G=(65200,65200,65200)/X=1.20/Y=1.50/T=20 num2str(Fscale)+" nN"
	//ModifyGraph width={Aspect,Astr[3]/Astr[4]}
//	ModifyGraph arrowMarker(makeshiftY)={lengthandangle,5,6,2,0}

	//ModifyGraph arrowMarker(makeshiftY)={lengthandangle,3+Afact,4+Afact,2,0} // 9/7/10 CC, below modifies this line for changing arrow appearances*************************
	ModifyGraph arrowMarker(makeshiftY)={lengthandangle,Afact,2+Afact,2,0}
	//ModifyGraph rgb(makeshiftY)=(65280,0,0)//  9/7/10 CC, below modifies this line for changing arrow appearances
	ModifyGraph rgb(makeshiftY)=(65535,0,0)
	
	ModifyGraph mode=3,marker=19,msize=1
	ModifyGraph tick=1,nticks=10,noLabel=0
	
	//added 6/1/2016
	AppendToGraph/T  makeshiftY_c2 vs makeshiftX_c2
	ModifyGraph arrowMarker(makeshiftY_c2)={lengthandangle_c2,Afact,2+Afact,2,0}
	//ModifyGraph rgb(makeshiftY)=(65280,0,0)//  9/7/10 CC, below modifies this line for changing arrow appearances
	ModifyGraph rgb(makeshiftY_c2)=(0,65535,0)
	
	ModifyGraph mode=3,marker=19,msize=1
	ModifyGraph tick=1,nticks=10,noLabel=0	
	//added 6/1/2016
	
	
//	BlueDisp[j]=sum(bluelengthandangle, 0,TotalNumPosts-1)/(MFact*BlueOccupado)
//8/28/06 TOTAL DISPLACEMENT NOT AVG
       BlueDisp[j] = sumlenblue
       BlueErr[j] = sqrt(sumlenerrblue)
//       BlueTotalE[j] = 0.5 * kspring * energyblue * (lattice /ASTr[5])^2   // 8/28/06 "U = 0.5 K x^2"  10.7 = nN/um , 4 = lattice cte in um, ASTr[5] = lattice cte in pixel  :  gives energy in attoJoules 
 // 2/15/12 bug fix - Craig and Dan      BlueTotalEerr[j] =  0.5 * kspring  * (lattice /ASTr[5])^2 * sqrt(energyerrblue)
       //showing total force  5-5-2016
       BlueTotalE[j] =   kspring * sqrt(energyblue) * (lattice /ASTr[5])
       
       BlueTotalEerr[j] =  0.5 * kspring  * (lattice /ASTr[5])^2 * energyerrblue

//put all this back in, this is for making nice looking plots 9/21/09 CMK

//Craig flags added:
if(bArrowflag==1)
	AppendToGraph/T blueshiftY vs blueshiftX 
	ModifyGraph arrowMarker(blueshiftY)={bluelengthandangle,5,6,2,0}
	ModifyGraph arrowMarker(blueshiftY)={bluelengthandangle,Afact,2+Afact,2,0}
	ModifyGraph rgb(blueshiftY)=(0,0,65280)
	ModifyGraph mode=3,marker=19,msize=1
	ModifyGraph tick=0,nticks=20,noLabel=2
endif
If(interplineflag==1)
	Make/o/n=(totalnumposts) Basex,Basey
	Basex = fitresultsbase[p][0][j]
	Basey = fitresultsbase[p][1][j]
	Appendtograph/T basey vs basex
	Modifygraph mode(basey)=3,marker(basey)=19,msize(basey)=2,rgb(basey)=(0,65535,0)
endif


	SetAxis/R left (ASTR[1]+ASTR[5]*(ASTR[4]-0.5)+ShiftWave[j][13]),(ASTR[1]-ASTR[5]*0.5+ShiftWave[j][13]);  // Fixed ROI, no shift!!!
	SetAxis top (ASTR[0]-ASTR[5]*0.5+ShiftWave[j][12]),(ASTR[0]+ASTR[5]*(ASTR[3]-0.5)+ShiftWave[j][12])

	
	//adding indices CMk 4/23/09
	if(indexflag==1)
		showindices(2,j)
 	endif
 	
 	//1/25/2014  Yu Shi
 	if (mod(j,freq)==0) 
	      AddMovieFrame
	endif

	endfor //j

     CloseMovie

    DoWindow/K AverageDisp
        AvDisp /= CellOccupado
        BlueDisp /= BlueOccupado
    if (timeflag == 1)
        display/K=1/N=AveragDisp AvDisp vs vTime
        appendtograph BlueDisp vs vTime
        Label Bottom "Time (min)"
    else

        display/K=1/N=AverageDisp  AvDisp,BlueDisp 
    endif
    ModifyGraph rgb(BlueDisp)=(0,0,65280)
    ModifyGraph mode(AvDisp)=4,marker(AvDisp)=19;DelayUpdate
    ErrorBars AvDisp Y,wave=(AvErr,AvErr)
    ModifyGraph mode(BlueDisp)=4,marker(BlueDisp)=19;DelayUpdate
    ErrorBars BlueDisp Y,wave=(BlueErr,BlueErr)
    TextBox/C/N=text1/F=0/M/A=RB ("Occupado "+num2str(CellOccupado)+ " Blue "+ num2str(BlueOccupado))
    Label Left "Avg Disp (pix)"
     ///Modified Jul. 27 2006
    SaveGraphCopy/O/W=AverageDisp/P=ToPractice as "AvDisp"+RootString+".pxp"  //Saving the graphs added  July 27, 2006
    SavePICT/O/P=ToPractice/E=-7/B=144 as  "AvDisp"+RootString+".tif"
//**********Making force calculation "totalF" will actually be an averge to stay consistent with "totalE"****
  variable TotalF1, totalF2
//  TotalF1 = kspring*avgdisp*celloccupado*lattice/Astr[5]*mfact
  duplicate/o/r=[][0] lengthandangle, justlength
  TotalF2 = sum(justlength)/lengthandangle[0][0]*10
 // print totalF
  print totalf2
  //***************************************
// 8/28/06 for Total Energy Plot
DoWindow/K TotalEnergy
     BlueTotalE /=BlueOccupado
    // print "Cell Occupado = , CellOccupado
   //print "Total E =", TotalE
//     TotalE /= CellOccupado   changed during 5/2/2016, changed to total energy
//     TotalE2  /= CellOccupado  changed during 5/2/2016, changed to total energy
     TotalEerr/=CellOccupado
     BlueTotalEerr/=BlueOccupado   // 2/15/12  bug fix - Dan and Craig (this line was absent)
     
     wavestats/q totalE
     avgE = V_avg
     wavestats/q totalF
     avgF = V_avg
     
     
     
if (timeflag == 1)
wave timewave
    
  
    display/K=1/N=TotalEnergy TotalE vs vTime 
    appendtograph BlueTotalE vs vTime
    Label Bottom "Time (min)"
    
else
    display/K=1/N=TotalEnergyChanges TotalE2 
    sprintf s10, "%s",RootString
    TextBox/C/B=1/N=text0/A=LT s10
    
    
    display/K=1/N=TotalEnergy  TotalE,BlueTotalE
    TextBox/C/B=1/N=text0/A=LT s10
endif
//Label Left "Avg E (fJ)"
Label Left "Avg F (nN)"
ModifyGraph rgb(BlueTotalE)=(0,0,65280)
ModifyGraph mode(TotalE)=4,marker(TotalE)=19, msize(TotalE)=1;DelayUpdate
ErrorBars TotalE Y,wave=(TotalEerr,TotalEerr)
ModifyGraph mode(BlueTotalE)=4,marker(BlueTotalE)=19, msize(BlueTotalE)=1; DelayUpdate
ErrorBars BlueTotalE Y,wave=(BlueTotalEerr,BlueTotalEerr)
//TextBox/C/N=text1/F=0/M/A=RB ("Occupado "+num2str(CellOccupado)+ " Blue "+ num2str(BlueOccupado))
TextBox/C/N=text1/F=0/M/A=RB ("TotalF "+num2str(avgF)+" Occupado "+num2str(CellOccupado)+ " Blue "+ num2str(BlueOccupado))
///Modified Jul. 27 2006
SaveGraphCopy/O/W=TotalEnergy/P=ToPractice as "TotalE"+RootString+".pxp"  //Saving the graphs added  July 27, 2006
SavePICT/O/P=ToPractice/E=-7/B=144 as  "TotalE"+RootString+".tif"


// *********Craig adding force vector summing*****************************
ForceSum()
//****************************************************************
                   
// pre-7/10/06 code			
//				if (Cellmatrix[i]==kCell)
//					makeshiftX[i]=FRS[i][kFX0][j]+ShiftWave[j][12]; makeshiftY[i]=FRS[i][kFY0][j]+ShiftWave[j][13]
//					lengthandangle[i][0]=MFact*sqrt((FRS[i][kFX0][j]-fitresultsbase[i][kFX0][j])^2+(FRS[i][kFY0][j]-fitresultsbase[i][kFY0][j] )^2)
//					lengthandangle[i][1]=-angle((FRS[i][kFX0][j]-fitresultsbase[i][kFX0 ][j]), (FRS[i][kFY0][j]-fitresultsbase[i][kFY0][j]))  
//					CellOccupado=CellOccupado+1    
//				elseif (Cellmatrix[i]==0)
//					blueshiftX[i]=FRS[i][kFX0][j]+ShiftWave[j][12]; blueshiftY[i]=FRS[i][kFY0][j]+ShiftWave[j][13]
//					bluelengthandangle[i][0]=MFact*sqrt((FRS[i][kFX0][j]-fitresultsbase[i][kFX0 ][j])^2+(FRS[i][kFY0][j]-fitresultsbase[i][kFY0][j] )^2)
//					bluelengthandangle[i][1]=-angle((FRS[i][kFX0][j]-fitresultsbase[i][kFX0 ][j]), (FRS[i][kFY0][j]-fitresultsbase[i][kFY0][j] ))   
//					BlueOccupado+=1	
//				endif
//			endif
//		endfor
//	AvDisp[j]=sum(lengthandangle, 0,TotalNumPosts-1)/(MFact*CellOccupado)
//	AppendToGraph/T  makeshiftY vs makeshiftX
//	TextBox/C/N=text0/A=RT/X=1.00/Y=1.00 num2str(j)
//	//ModifyGraph width={Aspect,Astr[3]/Astr[4]}
//	ModifyGraph arrowMarker(makeshiftY)={lengthandangle,5,6,2,0}
//	ModifyGraph rgb(makeshiftY)=(65280,0,0)
//	ModifyGraph mode=3,marker=19,msize=1
//	ModifyGraph tick=3,nticks=5,noLabel=2
//	
//	BlueDisp[j]=sum(bluelengthandangle, 0,TotalNumPosts-1)/(MFact*BlueOccupado)
//	AppendToGraph/T blueshiftY vs blueshiftX
//	ModifyGraph arrowMarker(blueshiftY)={bluelengthandangle,5,6,2,0}
//	ModifyGraph rgb(blueshiftY)=(0,0,65280)
//	ModifyGraph mode=3,marker=19,msize=1
//	ModifyGraph tick=3,nticks=5,noLabel=2
		
//	SetAxis/R left (ASTR[1]+ASTR[5]*(ASTR[4]-0.5)+ShiftWave[j][13]),(ASTR[1]-ASTR[5]*0.5+ShiftWave[j][13]);  // Fixed ROI, no shift!!!
//	SetAxis top (ASTR[0]-ASTR[5]*0.5+ShiftWave[j][12]),(ASTR[0]+ASTR[5]*(ASTR[3]-0.5)+ShiftWave[j][12])
//	AddMovieFrame
//	endfor //j
//CloseMovie
//display/K=1 AvDisp,BlueDisp; ModifyGraph rgb(BlueDisp)=(0,0,65280); 
//TextBox/C/N=text1/F=0/M/A=RB ("Occupado "+num2str(CellOccupado)+ " Blue "+ num2str(BlueOccupado))

end

///------------------------ View Jumps


function ViewJump(ctrlName): ButtonControl 
string CtrlName   
	variable f1j=39,f2j=40
	Prompt f1j, "Enter First File Index"     
	Prompt f2j, "Enter Last File Index"        
      	DoPrompt "Indexes", f1j, f2j  
	if (V_Flag==1)
		return 0
	endif          

SVAR FitRes_String, FitErr_String,ShiftString, AttributeString,MaskString,FitBottom,FitErrBottom

Duplicate/O $AttributeString ASTr
Duplicate/O $MaskString maskwave    
duplicate/O $FitRes_String  fitres_plotAlot
duplicate/O $FitErr_String  fiterr_plotAlot    // 6/02/06
// these next three lines for looking at zero positions
//wave FB_cell28, FBerr_cell28
//duplicate/O FB_cell28 fitres_plotAlot
//duplicate/O FBerr_cell28 fiterr_plotAlot    // 6/02/06
variable f1,f2
f1 = dimoffset(fitres_plotAlot,2)
f2 = dimsize(fitres_plotAlot,2) + dimoffset(fitres_plotAlot,2) -1
    calcdispGridEach(f1,f2,ASTr[3],ASTr[4],fitres_plotAlot,fiterr_plotAlot,$FitBottom,$FitErrBottom)
    //(f1,f2,ASTr[3],ASTr[4],fitres_plotAlot,fiterr_plotAlot,$ShiftString) //CMK removing 10/24/08
 
variable Nx = dimsize($MaskString,0)
variable Ny = dimsize($MaskString,1)

variable i,j
make/o/n=(Nx,Ny) testjumpx,testjumpy,testjumpr,testjumprelx,testjumprely, testjumprel
testjumpr = 0; testjumpx = 0; testjumpy = 0; testjumprelx = 0; testjumprely = 0
for (i=0; i < Nx; i+=1)
    for (j=0; j< Ny; j +=1)
        if (maskwave[i][j] !=kIgnore)
        testjumpx[i][j] = fitres_plotalot[i*Ny+j][kfX0][f2j]-fitres_plotalot[i*Ny+j][kfX0][f1j]
        testjumpy[i][j] = fitres_plotalot[i*Ny+j][kfY0][f2j]-fitres_plotalot[i*Ny+j][kfY0][f1j]
        testjumprelx[i][j] = testjumpx[i][j]/fitres_plotalot[i*Ny+j][kfX0][f1j]
        testjumprely[i][j] = testjumpy[i][j]/fitres_plotalot[i*Ny+j][kfY0][f1j]
        endif
     
    endfor
endfor
testjumpr = sqrt(testjumpx^2 + testjumpy^2)
testjumprel = sqrt(testjumprelx^2 + testjumprely^2)
// find cell only
variable ncp = 0
for (i=0; i < Nx; i+=1)
    for (j=0; j< Ny; j +=1)
        if (maskwave[i][j] == kCell)
            ncp +=1
        endif
    endfor
endfor
make/O/N=(ncp) jumprc,jumprcrel, cellx,celly
variable ic = 0
for (i=0; i < Nx; i+=1)
    for (j=0; j< Ny; j +=1)
        if (maskwave[i][j] == kCell)
            jumprc[ic] = testjumpr[i][j]
            jumprcrel[ic] = testjumprel[i][j]
            cellx[ic] = i
            celly[ic] = j
            ic +=1
        endif
    endfor
endfor


DoWindow/K Xjump
NewImage/N=Xjump testjumpx
appendtograph/t celly vs cellx
ModifyGraph mode=3

DoWindow/K Yjump
NewImage/N=Yjump testjumpy
appendtograph/t celly vs cellx
ModifyGraph mode=3

DoWindow/K Rjump
NewImage/N=Rjump testjumpr
appendtograph/t celly vs cellx
ModifyGraph mode=3


// histograms
variable maxv = wavemax(testjumpr)
Make/N=100/O W_HistR;DelayUpdate
Histogram/B={0,maxv/100,100} testjumpr,W_HistR
Make/N=100/O W_HistRC;DelayUpdate
Histogram/B={0,maxv/100,100} jumprc,W_HistRC
DoWindow/K HistR
display/N=HistR w_histr
appendtograph w_histrc
ModifyGraph rgb(W_HistRC)=(65535,0,0)
ModifyGraph rgb(W_HistR)=(0,0,65535)
end

///------------------------ Avg Jumps

// 4/20/07
function AvgJump(ctrlName): ButtonControl 
string CtrlName   
	variable f1L=35,f1R=39
	variable f2L = 76, f2R=80
	Prompt f1L, "First Window Left"     
	Prompt f1R, "First Window Right"     
	Prompt f2L, "Second Window Left"     
	Prompt f2R, "Second Window Right"     
      	DoPrompt "Indexes", f1L, F1R,F2L,F2R
	if (V_Flag==1)
		return 0
	endif          

SVAR FitRes_String, FitErr_String,ShiftString, AttributeString,MaskString,FitBottom,FitErrBottom

Duplicate/O $AttributeString ASTr
Duplicate/O $MaskString maskwave    
duplicate/O $FitRes_String  fitres_plotAlot
duplicate/O $FitErr_String  fiterr_plotAlot    // 6/02/06
// these next three lines for looking at zero positions
//wave FB_cell28, FBerr_cell28
//duplicate/O FB_cell28 fitres_plotAlot
//duplicate/O FBerr_cell28 fiterr_plotAlot    // 6/02/06
variable f1,f2
f1 = dimoffset(fitres_plotAlot,2)
f2 = dimsize(fitres_plotAlot,2) + dimoffset(fitres_plotAlot,2) -1
    calcdispGridEach(f1,f2,ASTr[3],ASTr[4],fitres_plotAlot,fiterr_plotAlot,$FitBottom,$FitErrBottom)
    //ApplyRotation(f1,f2,ASTr[3],ASTr[4],fitres_plotAlot,fiterr_plotAlot,$ShiftString) //CMK removing 10/24/08

variable Nx = dimsize($MaskString,0)
variable Ny = dimsize($MaskString,1)

variable i,j, iw, wx1,wy1,wrx1,wry1,wx2,wy2,wrx2,wry2
make/o/n=(Nx,Ny) testjumpx,testjumpy,testjumpr,testjumprelx,testjumprely, testjumprel
testjumpr = 0; testjumpx = 0; testjumpy = 0; testjumprelx = 0; testjumprely = 0
variable Nwin1 = f1R - F1L + 1
variable Nwin2 = f2R - F2L +1
for (i=0; i < Nx; i+=1)
    for (j=0; j< Ny; j +=1)
        if (maskwave[i][j] !=kIgnore)
            wx1=0; wy1=0; wx2=0; wy2=0
            for (iw = f1L; iw <= f1R; iw +=1)  // loop over first window
                wx1 += fitres_plotalot[i*Ny+j][kfX0][iw]
                wy1 += fitres_plotalot[i*Ny+j][kfY0][iw]
            endfor
            for (iw = f2L; iw <= f2R; iw +=1)   // loop over second window
                wx2 += fitres_plotalot[i*Ny+j][kfX0][iw]
                wy2 += fitres_plotalot[i*Ny+j][kfY0][iw]
            endfor
            
            testjumpx[i][j] = wx2/Nwin2 - wx1/Nwin1
            testjumpy[i][j] = wy2/Nwin2 - wy1/Nwin1
 //           testjumprelx[i][j] = testjumpx[i][j]/fitres_plotalot[i*Ny+j][kfX0][f1j]
 //           testjumprely[i][j] = testjumpy[i][j]/fitres_plotalot[i*Ny+j][kfY0][f1j]
        endif
     
    endfor
endfor
testjumpr = sqrt(testjumpx^2 + testjumpy^2)
//testjumprel = sqrt(testjumprelx^2 + testjumprely^2)
// find cell only
variable ncp = 0
for (i=0; i < Nx; i+=1)
    for (j=0; j< Ny; j +=1)
        if (maskwave[i][j] == kCell)
            ncp +=1
        endif
    endfor
endfor
make/O/N=(ncp) jumprc,jumprcrel, cellx,celly
variable ic = 0
for (i=0; i < Nx; i+=1)
    for (j=0; j< Ny; j +=1)
        if (maskwave[i][j] == kCell)
            jumprc[ic] = testjumpr[i][j]
//            jumprcrel[ic] = testjumprel[i][j]
            cellx[ic] = i
            celly[ic] = j
            ic +=1
        endif
    endfor
endfor

//save info to .ibw for rescaling in another program 

Save/O/P=ToPractice testjumpr as "testjumpr"+FitRes_String+".ibw"
Save/O/P=ToPractice cellx as "cellx"+FitRes_String+".ibw"
Save/O/P=ToPractice celly as "celly"+FitRes_String+".ibw"

DoWindow/K Xjump
NewImage/N=Xjump testjumpx
appendtograph/t celly vs cellx
ModifyGraph mode=3


DoWindow/K Yjump
NewImage/N=Yjump testjumpy
appendtograph/t celly vs cellx
ModifyGraph mode=3

DoWindow/K Rjump
NewImage/N=Rjump testjumpr
appendtograph/t celly vs cellx
ModifyGraph mode=3


// histograms
variable maxv = wavemax(testjumpr)
Make/N=100/O W_HistR;DelayUpdate
Histogram/B={0,maxv/100,100} testjumpr,W_HistR
Make/N=100/O W_HistRC;DelayUpdate
Histogram/B={0,maxv/100,100} jumprc,W_HistRC
DoWindow/K HistR
display/N=HistR w_histr
appendtograph w_histrc
ModifyGraph rgb(W_HistRC)=(65535,0,0)
ModifyGraph rgb(W_HistR)=(0,0,65535)

//CL add- put plots directly into notebook

DoWindow/K NB30
NewNotebook/F=1/N=NBAvgJump

string Nbinfo
sprintf Nbinfo, "%s- avg (%02d ,%02d), (%02d ,%02d)\r", Fitres_string,F1L, F1R, F2L, F2R
Notebook NBAvgJump text=Nbinfo
Notebook NBAvgJump scaling= {90,90}, picture = {Xjump, -1,1}
Notebook NBAvgJump text = "\r"
Notebook NBAvgJump scaling= {90,90}, picture = {Yjump, -1,1}
Notebook NBAvgJump text = "\r"
Notebook NBAvgJump scaling= {90,90}, picture = {Rjump, -1,1}
Notebook NBAvgJump text = "\r"
Notebook NBAvgJump scaling= {90,90}, picture = {HistR, -1,1}


end


function runcorrelation(ctrlName) : ButtonControl
string ctrlname

wave corrwave
variable f1, f2
Prompt f1, "Enter First File Index"     
Prompt f2, "Enter Last File Index"        
DoPrompt "Indexes", f1, f2


correlation(corrwave, f1, f2)

end




//********For drawing outline of cell mask************************************
  	Function Drawbutton(ctrlname):buttoncontrol
	  	string ctrlname
	  	Drawfunct()
	END
	  	
	 Function drawfunct()
	  	showtools/w=rotatedimage
		SetDrawLayer/w=rotatedimage ProgFront
		SetDrawEnv/w=rotatedimage fillpat=0, XCOORD=TOP, YCOORD=LEFT
		Button Drawit, title="Done Drawing", proc = acceptroi
	 END
		
	 Function acceptroi(ctrlname): buttoncontrol
		string ctrlname
		
		make/o/n=1 nodrawflag
		imagegenerateROImask/w=rotatedimage rotatedorig1
		if(v_flag==1)
			wave M_roimask
			wave rotatedorig1
			rotatedorig1 = rotatedorig1*M_roimask
			nodrawflag[0]=0
		else
			nodrawflag[0] = 1
	      endif
	      dowindow/k drawroipanel
	  	dowindow/k rotatedimage

	 END
//*******************************************************************	 
	 
// *********Craig adding force vector summing*****************************
// This runs automatically at the end of MakeItQuiverG, and also during force summing auto
Function ForceSum()

variable i,Fytot=0,Fxtot=0,Angletot=0,FmagXtot=0,FmagYtot=0,fractionleftx,fractionlefty,Frtot,fractionleftr,FmagRtot=0,totalnumposts
wave linmaskval, astr,lengthandangle

TotalNumPosts=ASTR[3]*ASTR[4]
 
for (i=1; i<=(TotalNumPosts); i+=1)  //Ignore post 0, which is used for the scale bar
	if(linmaskval[i]==1)
 		Fytot = Fytot+(lengthandangle[i][0]*sin(lengthandangle[i][1]))
 		Fxtot = Fxtot+(lengthandangle[i][0]*cos(lengthandangle[i][1]))   
 		FmagXtot = FmagXtot + abs(lengthandangle[i][0]*cos(lengthandangle[i][1]))   
 		FmagYtot = FmagYtot + abs(lengthandangle[i][0]*sin(lengthandangle[i][1])) 	
 		FmagRtot = FmagRtot + Lengthandangle[i][0]             
 	endif
 	Frtot = sqrt(fxtot^2+Fytot^2)
 	fractionleftx = abs(Fxtot/FmagXtot)
 	fractionlefty = abs(Fytot/FmagYtot)
 	fractionleftr = Frtot/FmagRtot
 endfor
angletot = angle(fxtot,fytot)
//Create output wave storing all of this information
Make/o/n=(1,7,0,0) ForceSumStats
ForceSumStats[0][0] = Fxtot 
ForceSumStats[0][1] = fractionleftx
ForceSumStats[0][2] = Fytot 
ForceSumStats[0][3] = Fractionlefty
ForceSumStats[0][4] = Frtot
ForceSumStats[0][5] = Fractionleftr
ForceSumStats[0][6] = angletot
print "FxTotal = ",Fxtot,"(",fractionleftx,")", "FyTotal = ", Fytot, "(",fractionlefty,")","Fr = ",Frtot,"(",fractionleftr,")","Resultant angle = ",angletot

END
//****************************************************************


//*********************for filtering low deflection posts*******************************
function ForceFilter(ctrlName): ButtonControl
    string CtrlName
    
NewPanel/K=1/N=ForceFiltering/W=(100,550,390,710)
Button BackupMask,pos={35,03},size={225,20},proc=maskbackup,title="Backup Mask Vals"
Button ShowHist,pos={35,33},size={225,20},proc=showhist,title="See Histogram"
Button filterforces,pos={35,63},size={225,20},proc=filterforces,title="Filter Forces"
Button RestoreMask,pos={35,93},size={225,20},proc=restoremaskval,title="Restore Mask"    
Button AutoFind, pos={35,123}, size={225,20},proc=FindForceCutoff, title="Autofind Cutoff"
END
//********************************************************
Function filterforces(ctrlname):buttoncontrol
string ctrlname

variable forcecutoff
variable i,j
wave linmaskval
wave lengthandangle

Prompt forcecutoff, "Enter value"
DoPrompt "Force Filter Panel", forcecutoff

duplicate/o/r=[0,*][0][0] lengthandangle forfilter

for(i=0;i<(dimsize(forfilter,0));i+=1)
	if(forfilter[i]<forcecutoff && linmaskval[i] !=0)
		linmaskval[i] = 3 //ignore it
	endif
endfor
print "Force cutoff =", forcecutoff

END
//******************backup original mask*************************
Function maskbackup(ctrlname):buttoncontrol
string ctrlname

wave linmaskval
duplicate/o linmaskval linmaskvalbackup

END
//*****************for restoring linmaskval after a force filter**********
Function restoremaskval(ctrlname):buttoncontrol
string ctrlname

wave linmaskval,linmaskvalbackup

linmaskval = linmaskvalbackup

END
//***********Show Histogram of vectors*****************************
Function showhist(ctrlname):buttoncontrol
string ctrlname

wave lengthandangle
variable numbins
	Prompt numbins,"How many bins?"
	DoPrompt "Set Histogram Bin Number",numbins
		
	duplicate/o/r=[0,*][0][0] lengthandangle forfilter
	forfilter[0] = 0 //do not include scale bar
	Make/N=(numbins)/O forfilter_Hist
	Histogram/C/B=1 forfilter,forfilter_Hist
	display forfilter_hist
	ModifyGraph mode=8,mrkThick=2
	SetAxis left 0,500
	
END	
//************************************************************************************************************************
Function FindForceCutoff(ctrlname):buttoncontrol
string ctrlname


wave lengthandangle,linmaskval, linmaskvalbackup

//The below needs to not be a user input and some automated value
variable cutoffmax,stepsize,i,j, numcutoffs,firsttimeflag=1
Prompt cutoffmax,"Enter cutoff max"
Prompt stepsize, "Enter step size"
DoPrompt "Set cutoff max and step size", cutoffmax, stepsize

Duplicate/o/r=[0,*][0][0] lengthandangle forfilter
forfilter[0]=0 // wipe out scale bar

numcutoffs = (ceil(cutoffmax/stepsize)+1) // +1 because we want to see zero
Make/O/n=(numcutoffs) CutOffValues
For(i=0; i< numcutoffs; i+=1) //Fill wave with cutoff values to test
	Cutoffvalues[i] = i*stepsize
endfor	

For(j=0; j< numcutoffs; j+=1)	
	for(i=0;i<(dimsize(forfilter,0));i+=1)
		if(forfilter[i]< Cutoffvalues[j] && linmaskval[i] ==1)
			linmaskval[i] = 3 //ignore it
		endif
	endfor
	//Sum the forces with filter in place and record the data
	ForceSum()
	wave ForceSumStats
	if(firsttimeflag == 1)
		Duplicate/O ForceSumStats ForceSumStatsAuto
		firsttimeflag = 0
	else
// COMMENTED OUT DEC. 9 2013 to get a successful compile - ASK CRAIG
//		Concatenate/NP=0 {ForceSumStats},ForceSumStatsAuto 
 	endif	
 endfor  
 Display/N=X ForceSumStatsAuto[][1] vs CutOffValues  //x%
 Display/N=Y ForceSumStatsAuto[][3] vs CutOffValues  //y%
 Display/N=R ForceSumStatsAuto[][5] vs CutOffValues  //r%

//Next, pick the best cutoff value and display to user

//Remove this when the code is correctly selecting the best value
linmaskval = linmaskvalbackup // make sure to restore to original mask
END
//************************************************************************************************************************************

function Calcbgshift(ctrlName): ButtonControl 
string ctrlName
SVAR FitRes_String, AttributeString,MaskString, bgshiftstring

variable iframe1, iframe2
string mask_wave_to_use

LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

//set start and ending frame
iframe1=dimoffset($FitRes_String,2)
iframe2 = dimoffset($FitRes_String,2)+dimsize($FitRes_String, 2)

Make/O/N=(dimsize($FitRes_String,2),2) $bgshiftstring

prompt mask_wave_to_use ,"maskwave" ,popup "manual;automatic"
doprompt "choose mask" mask_wave_to_use


if (stringmatch(mask_wave_to_use,"manual")==1)
	calcbackgroundshift(astr,linmaskval, $FitRes_string, $bgshiftstring,iframe1, iframe2)
else
	calcbackgroundshift(astr,linmaskval_intercept_int, $FitRes_string, $bgshiftstring,iframe1, iframe2)
endif


Save/O/P=ToPractice $bgshiftString  as bgshiftString+".ibw"	

end


//--------button  functionto calculate MSDs
//12/5/21  YS changed to remove plotting.  Plotting moved to PlotALot.
//added to do frequency sweep over different videos     3/1/2016
function computeMSD(ctrlName): ButtonControl 
string ctrlName
SVAR FitRes_String, AttributeString,MaskString, bgshiftstring
variable cellpostflag=0
variable framerate=10
variable spikeflag=0


wave zeropos


LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr


Prompt framerate, "type framerate"
Prompt spikeflag, " type 1 to remove spikes"
DoPrompt "Input Parameters for MSD calculation",framerate, spikeflag

if (V_flag==1)	
		return 0
endif




duplicate/O $fitres_string fitres_plotalot    

calculateMSD(fitres_plotalot, $bgshiftstring, Astr, $Maskstring, framerate, spikeflag)


end

//added 7/14/14
function calcheatmap(ctrlName): ButtonControl 
String ctrlname
SVAR FitRes_String, AttributeString,MaskString, bgshiftstring

LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

variable framerate=100
variable freq=1
variable corf, valf
variable heatflag
string corstring, valstring


Prompt framerate, "framerate of the video"
Prompt freq, "frequency of the signal"
Prompt corstring,"direction",  popup "X;Y"
Prompt valstring, "value", popup "Magnitude; Phase"

doprompt "input parameters",  framerate, freq, corstring, valstring

if ( stringmatch(corstring, "x")==1)
   corf=0
else 
    corf=2
endif

if ( stringmatch(valstring, "Magnitude") ==1 )
    valf=0
else 
   valf=1
endif

heatflag = corf + valf






calcallDLI(ASTR, $fitres_string, $bgshiftstring, $maskstring, framerate, freq, heatflag)



end

function calcMSDheatmap(ctrlname): ButtonControl
String ctrlname
SVAR AttributeString,MaskString

variable mapflag=0

prompt mapflag , "0 for max, 1 for slope"
doprompt "input parameters" , mapflag


LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

MSDstatsmap(MSD, $maskstring, Astr, mapflag)

end

function clacMSDstats(ctrlname): ButtonControl
String ctrlname
SVAR AttributeString, Maskstring

LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

MSDsumary(MSDstats, astr, $maskstring)
end


function calcDFDLIheatmap(ctrlname): ButtonControl
String ctrlname
SVAR FitRes_String, AttributeString,MaskString, bgshiftstring

LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr


variable freq_video = 1.0
variable framerate=100
variable freq_input=1, reffreq = 3
variable corf, valf
variable heatflag
string corstring, valstring
variable initialshift = 0  //for single video set as 0 temporarily.

string fitname,bgname,fitaffix


Prompt freq_video "Driving frequency of the movie to analysis"
Prompt framerate, "framerate of the video"
Prompt freq_input, "frequency of the signal"
Prompt reffreq, "frequency of the reference signal"
Prompt corstring,"direction",  popup "X;Y"
Prompt valstring, "value", popup "Magnitude; Phase"

doprompt "input parameters",  freq_video,framerate, freq_input,  reffreq, corstring, valstring

if ( stringmatch(corstring, "x")==1)
   corf=0
else 
    corf=2
endif

if ( stringmatch(valstring, "Magnitude") ==1 )
    valf=0
else 
   valf=1
endif

heatflag = corf + valf

print fitres_string

wave freq

findValue/V = (freq_video) freq


if (V_value < 0)
	print "Not a valid frequency"

else
	
	if (freq[V_value] <1)
	 	sprintf fitaffix "_0_%dhz"  (freq[V_value]*10)
	else 
	 	sprintf fitaffix "_%dhz"  freq[V_value]
	endif
	fitname = FitRes_String+fitaffix
	bgname = bgshiftstring+fitaffix
	
	
	DFcalcDLI(ASTR, $fitname, $bgname, $maskstring, framerate, freq_input, reffreq, heatflag, initialshift)
	
	
	dowindow/K HeatMapdisp
	newimage/N=HeatMapdisp M_rotatedimage
	appendtograph/T/W=Heatmapdisp linmasky vs linmaskx
	modifygraph mode=3, marker=0
	modifygraph zcolor(linmasky)={heatmap[][heatflag],DLIHeatmap_low,DLIHeatmap_high, Rainbow, 0}	
	DoUpdate
	ColorScale/C/N=text0/F=0/M/A=MC trace=LinMaskY;DelayUpdate
	ColorScale/C/N=text0/A=RC/X=33.49/Y=5.03 width=15,heightPct=50,lblMargin=20;DelayUpdate
	ColorScale/C/N=text0 "\\Z18Deflection Magnitude (\\F'GreekS'm\\F'Arial'm)"	
	ColorScale/C/N=text0/X=5.00/Y=5.03
	
endif


end



function calccor(ctrlname):buttoncontrol
String ctrlname
variable recalccorflag = 1
variable postnum=0
variable direction=0
variable absflag=0
variable lagtime = 1
variable bgflag = 1

wave fitres_ed


SVAR   AttributeString,MaskString

LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

Prompt recalccorflag, "0 if want to recalculate correaltion matrix"
Prompt postnum, "post number for calculating correlation"
Prompt direction, "0 for x , 1 for y"
Prompt absflag, "0 for signed, 1 for absolute"
Prompt  lagtime, "lagtime of map"
Prompt  bgflag, "1 for plot background posts"

doprompt "input parameters" , recalccorflag, postnum, direction, absflag, lagtime, bgflag

if (recalccorflag==0)
//   calculatecorrelation(fitres_ed, astr)
     calccorrelation(fitres_ed, astr, $maskstring)

endif


//corrheatmap_new(corrwave, $maskstring, astr, postnum,lagtime,idcell, bgflag)

end

//----------- COUNTGFitAll()-------------
//fitting posts center for all frequencies for magnetic microrheology
// 12/12/21 modified to ask for movie type:  .avi .MOV etc
//12/8/21 modified for optional showing each image file and printing frequencies to monitor progress
function COUNTGFitAll(ctrlName): ButtonControl
    string CtrlName
    
    NVAR NmovieFrames_G  // total number of movie frames
    string s1,s2
    variable flast = NmovieFrames_G-1
    variable showimageflag = 0   // show each image in fit?
    variable printfrequencyflag = 1  // flag for printing frequencies to monitor progress
    variable printeveryNframes

    variable f1=0,f2=1 // will get set inside COUNTcentroidgfit
    string movietypestring
    movietypestring = "avi"
    Prompt showimageflag, "Set = 1 to show image of each frame during fit (will be slower!)" 
    Prompt printfrequencyflag, " Set = 1 to print frequencies to monitor progress (0 to suppress): " 
    printeveryNframes = max(floor(NmovieFrames_G/10),1)
    Prompt printeveryNframes, " Print frame number every N frames to track progress (= -1 to suppress). N = "   
    prompt movietypestring, " Enter suffix of movie files (avi, MOV, etc): "
    DoPrompt "Fit frequency sweep data", showimageflag, printfrequencyflag, printeveryNframes, movietypestring
    
    SVAR MaskString, FitRes_String, FitErr_String, ShiftString      
    SVAR  AttributeString
    LoadWave/O/P=ToPractice AttributeString+".ibw"
    LoadWave/O/P=ToPractice ShiftString+".ibw"
    Duplicate/O $AttributeString ASTr
 	

   COUNTcentroidgfit2(f1,f2,ASTr[0],ASTr[1],ASTr[3],ASTr[4],ASTr[6],ASTr[7],$MaskString,$ShiftString,showimageflag, printfrequencyflag,printeveryNframes,movietypestring) 


	
end 


//calculate background shift under different frequency
function COUNTCalcbgshift(ctrlName): ButtonControl 
string ctrlName
SVAR FitRes_String, AttributeString,MaskString, bgshiftstring

variable iframe1, iframe2
variable i
wave freq
string fitaffix, fitname, bgname,fit_edname
variable nofreq = dimsize(freq,0)
wave fitres_ed

LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

//set start and ending frame




for (i = 0; i<nofreq; i+=1)
          if (freq[i] <1)
             sprintf fitaffix "_0_%dhz"  (freq[i]*10)
          else 
             sprintf fitaffix "_%dhz"  freq[i]
          endif
          fitname = FitRes_String+fitaffix
          bgname = bgshiftstring+fitaffix
          fit_edname = FitRes_String+"_ed"+fitaffix
         iframe1=dimoffset($fitname,2)
         iframe2 = dimoffset($fitname,2)+dimsize($fitname, 2)
         Make/O/N=(dimsize($fitname,2),2) $bgname
         calcbackgroundshift(astr, $MaskString, $Fitname, $bgname,iframe1, iframe2)
         WAVE fitres_ed  //  output of calcbackgroundshift
          duplicate/O fitres_ed $fit_edname
         Save/O/P=ToPractice $bgname  as bgname+".ibw"	

endfor



end


function COUNTcalcDFDLIheatmap(ctrlname): ButtonControl
String ctrlname
SVAR FitRes_String, AttributeString,MaskString, bgshiftstring

LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

variable framerate=100
variable reffreq = 3
variable corf, valf
variable heatflag
string corstring, valstring
wave freq

variable i 
variable freqthresh
variable nofreq = dimsize(freq,0)
variable initialshift
string fitaffix, fitname, bgname
String heatmapS , heatmapSE, heatmapREFSE, heatmapREFS
string hallS

wave wave2

variable npost = astr[3]*astr[4]
make/o/n=(npost, 2,nofreq)  relativeerr_n
make/o/n=(npost, 18, 17) heatmaprelative1=0

Prompt framerate, "framerate of the video"
//Prompt freq, "frequency of the signal"
Prompt reffreq, "frequency of the reference signal"
Prompt corstring,"direction",  popup "X;Y"
Prompt valstring, "value", popup "Magnitude; Phase"

doprompt "input parameters",  framerate,  reffreq, corstring, valstring

if ( stringmatch(corstring, "x")==1)
   corf=0
else 
    corf=2
endif

if ( stringmatch(valstring, "Magnitude") ==1 )
    valf=0
else 
   valf=1
endif

heatflag = corf + valf
freqthresh = framerate/2

wave relativeerrsingle

for (i = 0; i<nofreq; i+=1)



          if (freq[i] <1)
             sprintf fitaffix "_0_%dhz"  (freq[i]*10)
             sprintf heatmapS "DFDLIheatmap0_%d",(freq[i]*10)
             sprintf heatmapSE "DFDLIheatmaper0_%d", (freq[i]*10)
             sprintf heatmapREFSE "RefDFDLIheatmaper0_%d", (freq[i]*10)
             sprintf heatmapREFS "RefDFDLIheatmap0_%d", (freq[i]*10)
          else 
             sprintf fitaffix "_%dhz"  freq[i]
             sprintf heatmapS "DFDLIheatmap_%d", freq[i]
             sprintf heatmapSE "DFDLIheatmaper_%d", freq[i]
             sprintf heatmapREFSE "RefDFDLIheatmaper_%d", freq[i]
             sprintf heatmapREFS "RefDFDLIheatmap_%d", freq[i] 
          endif
          fitname = FitRes_String+fitaffix
          bgname = bgshiftstring+fitaffix
          
          
          sprintf HallS "magread_%1f.txt", freq[i]
          loadwave/N/O/P=ToPractice/G HallS
    
           initialshift  =  dimsize($fitname,2) - wave2[0]

          if (initialshift<0)
              initialshift = 3
           endif
    
          
          if (freq[i] <freqthresh)
          
          DFcalcDLI(ASTR, $fitname, $bgname, $maskstring, framerate, freq[i], reffreq, heatflag, initialshift)
                    
          else 
          DFcalcDLI(ASTR, $fitname, $bgname, $maskstring, framerate, abs(framerate - freq[i]), reffreq, heatflag, initialshift) 
          
          endif      
          
            
          
          Save/O/P=ToPractice heatmap as heatmapS+".ibw"
          Save/O/P=ToPractice heatmaperror as heatmapSE+".ibw"
          Save/O/P=ToPractice heatmapref as heatmapREFS+".ibw"
          Save/O/P=ToPractice heatmaperrorref as heatmapREFSE+".ibw"
          
          
          if (freq[i] <freqthresh)
          
          DFcalcDLI_break(ASTR, $fitname, $bgname, $maskstring, framerate, freq[i], reffreq, heatflag,i)
          
          else 
          DFcalcDLI_break(ASTR, $fitname, $bgname, $maskstring, framerate, abs(framerate - freq[i]), reffreq, heatflag,i) 
          
          endif      
          
          
          Save/O/P=ToPractice heatmap_bk as heatmapS+"_bk.ibw"
          Save/O/P=ToPractice heatmaperror_bk as heatmapSE+"_bk.ibw"
          Save/O/P=ToPractice heatmapref_bk as heatmapREFS+"_bk.ibw"
          Save/O/P=ToPractice heatmaperrorref_bk as heatmapREFSE+"_bk.ibw"
          
			if (i==0)   //plot DLI heatmap for f = 0.1 Hz only
				dowindow/K HeatMapdisp
				newimage/N=HeatMapdisp M_rotatedimage
				appendtograph/T/W=Heatmapdisp linmasky vs linmaskx
				modifygraph mode=3, marker=0
				modifygraph zcolor(linmasky)={heatmap[][heatflag],DLIHeatmap_low,DLIHeatmap_high, Rainbow, 0}
				DoUpdate 
				ColorScale/C/N=text0/F=0/M/A=MC trace=LinMaskY;DelayUpdate
				ColorScale/C/N=text0/A=RC/X=33.49/Y=5.03 width=15,heightPct=50,lblMargin=20;DelayUpdate
				ColorScale/C/N=text0 "\\Z18Deflection Magnitude (\\F'GreekS'm\\F'Arial'm)"	
				ColorScale/C/N=text0/X=5.00/Y=5.03         
          endif

endfor


          

end




function findanisotropy(datawave,astr)
wave  datawave, astr

wave linmaskval_intercept_int

variable npost = dimsize(datawave,0)
variable movielength = dimsize(datawave, 2)
variable selectratio = 3
variable cloudlength = movielength/selectratio
variable nxpost = astr[3]
variable nypost = astr[4]

variable ipost, i, j=0



make/o/n=(npost) cposanisotropy, cposanisotropy_raw
make/o/n=(cloudlength,2) poscloud

//alternative method of calculating anisotropy index
make/o/n=(movielength) poscloudx, poscloudy
make/o/n=(movielength) poscloudxx, poscloudyy, poscloudxy
variable covxx, covxy, covyy
make/o/n=(npost) cposanisotropy_new

wave W_eigen

variable icell=0


for ( ipost = 0; ipost < npost; ipost+=1)

      if(linmaskval_intercept_int[ipost] == 1)
       j=0
       for(i = 0; i < movielength; i+=selectratio)
       poscloud[j][0] = datawave[ipost][0][i]
       poscloud[j][1] = datawave[ipost][1][i]
       j+=1
       endfor
       matrixop/o poscloud = subtractmean(poscloud,1)      
       
       PCA/SEVC poscloud
       cposanisotropy[icell] = (W_eigen[0] - W_eigen[1])/(W_eigen[0] + W_eigen[1]) 
       cposanisotropy_raw[icell] = abs(W_eigen[0]/W_eigen[1])
      poscloudx = datawave[ipost][0][p]
      poscloudy = datawave[ipost][1][p]
      matrixop/o poscloudx = subtractmean(poscloudx,0)
      matrixop/o poscloudy = subtractmean(poscloudy,0)
      poscloudxx = poscloudx*poscloudx
      poscloudyy = poscloudy*poscloudy
      poscloudxy = poscloudx*poscloudy
      covxx = mean(poscloudxx)
      covyy = mean(poscloudyy)
      covxy = mean(poscloudxy)
      cposanisotropy_new[icell] = sqrt((covxx-covyy)^2 + 4*covxy^2)/(covxx+covyy) 
      icell+=1
endif

       
endfor

//duplicate/o M_rotatedimage M_rotatedimage_anisotropy

//dowindow/K anisotropy_heatmap
//newimage/N=anisotropy_heatmap M_rotatedimage_anisotropy
//appendtograph/T/W=anisotropy_heatmap linmasky vs linmaskx
//modifygraph mode =3, marker =0
//modifygraph zcolor(linmaskY)={cposanisotropy,0,1,Rainbow, 0}

redimension/n=(icell) cposanisotropy, cposanisotropy_raw, cposanisotropy_new
       
end


function findanisotropy_seg(datawave,astr, seglength, posttype)
wave  datawave, astr
variable seglength , posttype

wave linmaskval_intercept_int
wave fitres_ed


variable lamdathresh = 5

variable npost = dimsize(datawave,0)
variable movielength = dimsize(datawave, 2)
variable cloudlength = seglength
variable nseg = floor(movielength/seglength)
variable nxpost = astr[3]
variable nypost = astr[4]

variable ipost, i, j=0

make/o/n=(npost) cpostid

variable xdisp, ydisp
make/O/N=(npost, 2)   directw
variable CX, CY

variable ncell = 0
for(ipost=0; ipost<npost; ipost+=1)
         if( (linmaskval_intercept_int[ipost]>0) && (linmaskval_intercept_int[ipost]!=3) ) 	
            CX+=fitres_ed[ipost][0][0]
            CY+=fitres_ed[ipost][1][0]
            ncell+=1   
         endif
endfor         

CX=CX/ncell
CY=CY/ncell     





make/o/n=(npost*nseg) cposanisotropy, cposanisotropy_raw, angledif
make/o/n=(npost, nseg) cposanisotropy_raw_mat=0, angledif_mat=0
make/o/n=(cloudlength,2) poscloud

//alternative method of calculating anisotropy index
make/o/n=(cloudlength) poscloudx, poscloudy
make/o/n=(cloudlength) poscloudxx, poscloudyy, poscloudxy
variable covxx, covxy, covyy
make/o/n=(npost) cposanisotropy_num=0

wave W_eigen, M_C

variable icell=0

variable alpha1, beta1


variable iangle =  0



for ( ipost = 0; ipost < npost; ipost+=1)

   if(linmaskval_intercept_int[ipost] == posttype)
//     if( (linmaskval_intercept_int[ipost] >0) && (linmaskval_intercept_int[ipost]!=3) )
      cpostid[icell] = ipost
     for(j=0; j<nseg; j+=1)
       for(i = 0; i <seglength ; i+=1)
       poscloud[i][0] = datawave[ipost][0][i+j*seglength]
       poscloud[i][1] = datawave[ipost][1][i+j*seglength]
       endfor
       matrixop/o poscloud = subtractmean(poscloud,1)      
       
       PCA/SEVC/SCMT poscloud
       cposanisotropy[icell*nseg+j] = (W_eigen[0] - W_eigen[1])/(W_eigen[0] + W_eigen[1]) 
       cposanisotropy_raw[icell*nseg+j] = abs(W_eigen[0]/W_eigen[1])
       
       cposanisotropy_raw_mat[icell][j] =  abs(W_eigen[0]/W_eigen[1])     
         
      poscloudx = datawave[ipost][0][p]
      poscloudy = datawave[ipost][1][p]
   //   matrixop/o poscloudx = subtractmean(poscloudx,0)
  //    matrixop/o poscloudy = subtractmean(poscloudy,0)
  //    poscloudxx = poscloudx*poscloudx
  //    poscloudyy = poscloudy*poscloudy
  //    poscloudxy = poscloudx*poscloudy
  //    covxx = mean(poscloudxx)
  //    covyy = mean(poscloudyy)
  //    covxy = mean(poscloudxy)
  //    cposanisotropy_new[icell*nseg+j] = sqrt((covxx-covyy)^2 + 4*covxy^2)/(covxx+covyy) 
      
      

      //         xdisp = fitres_ed[ipost][0][0] - CX
//         ydisp = fitres_ed[ipost][1][0] - CY

           xdisp = datawave[ipost][0][j*seglength]
           ydisp = datawave[ipost][1][j*seglength]
         
       directw[ipost][0] = xdisp/sqrt(xdisp^2+ydisp^2)
       directw[ipost][1] = ydisp/sqrt(xdisp^2+ydisp^2)        
       
       alpha1 = atan(ydisp/xdisp)
       beta1 = atan(M_C[0][1]/M_C[0][0])
       
     
    
//    angledif[iangle] = asin(directw[ipost][1]*M_c[0][0] - directw[ipost][0]*M_c[0][1])
      if ( (alpha1 - beta1)>-pi/2 && (alpha1 - beta1)<pi/2 )
             angledif[iangle] = alpha1 - beta1
      elseif ( (alpha1 - beta1)<-pi/2)
               angledif[iangle] = pi + (alpha1 - beta1)
       elseif (  (alpha1 - beta1)>pi/2)      
               angledif[iangle] =-pi + ( alpha1 - beta1)
        endif       
        angledif_mat[icell][j] = angledif[iangle]
    iangle+=1

    if(abs(W_eigen[0]/W_eigen[1])>lamdathresh)
         cposanisotropy_num[icell] +=1
    
    endif      
      
      
      
   endfor
   
   icell+=1

endif

       
endfor

//duplicate/o M_rotatedimage M_rotatedimage_anisotropy

//dowindow/K anisotropy_heatmap
//newimage/N=anisotropy_heatmap M_rotatedimage_anisotropy
//appendtograph/T/W=anisotropy_heatmap linmasky vs linmaskx
//modifygraph mode =3, marker =0
//modifygraph zcolor(linmaskY)={cposanisotropy,0,1,Rainbow, 0}

redimension/n=(icell*nseg) cposanisotropy, cposanisotropy_raw
redimension/n=(iangle) angledif
redimension/n=(icell,nseg) cposanisotropy_raw_mat, angledif_mat
redimension/n=(icell) cposanisotropy_num, cpostid


       
end





function  calchistogram(ctrlName): ButtonControl 
string ctrlName
SVAR FitRes_String, AttributeString,MaskString, bgshiftstring
variable lagtime =1 , freq =1, startbin = 0.2 , posttype = 1
wave zeropos


LoadWave/O/P=ToPractice AttributeString+".ibw"
Duplicate/O $AttributeString ASTr

string fitname, bgname
string fitaffix


//string    filename="cellong2_5"
//string    f1="fitres_"+filename
prompt  lagtime, "enter lagtime"
prompt  freq, "frequency of trace"
prompt  startbin, "smallest bin"
prompt posttype, " type of posts"
//prompt  filename, "enter file folder name
doprompt "input all the stuff necessary",  freq,  lagtime, startbin,posttype 

    if (freq ==0)
           fitname = FitRes_String
           bgname = bgshiftstring
     elseif (freq !=0)      
           if (freq<1)
             sprintf fitaffix "_0_%dhz"  (freq*10)
          else 
             sprintf fitaffix "_%dhz"  freq
          endif

          fitname = FitRes_String+fitaffix
          
          bgname = bgshiftstring+fitaffix
     
     endif
drawhisto(ASTR, $fitname, $bgname, $maskstring, zeropos, lagtime, -startbin,posttype)

//drawhisto_cellref(ASTR, $fitname, $bgname, $maskstring, zeropos, lagtime, -startbin)

end



function calcvecfsum(datawave, maskwave)
wave datawave
wave maskwave

variable npost = dimsize(maskwave,0)
variable movielength = dimsize(datawave, 2)

variable ipost, iframe

variable kspring = 22.3
variable ncp1=0, ncp2=0


make/o/n=(movielength) vfsum_c1x=0, vfsum_c1y=0, vfsum_c2x=0, vfsum_c2y=0
make/o/n=(movielength) ftempx, ftempy
make/o/n=(movielength) vfsum_c1m, vfsum_c2m, vfsum_tm, vfsum_tx, vfsum_ty


for(ipost=0; ipost<npost; ipost+=1) 
        if(maskwave[ipost] == 1)
          ftempx = datawave[ipost][0][p]*kspring
          ftempy = datawave[ipost][1][p]*kspring
          vfsum_c1x = vfsum_c1x + ftempx
          vfsum_c1y = vfsum_c1y+ ftempy
          ncp1+=1
        elseif(maskwave[ipost] == 4)
          ftempx = datawave[ipost][0][p]*kspring
          ftempy = datawave[ipost][1][p]*kspring        
          vfsum_c2x = vfsum_c2x + ftempx          
          vfsum_c2y = vfsum_c2y + ftempy
         ncp2+=1
         endif
endfor

         
vfsum_c1m = sqrt(vfsum_c1x^2+vfsum_c1y^2)
vfsum_c2m = sqrt(vfsum_c2x^2+vfsum_c2y^2)
vfsum_tx = vfsum_c1x + vfsum_c2x
vfsum_ty = vfsum_c1y + vfsum_c2y
vfsum_tm = sqrt(vfsum_tx^2+vfsum_ty^2)



end




function calcvecfsum_2(datawave, datawave2,maskwave)
wave datawave
wave datawave2
wave maskwave


variable npost = dimsize(maskwave,0)

variable ipost, iframe

variable kspring = 22.3
variable ncp1=0, ncp2=0


variable vfsum_c1x=0, vfsum_c1y=0, vfsum_c2x=0, vfsum_c2y=0
variable ftempx, ftempy
variable vfsum_c1m, vfsum_c2m, vfsum_tm, vfsum_tx, vfsum_ty


for(ipost=0; ipost<npost; ipost+=1) 
        if(maskwave[ipost] == 1)
          ftempx = datawave[ipost][0]*cos(datawave[ipost][1])
          ftempy = datawave[ipost][1]*sin(datawave[ipost][1])
          vfsum_c1x = vfsum_c1x + ftempx
          vfsum_c1y = vfsum_c1y+ ftempy
          ncp1+=1
        elseif(maskwave[ipost] == 4)
          ftempx = datawave2[ipost][0]*cos(datawave2[ipost][1])
          ftempy = datawave2[ipost][1]*sin(datawave2[ipost][1])    
          vfsum_c2x = vfsum_c2x + ftempx          
          vfsum_c2y = vfsum_c2y + ftempy
         ncp2+=1
         endif
endfor

         
vfsum_c1m = sqrt(vfsum_c1x^2+vfsum_c1y^2)
vfsum_c2m = sqrt(vfsum_c2x^2+vfsum_c2y^2)
vfsum_tx = vfsum_c1x + vfsum_c2x
vfsum_ty = vfsum_c1y + vfsum_c2y
vfsum_tm = sqrt(vfsum_tx^2+vfsum_ty^2)


string windowname
for(iframe=2;iframe<347; iframe+=1)
     sprintf windowname, "bonddist%d", iframe
     killwindow $windowname
endfor
     



end





//***** Have some misfit blue posts you can't seem to find and ignore? Run this*****

Function fixblue()

wave linmaskval, bluelengthandangle
variable i

for(i=0;i<(dimsize(bluelengthandangle,0));i+=1)
	if(bluelengthandangle[i] > 5)
		linmaskval[i]=3
	endif
endfor

END
//************************************************************

