#pragma rtGlobals=1		// Use modern global access method.



// *************MODIFICATION HISTORY  - EFFECTIVE MAY 19, 2006
//11/3/14    added allowing analysis videos continuously

//6/16/14 use ROI for image processing

//6/5/14  change Gaussian fitting into centroid mask method

//12/18/13  use bandpass to do preprocessing

// 12/10/13 modified to read from movie

// 4/13/07  permits look at slopes of interp fits

// 4/5/07  COntains error checking for Version 6
   //    x0,y0 in ROI,    SigMin < abs(sigX), abs(sigY) < SigMax, abs(corr) < CorrMax, V_FitError !=0
   
//   3/24/07 NEW VERSION for Interp all frames

// 6/28/06:  Error bar propagation and CCD error option added
//       doGfit():   CCD error option
//    gfitforshifts()

//5/31/06  doPSFfit
//        Fit range could go out of bounds for posts on edge of image.  It now stays 5 pixels away from any boundary of image  (hopefully enough!
//5/31/06 PSFfitBaseinterpG()
//     fitresultsbase, fiterrorsbase have  fits for guideposts, interp X0,Y0 for all others

//static strconstant  ksRhodFileRoot = "RHOD"  // for 2-digit names
static strconstant  ksRhodFileRoot = "RHOD"  // for 3-digit names




// maskwave constants

static constant kEmpty = 0
static constant kCell = 1
static constant kGuide = 2
static constant kIgnore = 3
static constant kWirepost = 4
static constant kFitAll = -1
static constant kInterp = 5     // 5/24/06  new constant interp base fit

// indices of Shiftwave
// note that primary index of Shiftwave must correspond to file number
static constant kiXc = 0   // 4/15/07  Xc
static constant kiYc = 1  // Yc
static constant kisigXc = 2  // sigXc
static constant kisigYc = 3  // sigXy
static constant kiTheta = 4  // shift angle
static constant kisigXsf = 5   // 6/28/06 for error
static constant kisigYsf = 6   
static constant kiAngle = 7
static constant kiA = 8
static constant kiB = 9
static constant kiXsi  =10
static constant kiYsi = 11
static constant kiXsf = 12
static constant kiYsf =13

// indices for fitparams
static constant kNfitParamsG  = 9  // includes columns for chisq and FitErrorReporting

static constant kFbkgnd = 0
static constant kFAmpl = 1
static constant kFX0 = 0     //temporary change    7/11/14
static constant kFSigX= 3
static constant kFY0 = 1    //temporary change     7/11/14
static constant kFSigY = 5
static constant kFCorr = 6
static constant kFChisqG = 7
static constant kFerrorFlag = 8

// these used for G + quadratic for basefit
static constant kNfitParamsBase  = 9  // includes column for chisq
static constant kFB = 7       // ampl. of quadratic bit
static constant kFChisqB = 8

//static constant kfitrange  = 20   //+- range for fitting in pixels

// constants needed by fitrangewave
static constant kxmin = 0
static constant kxmax = 1
static constant kymin = 2
static constant kymax = 3

 // 5/19/06   reduction factor for reduced fitrange for refit.
 static constant kar = 0.75      

static constant deltaROI  = 1   // number of pixels to reduce ROI in each direction for refit-on-error
static constant SigMin = 1      //  minimum SigmaX and SigmaY for Gaussian fit
static constant SigMax =15
static constant CorrMax = 1.0

//------------------ GfitAllPosts() ---------------------

// New3/24/07

 function GfitAllPosts(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,autorefit,maskwave,shiftwave,InitGuessWave,pathstr)
    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    variable zmin,zmax      // min and max grayscale levels for display during fitting
    variable autorefit        //   = 1 for autorefit - fits each  post twice - second time with kfitrange reduced to 0.3 from 0.4
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
    wave InitGuessWave   // Initial guesses for fit to first post.   These will be copied into w_coef1[]
                                        //  Should have same dimensions an definitions as w_coef1.  As of 2/10/06 this is length = 7
                                        //  Initial guesses for X0 and Y0 (post center coords) and post bottoms need not be supplied as initial guesses
                                        //  Initial guesses needed:
                                        //  InitGuessWave[0] =   Background
                                        //  InitGuessWave[1] =   Ampl
                                         //  InitGuessWave[2] =  Post Radius (in um).  This is a fixed value, presumably determined 
                                         // by a prior call to PSFfitForShifts followed by suitable averaging.
    string pathstr            // path to data files
	nvar hexcp
	

  
  	 Variable/G kfitrange=floor(shiftwave[0][kIA]*0.4)    // added 6/28/06. Otherwise defined in doGr2FitBase 
   
   
   make/O/N=4 fitrangewave
 // fitrangewave = kfitrange //commented out 9/11/08 by CMK //have to think about this...
   if(hexcp==1)
	   fitrangewave[kxmin]=floor(shiftwave[0][kiA])
  	   fitrangewave[kxmax]= fitrangewave[kxmin]
 	   fitrangewave[kymin]=floor( shiftwave[0][kiB])
  	   fitrangewave[kymax] = fitrangewave[kymin]
  	else
  		 fitrangewave[kxmin]=.75*floor(shiftwave[0][kiA])
   		 fitrangewave[kxmax]= fitrangewave[kxmin]
  		 fitrangewave[kymin]=.75*floor( shiftwave[0][kiB])
  		 fitrangewave[kymax] = fitrangewave[kymin]
      endif
   make/O/n=(XPosts*YPosts,kNFitParamsG,f2-f1+1) fitresults // always created. calling routine must copy it as needed
   SetScale/I z f1,f2,"", fitresults                   
   fitresults = 0
   make/O/n=(XPosts*YPosts,kNFitParamsG-2,f2-f1+1) fiterrors // always created. two less columns - no chisq column or errorflag column needed
   SetScale/I z f1,f2,"", fiterrors                   
   fiterrors = 0
  
  make/o/n=(1,2) refitwave //stores post number and file number of major screwups for use in refitting

//fit Cell posts
    doGfitAll(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,autorefit,fitrangewave,maskwave,shiftwave,InitGuessWave,fitresults,fiterrors,pathstr)

end

//CMK added function buildhexmask in order to design a mask to maximize fitting area
function buildhexmask(xrmin, xrmax, yrmin, yrmax,buffer)

variable xrmin, xrmax, yrmin, yrmax,buffer


make/D/O/N=(dimsize(RotatedOrig1,0),dimsize(RotatedOrig1,1)) hexmask=0

variable p, q //indices for filling hex lattice
variable m =0,l= 0// used to increase/decrease number of columns = 1 to build hexagon 
variable xdim, ydim

xdim= xrmax-xrmin
ydim= yrmax-yrmin

//printf "xrmin = %d\t yrmin = %d\r", xrmin, yrmin 
//printf "xdim = %d\t ydim = %d\r", xdim, ydim  

   		   
for(p =(yrmin+buffer); p<(yrmin+(.35*ydim));p+=1)
	
		for(l=0; l<(p-yrmin-buffer); l+=1)
		
			hexmask[xrmin + xdim/2+l][p]=1
	    	       hexmask[xrmin+ xdim/2-l][p]=1
	              
	    	 endfor
		m+=1
endfor

for(p=(yrmin+(.35*ydim));p<(yrmin+(.65*ydim));p+=1)

	       for(l=0; l<(m); l+=1)
		
			hexmask[xrmin + xdim/2+l][p]=1
	    	       hexmask[xrmin + xdim/2-l][p]=1
	       
	       endfor
endfor
	
for(p=(yrmin+ (.65*ydim));p<(yrmax-buffer);p+=1)
		    
	       for(l=0; l<(m); l+=1)
		
			hexmask[xrmin+ xdim/2+l][p]=1
	    	       hexmask[xrmin+ xdim/2-l][p]=1
	       
	       endfor
	       m-=1
endfor	

end	

//*****************************************************
//   doGfitAll

// 12/10/13 modified to read from movie

// 3/24/07 modified from DoGFit2 - FITS ALL POSTS

//6/28/06 ccd error option added
//5/31/06  Image Boundary checking on Fitting ROI implemented
//5/26/06  Structure taken from doPSFfit()

// 5/24/06 - Fitting for interp. "base fit" added
//5/23/06  Bug fix - gets shifts right
// 5/22/06 New initial guess algorithm for cell/empty posts
//    first image:  ampl,bkgnd from previous fit post in first image
//                         post top center:  from base position
//    all other images:   all params from fit to current post in previous image
// 5/21/06  For initial guess in fits to posts in first image, count back to find previous non-Ignored post (was only looking one post back previously)
 // 5/19/06  Fits empty, cell and wire posts in LUT fit to enable Sasha's method of comparing
 //                   empty and cell posts for background.
 // 5/19/06  autorefit option added  (can do it for any post) 
 
//   2/10/06    General fitting routine using PSF fitting method.
//    Structure taken from doGfit()
 
//  OUTPUT
//  always generates fitresults wave which must be copied to avoid overwriting
//  This version of the fitresults wave uses its "z" index for the images, and has
//   its initial index value set to the first file number fit.  (Previous versions were zero-offset always)

//  also records Igor's estimate of errors of the parameters in wave "fiterrors" which has same dimensions
//  as fitresults




function doGfitAll(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,autorefit,fitrangewave,maskwave,shiftwave,InitGuessWave,FRwave,FEwave,pathstr)

    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    variable zmin,zmax      // min and max grayscale levels for display during fitting
    variable autorefit          // 5/19/06 = 1 for autorefit of each post.
    wave fitrangewave        // 3/8/06 has fitrange increments - default = kfitrange
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
    wave InitGuessWave   // Initial guesses for fit to first post.   These will be copied into w_coef1[]
                                        //  Should have same dimensions an definitions as w_coef1. 
                                        //  Initial guesses for X0 and Y0 (post center coords) need not be supplied 
       					//10/15/08 CMK adding in initial guesses for xo, yo	
       					//  Initial guesses needed:
                                        //  InitGuessWave[0] =   Background
                                        //  InitGuessWave[1] =   Ampl
  wave FRwave			// fit results wave.  Will be updated appropriately
  wave FEwave			// fit errors wave.  Will be updated appropriately
  string pathstr            // path to data files
    
   wave hexmask  //file containing 1's in hexagonal shape for maximum fitting area 9/11/08
    wave refitmask //stores badly fit posts and the file number
    
    
   variable iprev, jprev // 5/21/06 for fits in first image
   variable xrmin,xrmax,yrmin,yrmax   // 5/19/06  new for refit range
   string s1
   variable i, j, ifile, fitthispost
   variable xcenter, ycenter, Xsp, Ysp
   variable V_chisq
   variable AngleOfRotation
   variable dX,dY   // local variables for lattice constants of post array
   NVAR kfitrange
   NVAR CCDflag, CCDreadnoise, CCDeperADU   // 6/28/06
   nvar hexcp
   
   variable npost1,varkfx0,varkfy0,varkfsigx,varkfsigy,varkfbkgnd,varkfampl,varkfcorr //for potential refit
   
	//global variables needed for fitting with constraints
	variable/g g_xrmin, g_xrmax
	variable/g g_yrmin, g_yrmax
   

   
   variable hexmaskflag = 1 //setting for use of hexmask 9/11/08 
   
      // make w_coef1 and insert initial guesses  - positions obtained elsewhere.
      //  these are used once, for first fit.  After that, the previous fit is used as initial guess
 	make/D/O/n=(kNfitParamsG -2) w_coef1   // 4/4/07
 	make/D/O/n=(kNfitParamsG -2) W_Sigma
 
	variable V_FitError   //  4/4/07 - for curvefit error reporting
 	variable ErrorReporter, ErrorReporter2
 
 	wave xt2, yt2 //duplicates of the initial 3 locations used in roiselect, here used to determine angles
 
    variable zscale = 1 // 6/28/06 for ccd
    variable zoffset = 0
    variable zmin1 =zmin
    variable  zmax1 = zmax
    variable slopey =  0//for use instead of rotate image
    variable slopex = 0
    variable counter =0
    variable counter2=0
    
    
  if(hexcp==1)   
  
  	slopex=(YT2[1]-YT2[0])/(XT2[1]-XT2[0])
      if(mod(yposts,2)==0)
      slopey =  (XT2[2]-0.5*shiftwave[0][kiA]-XT2[0])/(YT2[2]-YT2[0])
      else
      slopey =  (XT2[2]-XT2[0])/(YT2[2]-YT2[0])
    endif
  else
  	
    slopey =  (XT2[2]-XT2[0])/(YT2[2]-YT2[0])
     slopex = (YT2[1]-YT2[0])/(XT2[1]-XT2[0])
   endif
// new 12/10/13 for movie input
   NVAR InputMovieID_G, NmovieFrames_G   //12/10/13
   WAVE M_MovieFrame,M_rgb2gray
   PlayMovieAction setFrontMovie = inputmovieID_G
   PlayMovieAction stop, frame = f1                // go to first frame to be analyzed
   for (ifile=f1; ifile <=f2; ifile +=1)   // loop over raw data files
       PlayMovieAction extract
       if (ifile < f2)
          PlayMovieAction step = 1   //  move to next frame to be ready for next pass through loop
       endif
// new 12/10/13 for movie input
       ImageTransform rgb2gray M_MovieFrame   // convert to grayscale in M_rgb2gray
       Duplicate/O M_rgb2gray Image1

// 12/10/13        sprintf s1, "%s%03d.tif", ksRhodFileRoot,ifile   // 4/3/07
       AngleOfRotation = shiftwave[ifile][kiAngle]
//  12/10/13   	ImageLoad/P=MovieFolder/T=tiff/O/N=image1 s1    // load, copy, and rotate the image
       Duplicate/O image1 OrigImage1
       if (CCDflag ==1) // use proper statistical errors from CCD  6/28/06
          origImage1 = OrigImage1 - CCDBias(x,y)
       endif
       ImageRotate/A=(AngleOfRotation) OrigImage1
       Duplicate/O M_RotatedImage RotatedOrig1
       if (CCDflag == 1)  //6/28/06
           RotatedOrig1 *= CCDeperADU  // convert to e-
           Duplicate/O RotatedOrig1 RotatedOrig1Err
           RotatedOrig1Err = sqrt(RotatedOrig1Err + CCDreadnoise^2)  
           zscale = CCDeperADU // conversion factor to e-
           zmin1 = CCDeperADU*(zmin - CCDavgBias())
           zmax1 = CCDeperADU*(zmax - CCDavgBias())
           zoffset =  CCDeperADU* CCDavgBias()
       endif
       DoWindow/K RotatedImage
       dowindow/k thresh //for testing
       newimage/n=thresh m_imagethresh
       NewImage/N=RotatedImage RotatedOrig1
       ModifyImage RotatedOrig1 ctab= {zmin1,zmax1,Grays,0}
       
       
       //use bandpass to do preprocess   Yu Shi
       bpass(20,3, rotatedorig1)
       duplicate/o mfinal rotatedorig1
       
       
       //include smoothing commands from Yu Shi
//	   MatrixFilter/N=(3)/P=1 point rotatedorig1
//	   MatrixFilter/N=(9)/P=1 gauss rotatedorig1
          
	//smooth/B/E=3 10,RotatedOrig1
	 //******preprocessing goes here******
//12/10/13		MatrixFilter/N=(3)/P=3 avg rotatedorig1
		//below was sued for fluroescence
		//MatrixFilter/N=(3)/P=5 min rotatedorig1
		//MatrixFilter/N=(3)/P=3 avg rotatedorig1
	 //*****************************          
	           
       // Use Initial Shifts
       xSp = kiXsi
       ySp = kiYsi
       
       // post array lattice constants
       dX = shiftwave[ifile][kiA]
       dY = shiftwave[ifile][kiB] //now defined in define array as sqrt(3)/2 lat b if hexcp ==1 or just kia if not
       
       // Draw a box around ROI
       SetDrawLayer ProgFront
       SetDrawEnv linefgc= (65535,65535,0),fillpat= 0,xcoord= top,ycoord= left, save 
       
    
       Drawrect XTop + shiftwave[ifile][xSp] - dX/2, YTop +  shiftwave[ifile][ySp] - dY/2 ,XTop +  shiftwave[ifile][xSp] + (XPosts- 0.5)*dX, YTop +  shiftwave[ifile][ySp] +(YPosts- 0.5)*dY
       SetAxis left  YTop + shiftwave[ifile][ySp] + YPosts*dY, YTop + shiftwave[ifile][ySp]- dY
       SetAxis top XTop + shiftwave[ifile][xSp] - dX, XTop + shiftwave[ifile][xSp] + XPosts*dX
           
       // Fit the posts
      
      make/D/O/n=(Xposts, Yposts) G_x //global variable for xcenter for constraints
       make/D/O/n=(Xposts, Yposts) G_y //global variable for xcenter for constraints
    
       for (i=0;i<=(Xposts-1); i+=1)
	    for (j=0;j<=(YPosts-1); j+=1)
	       // determine whether or not to fit this post
	       if (maskwave[i][j] == KIgnore)
	           fitthispost = 0
	       else
	           fitthispost = 1
	      endif
	      if (fitthispost == 1)
// Initial guess algorithm:
//        First frame:  use initial guesses
//        all other frames, use fit to current post in previous image
	           if (ifile == f1)  // first frame
                     
                     //if(mod(j,2)==0) //even rows initial guesses-- odd rows guesses are just shifted to the right by 0.5 lattice
                	
                	   	variable xoffset
                      w_coef1[kFbkgnd]  = InitGuessWave[kFbkgnd]  * zscale // 6/28/06
                      w_coef1[kFAmpl]  = InitGuessWave[kFAmpl] * zscale    // 6/28/06
                     
                      w_coef1[kFSigX]  = InitGuessWave[kFSigX]      //  
              
                      w_coef1[kFSigY]  = InitGuessWave[kFSigY]      //  
                      w_coef1[kFCorr]  = InitGuessWave[kFCorr]      //  
                  	
                       if(hexcp==1)
                       	  if(mod(j,2)!=0)
                       		xoffset = 0.5 *dx
                       	  else 
                       	  	xoffset = 0
                       	 
                       	  endif
                       	else
                       		xoffset = 0
                       	endif	
                  	  
                     xcenter = XTop + shiftwave[ifile][xSp] + i*dX + xoffset  + j * dx*slopey  // center of ROI, shifted by guess -- for odd rows, shifted over by 0.5 lattice spacing 
	                ycenter = YTop + shiftwave[ifile][ySp] + j*dY + i  * dy* slopex
	                w_coef1[kFX0] = xcenter    // x center of top of post - initial guess is zero displacement
	                w_coef1[kFY0] = ycenter    // ycenter of top of post - initial guess is zero displacement
                       
                       //need global variable to put into constraints wave -- these should remain constant, as opposed to xrmin, xrmax, which are adjustable 
   		           //adjusted by code. These are designed to be slightly larger (1.2*) in range than the "fitrangewave " just for safety to not miss the desired 
                        //post but should avoid any fitting of the wrong post
   		       
   		       	g_xrmin =xcenter-  fitrangewave[kxmin]
   				g_xrmax = xcenter + fitrangewave[kxmax]
   				g_yrmin =  ycenter-  fitrangewave[kymin]
   				g_yrmax = ycenter+  fitrangewave[kymax] 
   				
   				g_x[i][j] = xcenter
   				g_y[i][j] = ycenter
                       
                      else//ifile !=1
                      w_coef1[kFbkgnd]  = InitGuessWave[kFbkgnd]  * zscale  // 6/28/06
                      w_coef1[kFAmpl]  = InitGuessWave[kFAmpl] * zscale    // 6/28/06
                      w_coef1[kFSigX]  = InitGuessWave[kFSigX]      //  
                      w_coef1[kFSigY]  = InitGuessWave[kFSigY]      //  
                      w_coef1[kFCorr]  = InitGuessWave[kFCorr]      //  
                       	
                     
                       	
                       	if(hexcp==1)
                       	  if(mod(j,2)!=0)
                       		xoffset = 0.5 *dx
                       	   else 
                       	   	xoffset =0
                       	  endif
                       	else
                       		xoffset = 0
                       	endif
                       	
                      xcenter = XTop + shiftwave[ifile][xSp] + i*dX + xoffset  + j * dx*slopey  // center of ROI, shifted by guess -- for odd rows, shifted over by 0.5 lattice spacing 
	                ycenter = YTop + shiftwave[ifile][ySp] + j*dY + i  * dy* slopex
	                w_coef1[kFX0] = xcenter    // x center of top of post - initial guess is zero displacement
	                w_coef1[kFY0] = ycenter    // ycenter of top of post - initial guess is zero displacement 	
                       
                       //need global variable to put into constraints wave -- these should remain constant, as opposed to xrmin, xrmax, which are 
                       //adjusted by code. These are designed to be slightly larger (2*) in range than the "fitrangewave " just for safety to not miss the desired 
                        //post but should avoid any fitting of the wrong post
                       	g_xrmin =xcenter-fitrangewave[kxmin]
   				g_xrmax =xcenter+ fitrangewave[kxmax]
   				g_yrmin =  ycenter-fitrangewave[kymin]
   				g_yrmax = ycenter + fitrangewave[kymax] 
   				
   				g_x[i][j] = xcenter
   				g_y[i][j] = ycenter
   				
                       endif
                       
//                  else    // all other frames  -- 5/22/06 initial guess from fit to this post in previous image
//                  
//                  
//                        g_xrmin = g_x[i][j]+shiftwave[ifile][xsp]-  fitrangewave[kxmin]
//   				g_xrmax = g_x[i][j] +shiftwave[ifile][xsp]+ fitrangewave[kxmax]
//   				g_yrmin = g_y[i][j] + shiftwave[ifile][ySp] - fitrangewave[kymin]
//   				g_yrmax = g_y[i][j]+ shiftwave[ifile][ySp]+  fitrangewave[kymax] 
//                  
//                       w_coef1[kFbkgnd]  =  FRwave[i* YPosts + j][kFbkgnd][ifile-1]*zscale - zoffset    // bggnd  from previous fit to this post  
//                       w_coef1[kFAmpl]  = FRwave[i* YPosts + j][kFAmpl][ifile-1] *zscale  // Amplitude from previous fit to this post
//                       w_coef1[kFSigX]  = FRwave[i* YPosts + j][kFSigX][ifile-1]   // Amplitude from previous fit to this post
//                       w_coef1[kFSigY]  = FRwave[i* YPosts + j][kFSigY][ifile-1]   // Amplitude from previous fit to this post
//                       w_coef1[kFCorr]  = FRwave[i* YPosts + j][kFCorr][ifile-1]   // Amplitude from previous fit to this post
// 5/23/6   Bug fix: in next two lines, have to subtract off previous shift, and add the current shift
//	                 w_coef1[kFX0] = FRwave[i* YPosts + j][kFX0][ifile-1]  - shiftwave[ifile-1][xSp] + shiftwave[ifile][xSp]    // x center of top of post 
//	                 w_coef1[kFY0] = FRwave[i* YPosts + j][kFY0][ifile-1]  - shiftwave[ifile-1][ySp] + shiftwave[ifile][ySp]   // ycenter of top of post - 
//                  endif
                  duplicate/O w_coef1 w_coefInit
                  
                  xrmin = max(dimoffset(RotatedOrig1,0)+5,w_coef1[kFX0] - fitrangewave[kxmin])   //5/31/06  added to keep fitting range within bounds of image...Setting coordinates of a box within which to fit, really just equal to initial guess + calculated kfitrange, or edge of image
                  xrmax = min(dimoffset(RotatedOrig1,0) + dimsize(RotatedOrig1,0)-5, w_coef1[kFX0] + fitrangewave[kxmax])   // 5 seems a reasonable "border" offset
                  yrmin = max(dimoffset(RotatedOrig1,1)+5, w_coef1[kFY0] - fitrangewave[kymin])
                  yrmax = min(dimoffset(RotatedOrig1,1) + dimsize(RotatedOrig1,1)-5, w_coef1[kFY0] + fitrangewave[kymax])
   
   		  	if(hexcp ==1) //only use hexmask for hexcp arrays
   		  		 buildhexmask(xrmin,xrmax, yrmin, yrmax, 0) //calls function defined above to make a mask of ones in the shape of a hexagon, used xrmin, yrmin, and a "buffer" region of 10-- chosen somewhat arbitrarily
   		       endif
   		    // printf "xdim = %f\r", xrmax-xrmin
   			
   				
   		 //   make/o/t/n=2 T_Constraints={"k2< g_xrmax",  "k2 > g_xrmin", "k4> g_yrmin", "k4< g_yrmax"}  //not currently used 
   		    
   		     
   			
   		    //print  T_constraints   
   		   // print g_xrmax                                                                                                                                                                       				
	           V_FitError = 0
	           ErrorReporter = 0
	           
	           //choosing average of roi as background
	       
	               //**** changing roi to work with phase images******
	                xrmin = xrmin
	                 xrmax = xrmax
	                yrmin = yrmin
	                 yrmax = yrmax
	                //***********************************
	       //***for use with fluorescence***
//	                 xrmin = xrmin+15
//	                 xrmax = xrmax-15
//	                yrmin = yrmin+15
//	                 yrmax = yrmax-10
//	                 
	             //********
	       
	       
	           duplicate/o/R=[(xrmin), (xrmax)][(yrmin),( yrmax)] RotatedOrig1 ROIAvg
	           wavestats/q ROIavg
	          
	          //print "Original Background guess =", w_coef1[kFBkgnd]
	           w_coef1[kFbkgnd]= v_avg
	           
	          // print "New Average Background guess = ", w_coef1[kfBkgnd]
	           
	        
	         
	         setdrawenv gstart,gname=rects
	         
	         setdrawenv linefgc=(0,655350,0),linethick=5
	         DrawRect xrmin, yrmin, xrmax, yrmax
	          setdrawenv linefgc=(65535,0,0),linethick=1
	           drawrect g_xrmin, g_yrmin, g_xrmax, g_yrmax
	        setdrawenv gstop

	      	
	       //drawaction endinsert
		


                  if (CCDflag == 1)   // 6/28/06 fit with error bars
	               CurveFit/Q/N Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] /I=1 /W= RotatedOrig1Err
	           elseif(hexmaskflag==1)  	                                 // fit using hexmask  // fit with no weighting
	               if(ifile==f1) //show fits for first file, then go silent (remove /D)
	               
	               	if(hexcp==1) //don't use hexmask if not hexcp array
	               		wave m_particleperimeter
	               		wave m_imagethresh
	               		//some stuff to make sure hexmask is working-- can be removed later
					 appendimage/t m_imagethresh

	               		imagethreshold/i/t=0.5 hexmask
	               		imageanalyzeparticles/q/m=1 stats m_imagethresh
	               		imagethreshold/i/t=1 m_particleperimeter 
	            	   	//newimage/n=thresh m_imagethresh
					
					ModifyImage m_imagethresh ctab= {1,255,BlueRedGreen,0}
					ModifyImage m_imagethresh minRGB=NaN,maxRGB=0
				      
				      //*********Diagnostics******************************************
//				      if(i* YPosts + j == 3)//1199)
//						duplicate/o/r=[(xrmin),(xrmax)][(yrmin),(yrmax)] rotatedorig1 test
//						newimage test
//						//CurveFit/q Gauss2D kwCWave = w_coef1   test /D  // /C=T_Constraints
//					endif				            		
	              		//************************************************************
	              		CurveFit/q/n Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] /D  /m=hexmask // /C=T_Constraints
	              		removeimage m_imagethresh
	              	
	              	else    //shut down the windows of fitting result to increase speed   12/12/13  Yu Shi
	              		 CurveFit/q/n/W=2 Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] /D 
	                   endif
	               else
	               	if(hexcp==1)
	               		CurveFit/q/n/W=2 Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)]   /m=hexmask
	               	else
	               		CurveFit/q/n/W=2 Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)]
	               	endif
	               endif
	           endif
	           if (V_FitError != 0)  //  curvefit returned error?
	               ErrorReporter +=1
	               //print "V_fitError = ", v_fiterror, i,j,ifile
	           endif
	          
	           if ((w_coef1[kFX0] < g_xrmin) || (w_coef1[kFX0] > g_xrmax) || (w_coef1[kFY0] < g_yrmin) || (w_coef1[kFY0] > g_yrmax))  // post center out of range? changed these to global constants, not updated due to fit
	               ErrorReporter += 10
	               //print "Center out of range ", i,j,ifile
	           endif
              
	          if ((abs(w_coef1[kFSigX]) < SigMin) || (abs(w_coef1[kFSigY]) < SigMin) || (abs(w_coef1[kFSigX]) > SigMax) || (abs(w_coef1[kFSigY]) > SigMax))
	               //print "Sig out of range" , i,j,ifile
	               ErrorReporter += 100
	           endif
	          
	          if ((abs(w_coef1[kFCorr]) > CorrMax) )
	               //print "Corr out of range" , i,j,ifile
	               ErrorReporter += 1000
	           endif
	           
	           
	           //REFIT upon bad results
	           counter=0
	           do 
	              //break //testing
	              counter+=1
	              if(ErrorReporter == 0)
	              	break
	               endif
	              V_FitError = 0
	              xrmin += deltaROI
	              xrmax -=deltaROI 
	              yrmin+=deltaROI
	              yrmax -=deltaROI
	              //shifting around roi to try to get this straightened out
	              
	            print "REFITTING", i*yposts + j, ifile
	              //print "Initial Fit Results: ", w_coef1
	               
	           //     if ((abs(w_coef1[kFSigX]) < SigMin) )
	            //    	w_coef1=w_coefinit
	           //     	w_coef1[kfsigx]+=2
	                
	          //      elseif (abs(w_coef1[kFSigY]) < SigMin) 
	           //         w_coef1=w_coefinit
	            //    	w_coef1[kfsigy]+=2
	                	
	              //  elseif (abs(w_coef1[kFSigX]) > SigMax) 
	              //      w_coef1=w_coefinit
	              //  	w_coef1[kfsigx]-=2
	                	
	              //  elseif(abs(w_coef1[kFSigY]) > SigMax)
	              //  	w_coef1=w_coefinit
	              //	w_coef1[kfsigy]-=2 // last ditch attempt to fix a particular problem
	              	
	             //   else
	               // 	w_coef1=w_coefinit
	            //    endif
	             w_coef1=w_coefinit
	             
	             w_coef1[kFBkgnd]=v_avg
	             
                     if (CCDflag == 1)   // 6/28/06 fit with error bars
	                CurveFit/Q/n Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] /I=1 /W= RotatedOrig1Err
	                   errorreporter=0
	               else  // fit with no weighting
	                  if(ifile==f1)                  //
	                  	CurveFit/n/q/W=2 Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)]  /D /m=hexmask 
	                  else
	                   	CurveFit/n/q/W=2 Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)]   /m=hexmask ///c=T_Constraints
	                  endif
	                  errorreporter=0
	               endif
	              
	              //After first refit, recalculate error reporter
	              	if (V_FitError != 0)  //  curvefit returned error?
	             		 ErrorReporter +=1
	             		  //print "V_fitError = ", v_fiterror, i,j,ifile
	        		endif
	          
	         		  if ((w_coef1[kFX0] < g_xrmin) || (w_coef1[kFX0] > g_xrmax) || (w_coef1[kFY0] < g_yrmin) || (w_coef1[kFY0] > g_yrmax))  // post center out of range?
	           		    ErrorReporter += 10
	            		   //print "Center out of range ", i,j,ifile
	         
                 		endif
              
	        		  if ((abs(w_coef1[kFSigX]) < SigMin) || (abs(w_coef1[kFSigY]) < SigMin) || (abs(w_coef1[kFSigX]) > SigMax) || (abs(w_coef1[kFSigY]) > SigMax))
	           		  	  //print "Sig out of range" , i,j,ifile, sigmin, sigmax
	           	   		 ErrorReporter += 100
	         	 	endif
	         		 if ((abs(w_coef1[kFCorr]) > CorrMax) )
	               		//print "Corr out of range" , i,j,ifile
	             	 		ErrorReporter += 1000
	          		 endif
	              	 //if there are no errors after first refit, break
	              	if(ErrorReporter == 0)
	              		print "Successful Refit: Sigx =", w_coef1[kfsigx] ,", Sigy =", w_coef1[kfsigy], ", XO =", w_coef1[Kfx0], ", Y0=", w_coef1[kfy0]
	               		break
	             	 endif
	               
	            	if ((V_FitError != 0) || (abs(w_coef1[kFCorr]) > CorrMax) || ((abs(w_coef1[kFSigX]) < SigMin) || (abs(w_coef1[kFSigY]) < SigMin) || (abs(w_coef1[kFSigX]) > SigMax) || (abs(w_coef1[kFSigY]) > SigMax)) || ((w_coef1[kFX0] < xrmin) || (w_coef1[kFX0] > xrmax) || (w_coef1[kFY0] < yrmin) || (w_coef1[kFY0] > yrmax))) //  curvefit returned error?
	                   		ErrorReporter += 10000
	 	 		  endif
	           
	               
	        //After refitting, check again, print warning, go back to fitting again
	           if (V_FitError != 0)  //  curvefit returned error?
	          	  	//ErrorReporter +=1
	               	print "V_fitError = ", v_fiterror, i,j,ifile
	          endif
	           	if ((w_coef1[kFX0] < g_xrmin) || (w_coef1[kFX0] > g_xrmax) || (w_coef1[kFY0] < g_yrmin) || (w_coef1[kFY0] > g_yrmax))  // post center out of range?
	               	//ErrorReporter += 10
	              	print "Center out of range ", i,j,ifile, w_coef1[kfx0], w_coef1[kfy0]
                 endif
	         	 if ((abs(w_coef1[kFSigX]) < SigMin) || (abs(w_coef1[kFSigY]) < SigMin) || (abs(w_coef1[kFSigX]) > SigMax) || (abs(w_coef1[kFSigY]) > SigMax))
	              	print "Sig out of range" , i,j,ifile, w_coef1[kfsigx], w_coef1[kfsigy]
	              	
	             		 // ErrorReporter += 100
	          	 endif
	         	if ((abs(w_coef1[kFCorr]) > CorrMax) )
	               	print "Corr out of range" , i,j,ifile, w_coef1[kfcorr]
	               	
	              	//ErrorReporter += 1000
	          	endif
	          	if(counter ==  10)
		
			      wave refitwave     	
			      print "Counter = ", counter
			      refitwave[counter2][0]={i*yPosts+j}
			      refitwave[counter2][1]={ifile}
			      counter2+=1
	  	      endif
	      while(counter<=10)
	      //    while(0)     //not using this while trying to "fit with constraints"
// autorefit added 5/19/06
                  if (autorefit == 1)
 // 5/31/06 boundary check added to autorefit as well.
                       xrmin = max(dimoffset(RotatedOrig1,0)+5,w_coef1[kFX0] - kar*fitrangewave[kxmin])   //5/31/06  added to keep fitting range within bounds of image
                       xrmax = min(dimoffset(RotatedOrig1,0) + dimsize(RotatedOrig1,0)-5, w_coef1[kFX0] + kar*fitrangewave[kxmax])   // 5 seems a reasonable "border" offset
                       yrmin = max(dimoffset(RotatedOrig1,1)+5, w_coef1[kFY0] - kar*fitrangewave[kymin])
                       yrmax = min(dimoffset(RotatedOrig1,1) + dimsize(RotatedOrig1,1)-5, w_coef1[kFY0] + kar*fitrangewave[kymax])
                       V_FitError = 0
                       ErrorReporter2 = 0
                       w_coefInit = w_coef1
                       if (CCDflag == 1)   // 6/28/06 fit with error bars
	                    CurveFit/Q/N Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] /I=1 /W= RotatedOrig1Err/D 
	                 else                       //shut down the windows of fitting result to increase speed   12/12/13  Yu Shi
	               		if(ifile==f1)
	               		CurveFit/Q/N/W=2 Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] /D
	               		else
	               		CurveFit/Q/N/W=2 Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] 
	               		endif
	                 endif
	                 if (V_FitError != 0)  //  curvefit returned error?
	                     ErrorReporter2 +=2
	                     print "Refit V_fitError = ", v_fiterror, i,j,ifile
	                 endif
	                 if ((w_coef1[kFX0] < xrmin) || (w_coef1[kFX0] > xrmax) || (w_coef1[kFY0] < yrmin) || (w_coef1[kFY0] > yrmax))  // post center out of range?
	                     ErrorReporter2 += 20
	                     print "Refit Center out of range ", i,j,ifile
                        endif
	                 if ((abs(w_coef1[kFSigX]) < SigMin) || (abs(w_coef1[kFSigY]) < SigMin) || (abs(w_coef1[kFSigX]) > SigMax) || (abs(w_coef1[kFSigY]) > SigMax))
	                     print "Refit Sig out of range" , i,j,ifile
	                     ErrorReporter2 += 200
	                 endif
	                 if ((abs(w_coef1[kFCorr]) > CorrMax) )
	                      print "Refit Corr out of range" , i,j,ifile
	                      ErrorReporter2 += 2000
	                 endif
	                 if (ErrorReporter2 != 0)
	                     V_FitError = 0
	                     xrmin += deltaROI
	                     xrmax -=deltaROI 
	                     yrmin+=deltaROI
	                     yrmax -=deltaROI
	                     w_coef1 = w_coefInit
                            if (CCDflag == 1)   // 6/28/06 fit with error bars
	                         CurveFit/Q/N Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] /I=1 /W= RotatedOrig1Err
	                      else  // fit with no weighting
	                         if(ifile==f1)                    //shut down the windows of fitting result to increase speed   12/12/13  Yu Shi
	                         	 CurveFit/N Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] /D
	                         else
	                     		CurveFit/N Gauss2D kwCWave = w_coef1   RotatedOrig1[(xrmin),(xrmax)][(yrmin),(yrmax)] 
	                          endif
	                      endif
	                      if ((V_FitError != 0) || (abs(w_coef1[kFCorr]) > CorrMax)  || ((abs(w_coef1[kFSigX]) < SigMin) || (abs(w_coef1[kFSigY]) < SigMin)) || ((w_coef1[kFX0] < xrmin) || (w_coef1[kFX0] > xrmax) || (w_coef1[kFY0] < yrmin) || (w_coef1[kFY0] > yrmax))) //returned error?
	                          ErrorReporter2 += 20000
	                      endif
	                 endif

                  endif
             endif
               
               //******************
//               if(i* YPosts + j == 533)
//			print w_coef1[kfx0]
//			print w_coef1[kfy0]
//		   endif
             //********************
             //****** 4/15/13 Ellipse foci calculations should go here, before the results from the 2dgauss fit are recorded********
             //**************************************************************************************  
             
             
                    // record fit results for parameters
	      FRwave[i* YPosts + j][kFbkgnd][ifile]=(w_coef1[kFbkgnd]  + zoffset) /zscale //6/28/06  bggnd
             FRwave[i* YPosts + j][kFAmpl][ifile]=w_coef1[kFAmpl]/zscale // A
	      FRwave[i* YPosts + j][kFX0][ifile]=w_coef1[kFX0]   // x of post center top
	      FRwave[i* YPosts + j][kFSigX][ifile]=w_coef1[kFSigX]   // x of post center top
	      FRwave[i* YPosts + j][kFY0][ifile]=w_coef1[kFY0]   // x of post center top
	      FRwave[i* YPosts + j][kFSigY][ifile]=w_coef1[kFSigY]   // x of post center top
	      FRwave[i* YPosts + j][kFCorr][ifile]=w_coef1[kFCorr]   // x of post center top
	      FRwave[i* YPosts + j][kFChisqG][ifile]=V_Chisq   // Chisq
	      if (counter>=10) //only report those which have failed the refits CMK 6/2/09
	      FRwave[i* YPosts + j][kFerrorFlag][ifile]=ErrorReporter + ErrorReporter2  // 4/4/07  = 1 for fit eror, 10 for (x0,y0) out of range, 11 for both
	      endif
	      if (ErrorReporter + ErrorReporter2 >= 20000)  // All refits failed
	          print "All Refits failed. Error flag: ",ErrorReporter+ ErrorReporter2, i, j, ifile
	      endif
	            //now record error estimates
	       
	      FEwave[i* YPosts + j][kFbkgnd][ifile]=W_sigma[kFbkgnd] /zscale  // bggnd
             FEwave[i* YPosts + j][kFAmpl][ifile]=W_sigma[kFAmpl] /zscale  // A
	      FEwave[i* YPosts + j][kFX0][ifile]=W_sigma[kFX0]   // x0
	      FEwave[i* YPosts + j][kFSigX][ifile]=W_sigma[kFSigX]   // SigmaX
	      FEwave[i* YPosts + j][kFY0][ifile]=W_sigma[kFY0]   // y0
	      FEwave[i* YPosts + j][kFSigY][ifile]=W_sigma[kFSigY]   // SigmaY
	      FEwave[i* YPosts + j][kFCorr][ifile]=W_sigma[kFCorr]   // corr

		drawaction getgroup=rects, delete
       endfor //j
    endFor //i

   endfor  // loop over files
end

//*********************************************************
//CalcInterpPostCenters()

//7/11/14   changed accordingly for centroid method
//
// 4/13/07 this version makes files that record line fit params
// 3/25/07 copied from FitGuidePostsPSF()
// 7/5/06 Bugfix: x vs y was not mistakenly fit as y vs x, and "j" should have been "i" in error recording
// 6/4/06  Error propagation added for linear fits
// 5/24/06   Implements interpolation scheme as alternative to baseimage fit for post center positions
//     data are passed in a PSF fitresults-type file
//Modified to use hexcp flag to determine expected array

function  CalcInterpPostCenters(f1,f2,XTop,YTop,XPosts,YPosts,maskwave,datawave,X0p,Y0p)
    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    wave maskwave		// mask that IDs posts
    wave datawave           // contains results of fits for guideposts and empty posts when called
                                      // will have cell post positions written into it 
                                      //DATAWAVE of x0p and y0p is fitbottom/zeropos/zerow
//    wave errorwave       // fiterrors for guideposts and empty posts when called.
                                      // will have post position errors written into it.
    variable X0p,Y0p      // array indices for x,y position of top of post                                 
    wave W_coef, W_sigma
    variable i, j, ifile,nfitposts,ifd,a,b,c,d,sa,sb,sc,sd,ac1
    variable k //dummy variable
    nvar hexcp
    
    make/O/n=(2*XPosts,2)  vertlineparams   // waves that will hold coefs of linear fits
    make/O/n=(YPosts,2)  horizlineparams
    make/O/n=(2*XPosts,2)  vertlineerrs   // waves that will hold error of linear fit coefs  /6/4/06
    make/O/n=(YPosts,2)  horizlineerrs
    make/O/N=2 W_coef
    make/O/N=2 W_sigma
    
    // 4/13/07
    make/O/n=(2*Xposts,2,f2-f1+1) vertlineall
    make/O/n=(2*Xposts,2,f2-f1+1) vertlineallErr
    make/O/n=(Yposts,2,f2-f1+1) horizlineall
    make/O/n=(Yposts,2,f2-f1+1) horizlineallErr
    
   for (ifile=f1; ifile <=f2; ifile +=1)   // loop over images
        // do horizontal fits first
        for (j=0; j < YPosts; j+=1)  // loop over rows of posts
            // count guideposts in this row
            nfitposts =0
            for (i = 0; i < XPosts; i +=1)
                if ((maskwave[i][j] == kEmpty))  // 5/24/06 include both in fit
                    nfitposts +=1
                endif
            endfor
             if(nfitposts<3)
                	print "Horiz nfitposts =", nfitposts
                	
               endif
            // make and fill data waves for fit to this row
            make/O/N=(nfitposts) xfitdata, yfitdata, yfiterr   //6/4/06  add error wave
            
            ifd = 0
            for (i = 0; i < XPosts; i +=1)
                if ((maskwave[i][j] == kEmpty))  // 5/24/06 include both in fit
                 	  
                 	  
                    xfitdata[ifd] = datawave[i* YPosts + j][X0p][ifile]
                    yfitdata[ifd] = datawave[i* YPosts + j][Y0p][ifile]
//7/11/14                    yfiterr[ifd] = sqrt(errorwave[i* YPosts + j][X0p][ifile]^2 + errorwave[i* YPosts + j][Y0p][ifile]^2 )  // add x,y errors in quadrature. (see Bevington)
                    ifd +=1
                    
             
                  endif
            endfor
            // do Linear fit to y vs x, and transfer fit coefs to horizlineparams
            // note - fit is w[0] + w[1]*x
           //7/11/14  
            curvefit/Q Line yfitdata /X=xfitdata /I=1//7/11/14 /W=yfiterr           // 6/4/06  weighting by error bars
           
            horizlineparams[j][0] = W_coef[0]
            horizlineparams[j][1] = W_coef[1]
            horizlineerrs[j][0] = W_Sigma[0]        // 6/4/06 record error bars on fit params
            horizlineerrs[j][1] = W_Sigma[1]
            // 4/13/07
            horizlineall[j][0][ifile] = W_coef[0]
            horizlineall[j][1][ifile] = W_coef[1]
            horizlineallErr[j][0][ifile] = W_Sigma[0]        
            horizlineallErr[j][1][ifile] = W_Sigma[1]
            
            if(abs(w_coef[1]) > 1 )
       	     //print "horiz slope =", w_coef[1], i, j, ifile
        	    //print "horiz intercept =", w_coef[0]
              endif
        endfor //finishes j, rows of post
      
      
      //CMK rewriting to do 2 different lines-- hex code
      
      //NOTE:  fit column "one"  (i= 0) which is made up of two columns, (k=0, and k=1) 
        // Now do vertical fits 
        k=0
      if(hexcp ==1)
        for (i=0; i <  XPosts; i+=1)  // loop over columns of posts-- each column is actually made up of 2 lines
            // count guideposts in this column
            nfitposts =0
            for (j = 0; j < YPosts; j +=2) //only count half of each column
                if ((maskwave[i][j] == kEmpty))  // 5/24/06 include both in fit
                    nfitposts +=1
                endif
            endfor
              if(nfitposts<3)
          		//  print "vert A nfitposts = " , nfitposts
          	endif
            // make and fill data waves for fit to this column
            make/O/N=(nfitposts) xfitdata, yfitdata, xfiterr   //6/4/06 add error wave
            ifd = 0
            for (j = 0; j < YPosts; j+=2) //requires even number of rows!!!
            
          
                if ((maskwave[i][j] == kEmpty))  // 5/24/06 include both in fit
                    xfitdata[ifd] = datawave[i* YPosts + j][X0p][ifile]
                    yfitdata[ifd] = datawave[i* YPosts + j][Y0p][ifile]
//7/11/14                    xfiterr[ifd] = sqrt(errorwave[i* YPosts + j][X0p][ifile]^2 + errorwave[i* YPosts + j][Y0p][ifile]^2 )  // add x,y errors in quadrature. (see Bevington)
                    ifd += 1
                 endif
            endfor
            // do Linear fit to x vs y(!!), and transfer fit coefs to vertlineparams
            //  note - fitting x vs y is to avoid nearly infinite slopes...
           // 7/5/06 Bugfix: was curvefit/Q Line yfitdata /X=xfitdata /I=1 /W=xfiterr           // 6/4/06  weighting by error bars
           
           curvefit/Q Line xfitdata /X=yfitdata /I=1  //7/11/14 /W=xfiterr           // 6/4/06  weighting by error bars
           
            vertlineparams[k][0] = W_coef[0] //replaced i's by dummy variable k
            vertlineparams[k][1] = W_coef[1]
             // 7/5/06 Bugfix: index was mistakenly j
            vertlineerrs[k][0] = W_Sigma[0]        // 6/4/06 record error bars on fit params
            vertlineerrs[k][1] = W_Sigma[1]
            // 4/13/07
            vertlineall[k][0][ifile] = W_coef[0]
            vertlineall[k][1][ifile] = W_coef[1]
            vertlineallErr[k][0][ifile] = W_Sigma[0]        
            vertlineallErr[k][1][ifile] = W_Sigma[1]
            
            if(abs(w_coef[1]) > 1 )
            	//print "vert a slope =", w_coef[1], i, j, ifile
           	//print "vert a intercept =", w_coef[0]
            endif
            
            k+=1
              
       
            // count guideposts in the "second" column -- starting with post 1
            nfitposts =0
            for (j = 1; j < YPosts; j +=2)
                if ((maskwave[i][j] == kEmpty))  // 5/24/06 include both in fit
                    nfitposts +=1
                endif
            endfor
            
            if(nfitposts<3)
            
            		//print "vertB nfitposts = ", nfitposts
            	endif
            // make and fill data waves for fit to this column
            make/O/N=(nfitposts) xfitdata, yfitdata, xfiterr   //6/4/06 add error wave
            ifd = 0
            for (j = 1; j < YPosts; j +=2) 
                if ((maskwave[i][j] == kEmpty))  // 5/24/06 include both in fit
                    xfitdata[ifd] = datawave[i* YPosts + j][X0p][ifile]
                    yfitdata[ifd] = datawave[i* YPosts + j][Y0p][ifile]
//7/11/14                    xfiterr[ifd] = sqrt(errorwave[i* YPosts + j][X0p][ifile]^2 + errorwave[i* YPosts + j][Y0p][ifile]^2 )  // add x,y errors in quadrature. (see Bevington)
                    ifd += 1
                 endif
            endfor
            // do Linear fit to x vs y(!!), and transfer fit coefs to vertlineparams
            //  note - fitting x vs y is to avoid nearly infinite slopes...
           // 7/5/06 Bugfix: was curvefit/Q Line yfitdata /X=xfitdata /I=1 /W=xfiterr           // 6/4/06  weighting by error bars
            
            curvefit/Q Line xfitdata /X=yfitdata /I=1 //7/11/14 /W=xfiterr           // 6/4/06  weighting by error bars
            
            vertlineparams[k][0] = W_coef[0]
            vertlineparams[k][1] = W_coef[1]
             // 7/5/06 Bugfix: index was mistakenly j
            vertlineerrs[k][0] = W_Sigma[0]        // 6/4/06 record error bars on fit params
            vertlineerrs[k][1] = W_Sigma[1]
            // 4/13/07
            vertlineall[k][0][ifile] = W_coef[0]
            vertlineall[k][1][ifile] = W_coef[1]
            vertlineallErr[k][0][ifile] = W_Sigma[0]        
            vertlineallErr[k][1][ifile] = W_Sigma[1]
            
            if(abs(w_coef[1]) > 1 )
          	  //print "vert b slope =", w_coef[1], i, j, ifile
           	 //print "vert b intercept =", w_coef[0]
            endif
            k+=1
            
        endfor //ends each full column  consisting of two columns of posts
        
        else //if hexcp ==0 and it's a square array
         for (i=0; i <  XPosts; i+=1)  // loop over columns of posts--
           for (j = 1; j < YPosts; j +=1) //
                if ((maskwave[i][j] == kEmpty))  // 5/24/06 include both in fit
                    nfitposts +=1
                endif
            endfor
            
            if(nfitposts<3)
            
            		//print "vertB nfitposts = ", nfitposts
            endif
            // make and fill data waves for fit to this column
            make/O/N=(nfitposts) xfitdata, yfitdata, xfiterr   //6/4/06 add error wave
            ifd = 0
            for (j = 1; j < YPosts; j +=1) 
                if ((maskwave[i][j] == kEmpty))  // 5/24/06 include both in fit
                    xfitdata[ifd] = datawave[i* YPosts + j][X0p][ifile]
                    yfitdata[ifd] = datawave[i* YPosts + j][Y0p][ifile]
//7/11/14                    xfiterr[ifd] = sqrt(errorwave[i* YPosts + j][X0p][ifile]^2 + errorwave[i* YPosts + j][Y0p][ifile]^2 )  // add x,y errors in quadrature. (see Bevington)
                    ifd += 1
                 endif
            endfor
            
            // do Linear fit to x vs y(!!), and transfer fit coefs to vertlineparams
            //  note - fitting x vs y is to avoid nearly infinite slopes...
           // 7/5/06 Bugfix: was curvefit/Q Line yfitdata /X=xfitdata /I=1 /W=xfiterr           // 6/4/06  weighting by error bars
            
            curvefit/Q Line xfitdata /X=yfitdata /I=1 //   7/11/14/W=xfiterr           // 6/4/06  weighting by error bars
            
            vertlineparams[k][0] = W_coef[0]
            vertlineparams[k][1] = W_coef[1]
             // 7/5/06 Bugfix: index was mistakenly j
            vertlineerrs[k][0] = W_Sigma[0]        // 6/4/06 record error bars on fit params
            vertlineerrs[k][1] = W_Sigma[1]
            // 4/13/07
            vertlineall[k][0][ifile] = W_coef[0]
            vertlineall[k][1][ifile] = W_coef[1]
            vertlineallErr[k][0][ifile] = W_Sigma[0]        
            vertlineallErr[k][1][ifile] = W_Sigma[1]
            
             if(abs(w_coef[1]) > 1 )
          	  //print "vert b slope =", w_coef[1], i, j, ifile
           	 //print "vert b intercept =", w_coef[0]
            endif
            k+=1
            
             
       endfor
            if(abs(w_coef[1]) > 1 )
          	  //print "vert b slope =", w_coef[1], i, j, ifile
           	 //print "vert b intercept =", w_coef[0]
            endif
            k+=1
            
  	endif
              // Compute centers of cell posts as intersection of these fit lines
      
       k=0
        for (i = 0; i <  XPosts; i +=1)
         if(hexcp==1)
            for (j = 0; j < YPosts; j +=2)
                    a = horizlineparams[j][1]    // See calc DHR Book VII, p. 39
                    b = horizlineparams[j][0]    //  that calc uses y = ax+b
                    c = vertlineparams[k][1]       // and x = cy+d
                    d = vertlineparams[k][0]
                    sa = horizlineerrs[j][1]
                    sb = horizlineerrs[j][0]
                    sc = vertlineerrs[k][1]
                    sd = vertlineerrs[k][0]
                    datawave[i* YPosts + j][X0p][ifile] = (c*b +d)/(1 - a*c)
                    datawave[i* YPosts + j][Y0p][ifile]  = (a*d + b)/(1- a*c)
                    // 6/4/06 error propagation - see DHR Book VII, p. 74
                    ac1 = 1/(1-a*c)^2
//7/11/14                    errorwave[i* YPosts + j][X0p][ifile] = ac1*(ac1*c^2*(c*b+d)^2*sa^2 + c^2*sb^2 +ac1*(a*d+b)^2*sc^2 + sd^2)
//                    errorwave[i* YPosts + j][Y0p][ifile] = ac1*(ac1*(d+b*c)^2*sa^2 + sb^2 + ac1*a^2*(a*d+b)^2*sc^2 + a^2*sd^2)                    
            
           
            endfor
           k+= 1  //goes to second line of each column
    
       
          	  for (j = 1; j < YPosts; j +=2)
                    a = horizlineparams[j][1]    // See calc DHR Book VII, p. 39
                    b = horizlineparams[j][0]    //  that calc uses y = ax+b
                    c = vertlineparams[k][1]       // and x = cy+d
                    d = vertlineparams[k][0]
                    sa = horizlineerrs[j][1]
                    sb = horizlineerrs[j][0]
                    sc = vertlineerrs[k][1]
                    sd = vertlineerrs[k][0]
                    datawave[i* YPosts + j][X0p][ifile] = (c*b +d)/(1 - a*c)
                    datawave[i* YPosts + j][Y0p][ifile]  = (a*d + b)/(1- a*c)
                    // 6/4/06 error propagation - see DHR Book VII, p. 74
                    ac1 = 1/(1-a*c)^2
//7/11/14                    errorwave[i* YPosts + j][X0p][ifile] = ac1*(ac1*c^2*(c*b+d)^2*sa^2 + c^2*sb^2 +ac1*(a*d+b)^2*sc^2 + sd^2)
//7/11/14                    errorwave[i* YPosts + j][Y0p][ifile] = ac1*(ac1*(d+b*c)^2*sa^2 + sb^2 + ac1*a^2*(a*d+b)^2*sc^2 + a^2*sd^2)                    
        	    endfor
             k += 1
             

            
       
       else

       
       	 for (j = 1; j < YPosts; j +=1)
                    a = horizlineparams[j][1]    // See calc DHR Book VII, p. 39
                    b = horizlineparams[j][0]    //  that calc uses y = ax+b
                    c = vertlineparams[k][1]       // and x = cy+d
                    d = vertlineparams[k][0]
                    sa = horizlineerrs[j][1]
                    sb = horizlineerrs[j][0]
                    sc = vertlineerrs[k][1]
                    sd = vertlineerrs[k][0]
                    datawave[i* YPosts + j][X0p][ifile] = (c*b +d)/(1 - a*c)
                    datawave[i* YPosts + j][Y0p][ifile]  = (a*d + b)/(1- a*c)
                    // 6/4/06 error propagation - see DHR Book VII, p. 74
                    ac1 = 1/(1-a*c)^2
//7/11/14                    errorwave[i* YPosts + j][X0p][ifile] = ac1*(ac1*c^2*(c*b+d)^2*sa^2 + c^2*sb^2 +ac1*(a*d+b)^2*sc^2 + sd^2)
//                    errorwave[i* YPosts + j][Y0p][ifile] = ac1*(ac1*(d+b*c)^2*sa^2 + sb^2 + ac1*a^2*(a*d+b)^2*sc^2 + a^2*sd^2)                    
            endfor
          k += 1
    
       endif
       
       endfor //i
    endfor   // ifile

end
 



//*********************************************************
// CalcContractility
// 1/24/06
// Calculates various sums
// displacement wave:   sum(xi), sum(yi) to check Newton's Law....
//                                    sqrt(sum(xi^2 + yi^2) is contractility
//  difference from initial frame:    deltar sqrt(sum((xi-xi0)^2 + (yi-yi0)^2)
// makes a wave contractility[Nimages][8]   =(sumxi,sumyi,contr,deltar)_cell,(sumxi,sumyi,contr,deltar)_noncellposts

function CalcContractility(f1,f2,XPosts,YPosts,maskwave,dispwave,sdwave)
   variable f1,f2  	// beginning and ending file numbers
    variable XPosts,YPosts  // number of posts in x,y in ROI
    wave maskwave		// mask that IDs posts
    wave dispwave      // displacement data
    wave sdwave  // shifted position data
    
    variable i,j,k,ifile,ncp,nbp,wlen
    wlen = f2-f1 +1
    make/O/n=(wlen,8) contractility
    SetScale/I x f1,f2,"", contractility                   
    contractility = 0
    
    for (ifile=f1; ifile <=f2; ifile +=1)   // loop over images
        ncp = 0
        nbp = 0
        for (i=0;i<=(Xposts-1); i+=1)    // loop over posts in this image
	     for (j=0;j<=(YPosts-1); j+=1)
                if (maskwave[i][j] != kIgnore)
                   if (maskwave[i][j] == kCell)     // posts under cell
                       contractility[ifile][0] += dispwave[i* YPosts + j][kFX0][ifile]
                       contractility[ifile][1] += dispwave[i* YPosts + j][kFY0][ifile]
                       contractility[ifile][2] += dispwave[i* YPosts + j][kFX0][ifile]^2 + dispwave[i* YPosts + j][kFY0][ifile]^2
                       contractility[ifile][3] += (sdwave[i* YPosts + j][kFX0][ifile]- sdwave[i* YPosts + j][kFX0][f1])^2 + (sdwave[i* YPosts + j][kFY0][ifile]- sdwave[i* YPosts + j][kFY0][f1])^2
                       ncp +=1
                    elseif  ((maskwave[i][j] == kEmpty))   // posts for background calc
                       contractility[ifile][4] += dispwave[i* YPosts + j][kFX0][ifile]
                       contractility[ifile][5] += dispwave[i* YPosts + j][kFY0][ifile]
                       contractility[ifile][6] += dispwave[i* YPosts + j][kFX0][ifile]^2 + dispwave[i* YPosts + j][kFY0][ifile]^2
                       contractility[ifile][7] += (sdwave[i* YPosts + j][kFX0][ifile]- sdwave[i* YPosts + j][kFX0][f1])^2 + (sdwave[i* YPosts + j][kFY0][ifile]- sdwave[i* YPosts + j][kFY0][f1])^2
                       nbp +=1
                    endif
                endif
            endfor
        endfor
        contractility[ifile][0] /= ncp
        contractility[ifile][1] /= ncp
        contractility[ifile][2] = sqrt(contractility[ifile][2]/ncp)
        contractility[ifile][3] = sqrt(contractility[ifile][3]/ncp)
        contractility[ifile][4] /= nbp
        contractility[ifile][5] /= nbp
        contractility[ifile][6] = sqrt(contractility[ifile][6]/nbp)
        contractility[ifile][7] = sqrt(contractility[ifile][7]/nbp)
    endfor
end

   
//****************************************
// calcShiftsGrelF1rot

// 4/15/07 calculation of rotation added
//     includes both guide and empty posts in shift calculation

// 3/25/07 Calculates shifts relative to first frame

// 6/28/06  error weighting added (see Bevington p.70), and error of shift calculated
//   shift error stored in columns 5 and 6 of shiftwave

//5/26/06  Calculates shifts relative to results from base fit - probably first RHOD image in the gaussian fit scheme
//1/19/06  calculates (does not apply) shifts from results of a fit using the
// guide posts specified in maskwave

function  calcShiftsGrelF1rot(f1,f2,XPosts,YPosts,maskwave,shiftwave,fitres)
    variable f1,f2  	// beginning and ending file numbers
    variable XPosts,YPosts  // number of posts in x,y in ROI
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
    wave fitres                 // results from fit
//    wave fiterr                 //6/28/06 errors in fit
    
    variable i,j, ifile
    variable sumx,sumy,sumxs,sumys,sigx2,sigy2, nposts
    variable dx,dy,dxb,dyb
    for (ifile=f1; ifile <=f2; ifile +=1)   // loop over images
        sumx = 0; sumy = 0; sumxs = 0; sumys = 0; nposts = 0
         for (i=0;i<=(Xposts-1); i+=1)    // loop over posts in this image
	    for (j=0;j<=(YPosts-1); j+=1)
               if ((maskwave[i][j] == kEmpty))
                   sumx +=  fitres[i* YPosts + j][kFX0][ifile] 
                   sumy += fitres[i* YPosts + j][kFY0][ifile] 
//7/11/14                   sumxs += fiterr[i* YPosts + j][kFX0][ifile]^2
//7/11/14                   sumys += fiterr[i* YPosts + j][kFY0][ifile]^2
                   nposts +=1
               endif
          endfor
       endfor
       // write centers of arrays into shiftwave
       shiftwave[ifile][kiXc] = sumx/nposts
       shiftwave[ifile][kiYc] = sumy/nposts
       // write errors-of-mean of centers of arrays into shiftwave
       shiftwave[ifile][kisigXc] = sqrt(sumxs)/nposts
       shiftwave[ifile][kisigYc] = sqrt(sumys)/nposts
       // compute shifts and write to shiftwave
      if (ifile ==f1)
           shiftwave[ifile][kiXsf] = 0
           shiftwave[ifile][kiYsf] = 0
           shiftwave[ifile][kisigXsf] = 0   // error bars
           shiftwave[ifile][kisigYsf] = 0
       else
           shiftwave[ifile][kiXsf] = shiftwave[ifile][kiXc] - shiftwave[f1][kiXc] 
           shiftwave[ifile][kiYsf] = shiftwave[ifile][kiYc] - shiftwave[f1][kiYc] 
           shiftwave[ifile][kisigXsf] = sqrt(shiftwave[ifile][kisigXc]^2 + shiftwave[f1][kisigXc]^2)    // error bars
           shiftwave[ifile][kisigYsf] = sqrt(shiftwave[ifile][kisigYc]^2 + shiftwave[f1][kisigYc]^2)
       endif 
   endfor
   // now compute rotation angle for each frame
    for (ifile=f1; ifile <=f2; ifile +=1)   // loop over images
        sumx = 0; sumy = 0; sumxs = 0; sumys = 0; nposts = 0
         for (i=0;i<=(Xposts-1); i+=1)    // loop over posts in this image
	    for (j=0;j<=(YPosts-1); j+=1)
               if ((maskwave[i][j] == kEmpty))
                   dx = fitres[i* YPosts + j][kFX0][ifile] - shiftwave[ifile][kiXc]
                   dy = fitres[i* YPosts + j][kFY0][ifile] - shiftwave[ifile][kiYc]
                   dxb = fitres[i* YPosts + j][kFX0][f1] - shiftwave[f1][kiXc]
                   dyb = fitres[i* YPosts + j][kFY0][f1] - shiftwave[f1][kiYc]
                   sumy += dx*dyb - dy*dxb
                   sumx += dx*dxb + dy*dyb
               endif
          endfor
       endfor
       // compute angle
       shiftwave[ifile][kitheta] = atan(sumy/sumx)
    endfor   
end


// ApplyShiftsG

//7/7/06  propagates errorbars.  - replaces ApplyShifts()
// used by plotalot

function  ApplyShiftsG(f1,f2,XPosts,YPosts,shiftwave,datawave,errorwave)

    variable f1,f2  	// beginning and ending file numbers
    variable XPosts,YPosts  // number of posts in x,y in ROI
    wave shiftwave           // shifts for each frame, angle etc
    wave datawave                 // copy of results from fit to be shifted
    wave errorwave                 // copy of errors from fit to be shifted

    variable i,j, Nguideposts,valx,valy,ifile
    
    for (ifile=f1; ifile <=f2; ifile +=1)   // loop over images
        for (i=0;i<=(Xposts-1); i+=1)    // loop over posts in this image
	    for (j=0;j<=(YPosts-1); j+=1)
                   datawave[i* YPosts + j][kFX0][ifile] -= shiftwave[ifile][kiXsf]
                   datawave[i* YPosts + j][kFY0][ifile] -= shiftwave[ifile][kiYsf]
                   errorwave[i* YPosts + j][kFX0][ifile]  = sqrt( errorwave[i* YPosts + j][kFX0][ifile]^2 + shiftwave[ifile][kisigXsf]^2)
                   errorwave[i* YPosts + j][kFY0][ifile]  = sqrt( errorwave[i* YPosts + j][kFY0][ifile]^2 + shiftwave[ifile][kisigYsf]^2)
          endfor
     endfor
  endfor
end


// calcdispGridEach

// 3/25/07
// used in plotalot

function  calcdispGridEach(f1,f2,XPosts,YPosts,datawave,errorwave,zerow,zeroerr)

    variable f1,f2  	// beginning and ending file numbers
    variable XPosts,YPosts  // number of posts in x,y in ROI
    wave datawave                 // copy of results from fit to be shifted
    wave errorwave                 // copy of errors from fit to be shifted
    wave zerow       // contains gridded undeflected post positions
    wave zeroerr    // contains
   variable npost
   svar RootString
   
    variable i,j, Nguideposts,valx,valy,ifile
    
    for (ifile=f1; ifile <=f2; ifile +=1)   // loop over images
        for (i=0;i<=(Xposts-1); i+=1)    // loop over posts in this image
	    for (j=0;j<=(YPosts-1); j+=1)
                   datawave[i* YPosts + j][kFX0][ifile] -= zerow[i* YPosts + j][kFX0][ifile]
                   datawave[i* YPosts + j][kFY0][ifile] -= zerow[i* YPosts + j][kFY0][ifile]
//                   if(datawave[i*yposts+j][kfx0][ifile] >= 100)
//                   	npost=i*yposts+j
//                   	printf "Error: Post %d has an x displacement over 20 pix\r", i*yposts+j
//                   endif
//                   if(datawave[i*yposts+j][kfy0][ifile] >= 100)
//                         npost=i*yposts+j
//                   	printf "Error: Post %d has a displacement over 20 pix\r", i*yposts+j
//                   endif
                   errorwave[i* YPosts + j][kFX0][ifile]  = sqrt( errorwave[i* YPosts + j][kFX0][ifile]^2 + zeroerr[i* YPosts + j][kFX0][ifile]^2)
                   errorwave[i* YPosts + j][kFY0][ifile]  = sqrt( errorwave[i* YPosts + j][kFY0][ifile]^2 + zeroerr[i* YPosts + j][kFY0][ifile]^2)
          endfor
     endfor
  endfor

	

end


// ----------- ApplyRotation --------
//4/15/07  applies rotation to displacement wave

function  ApplyRotation(f1,f2,XPosts,YPosts,datawave,errorwave,shiftwave)
    variable f1,f2  	// beginning and ending file numbers
    variable XPosts,YPosts  // number of posts in x,y in ROI
    wave datawave                 // copy of results from fit to be shifted
    wave errorwave                 // copy of errors from fit to be shifted
    wave shiftwave
    
    variable i,j, theta,ifile,dxrot,dyrot
    ///CMK removing apply rotation, not relevant
    for (ifile=f1; ifile <=f2; ifile +=1)   // loop over images
        for (i=0;i<=(Xposts-1); i+=1)    // loop over posts in this image
	 for (j=0;j<=(YPosts-1); j+=1)
               dxrot = cos(theta)*datawave[i* YPosts + j][kFX0][ifile]  - sin(theta)*datawave[i* YPosts + j][kFY0][ifile] 
               dyrot = sin(theta)*datawave[i* YPosts + j][kFX0][ifile]  + cos(theta)*datawave[i* YPosts + j][kFY0][ifile] 
               datawave[i* YPosts + j][kFX0][ifile] = dxrot
               datawave[i* YPosts + j][kFY0][ifile] = dyrot
               // now do errors
               dxrot = cos(theta)^2*errorwave[i* YPosts + j][kFX0][ifile]^2 + sin(theta)^2*errorwave[i* YPosts + j][kFY0][ifile] ^2
               dyrot = sin(theta)^2*errorwave[i* YPosts + j][kFX0][ifile]^2  + cos(theta)^2*errorwave[i* YPosts + j][kFY0][ifile]^2
              errorwave[i* YPosts + j][kFX0][ifile] = sqrt(dxrot)
               errorwave[i* YPosts + j][kFY0][ifile] = sqrt(dyrot)
          endfor
      endfor
   endfor
end

//---------------Function RefitBadlyFitPosts------------------
//CMK 7-16-09
// 12/10/13 modified to read from movie

//Posts which cannot be refit with 10 tries are stored in "Refitwave" then this function runs "updategfit" on them
//
Function RefitBadlyFitPosts(ctrlname): ButtonControl


string ctrlname
wave refitwave
variable postnum
variable filenum

	wave ASTR,W_Sigma
	NVAR PInd, F1_Ind, F2_Ind 
      NVAR CCDflag   //6/28/06
	Variable npost=PInd, f1=F1_Ind, f2=F2_Ind, ifile, AngleofRotation=0
	SVAR FitRes_String, FitErr_String, ShiftString
	SVAR MaskString,  AttributeString
	Duplicate/O $FitRes_String FRS
	Duplicate/O $FitErr_String FRS_Err
	Duplicate/O $ShiftString ShiftWave
	//Prompt npost, "N of post"
	//Prompt f1, "Start File #"
	//Prompt f2, "End File #"
	//DoPrompt "File info", npost, f1, f2
	//if (V_flag==1)
	//	return 0 
	//endif 

     // ControlInfo CCDErrors
	 //   if (V_Value == 1)
	//        CCDflag = 1  // use ccd error weighting
	  //  else
	        CCDflag = 0
	//    endif	
	// temporary - 6/28/06
	//if (CCDflag == 1)
	 //   Print "ChangeGfit not implemented for CCD errors"
	//    return(0)
	//endif
	
	//first has to read refitwave
	
	variable k
	variable endpost = dimsize(refitwave,0)
	
	for(k = 0; k<=endpost;k+=1)
	
		
		npost = refitwave[k][0]
		
		if(refitwave[k][1]-1>0)
			f1= refitwave[k][1]-1 //read refitwave and fit the post starting with the file before the one where it went bad
			f2 = refitwave[k][1]
			
				
			variable v0=FRS[npost][kFBkgnd][f1], 	v1=FRS[npost][kFAmpl][f1],  v2=FRS[npost][kFX0][0] //changed f1 to zero in v4, v2
			variable v3=FRS[npost][kFSigX][f1], v4=FRS[npost][kFY0][0], 	v5=FRS[npost][kFSigY][f1]
			variable v6=FRS[npost][kFCorr][f1]
			variable Fctr=0.3
			variable MulLTX=Fctr, MulLTY=Fctr, MulRBX=Fctr, MulRBY=Fctr   //Multiplication factors for ROI
			
			//Prompt v0, "~Background";	Prompt v1, "~Amplitude"
			//Prompt v2, "X0";	Prompt v3, "SigX"; Prompt v4, "Y0"
			//Prompt v5, "SigY";	Prompt v6, "Corr";	
			//DoPrompt "Take a look at the parameters", v0, v1,v2, v3, v4, v5, v6
			//if (V_flag==1)
			//	return 0
			//endif
			FRS[npost][kFBkgnd][f1]=v0; 	FRS[npost][kFAmpl][f1]=v1;		FRS[npost][kFX0][f1]=v2;   //FRS[npost][7][f1]=V_chisq;	//changed f1 to zero in v4, v2
			FRS[npost][kFSigX][f1]=v3;	FRS[npost][kFY0][f1]=v4; 	FRS[npost][kFSigY][f1]=v5	;	FRS[npost][kFCorr][f1]=v6; 
		
			SaveGfit()
			//Prompt MulLTX, "LTX"	;Prompt MulLTY, "LTY"
			//Prompt MulRBX, "RBX";	Prompt  MulRBY, "RBY" 
			//DoPrompt "Multiplication factors ROI, Lattice "+num2str(ASTR[5]), MulLTX, MulLTY, MulRBX, MulRBY
			//if (V_flag==1)
			//	return 0
			//endif
			string s1
// new 12/10/13 for movie input
                  NVAR InputMovieID_G, NmovieFrames_G   //12/10/13
                  WAVE M_MovieFrame,M_rgb2gray
                  PlayMovieAction setFrontMovie = inputmovieID_G
                  PlayMovieAction stop, frame = f1                // go to first frame to be analyzed
                  for (ifile=f1; ifile <=f2; ifile +=1)   // loop over raw data files
                       PlayMovieAction extract
                       if (ifile < f2)
                          PlayMovieAction step = 1   //  move to next frame to be ready for next pass through loop
                        endif


// 12/10/13        sprintf s1, "%s%03d.tif", ksRhodFileRoot,ifile   // 4/3/07
                        AngleOfRotation = shiftwave[ifile][kiAngle]
//  12/10/13   	ImageLoad/P=MovieFolder/T=tiff/O/N=image1 s1    // load, copy, and rotate the image
                        Duplicate/O image1 OrigImage1

		       	ImageRotate/A=(AngleOfRotation) OrigImage1
		       	Duplicate/O M_RotatedImage RotatedOrig1
		       	DoWindow/K RotatedImage
		       	NewImage/K=1/N=RotatedImage RotatedOrig1
		      		ModifyImage RotatedOrig1 ctab= {AsTR[6],AStr[7],Grays,0}
		
		
				//NOTE: Astr[5] is lattice ct. in pixels
		      		SetAxis/R left (v4+ASTR[5]+ShiftWave[ifile][11]),(v4-ASTR[5]+ShiftWave[ifile][11]);  // Fixed ROI, no shift!!! //changed this to the shifts put in by hand for cases when the calculated shift is wrong, due to a misfit
			      SetAxis top (v2-ASTR[5]+ShiftWave[ifile][10]),(v2+ASTR[5]+ShiftWave[ifile][10])
			      	TextBox/C/N=text0/A=RT/X=1.00/Y=1.00 num2str(ifile)
			      	TextBox/C/N=text1/A=RT/X=5.00/Y=1.00 num2str(npost)
			      	
		       ///------------------------------------
				Make/O/N=7 Inits
				Inits[kFBkgnd]=FRS[npost][kFBkgnd][ifile]; 		Inits[kFAmpl]=FRS[npost][kFAmpl][ifile];		
				Inits[kFX0]=FRS[npost][kFX0][ifile];   	Inits[kFSigX]=FRS[npost][kFSigX][ifile];		
				Inits[kFY0]=FRS[npost][kFY0][ifile]; 		Inits[kFSigY]=FRS[npost][kFSigY][ifile];		
				Inits[kFCorr]=FRS[npost][kFCorr][ifile]; 	
				//FuncFitMD gauss2dquad w_coef  RotatedOrig1[pcsr(A),pcsr(B)][qcsr(A),qcsr(B)]/D 
				// For Gauss fit v3 was replaced with v2	
				//if(ifile==f1)
					CurveFit/Q Gauss2d  kwCWave=Inits RotatedOrig1[(v2-ASTR[5]*MulLTX+shiftwave[ifile][10]),(v2+MulRBX*ASTR[5])+shiftwave[ifile][10]][(v4-ASTR[5]*MulLTY+shiftwave[ifile][11]),(v4+ASTR[5]*MulRBY+shiftwave[ifile][11])] /D ///m=hexmask
				//else
					//CurveFit/Q Gauss2D  kwCWave=Inits RotatedOrig1[(v2-ASTR[5]*MulLTX),(v2+MulRBX*ASTR[5])][(v4-ASTR[5]*MulLTY),(v4+ASTR[5]*MulRBY)] /D ///m=hexmask
				//endif
					FRS[npost][kFBkgnd][ifile]=Inits[kFBkgnd]; 	FRS[npost][kFAmpl][ifile]=Inits[kFAmpl];	
					FRS[npost][kFX0][ifile]=Inits[kFX0];   FRS[npost][kFSigX][ifile]=Inits[kFSigX];	
					FRS[npost][kFY0][ifile]=Inits[kFY0]; 	FRS[npost][kFSigY][ifile]=Inits[kFSigY];		
					FRS[npost][kFCorr][ifile]=Inits[kFCorr]; 
					FRS[npost][kFChisqG][ifile]=V_chisq;	
			            //now record error estimates
			             FRS_Err[npost][kFbkgnd][ifile]=W_sigma[kFbkgnd]   // bggnd
			             FRS_Err[npost][kFAmpl][ifile]=W_sigma[kFAmpl]   // A
			             FRS_Err[npost][kFX0][ifile]=W_sigma[kFX0]   // x0
			             FRS_Err[npost][kFSigX][ifile]=W_sigma[kFSigX]   // SigmaX
			             FRS_Err[npost][kFY0][ifile]=W_sigma[kFY0]   // y0
			             FRS_Err[npost][kFSigY][ifile]=W_sigma[kFSigY]   // SigmaY
			             FRS_Err[npost][kFCorr][ifile]=W_sigma[kFCorr]   // corr
					SaveGfit()
				 
				endfor
			else
			print "I can't fit the -1th frame"
			endif
		
	endfor
end



function  bpass(lnoise, lobject, imagewave)
wave imagewave
variable lnoise, lobject

variable delta
variable thresh

variable i,j
variable nx, ny


nx=dimsize(imagewave, 0)
ny=dimsize(imagewave, 1)

duplicate/o imagewave mgaussian
duplicate/o imagewave mboxcar

matrixfilter/N=(lnoise)/P=1  avg mboxcar

matrixfilter/N=(lobject)/P=1  gauss mgaussian

make/O/B/U/N=(nx,ny) mfinal

//make/O/N=(nx,ny) mfinalt

//mfinalt=mgaussian-mboxcar

//thresh=0.1*wavemax(mfinalt)
//temp trial by setting thresh to zero
thresh = 0
//killwaves mfinalt


for(i=0; i<nx; i+=1)
    for(j=0; j<ny; j+=1)
        delta=mgaussian[i][j]-mboxcar[i][j]
        if (delta>thresh)
           mfinal[i][j]=delta
           else 
           mfinal[i][j]=0
        endif
        
   endfor
endfor

end

//added on 6/16/2014, use centroid method to find posts center
        
        
function centroidgfit(f1,f2,XTop,YTop,XPosts,YPosts,zmin, zmax,maskwave,shiftwave,InitGuessWave)
    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    variable zmin,zmax      // min and max grayscale levels for display during fitting
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
    wave InitGuessWave   // Initial guesses for fit to first post.   These will be copied into w_coef1[]
                                        //  Should have same dimensions an definitions as w_coef1.  As of 2/10/06 this is length = 7
                                        //  Initial guesses for X0 and Y0 (post center coords) and post bottoms need not be supplied as initial guesses
                                        //  Initial guesses needed:
                                        //  InitGuessWave[0] =   Background
                                        //  InitGuessWave[1] =   Ampl
                                         //  InitGuessWave[2] =  Post Radius (in um).  This is a fixed value, presumably determined 
                                         // by a prior call to PSFfitForShifts followed by suitable averaging.
//    string pathstr            // path to data files

   make/O/n=(XPosts*YPosts,2,f2-f1+1) fitresults // always created. calling routine must copy it as needed
   SetScale/I z f1,f2,"", fitresults    
  // make/O/n=(XPosts*YPosts,k2,f2-f1+1) fiterrors // always created. two less columns - no chisq column or errorflag column needed
//   SetScale/I z f1,f2,"", fiterrors 
 //  fiterrors=0
 //  make/O/N=4 fitrangewave

    docentoridgfit(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,maskwave,shiftwave, fitresults)

end   
       
       
function docentoridgfit(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,maskwave, shiftwave, FRwave)      
    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    variable zmin,zmax      // min and max grayscale levels for display during fitting
//    variable autorefit          // 5/19/06 = 1 for autorefit of each post.
//    wave fitrangewave        // 3/8/06 has fitrange increments - default = kfitrange
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
//    wave InitGuessWave   // Initial guesses for fit to first post.   These will be copied into w_coef1[]
                                        //  Should have same dimensions an definitions as w_coef1. 
                                        //  Initial guesses for X0 and Y0 (post center coords) need not be supplied 
       					//10/15/08 CMK adding in initial guesses for xo, yo	
       					//  Initial guesses needed:
                                        //  InitGuessWave[0] =   Background
                                        //  InitGuessWave[1] =   Ampl
  wave FRwave			// fit results wave.  Will be updated appropriately
//  wave FEwave			// fit errors wave.  Will be updated appropriately
//  string pathstr            // path to data files
  
//   variable xrmin,xrmax,yrmin,yrmax   // 5/19/06  new for refit range
//   string s1
   variable i, j, ifile, fitthispost
   variable ipost
   variable xcenter, ycenter, Xsp, Ysp
   variable V_chisq
   variable AngleOfRotation
   variable dX,dY   // local variables for lattice constants of post array
   variable centroidrad=7  // radius of centorid mask\
   variable xcmax, xcmin, ycmax, ycmin  // max and min coord in x and y direction
   variable ix, iy
   variable centx, centy, centval
   variable xdis, ydis
   variable epsx, epsy
   variable xoffset
   variable postrad = 2*centroidrad +1

 	wave xt2, yt2 //duplicates of the initial 3 locations used in roiselect, here used to determine angles
 
//    variable zscale = 1 // 6/28/06 for ccd
//    variable zoffset = 0
    variable zmin1 =zmin
    variable  zmax1 = zmax
    variable slopey =  0//for use instead of rotate image
    variable slopex = 0
    variable counter =0
    variable counter2=0
    
    variable xroimin, xroimax, yroimin, yroimax
    
    xroimin= max(floor(XT2[0]-50), 0)
    xroimax= ceil(XT2[1]+50) 
    yroimin = max(floor(YT2[0]-50), 0)
    yroimax = ceil(YT2[2]+50)
  
  
    	slopex=(YT2[1]-YT2[0])/(XT2[1]-XT2[0])
      if(mod(yposts,2)==0)
      slopey =  (XT2[2]-0.5*shiftwave[0][kiA]-XT2[0])/(YT2[2]-YT2[0])
      else
      slopey =  (XT2[2]-XT2[0])/(YT2[2]-YT2[0]) 
      endif
   
   
   // new 12/10/13 for movie input
   NVAR InputMovieID_G, NmovieFrames_G   //12/10/13
   WAVE M_MovieFrame,M_rgb2gray
   PlayMovieAction setFrontMovie = inputmovieID_G
   PlayMovieAction stop, frame = f1                // go to first frame to be analyzed
   for (ifile=f1; ifile <=f2; ifile +=1)   // loop over raw data files
       PlayMovieAction extract
       if (ifile < f2)
          PlayMovieAction step = 1   //  move to next frame to be ready for next pass through loop
       endif
// new 12/10/13 for movie input
       ImageTransform rgb2gray M_MovieFrame   // convert to grayscale in M_rgb2gray
       Duplicate/O M_rgb2gray Image1

// 12/10/13        sprintf s1, "%s%03d.tif", ksRhodFileRoot,ifile   // 4/3/07
       AngleOfRotation = shiftwave[ifile][kiAngle]
//  12/10/13   	ImageLoad/P=MovieFolder/T=tiff/O/N=image1 s1    // load, copy, and rotate the image
       Duplicate/O image1 OrigImage1
//       if (CCDflag ==1) // use proper statistical errors from CCD  6/28/06
//          origImage1 = OrigImage1 - CCDBias(x,y)
//       endif
       ImageRotate/A=(AngleOfRotation) OrigImage1
       Duplicate/O M_RotatedImage RotatedOrig1
       
       DoWindow/K RotatedImage
//       dowindow/k thresh //for testing
//       newimage/n=thresh m_imagethresh
       NewImage/N=RotatedImage RotatedOrig1
       ModifyImage RotatedOrig1 ctab= {zmin1,zmax1,Grays,0}
        //use bandpass to do preprocess   Yu Shi
        
       //use ROI
       duplicate/O/R=(xroimin, xroimax)(yroimin, yroimax)  rotatedorig1 rotatedROI
       bpass(postrad,3, rotatedROI)
       duplicate/o mfinal rotatedorig1
       
//       setscale/I x, xroimin, xroimax, rotatedorig1
//       setscale/I y, yroimin, yroimax, rotatedorig1
       
              // post array lattice constants
       dX = shiftwave[ifile][kiA]
       dY = shiftwave[ifile][kiB] //now defined in define array as sqrt(3)/2 lat b if hexcp ==1 or just kia if not
       

          for (i=0; i<Xposts; i+=1)
              for (j=0; j<Yposts; j+=1)
	       // determine whether or not to fit this post
	       if (maskwave[i][j] == KIgnore)
	           fitthispost = 0
	       else
	           fitthispost = 1
	      endif
	      ipost=i*yposts+j
	      if (fitthispost == 1)
                   if(ifile==f1)  //for the first frame	
                   	    if(mod(j,2)!=0)
                               xoffset = 0.5 *dx
                       	  else 
                       	  	xoffset = 0
                        endif
                            
                     xcenter = XTop + i*dX + xoffset  + j * dx*slopey - xroimin  // center of ROI, shifted by guess -- for odd rows, shifted over by 0.5 lattice spacing 
	               ycenter = YTop + j*dY + i  * dy* slopex - yroimin
	             else
	               xcenter=FRwave[ipost][0][ifile-1] - xroimin
	               ycenter=FRwave[ipost][1][ifile-1] - yroimin
	        endif


//trying to fit a possible bug	        
//	       xcmax = ceil( xcenter + centroidrad )
//	       xcmin = floor ( xcenter - centroidrad )
//	       ycmax =ceil( ycenter + centroidrad )
//	       ycmin =floor( ycenter - centroidrad )
//maximum number of loop now is 100

variable loopcounter = 0	       
	       
	      do
	      
	      loopcounter = loopcounter+1
	      
	       xcmax = ceil( xcenter + centroidrad )
	       xcmin = floor ( xcenter - centroidrad )
	       ycmax =ceil( ycenter + centroidrad )
	       ycmin =floor( ycenter - centroidrad )	      
	      
	      
	      centx=0
	      centy=0
	      centval=0
	      
	       for (ix=xcmin; ix<=xcmax; ix+=1)
	           for (iy=ycmin; iy<=ycmax; iy+=1)
	              xdis = ix - xcenter
	              ydis= iy - ycenter
	              if (sqrt( xdis^2+ydis^2) < centroidrad )
	                 centx += rotatedorig1[ix][iy] * ix
	                 centy += rotatedorig1[ix][iy] * iy
	                 centval += rotatedorig1[ix][iy]
	                 endif
	            endfor
	         endfor
	        
	       centx = centx / centval
	       centy = centy / centval
	       
	       epsx=centx-xcenter
	       epsy=centy-ycenter
	      
	      //if off center is larger than 0.6 pixel over 1 frame, redifine mask
               if (abs(epsx)>0.6)
                  xcenter=xcenter+1*sign(epsx)
               else
                  xcenter=centx
                  endif
                  
               if(abs(epsy)>0.6)
                   ycenter=ycenter+1*sign(epsy)
                   else
                   ycenter=centy
                   endif
                
                if( (abs(epsx)<=0.6) && (abs(epsy)<=0.6) ) 
                    FRwave[ipost][0][ifile]=centx+xroimin
                    FRwave[ipost][1][ifile]=centy +yroimin
                    break
                    
                    endif
                    
                if (loopcounter > 10000)
                    FRwave[ipost][0][ifile] = FRwave[ipost][0][ifile-1]
                    FRwave[ipost][1][ifile] = FRwave[ipost][1][ifile-1]
                    break
                  endif     
                   
                   
                   while (1)	       

	      

            endif
                
       endfor     
     endfor
    
    endfor
    
    
    
    
    end
                
	          	         
function COUNTcentroidgfit(f1,f2,XTop,YTop,XPosts,YPosts,zmin, zmax,maskwave,shiftwave,InitGuessWave)
    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    variable zmin,zmax      // min and max grayscale levels for display during fitting
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
    wave InitGuessWave   // Initial guesses for fit to first post.   These will be copied into w_coef1[]
                                        //  Should have same dimensions an definitions as w_coef1.  As of 2/10/06 this is length = 7
                                        //  Initial guesses for X0 and Y0 (post center coords) and post bottoms need not be supplied as initial guesses
                                        //  Initial guesses needed:
                                        //  InitGuessWave[0] =   Background
                                        //  InitGuessWave[1] =   Ampl
                                         //  InitGuessWave[2] =  Post Radius (in um).  This is a fixed value, presumably determined 
                                         // by a prior call to PSFfitForShifts followed by suitable averaging.
//    string pathstr            // path to data files

   variable i
   SVAR RootString, FitRes_String
   NVAR inputmovieID_G
   string fitaffix, fitname, moviename, movieaffix
   
   
   
   make/N=10/O/D freq={0.1,0.2,0.5,0.8,1.0,2.0,4.0,5.0,8.0,10.0,20.0,35.0,55.0,80.0,95.0,115.0,135.0}
   
   
   variable freqlen = dimsize(freq,0)
//   make/O/n=(XPosts*YPosts,2,f2-f1+1) fitresults // always created. calling routine must copy it as needed
//   SetScale/I z f1,f2,"", fitresults    




  // make/O/n=(XPosts*YPosts,k2,f2-f1+1) fiterrors // always created. two less columns - no chisq column or errorflag column needed
//   SetScale/I z f1,f2,"", fiterrors 
 //  fiterrors=0
 //  make/O/N=4 fitrangewave
     
          newpath videofolder    ///for cases when experiement and movies are stored at seperate places   
     
     for (i=0; i<freqlen;i+=1)
          if (freq[i] <1)
             sprintf fitaffix "_0_%dhz"  (freq[i]*10)
          else 
             sprintf fitaffix "_%dhz"  freq[i]
          ENDIF
          sprintf movieaffix "_%d.avi"  (i+1)
          moviename = Rootstring+movieaffix
          fitname =Fitres_string + fitaffix
          playmovie/p=videofolder as moviename 
          playmovieaction getID
          inputmovieID_G = V_Value
  // temporary change 11/13/2014         
           if( freq[i] < 1.0 )
             f2 = 18000
           elseif (freq[i]>=1 &&freq[i]<=10)
             f2 = 6000
           elseif (freq[i]>10)
             f2 = 3000  
           endif
           
           wave fitresults
           
           if (i>0)
           XTOP = fitresults[0][0][f2]
           YTOP = fitresults[0][1][f2]
           endif
           
           make/O/n=(XPosts*YPosts,2,f2-f1+1) fitresults // always created. calling routine must copy it as needed
           SetScale/I z f1,f2,"", fitresults    
           
           docentoridgfit(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,maskwave,shiftwave, fitresults)
           
           playmovieaction kill
           duplicate/O fitresults $fitname
           Save/O/P=topractice fitresults as fitname+".ibw"
           
    endfor
    
end   	          	         