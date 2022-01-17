#pragma rtGlobals=1		// Use modern global access method.// maskwave constantsstatic constant kEmpty = 0static constant kCell = 1static constant kGuide = 2static constant kIgnore = 3static constant kWirepost = 4static constant kFitAll = -1// indices of Shiftwave// note that primary index of Shiftwave must correspond to file numberstatic constant kiAngle = 7static constant kiA = 8static constant kiB = 9static constant kiXsi  =10static constant kiYsi = 11static constant kiXsf = 12static constant kiYsf =13// indices for fitparams for Gauss + quadratic fitstatic constant kNfitParamsBase  = 10  // includes column for chisq +1 for gauss2dquadstatic constant kFbkgnd = 0static constant kFAmpl = 1static constant kFX0 = 2static constant kFSigX= 3static constant kFY0 = 4static constant kFSigY = 5static constant kFCorr = 6static constant kFB = 9       // ampl. of quadratic bitstatic constant kFChisq = 8//static constant kfitrange  = 20   //+- range for fitting in pixels  //*****************************************************//   doGr2FitBase// 6/7/06 ccd error option added//   2/6/06    modified version of  doGfit routine.   // Uses modified 2D gaussian to fit a single BASE image of a 60x oil data set.//  Goal:  obtain post center positions for use in PSF fits // OUTPUT://   wave fitresultsbase      -  a 2D fit-results style file.  Note that there is an extra parameter compared//                                                 to the usual 2D gaussian -  the "B*r^2" term in the function //   wave fiterrorsbase      - Igor's estimate of errors of the fit parameters;  // Post bottom image is assumed to have name such as RHODBASE// Fits all posts except those flagged in mask as Ignorefunction doGr2FitBase(XTop,YTop,XPosts,YPosts,AngleOfRotation,Abase,Bbase,zmin,zmax,maskwave,InitGuessWave,pathstr)    variable XTop,YTop         //   position of upper left post    variable XPosts,YPosts   // number of posts in x,y in ROI    variable AngleOfRotation // Rotation angle for base image    variable Abase,Bbase       // lattice constants of post array    variable zmin,zmax      // min and max grayscale levels for display during fitting    wave maskwave		// mask that IDs posts    wave InitGuessWave   // Initial guesses for fit to first post.   These will be copied into W_coef[]                                        //  Should have same dimensions an definitions as W_coef.  As of 2/10/06 this is length = 8                                        //  Initial guesses for X0 and Y0 (post center coords) need not be supplied                                        //  Initial guesses needed:                                        //  InitGuessWave[0] =   Background                                        //  InitGuessWave[1] =   Ampl                                        //  InitGuessWave[3] =   Sig X   (in range 6-10?)                                        //  InitGuessWave[5] =   Sig Y                                        //  InitGuessWave[6] =   Corr   (0.01 is ok)                                        //  InitGuessWave[7] =   B        (approx 0.5 *Ampl?)                                            string pathstr            // path to data files       string s1   variable i, j, ifile, fitthispost   variable xcenter, ycenter, Xsp, Ysp   variable V_chisq   variable dX,dY   // local variables for lattice constants of post array   nvar hexcp     	Variable/G  kfitrange=floor(Abase*0.4)   NVAR CCDflag, CCDreadnoise, CCDeperADU   // 6/7/06    // make output waves     make/O/n=(XPosts*YPosts,kNFitParamsBase) fitresultsbase // always created. calling routine must copy it as needed   fitresultsbase = 0   make/O/n=(XPosts*YPosts,kNFitParamsBase-1) fiterrorsbase // always created. one less column - no chisq column needed   fiterrorsbase = 0      // sprintf s1, "%s_BASE.tif",pathstr   variable zscale = 1 // 6/7/06 for ccd     s1="Rhod000.tif"   //ImageLoad/P=toPractice/T=tiff/O/N=image1 s1    // load, copy, and rotate the image   ImageLoad/P=MovieFolder/T=tiff/O/N=image1 s1    // load, copy, and rotate the image   Duplicate/O image1 OrigImage1   if (CCDflag ==1) // use proper statistical errors from CCD  6/7/06        origImage1 = OrigImage1 - CCDBias(x,y)   endif   ImageRotate/A=(AngleOfRotation) OrigImage1   Duplicate/O M_RotatedImage RotatedBase   if (CCDflag == 1)       RotatedBase *= CCDeperADU  // convert to e-       Duplicate/O RotatedBase RotatedBaseErr       RotatedBaseErr = sqrt(RotatedBaseErr + CCDreadnoise^2)         zscale = CCDeperADU // conversion factor to e-       zmin = CCDeperADU*(zmin - CCDavgBias())       zmax = CCDeperADU*(zmax - CCDavgBias())   endif   DoWindow/K RotatedImage   NewImage/K=1/N=RotatedImage RotatedBase   ModifyImage RotatedBase ctab= {zmin,zmax,Grays,0}  // local variables for lattice constants of post array   dX = Abase   dY = Bbase          // Draw a box around ROI   SetDrawLayer ProgFront   SetDrawEnv linefgc= (65535,65535,0),fillpat= 0,xcoord= top,ycoord= left, save    DrawRect XTop  - dX/2, YTop  - dY/2 ,XTop  + (XPosts- 0.5)*dX, YTop  +(YPosts- 0.5)*dY   SetAxis left  YTop  + YPosts*dY, YTop - dY   SetAxis top XTop - dX, XTop  + XPosts*dX       // make w_coef and insert initial guesses   // zscale converts to e- counts if needed   make/O/n=(kNfitParamsBase -1) W_coef   W_coef[kFbkgnd]  = InitGuessWave[kFbkgnd]  * zscale // bggnd   W_coef[kFAmpl]  = InitGuessWave[kFAmpl]   * zscale // A   W_coef[kFSigX]  = InitGuessWave[kFSigX]      // SigmaX   W_coef[kFSigY]  = InitGuessWave[kFSigY]     // SigmaY   W_coef[kFCorr] = InitGuessWave[kFCorr]      // corr   W_coef[kFB]   = InitGuessWave[kFB]     * zscale       // B          // Fit the modified gaussians       for (i=0;i<=(Xposts-1); i+=1)	for (j=0;j<=(YPosts-1); j+=1)	       // determine whether or not to fit this post	   if (maskwave[i][j] == KIgnore)              fitthispost = 0	   else    // fit all posts except ignored ones	       fitthispost = 1	   endif	   if (fitthispost == 1)	       xcenter = XTop + i*dX 	       ycenter = YTop + j*dY 	      // initial guesses  for post center  (rest given by initial guess, or from fit to previous post)	     W_coef[kFX0] = xcenter	     W_coef[kFY0] = ycenter	     if (CCDflag == 1)   // 6/7/06 fit with error bars	         FuncFitMD/Q Gauss2Dquad  W_coef RotatedBase[(xcenter - kfitrange),(xcenter + kfitrange)][(ycenter - kfitrange),(ycenter + kfitrange)] /D /I=1 /W=RotatedBaseErr ///C=T_Constraints 	     else                        // fit with no weighting	         FuncFitMD/Q Gauss2Dquad  W_coef RotatedBase[(xcenter - kfitrange),(xcenter + kfitrange)][(ycenter - kfitrange),(ycenter + kfitrange)] /D  ///C=T_Constraints 	     endif	     wave w_sigma	     // convert fit params back from e- to ADU as needed	     fitresultsbase[i* YPosts + j][kFbkgnd] = W_coef[kFbkgnd]  / zscale  // bggnd	     fitresultsbase[i* YPosts + j][kFAmpl]  = W_coef[kFAmpl] / zscale  // A	     fitresultsbase[i* YPosts + j][kFX0]      = W_coef[kFX0]   // x0	     fitresultsbase[i* YPosts + j][kFSigX]   = W_coef[kFSigX]   // SigmaX	     fitresultsbase[i* YPosts + j][kFY0]      = W_coef[kFY0]   // y0	     fitresultsbase[i* YPosts + j][kFSigY]   = W_coef[kFSigY]   // SigmaY	     fitresultsbase[i* YPosts + j][kFCorr]  = W_coef[kFCorr]   // corr	     fitresultsbase[i* YPosts + j][kFB]       = W_coef[kFB]  / zscale // B	     fitresultsbase[i* YPosts + j][kFChisq] = V_Chisq   // Chisq	            //now record error estimates	     fiterrorsbase[i* YPosts + j][kFbkgnd] = W_sigma[kFbkgnd]  / zscale // bggnd	     fiterrorsbase[i* YPosts + j][kFAmpl]  = W_sigma[kFAmpl]/ zscale  // A	     fiterrorsbase[i* YPosts + j][kFX0]      = W_sigma[kFX0]   // x0	     fiterrorsbase[i* YPosts + j][kFSigX]   = W_sigma[kFSigX]   // SigmaX	     fiterrorsbase[i* YPosts + j][kFY0]      =  W_sigma[kFY0]   // y0	     fiterrorsbase[i* YPosts + j][kFSigY]   = W_sigma[kFSigY]   // SigmaY	     fiterrorsbase[i* YPosts + j][kFCorr]  = W_sigma[kFCorr]   // corr	     fiterrorsbase[i* YPosts + j][kFB]       = W_sigma[kFB] / zscale // corr	  endif       endfor //j    endFor //iend//************************************************// gaus2dquad()// Gaussian plus quadratic to model bottom of posts     function gauss2dquad(w,x,y) : fitfunc  wave w variable x,y variable bkgnd,A,X0,sx,Y0,sy,cor,B,val,r2  bkgnd = w[kFbkgnd] A = w[kFAmpl] X0 = w[kFX0] sx = w[kFSigX] Y0 = w[kFY0] sy = w[kFSigY] cor = w[kFCorr] B = w[kFB]  r2 = ((x-X0)/sx)^2 + ((y-Y0)/sy)^2  val = bkgnd + (A+B*r2)*exp( (-1/(2*(1-cor^2)))*( ((x-X0)/sx)^2 + ((y-Y0)/sy)^2 - 2*cor*(x-X0)*(y-Y0)/(sx*sy) ) ) return valendfunction CCDbias(x,y)    variable x,y        return 195.5    endfunction CCDavgBias()    return 195.5end