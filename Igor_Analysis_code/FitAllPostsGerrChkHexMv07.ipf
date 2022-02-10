#pragma rtGlobals=1		// Use modern global access method.



// *************MODIFICATION HISTORY  - EFFECTIVE MAY 19, 2006
//12/12/21 DHR COUNTcentroidgfit2:  passed movie suffix; other updates to make magnetic analysis run in Igor9
//12/8/21 Change centroidgfit2() to allow suppression of images and diagnostic printout.
//12/6/21  DHR:  docentroidgfit()  modified for Igor9 playmovie action
//12/7/21 modified to open/close movie
//       added wave references 
//12/3/21  DHR removed doGfitAll() and GfitAllposts()
//11/3/14    added allowing analysis videos continuously

//6/16/14 use ROI for image processing

//6/5/14  change Gaussian fitting into centroid mask method

//12/18/13  use bandpass to do preprocessing

// 12/10/13 modified to read from movie

// 4/13/07  permits look at slopes of interp fits

// 4/5/07  COntains error checking for Version 6
   //    x0,y0 in ROI,    SigMin < abs(sigX), abs(sigY) < SigMax, abs(corr) < CorrMax, V_FitError !=0
   
//   3/24/07 NEW VERSION for Interp all frames

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

//--------bpass()

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

//--------- centroidgfit2()
//added on 6/16/2014, use centroid method to find posts center        
//12/8/21 Change  to allow suppression of images and diagnostic printout.
        
function centroidgfit2(f1,f2,XTop,YTop,XPosts,YPosts,zmin, zmax,maskwave,shiftwave,showimageflag, printeveryNframes)
    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    variable zmin,zmax      // min and max grayscale levels for display during fitting
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
    variable showimageflag   //12/7/21 =0 default;  = 1 to show image each frame (slower! - old approach)
    variable  printeveryNframes // diagnostic printout in command window every N frames to track progress; = -1 to suppress
        
   make/O/n=(XPosts*YPosts,2,f2-f1+1) fitresults // always created. calling routine must copy it as needed
   SetScale/I z f1,f2,"", fitresults    

    docentoridgfit2(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,maskwave,shiftwave,showimageflag, printeveryNframes, fitresults)

end   
       
//-------docentroidgfit2() -----------------------
//12/8/21 Change  to allow suppression of images and diagnostic printout.
//12/7/21 modified to open/close movie
// 12/6/21 modified for Igor9 movie actions
//  Fixed for Igor 9 - have to have wave assignments after call that creates them

       
function docentoridgfit2(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,maskwave, shiftwave,showimageflag, printeveryNframes, FRwave)      
    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    variable zmin,zmax      // min and max grayscale levels for display during fitting
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
    variable showimageflag   //12/7/21 =0 default;  = 1 to show image each frame (slower! - old approach)
    variable  printeveryNframes // diagnostic printout in command window every N frames to track progress; = -1 to suppress
    wave FRwave			// fit results wave.  Will be updated appropriately

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
   SVAR MoviePath    //12/7/21

   PlayMovieAction Open = MoviePath  //12/7/21
   PlayMovieAction frame = f1                // go to first frame to be analyzed
   for (ifile=f1; ifile <=f2; ifile +=1)   // loop over raw data files
       if ((mod(ifile,printeveryNframes)==0) && printeveryNframes >0)
           printf "Ifile = %d\n",ifile
       endif
       PlayMovieAction  extract
       if (ifile < f2)
          PlayMovieAction step = 1   //  move to next frame to be ready for next pass through loop
       endif
// new 12/10/13 for movie input
       WAVE M_MovieFrame
       ImageTransform rgb2gray M_MovieFrame   // convert to grayscale in M_rgb2gray
       WAVE M_rgb2gray       
       Duplicate/O M_rgb2gray Image1

       AngleOfRotation = shiftwave[ifile][kiAngle]
       Duplicate/O image1 OrigImage1
       ImageRotate/A=(AngleOfRotation) OrigImage1
       WAVE M_RotatedImage
       Duplicate/O M_RotatedImage RotatedOrig1
       
       DoWindow/K RotatedImage
       if (showimageflag == 1)  // 12/8/21
           NewImage/N=RotatedImage RotatedOrig1
           ModifyImage RotatedOrig1 ctab= {zmin1,zmax1,Grays,0}
       endif
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
    PlayMovieAction kill   //12/7/21 to close movie 
    end
                
//-------	    COUNTcentroidgfit2() --------
// does fits on multiple movie files for Magnetic rheology expts.
//12/12/21  passed moviesuffix.
//12/12/21   Changed to work with Igor9-compatible movie handling
//12/8/21 modified for optional showing each image file and printing frequencies to monitor progress
   	         
function COUNTcentroidgfit2(f1,f2,XTop,YTop,XPosts,YPosts,zmin, zmax,maskwave,shiftwave,showimageflag, printfrequencyflag,printeveryNframes,movietypestring)
    variable f1,f2  	// beginning and ending file numbers
    variable XTop,YTop//   position of upper left post
    variable XPosts,YPosts  // number of posts in x,y in ROI
    variable zmin,zmax      // min and max grayscale levels for display during fitting
    wave maskwave		// mask that IDs posts
    wave shiftwave           // shifts for each frame, angle etc
    variable showimageflag   // 12/8/21  = 1 to show image of each frame (old, slower way)
    variable printfrequencyflag  // = 1 to print out each frequency to monitor progress
    variable printEveryNframes   //12/8/21 flag print Iframe counter to monitor progress  (-1 to suppress
    string movietypestring // 12/12/21  for different movie types
    
   variable i
   SVAR RootString, FitRes_String
   NVAR inputmovieID_G
   SVAR MoviePath
   
   string fitaffix, fitname, moviename, movieaffix
   string MoviePathSave
   string videofolderpath, videofolderpathstring
   
   // list of frequencies measured
   make/N=10/O/D freq={0.1,0.2,0.5,0.8,1.0,2.0,4.0,5.0,8.0,10.0,20.0,35.0,55.0,80.0,95.0,115.0,135.0}
   
   
   variable freqlen = dimsize(freq,0)
     
     newpath/C/O/M="Open folder containing the movies." videofolderpath    ///for cases when experiment and movies are stored at separate places   
     pathinfo videofolderpath
     videofolderpathstring = S_path
          
     for (i=0; i<freqlen;i+=1)
          if (printfrequencyflag ==1)  // 12/8/21
              printf "frequency = %f Hz \n", freq[i]
          endif
          if (freq[i] <1)
             sprintf fitaffix "_0_%dhz"  (freq[i]*10)
          else 
             sprintf fitaffix "_%dhz"  freq[i]
          ENDIF
//12/8/21 Several things will need to change
    // movie name suffix must be flexible for PC and MAC
    // Set the current movie's path in the variable "MoviePath" as this is what docentroidgfit2() wants
    //get rid of code down to "temporary change"
    // delete "kill" below
          sprintf movieaffix "_%d.%s"  (i+1), movietypestring   //12/12/21
          moviename = Rootstring+movieaffix
          fitname =Fitres_string + fitaffix
          MoviePathSave = MoviePath // save to restore at end of routine
          MoviePath =  videofolderpathstring + moviename
          
  // code to provide number of frames in movies.  This changes with frequency        
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
           
           docentoridgfit2(f1,f2,XTop,YTop,XPosts,YPosts,zmin,zmax,maskwave,shiftwave, showimageflag, printeveryNframes, fitresults)
           
           MoviePath = MoviePathSave // 12/12/21 restore just in case           
           duplicate/O fitresults $fitname
           Save/O/P=topractice fitresults as fitname+".ibw"
           
    endfor
    
end   	          	         