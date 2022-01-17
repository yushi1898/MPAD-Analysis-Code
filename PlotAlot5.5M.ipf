#pragma rtGlobals=1		// Use modern global access method.
#pragma rtGlobals=1		// Use modern global access method.

// *************MODIFICATION HISTORY  - EFFECTIVE JUNE 2, 2006

//5/7/14  remove error bar, add timescale and pixel micron ratio

// 1/3/14   added bgshift option to substract bgshift from deflection    Yu Shi
// 4/23/10 Version 5.3   forcing scale in plotalots to min,max of data not errorbars
//7/13/06:  adds (i,j) to text box; automatic save has been added
// 6/2/06 includes errorbars - working on it...!!!
// 3/3/06 and printing option




// maskwave constants

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
static constant kFX0 = 2
static constant kFSigX= 3
static constant kFY0 = 4
static constant kFSigY = 5
static constant kFCorr = 6
static constant kFChisqG = 7
static constant kFerrorFlag = 8


function PlotALot4EB(posti,postf,param1,param2,p1string,p2string,astrw,maskwave,datawave,bgshiftwave,cellpostflag,notebookflag,timeflag,bgshiftflag, framerate, spikeflag)
														   
    variable posti,postf              //posti==first post (0),postf==final post(55)
    variable param1,param2    // number of parameter to be plotted  param1==coef,param2==coef
    string p1string,p2string    // names of parameters for y-axes
    wave astrw                            // astr
    wave maskwave
    wave datawave                      // fitresults to be plotted
//    wave errorwave			   //  error bars (from Curvefit?)
    wave  bgshiftwave            //added by Yu Shi 1/3/14
    variable cellpostflag                // == 1 to plot only cellposts
   variable notebookflag  //Added by DHR 3/9/06
   variable bgshiftflag      //added by YuShi   1/3/14
     variable timeflag        //  = 1 for vs time  plots
     variable framerate
     variable spikeflag
    SVAR RootString
    
    
    wave linmaskval_begin, linmaskval_end
    wave linmaskval_begin_auto, linmaskval_end_auto
    string xstr,ystr,layoutstr
    variable ipost, ipcount, i,j,nxpost,nypost
    variable plotsperpage = 7
    variable plotsoncurrentpage = 0

    wave errorwave
	
    string nbname
    nxpost = astrw[3]
    nypost = astrw[4]
    
if (notebookflag == 1)
     //sprintf NBname, "NBPlotALot_"+RootString+"_"+"%d%d", param1,param2
     sprintf NBname, "NBPlotALot"
    DoWindow/K $NBname
    NewNotebook/K=1/F=1/N=$NBname
endif

    ipcount =0

    for (ipost = posti; ipost <= postf; ipost = ipost + 1)
       // get indices for maskwave
         i = floor(ipost/nypost)
         j = mod(ipost,nypost)
       
        if (mod(ipcount, plotsperpage) == 0)
            plotsoncurrentpage = 0
            DoWindow/K LayoutXY
            NewLayout/K=1/N=LayoutXY
        endif
        sprintf xstr,"Plot%dx",plotsoncurrentpage
        sprintf ystr,"Plot%dy",plotsoncurrentpage
    
 	if(cellpostflag==1)
 		
 	//	if(maskwave[i][j]==Kcell)
 	   if (linmaskval_begin_auto[ipost] ==1 || linmaskval_end_auto[ipost] ==1)	   
      		  dispfitres4EB(datawave,ipost,param1, p1string,xstr,errorwave,bgshiftwave,plotsoncurrentpage,maskwave[i][j], i,j,timeflag,cellpostflag,bgshiftflag, framerate, spikeflag)
        	  dispfitres4EB(datawave,ipost,param2,p2string,ystr,errorwave,bgshiftwave,plotsoncurrentpage+plotsperpage,maskwave[i][j],i,j,timeflag,cellpostflag,bgshiftflag, framerate, spikeflag)
    
   	       AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $xstr
      	       AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $ystr
              ipcount = ipcount+1
             plotsoncurrentpage =  plotsoncurrentpage + 1
                if ((plotsoncurrentpage == 7) || (ipost == postf))
                  DoWindow/F LayoutXY
             //if (printflag == 1)
                   //PrintLayout LayoutXY  //Modified Jul 27
                    if (notebookflag == 1)
                       Notebook $NBname scaling = {90,90}, picture = {LayoutXY, -1,1}
                       Notebook $NBname text = "\r"
                    else
                       SavePICT/E=-6/B=288
                    endif
                 endif
               endif
         else
          
               dispfitres4EB(datawave,ipost,param1, p1string,xstr,errorwave,bgshiftwave,plotsoncurrentpage,maskwave[i][j], i,j,timeflag,cellpostflag,bgshiftflag,framerate, spikeflag)
        	  dispfitres4EB(datawave,ipost,param2,p2string,ystr,errorwave,bgshiftwave,plotsoncurrentpage+plotsperpage,maskwave[i][j],i,j,timeflag,cellpostflag,bgshiftflag,framerate, spikeflag)
    
   	       AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $xstr
      	       AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $ystr
              ipcount = ipcount+1
             plotsoncurrentpage =  plotsoncurrentpage + 1
                if ((plotsoncurrentpage == 7) || (ipost == postf))
                  DoWindow/F LayoutXY
             //if (printflag == 1)
                   //PrintLayout LayoutXY  //Modified Jul 27
                    if (notebookflag == 1)
                       Notebook $NBname scaling = {90,90}, picture = {LayoutXY, -1,1}
                       Notebook $NBname text = "\r"
                    else
                       SavePICT/E=-6/B=288
                    endif
 	          endif
              endif
         endfor
    if (notebookflag==1)
	   SaveNotebook/P=ToPractice $NBname        //added by sasha, July 13  '06
	   //if (printFlag==1)
	  	//	PrintNotebook $NBname  //Modified Jul 27
	  // endif
     endif
end



function dispfitres4EB(frwave,npost,nparam,pnamestr,swname,errorwave,bgshiftwave,nplot,posttype,i,j,timeflag,cellpostflag,bgshiftflag,framerate, spikeflag)

    wave frwave          // fitresults wave
    wave  bgshiftwave          // background shift wave 
    variable npost       // number of post to be plotted
    variable nparam    // number id of fit parameter to be plotted
    string pnamestr    // y axis label
    string swname       // string of Plot (graph) name
    wave errorwave    // wave of errorbars
    variable nplot       // number of plot
    variable posttype // type of post for coloring graph
    variable i,j           // matrix indices of post 
    variable timeflag  // =1 for plots vs time
    variable cellpostflag
    variable bgshiftflag    // choose whether to use bgshift wave 
    variable framerate
    variable spikeflag
    
    string s1, s2,dws,ews
    wave timewave
    
//    wave linmaskval_begin, linmaskval_end
//    wave linmaskval_begin_auto, linmaskval_end_auto

wave linmaskval_intercept_int
wave linmaskval_intercept_int_sc
wave zeropos
    
    if (linmaskval_intercept_int[npost] == 6)
    posttype= 6
    elseif(linmaskval_intercept_int[npost] ==7)
    posttype = 7
    elseif(linmaskval_intercept_int[npost]==-1 )
    posttype=-1
    elseif(linmaskval_intercept_int[npost] ==-2 )
    posttype=-2
    elseif(linmaskval_intercept_int[npost] ==1)
    posttype = 1
    elseif(linmaskval_intercept_int[npost] ==3)
    posttype = 3
    elseif(linmaskval_intercept_int[npost] ==4)
    posttype = 3    
    endif
    
    sprintf dws, "dw%d",nplot
    sprintf ews, "ew%d",nplot
    variable wlen,woff
    variable psym
    variable rval = 0, gval = 0, bval = 0
    
    variable iframe
    
    variable pixratio=0.125   // pixel to micron ratio
    variable dstart
    
    woff = DimOffset(frwave,2)
    wlen = DimSize(frwave,2)
    make/O/N=(wlen) dwtemp,ewtemp //$dws,$ews, 
    SetScale/P x woff,1,"", dwtemp,ewtemp//$dws,$ews
    
    if (nparam==2)
    
    dwtemp = frwave[npost][0][p]                
    ewtemp = errorwave[npost][0][p]
    
    else
    
    dwtemp = frwave[npost][1][p]
    ewtemp = errorwave[npost][1][p]    
    endif
   //TEMPORARY! TO USE PLOT ALOT TO PLOT INDIV POST ENERGY
    //dwtemp*=dwtemp*.5*26.4*.109*.109 //convert displacement to energy 26.4 is spring constant , .109 pix to micron conversion


//changed it to deflection zeros at resting position of posts    


if(bgshiftflag==1)    
    if ( nparam== 2)
        for ( iframe=woff; iframe<(woff+wlen); iframe+=1) 
              dwtemp[iframe]-=bgshiftwave[iframe][0]
              dstart = zeropos[npost][0][0]
         endfor
    elseif (nparam==4)
      for ( iframe=woff; iframe<(woff+wlen); iframe+=1) 
               dwtemp[iframe]-=bgshiftwave[iframe][1]
               dstart = zeropos[npost][1][0]
               endfor
      endif
      
 endif

//dstart=dwtemp[0]   
    dwtemp=dwtemp-dstart

          
    
    
    
    dwtemp=dwtemp*pixratio
    

    
    if (spikeflag==1)
       wavestats/Q dwtemp
       
//          for (iframe=woff; iframe<(wlen+woff); iframe+=1)                 
//                 if ( abs(dwtemp[iframe] -V_avg) > 3 * V_sdev )
//                     if (iframe==woff)
//                     dwtemp[iframe] = 2* dwtemp[iframe+1] - dwtemp[iframe+2]
//                      elseif (iframe==(wlen+woff-1))
//                        dwtemp[iframe] = 2* dwtemp[iframe-1] - dwtemp[iframe-2]
//                      else
//                       dwtemp[iframe] = (dwtemp[iframe-1]+dwtemp[iframe+1] ) /2
//                 endif
//                 endif
                 
//                 endfor
             for (iframe=woff+1; iframe<(wlen+woff); iframe+=1)  
                  if ( abs(dwtemp[iframe]-dwtemp[iframe-1]) > 2 * V_sdev ) 
//                     if (iframe==woff+1)
                       dwtemp[iframe]=dwtemp[iframe-1]  
  //                  else  
  //                    dwtemp[iframe] = 2* dwtemp[iframe-1] - dwtemp[iframe-2]
  //                   endif
                   endif
                   
                endfor 
                  
                 
   endif
   
   setscale/P x 0, (1/framerate),"",dwtemp
     
    duplicate/O dwtemp $dws

   
    
  //  ewtemp*=.109 //convert pixels to microns
  //  ewtemp*=(ewtemp)
    
    
    duplicate/O ewtemp $ews
    
    // set up for different symbols
    if (posttype == kGuide)
        psym = 19
    elseif (posttype == kCell)
        psym = 8
        gval = 50000
    elseif (posttype == kEmpty)
        psym = 6
//         bval = 50000
        rval = 65280
        gval = 65280
    elseif (posttype == kIgnore)
        psym = 16
//        gval = 50000
    elseif (posttype == kWirepost)
        psym = 21
        rval = 30000
        bval =  30000
     elseif (posttype == 6)
       bval =50000
    elseif (posttype == -1)
    rval = 65280
    gval = 43520
    elseif ( posttype ==-2)    
    rval = 50000
//    gval = 25000
     elseif (posttype == 7)
     rval = 36864
     gval = 14592
     bval = 58880
    endif
    DoWindow/K $swname
    if (timeflag == 1)
        display/N=$swname $dws vs timewave
    else
    	 
    	 //putting this in temporarily to calculate energy
    	
//    	 variable kspring = 26.4
//    	 duplicate $dws energywave
//    	 energywave*=energywave*0.5*kspring*(.109)^2
//    	 
//    	 display/K=1/N = $swname energywave
    	 
    	 
        display/K=1/N = $swname $dws
       
    endif
    //set the scale to be min and max of data for tiral   12/11/13  Yu
    SetAxis left, wavemin($dws),wavemax($dws) // 4/23/10 forcing scale in plotalots to min,max of data not errorbars
  // SetAxis left, wavemin($ews),wavemax($ews)
//   ErrorBars $dws Y,wave=($ews,$ews)       disabled error bar   5/7/2014

    sprintf s1, "Post %d (%d,%d)\rType %d",npost,i,j,posttype
    TextBox/C/B=1/N=text0/A=LT s1
    Label left pnamestr
    ModifyGraph tick=2,mirror=1,standoff=0
    ModifyGraph mode=0,rgb=(rval,gval,bval)    //changed plotstyle 5/7/2014
    
    ModifyGraph margin(left)=50,margin(bottom)=23,margin(top)=2,margin(right)=5;
    ModifyGraph width=165,height=75
   
end

function dispfitres2(frwave,npost,nparam,pnamestr,swname)

    wave frwave          // fitresults wave
    variable npost       // number of post to be plotted
    variable nparam    // number id of fit parameter to be plotted
    string pnamestr    // y axis label
    string swname       // string of Plot (graph) name

    string s1, s2

    DoWindow/K $swname
    display/K=1/N = $swname frwave[npost][nparam][] 
    sprintf s1, "Post %d",npost
    TextBox/C/B=1/N=text0/A=LT s1
    Label left pnamestr
    ModifyGraph tick=2,mirror=1,standoff=0
    ModifyGraph mode=4,marker=19,rgb=(0,0,0)

    ModifyGraph margin(left)=50,margin(bottom)=23,margin(top)=2,margin(right)=5;
    ModifyGraph width=165,height=75
end

function comparenn()  //comparison of nearest neighbors, to determine correlation of response

wave frwave          // fitresults wave
    variable npost       // number of post to be plotted
    variable nparam    // number id of fit parameter to be plotted
    string pnamestr    // y axis label
    string swname       // string of Plot (graph) name


end


