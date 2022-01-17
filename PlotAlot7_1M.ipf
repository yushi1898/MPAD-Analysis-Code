#pragma rtGlobals=1		// Use modern global access method.
#pragma rtGlobals=1		// Use modern global access method.

// *************MODIFICATION HISTORY  - EFFECTIVE JUNE 2, 2006


//12/8/21  YS modified manual mask to gray for ignored posts
//12/8/21  YS modified MSD to be log scale, and to add units in both y and x axis
//12/7/21  DHR modified to make work in Igor 9:  added "DoUPdate" in a couple of key places.
// 12/5/21 YS modified to include plotting MSD items
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

static constant kFX0 = 2
static constant kFY0 = 4
static constant KFMSD = 1


function PlotALot4EB(posti,postf,param1,param2,p1string,p2string,astrw,maskwave,datawave,bgshiftwave,cellpostflag,notebookflag,bgshiftflag, framerate, spikeflag)
														   
    variable posti,postf              //posti==first post (0),postf==final post(55)
    variable param1,param2    // number of parameter to be plotted  
    string p1string,p2string    // names of parameters for y-axes
    wave astrw                            // astr
    wave maskwave
    wave datawave                      // fitresults to be plotted
//    wave errorwave			   //  error bars (from Curvefit?)
    wave  bgshiftwave            //added by Yu Shi 1/3/14
    variable cellpostflag                // == 1 to plot only cellposts
   variable notebookflag  //Added by DHR 3/9/06
   variable bgshiftflag      //added by YuShi   1/3/14
     variable framerate
     variable spikeflag
    SVAR RootString
    
    
    string xstr,ystr,layoutstr
    variable ipost, ipcount, i,j,nxpost,nypost
    variable plotsperpage = 7
    variable plotsoncurrentpage = 0
	 variable posttype
	
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
    
	    if (maskwave[ipost] == 6)
	    	posttype= 6
	    elseif(maskwave[ipost] ==7)
	    	posttype = 7
	    elseif(maskwave[ipost]==-1 )
	    	posttype=-1
	    elseif(maskwave[ipost] ==-2 )
	    	posttype=-2
	    elseif(maskwave[ipost] ==1)
	    	posttype = 1
	    elseif(maskwave[ipost] ==3)
	    	posttype = 3
	    elseif(maskwave[ipost] ==4)
	    	posttype = 3 
	    elseif(maskwave[ipost] == 0)
	    	posttype = 0   
	    endif    
    
    
    
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
 	       if (posttype == 1 || posttype == 6 || posttype ==7)	   
      	     dispfitres4EB(datawave,ipost,param1, p1string,xstr,bgshiftwave,plotsoncurrentpage,posttype, i,j,cellpostflag,bgshiftflag, framerate, spikeflag)
        	     dispfitres4EB(datawave,ipost,param2,p2string,ystr,bgshiftwave,plotsoncurrentpage+plotsperpage,posttype,i,j,cellpostflag,bgshiftflag, framerate, spikeflag)
    
   	        AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $xstr
      	     AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $ystr
              ipcount = ipcount+1
              plotsoncurrentpage =  plotsoncurrentpage + 1
              if ((plotsoncurrentpage == 7) || (ipost == postf))
                   DoUpdate  //12/7/21  added to make work in Igor9.  
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
          
          dispfitres4EB(datawave,ipost,param1, p1string,xstr,bgshiftwave,plotsoncurrentpage,posttype, i,j,cellpostflag,bgshiftflag,framerate, spikeflag)
        	 dispfitres4EB(datawave,ipost,param2,p2string,ystr,bgshiftwave,plotsoncurrentpage+plotsperpage,posttype,i,j,cellpostflag,bgshiftflag,framerate, spikeflag)
    
   	    AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $xstr
      	 AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $ystr
          ipcount = ipcount+1
          plotsoncurrentpage =  plotsoncurrentpage + 1
          if ((plotsoncurrentpage == 7) || (ipost == postf))
              DoUpdate  //12/7/21  added to make work in Igor9.  
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



function dispfitres4EB(frwave,npost,nparam,pnamestr,swname,bgshiftwave,nplot,posttype,i,j,cellpostflag,bgshiftflag,framerate, spikeflag)

    wave frwave          // fitresults wave
    wave  bgshiftwave          // background shift wave 
    variable npost       // number of post to be plotted
    variable nparam    // number id of fit parameter to be plotted
    string pnamestr    // y axis label
    string swname       // string of Plot (graph) name
    variable nplot       // number of plot
    variable posttype // type of post for coloring graph
    variable i,j           // matrix indices of post 
    variable cellpostflag
    variable bgshiftflag    // choose whether to use bgshift wave 
    variable framerate
    variable spikeflag
    
    string s1, s2,dws,ews
    wave timewave


wave zeropos
    
wave MSD_ed, MSDxwave
    
    sprintf dws, "dw%d",nplot
    sprintf ews, "ew%d",nplot
    variable wlen,woff
    variable psym
    variable rval = 0, gval = 0, bval = 0
    
    variable iframe
    
    variable pixratio=0.125   // pixel to micron ratio
    variable dstart = 0
    
    woff = DimOffset(frwave,2)
    wlen = DimSize(frwave,2)
    make/O/N=(wlen) dwtemp,ewtemp //$dws,$ews, 
    SetScale/P x woff,1,"", dwtemp,ewtemp//$dws,$ews
    
    if (nparam==2)
    
    dwtemp = frwave[npost][0][p]                
    
    elseif (nparam == 4)
    
    dwtemp = frwave[npost][1][p]
    
    elseif (nparam == 1) 
    dwtemp = MSD_ed[p][npost]  
    endif
   //TEMPORARY! TO USE PLOT ALOT TO PLOT INDIV POST ENERGY
    //dwtemp*=dwtemp*.5*26.4*.109*.109 //convert displacement to energy 26.4 is spring constant , .109 pix to micron conversion


//changed it to deflection zeros at resting position of posts    


if(bgshiftflag!=0)    
    if ( nparam== 2)
        for ( iframe=woff; iframe<(woff+wlen); iframe+=1) 
              dwtemp[iframe]-=bgshiftwave[iframe][0]
              
         endfor
         
     		if (bgshiftflag==2)
         
       	 	 dstart = zeropos[npost][0][0]
      	   dwtemp=dwtemp-dstart
      	   dwtemp=dwtemp*pixratio
         endif
     
     
    elseif (nparam==4)
      	for ( iframe=woff; iframe<(woff+wlen); iframe+=1) 
     	          dwtemp[iframe]-=bgshiftwave[iframe][1]
               
     	 	endfor
      	if (bgshiftflag==2)
      		dstart = zeropos[npost][1][0]
      		dwtemp=dwtemp-dstart
      		dwtemp=dwtemp*pixratio
     		endif
     endif
      
 endif
   
    
    
    
    

    
    if (spikeflag==1)
       wavestats/Q dwtemp
       
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

    
    // set up for different symbols
    if (posttype == kGuide)
        psym = 19
    elseif (posttype == 1)
        psym = 8
        gval = 50000
    elseif (posttype == 0)
        psym = 6
        bval = 0
        rval = 0
        gval = 0
    elseif (posttype == 3)
        psym = 16
        rval = 39321
        bval =  39321       
        gval = 39321
    elseif (posttype == 4)
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

	 if (nparam ==1)
	 	display/K=1/N = $swname $dws vs MSDxwave 
	 	ModifyGraph log=1
	 	Label bottom "Lag time (s)"
	 	Label left pnamestr+" (nm\S2\M)"
	 else
    	display/K=1/N = $swname $dws 
    	Label bottom "Time (s)"
    	if(bgshiftflag!=2) 
    		Label left pnamestr+" (px)"
    	else
    		Label left pnamestr+" (\\F'GreekC'm\\F'Arial'm)"
    	endif
    endif
    

    //set the scale to be min and max of data for tiral   12/11/13  Yu
    SetAxis left, wavemin($dws),wavemax($dws) // 4/23/10 forcing scale in plotalots to min,max of data not errorbars


    sprintf s1, "Post %d (%d,%d)\rType %d",npost,i,j,posttype
    TextBox/C/B=1/N=text0/A=LT s1
//    Label left pnamestr
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


