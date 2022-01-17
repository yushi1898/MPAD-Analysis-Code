 #pragma rtGlobals=1		// Use modern global access method and strict wave access.
// *************MODIFICATION HISTORY  - EFFECTIVE Dec 12, 2013

 
// 9/7/21 compilation bug in function makeNNunit()  fixed by DHR 

//velocity correlation function different    


//9/20/17   modified erf fitting regim, depend on sigma value
//3/1/2016      calculate the frequency dependence for MSD curve

//10/26/2015   enable MSD calculation based on local vector coordinate



//10/20/2015    correlation function based on local displacement vector    


//6/11/2015   normalized correlation function
//4/30/2015   added a new way of calculating correlation function(trial)
//12/3/2014 change the MDS calculation procedure to improve speed
//11/4/2014  added a procedure to calculate background shift under all frequency
//add function to remove spikes
//calcbacgroundshift and calcMSD function in this procedure
//change the bgshift calculation, fixed the bug in plot     1/3/14
//added function to calculate MSDstatistics
//added timescale adjustment

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

static constant kfx0 = 0
static constant kFy0 = 1




//12/11/13    to calculate background shift    Yu Shi
function calcbackgroundshift(astrw,maskwave,datawave,bgshift, iframe1, iframe2)
wave maskwave
wave datawave
wave astrw
wave bgshift
variable iframe1,iframe2


Variable npost
variable ipost
variable i,j,k
variable nxpost, nypost
variable iframe
variable nbgpost=0

variable navgf = 10
variable nframe_s = floor((iframe2-iframe1)/navgf)






//set row and colum number of post
nxpost=astrw[3]
nypost=astrw[4]
npost=nxpost*nypost


for(i=0; i<nxpost; i+=1)
    for(j=0; j<nypost; j+=1) 
        ipost=i*nypost+j
         if(maskwave[ipost]==kEmpty)
         nbgpost+=1
         endif
      endfor
 endfor


for(iframe=iframe1; iframe<iframe2; iframe+=1)
    bgshift[iframe][0]=0
    bgshift[iframe][1]=0
     for(i=0; i<nxpost; i+=1)
          for(j=0; j<nypost; j+=1) 
                ipost=i*nypost+j
                if(iframe>iframe1)         //set the shift of the first frame 0
                   if(maskwave[ipost]==kEmpty)
                       bgshift[iframe][0]+=datawave[ipost][kFx0][iframe]-datawave[ipost][kFx0][iframe1]
                       bgshift[iframe][1]+=datawave[ipost][kfy0][iframe]-datawave[ipost][kfy0][iframe1]
                   endif
                 endif
            endfor
       endfor
       bgshift[iframe][0]/=nbgpost
       bgshift[iframe][1]/=nbgpost
 endfor
 
 make/O/N=(npost,2,(iframe2-iframe1)) fitres_ed
 make/o/n=(npost,2,nframe_s) fitres_ed_s
 make/o/n=(navgf)   xtemp_avgseg, ytemp_avgseg
 
 for (iframe = iframe1; iframe<iframe2; iframe+=1)
      for(i=0; i<nxpost; i+=1)
            for(j=0; j<nypost; j+=1) 
                    ipost=i*nypost+j
                    if(maskwave[ipost]!=kignore) 
                      fitres_ed[ipost][kfx0][iframe] = datawave[ipost][kfx0][iframe] - bgshift[iframe][0]
                      fitres_ed[ipost][kfy0][iframe] = datawave[ipost][kfy0][iframe] - bgshift[iframe][1]
                    endif
             endfor
        endfor
 endfor
 
 
 
  for(ipost = 0; ipost<npost; ipost+=1)       
  if (maskwave[ipost] !=kignore)
        for(iframe=0;iframe<iframe2; iframe+=1)
              j = mod(iframe,navgf)
              k=floor(iframe/navgf)
              xtemp_avgseg[j] = datawave[ipost][0][iframe] - bgshift[iframe][0]
              ytemp_avgseg[j] = datawave[ipost][1][iframe] - bgshift[iframe][1]
              if(j == (navgf-1) )
                 wavestats/q xtemp_avgseg         
                  fitres_ed_s[ipost][0][k] = V_avg
                wavestats/q ytemp_avgseg  
                  fitres_ed_s[ipost][1][k] = V_avg
              endif
       endfor   
  endif     
 endfor       
 
 
end
                  
function PlotMSDM(posti,postf,astrw, maskwave,datawave,cellpostflag,notebookflag,timeflag,maxt, slopet, framerate, spikeflag, directionflag)
  variable posti,postf              //posti==first post (0),postf==final post(55)
  wave astrw                            // astr
  wave maskwave
   wave datawave                      // fitresults to be plotted  
   variable cellpostflag                // == 1 to plot only cellposts
   variable notebookflag  //Added by DHR 3/9/06
   variable timeflag        //  = 1 for vs time  plots
   variable maxt
   variable slopet
   variable framerate
   variable spikeflag
   variable directionflag     // =1 for ploting in local vector frame
   
   
   
    SVAR RootString 
    SVAR  fiterr_String,bgshiftstring
    variable ipost, ipcount, i,j,nxpost,nypost, npost
    variable plotsperpage = 7
    variable plotsoncurrentpage = 0
    variable bgshiftflag=1
    string nbname
    string xstr, MSDstr, MSDrstr, MSDtstr
    nxpost = astrw[3]
    nypost = astrw[4]
    npost = nxpost* nypost
    
    variable imax,maxflag
    variable islope,slopeflag
    variable L
    variable despikecoef = 4    //use this coef to determine how far the value from avg then it's called a spike  6/2/2014
    
    
    
    
    variable iframe1, iframe2, iframe
    iframe1=dimoffset(datawave,2)
    iframe2=dimsize(datawave, 2)+iframe1
    L=floor((iframe2-iframe1)/5)
    

    
    if (spikeflag==1)
 
       make/O/N=(iframe2) despiketemp   
       duplicate/O datawave despikewave
       
       for (ipost=0; ipost<npost; ipost+=1)
            despiketemp=datawave[ipost][kfx0][p]
            wavestats/Q despiketemp
            for (iframe=iframe1; iframe<iframe2; iframe+=1)                 
                 if ( (despiketemp[iframe] -V_avg) > despikecoef * V_sdev )
                     if (iframe==iframe1)
                     despikewave[ipost][kfx0][iframe] = 2* datawave[ipost][kfx0][iframe+1] - datawave[ipost][kfx0][iframe+2]
                      elseif (iframe==(iframe2-1))
                        despikewave[ipost][kfx0][iframe] = 2* datawave[ipost][kfx0][iframe-1] - datawave[ipost][kfx0][iframe-2]
                      else
                       despikewave[ipost][kfx0][iframe] = (datawave[ipost][kfx0][iframe-1]+datawave[ipost][kfx0][iframe+1] ) /2
                 endif
                 endif
                 
                 endfor
                 endfor
//                 
        for (ipost=0; ipost<npost; ipost+=1)
            despiketemp=datawave[ipost][kfy0][p]
            wavestats/Q despiketemp
            for (iframe=iframe1; iframe<iframe2; iframe+=1)
                 
                 if ( (despiketemp[iframe] -V_avg) > despikecoef * V_sdev )
                     if (iframe==iframe1)
                     despikewave[ipost][kfy0][iframe] = 2* datawave[ipost][kfy0][iframe+1] - datawave[ipost][kfy0][iframe+2]
                      elseif (iframe==(iframe2-1))
                        despikewave[ipost][kfy0][iframe] = 2* datawave[ipost][kfy0][iframe-1] - datawave[ipost][kfy0][iframe-2]
                      else
                       despikewave[ipost][kfy0][iframe] = (datawave[ipost][kfy0][iframe-1]+datawave[ipost][kfy0][iframe+1] ) /2
                 endif
                 endif
                 
                 endfor
                 endfor                
                 
                 duplicate/O despikewave datawave
                 killwaves despikewave
                 
      endif
       
    
    
    if (notebookflag == 1)
     //sprintf NBname, "NBPlotALot_"+RootString+"_"+"%d%d", param1,param2
     sprintf NBname, "NBPlotMSD"
    DoWindow/K $NBname
    NewNotebook/K=1/F=1/N=$NBname
endif
    ipcount =0
    
    make/N=((postf-posti),2)/O MSDstats
    
    
    wave MSD, MSDr, MSDt 
    wave W_coef 
    variable MSDmax
    
    imax=0
    islope=0
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
        sprintf MSDstr,"Plot%dMSD",plotsoncurrentpage
        sprintf MSDrstr, "Plot%dMSDr", plotsoncurrentpage
        sprintf MSDtstr, "Plot%dMSDt", plotsoncurrentpage
    
 	if(cellpostflag==1)
 		
 		if(maskwave[i][j]==Kcell)
               
              if (directionflag == 0)  		
 			   
      		  dispfitres4EB(datawave,ipost,2, "X",xstr, $Fiterr_String,$bgshiftstring,plotsoncurrentpage,maskwave[i][j], i,j,timeflag,cellpostflag,bgshiftflag,framerate, spikeflag)
        	  calcMSD(datawave,$bgshiftstring,ipost,MSDstr,plotsoncurrentpage+plotsperpage,maskwave[i][j],i,j,timeflag,cellpostflag,maxt, slopet, framerate, spikeflag, 0 )             

                                
   	       AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $xstr
      	       AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $MSDstr

             elseif (directionflag ==1)
                
                calcMSD(datawave,$bgshiftstring,ipost,MSDrstr,plotsoncurrentpage,maskwave[i][j],i,j,timeflag,cellpostflag,maxt, slopet, framerate, spikeflag, 1)
                calcMSD(datawave,$bgshiftstring,ipost,MSDtstr,plotsoncurrentpage+plotsperpage,maskwave[i][j],i,j,timeflag,cellpostflag,maxt, slopet, framerate, spikeflag,2)   	       

   	       AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $MSDrstr
      	       AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $MSDtstr
 
             endif
                   	       
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

      //temporary try by changing kfx0 to 2    

         if (directionflag == 0)


               dispfitres4EB(datawave,ipost,2, "X",xstr,$Fiterr_String,$bgshiftstring,plotsoncurrentpage,maskwave[i][j], i,j,timeflag,cellpostflag,bgshiftflag,framerate, spikeflag)
        	 calcMSD(datawave,$bgshiftstring, ipost,MSDstr,plotsoncurrentpage+plotsperpage,maskwave[i][j],i,j,timeflag,cellpostflag,maxt,slopet,framerate, spikeflag, 0)

    
   	       AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $xstr      	       
      	       AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $MSDstr



          elseif (directionflag ==1)
                
              calcMSD(datawave,$bgshiftstring,ipost,MSDrstr,plotsoncurrentpage,maskwave[i][j],i,j,timeflag,cellpostflag,maxt, slopet, framerate, spikeflag,1 )
              calcMSD(datawave,$bgshiftstring,ipost,MSDtstr,plotsoncurrentpage+plotsperpage,maskwave[i][j],i,j,timeflag,cellpostflag,maxt, slopet, framerate, spikeflag,2)   	       

   	       AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $MSDrstr
      	       AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $MSDtstr
 
              endif

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
         
//    redimension/N=(max(islope, imax),2) MSDstats    
    
//    make/O/N=(imax,2) MSDmask  
//    for (i=0; i<imax; i+=1)
//         MSDmask[i][0]=datawave[MSDstats[i]][kfx0][0]
//         MSDmask[i][1]=datawave[MSDstats[i]][kfy0][0]
//         endfor
    
    
    if (notebookflag==1)
	   SaveNotebook/P=ToPractice $NBname        //added by sasha, July 13  '06
	   //if (printFlag==1)
	  	//	PrintNotebook $NBname  //Modified Jul 27
	  // endif
     endif
     
     
  end
  
  
  
  
  
  
function calcMSD(datawave,bgshiftwave, npost,swname,nplot,posttype, m,n ,timeflag,cellpostflag, maxt, slopet,framerate, spikeflag, directionindex)
wave datawave          // fitresults wave
wave bgshiftwave     //background shift
variable npost       // number of post to be plotted
variable nplot       // number of plot
variable posttype // type of post for coloring graph
variable m,n          // matrix indices of post
variable timeflag  // =1 for plots vs time
variable cellpostflag
string swname       // string of Plot (graph) name
variable maxt
variable slopet
variable framerate
variable spikeflag
variable directionindex

variable i,j
variable pixratio=0.125


string s1, s2,dws
 wave timewave
 
 sprintf dws, "dw%d",nplot
 
 variable psym
 variable rval = 0, gval = 0, bval = 0
variable  DX, DY
variable  L, DataR
variable maxflag
variable slopeflag


DataR=dimsize(datawave,2)
//L=floor(DataR/5)
L=floor(DataR/4)

wave MSD, MSDxwave, MSDr, MSDt
//change made on 9/20/2015, rescale MSD
wave MSD_ed
    
make/O/N=(L-1) MSDtemp
    
if (directionindex == 0 )
    MSDtemp=MSD_ed[p][npost]
elseif (directionindex ==1 ) 
    MSDtemp = MSDr[p][npost]
elseif (directionindex ==2)
     MSDtemp = MSDt[p][npost]
endif

     
setscale/P x 0, (1/framerate), "", MSDtemp

duplicate/O MSDtemp $dws



//DataR=dimsize(datawave,2)
//DataR=6000
//L=DataR/5
//make/O/N=(L-1) MSD


//for(i=0;i<=L-2;i+=1)
 //    MSD[i]=0
  //  for(j=0;j+i+1<DataR;j+=1)
   //     DX=(datawave[npost][kfx0][j+i+1]-bgshiftwave[j+i+1][0]) - (datawave[npost][kfx0][j]-bgshiftwave[j][0])
  //      DY=(datawave[npost][kfy0][j+i+1]-bgshiftwave[j+i+1][1]) - (datawave[npost][kfy0][j]-bgshiftwave[j][1])
   //     MSD[i]=MSD[i]+DX^2+DY^2
    //    EndFor
     //if(j!=0)
  //   MSD[i]=MSD[i]/j
  //   endif

     
//EndFor     

//     MSD=MSD*(pixratio^2)

//setscale/P x 0, (1/framerate),"",MSD   
     
//duplicate/O MSD  $dws     

    if (posttype == kGuide)
        psym = 19
    elseif (posttype == kCell)
        psym = 8
        rval = 50000
    elseif (posttype == kEmpty)
        psym = 6
         bval = 50000
    elseif (posttype == kIgnore)
        psym = 16
        gval = 50000
    elseif (posttype == kWirepost)
        psym = 21
        rval = 30000
        bval =  30000
    endif
    
    dowindow/K $swname

    if (timeflag == 1)
        display/N=$swname $dws vs timewave
    else
         display/K=1/N = $swname $dws vs MSDxwave  	
         
    endif 
    
    SetAxis left, wavemin($dws),wavemax($dws)

    sprintf s1, "Post %d (%d,%d)\rType %d",npost,m,n,posttype
    TextBox/C/B=1/N=text0/A=LT s1
    ModifyGraph tick=2,mirror=1,standoff=0
    ModifyGraph mode=4,marker=psym,rgb=(rval,gval,bval)
    ModifyGraph log=1,tick=2,mirror=1,standoff=0
    ModifyGraph margin(left)=50,margin(bottom)=23,margin(top)=2,margin(right)=5;
    ModifyGraph width=165,height=75     


end


//currently in use
function calculateMSD(datawave, bgshiftwave, astrw, maskwave, cellflag, framerate, spikeflag, filterflag, freq)
wave datawave
wave bgshiftwave
wave astrw
wave maskwave
//wave zerow     // store undeflected post positions
variable cellflag
variable framerate
variable spikeflag
variable filterflag
variable freq

variable nxpost, nypost, npost, ipost
variable iframe1, iframe2, iframe
variable i, j 
variable dstart, pixratio=0.125
variable DataR, L
variable DX, DY
variable xdisp, ydisp     // displacement from untdeflected postions
variable directx, directy   // direction vector of local vector coordinates


nxpost=astrw[3]
nypost=astrw[4]
npost = nxpost*nypost

iframe1=dimoffset(datawave,2)
iframe2=dimsize(datawave, 2)+ iframe1

DataR=iframe2-iframe1
//L=floor(DataR/5)
L=floor(DataR/4)


if (filterflag==1)
   duplicate/o MSD MSD_orig
endif

make/O/N=(L-1, npost) MSD, MSD_ed , MSDr, MSDt

make/O/N=(L-1)  MSDxwave

for (ipost=0; ipost<npost; ipost+=1)


         
         make/O/N=(DataR)  dwtempx, dwtempy
         dwtempx=datawave[ipost][kfx0][p]
         dwtempy=datawave[ipost][kfy0][p]
        
    
         for( iframe=iframe1; iframe<iframe2; iframe+=1)
              dwtempx[iframe]-=bgshiftwave[iframe][0]
              dwtempy[iframe]-=bgshiftwave[iframe][1]
              endfor
              
              
//         xdisp = mean(dwtempx) - zerow[ipost][0][0]
//         ydisp = mean(dwtempy) - zerow[ipost][1][0]
   
//        directx = xdisp/sqrt(xdisp^2+ydisp^2)
//        directy = ydisp/sqrt(xdisp^2+ydisp^2)         
         
              
         dstart=dwtempx[0]
         dwtempx=dwtempx-dstart
         dstart=dwtempy[0]
         dwtempy=dwtempy- dstart
         dwtempx=dwtempx*pixratio
         dwtempy=dwtempy*pixratio 
         
     if (spikeflag==1)        
         wavestats/Q dwtempx
//         for( iframe=iframe1+1; iframe<iframe2; iframe+=1) 
//            if ( abs(dwtempx[iframe] -V_avg]) > 3 * V_sdev )
//                     if (iframe==iframe1)
//                     dwtempx[iframe] = 2* dwtempx[iframe+1] - dwtempx[iframe+2]
//                      elseif (iframe==(iframe2-1))
//                        dwtempx[iframe] = 2* dwtempx[iframe-1] - dwtempx[iframe-2]
//                      else
 //                      dwtempx[iframe] = (dwtempx[iframe-1]+dwtempx[iframe+1] ) /2
 //                endif
 //                endif
                 
 //                endfor              
             for (iframe=iframe1+1; iframe<iframe2; iframe+=1)  
                  if ( abs(dwtempx[iframe]-dwtempx[iframe-1]) > 3 * V_sdev ) 
                       dwtempx[iframe]=dwtempx[iframe-1]  
                   endif
                   endfor


         wavestats/Q dwtempy
//         for( iframe=iframe1; iframe<iframe2; iframe+=1) 
//            if ( abs(dwtempy[iframe] -V_avg) > 3 * V_sdev )
//                     if (iframe==iframe1)
//                     dwtempy[iframe] = 2* dwtempy[iframe+1] - dwtempy[iframe+2]
//                      elseif (iframe==(iframe2-1))
//                        dwtempy[iframe] = 2* dwtempy[iframe-1] - dwtempy[iframe-2]
//                      else
//                       dwtempy[iframe] = (dwtempy[iframe-1]+dwtempy[iframe+1] ) /2
//                 endif
//                 endif
                 
//                 endfor   

             for (iframe=iframe1+1; iframe<iframe2; iframe+=1)  
                  if ( abs(dwtempy[iframe]-dwtempy[iframe-1]) > 1.5 * V_sdev ) 
                       dwtempy[iframe]=dwtempy[iframe-1]  
                   endif
                   endfor
                   
   endif
                 
                          
   if (filterflag==1)
      notchfilter(dwtempx, freq)
      notchfilter(dwtempy, freq)

   endif
     
     
    
     
//for(i=0;i<=L-2;i+=1)
//     MSD[i][ipost]=0
//    for(j=0;j+i+1<DataR;j+=1)
//        DX=dwtempx[i+j+1]-dwtempx[j]
//        DY=dwtempy[i+j+1]-dwtempy[j]
//        MSD[i][ipost]=MSD[i][ipost]+DX^2+DY^2
//        EndFor
//     if(j!=0)
//     MSD[i][ipost]=MSD[i][ipost]/j
//     endif
//endfor

//manipulating array directly instead of each elements(12/3/14)
for (i =1; i<100; i+=1)
     make/o/n=(DataR-i)  wmsd, wmsdr, wmsdt
     duplicate/o/r = [0,datar-i-1] dwtempx, wmsdx1
     duplicate/o/r = [0,datar-i-1] dwtempy, wmsdy1
     duplicate/o/r = [i, datar-1] dwtempx, wmsdx2
     duplicate/o/r = [i, datar-1] dwtempy, wmsdy2     
     wmsd = (wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
//     wmsdr =( (wmsdx2[p] - wmsdx1[p] ) * directx + (wmsdy2[p] - wmsdy1[p]) * directy )^2
//     wmsdt =( (wmsdx2[p] - wmsdx1[p] ) * directy - (wmsdy2[p] - wmsdy1[p]) * directx ) ^2         
     MSD[i-1][ipost] = mean(wmsd)
//     MSDr[i-1][ipost] = mean(wmsdr)
//     MSDt[i-1][ipost] = mean(wmsdt)          
     MSDxwave[i-1] = i/framerate
     endfor
     
for (i=100; i<L; i+=10)   // i is the lag time and j is the index
     j = 100+ (i-100)/10
     make/o/n=(DataR-i)  wmsd, wmsdr, wmsdt
     duplicate/o/r = [0,datar-i-1] dwtempx, wmsdx1
     duplicate/o/r = [0,datar-i-1] dwtempy, wmsdy1
     duplicate/o/r = [i, datar-1] dwtempx, wmsdx2
     duplicate/o/r = [i, datar-1] dwtempy, wmsdy2     
     wmsd = (wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
//     wmsdr =( (wmsdx2[p] - wmsdx1[p] ) * directx + (wmsdy2[p] - wmsdy1[p]) * directy )^2
//     wmsdt =( (wmsdx2[p] - wmsdx1[p] ) * directy - (wmsdy2[p] - wmsdy1[p]) * directx ) ^2          
     MSD[j-1][ipost] = mean(wmsd)
//     MSDr[j-1][ipost] = mean(wmsdr)
//     MSDt[j-1][ipost] = mean(wmsdt)
     MSDxwave[j-1] = i/framerate
     endfor          


     
     
     
endfor

redimension/N =(j) MSDxwave 
redimension/N=(j,-1) MSD, MSDr, MSDt     

//added rescale MSD here   9/20/2015
//wave MSDfloor
variable MSDfloorV

make/N=(j)/O   MSDfloor

for (ipost =0; ipost<npost; ipost+=1)
    MSDfloor = MSD[p][ipost]
    MSDfloorV = mean(MSDfloor, 0, 10)
    MSD_ed[][ipost] = MSD[p][ipost]/MSDfloorV
    
endfor 




setscale/P x 0, (1/framerate), "", MSD

end



function calculateMSD_vshort(datawave, astrw, maskwave, cellflag, framerate, spikeflag, filterflag, freq)
wave datawave
//wave bgshiftwave
wave astrw
wave maskwave
//wave zerow     // store undeflected post positions
variable cellflag
variable framerate
variable spikeflag
variable filterflag
variable freq

variable nxpost, nypost, npost, ipost
variable iframe1, iframe2, iframe
variable i, j 
variable dstart, pixratio=0.125
variable DataR, L
variable DX, DY
variable xdisp, ydisp     // displacement from untdeflected postions
variable directx, directy   // direction vector of local vector coordinates


nxpost=astrw[3]
nypost=astrw[4]
npost = nxpost*nypost

iframe1=dimoffset(datawave,2)
iframe2=dimsize(datawave, 2)+ iframe1

DataR=iframe2-iframe1
//L=floor(DataR/5)
L=floor(DataR/4)


if (filterflag==1)
   duplicate/o MSD MSD_orig
endif

make/O/N=(L-1, npost) MSD, MSD_ed , MSDr, MSDt

make/O/N=(L-1)  MSDxwave

for (ipost=0; ipost<npost; ipost+=1)


         
         make/O/N=(DataR)  dwtempx, dwtempy
         dwtempx=datawave[ipost][kfx0][p]
         dwtempy=datawave[ipost][kfy0][p]
        
    
//         for( iframe=iframe1; iframe<iframe2; iframe+=1)
//              dwtempx[iframe]-=bgshiftwave[iframe][0]
//              dwtempy[iframe]-=bgshiftwave[iframe][1]
//              endfor
              
              
//         xdisp = mean(dwtempx) - zerow[ipost][0][0]
//         ydisp = mean(dwtempy) - zerow[ipost][1][0]
   
//        directx = xdisp/sqrt(xdisp^2+ydisp^2)
//        directy = ydisp/sqrt(xdisp^2+ydisp^2)         
         
              
         dstart=dwtempx[0]
         dwtempx=dwtempx-dstart
         dstart=dwtempy[0]
         dwtempy=dwtempy- dstart
         dwtempx=dwtempx*pixratio
         dwtempy=dwtempy*pixratio 
         
     if (spikeflag==1)        
         wavestats/Q dwtempx
//         for( iframe=iframe1+1; iframe<iframe2; iframe+=1) 
//            if ( abs(dwtempx[iframe] -V_avg]) > 3 * V_sdev )
//                     if (iframe==iframe1)
//                     dwtempx[iframe] = 2* dwtempx[iframe+1] - dwtempx[iframe+2]
//                      elseif (iframe==(iframe2-1))
//                        dwtempx[iframe] = 2* dwtempx[iframe-1] - dwtempx[iframe-2]
//                      else
 //                      dwtempx[iframe] = (dwtempx[iframe-1]+dwtempx[iframe+1] ) /2
 //                endif
 //                endif
                 
 //                endfor              
             for (iframe=iframe1+1; iframe<iframe2; iframe+=1)  
                  if ( abs(dwtempx[iframe]-dwtempx[iframe-1]) > 3 * V_sdev ) 
                       dwtempx[iframe]=dwtempx[iframe-1]  
                   endif
                   endfor


         wavestats/Q dwtempy
//         for( iframe=iframe1; iframe<iframe2; iframe+=1) 
//            if ( abs(dwtempy[iframe] -V_avg) > 3 * V_sdev )
//                     if (iframe==iframe1)
//                     dwtempy[iframe] = 2* dwtempy[iframe+1] - dwtempy[iframe+2]
//                      elseif (iframe==(iframe2-1))
//                        dwtempy[iframe] = 2* dwtempy[iframe-1] - dwtempy[iframe-2]
//                      else
//                       dwtempy[iframe] = (dwtempy[iframe-1]+dwtempy[iframe+1] ) /2
//                 endif
//                 endif
                 
//                 endfor   

             for (iframe=iframe1+1; iframe<iframe2; iframe+=1)  
                  if ( abs(dwtempy[iframe]-dwtempy[iframe-1]) > 1.5 * V_sdev ) 
                       dwtempy[iframe]=dwtempy[iframe-1]  
                   endif
                   endfor
                   
   endif
                 
                          
   if (filterflag==1)
      notchfilter(dwtempx, freq)
      notchfilter(dwtempy, freq)

   endif
     
     
    
     
//for(i=0;i<=L-2;i+=1)
//     MSD[i][ipost]=0
//    for(j=0;j+i+1<DataR;j+=1)
//        DX=dwtempx[i+j+1]-dwtempx[j]
//        DY=dwtempy[i+j+1]-dwtempy[j]
//        MSD[i][ipost]=MSD[i][ipost]+DX^2+DY^2
//        EndFor
//     if(j!=0)
//     MSD[i][ipost]=MSD[i][ipost]/j
//     endif
//endfor

//manipulating array directly instead of each elements(12/3/14)
for (i =1; i<100; i+=1)
     make/o/n=(DataR-i)  wmsd, wmsdr, wmsdt
     duplicate/o/r = [0,datar-i-1] dwtempx, wmsdx1
     duplicate/o/r = [0,datar-i-1] dwtempy, wmsdy1
     duplicate/o/r = [i, datar-1] dwtempx, wmsdx2
     duplicate/o/r = [i, datar-1] dwtempy, wmsdy2     
     wmsd = (wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
//     wmsdr =( (wmsdx2[p] - wmsdx1[p] ) * directx + (wmsdy2[p] - wmsdy1[p]) * directy )^2
//     wmsdt =( (wmsdx2[p] - wmsdx1[p] ) * directy - (wmsdy2[p] - wmsdy1[p]) * directx ) ^2         
     MSD[i-1][ipost] = mean(wmsd)
//     MSDr[i-1][ipost] = mean(wmsdr)
//     MSDt[i-1][ipost] = mean(wmsdt)          
     MSDxwave[i-1] = i/framerate
     endfor
     
for (i=100; i<L; i+=10)   // i is the lag time and j is the index
     j = 100+ (i-100)/10
     make/o/n=(DataR-i)  wmsd, wmsdr, wmsdt
     duplicate/o/r = [0,datar-i-1] dwtempx, wmsdx1
     duplicate/o/r = [0,datar-i-1] dwtempy, wmsdy1
     duplicate/o/r = [i, datar-1] dwtempx, wmsdx2
     duplicate/o/r = [i, datar-1] dwtempy, wmsdy2     
     wmsd = (wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
//     wmsdr =( (wmsdx2[p] - wmsdx1[p] ) * directx + (wmsdy2[p] - wmsdy1[p]) * directy )^2
//     wmsdt =( (wmsdx2[p] - wmsdx1[p] ) * directy - (wmsdy2[p] - wmsdy1[p]) * directx ) ^2          
     MSD[j-1][ipost] = mean(wmsd)
//     MSDr[j-1][ipost] = mean(wmsdr)
//     MSDt[j-1][ipost] = mean(wmsdt)
     MSDxwave[j-1] = i/framerate
     endfor          


     
     
     
endfor

redimension/N =(j) MSDxwave 
redimension/N=(j,-1) MSD, MSDr, MSDt     

//added rescale MSD here   9/20/2015
//wave MSDfloor
variable MSDfloorV

make/N=(j)/O   MSDfloor

for (ipost =0; ipost<npost; ipost+=1)
    MSDfloor = MSD[p][ipost]
    MSDfloorV = mean(MSDfloor, 0, 10)
    MSD_ed[][ipost] = MSD[p][ipost]/MSDfloorV
    
endfor 




setscale/P x 0, (1/framerate), "", MSD

end






function  displayallmsd(MSD)
wave MSD
wave idcell

variable L
variable npost,ipost
variable h=0
string MSDdisstr

wave msdxwave, linmaskval_intercept_int

L=dimsize(MSD,0)
npost=dimsize(MSD,1)

dowindow/K msdall

make/O/N=(L) displaytemp, displaytemp1
displaytemp=MSD[p][0]

display/N=msdall displaytemp vs msdxwave
modifygraph rgb(displaytemp) = (0,0,50000)


for (ipost=1; ipost<npost; ipost+=1)
//     if (idcell[h]==ipost)
     sprintf MSDdisstr,  "MSD%d", ipost
//     make/O/N=(L) $MSDdisstr
//     displaytemp1=MSD[p][ipost]
//     duplicate/O displaytemp1 $MSDdisstr
//     appendtograph $MSDdisstr
    if (linmaskval_intercept_int[ipost] ==1)
     appendtograph/W=msdall  MSD[][ipost]/TN = $MSDdisstr vs msdxwave

       modifygraph rgb($MSDdisstr) = (500000,0,0)
     h+=1
//     else
//        modifygraph rgb($MSDdisstr) = (0,0,50000)
     endif
     
 endfor
 
 end
 
 
 //using noise floor subtracted MSD for calculation
 //MSD statistics
 function MSDstatsmap(MSD, maskwave, astrw, mapflag)
 wave MSD
 wave maskwave
 wave astrw
 variable mapflag
 
variable nxpost, nypost
variable i, j 
variable npost 
variable ipost
variable MSDmax
variable slope
variable MSDlength
npost = dimsize(MSD, 1)
MSDlength = dimsize(MSD,0)

wave W_coef

nxpost = astrw[3] 
nypost = astrw[4]

make/N=(2,npost)/O MSDstats
make/N=(MSDlength)/O MSDt, MSDtlog, timelog

wave MSDxwave

timelog = log(MSDxwave)

for (ipost=0; ipost<npost; ipost+=1)
     i = floor(ipost/nypost)
     j = mod(ipost, nypost)
     MSDt = MSD[p][ipost] 
     MSDt = MSD[p][ipost] - mean(MSDt,0,4)
     if (maskwave[i][j] != 3)
        MSDstats[0][ipost] = wavemax(MSDt)
        MSDtlog = log(MSDt)
     //trying new fitting range
        CurveFit/Q /W=2/NTHR=0 line  MSDtlog[20,100]  /X=timelog /D 
        slope = W_coef[1]
//     if (slope<0)
//        slope = 0 
//     endif
        MSDstats[1][ipost] = slope
      else 
         MSDstats[0][ipost] = 0
         MSDstats[1][ipost] = 0
      endif
     
     
endfor

dowindow/K MSD_heatmap_disp
newimage/N=MSD_heatmap_disp M_rotatedimage
appendtograph/T/W=MSD_heatmap_disp linmasky vs linmaskx
modifygraph mode=3, marker=0
modifygraph zcolor(linmasky) = {MSDstats[mapflag][], *, *, Rainbow, 0}

end


function MSDsumary(MSDstats, astrw, maskwave)
wave MSDstats, astrw, maskwave

variable npost, nxpost, nypost
variable postcountc=0, postcount=0
variable sumslopec=0, summaxc=0, sumslope=0, summax=0
variable i, j 
nxpost =astrw[3]
nypost =astrw[4]
npost = nxpost*nypost

make/N=4/O MSDsummary

variable ipost

for (ipost=1; ipost<npost; ipost+=1) 
     i = floor(ipost/nypost)
     j = mod(ipost, nypost)
         if (maskwave[i][j] ==1)
            postcountc+=1
            sumslopec+= MSDstats[0][ipost]
            summaxc+=MSDstats[1][ipost]
         elseif(maskwave[i][j]==0)
           postcount+=1
           sumslope+=MSDstats[0][ipost]
           summax+=MSDstats[1][ipost]
      endif
endfor

MSDsummary[0] = sumslopec/postcountc
MSDsummary[1] = summaxc/postcountc
MSDsummary[2] = sumslope/postcount
MSDsummary[3] = summax/postcount

end


function singlemsd(dataw)
wave dataw

variable i, j 
variable L
variable len = dimsize(dataw, 0)
variable DX, DY
L = floor(len/5)
make/N=(L)/O   datawMSD=0

for (i = 0; i<L; i+=1)
     datawMSD[i] = 0
     for (j= 0; j+i+1<len; j+=1)
          DX = dataw[i+j+1] - dataw[j]
//          DY = dataw[i+j+1] - dataw[j]
          datawMSD[i] = datawMSD[i] + DX^2 //+DY^2
      endfor
      
   if (j!=0)
   datawMSD[i] = datawMSD[i]/j
   endif
endfor

end.







function rescaleMSD(MSD)
wave MSD

variable L = dimsize(MSD,0)
variable npost = dimsize(MSD,1)
variable ipost, i 
variable MSDmin

make/N=(L)/O MSDtemp



for (ipost = 0; ipost<npost; ipost+=1)
     MSDtemp = MSD[p][ipost]
//     MSDmin = wavemin(MSDtemp)
     MSdmin = MSDtemp[0]
     for ( i =0; i<L; i+=1)
          MSD[i][ipost] = MSD[i][ipost] - MSDmin 
      endfor
      
endfor

end

function rescalesingleMSD(MSD,ipost)
wave MSD
variable ipost

variable L = dimsize(MSD,0)
variable npost = dimsize(MSD,1)
variable  i 
variable MSDmin

make/N=(L)/O MSDtemp



     MSDtemp = MSD[p][ipost]
     MSDmin = wavemin(MSDtemp)
//     MSDmin = MSDtemp[0]
     for ( i =0; i<L; i+=1)
          MSD[i][ipost] = MSD[i][ipost] - MSDmin 
      endfor
      


end


//added 3/7/2017 for calculating individual magnetic posts MSD at low frequency
function calcactiveavgMSD(freq,postnum)
wave freq
variable postnum
SVAR FitRes_String
//wave datawave

variable numofavg=4    //number of averaged frequency, set to be 0.1, 0.2, 0.5,0.8 by default
variable ifreq
variable i,j
string fitname, fitaffix
variable movielength 
variable MSDlength
wave datawMSD


for ( ifreq = 0 ; ifreq<numofavg; ifreq+=1)
       sprintf fitaffix "_ed_0_%dhz" (freq[ifreq]*10) 
       fitname = fitres_string+fitaffix
       movielength = dimsize($fitname,2)
       make/o/n=(movielength)  temptracewave=0 
       duplicate/o $fitname datawave
       temptracewave = datawave[postnum][0][p] - datawave[postnum][0][0]
       temptracewave = temptracewave*125
       notchfilter2(temptracewave,freq[ifreq])
       notchfilter3(temptracewave,7)
       singlemsd(temptracewave)
       if(ifreq==0)
           MSDlength = dimsize(datawMSD,0)
           make/o/n=(MSDlength) MSDavg =0
       endif    
           MSDavg += log(datawMSD)
endfor

MSDavg = 10^(MSDavg/numofavg)

end           
           

function MSDslopehisto()



wave MSDstats
wave linmaskval
variable ngap = 0.1

variable npost = dimsize(MSDstats,1)
variable i
make/o/n=(npost) MSDslope=0

wave MSDslope_histo1

for(i=0; i<npost; i+=1)
     if (linmaskval[i] ==1)
         MSDslope[i] = MSDstats[1][i]
     endif
 endfor
 
 wavestats/q MSDslope 
 variable slope_max = V_max
 variable binstart = 0.01 
 variable nobins = ceil(slope_max-binstart)/ngap
 
 
 histogram/B={binstart, ngap, nobins} MSDslope ,MSDslope_histo1
 
 end       
 
 
 
 
 
 //using color code for different cases
 
 function MSDslopevslag(avgend, subflag, avgflag)
 variable avgend     //end of points for noise floor average
 variable subflag     //choose whether display subtract noisefloor
 variable avgflag    //choose whether to use fitting method for calculating slope

 
 wave MSD
 wave linmaskval, linmaskval_intercept_int
 wave msdxwave
 
 duplicate/o msdxwave msdxwave_log
 msdxwave_log = log(msdxwave)
 
 duplicate/o MSD MSD_diff, MSD_sub_diff, msd_sub
 
// duplicate/o MSD MSD_sub
 variable npost = dimsize(MSD,1)
 variable msdlength = dimsize(MSD,0)
 make/o/n=(msdlength) msdtemp, msdtemp_sub,msdtemp_log, msdtemp_sub_log
 variable ipost
 variable ipcount = 0, plotsoncurrentpage = 0
 variable plotsperpage = 7
 wave msdxwave_diff
 
 variable startflag=0
 
 string  NBname, msdstr, slopestr
 
    sprintf NBname, "NBslopevslag"
    DoWindow/K $NBname
    NewNotebook/K=1/F=1/N=$NBname
 
 wave w_coef
 
 for (ipost=0; ipost<npost; ipost+=1)

 
         if (mod(ipcount, plotsperpage) == 0)
            plotsoncurrentpage = 0
            DoWindow/K LayoutXY
            NewLayout/K=1/N=LayoutXY
        endif
         sprintf msdstr,"Plot%dmsd",plotsoncurrentpage
         sprintf slopestr,"Plot%dslope",plotsoncurrentpage
 
        msdtemp = MSD[p][ipost]
        curvefit/NTHR=0/Q/W=2 power msdtemp[0,avgend] /X=MSDxwave /D
//        msdtemp_sub = msdtemp - mean(msdtemp,0,avgend)
        msdtemp_sub = msdtemp - w_coef[0]
        msdtemp_log = log(msdtemp)
        msdtemp_sub_log = log(msdtemp_sub)
        
        //using linear fitting to calcualte slope of small segment

        
        
        //using built-in function to calculate derivative
        differentiate msdtemp_log /x=msdxwave_log/D=msdtemp_diff
        differentiate msdtemp_sub_log /X=msdxwave_log/D=msdtemp_sub_diff

//trying alternative way of smoothing
//changed for new mask trial

//        if ( linmaskval[ipost]==1 && avgflag==1)
         if( linmaskval_intercept_int[ipost] !=0 && linmaskval_intercept_int[ipost] !=3 && avgflag ==1)
                           
         calcmsdslope(msdtemp_log, msdxwave_log)
         duplicate/o msdtempwave_diff msdtemp_diff
         calcmsdslope(msdtemp_sub_log, msdxwave_log)
        duplicate/o msdtempwave_diff msdtemp_sub_diff
        endif

//       if(linmaskval_intercept_int[ipost] !=0 && linmaskval_intercept_int[ipost] !=3 && avgflag ==2)      //changed temproaily 10/30/18
       if( linmaskval_intercept_int[ipost] !=3 && avgflag ==1)            
          smoothmsd(msdtemp_sub_log, msdxwave_log,startflag)
          duplicate/o msdtempwave_diff msdtemp_sub_diff
          duplicate/o msdtempwave_smooth  msdtemp_smooth
            startflag+=1
        endif
        
        
        
        MSD_diff[][ipost] = msdtemp_diff[p]
        MSD_sub_diff[][ipost] = msdtemp_sub_diff[p]
        msd_sub[][ipost] = msdtemp_sub[p]
        
//        if (linmaskval_intercept_int[ipost] !=0 && linmaskval_intercept_int[ipost] !=3 )           //changed temprarily 10/30/2018
          if (linmaskval_intercept_int[ipost] !=3 )
            if (subflag==0)
                  dispmsdslope(msdtemp,ipost, 0,msdstr, plotsoncurrentpage, msdxwave)
                  if (avgflag==1)
                  dispmsdslope(msdtemp_diff,ipost, 1, slopestr, plotsoncurrentpage+plotsperpage, msdxwave_diff)
                  else
                  dispmsdslope(msdtemp_diff,ipost, 1, slopestr, plotsoncurrentpage+plotsperpage, msdxwave)
                  endif
            else
                  dispmsdslope(msdtemp_sub,ipost, 0, msdstr, plotsoncurrentpage, msdxwave)
                  if (avgflag==1)
                  dispmsdslope(msdtemp_sub_diff,ipost, 1, slopestr, plotsoncurrentpage+plotsperpage, msdxwave_diff)
                  else
                  dispmsdslope(msdtemp_sub_diff,ipost, 1, slopestr, plotsoncurrentpage+plotsperpage, msdxwave)
                  endif                       
            endif
   	       AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $msdstr
      	       AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $slopestr            
              ipcount = ipcount+1
             plotsoncurrentpage =  plotsoncurrentpage + 1
                if ((plotsoncurrentpage == 7) || (ipost == npost))
                  DoWindow/F LayoutXY
             //if (printflag == 1)
                   //PrintLayout LayoutXY  //Modified Jul 27
                       Notebook $NBname scaling = {90,90}, picture = {LayoutXY, -1,1}
                       Notebook $NBname text = "\r"
               endif        
                       
           endif               
        
endfor


SaveNotebook/P=ToPractice $NBname



end


function cellmotdetect(datawave, fitlength,slopethresh,bg_thresh,recalc_msd_flag)
wave datawave
variable fitlength, slopethresh,bg_thresh,recalc_msd_flag

variable movielength = dimsize(datawave,2)
//variable MSDlength = fitlength/5
variable npost = dimsize(datawave,0)
//variable slopethresh = 0.5
variable magthresh = 1e-5*8^2
variable ipost, i

wave linmaskval, linmaskval_begin, linmaskval_end

duplicate/o linmaskval linmaskval_begin_auto, linmaskval_end_auto, linmaskval_intercept

wave datawmsd, w_coef, msdxwave_seg

if (recalc_msd_flag == 0)

make/o/n=(fitlength) begintempx, begintempy, endtempx, endtempy
make/o/n=(npost) beginslope, endslope, beginmag, endmag
for (ipost = 0; ipost<npost; ipost+=1)
      begintempx = datawave[ipost][0][p]
      begintempy = datawave[ipost][1][p]
      endtempx = datawave[ipost][0][movielength-fitlength+p]
      endtempy = datawave[ipost][1][movielength-fitlength+p]
      
      
      
      calcsinglemsdlong(begintempx, begintempy)
      duplicate/o datawmsd begintempmsd
      duplicate/o datawmsd datawmsd_log
      duplicate/o msdxwave_seg msdxwave_seglog
      datawmsd_log = log(begintempmsd)
      msdxwave_seglog = log(msdxwave_seg)
      CurveFit/Q /W=2/NTHR=0 line  datawmsd_log[120,140]  /X=msdxwave_seglog /D 
      beginslope[ipost] = w_coef[1]
      beginmag[ipost] = begintempmsd[110]
      
      
      
      
      
      calcsinglemsdlong(endtempx, endtempy)
      duplicate/o datawmsd endtempmsd
      datawmsd_log = log(endtempmsd)
      CurveFit/Q /W=2/NTHR=0 line  datawmsd_log[100,120]  /X=msdxwave_seglog /D 
      endslope[ipost] = w_coef[1]
      endmag[ipost] = endtempmsd[110]
      
      if(beginslope[ipost]>slopethresh && linmaskval[ipost]!=3)
//        if(beginmag[ipost] >magthresh && linmaskval[ipost]!=3)
          linmaskval_begin_auto[ipost] =1
       elseif(beginslope[ipost]<bg_thresh && linmaskval[ipost]!=3)
           linmaskval_begin_auto[ipost]=0
       else
         	linmaskval_end_auto[ipost] = 3           
        endif
        
        if(endslope[ipost]>slopethresh && linmaskval[ipost]!=3)
//        if(endmag[ipost] >magthresh && linmaskval[ipost]!=3)
            linmaskval_end_auto[ipost] =1
         elseif (endslope[ipost]<bg_thresh && linmaskval[ipost]!=3)
            linmaskval_end_auto[ipost] =0
         else
         	linmaskval_end_auto[ipost] = 3
         endif

      
endfor

elseif (recalc_msd_flag == 1)
	wave beginslope, endslope, beginmag, endmag
	
for (ipost = 0; ipost<npost; ipost+=1)
      if(beginslope[ipost]>slopethresh && linmaskval[ipost]!=3)
//        if(beginmag[ipost] >magthresh && linmaskval[ipost]!=3)
          linmaskval_begin_auto[ipost] =1
       elseif(beginslope[ipost]<bg_thresh && linmaskval[ipost]!=3)
           linmaskval_begin_auto[ipost]=0
       else
         	linmaskval_end_auto[ipost] = 3           
        endif
        
        if(endslope[ipost]>slopethresh && linmaskval[ipost]!=3)
//        if(endmag[ipost] >magthresh && linmaskval[ipost]!=3)
            linmaskval_end_auto[ipost] =1
         elseif (endslope[ipost]<bg_thresh && linmaskval[ipost]!=3)
            linmaskval_end_auto[ipost] =0
         else
         	linmaskval_end_auto[ipost] = 3
         endif
endfor

endif


linmaskval_intercept = linmaskval_begin_auto*linmaskval_end_auto

end


function cellmotdetect_2(datawave, fitlength)
wave datawave
variable fitlength

variable movielength = dimsize(datawave,2)
//variable MSDlength = fitlength/5
variable npost = dimsize(datawave,0)
variable slopethresh = 0.5
variable magthresh = 1e-5*8^2
variable ipost, i

variable startflag = 0

variable Rmin, Rmax


//wave linmaskval, linmaskval_begin, linmaskval_end

//duplicate/o linmaskval linmaskval_begin_auto, linmaskval_end_auto

make/o/n=(npost) beginslope_2, endslope_2, beginmag_2, endmag_2, beginptime_2, endptime_2

wave datawmsd, w_coef, msdxwave_seg

Rmin = 80
Rmax = dimsize(msdxwave_seg,0)

make/o/n=(fitlength) begintempx, begintempy, endtempx, endtempy

wave msdtempwave_diff, msdtempwave_smooth, msdxwave_diff

variable ptimethresh = 0.8

for (ipost = 0; ipost<npost; ipost+=1)
      begintempx = datawave[ipost][0][p]
      begintempy = datawave[ipost][1][p]
      endtempx = datawave[ipost][0][movielength-fitlength+p]
      endtempy = datawave[ipost][1][movielength-fitlength+p]
      
      
      
      calcsinglemsdlong(begintempx, begintempy)
      duplicate/o datawmsd begintempmsd
      duplicate/o datawmsd datawmsd_log
      duplicate/o msdxwave_seg msdxwave_seglog
      datawmsd_log = log(begintempmsd)
      msdxwave_seglog = log(msdxwave_seg)
      CurveFit/Q /W=2/NTHR=0 line  datawmsd_log[100,120]  /X=msdxwave_seglog /D 
      beginslope_2[ipost] = w_coef[1]
      beginmag_2[ipost] = begintempmsd[110]
     
     

      
          smoothmsd(datawmsd_log, msdxwave_seglog,startflag)
          duplicate/o msdtempwave_diff msdtemp_sub_diff
          duplicate/o msdtempwave_smooth  msdtemp_smooth
            startflag+=1    
              
          wavestats/q/R=(Rmin, Rmax)  msdtemp_sub_diff
  
  
         if(V_max > 0.5)
          findlevel/Q/edge=2/R=(Rmin, Rmax) msdtemp_sub_diff, (V_max*ptimethresh)
          beginptime_2[ipost] = msdxwave_diff[round(V_levelx)]
          else 
          beginptime_2[ipost] = 0
          endif    
      
      
      calcsinglemsdlong(endtempx, endtempy)
      duplicate/o datawmsd endtempmsd
      datawmsd_log = log(endtempmsd)
      CurveFit/Q /W=2/NTHR=0 line  datawmsd_log[100,120]  /X=msdxwave_seglog /D 
      endslope_2[ipost] = w_coef[1]
      endmag_2[ipost] = endtempmsd[110]
      
      
      
             smoothmsd(datawmsd_log, msdxwave_seglog,startflag)
          duplicate/o msdtempwave_diff msdtemp_sub_diff
          duplicate/o msdtempwave_smooth  msdtemp_smooth
            startflag+=1    
              
          wavestats/q/R=(Rmin, Rmax)  msdtemp_sub_diff
  
  
         if(V_max > 0.5)
          findlevel/Q/edge=2/R=(Rmin, Rmax) msdtemp_sub_diff, (V_max*ptimethresh)
          endptime_2[ipost] = msdxwave_diff[round(V_levelx)]
          else 
          endptime_2[ipost] = 0          
          endif       
      
      
endfor


end


function cellmotdetect_seg(datawave, fitlength, slopethresh)
wave datawave
variable fitlength
variable slopethresh 

wave linmaskval

variable movielength = dimsize(datawave,2)
//variable MSDlength = fitlength/5
variable npost = dimsize(datawave,0)

variable magthresh = 1e-5*8^2
variable ipost, i

variable startflag = 0

variable Rmin, Rmax


//wave linmaskval, linmaskval_begin, linmaskval_end

//duplicate/o linmaskval linmaskval_begin_auto, linmaskval_end_auto



wave datawmsd, w_coef, msdxwave_seg

Rmin = 80
Rmax = dimsize(msdxwave_seg,0)

variable segnum, j

segnum = floor(movielength/fitlength)

make/o/n=(fitlength) segtempx, segtempy

make/o/n=(npost,segnum) segslope,  segmag, segptime
//make 2D mask 
make/o/n=(npost, segnum) linmaskval_seg

wave msdtempwave_diff, msdtempwave_smooth, msdxwave_diff

variable ptimethresh = 0.8

for (ipost = 0; ipost<npost; ipost+=1)
      for(j=0; j<segnum; j+=1)
      segtempx = datawave[ipost][0][j*fitlength+p]
      segtempy = datawave[ipost][1][j*fitlength+p]
      
      
      
      calcsinglemsdlong(segtempx, segtempy)
      duplicate/o datawmsd segtempmsd
      duplicate/o datawmsd datawmsd_log
      duplicate/o msdxwave_seg msdxwave_seglog
      datawmsd_log = log(segtempmsd)
      msdxwave_seglog = log(msdxwave_seg)
      if(fitlength > 1000)
      CurveFit/Q /W=2/NTHR=0 line  datawmsd_log[100,120]  /X=msdxwave_seglog 
      else
      CurveFit/Q /W=2/NTHR=0 line  datawmsd_log[80,100]  /X=msdxwave_seglog
      endif
      segslope[ipost][j] = w_coef[1]
      segmag[ipost][j] = segtempmsd[110]
      
      if(segslope[ipost][j]>slopethresh && linmaskval[ipost]!=3)
//        if(beginmag[ipost] >magthresh && linmaskval[ipost]!=3)
          linmaskval_seg[ipost][j] =1
       else
           linmaskval_seg[ipost][j]=0
        endif
      
      
          smoothmsd(datawmsd_log, msdxwave_seglog,startflag)
          duplicate/o msdtempwave_diff msdtemp_sub_diff
          duplicate/o msdtempwave_smooth  msdtemp_smooth
            startflag+=1    
              
          wavestats/q/R=(Rmin, Rmax)  msdtemp_sub_diff
  
  
         if(V_max > 0.5)
          findlevel/Q/edge=2/R=(Rmin, Rmax) msdtemp_sub_diff, (V_max*ptimethresh)
         segptime[ipost][j] = msdxwave_diff[round(V_levelx)]
          else 
         segptime[ipost][j] = 0
          endif       

      
endfor
      
endfor


end


function calcMSD_part_seg()

wave segmag
wave linmaskval_intercept_int

variable npost = dimsize(linmaskval_intercept_int,0)
variable ipost
variable nlb = 0
variable nhb = 0
variable pixelratio = 8

variable segnum = dimsize(segmag,1)
variable j 

make/o/n=(segnum,4) MSD_bistats_seg2
make/o/n=(npost) msd_mag_lb_seg, msd_mag_hb_seg

for (j = 0; j<segnum; j+=1)
nlb=0
nhb=0
for (ipost =0; ipost < npost; ipost+=1) 
      if (linmaskval_intercept_int[ipost] ==1)
          msd_mag_lb_seg[nlb] = log(segmag[ipost][j]/pixelratio^2) 
          nlb+=1
      elseif (linmaskval_intercept_int[ipost] == 6)
          msd_mag_hb_seg[nhb] = log(segmag[ipost][j]/pixelratio^2)
          nhb+=1
      endif

endfor


redimension/n=(nlb) msd_mag_lb_seg
redimension/n=(nhb) msd_mag_hb_seg


wavestats/q msd_mag_lb_seg
MSD_bistats_seg2[j][0] = 10^V_avg
msd_bistats_seg2[j][1] = 10^V_sem
wavestats/q msd_mag_hb_seg
MSD_bistats_seg2[j][2] = 10^V_avg
msd_bistats_seg2[j][3] = 10^V_sem



msd_mag_lb_seg = 10^(msd_mag_lb_seg)
msd_mag_hb_seg = 10^(msd_mag_hb_seg)

endfor

end


function sortptime()
wave msd_ptime
variable npost = dimsize(msd_ptime,0)

variable ipost

for (ipost = 0 ; ipost<npost; ipost+=1)
     if (msd_ptime[ipost]<1)
         msd_ptime[ipost] = 0 
     endif
     
endfor

end


function sortbeginvsend(maskwave)
wave maskwave

variable npost = dimsize(maskwave, 0 )
variable ipost
wave beginslope_2, endslope_2, beginptime_2, endptime_2
make/o/n=(npost) beginslope_lb, beginslope_hb, endslope_lb, endslope_hb, beginptime_lb, beginptime_hb, endptime_lb, endptime_hb
variable nlb=0, nhb=0

for (ipost=0; ipost<npost; ipost+=1)
      if (maskwave[ipost] == 1)
          beginslope_lb[nlb] = beginslope_2[ipost]
          endslope_lb[nlb] = endslope_2[ipost]
          beginptime_lb[nlb] = beginptime_2[ipost]
          endptime_lb[nlb] = endptime_2[ipost]
          nlb+=1
      elseif (maskwave[ipost] == 6)
          beginslope_hb[nhb] = beginslope_2[ipost]
          endslope_hb[nhb] = endslope_2[ipost]
          beginptime_hb[nhb] = beginptime_2[ipost]
          endptime_hb[nhb] = endptime_2[ipost]
          nhb+=1          
       endif    
endfor

redimension/N=(nlb) beginslope_lb, endslope_lb, beginptime_lb, endptime_lb
redimension/N=(nhb) beginslope_hb,endslope_hb, beginptime_hb, endptime_hb





end
             






function calcsinglemsdlong(dwtempx1, dwtempy1)
wave dwtempx1, dwtempy1

variable i=0, j=0 
variable L
variable len = dimsize(dwtempx1, 0)
variable DX, DY
L = floor(len/5)
make/N=(L)/O   datawMSD=0, msdxwave_seg = 0


     
for (i = 1; i<100; i+=1)
     duplicate/o/r = [0,len-i-1] dwtempx1, wmsdx1
     duplicate/o/r = [0,len-i-1] dwtempy1, wmsdy1
     duplicate/o/r = [i, len-1] dwtempx1, wmsdx2
     duplicate/o/r = [i, len-1] dwtempy1, wmsdy2     
     make/o/n=(len-i)  wmsd
     wmsd = (wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
     datawmsd[j] = mean(wmsd)
     msdxwave_seg[j] = i
     j=j+1
      endfor
      
for(i=100; i<L;i+=10)
     duplicate/o/r = [0,len-i-1] dwtempx1, wmsdx1
     duplicate/o/r = [0,len-i-1] dwtempy1, wmsdy1
     duplicate/o/r = [i, len-1] dwtempx1, wmsdx2
     duplicate/o/r = [i, len-1] dwtempy1, wmsdy2     
     make/o/n=(len-i)  wmsd
     wmsd = (wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
     datawmsd[j] = mean(wmsd) 
     msdxwave_seg[j] = i+1
     j=j+1
endfor


        
redimension/N=(j)      datawmsd , msdxwave_seg



end.


function findinteriorposts(linmaskval_intercept, posttype)
wave linmaskval_intercept
variable posttype

wave linmaskval_end_auto, linmaskval_begin_auto

duplicate/o linmaskval_intercept linmaskval_intercept_int     //mask for interior posts

wave zeropos, fitres_ed
duplicate/o fitres_ed  fitres_ed_sub


variable npost = dimsize(linmaskval_intercept,0)
variable ipost, i
variable movielength = dimsize(fitres_ed,2)

make/o/n=(movielength) sforce=0
variable kspring, forcethresh

//variable kspring = 15.7     //for m6
//variable kspring = 5.7          // for m10
//variable kspring = 22.3       //for m4
//variable forcethresh = 3         //for m6
//variable forcethresh = 2       //for m10
//variable forcethresh = 10            //for M4


if(posttype ==10)
    kspring = 5.7
    forcethresh = 2
    elseif(posttype ==6)
    kspring =15.7
    forcethresh =  5
    elseif(posttype ==4)
    kspring = 22.3
    forcethresh = 10
endif    


variable forcethreshl = 1
variable pratio = 0.125
fitres_ed_sub = (fitres_ed[p][q][r] - zeropos[p][q][0])*pratio
variable V_max, V_min
variable noll     //number of points below threshold
variable pthresh= 16000    //number of points below threshold that forms plateau

for(ipost = 0; ipost<npost; ipost+=1)
       sforce = sqrt(fitres_ed_sub[ipost][0][p]^2 + fitres_ed_sub[ipost][1][p]^2)
       wavestats/q sforce
       noll=0
//modified categorization into always highly deflected, alwayls weakly deflected, ones that have transition
      for(i = 0; i<movielength; i+=1)
             if (sforce[i]>forcethreshl/kspring && V_max>forcethresh/kspring)
                  noll+=1
             endif
       endfor           

//       if (V_max > forcethresh/kspring && linmaskval_intercept[ipost] ==1)
//       if (noll >pthresh && linmaskval_intercept[ipost] ==1)
         if (V_min >forcethreshl/kspring && V_max>forcethresh/kspring && linmaskval_intercept[ipost] ==1)
            linmaskval_intercept_int[ipost] = 6
//       elseif (noll < pthresh && noll>0 && linmaskval_intercept[ipost] ==1)
         elseif (V_min < forcethreshl/kspring && V_max>forcethresh/kspring && linmaskval_intercept[ipost] ==1)
            linmaskval_intercept_int[ipost] =7
       elseif ( V_max<=forcethresh/kspring && linmaskval_intercept[ipost] ==1)
            linmaskval_intercept_int[ipost] = 1
       elseif (linmaskval_end_auto[ipost]==0 && linmaskval_begin_auto[ipost] ==1)
            linmaskval_intercept_int[ipost] =-1
       elseif (linmaskval_end_auto[ipost] ==1 && linmaskval_begin_auto[ipost] ==0)    
             linmaskval_intercept_int[ipost] =-2
       elseif( linmaskval_end_auto[ipost] ==3 || linmaskval_begin_auto[ipost] ==3)
             linmaskval_intercept_int[ipost] =3      
       endif     
endfor


end     


function findinteriorposts_dis(datawave,linmaskval_intercept,posttype)
wave linmaskval_intercept, datawave
variable posttype
wave linmaskval_end_auto, linmaskval_begin_auto

duplicate/o linmaskval_intercept linmaskval_intercept_int     //mask for interior posts

wave zeropos
duplicate/o datawave  fitres_ed_sub



variable npost = dimsize(linmaskval_intercept,0)
variable ipost, i
variable movielength = dimsize(datawave,2)

make/o/n=(movielength) sforce=0

//variable kspring = 15.7     //for m6
//variable kspring = 5.7          // for m10
//variable kspring = 22.3       //for m4
//variable forcethresh = 5         //for m6
//variable forcethresh = 1.5       //for m10
//variable forcethresh = 10            //for M4

variable kspring, forcethresh
variable forcethreshl = 3




if(posttype ==10)
    kspring = 5.7
    forcethresh = 3
    forcethreshl = 1
    
    
    elseif(posttype == 5)
    kspring = 18.16857
    forcethresh =8
    forcethreshl = 2.5 
    
    
//change to M6 temporarily    
//    elseif(posttype == 5)
//    kspring = 18.16857
//    forcethresh =5
//    forcethreshl = 2
       
    elseif(posttype ==6)
//    kspring =15.7
//    kspring = 22.3         //M4
//     kspring = 18.16857   //M5
    kspring = 11.5
    forcethresh =  5
    forcethreshl = 2
        
    elseif(posttype ==7)
    kspring =11.5
    forcethresh =  4
    forcethreshl = 2

//change to M6 temporarily
//    elseif(posttype ==7)
//    kspring =11.5
//    forcethresh =  5
//    forcethreshl = 2
    
    elseif(posttype ==4)
    kspring = 22.3
    forcethresh = 10
    forcethreshl = 4   


//change to M6 temporarily
//    elseif(posttype ==4)
//    kspring = 22.3
//    forcethresh = 5
//    forcethreshl = 2        
    
    
endif    




variable pratio = 0.125
fitres_ed_sub = (datawave[p][q][r] - zeropos[p][q][1])*pratio
variable V_max, V_min
variable noll     //number of points below threshold
variable pthresh= 16000    //number of points below threshold that forms plateau

make/o/n=(npost)  pforce_avg
make/o/n=(npost, movielength) pforce


for(ipost = 0; ipost<npost; ipost+=1)
       sforce = sqrt(fitres_ed_sub[ipost][0][p]^2 + fitres_ed_sub[ipost][1][p]^2)
       wavestats/q sforce
       noll=0
       pforce[ipost][] = sforce[q]*kspring
       pforce_avg[ipost] = V_avg*kspring
//modified categorization into always highly deflected, alwayls weakly deflected, ones that have transition
      for(i = 0; i<movielength; i+=1)
             if (sforce[i]>forcethreshl/kspring && V_max>forcethresh/kspring)
                  noll+=1
             endif
       endfor           

//       if (V_max > forcethresh/kspring && linmaskval_intercept[ipost] ==1)
//       if (noll >pthresh && linmaskval_intercept[ipost] ==1)
         if (V_avg >forcethresh/kspring && linmaskval_intercept[ipost] ==1)
            linmaskval_intercept_int[ipost] = 6
//       elseif (noll < pthresh && noll>0 && linmaskval_intercept[ipost] ==1)
         elseif (V_avg < forcethresh/kspring && V_max>forcethreshl/kspring && linmaskval_intercept[ipost] ==1)
            linmaskval_intercept_int[ipost] =7
       elseif ( V_max<=forcethreshl/kspring && linmaskval_intercept[ipost] ==1)
            linmaskval_intercept_int[ipost] = 1
       elseif (linmaskval_end_auto[ipost]!=1 && linmaskval_begin_auto[ipost] ==1)
            linmaskval_intercept_int[ipost] =-1
       elseif (linmaskval_end_auto[ipost] ==1 && linmaskval_begin_auto[ipost] !=1)    
             linmaskval_intercept_int[ipost] =-2
 
       elseif( linmaskval_end_auto[ipost] ==0 && linmaskval_begin_auto[ipost] ==0)
             linmaskval_intercept_int[ipost] =0 
        
        else
        	 	linmaskval_intercept_int[ipost] =3	     
        
       endif     
endfor


end     



function findinteriorposts_seg(dataw, linmaskval_seg,posttype)
wave dataw
wave linmaskval_seg
variable posttype
wave linmaskval_end_auto, linmaskval_begin_auto

duplicate/o linmaskval_seg linmaskval_seg_int     //mask for interior posts

wave zeropos
duplicate/o dataw  fitres_ed_sub


variable npost = dimsize(linmaskval_seg,0)
variable ipost, i, iseg
variable movielength = dimsize(dataw,2)
variable segnum = dimsize(linmaskval_seg, 1)
variable seglength = floor(movielength/segnum)

//change mask type to 2D matrix


make/o/n=(movielength) sforce=0

variable kspring, forcethresh,forcethreshl

if(posttype ==4)
    kspring = 22.3
    forcethresh = 10
    forcethreshl = 4
    elseif(posttype == 5)
    kspring = 18.16857
    forcethresh = 8
    forcethreshl = 2.5    
    elseif(posttype ==6)
    kspring =15.7
    forcethresh =  5
    forcethreshl = 2
    elseif(posttype ==7)
    kspring =11.5
    forcethresh =  4
    forcethreshl = 2
    elseif(posttype ==10)
    kspring = 5.7
    forcethresh = 3            //1.5 originally
    forcethreshl = 1

endif 


//variable kspring = 15.7     //for m6
//variable kspring = 5.7          // for m10
//variable kspring = 22.3       //for m4
//variable forcethresh = 5         //for m6
//variable forcethresh = 1.5       //for m10
//variable forcethresh = 10            //for M4
//variable forcethreshl = 2

variable pratio = 0.125
fitres_ed_sub = (dataw[p][q][r] - zeropos[p][q][0])*pratio
variable V_max, V_min
variable noll     //number of points below threshold
variable pthresh= movielength*0.8    //number of points below threshold that forms plateau

make/o/n=(npost)  pforce_avg
make/o/n=(npost, movielength) pforce



for(ipost = 0; ipost<npost; ipost+=1)
       sforce = sqrt(fitres_ed_sub[ipost][0][p]^2 + fitres_ed_sub[ipost][1][p]^2)
       pforce[ipost][] = sforce*kspring    
     for(iseg = 0; iseg<segnum; iseg+=1)       
       wavestats/q/r=(iseg*seglength, (iseg+1)*seglength-1) sforce
       pforce_avg[ipost] = V_avg*kspring
 
//modified categorization into always highly deflected, alwayls weakly deflected, ones that have transition     

       if (V_avg >forcethresh/kspring && linmaskval_seg[ipost][iseg] ==1)
            linmaskval_seg_int[ipost][iseg] = 6
        elseif (V_avg < forcethresh/kspring && V_max>forcethreshl/kspring && linmaskval_seg[ipost][iseg] ==1)
            linmaskval_seg_int[ipost][iseg] =7
       elseif ( V_max<=forcethreshl/kspring && linmaskval_seg[ipost][iseg] ==1)
            linmaskval_seg_int[ipost][iseg] = 1
       endif    
       
       endfor 
        
endfor


end  






function dispmsdslope(datawave, ipost,npara,  swname, nplot, msdxwave)
wave datawave,msdxwave
variable ipost,nplot, npara
string swname
string dws,s1

sprintf dws, "dw%d", nplot
duplicate/o datawave $dws

wave linmaskval_intercept_int

variable rval=0, bval=0, gval=0

if(linmaskval_intercept_int[ipost] == -1)
    rval = 65280
    gval = 43520
elseif (linmaskval_intercept_int[ipost] == -2)
    rval = 50000
elseif (linmaskval_intercept_int[ipost] == 1)
    gval =50000
elseif (linmaskval_intercept_int[ipost] == 6)
    bval = 50000  
elseif( linmaskval_intercept_int[ipost] == 7)   
     rval = 36864
     gval = 14592
     bval = 58880
endif



    dowindow/K $swname
    display/K=1/N=$swname $dws vs msdxwave
if( npara == 1) //&& wavemax($dws)>4)
    SetAxis left, -0.5,2.5
    endif 
if( npara ==0)
     SetAxis left, 10^-8, 0.1 
     endif   

    sprintf s1, "Post %d",ipost
    TextBox/C/B=1/N=text0/A=LT s1
    ModifyGraph tick=2,mirror=1,standoff=0
    ModifyGraph mode=0,rgb=(rval,gval,bval) 
    if (npara == 0)
    ModifyGraph log=1,tick=2,mirror=1,standoff=0
    else
    ModifyGraph log(bottom)=1, tick=2,mirror=1,standoff=0
    endif
    ModifyGraph margin(left)=50,margin(bottom)=23,margin(top)=2,margin(right)=5;
    ModifyGraph width=165,height=75     

end

function calcmsdslope(msdtempwave, msdxwave_log)
wave msdtempwave, msdxwave_log
variable halflength = 5
variable fitlength = 2*halflength+1
variable msdlength = dimsize(msdtempwave,0)
variable i, j
make/o/n=(msdlength - 2*halflength) msdtempwave_diff=0
wave W_coef, msdxwave
duplicate/o/R=(halflength, msdlength-halflength-1) msdxwave msdxwave_diff


for (i=halflength; i<msdlength-halflength; i+=1) 
       CurveFit/Q /W=2/NTHR=0 line  msdtempwave[i-halflength,i+halflength]  /X=msdxwave_log /D 
       msdtempwave_diff[i-halflength] = W_coef[1]
endfor

end


function smoothmsd(msdtempwave, msdxwave_log,startflag)
wave msdtempwave, msdxwave_log
variable startflag


variable smoothrange = 0.6
wavestats/q msdxwave_log
variable msdxmin = V_min
variable msdxmax = V_max

variable msdlength = dimsize(msdtempwave,0)
wave W_coef
variable istart, iend
variable i



//if (startflag ==0)
    findlevel/P/Q/edge=1 msdxwave_log (V_min+smoothrange/2)
    istart = round(V_levelx)
    findlevel/P/Q/edge=1 msdxwave_log (V_max - smoothrange/2)
    iend = round(V_levelx)
    wave msdxwave
    duplicate/o/R=(istart, iend-1) msdxwave msdxwave_diff
    make/o/n=(msdlength)  rangestart=0, rangeend=0
    make/o/n=(iend-istart) msdtempwave_diff, msdtempwave_smooth
    for(i=istart; i<iend;i+=1)
          findlevel/P/Q/edge=1 msdxwave_log (msdxwave_log[i] - smoothrange/2)
          rangestart[i] = round(V_levelx)
           findlevel/P/Q/edge=1 msdxwave_log (msdxwave_log[i] + smoothrange/2)
           rangeend[i] = round(V_levelx)
    endfor
//endif    


for (i=istart; i<iend; i+=1)
      CurveFit/W=2/Q/NTHR=0/K={msdxwave_log[i]} poly_XOffset 3,  msdtempwave[rangestart[i],rangeend[i]] /X=msdxwave_log /D
      msdtempwave_smooth[i-istart]= W_coef[0]
      msdtempwave_diff[i-istart] = W_coef[1]
      
endfor
      
end



function findpara(Rmin, Rmax, slopethresh, msdwave_diff)
variable Rmin, Rmax
wave MSDwave_diff
variable slopethresh   // threshold for MSD slope underwhich would be consider as start of turn over
wave msdxwave_diff, msdxwave



//variable msdvalue = 20                    //for 10 fps videos
variable msdvalue = 10                    //for 100 fps videos
variable gaplength = 100
wave linmaskval
wave linmaskval_intercept_int
wave msd



variable npost = dimsize(MSDwave_diff, 1)
variable msdlength = dimsize(MSDwave_diff, 0)
variable msd_fixedlagtime
make/o/n=(msdlength) msd_diff_temp
make/o/n=(npost)  msd_slope_max=0, msd_ptime=0, msd_slope_para=0, msd_mag_para=0, msd_pdisp = 0,msd_slope_err=0

findvalue/V=(msdvalue)  msdxwave

msd_fixedlagtime = V_value


variable ipost,i
for(ipost = 0; ipost<npost; ipost+=1)
      msd_diff_temp = msdwave_diff[p][ipost]
      wavestats/q/R=(Rmin,Rmax)   msd_diff_temp
//      if(linmaskval[ipost] ==1)
         if ( linmaskval_intercept_int[ipost] !=0 && linmaskval_intercept_int[ipost] !=3)
//           if (  linmaskval_intercept_int[ipost] !=3)            //changed temprorarily
    //       msd_slope_max[ipost] =  V_max

         msd_slope_para[ipost] = msd_diff_temp[msd_fixedlagtime]
         msd_mag_para[ipost] = msd[msd_fixedlagtime][ipost]
         msd_slope_max[ipost] = mean(msd_diff_temp,Rmin,Rmax)
         msd_slope_err[ipost] = V_sdev
//      for(i=Rmin; i<msdlength; i+=1)
//            if (msd_diff_temp[i] <slopethresh && msd_slope_max[ipost]>1)
//                 msd_ptime[ipost] = msdxwave_diff[i]
//            endif
//              if ( msd_diff_temp[i] < slopethresh*msd_slope_max[ipost])
//                   msd_ptime[ipost] = msdxwave_diff[i]
          if(V_max > 1)
          findlevel/Q/edge=2/R=(Rmin, 500) msd_diff_temp, (V_max*slopethresh)
          msd_ptime[ipost] = msdxwave_diff[round(V_levelx)]
         msd_pdisp[ipost ]= msd[round(V_levelx)][ipost]
          endif

             // endif     
     // endfor    
      endif 
endfor

end      



function calcavgmsd()
wave MSD
wave linmaskval
variable ipost
variable npost = dimsize(MSD,1)
variable msdlength = dimsize(MSD,0)
variable ncell=0
make/o/n=(msdlength)  MSDavg=0

for(ipost=0;ipost<npost;ipost+=1)
     if(linmaskval[ipost]==1)
        MSDavg+=log(MSD[p][ipost])
        ncell+=1
     endif
endfor 

MSDavg = 10^(MSDavg/ncell)

end



//find periphery posts
function findperipposts(linmaskval, astr)
wave linmaskval, astr
variable nxpost = astr[3]
variable nypost = astr[4]
variable ipost, jpost, postid
variable linindicator            //indicate even or odd number of line

duplicate/o linmaskval linmaskval_edge

for (ipost=1; ipost<nxpost-1; ipost+=1)
     for( jpost = 1; jpost<nypost-1; jpost+=1)
              postid = ipost*nypost+jpost
              if(linmaskval[postid]==1)
                  linindicator = -1 + 2*mod(jpost,2)
                 if (linmaskval[postid-nypost]*linmaskval[postid+nypost]*linmaskval[postid+1]*linmaskval[postid-1]*linmaskval[postid+linindicator*nypost+1]*linmaskval[postid+linindicator*nypost-1]==0)
                    linmaskval_edge[postid] = -1
                endif
          endif  
     endfor
endfor      

end    

//tried to automize picking highly and lowly bent posts

function disectmsdmean(msdwave, msdxwave, lengthandangle)
wave msdwave, msdxwave, lengthandangle

wave   linmaskval_intercept_int, linmaskval_edge
variable npost = dimsize(msdwave,1)
variable msdlength = dimsize(msdwave,0)
variable ipost
variable nedge =0 , ninte =0
variable mthresh = 5
variable nhb=0, nlb=0

make/o/n=(msdlength)  msdtemp, msd_edge=0, msd_inte=0, msd_hb=0, msd_lb=0

duplicate/o msdxwave msdxwave_log
msdxwave_log = log(msdxwave)

for(ipost =1; ipost<npost; ipost+=1)
     msdtemp = msdwave[p][ipost] //- msdwave[0][ipost] 
     if (linmaskval_edge[ipost] == 1)
        msd_inte = msd_inte+ log(msdtemp)
        ninte+=1
     elseif (linmaskval_edge[ipost] == -1)
        msd_edge = msd_edge + log(msdtemp)
        nedge+=1
     endif   
//     if(lengthandangle[ipost][0] >= mthresh)
        if (linmaskval_intercept_int[ipost] == 6)
         msd_hb = msd_hb+log(msdtemp)
         nhb +=1
//     elseif((lengthandangle[ipost][0] < mthresh) && (lengthandangle[ipost][0]>0))    
       elseif( linmaskval_intercept_int[ipost] ==1)
         msd_lb = msd_lb+log(msdtemp)
         nlb +=1
     endif             

endfor        

msd_inte = msd_inte/ninte
msd_edge = msd_edge/nedge

msd_hb = msd_hb/nhb
msd_lb = msd_lb/nlb




//CurveFit/NTHR=0 line  msd_inte[100,120] /X=MSDxwave_log /D 
CurveFit/NTHR=0 line  msd_hb[100,110] /X=MSDxwave_log /D 

//CurveFit/NTHR=0 line  msd_edge[100,120] /X=msdxwave_log /D 
CurveFit/NTHR=0 line  msd_lb[100,110] /X=msdxwave_log /D 


msd_inte = 10^(msd_inte)

msd_edge = 10^(msd_edge)

msd_hb = 10^(msd_hb)
msd_lb = 10^(msd_lb)


//duplicate/o msd_hb msd_hb_ed
//duplicate/o msd_lb  msd_lb_ed





end



function disectmsd()


variable forcethresh = 5
wave msd_slope_para
wave msd_mag_para
wave lengthandangle

wave linmaskval_intercept_int

wave msd_slope_max, msd_ptime, msd_pdisp


variable npost = dimsize(msd_slope_para,0)
variable ipost

make/o/n=(npost) msd_slope_para_hb=0, msd_slope_para_lb=0, msd_mag_para_hb=0, msd_mag_para_lb=0
make/o/n=(npost) msd_ptime_hb=0, msd_ptime_lb=0
make/o/n=(npost) msd_slope_para_mb=0, msd_mag_para_mb=0
make/o/n=(npost) msd_pdisp_hb=0, msd_pdisp_lb=0
make/o/n=(npost) msd_slope_para_bg=0, msd_mag_para_bg = 0


make/o/n=(4,2) msd_bistats

//make/o/n=(npost) slope_lb_avg=0, slope_hb_avg=0, mag_lb_avg=0, mag_hb_avg=0

variable nhb=0, nlb=0, nbg = 0


for( ipost=0; ipost<npost; ipost+=1)

    if(msd_slope_max[ipost]>0.5 && msd_slope_max[ipost]<2)

//       if( msd_slope_max[ipost]<0.5)
//       if(lengthandangle[ipost][0] > forcethresh)
         if (linmaskval_intercept_int[ipost] == 6)
//           msd_slope_para_hb[ipost] = msd_slope_max[ipost]
//           msd_mag_para_hb[ipost] = log(msd_mag_para[ipost])
          
           
           msd_slope_para_hb[nhb] = msd_slope_max[ipost]
           msd_mag_para_hb[nhb] = log(msd_mag_para[ipost])
           msd_ptime_hb[nhb] = msd_ptime[ipost]
           msd_pdisp_hb[nhb] = msd_pdisp[ipost]
           
           nhb+=1

           

//       elseif(lengthandangle[ipost][0]>0 && lengthandangle[ipost][0]<=forcethresh)
          elseif (linmaskval_intercept_int[ipost] == 1)
//           msd_slope_para_lb[ipost] = msd_slope_max[ipost]
//           msd_mag_para_lb[ipost] = log(msd_mag_para[ipost])       
           msd_slope_para_lb[nlb] = msd_slope_max[ipost]
           msd_mag_para_lb[nlb] = log(msd_mag_para[ipost])              
           msd_ptime_lb[nlb] = msd_ptime[ipost]
           msd_pdisp_lb[nlb] = msd_pdisp[ipost]           
           
           nlb+=1
           
           elseif (linmaskval_intercept_int[ipost] == 7)
           msd_slope_para_mb[ipost] = msd_slope_max[ipost]
           msd_mag_para_mb[ipost] = msd_mag_para[ipost]
           
//         added temporarily     10/31/2018           
           elseif (linmaskval_intercept_int[ipost] == 0)
           msd_slope_para_bg[nbg] = msd_slope_max[ipost]
           msd_mag_para_bg[nbg] = msd_mag_para[ipost]
           nbg+=1
           
           
           
       endif
    
    
    endif

       
endfor


redimension/n=(nhb) msd_slope_para_hb, msd_ptime_hb, msd_mag_para_hb
redimension/n=(nlb) msd_slope_para_lb, msd_ptime_lb, msd_mag_para_lb
redimension/n=(nbg) msd_slope_para_bg, msd_mag_para_bg

duplicate/o msd_ptime_hb keywave

//sort keywave msd_ptime_hb

duplicate/o msd_ptime_lb keywave

//sort keywave msd_ptime_lb

wavestats/q msd_slope_para_hb
msd_bistats[0][0] = V_avg
msd_bistats[0][1] = V_sdev
wavestats/q msd_mag_para_hb
msd_bistats[1][0] = 10^(V_avg)
msd_bistats[1][1] = V_sdev
wavestats/q msd_slope_para_lb
msd_bistats[2][0] = V_avg
msd_bistats[2][1] = V_sdev
wavestats/q msd_mag_para_lb
msd_bistats[3][0] = 10^(V_avg)
msd_bistats[3][1] = V_sdev

msd_mag_para_hb = 10^(msd_mag_para_hb)
msd_mag_para_lb = 10^(msd_mag_para_lb)


Make/N=20/O MSD_slope_para_hb_Hist;DelayUpdate
Histogram/B={0.1,0.1,20}  MSD_slope_para_hb, MSD_slope_para_hb_Hist

Make/N=20/O MSD_slope_para_lb_Hist;DelayUpdate
Histogram/B={0.1,0.1,20}  MSD_slope_para_lb, MSD_slope_para_lb_Hist

Make/N=6/O MSD_ptime_hb_Hist;DelayUpdate
Histogram/B={20,50,6}  MSD_ptime_hb, MSD_ptime_hb_Hist

//Make/N=6/O MSD_ptime_lb_Hist;DelayUpdate                   //m6
//Histogram/B={20,50,6}  MSD_ptime_lb, MSD_ptime_lb_Hist

Make/N=10/O MSD_ptime_lb_Hist;DelayUpdate                  //m10
Histogram/B={10,20,10}  MSD_ptime_lb, MSD_ptime_lb_Hist



end       




function msdseperate()
//variable directionflag         //1 for east/west, 2 for north/south



wave astr
wave msd_slope_max
wave linmaskval_intercept_int
variable npost = dimsize(msd_slope_max, 0)
variable nxpost = astr[3]
variable nypost = astr[4]
variable ipost, i, j 

wave msd_slope_para_lb, msd_slope_para_hb
variable nlb = dimsize(msd_slope_para_lb,0)
variable nhb = dimsize(msd_slope_para_hb, 0)

variable countlb=0, counthb=0
 
make/o/n=(nlb) indexlb=0, indexlb_even = 0
make/o/n=(nhb) indexhb=0 

make/o/n=(nlb)  lbx, lby
make/o/n=(nhb) hbx, hby

wave linmaskX, linmaskY


for (ipost=0; ipost<npost; ipost+=1)

       i = floor(ipost/nypost)
       j = mod(ipost, nypost)
       if ( linmaskval_intercept_int[ipost] == 1) 
               indexlb[countlb] = j
               indexlb_even[countlb] = mod(countlb,2) 
               lbx[countlb] = linmaskX[ipost]
               lby[countlb] = linmaskY[ipost]
               countlb +=1
       elseif ( linmaskval_intercept_int[ipost] == 6) 
                indexhb[counthb] = j
                              
                hbx[counthb] = linmaskX[ipost]
                hby[counthb] = linmaskY[ipost]
                counthb +=1  
       endif
       
       
     
endfor


          

variable nlbhalf = floor(nlb/2)
variable nhbhalf = floor(nhb/2)

variable tvalue, sdev_1, sdev_2

make/o/n=(3,5) sepstats
wavestats/q/R=(0,nlbhalf-1) msd_slope_para_lb

sepstats[0][0] = V_avg
sepstats[0][1] = V_sem
sdev_1 = V_sdev
wavestats/q/R=(nlbhalf,nlb-1) msd_slope_para_lb

sepstats[0][2] = V_avg
sepstats[0][3] = V_sem
sdev_2 = V_sdev

tvalue = t_value_calc( sepstats[0][0] ,sepstats[0][2], sdev_1, sdev_2,nlbhalf, nlb-nlbhalf)
sepstats[0][4] = studentA(tvalue, nlb-2)


duplicate/o/R=(0,nlbhalf-1) msd_slope_para_lb msd_slope_lb_es
duplicate/o/R=(nlbhalf,nlb-1) msd_slope_para_lb msd_slope_lb_we
duplicate/o/R=(0,nlbhalf-1) lbx, lbx_es
duplicate/o/R=(nlbhalf, nlb-1) lbx, lbx_we
duplicate/o/R=(0,nlbhalf-1) lby, lby_es
duplicate/o/R=(nlbhalf, nlb-1) lby, lby_we



duplicate/o msd_slope_para_lb msd_slope_para_lb_sort
duplicate/o lbx lbx_sort
duplicate/o lby lby_sort
sort indexlb msd_slope_para_lb_sort
sort indexlb lbx_sort
sort indexlb lby_sort
wavestats/q/R=(0,nlbhalf-1) msd_slope_para_lb_sort
sepstats[1][0] = V_avg
sepstats[1][1] = V_sem
sdev_1 = V_sdev
wavestats/q/R=(nlbhalf, nlb-1) msd_slope_para_lb_sort
sepstats[1][2] = V_avg
sepstats[1][3] = V_sem
sdev_2 = V_sdev

tvalue = t_value_calc( sepstats[1][0] ,sepstats[1][2], sdev_1, sdev_2,nlbhalf, nlb-nlbhalf)
sepstats[1][4] = studentA(tvalue, nlb-2)

duplicate/o/R=(0,nlbhalf-1) msd_slope_para_lb_sort msd_slope_lb_so
duplicate/o/R=(nlbhalf,nlb-1) msd_slope_para_lb_sort msd_slope_lb_no
duplicate/o/R=(0,nlbhalf-1) lbx_sort, lbx_so
duplicate/o/R=(nlbhalf, nlb-1) lbx_sort, lbx_no
duplicate/o/R=(0,nlbhalf-1) lby_sort, lby_so
duplicate/o/R=(nlbhalf, nlb-1) lby_sort, lby_no


sort indexlb_even msd_slope_para_lb_sort
sort indexlb_even lbx_sort
sort indexlb_even lby_sort
wavestats/q/R=(0,nlbhalf-1) msd_slope_para_lb_sort
sepstats[2][0] = V_avg
sepstats[2][1] = V_sem
sdev_1 = V_sdev
wavestats/q/R=(nlbhalf, nlb-1) msd_slope_para_lb_sort
sepstats[2][2] = V_avg
sepstats[2][3] = V_sem
sdev_2 = V_sdev

tvalue = t_value_calc( sepstats[2][0] ,sepstats[2][2], sdev_1, sdev_2,nlbhalf, nlb-nlbhalf)
sepstats[2][4] = studentA(tvalue, nlb-2)

duplicate/o/R=(0,nlbhalf-1) msd_slope_para_lb_sort msd_slope_lb_even
duplicate/o/R=(nlbhalf,nlb-1) msd_slope_para_lb_sort msd_slope_lb_odd
duplicate/o/R=(0,nlbhalf-1) lbx_sort, lbx_even
duplicate/o/R=(nlbhalf, nlb-1) lbx_sort, lbx_odd
duplicate/o/R=(0,nlbhalf-1) lby_sort, lby_even
duplicate/o/R=(nlbhalf, nlb-1) lby_sort, lby_odd


end



function/c t_value_calc(xavg, yavg, xsdev, ysdev,nx,ny)

variable xavg, yavg, xsdev, ysdev, nx, ny

variable tvalue

tvalue = (xavg-yavg)/(sqrt(((nx-1)*xsdev^2+(ny-1)*ysdev^2)/(nx+ny-2))*sqrt(1/nx+1/ny))

//print(tvalue)

return tvalue

end


function  selectcellposts()

wave zeropos
wave msd_slope_max
wave msd_mag_para

wave linmaskval_intercept_int
wave linmaskval_orig

variable npost = dimsize(zeropos,0)
variable i
variable nc=0

make/o/n=(npost,4)   slopecolordata

for(i =0; i<npost; i+=1)
         slopecolordata[i][0] = zeropos[i][0]
         slopecolordata[i][1] = zeropos[i][1]      
//      if (linmaskval_intercept_int[i] >0)          
         if (linmaskval_orig[i] >0)    
         slopecolordata[i][2] = msd_slope_max[i]
         slopecolordata[i][3] = msd_mag_para[i]
      elseif ( linmaskval_orig[i] ==-1)
         slopecolordata[i][2] = -1
          slopecolordata[i][3] = -1
      elseif (  linmaskval_orig[i]  ==-2)
        slopecolordata[i][2] = -2
         slopecolordata[i][3] = -2
      else
        slopecolordata[i][2] = 0  
         slopecolordata[i][3] = 0
      endif
      
endfor

//redimension/n=(nc,3) slopecolordata

end          


function calcMSD_part()

wave beginmag
wave linmaskval_intercept_int

variable npost = dimsize(linmaskval_intercept_int,0)
variable ipost
variable nlb = 0
variable nhb = 0
variable pixelratio = 8



make/o/n=(2,2) MSD_bistats_seg
make/o/n=(npost) msd_mag_lb_seg, msd_mag_hb_seg

for (ipost =0; ipost < npost; ipost+=1)
      if (linmaskval_intercept_int[ipost] ==1)
          msd_mag_lb_seg[nlb] = log(beginmag[ipost]/pixelratio^2) 
          nlb+=1
      elseif (linmaskval_intercept_int[ipost] == 6)
          msd_mag_hb_seg[nhb] = log(beginmag[ipost]/pixelratio^2)
          nhb+=1
      endif

endfor


redimension/n=(nlb) msd_mag_lb_seg
redimension/n=(nhb) msd_mag_hb_seg


wavestats/q msd_mag_lb_seg
MSD_bistats_seg[0][0] = 10^V_avg
msd_bistats_seg[0][1] = 10^V_sem
wavestats/q msd_mag_hb_seg
MSD_bistats_seg[1][0] = 10^V_avg
msd_bistats_seg[1][1] = 10^V_sem



msd_mag_lb_seg = 10^(msd_mag_lb_seg)
msd_mag_hb_seg = 10^(msd_mag_hb_seg)

end









//*********************************************** Correlation Function***************************************//

function calculatecorrelation(dataw, astrw)
wave dataw
wave astrw

variable npostx = astrw[3]
variable nposty = astrw[4]
variable npost = npostx*nposty
variable ipost, jpost
variable f1 = dimoffset(dataw,2)
variable f2 = dimsize(dataw, 2) + f1
variable iframe
variable diff1x, diff1y, diff2x, diff2y
variable corsumx, corsumy, corsum

wave cortemp
wave cortemp1x, cortemp2x, cortemp1y, cortemp2y
make/N=(npost, npost,2 )/O corrwave
corrwave = 0
make/N=(npost, 2, 2)/O datawstat
make/N=(f2-f1)/O cortemp

for (ipost = 0 ; ipost <npost; ipost+=1)
     cortemp = dataw[ipost][0][p]
     wavestats/Q cortemp
     datawstat[ipost][0][0] = V_avg
     datawstat[ipost][1][0] = V_sdev
     cortemp = dataw[ipost][1][p]
     wavestats/Q cortemp
     datawstat[ipost][0][1] = V_avg
     datawstat[ipost][1][1] = V_sdev

endfor

make/N=(f2-f1)/O cortemp1x, cortemp2x, cortemp1y, cortemp2y

for (ipost = 0 ; ipost <npost; ipost+=1)
      for(jpost = ipost; jpost<npost; jpost+=1)
           cortemp1x = dataw[ipost][0][p]
           cortemp1y = dataw[ipost][1][p]
           cortemp2x = dataw[jpost][0][p]
           cortemp2y = dataw[jpost][1][p]
           corsumx = 0 
           corsumy = 0
           for (iframe = f1; iframe<f2; iframe+=1)
                diff1x = cortemp1x[iframe] - datawstat[ipost][0][0]
                diff1y = cortemp1y[iframe] - datawstat[ipost][0][1]
                diff2x = cortemp2x[iframe] - datawstat[jpost][0][0]
                diff2y = cortemp2y[iframe] - datawstat[jpost][0][1]               
                corsumx+= diff1x*diff2x
                corsumy+= diff1y*diff2y
                
           endfor
           corrwave[ipost][jpost][0] = corsumx/(datawstat[ipost][1][0]*datawstat[jpost][1][0])/(f2-f1)
           corrwave[jpost][ipost][0] = corsumx/(datawstat[ipost][1][0]*datawstat[jpost][1][0])/(f2-f1)
           corrwave[ipost][jpost][1] = corsumy/(datawstat[ipost][1][1]*datawstat[jpost][1][1])/(f2-f1)
           corrwave[jpost][ipost][1] = corsumy/(datawstat[ipost][1][1]*datawstat[jpost][1][1])/(f2-f1)
       endfor
endfor

make/N=(npost, npost)/O corrwaveall

for (ipost = 0 ; ipost<npost; ipost+=1)
     for (jpost = ipost; jpost < npost; jpost+=1)
          corrwaveall[ipost][jpost] = sqrt((corrwave[ipost][jpost][0]^2 + corrwave[ipost][jpost][1]^2)/2)
          corrwaveall[jpost][ipost] = corrwaveall[ipost][jpost]
     endfor
endfor

end


function  corrheatmap(corrwave, maskwave, astrw, targetpost)
wave corrwave
wave maskwave
wave astrw
variable  targetpost

variable npostx = astrw[3]
variable nposty = astrw[4]
variable npost = npostx*nposty

variable ipost, jpost

wave linmasky 
wave linmaskx

//for (ipost = 0; ipost<npost; ipost+=1)
//     corrwave[ipost][ipost] = 1
//endfor

make/N=(npost)/O corrhm


         

corrhm = corrwave[targetpost][p][0]


dowindow/K HeatMapcorr
newimage/N=HeatMapcorr M_rotatedimage
appendtograph/T/W=Heatmapcorr linmasky vs linmaskx
modifygraph mode=3, marker=0
modifygraph zcolor(linmasky)={corrhm,*,*, Rainbow, 0}

end     








function   calccorrelation (dataw, astrw, maskwave)
wave dataw, astrw, maskwave


variable npostx = astrw[3]
variable nposty = astrw[4]
variable npost = npostx*nposty
variable ipost, jpost
variable f1 = dimoffset(dataw,2)
variable f2 = dimsize(dataw, 2) + f1
variable iframe
variable diff1x, diff1y, diff2x, diff2y
wave    corsumx, corsumy, corsum
variable  lagtime
variable  laglength = floor((f2 - f1)/6)
variable  i, j
variable ncell
variable postid
variable icell

//counting number of cell posts and their post numbers
make/N=1000/O idcell
ncell = 0
for(i=0; i<npostx; i+=1)
    for(j=0; j<nposty; j+=1) 
        ipost=i*nposty+j
         if(maskwave[i][j]==kCell || maskwave[i][j] == 6 ) 	
            idcell[ncell]=ipost
            ncell+=1
            endif
       endfor
   endfor
   
   
// ncell = 10

redimension/N=(ncell) idcell
 
wave cortemp1x, cortemp2x, cortemp1y, cortemp2y
wave  corix1, corix2, coriy1, coriy2, corjx1, corjx2, corjy1, corjy2
make/N=(ncell, ncell,(100+(laglength-100)/100),2 )/O corrwave
make/N=(ncell, ncell)/O corrtemp
corrwave = 0
make/O/N=(f2-f1)/O cortemp1x, cortemp2x, cortemp1y, cortemp2y
make/O/N=(laglength-1) corrxwave
make/O/N=(laglength-1) corravg
make/O/N=(ncell, f2, (100+(laglength-100)/100),2) lagwave


//calculate correlation for background post to reduce impact from residue background drifting  9/21/2015
make/O/N=(npost, npost, 2) bgcorrwave
make/O/N=(npost, f2,2)   bglagwave



//used for normalization
variable inormx, jnormx, inormy, jnormy , normx, normy
wave inormwx, jnormwx, inormwy, jnormwy


//calculate background correlation level
for (ipost = 0; ipost<npost; ipost+=1)
     
         cortemp1x = dataw[ipost][0][p]
         cortemp1y = dataw[ipost][1][p]
         lagtime = 100
               duplicate/o/r=[0,f2-lagtime-1] cortemp1x  corix1
               duplicate/o/r=[lagtime, f2-1] cortemp1x  corix2
               duplicate/o/r=[0,f2-lagtime-1] cortemp1y  coriy1
               duplicate/o/r=[lagtime, f2-1] cortemp1y  coriy2
               for (iframe = 0; iframe<f2-lagtime; iframe+=1)         
               bglagwave[ipost][iframe][0] = corix2[iframe]-corix1[iframe]
               bglagwave[ipost][iframe][1] = coriy2[iframe]-coriy1[iframe]
               endfor
endfor

for (ipost = 0; ipost<npost; ipost+=1)
     for (jpost = ipost; jpost<npost; jpost+=1)

               make/N=(f2-lagtime)/O corsumx, corsumy, inormwx, jnormwx, inormwy, jnormwy

                              
               corsumx = bglagwave[ipost][p][0]*bglagwave[jpost][p][0]
               corsumy = bglagwave[ipost][p][1]*bglagwave[jpost][p][1]
               
               //normalization
               inormwx = bglagwave[ipost][p][0]*bglagwave[ipost][p][0]
               jnormwx = bglagwave[jpost][p][0]*bglagwave[jpost][p][0]
               inormwy = bglagwave[ipost][p][1]*bglagwave[ipost][p][1]
               jnormwy = bglagwave[jpost][p][1]*bglagwave[jpost][p][1]
               
               inormx = sqrt( mean(inormwx) )
               jnormx = sqrt( mean(jnormwx) )
               inormy = sqrt( mean(inormwy) )
               jnormy = sqrt( mean(jnormwy) )
               normx = inormx * jnormx
               normy = inormy * jnormy

               
               bgcorrwave[ipost][jpost][0] = mean(corsumx)/normx
               bgcorrwave[jpost][ipost][0] = mean(corsumx)/normx
               bgcorrwave[ipost][jpost][1] = mean(corsumy)/normy
               bgcorrwave[jpost][ipost][1] = mean(corsumy)/normy 
endfor
endfor          
               

for (icell = 0; icell<ncell; icell+=1)  
  
         ipost = idcell[icell]
                
         cortemp1x = dataw[ipost][0][p]
         cortemp1y = dataw[ipost][1][p]
         
    for (lagtime = 1; lagtime <100; lagtime+=1)
               duplicate/o/r=[0,f2-lagtime-1] cortemp1x  corix1
               duplicate/o/r=[lagtime, f2-1] cortemp1x  corix2
               duplicate/o/r=[0,f2-lagtime-1] cortemp1y  coriy1
               duplicate/o/r=[lagtime, f2-1] cortemp1y  coriy2
               for (iframe = 0; iframe<f2-lagtime; iframe+=1)         
               lagwave[icell][iframe][lagtime-1][0] = corix2[iframe]-corix1[iframe]
               lagwave[icell][iframe][lagtime-1][1] = coriy2[iframe]-coriy1[iframe]
               endfor
     endfor
     
     

//making wave for displacement over different lagtime
     
     for(lagtime = 100; lagtime<laglength; lagtime+=100)
               j = 100+ (lagtime-100)/100     
               duplicate/o/r=[0,f2-lagtime-1] cortemp1x  corix1
               duplicate/o/r=[lagtime, f2-1] cortemp1x  corix2
               duplicate/o/r=[0,f2-lagtime-1] cortemp1y  coriy1
               duplicate/o/r=[lagtime, f2-1] cortemp1y  coriy2              
               for (iframe = 0; iframe<f2-lagtime; iframe+=1)         
               lagwave[icell][iframe][j-1][0] = corix2[iframe]-corix1[iframe]
               lagwave[icell][iframe][j-1][1] = coriy2[iframe]-coriy1[iframe]
               endfor
      endfor

endfor
      






for ( ipost= 0; ipost < ncell; ipost+=1)    
     for ( jpost = ipost ; jpost<ncell ; jpost+=1)
 //             ipsot = idcell[icell]
 //             jpost = idcell[jcell]
 
     
          for ( lagtime = 1; lagtime < 100; lagtime+=1)
               make/N=(f2-lagtime)/O corsumx, corsumy, inormwx, jnormwx, inormwy, jnormwy

                              
               corsumx = lagwave[ipost][p][lagtime-1][0]*lagwave[jpost][p][lagtime-1][0]
               corsumy = lagwave[ipost][p][lagtime-1][1]*lagwave[jpost][p][lagtime-1][1]
               
            //  corsumr = 
           //  corsumt =
               
               //normalization
               inormwx = lagwave[ipost][p][lagtime-1][0]*lagwave[ipost][p][lagtime-1][0]
               jnormwx = lagwave[jpost][p][lagtime-1][0]*lagwave[jpost][p][lagtime-1][0]
               inormwy = lagwave[ipost][p][lagtime-1][1]*lagwave[ipost][p][lagtime-1][1]
               jnormwy = lagwave[jpost][p][lagtime-1][1]*lagwave[jpost][p][lagtime-1][1]
               
               inormx = sqrt( mean(inormwx) )
               jnormx = sqrt( mean(jnormwx) )
               inormy = sqrt( mean(inormwy) )
               jnormy = sqrt( mean(jnormwy) )
               normx = inormx * jnormx
               normy = inormy * jnormy

               
               corrwave[ipost][jpost][lagtime-1][0] = mean(corsumx)/normx
               corrwave[jpost][ipost][lagtime-1][0] = mean(corsumx)/normx
               corrwave[ipost][jpost][lagtime-1][1] = mean(corsumy)/normy
               corrwave[jpost][ipost][lagtime-1][1] = mean(corsumy)/normy     
               corrxwave[lagtime-1] = lagtime * 0.01
               
               
         endfor 
         
          for ( lagtime = 100; lagtime < laglength; lagtime+=100)
               j = 100+ (lagtime-100)/100
               make/N=(f2-lagtime)/O corsumx, corsumy, inormwx, jnormwx, inormwy, jnormwy
                              
               corsumx = lagwave[ipost][p][j-1][0]*lagwave[jpost][p][j-1][0]
               corsumy = lagwave[ipost][p][j-1][1]*lagwave[jpost][p][j-1][1]
               
               
//               duplicate/o/r=[0,f2-lagtime-1] cortemp1x  corix1
//               duplicate/o/r=[lagtime, f2-1] cortemp1x  corix2
//               duplicate/o/r=[0,f2-lagtime-1] cortemp2x  corjx1
//               duplicate/o/r=[lagtime, f2-1] cortemp2x  corjx2
//               duplicate/o/r=[0,f2-lagtime-1] cortemp1y  coriy1
//               duplicate/o/r=[lagtime, f2-1] cortemp1y  coriy2
//               duplicate/o/r=[0,f2-lagtime-1] cortemp2y  corjy1
//               duplicate/o/r=[lagtime, f2-1] cortemp2y  corjy2
               
//               corsumx = (corix2[p] - corix1[p])*(corjx2[p] - corjx1[p])
//               corsumy = (coriy2[p] - coriy1[p])*(corjy2[p] - corjy1[p])


               inormwx = lagwave[ipost][p][j-1][0]*lagwave[ipost][p][j-1][0]
               jnormwx = lagwave[jpost][p][j-1][0]*lagwave[jpost][p][j-1][0]
               inormwy = lagwave[ipost][p][j-1][1]*lagwave[ipost][p][j-1][1]
               jnormwy = lagwave[jpost][p][j-1][1]*lagwave[jpost][p][j-1][1]
               
               inormx = sqrt( mean(inormwx) )
               jnormx = sqrt( mean(jnormwx) )
               inormy = sqrt( mean(inormwy) )
               jnormy = sqrt( mean(jnormwy) )
               normx = inormx * jnormx
               normy = inormy * jnormy

               
               corrwave[ipost][jpost][j-1][0] = mean(corsumx)/normx
               corrwave[ipost][jpost][j-1][1] = mean(corsumy)/normy
               corrwave[jpost][ipost][j-1][0] = mean(corsumx)/normx
               corrwave[jpost][ipost][j-1][1] = mean(corsumy)/normy                    
               corrxwave[j-1] = lagtime * 0.01
               
                              
         endfor 
         
endfor
endfor         

redimension/N =(j)   corrxwave, corravg
redimension/N=(-1,-1,j,-1)  corrwave     
 
 
variable nop=ncell*(ncell-1)


 
for ( lagtime = 0; lagtime<j; lagtime+=1) 
    corrtemp = corrwave[p][q][lagtime][0]
    matrixOP/O tracetemp = trace(corrtemp) 
    corravg[lagtime] = (sum(corrtemp) - tracetemp[0])/nop
endfor


setscale z, 0, 0.01,  corrwave 
 
end         
                   
//                     diff1x = cortemp1x[iframe+lagtime] - cortemp1x[iframe] 
//                     diff2x = cortemp2x[iframe+lagtime] - cortemp2x[iframe] 
//                     diff1y = cortemp1y[iframe+lagtime] - cortemp1y[iframe]
//                     diff2y = cortemp2y[iframe+lagtime] - cortemp2y[iframe]
                     
//                     corsumx = corsumx+diff1x*diff2x
//                     corsumy = corsumy+diff1y*diff2y
                     
                     
                     
function  corrheatmap_new(corrwave,bgcorrwave maskwave, astrw, targetpost, lagtime,idcell, bgflag)
wave corrwave
wave bgcorrwave
wave maskwave
wave astrw
wave idcell
variable  targetpost
variable  lagtime
variable bgflag
variable  j 
variable posti, postj



if( lagtime<1)
   j = lagtime*100
else 
   j = 100 + (lagtime*100 - 100)/100
endif


variable npostx = astrw[3]
variable nposty = astrw[4]
variable npost = npostx*nposty

variable ipost, jpost

variable bgcorrV=0, nbgpost=0

wave linmasky 
wave linmaskx

//for (ipost = 0; ipost<npost; ipost+=1)
//     corrwave[ipost][ipost] = 1
//endfor

make/N=(npost)/O corrhm = 0
variable ncell = dimsize(idcell,0)


//
if (bgflag ==0)         
for( ipost = 0; ipost <ncell; ipost+=1)
      corrhm[idcell[ipost]]  = corrwave[targetpost][ipost][j][0]
endfor       
//

else

//when plotting background posts    9/22/2015
for (ipost=0; ipost<npost; ipost+=1)
       posti = floor(ipost/nposty)
       postj = mod(ipost, nposty)
       if (maskwave[posti][postj] == kEmpty)
           bgcorrV += bgcorrwave[targetpost][ipost][0]
           nbgpost+=1
       endif
endfor

bgcorrV = bgcorrV/nbgpost

for (ipost=0; ipost<npost; ipost+=1)
      corrhm[ipost] = bgcorrwave[targetpost][ipost][0]
endfor

endif





dowindow/K HeatMapcorr
newimage/N=HeatMapcorr M_rotatedimage
appendtograph/T/W=Heatmapcorr linmasky vs linmaskx
modifygraph mode=3, marker=0
modifygraph zcolor(linmasky)={corrhm,*,*, Rainbow, 0}

end  









function calcMSDfreq(freq, ipost)
wave freq
variable ipost

variable i 
variable nofreq = dimsize(freq,0)
variable MSDmax
//variable MSDmin
variable MSDrange = 5
string MSDname , fitaffix

make/o/n=(nofreq)  MSDmagfreq


for (i = 0; i<nofreq; i+=1)
          if (freq[i] <1)
             sprintf fitaffix "_0_%dhz"  (freq[i]*10)
          else 
             sprintf fitaffix "_%dhz"  freq[i]
          endif
       MSDname = "MSD"+fitaffix
       
      duplicate/o $MSDname   MSDtemp
      
      MSDmax = dimsize(MSDtemp, 0)
      
      MSdmagfreq[i] = MSDtemp[msdmax][ipost] - MSDtemp[0][ipost]
      
 endfor
 
 
 end     
      
      




function bondcorr(seglength, datawave)
//variable posttype      // 1 for low traction force posts, 6 for high traction force posts
variable seglength    //units in frames   
wave datawave     

wave astr
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost
variable movielength = dimsize(datawave,2)

variable nseg = floor(movielength/seglength)

wave zeropos


wave linmaskval_intercept_int
variable i, j, k, l 
variable postid
make/o/n=3 postpairid
variable postpairid_cor
make/o/n=(npost*3, nseg) disp_corr, velocity_corr,disp_t_corr,disp_cov, disp_nfactor
make/o/n=(npost*3, 2) bondid
make/o/n=(npost*3,2) bonddisp           //displacement between 1st and last frames for all bonds
make/o/n=(npost*3) bond_x, bond_y
make/o/n=(movielength)  pairwave1x, pairwave1y, pairwave2x, pairwave2y
make/o/n=(movielength)  pairwave1p, pairwave1t, pairwave2p, pairwave2t
make/o/n=(seglength) pairwave1p_seg, pairwave2p_seg,pairwave1t_seg, pairwave2t_seg

make/o/n=(movielength,2,npost*3) pairdispwave, pairdispwave_t
make/o/n=(seglength, 2, npost*3, nseg) pairdispwave_seg, pairdispwave_t_seg
 

variable dx, dy, theta

variable diff_p1, diff_p2, diff_t1, diff_t2

variable pairid = 0 

variable alpha, v_covariance, v_nfactor

//removing posts at edge for the moment)
for (i=1; i<nxpost; i+=1)
      for(j=1; j<nypost; j+=1)
            postid = i*nypost+j
            if (linmaskval_intercept_int[postid] ==1 || linmaskval_intercept_int[postid] ==6 || linmaskval_intercept_int[postid] ==7)
                postpairid[0] = (i-1)*nypost+ j
                postpairid[1] = i*nypost+j-1
                postpairid[2] = (i+2*(mod(j,2)-0.5))*nypost + j-1
                for( k = 0; k<3;k+=1)
                       postpairid_cor = postpairid[k] 
                       if (linmaskval_intercept_int[postid] ==1 || linmaskval_intercept_int[postid] ==6 || linmaskval_intercept_int[postid] ==7)
                           bondid[pairid][0] = postid
                           bondid[pairid][1] = postpairid_cor
                           bond_x[pairid] = 0.5*(zeropos[postid][0]+zeropos[postpairid_cor][0])
                           bond_y[pairid] = 0.5*(zeropos[postid][1]+zeropos[postpairid_cor][1])
                           dx = zeropos[postid][0] - zeropos[postpairid_cor][0]
                           dy = zeropos[postid][1] - zeropos[postpairid_cor][1]
                           theta = acos (dx/sqrt(dx^2+dy^2))
                           
                           pairwave1x = datawave[postid][0][p]
                           pairwave2x = datawave[postpairid_cor][0][p] 
                           pairwave1y = datawave[postid][1][p]
                           pairwave2y = datawave[postpairid_cor][1][p]
                           pairwave1p = pairwave1x*cos(theta)+pairwave1y*sin(theta)
                           pairwave1t = -pairwave1x*sin(theta)+pairwave1y*cos(theta)
                           pairwave2p = pairwave2x*cos(theta)+pairwave2y*sin(theta)
                           pairwave2t =  -pairwave2x*sin(theta)+pairwave2y*cos(theta)    
                           
                             
                           pairdispwave[][0][pairid] = pairwave1p[p]                            
                           pairdispwave[][1][pairid] = pairwave2p[p]       
                           pairdispwave_t[][0][pairid] = pairwave1t[p]                            
                           pairdispwave_t[][1][pairid] = pairwave2t[p]   
                           
                           bonddisp[pairid][0] = pairwave1p[movielength-1] - pairwave1p[0]
                           bonddisp[pairid][1] = pairwave2p[movielength-1] - pairwave2p[0]                            
                           
                           for (l=0; l<nseg; l+=1)
                                 pairwave1p_seg = pairwave1p[p+l*seglength]
                                 pairwave2p_seg = pairwave2p[p+l*seglength]
                                 pairwave1t_seg =  pairwave1t[p+l*seglength]
                                 pairwave2t_seg =  pairwave2t[p+l*seglength]
                                 
                                 alpha = statscorrelation(pairwave1p_seg, pairwave2p_seg)
                                 v_covariance = calccov(pairwave1p_seg, pairwave2p_seg)
                                 v_nfactor = calcnfactor(pairwave1p_seg, pairwave2p_seg)
                                 
                                                                
//                                 alpha = altcorrelation(pairwave1p_seg, pairwave2p_seg)
                                 disp_corr[pairid][l] = alpha
                                 disp_cov[pairid][l] = v_covariance
                                 disp_nfactor[pairid][l] = v_nfactor
                                 
                                 
                                 alpha = statscorrelation(pairwave1t_seg, pairwave2t_seg)
//                                 alpha = altcorrelation(pairwave1t_seg, pairwave2t_seg)
                                 disp_t_corr[pairid][l] = alpha
                                 
                            endfor  
                               
                           pairid+=1
                       endif
                 endfor    
               endif
          endfor                 
endfor                           


redimension/n=(pairid, nseg) disp_corr, disp_t_corr,disp_cov, disp_nfactor
redimension/n=(pairid,2) bondid, bonddisp
redimension/n=(movielength, 2, pairid) pairdispwave, pairdispwave_t
redimension/n=(pairid) bond_x, bond_y

dowindow/k bondcorr0
newimage/N=bondcorr0 rotatedorig1
appendtograph/T/W=bondcorr0 bond_y vs bond_x
ModifyGraph mode=3
ModifyGraph zColor(bond_y)={disp_corr[*][0],-1,1,Rainbow,0}
ModifyGraph marker=19,msize=2



end



function bondcorr_seg(posttype, seglength, datawave)
variable posttype      // 1 for low traction force posts, 6 for high traction force posts
variable seglength    //units in frames   
wave datawave     

wave astr
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost
variable movielength = dimsize(datawave,2)

variable nseg = floor(movielength/seglength)

wave zeropos


wave linmaskval_seg_int
variable i, j, k, l 
variable postid
make/o/n=3 postpairid
variable postpairid_cor
make/o/n=(npost*3*nseg,1) disp_corr_seg, disp_t_corr_seg//,disp_cov, disp_nfactor
make/o/n=(npost*3*nseg, 3) bondid
make/o/n=(npost*3*nseg,2) bonddisp           //displacement between 1st and last frames for all bonds
make/o/n=(npost*3) bond_x, bond_y
make/o/n=(seglength)  pairwave1x, pairwave1y, pairwave2x, pairwave2y
make/o/n=(seglength)  pairwave1p, pairwave1t, pairwave2p, pairwave2t


//for extended waves

make/o/n=(seglength*2)  pairwave1x_ext, pairwave1y_ext, pairwave2x_ext, pairwave2y_ext
make/o/n=(seglength*2)  pairwave1p_ext, pairwave1t_ext, pairwave2p_ext, pairwave2t_ext


make/o/n=(seglength,2,npost*3) pairdispwave_seg, pairdispwave_seg_t
make/o/n=(2*seglength,2,npost*3) pairdispwave_seg_ext
 

variable dx, dy, theta

variable diff_p1, diff_p2, diff_t1, diff_t2

variable pairid = 0 

variable alpha, v_covariance, v_nfactor

//removing posts at edge for the moment)
for (i=1; i<nxpost; i+=1)
      for(j=1; j<nypost; j+=1)
            postid = i*nypost+j
            for (l=0; l<nseg; l+=1)
              if (linmaskval_seg_int[postid][l] == posttype)
                  postpairid[0] = (i-1)*nypost+ j
                  postpairid[1] = i*nypost+j-1
                  postpairid[2] = (i+2*(mod(j,2)-0.5))*nypost + j-1
                for( k = 0; k<3;k+=1)
                       postpairid_cor = postpairid[k] 
                       if (linmaskval_seg_int[postpairid_cor][l] == posttype)
                           bondid[pairid][0] = postid
                           bondid[pairid][1] = postpairid_cor
                           bondid[pairid][2] = l
                           bond_x[pairid] = 0.5*(zeropos[postid][0]+zeropos[postpairid_cor][0])
                           bond_y[pairid] = 0.5*(zeropos[postid][1]+zeropos[postpairid_cor][1])
                           dx = zeropos[postid][0] - zeropos[postpairid_cor][0]
                           dy = zeropos[postid][1] - zeropos[postpairid_cor][1]
                           theta = acos (dx/sqrt(dx^2+dy^2))
                           
                           pairwave1x = datawave[postid][0][p+l*seglength]
                           pairwave2x = datawave[postpairid_cor][0][p+l*seglength] 
                           pairwave1y = datawave[postid][1][p+l*seglength]
                           pairwave2y = datawave[postpairid_cor][1][p+l*seglength]
                           pairwave1p = pairwave1x*cos(theta)+pairwave1y*sin(theta)
                           pairwave1t = -pairwave1x*sin(theta)+pairwave1y*cos(theta)
                           pairwave2p = pairwave2x*cos(theta)+pairwave2y*sin(theta)
                           pairwave2t =  -pairwave2x*sin(theta)+pairwave2y*cos(theta)    
                           
                             
                           pairdispwave_seg[][0][pairid] = pairwave1p[p]                            
                           pairdispwave_seg[][1][pairid] = pairwave2p[p]     
                           pairdispwave_seg_t[][0][pairid] = pairwave1t[p]                            
                           pairdispwave_seg_t[][1][pairid] = pairwave2t[p]   
                                                    

                          alpha = statscorrelation(pairwave1p, pairwave2p)
 //                          v_covariance = calccov(pairwave1p_seg, pairwave2p_seg)
 //                          v_nfactor = calcnfactor(pairwave1p_seg, pairwave2p_seg)
                                
                                                                
//                                 alpha = altcorrelation(pairwave1p_seg, pairwave2p_seg)
                                 disp_corr_seg[pairid][0] = alpha
//                                 disp_cov[pairid][l] = v_covariance
//                                 disp_nfactor[pairid][l] = v_nfactor
                                 
                                 
                                 alpha = statscorrelation(pairwave1t, pairwave2t)
//                                 alpha = altcorrelation(pairwave1t_seg, pairwave2t_seg)
                                 disp_t_corr_seg[pairid][0] = alpha
                                 
                                 
//for extending those pairs at the edge                           
                           pairwave1x_ext = datawave[postid][0][p+(l-0.5)*seglength]
                           pairwave2x_ext = datawave[postpairid_cor][0][p+(l-0.5)*seglength] 
                           pairwave1y_ext = datawave[postid][1][p+(l-0.5)*seglength]
                           pairwave2y_ext = datawave[postpairid_cor][1][p+(l-0.5)*seglength]
                           pairwave1p_ext = pairwave1x_ext*cos(theta)+pairwave1y_ext*sin(theta)
                           pairwave1t_ext = -pairwave1x_ext*sin(theta)+pairwave1y_ext*cos(theta)
                           pairwave2p_ext = pairwave2x_ext*cos(theta)+pairwave2y_ext*sin(theta)
                           pairwave2t_ext =  -pairwave2x_ext*sin(theta)+pairwave2y_ext*cos(theta)  
                                                            
                           pairdispwave_seg_ext[][0][pairid] = pairwave1p_ext[p]                            
                           pairdispwave_seg_ext[][1][pairid] = pairwave2p_ext[p]                                    
                                 
                                                                
                           pairid+=1
                       endif
                 endfor    
               endif
              endfor 
          endfor                 
endfor                           


redimension/n=(pairid,1) disp_corr_seg, disp_t_corr_seg//,disp_cov, disp_nfactor
redimension/n=(pairid,3) bondid//, bonddisp
redimension/n=(seglength, 2, pairid) pairdispwave_seg, pairdispwave_seg_t
redimension/n=(seglength*2, 2, pairid) pairdispwave_seg_ext
redimension/n=(pairid) bond_x, bond_y



end





function bondcorr_2ndneighbor(posttype, seglength, datawave)
variable posttype      // 1 for low traction force posts, 6 for high traction force posts
variable seglength    //units in frames   
wave datawave     



wave astr
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost
variable movielength = dimsize(datawave,2)

variable nseg = floor(movielength/seglength)

wave zeropos


wave linmaskval_intercept_int
variable i, j, k, l 
variable postid
make/o/n=3 postpairid
variable postpairid_cor
make/o/n=(npost*3, nseg) disp_corr_2nd, velocity_corr, disp_cov_2nd, disp_nfactor_2nd
make/o/n=(npost*3) bond_x, bond_y
make/o/n=(npost*3, 2) bonddisp_2nd
make/o/n=(movielength)  pairwave1x, pairwave1y, pairwave2x, pairwave2y
make/o/n=(movielength)  pairwave1p, pairwave1t, pairwave2p, pairwave2t
make/o/n=(seglength) pairwave1p_seg, pairwave2p_seg

make/o/n=(movielength,2,npost*3) pairdispwave
 

variable dx, dy, theta

variable pairid = 0 

variable alpha, v_covariance, v_nfactor

//removing posts at edge for the moment)
for (i=2; i<nxpost; i+=1)
      for(j=1; j<nypost; j+=1)
            postid = i*nypost+j
            if (linmaskval_intercept_int[postid] == posttype)
                postpairid[0] = (i+(mod(j,2)-2))*nypost+ j +1
                postpairid[1] = i*nypost+j-2
                postpairid[2] = (i+(mod(j,2)-2))*nypost + j - 1
                for( k = 0; k<3;k+=1)
                       postpairid_cor = postpairid[k] 
                       if (linmaskval_intercept_int[postpairid_cor] == posttype)
                           bond_x[pairid] = 0.5*(zeropos[postid][0]+zeropos[postpairid_cor][0])
                           bond_y[pairid] = 0.5*(zeropos[postid][1]+zeropos[postpairid_cor][1])
                           dx = zeropos[postid][0] - zeropos[postpairid_cor][0]
                           dy = zeropos[postid][1] - zeropos[postpairid_cor][1]
                           theta = acos (dx/sqrt(dx^2+dy^2))
                           
                           pairwave1x = datawave[postid][0][p]
                           pairwave2x = datawave[postpairid_cor][0][p] 
                           pairwave1y = datawave[postid][1][p]
                           pairwave2y = datawave[postpairid_cor][1][p]
                           pairwave1p = pairwave1x*cos(theta)+pairwave1y*sin(theta)
                           pairwave1t = -pairwave1x*sin(theta)+pairwave1y*cos(theta)
                           pairwave2p = pairwave2x*cos(theta)+pairwave2y*sin(theta)
                           pairwave2t =  -pairwave2x*sin(theta)+pairwave2y*cos(theta)    
                              
                           pairdispwave[][0][pairid] = pairwave1p[p]                            
                           pairdispwave[][1][pairid] = pairwave2p[p]       

                           bonddisp_2nd[pairid][0] = pairwave1p[movielength-1] - pairwave1p[0]
                           bonddisp_2nd[pairid][1] = pairwave2p[movielength-1] - pairwave2p[0]    
                           
                           
                           for (l=0; l<nseg; l+=1)
                                 pairwave1p_seg = pairwave1p[p+l*seglength]
                                 pairwave2p_seg = pairwave2p[p+l*seglength]
                                 
                                 

                                 
                                 
                                 alpha = statscorrelation(pairwave1p_seg, pairwave2p_seg)
                                 v_covariance = calccov(pairwave1p_seg, pairwave2p_seg)
                                 v_nfactor = calcnfactor(pairwave1p_seg, pairwave2p_seg)
                                                                  
//                                 alpha = altcorrelation(pairwave1p_seg, pairwave2p_seg)
                                 disp_corr_2nd[pairid][l] = alpha
                                 disp_cov_2nd[pairid][l] = v_covariance
                                 disp_nfactor_2nd[pairid][l] = v_nfactor
                          
//                                 velocity_corr[pairid][l] = velo_corr(pairwave1p_seg, pairwave2p_seg)

                                 
                                 
                            endfor  
                               
                           pairid+=1
                       endif
                 endfor    
               endif
          endfor                 
endfor                           


redimension/n=(pairid, nseg) disp_corr_2nd, disp_cov_2nd, disp_nfactor_2nd
redimension/n=(movielength, 2, pairid) pairdispwave
redimension/n=(pairid) bond_x, bond_y

end


function bondcorr_ran(posttype, seglength, datawave)
variable posttype      // 1 for low traction force posts, 6 for high traction force posts
variable seglength    //units in frames   
wave datawave     

wave astr
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost
variable movielength = dimsize(datawave,2)

variable nseg = floor(movielength/seglength)
wave disp_corr
variable numbond = dimsize(disp_corr,0)

wave zeropos


wave linmaskval_intercept_int
variable i, j, k, l 
variable postid
make/o/n=3 postpairid
variable postpairid_cor
make/o/n=(npost*3, nseg) disp_corr_ran, velocity_corr, disp_cov_ran, disp_nfactor_ran
make/o/n=(npost*3) bond_x, bond_y
make/o/n=(npost*3, 2) bonddisp_ran
make/o/n=(movielength)  pairwave1x, pairwave1y, pairwave2x, pairwave2y
make/o/n=(movielength)  pairwave1p, pairwave1t, pairwave2p, pairwave2t
make/o/n=(seglength) pairwave1p_seg, pairwave2p_seg,pairwave1t_seg, pairwave2t_seg

make/o/n=(movielength,2,npost*3) pairdispwave, pairdispwave_t

make/o/n=(npost) posttypeid
variable postcount=0

for( i =0; i<npost;i+=1)
      if (linmaskval_intercept_int[i] == posttype)
          posttypeid[postcount] = i
          postcount+=1
      endif
endfor

redimension/n=(postcount) posttypeid          

 

variable dx, dy, theta

variable pairid = 0 

variable alpha, v_covariance, v_nfactor

//removing posts at edge for the moment)
for (i=0; i<postcount; i+=1)
           postid = posttypeid[i]
           for(k=0;k<3;k+=1)
                          j=abs(floor(enoise(postcount-1)))
                          if(j!=i)
                          postpairid_cor = posttypeid[j]
       
                           bond_x[pairid] = 0.5*(zeropos[postid][0]+zeropos[postpairid_cor][0])
                           bond_y[pairid] = 0.5*(zeropos[postid][1]+zeropos[postpairid_cor][1])
                           dx = zeropos[postid][0] - zeropos[postpairid_cor][0]
                           dy = zeropos[postid][1] - zeropos[postpairid_cor][1]
                           theta = acos (dx/sqrt(dx^2+dy^2))
                           
                           pairwave1x = datawave[postid][0][p]
                           pairwave2x = datawave[postpairid_cor][0][p] 
                           pairwave1y = datawave[postid][1][p]
                           pairwave2y = datawave[postpairid_cor][1][p]
                           pairwave1p = pairwave1x*cos(theta)+pairwave1y*sin(theta)
                           pairwave1t = -pairwave1x*sin(theta)+pairwave1y*cos(theta)
                           pairwave2p = pairwave2x*cos(theta)+pairwave2y*sin(theta)
                           pairwave2t =  -pairwave2x*sin(theta)+pairwave2y*cos(theta)    
                              
                           pairdispwave[][0][pairid] = pairwave1p[p]                            
                           pairdispwave[][1][pairid] = pairwave2p[p]       
                           pairdispwave_t[][0][pairid] = pairwave1t[p]                            
                           pairdispwave_t[][1][pairid] = pairwave2t[p]        
                           
                           bonddisp_ran[pairid][0] = pairwave1p[movielength-1] - pairwave1p[0]
                           bonddisp_ran[pairid][1] = pairwave2p[movielength-1] - pairwave2p[0]                                                   
                           
                           for (l=0; l<nseg; l+=1)
                                 pairwave1p_seg = pairwave1p[p+l*seglength]
                                 pairwave2p_seg = pairwave2p[p+l*seglength]
                                 pairwave1t_seg =  pairwave2t[p+l*seglength]
                                 pairwave2t_seg =  pairwave2t[p+l*seglength]
                                 
                                 alpha = statscorrelation(pairwave1p_seg, pairwave2p_seg)
                                 v_covariance = calccov(pairwave1p_seg, pairwave2p_seg)
                                 v_nfactor = calcnfactor(pairwave1p_seg, pairwave2p_seg)      
                                                            
//                                 alpha = altcorrelation(pairwave1p_seg, pairwave2p_seg)
                                 disp_corr_ran[pairid][l] = alpha
                                 disp_cov_ran[pairid][l] = v_covariance
                                 disp_nfactor_ran[pairid][l] = v_nfactor                                 
//                                 alpha = statscorrelation(pairwave1t_seg, pairwave2t_seg)
//                                 disp_corr_ran[pairid][l] = alpha
 //                                velocity_corr[pairid][l] = velo_corr(pairwave1p_seg, pairwave2p_seg)

                                 
                                 
                            endfor  
                               
                           pairid+=1
                           endif
               endfor
               
endfor                           


redimension/n=(numbond, nseg) disp_corr_ran, disp_cov_ran ,disp_nfactor_ran
redimension/n=(movielength, 2, numbond) pairdispwave
redimension/n=(numbond, 2) bonddisp_ran
redimension/n=(pairid) bond_x, bond_y

end



function bondcorr_seg_ran(posttype, seglength, datawave)
variable posttype      // 1 for low traction force posts, 6 for high traction force posts
variable seglength    //units in frames   
wave datawave     
wave disp_corr_seg

wave astr
variable n_pair_gen = dimsize(disp_corr_seg,0)
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost
variable movielength = dimsize(datawave,2)

variable nseg = floor(movielength/seglength)

wave zeropos


wave linmaskval_seg_int
variable i, j, k, l 
variable postid
make/o/n=3 postpairid
variable postpairid_cor
make/o/n=(npost*3*nseg,1) disp_corr_ran, disp_t_corr_ran//,disp_cov, disp_nfactor
make/o/n=(npost*3*nseg, 3) bondid_ran
make/o/n=(npost*3*nseg,2) bonddisp           //displacement between 1st and last frames for all bonds
make/o/n=(npost*3) bond_x, bond_y
make/o/n=(seglength)  pairwave1x, pairwave1y, pairwave2x, pairwave2y
make/o/n=(seglength)  pairwave1p, pairwave1t, pairwave2p, pairwave2t


//for extended waves

make/o/n=(seglength*2)  pairwave1x_ext, pairwave1y_ext, pairwave2x_ext, pairwave2y_ext
make/o/n=(seglength*2)  pairwave1p_ext, pairwave1t_ext, pairwave2p_ext, pairwave2t_ext


make/o/n=(seglength,2,npost*3) pairdispwave_seg, pairdispwave_seg_t
make/o/n=(2*seglength,2,npost*3) pairdispwave_seg_ext
 

variable dx, dy, theta

variable diff_p1, diff_p2, diff_t1, diff_t2

variable pairid = 0 

variable alpha, v_covariance, v_nfactor

//removing posts at edge for the moment)
for (i=0; i<n_pair_gen; i+=1)
	do 
		postid =abs(floor(enoise(npost-1)))
		postpairid_cor = abs(floor(enoise(npost-1)))
//		k = abs(floor(enoise(nseg-1)))
		l = abs(floor(enoise(nseg-1)))
		k=l
		dx = zeropos[postid][0] - zeropos[postpairid_cor][0]
                dy = zeropos[postid][1] - zeropos[postpairid_cor][1]
	while((linmaskval_seg_int[postpairid_cor][l] != posttype) || (linmaskval_seg_int[postid][k] != posttype) || (sqrt(dx^2+dy^2)<100))
	
                           bondid_ran[pairid][0] = postid
                           bondid_ran[pairid][1] = postpairid_cor
                           bondid_ran[pairid][2] = l
                           bond_x[pairid] = 0.5*(zeropos[postid][0]+zeropos[postpairid_cor][0])
                           bond_y[pairid] = 0.5*(zeropos[postid][1]+zeropos[postpairid_cor][1])

                           theta = acos (dx/sqrt(dx^2+dy^2))
                           
                           pairwave1x = datawave[postid][0][p+k*seglength]
                           pairwave2x = datawave[postpairid_cor][0][p+l*seglength] 
                           pairwave1y = datawave[postid][1][p+k*seglength]
                           pairwave2y = datawave[postpairid_cor][1][p+l*seglength]
                           pairwave1p = pairwave1x*cos(theta)+pairwave1y*sin(theta)
                           pairwave1t = -pairwave1x*sin(theta)+pairwave1y*cos(theta)
                           pairwave2p = pairwave2x*cos(theta)+pairwave2y*sin(theta)
                           pairwave2t =  -pairwave2x*sin(theta)+pairwave2y*cos(theta)    
                           
                             
                           pairdispwave_seg[][0][i] = pairwave1p[p]                            
                           pairdispwave_seg[][1][i] = pairwave2p[p]     
                           pairdispwave_seg_t[][0][i] = pairwave1t[p]                            
                           pairdispwave_seg_t[][1][i] = pairwave2t[p]   
                                                    

                          alpha = statscorrelation(pairwave1p, pairwave2p)
 //                          v_covariance = calccov(pairwave1p_seg, pairwave2p_seg)
 //                          v_nfactor = calcnfactor(pairwave1p_seg, pairwave2p_seg)
                                
                                                                
//                                 alpha = altcorrelation(pairwave1p_seg, pairwave2p_seg)
                                 disp_corr_ran[pairid][0] = alpha
//                                 disp_cov[pairid][l] = v_covariance
//                                 disp_nfactor[pairid][l] = v_nfactor
                                 
                                 
                                 alpha = statscorrelation(pairwave1t, pairwave2t)
//                                 alpha = altcorrelation(pairwave1t_seg, pairwave2t_seg)
                                 disp_t_corr_ran[pairid][0] = alpha
                                 
                                 
//for extending those pairs at the edge                           
                           pairwave1x_ext = datawave[postid][0][p+(k-0.5)*seglength]
                           pairwave2x_ext = datawave[postpairid_cor][0][p+(l-0.5)*seglength] 
                           pairwave1y_ext = datawave[postid][1][p+(k-0.5)*seglength]
                           pairwave2y_ext = datawave[postpairid_cor][1][p+(l-0.5)*seglength]
                           pairwave1p_ext = pairwave1x_ext*cos(theta)+pairwave1y_ext*sin(theta)
                           pairwave1t_ext = -pairwave1x_ext*sin(theta)+pairwave1y_ext*cos(theta)
                           pairwave2p_ext = pairwave2x_ext*cos(theta)+pairwave2y_ext*sin(theta)
                           pairwave2t_ext =  -pairwave2x_ext*sin(theta)+pairwave2y_ext*cos(theta)  
                                                            
                           pairdispwave_seg_ext[][0][pairid] = pairwave1p_ext[p]                            
                           pairdispwave_seg_ext[][1][pairid] = pairwave2p_ext[p]                                    
                                 
                                                                
                           pairid+=1
               
endfor                           


redimension/n=(pairid,1) disp_corr_ran, disp_t_corr_ran//,disp_cov, disp_nfactor
redimension/n=(pairid,3) bondid_ran//, bonddisp
redimension/n=(seglength, 2, pairid) pairdispwave_seg, pairdispwave_seg_t
redimension/n=(seglength*2, 2, pairid) pairdispwave_seg_ext
redimension/n=(pairid) bond_x, bond_y



end





function bondcorr_other(posttype, seglength, datawave)
variable posttype      // 1 for low traction force posts, 6 for high traction force posts
variable seglength    //units in frames   
wave datawave     


wave astr
variable postdis = astr[5]


variable movielength = dimsize(datawave,2)

variable nseg = floor(movielength/seglength)

wave zeropos
wave posttypeid

variable npost = dimsize(posttypeid,0)
variable ipost, jpost

variable i, j, k, l 
variable postid
make/o/n=3 postpairid
variable postpairid_cor
make/o/n=(npost^2, nseg) disp_corr_other,disp_t_corr_other
make/o/n=(npost^2, 2) bonddisp_other
make/o/n=(movielength)  pairwave1x, pairwave1y, pairwave2x, pairwave2y
make/o/n=(movielength)  pairwave1p, pairwave1t, pairwave2p, pairwave2t
make/o/n=(seglength) pairwave1p_seg, pairwave2p_seg,pairwave1t_seg, pairwave2t_seg

make/o/n=(movielength,2,npost*3) pairdispwave, pairdispwave_t
 

variable dx, dy, theta

variable pairid = 0 

variable alpha

//removing posts at edge for the moment)
for (ipost=0; ipost<npost; ipost+=1)
      i = posttypeid[ipost]
      for(jpost=ipost+1; jpost<npost; jpost+=1)
      j = posttypeid[jpost]
                           dx = zeropos[i][0] - zeropos[j][0]
                           dy = zeropos[i][1] - zeropos[j][1]
                     if (sqrt(dx^2+dy^2)>3*postdis)      
                           theta = acos (dx/sqrt(dx^2+dy^2))
                           
                           pairwave1x = datawave[i][0][p]
                           pairwave2x = datawave[j][0][p] 
                           pairwave1y = datawave[i][1][p]
                           pairwave2y = datawave[j][1][p]
                           pairwave1p = pairwave1x*cos(theta)+pairwave1y*sin(theta)
                           pairwave1t = -pairwave1x*sin(theta)+pairwave1y*cos(theta)
                           pairwave2p = pairwave2x*cos(theta)+pairwave2y*sin(theta)
                           pairwave2t =  -pairwave2x*sin(theta)+pairwave2y*cos(theta)     
                           
                                       
                           bonddisp_other[pairid][0] = pairwave1p[movielength-1] - pairwave1p[0]
                           bonddisp_other[pairid][1] = pairwave2p[movielength-1] - pairwave2p[0]                               
                           
                                          
                           
                           for (l=0; l<nseg; l+=1)
                                 pairwave1p_seg = pairwave1p[p+l*seglength]
                                 pairwave2p_seg = pairwave2p[p+l*seglength]
                                 pairwave1t_seg =  pairwave1t[p+l*seglength]
                                 pairwave2t_seg =  pairwave2t[p+l*seglength]
                                 
                                 alpha = statscorrelation(pairwave1p_seg, pairwave2p_seg)
//                                 alpha = altcorrelation(pairwave1p_seg, pairwave2p_seg)
                                 disp_corr_other[pairid][l] = alpha
                                 alpha = statscorrelation(pairwave1t_seg, pairwave2t_seg)
                                 disp_t_corr_other[pairid][l] = alpha
                                 
                            endfor  
                               
                           pairid+=1
                       endif
          endfor                 
endfor                           


redimension/n=(pairid, nseg) disp_corr_other, disp_t_corr_other
redimension/n=(pairid,2) bondid, bonddisp_other
redimension/n=(movielength, 2, pairid) pairdispwave, pairdispwave_t
redimension/n=(pairid) bond_x, bond_y



end







function velo_corr(wave1, wave2)

wave wave1, wave2

variable boxcarsize =2000

duplicate/o wave1 wave1_smooth, wave2_smooth

boxsmooth boxcarsize, wave1, wave1_smooth
boxsmooth boxcarsize, wave2, wave2_smooth

variable movielength = dimsize(wave1_smooth, 0)
variable i 

variable j=0

wave W_coef


make/o/n=(movielength-boxcarsize) wave1_diff, wave2_diff


for (i=boxcarsize/2; i<movielength-boxcarsize/2; i+=boxcarsize/10)
     curvefit/w=2/q=1 line wave1_smooth[i-boxcarsize/2, i+boxcarsize/2]
     wave1_diff[j] = w_coef[1]
     curvefit/w=2/q=1 line wave2_smooth[i-boxcarsize/2, i+boxcarsize/2] 
     wave2_diff[j] = w_coef[1]    
     j+=1   
endfor

redimension/n=(j) wave1_diff, wave2_diff

duplicate/o wave1_diff productwave_diff, wave1_sq_diff, wave2_sq_diff

productwave_diff =wave1_diff*wave2_diff
wave1_sq_diff = wave1_diff*wave1_diff
wave2_sq_diff = wave2_diff*wave2_diff

return sum(productwave_diff)/sqrt(sum(wave1_sq_diff)*sum(wave2_sq_diff))

end


function velo_pair(datawave,boxcarsize)
wave datawave
variable boxcarsize           //for short movies



variable npair = dimsize(datawave, 2)
variable videolength = dimsize(datawave, 0)

//variable boxcarsize = 99

variable i

make/o/n=(videolength) tempdisp_1, tempdisp_2, tempdisp_1_smt, tempdisp_2_smt

make/o/n=(videolength-boxcarsize-1,2,npair) vpair

for( i =0; i<npair; i+=1)

       tempdisp_1 = datawave[p][0][i]
       tempdisp_2 = datawave[p][1][i]
       
       boxsmooth boxcarsize, tempdisp_1, tempdisp_1_smt
       boxsmooth boxcarsize, tempdisp_2, tempdisp_2_smt
       
       vpair[][0][i] = tempdisp_1_smt[p+boxcarsize+1] - tempdisp_1_smt[p]
       vpair[][1][i] = tempdisp_2_smt[p+boxcarsize+1] - tempdisp_2_smt[p]
       
endfor




//vpair = datawave[p+boxcarsize+1][q][r] - datawave[p][q][r]       


       
       
end


function corrpeak(v_pair, peakdifthresh, boxcarsize)
wave v_pair
//variable peakdifthresh = 50             
variable peakdifthresh        // for short movies
variable boxcarsize       //same as velo_pair function


variable npair = dimsize(v_pair,2)
variable movielength = dimsize(v_pair, 0)
//variable peakthresh = 0.008
variable peakthresh = 0.004

variable i,j,k,l,m
variable Rend = movielength
variable Rstart = 0 
variable npeak = 50

variable sig1=0, sig2=0

variable counter = 0 

variable peakmax_loc1, peakmax_loc2

make/o/n=(npair) steppair=0, steploc1=0, steploc2=0

make/o/n=(movielength) velo_temp1, velo_temp2
make/o/n=(npeak)  peakpos_1, peakpos_2
make/o/n=(npair) peakmaxdif


wave peakposmax, peakpos_alt, peakpos

l=0
m=0
for (i=0; i<npair; i+=1)
      velo_temp1 = v_pair[p][0][i]
      velo_temp2 = v_pair[p][1][i]
      
      peakpos_1 = 0
      peakpos_2 = 0
      
      findmultipeak(velo_temp1,peakthresh, peakdifthresh)
//    duplicate/o peakpos peakpos_1
      duplicate/o peakpos_alt peakpos_1
      peakmax_loc1 = peakposmax[0]
      findmultipeak(velo_temp2,peakthresh, peakdifthresh)
//      duplicate/o peakpos peakpos_2
      duplicate/o peakpos_alt peakpos_2    
      peakmax_loc2 = peakposmax[0]  

      
      if(peakmax_loc1*peakmax_loc2<0)
         peakmaxdif[m] = peakmax_loc1+peakmax_loc2
         m+=1
      endif    
      
      
      
      for ( j=0; j<npeak;j+=1)
          if ( peakpos_1[j] !=0)
              for ( k=0; k<npeak; k+=1)
                   if ( (abs(peakpos_1[j]+peakpos_2[k])<peakdifthresh)&&(peakpos_1[j]*peakpos_2[k])<0)
                       steppair[l] = i
                       steploc1[l] = peakpos_1[j]+sign(peakpos_1[j])*floor(boxcarsize/2)                       
                       steploc2[l] = peakpos_2[k]+sign(peakpos_2[k])*floor(boxcarsize/2)
                       l+=1                                 //might result in curves with multiple pairs get repeated
                   endif
               endfor
          endif
      endfor             
                   
                    
      
      
      
      
endfor




redimension/n=(l)   steppair, steploc1, steploc2
redimension/n=(m) peakmaxdif

end





function matchseg()

variable seglength =200

wave steppair, steploc1, steploc2, selectpair_seg, pairdispwave, selectpair_segnumber

variable npair = dimsize(selectpair_seg,0)
variable nstep = dimsize(steppair,0)
variable movielength = dimsize(pairdispwave,0)


variable i,j
make/o/n=(seglength*2, 2, nstep) pairstepdisp
make/o/n=(seglength*2, 2 ,nstep) pairstepdisp_2nddif

make/o/n=(nstep) peakdif_pair

make/o/n=(movielength)  p1_disp, p2_disp

variable segstart1, segstart2
variable pairnum
variable segnum1, segnum2
variable segstart_2nddif1, segstart_2nddif2


//generate scaled curve
make/o/n=(seglength*3) p1_scaled, p2_scaled, p1_scaledx, p2_scaledx
make/o/n=(seglength*3, 2, nstep) pairstepscaled, pairstepscaled_x
make/o/n=(nstep,2) V_tau, V_height, V_terr
variable l
variable TC = 100


string  NBname, p1str, p2str
 
            sprintf NBname, "postpairsteps"
            DoWindow/K $NBname
            NewNotebook/K=1/F=1/N=$NBname


variable count = 0
wave w_coef


//parameters for fitting
make/o/n=4 w_coef1
variable baseline, stepheight, stepwidth, stepcent


for(i=0; i<npair; i+=1)
    for(j=0; j<nstep;j+=1)
         if (steppair[j] == selectpair_seg[i]) 
             pairnum = selectpair_seg[i]
             segnum1 = floor(abs(steploc1[j]/900))
             segnum2 = floor(abs(steploc2[j]/900))
             if (segnum1 == selectpair_segnumber[i]|| segnum2==selectpair_segnumber[i])   
             p1_disp = pairdispwave[p][0][pairnum]
             p2_disp = pairdispwave[p][1][pairnum]
             
             
             calc2nddiff(p1_disp, seglength/4)
             duplicate/o dwave_2nddiff p1_2nddiff
             segstart_2nddif1 = find2nddif_peak(p1_2nddiff, steploc1[j],seglength/2, seglength/4) - seglength/2            //align at the beginning of the step
//             segstart_2nddif1 = find2nddif_peak(p1_2nddiff, -steploc1[j],seglength/2, seglength/4) - 3*seglength/2       //align at the end of the step      

//scaleing part
             curvefit/NTHR=0/K={segstart_2nddif1+seglength/2}/Q/W=2 exp_Xoffset p1_disp[segstart_2nddif1+seglength/2, segstart_2nddif1+2*seglength] /D
             V_tau[count][0] = W_coef[2]
             p1_scaled = (p1_disp[p+segstart_2nddif1]-W_coef[0])/-W_coef[1]
             for( l=0; l<3*seglength; l+=1)
                  if(l<seglength/2)
                       p1_scaledx[l] = l
                  else
                       p1_scaledx[l] = seglength/2+(l-seglength/2)/W_coef[2]*TC
                   endif          
              endfor     
             pairstepscaled[][0][count] = p1_scaled[p]
             pairstepscaled_x[][0][count] = p1_scaledx[p]


             
             calc2nddiff(p2_disp, seglength/4)
             duplicate/o dwave_2nddiff p2_2nddiff             
             segstart_2nddif2 = find2nddif_peak(p2_2nddiff, steploc2[j],seglength/2, seglength/4) - seglength/2    //for align at the starting point
//             segstart_2nddif2 = find2nddif_peak(p2_2nddiff, -steploc2[j],seglength/2, seglength/4) - 3*seglength/2     //for align at ending point             

             curvefit/NTHR=0/K={segstart_2nddif2+seglength/2}/Q/W=2 exp_Xoffset p2_disp[segstart_2nddif2+seglength/2, segstart_2nddif2+2*seglength] /D
             V_tau[count][1] = W_coef[2]
             p2_scaled = (p2_disp[p+segstart_2nddif2]-W_coef[0])/-W_coef[1]
             for( l=0; l<3*seglength; l+=1)
                  if(l<seglength/2)
                       p2_scaledx[l] = l
                  else
                       p2_scaledx[l] = seglength/2+(l-seglength/2)/W_coef[2]*TC
                   endif          
              endfor                   
             pairstepscaled[][1][count] = p2_scaled[p]
             pairstepscaled_x[][1][count] = p2_scaledx[p]

                       
             segstart1 = floor(abs(steploc1[j])) - seglength
             segstart2 = floor(abs(steploc2[j])) - seglength

             
//fitting with error function            
             baseline = p1_disp[segstart1] 
             stepheight = p1_disp[segstart1+2*seglength] - p1_disp[segstart1]
             stepwidth = 100
             stepcent = segstart1+seglength
             
             w_coef1 = {baseline, stepheight, stepcent, stepwidth}
             funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p1_disp[segstart1, segstart1+2*seglength] /D
             
             V_tau[count][0] = w_coef1[3]
             V_height[count][0] = w_coef1[1]

             baseline = p2_disp[segstart2] 
             stepheight = p2_disp[segstart2+2*seglength] - p2_disp[segstart2]
             stepwidth = 100
             stepcent = segstart2+seglength
             
             w_coef1 = {baseline, stepheight, stepcent, stepwidth}
             funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p2_disp[segstart2, segstart2+2*seglength] /D          
                
             V_tau[count][1] = w_coef1[3]
             V_height[count][1] = w_coef1[1]             
             



             
             pairstepdisp[][0][count] = pairdispwave[p+segstart1][0][pairnum]
             pairstepdisp[][1][count] = pairdispwave[p+segstart2][1][pairnum]
             
             pairstepdisp_2nddif[][0][count] = pairdispwave[p+segstart_2nddif1][0][pairnum]
             pairstepdisp_2nddif[][1][count] = pairdispwave[p+segstart_2nddif2][1][pairnum]
             
             
               peakdif_pair[count] = abs(steploc1[j]) - abs(steploc2[j]) 
             count+=1

             Dowindow/k Layoutpair
             Newlayout/K=1/N=layoutpair
             

             
             p1str = "pair1origw"
             p2str = "pair2origw"
             
             
             killwindow $p1str
             killwindow $p2str
             display/N=$p1str pairdispwave[][0][pairnum]
             SetAxis Bottom abs(steploc1[j]) - 450, abs(steploc1[j])+450
             wavestats/q/R=( abs(steploc1[j]) - 450, abs(steploc1[j])+450) p1_disp
             SetAxis left V_min, V_max
             appendtograph/W=$p1str fit_p1_disp
             ModifyGraph mode(fit_p1_disp)=4,rgb(fit_p1_disp)=(0,15872,65280)
             display/N=$p2str pairdispwave[][1][pairnum]
             SetAxis Bottom abs(steploc2[j]) - 450, abs(steploc2[j])+450
             wavestats/q/R=( abs(steploc2[j]) - 450, abs(steploc2[j])+450) p2_disp
             SetAxis left V_min, V_max             
             appendtograph/W=$p2str fit_p2_disp
             ModifyGraph mode(fit_p2_disp)=4,rgb(fit_p2_disp)=(0,15872,65280)

             appendlayoutobject/F=0/R=(60,21+100*1, 260, 120+100*1) graph $p1str
             appendlayoutobject/F=0/R=(281,21 + 100*1,500,120 + 100*1)  graph $p2str             
             

                  DoWindow/F Layoutpair
                  Notebook $NBname scaling = {90,90}, picture = {Layoutpair, -1,1}
                  Notebook $NBname text = "\r"
           
             endif
         endif
      endfor
endfor


redimension/n=(seglength*2,2,count) pairstepdisp, pairstepdisp_2nddif
redimension/n=(count) peakdif_pair


redimension/n=(count, 2) V_tau, V_height, V_terr
make/o/n=(seglength*3, 2, count) pairstepscaled, pairstepscaled_x


end          



function matchseg_short_alt()

wave steppair, steploc1, steploc2, selectpair, pairdispwave_seg, vpair
wave steploc1_old, steploc2_old

wave pairdispwave

variable npair = dimsize(selectpair,0)
variable nstep = dimsize(steppair,0)
variable seglength = 200
variable boxcarsize = 10
variable i,j, iframe

variable movielength = dimsize(pairdispwave_seg,0)

//wave fitres_ed_sub
//wave bondid, segid


make/o/n=(seglength*2, 2, nstep) pairstepdisp
make/o/n=(nstep,2) V_tau, V_height, V_terr
make/o/n=(seglength*2, 2 ,nstep) pairstepdisp_2nddif

make/o/n=(nstep) peakdif_pair, synpairid

variable segstart1, segstart2
variable pairnum

variable boxcarsmooth = 100


make/o/n=(boxcarsize) dispseg
make/o/n=(movielength) p1_smth, p2_smth, p1_disp, p2_disp
make/o/n=(movielength-boxcarsmooth) p1_v, p2_v

//generate scaled curve
make/o/n=(seglength*3) p1_scaled, p2_scaled, p1_scaledx, p2_scaledx
make/o/n=(seglength*3, 2, nstep) pairstepscaled, pairstepscaled_x

//parameters for fitting
make/o/n=4 w_coef1
variable baseline, stepheight, stepwidth, stepcent
variable tau_thresh = seglength*0.6
variable errthresh = 0.009

variable l
variable TC = 100


//variable heightthresh = 0.006
//variable heightratio = 0.7
variable heightthresh = 0.003
variable heightratio = 0.1
variable tauthresh = 600
variable h_p1, h_p2




variable count = 0
variable segstart_2nddif1, segstart_2nddif2

string  NBname, p1str, p2str
 
            sprintf NBname, "postpairsteps"
            DoWindow/K $NBname
            NewNotebook/K=1/F=1/N=$NBname
            
            
wave dwave_2nddiff            
wave w_coef


for(i=0; i<npair; i+=1)
    for(j=0; j<nstep;j+=1)
         if (steppair[j] == selectpair[i]) 
             pairnum = selectpair[i]          
           
             
             p1_disp = pairdispwave_seg[p][0][pairnum]
             p2_disp = pairdispwave_seg[p][1][pairnum] 
                         
             calc2nddiff(p1_disp, seglength/4)
             duplicate/o dwave_2nddiff p1_2nddiff
             segstart_2nddif1 = find2nddif_peak(p1_2nddiff, steploc1[j],seglength/2, seglength/4) - seglength/2        //for align at starting point


             
             calc2nddiff(p2_disp, seglength/4)
             duplicate/o dwave_2nddiff p2_2nddiff              
             segstart_2nddif2 = find2nddif_peak(p2_2nddiff, steploc2[j],seglength/2, seglength/4) - seglength/2        //for align at starting point
         
           
           
             segstart1 = floor(abs(steploc1[j])) - seglength
             segstart2 = floor(abs(steploc2[j])) - seglength
             
//fitting with error function            
             baseline = p1_disp[segstart_2nddif1] 
             stepheight = p1_disp[segstart_2nddif1+2*seglength] - p1_disp[segstart_2nddif1]
             stepwidth = tc
             stepcent = segstart1+seglength
             
             w_coef1 = {baseline, stepheight, stepcent, stepwidth}
             funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p1_disp[segstart1, segstart1+2*seglength] /D
                       
             

             //added 9/20/17   modify fitting range based on tau value
             if(w_coef1[3] > tau_thresh)
                   funcfit/NTHR=0/Q/W=2/H="1111" step_1, w_coef1 p1_disp[segstart1-seglength, segstart1+3*seglength] /D             
//                 funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p1_disp[segstart1-0.5*seglength, segstart1+2.5*seglength] /D     
             elseif(w_coef1[3] <0)     
                 w_coef1 = {baseline, stepheight, stepcent, stepwidth}       
                 funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p1_disp[segstart1+0.5*seglength, segstart1+1.5*seglength] /D   
             endif    
 
              V_tau[count][0] = w_coef1[3]
             V_height[count][0] = w_coef1[1]
             h_p1 = abs(w_coef1[1])
             
             V_terr[count][0] = V_chisq     
             
             if ( w_coef1[3] <5 || V_chisq>errthresh || w_coef1[3]>tauthresh || h_p1 < heightthresh)   
                continue
             endif       

             baseline = p2_disp[segstart_2nddif2] 
             stepheight = p2_disp[segstart_2nddif2+2*seglength] - p2_disp[segstart_2nddif2]
             stepwidth = tc
             stepcent = segstart2+seglength
             
             w_coef1 = {baseline, stepheight, stepcent, stepwidth}
             funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p2_disp[segstart2, segstart2+2*seglength]  /D            


             //added 9/20/17   modify fitting range based on tau value
             if(w_coef1[3] > tau_thresh)
                 funcfit/NTHR=0/Q/W=2/H="1111" step_1, w_coef1 p2_disp[segstart2-seglength, segstart2+3*seglength] /D             
//                 funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p2_disp[segstart2-0.5*seglength, segstart2+2.5*seglength] /D  
             elseif(w_coef1[3] <0)       
                 w_coef1 = {baseline, stepheight, stepcent, stepwidth}     
                 funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p2_disp[segstart2+0.5*seglength, segstart2+1.5*seglength] /D    
             endif                                   

             V_tau[count][1] = w_coef1[3]
             V_height[count][1] = w_coef1[1]  
              h_p2 = abs(w_coef1[1])
                        
              V_terr[count][1] = V_chisq 
              
             if ( w_coef1[3] <5 || V_chisq>errthresh || w_coef1[3]>tauthresh || h_p2 < heightthresh)   
                continue
             endif                

//           

            if( (min(h_p1,h_p2)/max(h_p1,h_p2)) < heightratio)
                 continue
             endif
             
                 
             peakdif_pair[count] = abs(steploc1[j]) - abs(steploc2[j])
              
              synpairid[count] = pairnum
 
             
             pairstepdisp[][0][count] = pairdispwave_seg[p+segstart1][0][pairnum]
             pairstepdisp[][1][count] = pairdispwave_seg[p+segstart2][1][pairnum]
             pairstepdisp_2nddif[][0][count] = pairdispwave_seg[p+segstart_2nddif1][0][pairnum]
             pairstepdisp_2nddif[][1][count] = pairdispwave_seg[p+segstart_2nddif2][1][pairnum]             
             
             
             count+=1
             
             Dowindow/k Layoutpair
             Newlayout/K=1/N=layoutpair
             

             
             p1str = "pair1origw"
             p2str = "pair2origw"
             
             
             killwindow $p1str
             killwindow $p2str
             display/N=$p1str pairdispwave_seg[][0][pairnum]
             appendtograph/W=$p1str fit_p1_disp
             ModifyGraph mode(fit_p1_disp)=4,rgb(fit_p1_disp)=(0,15872,65280)
             display/N=$p2str pairdispwave_seg[][1][pairnum]
             appendtograph/W=$p2str fit_p2_disp
             ModifyGraph mode(fit_p2_disp)=4,rgb(fit_p2_disp)=(0,15872,65280)

             appendlayoutobject/F=0/R=(60,21+100*1, 260, 120+100*1) graph $p1str
             appendlayoutobject/F=0/R=(281,21 + 100*1,500,120 + 100*1)  graph $p2str             
             
             
                  DoWindow/F Layoutpair
                  Notebook $NBname scaling = {90,90}, picture = {Layoutpair, -1,1}
                  Notebook $NBname text = "\r"             
             
             
         endif
      endfor
endfor


redimension/n=(seglength*2,2,count) pairstepdisp
redimension/n=(seglength*2,2,count) pairstepdisp_2nddif
redimension/n=(count, 2) V_tau, v_height, V_terr
redimension/n=(count) peakdif_pair, synpairid

make/o/n=(seglength/boxcarsize*2, 2,count) pairstepdisp_avg, pairstepdisp_2nddifavg

//for scaling
make/o/n=(seglength*3, 2, count) pairstepscaled, pairstepscaled_x
make/o/n=(seglength*3/boxcarsize, 2, count) pairstepscaled_avg, pairstepscaled_x_avg
//




//

end             




   



function matchseg_short()

wave steppair, steploc1, steploc2, selectpair, pairdispwave, vpair
wave steploc1_old, steploc2_old

variable npair = dimsize(selectpair,0)
variable nstep = dimsize(steppair,0)
//variable seglength = 200
variable seglength = 2000
variable boxcarsize = 10
variable i,j, iframe

variable movielength = dimsize(pairdispwave,0)

make/o/n=(seglength*2, 2, nstep) pairstepdisp
make/o/n=(nstep,2) V_tau, V_height, V_terr
make/o/n=(seglength*2, 2 ,nstep) pairstepdisp_2nddif

make/o/n=(nstep) peakdif_pair

variable segstart1, segstart2
variable pairnum

variable boxcarsmooth = 1000

make/o/n=(boxcarsize) dispseg
make/o/n=(movielength) p1_smth, p2_smth, p1_disp, p2_disp
make/o/n=(movielength-boxcarsmooth) p1_v, p2_v

//generate scaled curve
make/o/n=(seglength*3) p1_scaled, p2_scaled, p1_scaledx, p2_scaledx
make/o/n=(seglength*3, 2, nstep) pairstepscaled, pairstepscaled_x

//parameters for fitting
make/o/n=4 w_coef1
variable baseline, stepheight, stepwidth, stepcent
variable tau_thresh = seglength*0.6
variable errthresh = 0.09

variable l
variable TC = 1000





variable count = 0
variable segstart_2nddif1, segstart_2nddif2

string  NBname, p1str, p2str
 
            sprintf NBname, "postpairsteps"
            DoWindow/K $NBname
            NewNotebook/K=1/F=1/N=$NBname
            
            
wave dwave_2nddiff            
wave w_coef


for(i=0; i<npair; i+=1)
    for(j=0; j<nstep;j+=1)
         if (steppair[j] == selectpair[i]) 
             pairnum = selectpair[i]          
             
             p1_disp = pairdispwave[p][0][pairnum]
             p2_disp = pairdispwave[p][1][pairnum] 
                         
             calc2nddiff(p1_disp, seglength/4)
             duplicate/o dwave_2nddiff p1_2nddiff
             segstart_2nddif1 = find2nddif_peak(p1_2nddiff, steploc1[j],seglength/2, seglength/4) - seglength/2        //for align at starting point
//             segstart_2nddif1 = find2nddif_peak(p1_2nddiff, -steploc1[j],seglength/2, seglength/4) - 3*seglength/2        //for align at ending point

//scaleing part
             curvefit/NTHR=0/K={segstart_2nddif1+seglength/2}/Q/W=2 exp_Xoffset p1_disp[segstart_2nddif1+seglength/2, segstart_2nddif1+2*seglength] /D
//             V_tau[count][0] = W_coef[2]
             p1_scaled = (p1_disp[p+segstart_2nddif1]-W_coef[0])/-W_coef[1]          
             for( l=0; l<3*seglength; l+=1)
                  if(l<seglength/2)
                       p1_scaledx[l] = l
                  else
                       p1_scaledx[l] = seglength/2+(l-seglength/2)/W_coef[2]*1000
                   endif          
              endfor     
             pairstepscaled[][0][count] = p1_scaled[p]
             pairstepscaled_x[][0][count] = p1_scaledx[p]


             
             calc2nddiff(p2_disp, seglength/4)
             duplicate/o dwave_2nddiff p2_2nddiff              
             segstart_2nddif2 = find2nddif_peak(p2_2nddiff, steploc2[j],seglength/2, seglength/4) - seglength/2        //for align at starting point
//             segstart_2nddif2 = find2nddif_peak(p2_2nddiff, -steploc2[j],seglength/2, seglength/4) - 3*seglength/2     //for align at ending point

//scaling part
             curvefit/NTHR=0/K={segstart_2nddif2+seglength/2}/Q/W=2 exp_Xoffset p2_disp[segstart_2nddif2+seglength/2, segstart_2nddif2+2*seglength] /D
//             V_tau[count][1] = W_coef[2]
             p2_scaled = (p2_disp[p+segstart_2nddif2]-W_coef[0])/-W_coef[1]
             for( l=0; l<3*seglength; l+=1)
                  if(l<seglength/2)
                       p2_scaledx[l] = l
                  else
                       p2_scaledx[l] = seglength/2+(l-seglength/2)/W_coef[2]*1000
                   endif          
              endfor                   
             pairstepscaled[][1][count] = p2_scaled[p]
             pairstepscaled_x[][1][count] = p2_scaledx[p]

           
           
           
             segstart1 = floor(abs(steploc1[j])) - seglength
             segstart2 = floor(abs(steploc2[j])) - seglength
             
//fitting with error function            
             baseline = p1_disp[segstart1] 
             stepheight = p1_disp[segstart1+2*seglength] - p1_disp[segstart1]
             stepwidth = 1000
             stepcent = segstart1+seglength
             
             w_coef1 = {baseline, stepheight, stepcent, stepwidth}
             funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p1_disp[segstart1, segstart1+2*seglength] /D
             
             //added 9/20/17   modify fitting range based on tau value  
             
             if(w_coef1[3] > tau_thresh)
                   funcfit/NTHR=0/Q/W=2/H="1111" step_1, w_coef1 p1_disp[segstart1-seglength, segstart1+3*seglength] /D             
//                 funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p1_disp[segstart1-0.5*seglength, segstart1+2.5*seglength] /D     
             elseif(w_coef1[3] <0)     
                 w_coef1 = {baseline, stepheight, stepcent, stepwidth}       
                 funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p1_disp[segstart1+0.5*seglength, segstart1+1.5*seglength] /D   
             endif    
 
              V_tau[count][0] = w_coef1[3]
             V_height[count][0] = w_coef1[1]             
             V_terr[count][0] = V_chisq    
              
             
             if ( w_coef1[3] <1 || V_chisq>errthresh)   
                continue
             endif                    
             

             baseline = p2_disp[segstart2] 
             stepheight = p2_disp[segstart2+2*seglength] - p2_disp[segstart2]
             stepwidth = 1000
             stepcent = segstart2+seglength
             
             w_coef1 = {baseline, stepheight, stepcent, stepwidth}
             funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p2_disp[segstart2, segstart2+2*seglength] /D             


             //added 9/20/17   modify fitting range based on tau value
             if(w_coef1[3] > tau_thresh)
                 funcfit/NTHR=0/Q/W=2/H="1111" step_1, w_coef1 p2_disp[segstart2-seglength, segstart2+3*seglength] /D             
//                 funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p2_disp[segstart2-0.5*seglength, segstart2+2.5*seglength] /D  
             elseif(w_coef1[3] <0)       
                 w_coef1 = {baseline, stepheight, stepcent, stepwidth}     
                 funcfit/NTHR=0/Q/W=2 step_1, w_coef1 p2_disp[segstart2+0.5*seglength, segstart2+1.5*seglength] /D    
             endif                                   

             V_tau[count][1] = w_coef1[3]
             V_height[count][1] = w_coef1[1]             
              V_terr[count][1] = V_chisq 
              
             if ( w_coef1[3] <1 || V_chisq>errthresh)   
                continue
             endif 
//             

//             
             peakdif_pair[count] = abs(steploc1[j]) - abs(steploc2[j])
             
             pairstepdisp[][0][count] = pairdispwave[p+segstart1][0][pairnum]
             pairstepdisp[][1][count] = pairdispwave[p+segstart2][1][pairnum]
             pairstepdisp_2nddif[][0][count] = pairdispwave[p+segstart_2nddif1][0][pairnum]
             pairstepdisp_2nddif[][1][count] = pairdispwave[p+segstart_2nddif2][1][pairnum]             
             
             
             count+=1
             
             Dowindow/k Layoutpair
             Newlayout/K=1/N=layoutpair
             

             
             p1str = "pair1origw"
             p2str = "pair2origw"
             
             
             killwindow $p1str
             killwindow $p2str
             display/N=$p1str pairdispwave[][0][pairnum]
             appendtograph/W=$p1str fit_p1_disp
             ModifyGraph mode(fit_p1_disp)=4,rgb(fit_p1_disp)=(0,15872,65280)
             display/N=$p2str pairdispwave[][1][pairnum]
             appendtograph/W=$p2str fit_p2_disp
             ModifyGraph mode(fit_p2_disp)=4,rgb(fit_p2_disp)=(0,15872,65280)

             appendlayoutobject/F=0/R=(60,21+100*1, 260, 120+100*1) graph $p1str
             appendlayoutobject/F=0/R=(281,21 + 100*1,500,120 + 100*1)  graph $p2str             
             
             boxsmooth boxcarsmooth, p1_disp, p1_smth
             boxsmooth boxcarsmooth, p2_disp, p2_smth
             
             p1str = "pair1smthw"
             p2str = "pair2smthw"


             killwindow $p1str
             killwindow $p2str             
             display/N=$p1str p1_smth
             display/N=$p2str p2_smth

             appendlayoutobject/F=0/R=(60,21+100*2, 260, 120+100*2) graph $p1str
             appendlayoutobject/F=0/R=(281,21 + 100*2,500,120 + 100*2)  graph $p2str     
             
             
             p1str = "pair1vw"
             p2str = "pair2vw"


             killwindow $p1str
             killwindow $p2str            
             p1_v = vpair[p][0][pairnum]
             p2_v = vpair[p][1][pairnum]
              
             display/N=$p1str p1_v
             SetDrawEnv/W=$p1str xcoord= bottom, linefgc= (0,15872,65280)
             DrawLine abs(steploc1_old[j])-boxcarsmooth/2,-0.05,abs(steploc1_old[j])-boxcarsmooth/2,1.08  
             SetDrawEnv/W=$p1str xcoord= bottom, linefgc= (0,0,0)             
             DrawLine abs(steploc1[j])-boxcarsmooth/2,-0.05,abs(steploc1[j])-boxcarsmooth/2,1.08  
             display/N=$p2str p2_v
             SetDrawEnv/W=$p2str xcoord= bottom, linefgc=  (0,15872,65280)
             DrawLine abs(steploc2_old[j])-boxcarsmooth/2,-0.05,abs(steploc2_old[j])-boxcarsmooth/2,1.08   
             SetDrawEnv/W=$p2str xcoord= bottom, linefgc=  (0,0,0)
             DrawLine abs(steploc2[j])-boxcarsmooth/2,-0.05,abs(steploc2[j])-boxcarsmooth/2,1.08   
                          
                                      
             appendlayoutobject/F=0/R=(60,21+100*3, 260, 120+100*3) graph $p1str
             appendlayoutobject/F=0/R=(281,21 + 100*3,500,120 + 100*3)  graph $p2str                           
             
                  DoWindow/F Layoutpair
                  Notebook $NBname scaling = {90,90}, picture = {Layoutpair, -1,1}
                  Notebook $NBname text = "\r"             
             
             
         endif
      endfor
endfor


redimension/n=(seglength*2,2,count) pairstepdisp
redimension/n=(seglength*2,2,count) pairstepdisp_2nddif
redimension/n=(count, 2) V_tau, v_height, V_terr
redimension/n=(count) peakdif_pair

make/o/n=(seglength/boxcarsize*2, 2,count) pairstepdisp_avg, pairstepdisp_2nddifavg

//for scaling
make/o/n=(seglength*3, 2, count) pairstepscaled, pairstepscaled_x
make/o/n=(seglength*3/boxcarsize, 2, count) pairstepscaled_avg, pairstepscaled_x_avg
//



for(i=0; i<count; i+=1)
      j=0
      for ( iframe= 0; iframe<(2*seglength); iframe+=boxcarsize)
            dispseg = pairstepdisp[iframe+p][0][i]
            pairstepdisp_avg[j][0][i] = mean(dispseg)
            dispseg = pairstepdisp[iframe+p][1][i]
            pairstepdisp_avg[j][1][i] = mean(dispseg)
            dispseg = pairstepdisp_2nddif[iframe+p][0][i]
            pairstepdisp_2nddifavg[j][0][i] = mean(dispseg)
            dispseg = pairstepdisp_2nddif[iframe+p][1][i]
            pairstepdisp_2nddifavg[j][1][i] = mean(dispseg)    
          
            j+=1
       endfor
endfor            



//for scaling 

for(i=0; i<count; i+=1)
      j=0
      for ( iframe= 0; iframe<(3*seglength); iframe+=boxcarsize)
            dispseg = pairstepscaled[iframe+p][0][i]
            pairstepscaled_avg[j][0][i] = mean(dispseg)
            dispseg = pairstepscaled[iframe+p][1][i]
            pairstepscaled_avg[j][1][i] = mean(dispseg)
            dispseg = pairstepscaled_x[iframe+p][0][i]
            pairstepscaled_x_avg[j][0][i] = mean(dispseg)/10
            dispseg =  pairstepscaled_x[iframe+p][1][i]
            pairstepscaled_x_avg[j][1][i] = mean(dispseg)/10    
          
            j+=1
       endfor
endfor   

//

end             












function bondcorr_plotalot(datawave, bondpairwave)

wave datawave
wave bondpairwave

variable ipair
variable npair = dimsize(datawave,0)
variable wavelength = dimsize(bondpairwave, 0)

variable corval


 variable ipcount = 0, plotsoncurrentpage = 0
 variable plotsperpage = 7
 wave msdxwave_diff
 
 
 string  NBname, p1str, p2str
 
            sprintf NBname, "postpairplotalot"
            DoWindow/K $NBname
            NewNotebook/K=1/F=1/N=$NBname

make/o/n=(wavelength)  pairwave1, pairwave2
make/o/n=(wavelength, 2*plotsperpage) colorwave            
            
for ( ipair = 0; ipair<npair; ipair+=1)             

          pairwave1 = bondpairwave[p][0][ipair]
          pairwave2 = bondpairwave[p][1][ipair]
          

          
          SetScale/P x 0,0.01,"",pairwave1, pairwave2
            
         if (mod(ipcount, plotsperpage) == 0)
            plotsoncurrentpage = 0
            DoWindow/K LayoutXY
            NewLayout/K=1/N=LayoutXY
        endif    


          corval = datawave[ipair][0]
          colorwave[][plotsoncurrentpage] = corval
          colorwave[][plotsoncurrentpage+plotsperpage] = corval
      
    
            sprintf p1str,"Plot%dp1",plotsoncurrentpage
            sprintf p2str,"Plot%dp2",plotsoncurrentpage
            
            dispwavepair(pairwave1,ipair,p1str, plotsoncurrentpage,corval, colorwave)
            
            dispwavepair(pairwave2,ipair, p2str, plotsoncurrentpage+plotsperpage,corval,colorwave)


   	    AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $p1str
      	    AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $p2str    
              ipcount = ipcount+1
             plotsoncurrentpage =  plotsoncurrentpage + 1
                if ((plotsoncurrentpage == 7) || (ipair == npair))
                  DoWindow/F LayoutXY
             //if (printflag == 1)
                   //PrintLayout LayoutXY  //Modified Jul 27
                       Notebook $NBname scaling = {90,90}, picture = {LayoutXY, -1,1}
                       Notebook $NBname text = "\r"
               endif        
               

endfor               
               
               
               
               
SaveNotebook/P=ToPractice $NBname


end





function dispwavepair(datawavep, ipost,  swname, nplot, corval,cwave)
wave datawavep, cwave
variable ipost,nplot
string swname
variable corval

string dws,s1

sprintf dws, "dw%d", nplot
duplicate/o datawavep $dws





    dowindow/K $swname
    
    display/K=1/N=$swname $dws 
    SetAxis left, wavemin($dws),wavemax($dws)

    sprintf s1, "pairid=%d, corval = %g",ipost, corval
    TextBox/C/B=1/N=text0/A=LT s1
    ModifyGraph tick=2,mirror=1,standoff=0
    ModifyGraph mode=0, zColor={cwave[][nplot],-1,1,Rainbow,0}
    ModifyGraph margin(left)=50,margin(bottom)=23,margin(top)=2,margin(right)=5;
    ModifyGraph width=165,height=75     

end


function altcorrelation(datawave1, datawave2)
wave datawave1, datawave2

variable wavelength = dimsize(datawave1,0)
variable  tlag=3000
variable i

variable corrcoef

make/o/n=(wavelength-tlag) dwave1_bg, dwave1_ed, dwave2_bg, dwave2_ed, dwave1_diff, dwave2_diff
make/o/n=(wavelength-tlag) dwave1_diff_sq, dwave2_diff_sq, dwave12_diff

duplicate/o/R=(0, wavelength-tlag-1) datawave1, dwave1_bg
duplicate/o/R=(0,wavelength-tlag-1) datawave2, dwave2_bg
duplicate/o/R=(tlag, wavelength-1) datawave1, dwave1_ed
duplicate/o/R=(tlag, wavelength-1) datawave2, dwave2_ed

dwave1_diff = dwave1_ed - dwave1_bg
dwave2_diff = dwave2_ed - dwave2_bg

dwave1_diff_sq = dwave1_diff*dwave1_diff
dwave2_diff_sq = dwave2_diff*dwave2_diff
dwave12_diff = dwave1_diff*dwave2_diff

      
//corrcoef = sum(dwave12_diff)/sqrt(sum(dwave1_diff_sq)*sum(dwave2_diff_sq))

corrcoef = statscorrelation(dwave1_diff, dwave2_diff)

return(corrcoef)

end


function calccov(dwave1, dwave2)
wave dwave1, dwave2
variable cov


duplicate/o dwave1 dwave1_ed
duplicate/o dwave2 dwave2_ed

dwave1_ed = dwave1 - mean(dwave1)
dwave2_ed = dwave2 - mean(dwave2)

duplicate/o dwave1 dwavecross

dwavecross = dwave1_ed*dwave2_ed

cov = sum(dwavecross)

return(cov)

end


function calcnfactor(dwave1, dwave2)
wave dwave1 , dwave2
variable nfactor

duplicate/o dwave1 dwave1_ed
duplicate/o dwave2 dwave2_ed

dwave1_ed = dwave1 - mean(dwave1)
dwave2_ed = dwave2 - mean(dwave2)

duplicate/o dwave1 dwave1sq, dwave2sq

dwave1sq = dwave1_ed^2
dwave2sq = dwave2_ed^2

nfactor = sqrt(sum(dwave1sq)*sum(dwave2sq))

return(nfactor)

end





function selectacpair(datawave, bondpairwave_p, bondpairwave_t, corthresh)
wave datawave, bondpairwave_p, bondpairwave_t
variable corthresh
variable ratiothresh = 10

variable ipair
variable npair = dimsize(datawave,0)
variable movielength = dimsize(bondpairwave_p, 0)
variable paircount = 0

make/o/n=(movielength) temppairwave
make/o/n=(npair) selectpair, pairratio

variable diff_p1, diff_p2, diff_t1, diff_t2
variable pthresh,tthresh

pthresh = 0.02
tthresh = 0.03

for (ipair=0; ipair<npair; ipair+=1)
//       if(abs(datawave[ipair][0]) > corthresh)

          temppairwave = bondpairwave_p[p][0][ipair]
          wavestats/q temppairwave
          diff_p1 = V_max - V_min
          temppairwave = bondpairwave_p[p][1][ipair]
          wavestats/q temppairwave
          diff_p2 = V_max - V_min
          
                         
          pairratio[ipair] = max(diff_p1/diff_p2, diff_p2/diff_p1)        
         if(datawave[ipair][0]<-corthresh && pairratio[ipair]<ratiothresh)
          if ((diff_p1>pthresh) && (diff_p2>pthresh))
               temppairwave = bondpairwave_t[p][0][ipair]
               wavestats/q temppairwave
               diff_t1 = V_max - V_min
               temppairwave = bondpairwave_t[p][1][ipair]
               wavestats/q temppairwave
               diff_t2 = V_max - V_min
               if ((diff_t1<tthresh)&& (diff_t2<tthresh))
                     selectpair[paircount] = ipair
                      paircount+=1
               endif
          endif
     
     
     endif

endfor          


redimension/n=(paircount)  selectpair

end



function selectacpair_cat(datawave, datawave_t, bondpairwave_p, bondpairwave_t, corthresh)
wave datawave, datawave_t, bondpairwave_p, bondpairwave_t
variable corthresh
variable ratiothresh = 2

variable ipair
variable npair = dimsize(datawave,0)
variable movielength = dimsize(bondpairwave_p, 0)
variable paircount = 0
make/o/n=4 paircount_cat = 0

make/o/n=(movielength) temppairwave
make/o/n=(npair,4) selectpair_cat=0, pairratio

variable V_flag

variable diff_p1, diff_p2, diff_t1, diff_t2
variable pthresh,tthresh

pthresh = 0.02
tthresh = 0.03

for (ipair=0; ipair<npair; ipair+=1)
//       if(abs(datawave[ipair][0]) > corthresh)

          temppairwave = bondpairwave_p[p][0][ipair]
          wavestats/q temppairwave
          diff_p1 = V_max - V_min
          temppairwave = bondpairwave_p[p][1][ipair]
          wavestats/q temppairwave
          diff_p2 = V_max - V_min

          temppairwave = bondpairwave_t[p][0][ipair]
          wavestats/q temppairwave
          diff_t1 = V_max - V_min
          temppairwave = bondpairwave_t[p][1][ipair]
          wavestats/q temppairwave
          diff_t2 = V_max - V_min          
          
          
          V_flag = 5
                         
          pairratio[ipair] = max(diff_p1/diff_p2, diff_p2/diff_p1)        
      if( pairratio[ipair]<ratiothresh) 
          if ((diff_p1>pthresh) && (diff_p2>pthresh)&&(diff_t1<tthresh)&& (diff_t2<tthresh))
              if(datawave[ipair][0]<-corthresh)
                   V_flag = 0   
                   paircount = paircount_cat[V_flag]
                   selectpair_cat[paircount][V_flag] = ipair
                   paircount_cat[V_flag]+=1                    
               elseif(datawave[ipair][0]>corthresh)
                    V_flag = 1 
                   paircount = paircount_cat[V_flag]
                   selectpair_cat[paircount][V_flag] = ipair
                   paircount_cat[V_flag]+=1                         
               endif                     
           elseif ((diff_p1<pthresh) && (diff_p2<pthresh)&&(diff_t1>tthresh)&& (diff_t2>tthresh))   
              if(datawave_t[ipair][0]<-corthresh)
                   V_flag = 2   
                   paircount = paircount_cat[V_flag]
                   selectpair_cat[paircount][V_flag] = ipair
                   paircount_cat[V_flag]+=1                    
               elseif(datawave_t[ipair][0]>corthresh)
                    V_flag = 3     
                   paircount = paircount_cat[V_flag]
                   selectpair_cat[paircount][V_flag] = ipair
                   paircount_cat[V_flag]+=1                     
               endif                           
            endif     

                 
        endif          
endfor          

wavestats/q paircount_cat
redimension/n=(V_max,4)  selectpair_cat

end











function bifurcatetail(bonddispw, bondcorrw,option)
wave bonddispw, bondcorrw
variable option            //1 for only highly anticorrelated, 0 for all with large displacements

variable nbond = dimsize(bonddispw, 0)
variable i

if(option == 1)

for( i= 0; i<nbond; i+=1) 
       if (bondcorrw[i][0] > -0.8)
           bonddispw[i][0] = 0 
           bonddispw[i][1] = 0
       endif
endfor

elseif (option ==0)

for( i= 0; i<nbond; i+=1) 
       if (abs(bonddispw[i][0])<0.005 || abs(bonddispw[i][1])<0.005 )
           bonddispw[i][0] = 0 
           bonddispw[i][1] = 0
       endif
endfor


elseif (option==2)

for( i= 0; i<nbond; i+=1) 
       if (abs(bonddispw[i][0])<0.005 || abs(bonddispw[i][1])<0.005||bondcorrw[i][0] > -0.8)
           bonddispw[i][0] = 0 
           bonddispw[i][1] = 0
       endif
endfor



endif

end           
       



function quadstats(datawave)
wave datawave

variable np = dimsize(datawave,0)
variable i
make/o/n=4 quadcount=0


for(i=0; i<np; i+=1)
      if (datawave[i][0] >0 && datawave[i][1] >0)
            quadcount[0] +=1
      elseif (datawave[i][0]<0 && datawave[i][1]>0)
             quadcount[1] +=1
      elseif (datawave[i][0]<0 && datawave[i][1] <0)
              quadcount[2]+=1
      elseif( datawave[i][0]>0 && datawave[i][1] <0) 
               quadcount[3]+=1
      endif
endfor

end          


function overlaysp()
wave selectpair                              
wave disp_cov, disp_nfactor

variable npair = dimsize(disp_cov, 0)
make/o/n=(npair) disp_cov_select, disp_nfactor_select
variable i
variable j=0

for (i=0; i<npair; i+=1)
      if (i==selectpair[j])
          disp_cov_select[j] = disp_cov[i]
          disp_nfactor_select[j] = disp_nfactor[i]
          j+=1
      endif

endfor          

redimension/N=(j) disp_cov_select, disp_nfactor_select

end





function bonddist(posttype, datawave)

variable posttype      // 1 for low traction force posts, 6 for high traction force posts
wave datawave     


wave astr
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost
variable movielength = dimsize(datawave,2)
variable freq = 200
//variable freq = 1
wave astr

variable angleofrotation = astr[2]

NVAR  inputmovieID_G

wave zeropos


wave linmaskval_intercept_int
variable i, j, k, l 
variable postid
make/o/n=3 postpairid
variable postpairid_cor
make/o/n=(npost*3, movielength) disp_dist
make/o/n=(npost*3) bond_x, bond_y
make/o/n=(movielength)  pairwave1x, pairwave1y, pairwave2x, pairwave2y
make/o/n=(movielength)  pairwave1p, pairwave1t, pairwave2p, pairwave2t

variable iframe
variable nframe = floor(movielength/freq)
//variable nframe = floor(movielength/freq)

variable dx, dy, theta

variable pairid = 0 

variable alpha

//removing posts at edge for the moment)
for (i=1; i<nxpost; i+=1)
      for(j=1; j<nypost; j+=1)
            postid = i*nypost+j
            if (linmaskval_intercept_int[postid] == posttype)
                postpairid[0] = (i-1)*nypost+ j
                postpairid[1] = i*nypost+j-1
                postpairid[2] = (i+2*(mod(j,2)-0.5))*nypost + j-1
                for( k = 0; k<3;k+=1)
                       postpairid_cor = postpairid[k] 
                       if (linmaskval_intercept_int[postpairid_cor] == posttype)
                           bond_x[pairid] = 0.5*(zeropos[postid][0]+zeropos[postpairid_cor][0])
                           bond_y[pairid] = 0.5*(zeropos[postid][1]+zeropos[postpairid_cor][1])
                           dx = zeropos[postid][0] - zeropos[postpairid_cor][0]
                           dy = zeropos[postid][1] - zeropos[postpairid_cor][1]
                           theta = acos (dx/sqrt(dx^2+dy^2))
                           
                           pairwave1x = datawave[postid][0][p]
                           pairwave2x = datawave[postpairid_cor][0][p] 
                           pairwave1y = datawave[postid][1][p]
                           pairwave2y = datawave[postpairid_cor][1][p]
                           pairwave1p = pairwave1x*cos(theta)+pairwave1y*sin(theta)
                           pairwave1t = -pairwave1x*sin(theta)+pairwave1y*cos(theta)
                           pairwave2p = pairwave2x*cos(theta)+pairwave2y*sin(theta)
                           pairwave2t =  -pairwave2x*sin(theta)+pairwave2y*cos(theta)    
                              
                           

                            for(iframe=0; iframe<nframe; iframe+=1)     
                                  disp_dist[pairid][iframe] = pairwave1p[iframe*freq]-pairwave2p[iframe*freq]
                            endfor      
                                 

                           pairid+=1
                       endif
                 endfor    
               endif
          endfor                 
endfor                           


redimension/n=(pairid, nframe) disp_dist
redimension/n=(pairid) bond_x, bond_y

Newmovie/o/P=ToPractice/F=(10)/L as "bonddist.mov"


playmovieaction setfrontmovie= inputMovieID_G
playmovieaction  stop, frame=0




for(iframe =0;iframe<nframe;iframe+=1)
     playmovieaction frame = iframe*freq
     playmovieaction extract
     
             imagetransform rgb2gray M_movieframe
             duplicate/O M_rgb2gray image1 
		Duplicate/O image1 OrigImage1
		ImageRotate/A=(AngleOfRotation) OrigImage1
		Duplicate/O M_RotatedImage RotatedOrig1
		     
     dowindow/k bonddist0
     newimage/N=bonddist0 rotatedorig1
     appendtograph/T/W=bonddist0 bond_y vs bond_x
     ModifyGraph mode=3
     ModifyGraph zColor(bond_y)={disp_dist[*][iframe],-0.15,0.15,Rainbow,0}
     ModifyGraph marker=19,msize=1.5
     Addmovieframe
endfor   

closemovie  



end


function bonddist1(posttype, datawave)

variable posttype      // 1 for low traction force posts, 6 for high traction force posts
wave datawave     


wave astr
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost
variable movielength = dimsize(datawave,2)
variable freq = 1
wave astr

variable angleofrotation = astr[2]

NVAR  inputmovieID_G

wave zeropos


wave linmaskval_intercept_int_sc
variable i, j, k, l 
variable postid
make/o/n=3 postpairid
variable postpairid_cor
make/o/n=(npost*3, movielength) disp_dist
make/o/n=(movielength) disp_dist_avg=0
make/o/n=(npost*3) bond_x, bond_y
make/o/n=(movielength)  pairwave1x, pairwave1y, pairwave2x, pairwave2y
make/o/n=(movielength)  pairwave1p, pairwave1t, pairwave2p, pairwave2t

variable iframe
variable nframe = floor(movielength/freq)


variable dx, dy, theta

variable pairid = 0 

variable alpha

//removing posts at edge for the moment)
for (i=1; i<nxpost; i+=1)
      for(j=1; j<nypost; j+=1)
            postid = i*nypost+j
            if (linmaskval_intercept_int_sc[postid] == posttype)
                postpairid[0] = (i-1)*nypost+ j
                postpairid[1] = i*nypost+j-1
                postpairid[2] = (i+2*(mod(j,2)-0.5))*nypost + j-1
                for( k = 0; k<3;k+=1)
                       postpairid_cor = postpairid[k] 
                       if (linmaskval_intercept_int_sc[postpairid_cor] == posttype)
                           bond_x[pairid] = 0.5*(zeropos[postid][0]+zeropos[postpairid_cor][0])
                           bond_y[pairid] = 0.5*(zeropos[postid][1]+zeropos[postpairid_cor][1])
                           dx = zeropos[postid][0] - zeropos[postpairid_cor][0]
                           dy = zeropos[postid][1] - zeropos[postpairid_cor][1]
                           theta = acos (dx/sqrt(dx^2+dy^2))
                           
                           pairwave1x = datawave[postid][0][p]
                           pairwave2x = datawave[postpairid_cor][0][p] 
                           pairwave1y = datawave[postid][1][p]
                           pairwave2y = datawave[postpairid_cor][1][p]
                           pairwave1p = pairwave1x*cos(theta)+pairwave1y*sin(theta)
                           pairwave1t = -pairwave1x*sin(theta)+pairwave1y*cos(theta)
                           pairwave2p = pairwave2x*cos(theta)+pairwave2y*sin(theta)
                           pairwave2t =  -pairwave2x*sin(theta)+pairwave2y*cos(theta)    
                              
                           

                            for(iframe=0; iframe<nframe; iframe+=1)     
                                  disp_dist[pairid][iframe] = pairwave1p[iframe*freq]-pairwave2p[iframe*freq]
                                  disp_dist_avg[iframe] += pairwave1p[iframe]-pairwave2p[iframe]
                            endfor      
                                 

                           pairid+=1
                       endif
                 endfor    
               endif
          endfor                 
endfor                           


redimension/n=(pairid, nframe) disp_dist
redimension/n=(pairid) bond_x, bond_y

disp_dist_avg = disp_dist_avg/pairid


end



function bonddist1_edge(datawave)

wave datawave  

wave linmaskval_intercept_int_sc, bondid

variable posttype      // 1 for low traction force posts, 6 for high traction force posts
   
variable npair = dimsize(bondid, 0)
variable nframe = dimsize(datawave, 2)

variable ipair, postA, postB, iframe,count=0
make/o/n=(nframe) pairwave1p, pairwave2p

wave pairdispwave
//make/o/n=(npair, nframe) disp_dist_edge
make/o/n=(nframe) disp_dist_avg_edge

for (ipair = 0; ipair<npair; ipair+=1)
      postA = bondid[ipair][0]
      postB = bondid[ipair][1]
      if(linmaskval_intercept_int_sc[postA] != linmaskval_intercept_int_sc[postB])
           pairwave1p = pairdispwave[p][0][ipair]
           pairwave2p = pairdispwave[p][1][ipair]
           count+=1

                            for(iframe=0; iframe<nframe; iframe+=1)     
  //                                disp_dist_edge[ipair][iframe] = pairwave1p[iframe]-pairwave2p[iframe]
                                  disp_dist_avg_edge[iframe] += pairwave1p[iframe]-pairwave2p[iframe]
                            endfor               
       endif
endfor
           


disp_dist_avg_edge = disp_dist_avg_edge/count


end









function selectacpair_seg(datawave, bondpairwave_p, bondpairwave_t, corthresh)
wave datawave, bondpairwave_p, bondpairwave_t
variable corthresh
variable ratiothresh = 2

variable ipair, iseg
variable npair = dimsize(datawave,0)
variable nseg = dimsize(datawave,1)
variable movielength = dimsize(bondpairwave_p, 0)
variable segmentlength = floor(movielength/nseg)
variable paircount = 0

make/o/n=(segmentlength) temppairwave
make/o/n=(npair) selectpair_seg, pairratio
make/o/n=(npair) selectpair_segnumber


variable diff_p1, diff_p2, diff_t1, diff_t2
variable pthresh,tthresh

pthresh = 0.02
tthresh = 0.03

for (ipair=0; ipair<npair; ipair+=1)
//       if(abs(datawave[ipair][0]) > corthresh)
      for(iseg=0; iseg<nseg; iseg+=1)
      
          temppairwave = bondpairwave_p[p+iseg*segmentlength][0][ipair]
          wavestats/q temppairwave
          diff_p1 = V_max - V_min
          temppairwave = bondpairwave_p[p+iseg*segmentlength][1][ipair]
          wavestats/q temppairwave
          diff_p2 = V_max - V_min
          
                         
          pairratio[ipair] = max(diff_p1/diff_p2, diff_p2/diff_p1)        
         if(datawave[ipair][iseg]<-corthresh && pairratio[ipair]<ratiothresh)
          if ((diff_p1>pthresh) && (diff_p2>pthresh))
               temppairwave = bondpairwave_t[p+iseg*segmentlength][0][ipair]
               wavestats/q temppairwave
               diff_t1 = V_max - V_min
               temppairwave = bondpairwave_t[p+iseg*segmentlength][1][ipair]
               wavestats/q temppairwave
               diff_t2 = V_max - V_min
               if ((diff_t1<tthresh)&& (diff_t2<tthresh))
                     selectpair_seg[paircount] = ipair
                     selectpair_segnumber[paircount] = iseg
                      paircount+=1
                endif
          endif
     
     
       endif
     
     
     endfor

endfor          


redimension/n=(paircount)  selectpair_seg, selectpair_segnumber

end



function findmultipeak(dataw, peakthresh,peakdifthresh)
wave dataw
variable peakthresh, peakdifthresh


variable boxsize =peakdifthresh*0.2
variable movielength = dimsize(dataw,0)
variable sig1
variable i,j
variable Rstart = 0
variable SP_begin, SP_end
variable peakdis_multiplier, V_peakpos
variable counter=0
variable npeak=10
make/o/n=(npeak) peakpos=0,  peakpos_alt = 0
make/o/n=1 peakposmax=0
variable Vmax=0
variable fits, fitr, fite
variable peak_flag
peakdis_multiplier = 5

fitr = peakdifthresh*0.4
//make/o/n=(2*fitr+1) fit_dataw


      do
      
      findpeak/B=(boxsize)/M=(peakthresh)/q/R=(Rstart, movielength) dataw
      peak_flag = V_flag
      SP_begin = v_trailingedgeloc
      SP_end = v_trailingedgeloc + peakdis_multiplier*peakdifthresh
      V_peakpos = V_peakloc
      
      if(peak_flag ==0)
         sig1 = 1
         Rstart = V_trailingedgeloc+peakdifthresh
          findpeak/B=(boxsize)/M=(-peakthresh)/q/N/R=(SP_begin,SP_end) dataw     
          
          
          if(V_flag!=0 || (V_peakloc - V_peakpos)>peakdis_multiplier*peakdifthresh)
              peakpos[j] = sig1*V_peakpos   
               j+=1         
           endif    
      endif      
      
      counter+=1
      
     if(counter>50)
        break
     endif   
     
     
      
      while(peak_flag ==0 && Rstart<movielength) //&& numtype(Rstart==0))
      
  
  
    Rstart = 0 
    counter=0
    
    
    do
           findpeak/B=(boxsize)/M=(-peakthresh)/q/N/R=(Rstart,movielength) dataw
           
      peak_flag = V_flag
      SP_begin = v_trailingedgeloc
      SP_end = v_trailingedgeloc + peakdis_multiplier*peakdifthresh
      V_peakpos = V_peakloc           
           
         if(V_flag ==0)
            sig1 = -1
            Rstart = V_trailingedgeloc+peakdifthresh
          findpeak/B=(boxsize)/M=(peakthresh)/q/R=(SP_begin,SP_end) dataw        
          
          if(V_flag!=0 || (V_peakloc - V_peakpos)>peakdis_multiplier*peakdifthresh)
              peakpos[j] = sig1*V_peakpos
               j+=1         
           endif                       
            
                        
         endif       
       

       
      counter+=1
      
     if(counter>50)
        break
     endif   
      
      while(V_flag ==0 && Rstart<movielength) //&& numtype(Rstart==0)) 
 
 for(i=0;i<npeak; i+=1)
    if (peakpos[i]!=0)
      fits = abs(peakpos[i]) - fitr
      fite = abs(peakpos[i]) + fitr
      duplicate/o/R=(fits, fite) dataw dataw_r, fit_dataw
      curvefit/Q/W=2 poly 3 ,dataw_r /D=fit_dataw
      if ( peakpos[i]>0)
      findpeak/q fit_dataw
      peakpos_alt[i] = V_peakloc
      else 
      findpeak/n/q fit_dataw
      peakpos_alt[i] = -V_peakloc 
      endif
      
      if(abs(V_peakval) > vmax)
        vmax = abs(V_peakval)
        peakposmax[0] = peakpos_alt[i]  
      endif   
      
    endif  
      
 endfor
 
 
 
      
      
end      





function findmultipeak_old(dataw, peakthresh,peakdifthresh)
wave dataw
variable peakthresh, peakdifthresh


variable boxsize =100
variable movielength = dimsize(dataw,0)
variable sig1
variable i,j
variable Rstart = 0
variable counter=0
variable npeak=10
make/o/n=(npeak) peakpos=0,  peakpos_alt = 0
make/o/n=1 peakposmax=0
variable Vmax=0
variable fits, fitr, fite
fitr = 400
//make/o/n=(2*fitr+1) fit_dataw


      do
      
      findpeak/M=(peakthresh)/q/R=(Rstart, movielength) dataw
      if(V_flag ==0)
         sig1 = 1
         Rstart = V_trailingedgeloc+peakdifthresh
         peakpos[j] = sig1*V_peakloc       
         j+=1         
      endif      
      
      counter+=1
      
     if(counter>50)
        break
     endif   
      
      while(V_flag ==0) //&& numtype(Rstart==0))
      
  
  
    Rstart = 0 
    counter=0
    do
           findpeak/M=(-peakthresh)/q/N/R=(Rstart,movielength) dataw
         if(V_flag ==0)
            sig1 = -1
            Rstart = V_trailingedgeloc+peakdifthresh
            peakpos[j] = sig1*V_peakloc          
            j+=1   
         endif       
       

       
      counter+=1
      
     if(counter>50)
        break
     endif   
      
      while(V_flag ==0) //&& numtype(Rstart==0)) 
 
 
 
      
      
end      



function corrpeak_old(v_pair, peakdifthresh, boxcarsize)
wave v_pair
//variable peakdifthresh = 50             
variable peakdifthresh        // for short movies
variable boxcarsize       //same as velo_pair function


variable npair = dimsize(v_pair,2)
variable movielength = dimsize(v_pair, 0)
variable peakthresh = 0.008

variable i,j,k,l,m
variable Rend = movielength
variable Rstart = 0 
variable npeak = 10 

variable sig1=0, sig2=0

variable counter = 0 

variable peakmax_loc1, peakmax_loc2

make/o/n=(npair) steppair=0, steploc1_old=0, steploc2_old=0

make/o/n=(movielength) velo_temp1, velo_temp2
make/o/n=(npeak)  peakpos_1, peakpos_2
make/o/n=(npair) peakmaxdif


wave peakposmax, peakpos_alt, peakpos

l=0
m=0
for (i=0; i<npair; i+=1)
      velo_temp1 = v_pair[p][0][i]
      velo_temp2 = v_pair[p][1][i]
      
      peakpos_1 = 0
      peakpos_2 = 0
      
      findmultipeak_old(velo_temp1,peakthresh, peakdifthresh)
//    duplicate/o peakpos peakpos_1
      duplicate/o peakpos peakpos_1
      peakmax_loc1 = peakposmax[0]
      findmultipeak_old(velo_temp2,peakthresh, peakdifthresh)
//      duplicate/o peakpos peakpos_2
      duplicate/o peakpos peakpos_2    
      peakmax_loc2 = peakposmax[0]  

      
      if(peakmax_loc1*peakmax_loc2<0)
         peakmaxdif[m] = peakmax_loc1+peakmax_loc2
         m+=1
      endif    
      
      
      
      for ( j=0; j<npeak;j+=1)
          if ( peakpos_1[j] !=0)
              for ( k=0; k<npeak; k+=1)
                   if ( (abs(peakpos_1[j]+peakpos_2[k])<peakdifthresh)&&(peakpos_1[j]*peakpos_2[k])<0)
                       steppair[l] = i
                       steploc1_old[l] = peakpos_1[j]+sign(peakpos_1[j])*floor(boxcarsize/2)                       
                       steploc2_old[l] = peakpos_2[k]+sign(peakpos_2[k])*floor(boxcarsize/2)
                       l+=1                                 //might result in curves with multiple pairs get repeated
                   endif
               endfor
          endif
      endfor             
                   
                    
      
      
      
      
endfor




redimension/n=(l)   steppair, steploc1_old, steploc2_old
redimension/n=(m) peakmaxdif

end









function calcindcrosscorrelation(dwave1, dwave2, corrlength)
wave dwave1, dwave2
variable corrlength


variable movielength = dimsize(dwave1,0)
variable i, j, lagtime


//make/o/n=(2*corrlength+1) indcorrwave
make/o/n=(2*corrlength/10+1) indcorrwave
variable meanw1, meanw2, stdw1, stdw2
variable posindex, negindex

//wavestats/Q dwave1
//meanw1 = V_avg
//stdw1 = V_sdev

//wavestats/Q dwave2
//meanw2 = V_avg
//stdw2 = V_sdev


for (lagtime =0; lagtime<=corrlength; lagtime+=10)

     make/o/n=(movielength - lagtime) dwave1_temp=0, dwave2_temp=0, corrwave_temp=0, dwave1_temp_ed =0, dwave2_temp_ed = 0
     posindex = (corrlength+lagtime)/10
     negindex = (corrlength-lagtime)/10
     dwave1_temp = dwave1[p+lagtime]// - meanw1
     dwave2_temp = dwave2[p]// - meanw2

    wavestats/Q dwave1_temp
    meanw1 = V_avg
    stdw1 = V_sdev

    wavestats/Q dwave2_temp
    meanw2 = V_avg
    stdw2 = V_sdev
    
    dwave1_temp_ed = dwave1_temp - meanw1
    dwave2_temp_ed = dwave2_temp - meanw2
    
     
     corrwave_temp = (dwave1_temp_ed*dwave2_temp_ed) /(stdw1*stdw2)
     indcorrwave[posindex] = mean(corrwave_temp)
     dwave1_temp = dwave1[p]// - meanw1
     dwave2_temp = dwave2[p+lagtime]// - meanw2
     
     
     wavestats/Q dwave1_temp
     meanw1 = V_avg
     stdw1 = V_sdev

    wavestats/Q dwave2_temp
    meanw2 = V_avg
    stdw2 = V_sdev     

    dwave1_temp_ed = dwave1_temp - meanw1
    dwave2_temp_ed = dwave2_temp - meanw2
     
     corrwave_temp = (dwave1_temp_ed*dwave2_temp_ed) /(stdw1*stdw2)
     indcorrwave[negindex] = mean(corrwave_temp)     


endfor     

setscale x -corrlength, 10,"", indcorrwave




end



function calcallcrosscorrelation(pairwave)
wave pairwave


variable npair = dimsize(pairwave,2)
variable movielength = dimsize(pairwave, 0)
variable corrlength = floor(movielength-1)
variable j=0 

wave selectpair


make/o/n=(movielength) pairwavetemp1, pairwavetemp2, pairwavetemp1_smth, pairwavetemp2_smth
make/o/n=(2*corrlength/10+1,npair) crosscorr=0

wave indcorrwave

variable i 

for(i=0;i<npair; i+=1)
//       if (i==selectpair[j])
       pairwavetemp1 = pairwave[p][0][i]
       pairwavetemp2 = pairwave[p][1][i]
//       boxsmooth 10, pairwavetemp1, pairwavetemp1_smth
//       boxsmooth 10, pairwavetemp2, pairwavetemp2_smth       
       
       calcindcrosscorrelation(pairwavetemp1, pairwavetemp2, corrlength)
       
       crosscorr[][i]= indcorrwave[p]
//       j+=1
//       endif
       
       
endfor

//crosscorravg = crosscorravg/npair

end       


function selectcrosscorr(dwave, refwave)
wave dwave 
wave refwave

variable npair = dimsize(dwave,1)
variable corrlength = dimsize(dwave, 0)
variable i, counter
variable selectthresh= -0.7

counter= 0 
make/o/n=(corrlength,npair) crosscorr_select

for(i=0; i<npair;i+=1)
     if (refwave[i][0] < selectthresh)
         crosscorr_select[][counter] = dwave[p][i]
         counter+=1
     endif
endfor


redimension/n=(corrlength,counter) crosscorr_select
end     


function calccrosscorravg(dwave,methodflag)
wave dwave
variable methodflag    // 0 for simpel averaging, 1 for abs value averaging , 2 for seperating + and -

variable npair = dimsize(dwave, 1)
variable corrlength = dimsize(dwave,0)
variable i, j, V_corr
variable counter_p, counter_n
make/o/n=(corrlength) dwave_temp, crosscorr_avg=0, crosscorr_avgp=0, crosscorr_avgn=0



     if (methodflag ==0)
         for (i=0; i<npair; i+=1)
             dwave_temp = dwave[p][i]
             crosscorr_avg += dwave_temp 
         endfor   
         crosscorr_avg = crosscorr_avg/npair
      elseif(methodflag ==1)      
           for (i=0; i<npair; i+=1)
               dwave_temp = dwave[p][i]
               crosscorr_avg += abs(dwave_temp) 
           endfor   
           crosscorr_avg = crosscorr_avg/npair     
       elseif(methodflag==2)
            for (j=0; j<corrlength; j+=1)
                 counter_p=0
                 counter_n=0
                 for(i=0; i<npair; i+=1) 
                     V_corr = dwave[j][i]
                      if(V_corr>=0)
                         crosscorr_avgp[j]+=V_corr
                         counter_p+=1
                       else
                          crosscorr_avgn[j]+=V_corr  
                          counter_n+=1
                       endif            
                   endfor
                   crosscorr_avgp[j] = crosscorr_avgp[j]/counter_p
                   crosscorr_avgn[j] = crosscorr_avgn[j]/counter_n
               endfor    
     
     
     
       endif     



end



function calc2nddiff(dwave, anafreq)
wave dwave
variable anafreq

duplicate/o dwave dwave_smth


boxsmooth (anafreq/50), dwave, dwave_smth

variable movielength = dimsize(dwave_smth,0)

make/o/n=(movielength-2*anafreq) dwave_2nddiff

variable i

for(i=anafreq;i<movielength-anafreq;i+=1)
     dwave_2nddiff[i-anafreq] = dwave_smth[i+anafreq]+dwave_smth[i-anafreq]-2*dwave_smth[i]
endfor

end     


function find2nddif_peak(dwave, guess, guessrange, anafreq)
wave dwave
variable guess, guessrange, anafreq

variable guess_2
variable fits, fite
variable peakloc
fits = abs(guess) - guessrange - anafreq
fite = abs(guess) + guessrange - anafreq

wavestats/q/r=(fits, fite) dwave
if(guess<0)
   guess_2 = V_minloc
else
    guess_2 = v_maxloc
endif       

fits = guess_2 - guessrange/5
fite = guess_2 + guessrange/5

duplicate/o/R=(fits, fite) dwave dwave_r, fit_dataw
curvefit/Q/W=2 poly 3 ,dwave_r /D=fit_dataw

wavestats/q fit_dataw

if(guess<0)
  peakloc = V_minloc+ anafreq
  
else
  peakloc = V_maxloc + anafreq
  
endif

return(peakloc)    



end


Function step_1(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = a + b*erf((x-c)/d)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = a
	//CurveFitDialog/ w[1] = b
	//CurveFitDialog/ w[2] = c
	//CurveFitDialog/ w[3] = d

	return w[0] + w[1]*erf((x-w[2])/w[3])
End



Function polypeak(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) =  a+ k0*(x-x0)^2
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 3
	//CurveFitDialog/ w[0] = k0
	//CurveFitDialog/ w[1] = x0
	//CurveFitDialog/ w[2] = a

	return w[2] + w[0]*(x-w[1])^2
End



function peakreso(dwave, r1, r2,k0,a)
wave dwave
variable r1 ,r2, k0, a

variable x0 = (r1+r2)/2
variable v_chisqmin, v_chisqthresh

make/o/n=3  peak_coef

peak_coef[0] = k0
peak_coef[1] = x0
peak_coef[2] = a 

funcfit/NTHR=0/q/w=2 polypeak peak_coef v1_1[r1,r2] 
v_chisqmin = v_chisq
v_chisqthresh = v_chisqmin*(1+1/(v_npnts-2))


x0 = peak_coef[1]

variable trialnum = 50
variable trialrange = 10
variable gap = 2*trialrange/trialnum
make/o/n= (trialnum) peakpos_trial, chisq
variable i


for(i=0; i<trialnum; i+=1) 
     peakpos_trial[i] = (i-trialnum/2)*gap + x0
     peak_coef[0] = k0
     peak_coef[1] = peakpos_trial[i]
     peak_coef[2] = a
     funcfit/NTHR=0/H="010"/q/w=2 polypeak peak_coef v1_1[r1,r2]  
     chisq[i] = v_chisq
endfor

findlevels/q/D=peak_findlevel chisq, v_chisqthresh
  
  
end     





function sortratio()
wave v_tau, v_height, peakdif_pair

wave pairstepdisp

duplicate/o pairstepdisp pairstepdisp_select

variable npair = dimsize(v_tau, 0)
variable i
variable height_ratio, tau_ratio
variable ratiothresh = 0.7
variable height_val_1, height_val_2
variable tau_val_1, tau_val_2
variable pairdif_thresh = 30
variable heightthresh= 0.006
variable tauthresh =600

variable corr_tau, corr_height, t_tau, t_height

variable selectcount = 0

duplicate/o v_tau v_tau_select, v_height_select

make/o/n=(2*npair) pair_select_mask=0, tau_select_p1, tau_select_p2, height_select_p1, height_select_p2


for(i=0; i<npair; i+=1)
     height_val_1 = abs(v_height[i][0])
     height_val_2 = abs(v_height[i][1])
     height_ratio = min(height_val_1, height_val_2)/max(height_val_1, height_val_2)
     
     tau_val_1 = abs(v_tau[i][0])
     tau_val_2 = abs(v_tau[i][1])
     tau_ratio = min(tau_val_1,tau_val_2)/max(tau_val_1, tau_val_2)     
     

//       if (abs(peak_dif_pair[i]) < pairdif_thresh)    
//     if((height_ratio>ratiothresh)) && (abs(peak_dif_pair[i]) < pairdif_thresh))

     if((height_ratio>ratiothresh) && (height_val_1> heightthresh) && (height_val_2>heightthresh)&&(tau_val_1<tauthresh)&&(tau_val_2<tauthresh))      
//      if(tau_ratio>ratiothresh)
//        if(tau_val_1<tauthresh && tau_val_2<tauthresh)
        v_tau_select[selectcount][0] = v_tau[i][0]
        v_tau_select[selectcount][1] = v_tau[i][1]
        v_height_select[selectcount][0] = abs(v_height[i][0])
        v_height_select[selectcount][1] = abs(v_height[i][1])
        pairstepdisp_select[][0][selectcount] = pairstepdisp[p][0][i]
        pairstepdisp_select[][1][selectcount] = pairstepdisp[p][1][i]
        pair_select_mask[i] = 1
        pair_select_mask[i+npair]=1
        tau_select_p1[selectcount] = tau_val_1
        tau_select_p2[selectcount] = tau_val_2
        height_select_p1[selectcount] = height_val_1
        height_select_p2[selectcount] = height_val_2
     
        selectcount+=1
     endif
     
endfor



redimension/n=(selectcount,2) v_tau_select, v_height_select
redimension/n=(selectcount) tau_select_p1, tau_select_p2, height_select_p1, height_select_p2
redimension/n=(-1,2,selectcount) pairstepdisp_select

corr_tau = statscorrelation(tau_select_p1, tau_select_p2)
corr_height = statscorrelation(height_select_p1, height_select_p2)


tau_select_p1 = tau_select_p1*0.1
tau_select_p2 = tau_select_p2*0.1

t_tau = corr_tau*sqrt((selectcount-2)/(1-corr_tau^2))
print(corr_tau)
print(corr_tau/t_tau)
print(studentA(t_tau, selectcount-2))



end      



function fsum(ptype)

variable ptype

wave pforce
wave linmaskval_intercept_int_sc

variable npost = dimsize(pforce,0)
variable movielength = dimsize(pforce,1)
variable i, npt=0

make/o/n=(movielength) cellforce_avg=0

for(i=0; i<npost; i+=1)
    if (linmaskval_intercept_int_sc[i] == ptype)
        cellforce_avg += pforce[i][p]
        npt+=1
    endif
endfor

cellforce_avg = cellforce_avg/npt

end        
   
   
   
   
   
function calcmsdfloor(MSDwave,maskwave)

wave MSDwave
wave maskwave
variable npost = dimsize(MSDwave, 1)
variable floorsize =  20
make/o/n=(floorsize) MSDfloortemp
make/o/n=(npost) MSD_floor_bg=0, MSD_floor_lb=0
variable countbg = 0, countlb=0


variable i
for(i=0;i<npost;i+=1)
        MSDfloortemp = MSDwave[p][i]
        if(maskwave[i] == 0)
           MSD_floor_bg[countbg] = log(mean(MSDfloortemp))
           countbg+=1
        elseif(maskwave[i] == 1)
           MSD_floor_lb[countlb] = log(mean(MSDfloortemp))
           countlb+=1
        endif
        
endfor

redimension/n=(countbg) MSD_floor_bg
redimension/n=(countlb) MSD_floor_lb

wavestats/q MSD_floor_bg
make/o/n=(3) MSD_floor_bg_stats, MSD_floor_lb_stats
MSD_floor_bg_stats[0] = V_avg
MSD_floor_bg_stats[1] = V_sdev
MSD_floor_bg_stats[2] = 10^V_avg 


wavestats/q MSD_floor_lb
MSD_floor_lb_stats[0] = V_avg
MSD_floor_lb_stats[1] = V_sdev
MSD_floor_lb_stats[2] = 10^V_avg 

end

   
    
      
      



function detectstaticpost(datawave, maskwave,thresh)
wave datawave, maskwave
variable thresh

variable movielength = dimsize(datawave,2)
variable seglength = 3000
variable npost = dimsize(datawave,0)


duplicate/o maskwave linmaskval_intercept_int_static

variable segnum = floor(movielength/seglength)

variable i, j

wave msdxwave_seg, datawmsd

variable msdlength = dimsize(msdxwave_seg,0)

make/o/n=(seglength) disptempx, disptempx_abs,disptempy, disptempy_abs, disptempran, disptempran_abs

make/o/n=(msdlength, npost) MSD_static
make/o/n=(npost*segnum, 5) Fmax_cp                 //[][0] : Fxmax
                                                             //[][1]: Fymax
                                                             //[][2]:postid
                                                             //[][3]:segid
                                                             //[][4]:maskwave
                                                             
make/o/n=(npost*segnum,2) staticpostid                                                             

variable dispx_max, dispy_max

variable countcp = 0
variable countxstats = 0 
variable countystats = 0

string NBname
sprintf NBname, "statictraceplotalot"

Dowindow/K $NBname
NewNotebook/K=1/F=1/N=$NBname

variable countonpage = 0
variable plotperpage = 7

string pxstr, pystr

wave NNpostid


for ( i = 0; i<npost; i+=1)
   // if ((maskwave[i] !=3) || (maskwave[i] !=0))
      if ( (maskwave[i] ==1) || (maskwave[i] == 7)) 
          linmaskval_intercept_int_static[i] = 1
      for(j=0;j<segnum;j+=1)
          disptempy = datawave[i][1][p+seglength*j]
          disptempy_abs = abs(disptempy)
          dispy_max = wavemax(disptempy_abs)
          
          disptempx = datawave[i][0][p+seglength*j]
          disptempx_abs = abs(disptempx)  
          dispx_max = wavemax(disptempx_abs)
          
          
          Fmax_cp[countcp][0] = dispx_max
          Fmax_cp[countcp][1] = dispy_max
          Fmax_cp[countcp][2] = i
          Fmax_cp[countcp][3] = j
          Fmax_cp[countcp][4] = maskwave[i]
          countcp+=1
                  
//          if( mean(disptempy_abs) < thresh)
            if(wavemax(disptempy_abs)<thresh)
            countxstats+=1
            

            
//            if ( mean(disptempx_abs) < thresh)
              if( wavemax(disptempx_abs)<thresh)
                              
                calcsinglemsdlong_alt(disptempx, disptempy)
                msd_static[][countystats] = datawmsd[p]
                staticpostid[countystats][0] = i     
                staticpostid[countystats][1] = j     
                linmaskval_intercept_int_static[i] = -4      
                countystats+=1

                  if(countonpage == 0)          
                     DoWindow/K Layoutxytrace
                     Newlayout/K=1/N=Layoutxytrace
                  endif

            sprintf pxstr,"Plot%dx",countonpage
            sprintf pystr,"Plot%dy",countonpage
            
            disptrace(disptempx,i,pxstr, countonpage)
            
            disptrace(disptempy,i, pystr, countonpage+plotperpage)


   	    AppendLayoutObject/F=0/R=(41,21 + 100*countonpage,260,120 + 100*countonpage) graph $pxstr
      	    AppendLayoutObject/F=0/R=(281,21 + 100*countonpage,500,120 + 100*countonpage)  graph $pystr  

            countonpage+=1
            
                  if(countonpage == plotperpage)
                       countonpage=0
                       DoWindow/F layoutxytrace
                       Notebook $NBname scaling = {90,90}, picture = {Layoutxytrace, -1,1}
                       Notebook $NBname text = "\r"
                  endif

                                
            endif
            
            

            
         
         endif
     endfor
     
   endif
   
                     if((countonpage!=(plotperpage)) && (i == (npost-1)) )

                       DoWindow/F layoutxytrace
                       Notebook $NBname scaling = {90,90}, picture = {Layoutxytrace, -1,1}
                       Notebook $NBname text = "\r"
                  endif
   
   
   
   if(maskwave[i]<0)
      linmaskval_intercept_int_static[i] = -0.5
   endif
   
   
endfor

redimension/n=(msdlength, countystats) MSD_static
redimension/n=(countcp,5) Fmax_cp
redimension/n=(countystats,2) staticpostid
print(countxstats)
print(countystats)


display Fmax_cp[][0] vs Fmax_cp[][1]
ModifyGraph log(left)=1
ModifyGraph log=1
ModifyGraph mode=3
ModifyGraph zColor(Fmax_cp)={Fmax_cp[*][4],-2,7,Rainbow,0}



end      
      
      

function disptrace(datawave, ipost,  swname, nplot)
wave datawave
variable ipost,nplot
string swname

string dws,s1

sprintf dws, "dw%d", nplot
duplicate/o datawave $dws





    dowindow/K $swname
    
    display/K=1/N=$swname $dws 
    SetAxis left, wavemin($dws),wavemax($dws)

    sprintf s1, "pairid=%d",ipost
    TextBox/C/B=1/N=text0/A=LT s1
    ModifyGraph tick=2,mirror=1,standoff=0
    ModifyGraph margin(left)=50,margin(bottom)=23,margin(top)=2,margin(right)=5;
    ModifyGraph width=165,height=75     

end





function findneighborpost(nxpost, nypost, postid)
variable nxpost, nypost, postid

variable i = floor(postid/nypost)
variable j = mod(postid, nypost)

make/o/n=6 NNpostid
NNpostid[0] = (i-1)*nypost+j
NNpostid[1] = (i+1)*nypost+j
NNpostid[2] = i*nypost+j-1
NNpostid[3] = i*nypost+j+1
NNpostid[4] = (i+2*(mod(j,2)-0.5))*nypost + j-1
NNpostid[5] = (i+2*(mod(j,2)-0.5))*nypost + j+1


end


function calcneighbro(spostid,maskwave, astr)
wave spostid, maskwave, astr

variable npost = dimsize(maskwave,0)
variable NNnum = 6

variable i, j ,k

variable nxpost = astr[3]
variable nypost = astr[4]

wave NNpostid
variable neighborpost, neighbortype


variable ncortical = 0, count_nb=0

make/o/n=(npost) postcount

for( i=0; i<npost; i+=1)
     if (maskwave[i] == 1 || maskwave[i] == 7)

         count_nb=0
         findneighborpost(nxpost, nypost, i)
         for(j=0;j<NNnum;j+=1)
              neighborpost = NNpostid[j]
              neighbortype = maskwave[neighborpost]
             if( (neighbortype>0) && (neighbortype!=3) )
                 count_nb+=1
             endif
             
        endfor
        
        postcount[ncortical] = count_nb
        ncortical+=1
     endif
     
endfor


redimension/n=(ncortical) postcount


count_nb = 0
variable nspost = dimsize(spostid,0)

make/o/n=(nspost) postcount_s=0

for(i=0;i<nspost; i+=1)
        findneighborpost(nxpost, nypost,spostid[i])
        for(j=0;j<NNnum;j+=1)
             neighborpost = NNpostid[j]
             neighbortype = maskwave[neighborpost]
             if(neighbortype!=0)
                postcount_s[i]+=1
             endif
             
             for(k=0;k<nspost;k+=1)
                 if (spostid[k][0] == neighborpost && spostid[k][1] == spostid[i][1])
                    postcount_s[i] -= 1
                    break
                 endif
             endfor
             
        endfor


endfor

end    

function calcforce()
wave fitres_ed_sub
wave linmaskval_intercept_int
variable kspring = 15.7
variable movielength = dimsize(fitres_ed_sub,2)

variable npost= dimsize(fitres_ed_sub,0) 
variable i

make/o/n=(movielength) sforce
make/o/n=(npost) pforce_sd
variable avgforce = 0
variable ncpost = 0 

for(i=0;i<npost;i+=1)
     if(linmaskval_intercept_int[i]>0 && linmaskval_intercept_int[i]!=3)
//      if(linmaskval_intercept_int[i]==1)
         sforce = kspring*sqrt(fitres_ed_sub[i][0][p]^2+fitres_ed_sub[i][1][p]^2)
         avgforce+=mean(sforce)
         wavestats/q sforce
         pforce_sd[i] = V_sdev
         ncpost +=1
     endif
     
endfor

avgforce=avgforce/ncpost

print(avgforce)

end    



function calcforcesd()

wave pforce

variable npost = dimsize(pforce,0)
variable ipost
variable movielength = dimsize(pforce,1)

make/o/n=(movielength) sforce
make/o/n=(npost) pforce_sd


for(ipost=0;ipost<npost;ipost+=1)
    sforce = pforce[ipost][p]
    wavestats/q sforce
    pforce_sd[ipost] = v_sdev
endfor

end
    
    

function splitedmsd(datawave,maskwave)
wave datawave, maskwave


newpath/c/o lmaxfolder

variable postid

variable seglength = 3000

variable postnum = dimsize(datawave, 0)
variable movielength = dimsize(datawave,2)
variable nseg = floor(movielength/seglength)

variable i


make/o/n=(seglength) disptempx, disptempy

wave datawmsd , msdxwave_seg
variable msdlength = dimsize(msdxwave_seg, 0)
make/o/n=(msdlength, nseg) msd_sp_seg
make/o/n=(nseg) lmsd_single, lmax_single

variable tau = 200      //20s in lag time
variable spreadpara_lmsd, spreadpara_lmax

make/o/n=(seglength-tau)  wmsd

make/o/n=(postnum) lmsd_single_spread_lb, lmax_single_spread_lb, lmsd_single_spread_hb, lmax_single_spread_hb

variable nlb = 0
variable nhb = 0

for(postid=0;postid<postnum;postid+=1)
if(maskwave[postid] == 1 || maskwave[postid] == 6)
for (i=0;i<nseg;i+=1)
     disptempx = datawave[postid][0][p+i*seglength]*125
     disptempy = datawave[postid][1][p+i*seglength]*125
     duplicate/o/r=[0,seglength-tau-1] disptempx wmsdx1
     duplicate/o/r=[0,seglength-tau-1] disptempy wmsdy1
     duplicate/o/r=[tau,seglength-1] disptempx wmsdx2
     duplicate/o/r=[tau,seglength-1] disptempy wmsdy2
     wmsd =( wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
     lmax_single[i] = sqrt(wavemax(wmsd))
     lmsd_single[i] = sqrt(mean(wmsd))

     

     
//     lmax_single[i] = sqrt(calcsinglemsdlong_alt(disptempx, disptempy))
//     msd_sp_seg[][i] = datawmsd[p]
//     lmsd_single[i] = sqrt(datawmsd[189])            //msd at tau = 100s     
endfor     

     wavestats/q lmsd_single
     spreadpara_lmsd = V_sdev/V_avg
     wavestats/q lmax_single
     spreadpara_lmax = V_sdev/V_avg


     if(maskwave[postid] == 1)
        lmsd_single_spread_lb[nlb] = spreadpara_lmsd
        lmax_single_spread_lb[nlb] = spreadpara_lmax
        nlb+=1
     endif
     
     
     if(maskwave[postid] == 6)
        lmsd_single_spread_hb[nhb] = spreadpara_lmsd
        lmax_single_spread_hb[nhb] = spreadpara_lmax     
        nhb+=1
     endif
endif

endfor


redimension/n=(nlb) lmsd_single_spread_lb,lmax_single_spread_lb
redimension/n=(nhb) lmsd_single_spread_hb, lmax_single_spread_hb

string lstarstr_lb = "lmsd_single_spread_lb.txt"
string lmaxstr_lb = "lmax_single_spread_lb.txt"
string lstarstr_hb = "lmsd_single_spread_hb.txt"
string lmaxstr_hb = "lmax_single_spread_hb.txt"


save/G/o/A=2/p=lmaxfolder lmsd_single_spread_lb as lstarstr_lb
save/G/o/A=2/p=lmaxfolder lmax_single_spread_lb as lmaxstr_lb
save/G/o/A=2/p=lmaxfolder lmsd_single_spread_hb as lstarstr_hb
save/G/o/A=2/p=lmaxfolder lmax_single_spread_hb as lmaxstr_hb


end    





function calcsinglemsdlong_alt(dwtempx1, dwtempy1)
wave dwtempx1, dwtempy1

variable i=0, j=0 
variable L
variable len = dimsize(dwtempx1, 0)
variable DX, DY
L = floor(len/1)
make/N=(L)/O   datawMSD=0, msdxwave_seg = 0
variable v_max


     
for (i = 1; i<100; i+=1)
     duplicate/o/r = [0,len-i-1] dwtempx1, wmsdx1
     duplicate/o/r = [0,len-i-1] dwtempy1, wmsdy1
     duplicate/o/r = [i, len-1] dwtempx1, wmsdx2
     duplicate/o/r = [i, len-1] dwtempy1, wmsdy2     
     make/o/n=(len-i)  wmsd
     wmsd = (wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
     datawmsd[j] = mean(wmsd)
     msdxwave_seg[j] = i
     j=j+1
      endfor
      
for(i=100; i<L;i+=10)
     duplicate/o/r = [0,len-i-1] dwtempx1, wmsdx1
     duplicate/o/r = [0,len-i-1] dwtempy1, wmsdy1
     duplicate/o/r = [i, len-1] dwtempx1, wmsdx2
     duplicate/o/r = [i, len-1] dwtempy1, wmsdy2     
     make/o/n=(len-i)  wmsd
     wmsd = (wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
     datawmsd[j] = mean(wmsd) 
     if (j==189)
         v_max= wavemax(wmsd)
     endif
     msdxwave_seg[j] = i+1
     j=j+1
endfor


        
redimension/N=(j)      datawmsd , msdxwave_seg
msdxwave_seg = msdxwave_seg/10

//return v_max

end.

    
    
    
    
function splitedmsd_single(datawave,postid,maskwave)
wave datawave, maskwave
variable postid



variable seglength = 3000

variable postnum = dimsize(datawave, 0)
variable movielength = dimsize(datawave,2)
variable nseg = floor(movielength/seglength)

variable i


make/o/n=(seglength) disptempx, disptempy

wave datawmsd , msdxwave_seg
variable msdlength = dimsize(msdxwave_seg, 0)
make/o/n=(msdlength, nseg) msd_sp_seg
make/o/n=(nseg) lmsd_single, lmax_single

variable tau =1000        //100s in lag time
variable spreadpara_lmsd, spreadpara_lmax

make/o/n=(seglength-tau)  wmsd

make/o/n=(postnum) lmsd_single_spread_lb, lmax_single_spread_lb, lmsd_single_spread_hb, lmax_single_spread_hb

variable nlb = 0
variable nhb = 0

//for(postid=0;postid<postnum;postid+=1)
//if(maskwave[postid] == 1 || maskwave[postid] == 6)
for (i=0;i<nseg;i+=1)
     disptempx = datawave[postid][0][p+i*seglength]*125
     disptempy = datawave[postid][1][p+i*seglength]*125
     duplicate/o/r=[0,seglength-tau-1] disptempx wmsdx1
     duplicate/o/r=[0,seglength-tau-1] disptempy wmsdy1
     duplicate/o/r=[tau,seglength-1] disptempx wmsdx2
     duplicate/o/r=[tau,seglength-1] disptempy wmsdy2
     wmsd =( wmsdx2[p] - wmsdx1[p] )^2 +(wmsdy2[p] - wmsdy1[p])^2
     lmax_single[i] = sqrt(wavemax(wmsd))
     lmsd_single[i] = sqrt(mean(wmsd))

     

     
     lmax_single[i] = sqrt(calcsinglemsdlong_alt(disptempx, disptempy))
     msd_sp_seg[][i] = datawmsd[p]
//     lmsd_single[i] = sqrt(datawmsd[189])            //msd at tau = 100s     
endfor     
    
      wavestats/q lmsd_single
     print ( V_sdev/V_avg ) 
    
    
    
 end   
 
 
 
 
 function calccagesize(dataw, maskw)
 
 wave dataw, maskw
 
 variable npost = dimsize(dataw,0)
 variable movielength = dimsize(dataw,2)
 
 variable interval = 100
 
 variable ntime = floor(movielength/interval)
 
 make/o/n=(ntime*npost) cagesize_all, cagesize_cortical, cagesize_all_x, cagesize_all_y,cagesize_c_x, cagesize_c_y, cagesize_cs_x, cagesize_cs_y
 
 variable i,j
 variable icell = 0
 
 for(i=0;i<npost;i+=1)
      if (maskw[i] >0 && maskw[i]!=3)
//        if (maskw[i] ==1)
          for(j=0;j<ntime;j+=1)
               cagesize_all[icell*ntime+j] = sqrt(dataw[i][0][j*interval]^2+dataw[i][1][j*interval]^2)
               cagesize_all_x[icell*ntime+j] = dataw[i][0][j*interval]
               cagesize_all_y[icell*ntime+j] = dataw[i][1][j*interval]
          endfor
          icell+=1
       endif
endfor     


redimension/n=(ntime*icell) cagesize_all, cagesize_all_x, cagesize_all_y



 icell = 0
 
 for(i=0;i<npost;i+=1)
//      if (maskw[i] >0 && maskw[i]!=3)
        if (maskw[i] ==1)
          for(j=0;j<ntime;j+=1)
               cagesize_cortical[icell*ntime+j] = sqrt(dataw[i][0][j*interval]^2+dataw[i][1][j*interval]^2)
               cagesize_c_x[icell*ntime+j] = dataw[i][0][j*interval]
               cagesize_c_y[icell*ntime+j] = dataw[i][1][j*interval]
               cagesize_cs_x[icell] = dataw[i][0][9000]
               cagesize_cs_y[icell] = dataw[i][1][9000]
          endfor
          icell+=1
       endif
endfor     


redimension/n=(ntime*icell) cagesize_cortical, cagesize_c_x,cagesize_c_y
redimension/n=(icell) cagesize_cs_x, cagesize_cs_y




end




function makeNNandCCunit()   
wave fitres_ed_sub
wave pairdispwave_seg,pairdispwave_seg_t,disp_corr_seg
wave vpair
wave steppair, steploc1, steploc2


//make near neighbor pairs
bondcorr_seg(1,900,fitres_ed_sub)
//calculate smoothed velocity on post pairs
velo_pair(pairdispwave_seg,50)
//correlate peaks between pairs
corrpeak(vpair, 50,50)      

//create near neighbor contracting pairs
selectacpair(disp_corr_seg,pairdispwave_seg,pairdispwave_seg_t,0.7)

matchseg_short_alt()


end



function makeNNandCCunit_ran()   
wave fitres_ed_sub
wave pairdispwave_seg,pairdispwave_seg_t,disp_corr_seg
wave vpair
wave steppair, steploc1, steploc2


//make near neighbor pairs
bondcorr_seg_ran(1,900,fitres_ed_sub)
//calculate smoothed velocity on post pairs
velo_pair(pairdispwave_seg,50)
//correlate peaks between pairs
corrpeak(vpair, 50,50)      

//create near neighbor contracting pairs
selectacpair(disp_corr_ran,pairdispwave_seg,pairdispwave_seg_t,0.8)

matchseg_short_alt()


end



   
function makeNNunit()  
// 9/7/21 compilation bug fixed by DHR 
wave fitres_ed
wave zeropos
wave linmaskval
wave pairdispwave_seg,pairdispwave_seg_t,disp_corr_seg
wave vpair
wave steppair, steploc1, steploc2

variable npost = dimsize(linmaskval , 0)
variable movielength = dimsize(fitres_ed_sub, 2)
variable seglength = 1000
variable nseg = floor(movielength/seglength)
make/o/n=(npost, nseg) linmaskval_seg_int
//duplicate/o linmaskval linmaskval_intercept_int
wave linmaskval_intercept_int

duplicate/o fitres_ed fitres_ed_sub

variable pratio = 0.125
fitres_ed_sub = (fitres_ed[p][q][r] - zeropos[p][q][0])*pratio

variable i , j
for(i=0; i<nseg;i+=1)
    linmaskval_seg_int[][i] = linmaskval_intercept_int[p]
endfor



//make near neighbor pairs
bondcorr_seg(1,seglength,fitres_ed_sub)
//calculate smoothed velocity on post pairs
velo_pair(pairdispwave_seg,100)
//correlate peaks between pairs
corrpeak(vpair, 50,100)  

//create near neighbor contracting pairs
selectacpair(disp_corr_seg,pairdispwave_seg,pairdispwave_seg_t,0.6)
matchseg_short_alt()

bondcorr_edge(fitres_ed_sub)

wave synpairid, bond_edge, synpairid_edge
variable count = 0

make/o/n=(dimsize(synpairid,0)) synpairid_edge

for (i = 0; i<dimsize(synpairid,0); i+=1)
      for(j = 0; j<dimsize(bond_edge,0); j+=1)      //   9/7/21 Was j < bond_edge - didn't compile in Igor 9
            if(synpairid[i] == bond_edge[j])
                synpairid_edge[count] = synpairid[i] 
                count+=1
            endif
      endfor
endfor


//bondcorr(movielength,fitres_ed_sub)

end








function bondcorr_edge(datawave)

wave datawave  

wave disp_corr_seg

wave linmaskval_intercept_int_sc, bondid

variable posttype      // 1 for low traction force posts, 6 for high traction force posts
   
variable npair = dimsize(bondid, 0)
variable nframe = dimsize(datawave, 2)

variable ipair, postA, postB, iframe,count=0
make/o/n=(nframe) pairwave1p, pairwave2p

wave pairdispwave


make/o/n=(npair) disp_corr_edge=0, bond_edge=0

for (ipair = 0; ipair<npair; ipair+=1)
      postA = bondid[ipair][0]
      postB = bondid[ipair][1]
      if(linmaskval_intercept_int_sc[postA] != linmaskval_intercept_int_sc[postB])
           disp_corr_edge[count] = disp_corr_seg[ipair][0]        
           bond_edge[count] = ipair  
           count+=1 
       endif
endfor
           


//disp_dist_avg_edge = disp_dist_avg_edge/count

redimension/n=(count) disp_corr_edge, bond_edge


end








//plot pairs at the bonuadry between cells

function bonddist_edge_plotalot(datawave, bondpairwave)

wave datawave
wave bondpairwave

wave linmaskval_intercept_int_sc, bondid

variable ipair, iseg
variable npair = dimsize(datawave,0)
variable nseg = dimsize(datawave,1)
variable wavelength = dimsize(bondpairwave, 0)


variable corval


 variable ipcount = 0, plotsoncurrentpage = 0
 variable plotsperpage = 7
 wave msdxwave_diff
 
 variable postA, postB
 
 
 string  NBname, p1str, p2str
 
            sprintf NBname, "postpairplotalot"
            DoWindow/K $NBname
            NewNotebook/K=1/F=1/N=$NBname

make/o/n=(wavelength)  pairwave1, pairwave2
make/o/n=(wavelength, 2*plotsperpage) colorwave        

make/o/n= (npair*nseg) disp_corr_edge   

make/o/n=(wavelength/2+1, npair)  fftphase_diff_edge
            
for ( ipair = 0; ipair<npair; ipair+=1)             

          postA = bondid[ipair][0]
          postB = bondid[ipair][1]
           
  if(linmaskval_intercept_int_sc[postA] != linmaskval_intercept_int_sc[postB])

          pairwave1 = bondpairwave[p][0][ipair]
          pairwave2 = bondpairwave[p][1][ipair]
          
          FFT/OUT=5/DEST=pairwave1_FFT pairwave1
          FFT/OUT=5/DEST=pairwave2_FFT pairwave2
          
          duplicate/o pairwave1_FFT pairwave_diff_FFT
          
          pairwave_diff_FFT = mod(pairwave1_FFT - pairwave2_FFT, pi)
          
          fftphase_diff_edge[][ipcount] = pairwave_diff_FFT[p]
          
         for(iseg = 0; iseg<nseg; iseg+=1)  
               disp_corr_edge[ipcount*nseg+iseg] = datawave[ipair][iseg]              
           endfor    

          
          SetScale/P x 0,0.01,"",pairwave1, pairwave2
            
         if (mod(ipcount, plotsperpage) == 0)
            plotsoncurrentpage = 0
            DoWindow/K LayoutXY
            NewLayout/K=1/N=LayoutXY
        endif    


          corval = datawave[ipair][0]
          colorwave[][plotsoncurrentpage] = corval
          colorwave[][plotsoncurrentpage+plotsperpage] = corval
      
    
            sprintf p1str,"Plot%dp1",plotsoncurrentpage
            sprintf p2str,"Plot%dp2",plotsoncurrentpage
            
            dispwavepair(pairwave1,ipair,p1str, plotsoncurrentpage,corval, colorwave)
            
            dispwavepair(pairwave2,ipair, p2str, plotsoncurrentpage+plotsperpage,corval,colorwave)


   	    AppendLayoutObject/F=0/R=(41,21 + 100*plotsoncurrentpage,260,120 + 100*plotsoncurrentpage) graph $p1str
      	    AppendLayoutObject/F=0/R=(281,21 + 100*plotsoncurrentpage,500,120 + 100*plotsoncurrentpage)  graph $p2str    
              ipcount = ipcount+1
             plotsoncurrentpage =  plotsoncurrentpage + 1
                if ((plotsoncurrentpage == 7) || (ipair == npair))
                  DoWindow/F LayoutXY
             //if (printflag == 1)
                   //PrintLayout LayoutXY  //Modified Jul 27
                       Notebook $NBname scaling = {90,90}, picture = {LayoutXY, -1,1}
                       Notebook $NBname text = "\r"
               endif        
               

endif

endfor               
               
redimension/N=(ipcount*nseg) disp_corr_edge        
redimension/N=(wavelength/2+1, ipcount) fftphase_diff_edge
               
               
SaveNotebook/P=ToPractice $NBname


end


function bondid_unique(bondid)
wave bondid
variable num_pair = dimsize(bondid, 0)

variable npair_unique = 0 
variable i 

make/o/n=2 bondid_prev0, bondid_prev1, bondid_prev2
for (i = 0;i<num_pair;i+=1)
	if (i != 0)
		if ( (bondid[i][0] != bondid_prev0[0] || bondid[i][1] != bondid_prev0[1]) &&(bondid[i][0] != bondid_prev1[0] || bondid[i][1] != bondid_prev1[1])&&(bondid[i][0] != bondid_prev2[0] || bondid[i][1] != bondid_prev2[1]))
			npair_unique +=1
		endif
	endif
	bondid_prev0[0] = bondid[i][0]
	bondid_prev0[1] = bondid[i][1]
	
	if (i >0)
		bondid_prev1[0] = bondid[i-1][0]
		bondid_prev1[1] = bondid[i-1][1]
	endif
	
	if (i>1)
		bondid_prev2[0] = bondid[i-2][0]
		bondid_prev2[1] = bondid[i-2][1]
	endif	
	
endfor

print (npair_unique)

end



function make_cell_coord()
wave linmaskval_intercept_int, linmaskval_seg_int
wave fitres_ed_s
wave lengthandangle
wave pforce_avg

variable npost = dimsize(fitres_ed_s,0)
variable i
variable n_cellpost = 0


make/o/n=(npost ,2) cell_coord, force_vec

for (i=0;i<npost;i+=1)
	if (linmaskval_intercept_int[i] !=0 && linmaskval_intercept_int[i] !=3&& linmaskval_intercept_int[i] !=-2)
//	if (linmaskval_seg_int[i][19] !=0 && linmaskval_seg_int[i][19] !=3)
		cell_coord[n_cellpost][0] = fitres_ed_s[i][0]
		cell_coord[n_cellpost][1] = fitres_ed_s[i][1]
		force_vec[n_cellpost][0] = pforce_avg[i][0]
		force_vec[n_cellpost][1] = lengthandangle[i][1]
		n_cellpost+=1
	endif
endfor

redimension/n=(n_cellpost,2) cell_coord, force_vec

end



function disect_pforce()
wave pforce_avg
wave linmaskval_intercept_int

variable i
variable n_post = dimsize(pforce_avg,0)


variable n_lb = 0
variable n_hb = 0
variable n_cell = 0

variable pforce_avg_lb = 0 
variable pforce_avg_hb = 0
variable pforce_avg_whole = 0

for (i=0;i<n_post;i+=1)
	if (linmaskval_intercept_int[i] == 1)
		n_lb+=1
		n_cell+=1
		pforce_avg_lb+=pforce_avg[i]
		pforce_avg_whole +=pforce_avg[i]
	elseif (linmaskval_intercept_int[i] == 6)
		n_hb+=1
		n_cell+=1
		pforce_avg_hb+= pforce_avg[i]
		pforce_avg_whole +=pforce_avg[i]		
	elseif(linmaskval_intercept_int[i] == 7)
		n_cell+=1
		pforce_avg_whole +=pforce_avg[i]		
	endif
endfor

pforce_avg_lb = pforce_avg_lb/n_lb
pforce_avg_hb = pforce_avg_hb/n_hb


//print(n_lb)
//print(pforce_avg_lb)
//print(n_hb)
//print(pforce_avg_hb)
print(n_cell)
print(pforce_avg_whole)
print(pforce_avg_whole/n_cell)
end



