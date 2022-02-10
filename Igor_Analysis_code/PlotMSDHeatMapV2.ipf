#pragma TextEncoding = "UTF-8"
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later
#pragma rtGlobals=1		// Use modern global access method and strict wave access.




//12/29/21 Yu: PlotOneFrameHeat03 generate hexgon heatmap plot, MSD scatter plot generate scatter plot in PNAS fig2


// June 9, 2017
// version 02:
// Works from filename roots
// Files:   rootx,rooty,roota,rootmsd   a is exponent, msd is MSD at fixed lag time, these = 0 for bkgnd, -1 for initially engaged only, -2 for finally engaged only
// Edit these by hand to remove any bad points - need routine where these can be plotted to find bad points (pull from old code)
// run sort to make:
//     rootxc,rootyc,rootac, rootmsdc - for cell posts that will be shown in color
//     rootxb, rootyb - "background" posts' positions
//     rootxi, rootyi - positions of initially- engaged only posts 
//     rootxf, rootyf - positions of finally-engaged only posts

// Code to make a hexagon-based heat map 
// Used for Yu Shi data  2017
// Code of origin:  PlotStack03.ipf  PlotOneFrameStack() from StretchAnalysisCode from Craig Copeland Thesis

// removes data of form (x,y,0) from x,y,z waves
 function SortThePosts(croot)
    string croot  // root of cell name
    
    
    wave fitres_ed
    wave msd_slope_mean
    wave msd_mag_para
    wave linmaskval_intercept_int
    
    variable i
    
    // make name strings for the input waves
    variable npost = dimsize(fitres_ed,0)
    make/o/n=(npost) xin, yin
	xin = fitres_ed[p][0][0]
	yin = fitres_ed[p][1][0]
	
	duplicate/o xin $(croot+"x")
	duplicate/o yin $(croot+"y")
	
    WAVE ain = msd_slope_mean
    WAVE msdin = msd_mag_para
    WAVE forcein = pforce_plot
    WAVE maskin = linmaskval_intercept_int

    // make the "cell" waves
    duplicate/O xin $(croot+"xc")
    WAVE xc = $(croot+"xc")
    duplicate/O yin $(croot+"yc")
    WAVE yc = $(croot+"yc")
    duplicate/O ain $(croot+"ac")
    WAVE ac = $(croot+"ac")
    duplicate/O msdin $(croot+"msdc")
    WAVE msdc = $(croot+"msdc")
    duplicate/O forcein $(croot+"forcec")
    WAVE forcec = $(croot+"forcec")
    duplicate/o maskin $(croot+"maskc")
    wave maskc = $(croot+"maskc")
       

    
    i = 0      
    do
        if (maskc[i] <=0 || maskc[i] == 3)   // keep only posts with positive slope values
            DeletePoints i,1,xc,yc,ac,msdc,forcec, maskc
        else
            i +=1
        endif
    while (i < dimsize(xc,0))
    
    // make the background position waves - z wave will be a dummy
    duplicate/O xin $(croot+"xb")
    WAVE xb = $(croot+"xb")
    duplicate/O yin $(croot+"yb")
    WAVE yb = $(croot+"yb")
    duplicate/O maskin $(croot+"dummy")
    WAVE maskb = $(croot+"dummy")

    i = 0
    do
        if (maskb[i] !=0)   // keep only posts with zero "slope"
            DeletePoints i,1,xb,yb,maskb
        else
            i +=1
        endif
    while (i < dimsize(maskb,0))
 
    // make the  initially engaged position waves - z wave will be a dummy
    duplicate/O xin $(croot+"xi")
    WAVE xi = $(croot+"xi")
    duplicate/O yin $(croot+"yi")
    WAVE yi = $(croot+"yi")
    duplicate/O maskin $(croot+"dummy")
    WAVE maski = $(croot+"dummy")

    i = 0
    do
        if (maski[i] != -1)   // keep only posts with -1 initially engaged
            DeletePoints i,1,xi,yi,maski
        else
            i +=1
        endif
    while (i < dimsize(xi,0))
 
    // make the finally engaged position waves - z wave will be a dummy
    duplicate/O xin $(croot+"xf")
    WAVE xf = $(croot+"xf")
    duplicate/O yin $(croot+"yf")
    WAVE yf = $(croot+"yf")
    duplicate/O maskin $(croot+"dummy")
    WAVE maskf = $(croot+"dummy")

    i = 0
    do
        if (maskf[i] != -2)   // keep only posts with -2 - not initially engaged
            DeletePoints i,1,xf,yf,maskf
        else
            i +=1
        endif
    while (i < dimsize(xf,0))
 
   
end

//--------------------PlotOneFrameHeat03() ------------------
// 1/27/14
//  plots one frame of some scalar already calculated.  linmaskwaves must be peviously made
// Renamed and modified 5/30/17
// 6/12/17 re-done to work from file name
//
// Expects there to be
//     rootxc,rootyc,rootac, rootmsdc - for cell posts that will be shown in color
//     rootxb, rootyb - "background" posts' positions
//     rootxi, rootyi - positions of initially- engaged only posts 
//     rootxf, rootyf - positions of finally-engaged only posts

function PlotOneFrameHeat03(PlotNameStr,croot,aormsdflag,scaleflag,ticflag,bkgndflag,initfinalflag,markersize)
    string PlotNameStr   // Name of window for this plot 
    string  croot    // root of wave names
    variable  aormsdflag  //  =0 for exponent (alpha) = 1 for MSD at some fixed time  = 2 for force vector map
    variable scaleflag  // = 0 to autoscale this frame; = 1 to scale from zero to max(zwave)
    variable ticflag  // = 0 for no tics, = 1  no tics
    variable bkgndflag // = 0 for no bkgnd posts, = 1 for bkgndposts
    variable initfinalflag  // = 1 to plot initial only and final only cell posts  (color/shape currently fixed to gray in code  
    variable markersize // size of hexagons
    
    
    variable xmin //  these 4 set the plotting range (set xmin =-1 for default)
    variable ymin
    variable xmax
    variable ymax
    
    variable imageXsize = 1024 // default size of image  (1024x1024 for fastcamera, 1350x1000 for Roper)
    variable imageYsize = 1024
    
    variable zmin1 = 0.65     // these seem to be good alpha values for Cells 1-3   - reset from 0.5 to 0.65 6/12/17
    variable zmax1 = 1.85
    variable zminMSD = 7e-6  // seem to be good MSD values for Cells1 -3
    variable zmaxMSD = 0.0075
    variable i,j
    string TBstr,backgroundstr
    
    
    
    make/o/n=(8,3) directcolorwave
	directcolorwave[1][1] = 65280
	directcolorwave[6][2] = 65280
	directcolorwave[7][0] = 57000
	directcolorwave[7][1] = 57000	
	
//    string nameCellw = NameOfWave(linmaskstacky)   // for use with modifygraph
//    string nameBkdw = NameOfWave(by)
  
    // make pointers to the various waves
    
    //autoscale image
    wave xw = $(croot+"x")
    wave yw = $(croot+"y")
    
    variable npost = dimsize(xw,0)
    variable ipost
    variable xmin_temp,ymin_temp
    xmin = 1E5
    ymin = 1E5
    for (ipost=0;ipost<npost;ipost+=1) 
    	xmin_temp = xw[ipost]
    	ymin_temp = yw[ipost]
    	if (xmin > xmin_temp && xmin_temp >0)
    		xmin = xmin_temp
    	endif

    	if (ymin > ymin_temp && ymin_temp >0)
    		ymin = ymin_temp
    	endif    
    	
    endfor	
    xmin = xmin - 20
    xmax = wavemax(xw) + 20
    ymin = ymin - 20
    ymax = wavemax(yw) + 20
    
    make/o/n=(1,2) force_pos, force_scale
    force_scale[0][0] = 20
    force_pos[0][0] = xmin+20
    force_pos[0][1] = ymin+20
    
    
    // cell waves
    WAVE xc = $(croot+"xc")
    WAVE yc = $(croot+"yc")
    if (aormsdflag == 0)  //want to plot exponent
        WAVE zwave = $(croot+"ac")
    elseif (aormsdflag == 1) // want to plot MSD magnitude
        WAVE zwave = $(croot+"msdc")
    elseif (aormsdflag == 2) // want to plot force vector map
        wave zwave = $(croot+"maskc")
    endif
    // background waves if desired
    if (bkgndflag == 1)  
        WAVE xb = $(croot+"xb")
        WAVE yb = $(croot+"yb")
    endif
    // initially engaged or finally engaged waves if desired
    if (initfinalflag == 1)
        WAVE xi = $(croot+"xi")
        WAVE yi = $(croot+"yi")
        WAVE xf = $(croot+"xf")
        WAVE yf = $(croot+"yf")
    endif

//  make a dummy image
    sprintf backgroundstr, "%sbkd",PlotNameStr
    if (xmin >=0) // use user-given image dimensions
        make/O/N=(xmax-xmin+1,ymax-ymin+1)  $backgroundstr
        SetScale/I x xmin,xmax,"", $backgroundstr
        SetScale/I y ymin,ymax,"", $backgroundstr
     else // default
        make/O/N=(imageXsize,ImageYsize)  $backgroundstr
    endif
    WAVE bkdwave = $backgroundstr

    dowindow/k $PlotNameStr

    newimage/k=1/N=$PlotNameStr bkdwave
    if (xmin >=0)
        SetAxis top xmin - 0.5, xmax + 0.5
        SetAxis/R left ymax+0.5,ymin -0.5
     else
         SetAxis top  -0.5, imageXsize + 0.5
        SetAxis/R left imageYsize+0.5, -0.5
     endif
     
// append the wave we care about
  appendtograph/t yc vs xc
  
      
    if (scaleflag == 0) // autoscale this frame
        ModifyGraph zcolor($(croot+"yc")) = {zwave,*,*,rainbow,0}      
    else //  scale using values given in code - should work for multiple cells - currently only works for exponent
         if  (aormsdflag == 0)  //are  plotting exponent
             ModifyGraph zcolor($(croot+"yc")) = {zwave,zmin1,zmax1,rainbow,0}
         elseif (aormsdflag == 1)  // are plotting MSD mag - need log scale
             ModifyGraph logZColor($(croot+"yc"))=1
             ModifyGraph zcolor($(croot+"yc")) = {zwave,zminMSD,zmaxMSD,rainbow,0}     
         elseif (aormsdflag == 2) // plot force vector map
              ModifyGraph zColor($(croot+"yc"))={zwave,*,*,cindexRGB,0,directcolorwave}
         endif
   		
    endif   


    ModifyGraph mode=3,marker=55,msize=markersize
    if (aormsdflag == 2)
    	appendtograph/t yc/tn=forcetrace vs xc
    	modifygraph mode=3
   	 ModifyGraph mrkThick(forcetrace)=0, arrowMarker(forcetrace)={$(croot+"forcec"),2,4,2,0}
    endif



    if (ticflag == 0)
        ModifyGraph tick=3,noLabel=2
    endif
    
    if (bkgndflag != 0)   // if plotting backgrond posts
        appendtograph/t yb vs xb
        ModifyGraph mode=3 
//      ModifyGraph marker($nameBkdw)=55,msize($nameBkdw)=markersize,rgb($nameBkdw)=(56797,56797,56797)  // gray hexagons
//       ModifyGraph marker($nameBkdw)=55,msize($nameBkdw)=markersize,rgb($nameBkdw)=(49151,53155,65535) // blue hexagons
        ModifyGraph marker( $(croot+"yb"))=19,msize($(croot+"yb"))=markersize/2,rgb($(croot+"yb"))=(56797,56797,56797)  // gray circles
//       ModifyGraph marker($nameBkdw)=19,msize($nameBkdw)=8,rgb($nameBkdw)=(49151,53155,65535)  // blue circles
    endif   

    if (initfinalflag != 0)   // if plotting initially and finally engaged posts
        appendtograph/t yi vs xi
        ModifyGraph mode=3 
        ModifyGraph marker( $(croot+"yi"))=55,msize($(croot+"yi"))=markersize,rgb($(croot+"yi"))=(56797,56797,56797)  // gray hexagons
        appendtograph/t yf vs xf
        ModifyGraph mode=3 
        ModifyGraph marker( $(croot+"yf"))=55,msize($(croot+"yf"))=markersize,rgb($(croot+"yf"))=(56797,56797,56797)  // gray hexagons
    endif   

  ModifyImage $backgroundstr explicit=1
//  ModifyImage $backgroundstr eval={0,56797,56797,56797} //gray background fill
  ModifyImage $backgroundstr eval={0,49151,53155,65535} // light blue background fill
  ModifyGraph width={perUnit,0.5,top},height={perUnit,0.5,left} 
    doupdate
  
  if (aormsdflag == 0)  
  	ColorScale/C/N=text0/F=0/M/A=RC trace=$(croot+"yc");DelayUpdate 
  	ColorScale/C/N=text0 "\\Z16\\F'GreekS'a"
  elseif (aormsdflag ==1)
  	ColorScale/C/N=text0/F=0/M/A=RC trace=$(croot+"yc"),log=1,lblMargin=30;DelayUpdate 
  	ColorScale/C/N=text0 "\\Z16 MSD (\\F'GreekS't=10\\F'Arial's)"  	
  elseif (aormsdflag == 2)
	appendtograph/T force_pos[][1] vs force_pos[][0]
	ModifyGraph mode=3
	ModifyGraph arrowMarker(force_pos)=0
	ModifyGraph mrkThick(force_pos)=0,arrowMarker(force_pos)={force_scale,2,4,2,0}  
	TextBox/C/N=text0/F=0/M/A=LT "20 nN"
	TextBox/C/N=text0/B=1/X=(20/xmax*100)/Y=(20/ymax*100+2)
	Legend/C/N=text1/J/F=0/B=1/M/A=LT "\\K(0,65280,0)\\W1019 \\K(0,0,0)\\F'Arial'Cortical\r\\K(0,0,65280)\\W1019 \\K(0,0,0)\\F'Arial'Stress fiber";DelayUpdate
	AppendText "\\K(57000,57000,0)\\W1019 \\K(0,0,0)\\F'Arial'Not Assigned"
	Legend/C/N=text1/J/A=RB/X=5.00/Y=5.00
	Legend/C/N=text1/J/B=(49151,53155,65535)
	
  endif

    
   //  TextBox/C/N=text0/F=0/A=LT/X=5/Y= 0/B=1 TBstr
//	SavePICT/P=datapath/E=-6/B=72 as "EDiffPlots"+StringFromList(n,CellListStr)+num2str(iframe)

  end
  
  
  
  
function MSD_scatter_plot(PlotNameStr,croot)
    string PlotNameStr   // Name of window for this plot 
    string  croot    // root of wave names
    
    
    make/o/n=(8,3) directcolorwave
	directcolorwave[1][1] = 65280
	directcolorwave[6][2] = 65280
	directcolorwave[7][0] = 57000
	directcolorwave[7][1] = 57000

	wave msdmag = $(croot+"msdc")
	wave msda = $(croot+"ac")
	wave maskw = $(croot+"maskc")
	wave forcew = $(croot+"forcec")
	


	dowindow/k $(PlotNameStr+"msdvsforce")
	Display/k=1 /N=$(PlotNameStr+"msdvsforce")/W=(209.25,292.25,630.75,514.25) msdmag vs forcew[*][0]
	ModifyGraph width={Aspect,2.5},height=108
	ModifyGraph mode=3
	ModifyGraph marker=19,msize=3
	ModifyGraph zColor($(croot+"msdc"))={$(croot+"maskc"),*,*,cindexRGB,0,directcolorwave}
	ModifyGraph log(left)=1
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph fSize=18
	ModifyGraph lblMargin(left)=13,lblMargin(bottom)=6,margin(top)=144,lblLatPos(left)=-90
	ModifyGraph standoff=1
	Label left "MSD(\\F'Symbol't\\F'Arial' = 10 s) (\\F'Symbol'm\\F'Arial'm\\S2\\M)"
	Label bottom "Traction force (nN)"
	Legend/C/N=text1/J/F=0/B=1/M/A=LT "\\K(0,65280,0)\\W1019 \\K(0,0,0)\\F'Arial'Cortical\r\\K(0,0,65280)\\W1019 \\K(0,0,0)\\F'Arial'Stress fiber";DelayUpdate
	AppendText "\\K(57000,57000,0)\\W1019 \\K(0,0,0)\\F'Arial'Not Assigned"
	Legend/C/N=text1/J/A=RB/X=5.00/Y=5.00
	Legend/C/N=text1/J/B=2
//	SetAxis left 1.1766452e-05,0.005
//	SetAxis bottom 0,40

	
	dowindow/k $(PlotNameStr+"msdvsslope")
	Display/k=1 /N=$(PlotNameStr+"msdvsslope")/W=(355.5,155,727.5,338.75) msdmag vs msda
	ModifyGraph width={Aspect,2.5},height=108
	ModifyGraph mode=3
	ModifyGraph marker=19,msize=3
	ModifyGraph zColor($(croot+"msdc"))={$(croot+"maskc"),*,*,cindexRGB,0,directcolorwave}
	ModifyGraph log(left)=1
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph fSize=18
	ModifyGraph lblMargin(bottom)=4
	ModifyGraph standoff=1	
	Label bottom "MSD exponent"
	
	dowindow/k $PlotNameStr
	Newlayout/k=1/N=$PlotNameStr
	ModifyLayout units=0
	AppendLayoutObject/F=0/R=(59,130,468,442) graph $(PlotNameStr+"msdvsforce")
    AppendLayoutObject/F=0/R=(92,79,467 ,266)  graph $(PlotNameStr+"msdvsslope")
	ModifyLayout frame=0,trans=1
	
end	

//--------------------QuickPlot() ------------------
// 6/9/17 makes a quick heatmap plot for diagnostic purposes


function QuickPlot(PlotNameStr, linmaskstackx, linmaskstacky, zwave,markersize)
    string PlotNameStr   // Name of window for this plot 
    wave  linmaskstackx  // one column wave of x  coordinates of points to be plotted
    wave  linmaskstacky  // one column wave of y  coordinates of points to be plotted
    wave zwave        //  wave containing scalar value to be plotted
    variable markersize // size of hexagons
    
    variable imageXsize = 1024 // default size of image  (1024x1024 for fastcamera, 1350x1000 for Roper)
    variable imageYsize = 1024
    
    variable zmin1 = 0.5     // these seem to be good values for Cells 1-3  
    variable zmax1 = 1.85
    variable i,j
    string TBstr,backgroundstr
    string nameCellw = NameOfWave(linmaskstacky)   // for use with modifygraph

 
//  make a dummy image
    sprintf backgroundstr, "%sbkd",PlotNameStr
        make/O/N=(imageXsize,ImageYsize)  $backgroundstr
    WAVE bkdwave = $backgroundstr

    dowindow/k $PlotNameStr

    newimage/k=1/N=$PlotNameStr bkdwave
         SetAxis top  -0.5, imageXsize + 0.5
        SetAxis/R left imageYsize+0.5, -0.5
     
// append the wave we care about
  appendtograph/t linmaskstacky vs linmaskstackx
      
        ModifyGraph zcolor($nameCellw) = {zwave,*,*,rainbow,0}
         
        

    ModifyGraph mode=3,marker=55,msize=markersize

     
  ModifyImage $backgroundstr explicit=1
//  ModifyImage $backgroundstr eval={0,56797,56797,56797} //gray background fill
  ModifyImage $backgroundstr eval={0,49151,53155,65535} // light blue background fill

    doupdate
    
   //  TextBox/C/N=text0/F=0/A=LT/X=5/Y= 0/B=1 TBstr
//	SavePICT/P=datapath/E=-6/B=72 as "EDiffPlots"+StringFromList(n,CellListStr)+num2str(iframe)

  end




//-- UTILITY ROUTINES  

// removes data of form (x,y,0) from x,y,z waves
 function removezeros(wx,wz)

wave wx,wz
variable i = 0
do
    if (wz[i] <=0)
        DeletePoints i,1,wx,wz
    else
        i +=1
    endif
while (i < dimsize(wx,0))
   
end


// removes data of form (x,y,#) from x,y,z waves where # != 0
 function keepzeros(wx,wy,wz)

wave wx,wy,wz
variable i = 0
do
    if (wz[i] !=0)
        DeletePoints i,1,wx,wy,wz
    else
        i +=1
    endif
while (i < dimsize(wx,0))
   
end

//--------- OLD VERSION -----------

//--------------------PlotOneFrameHeat02() ------------------
// 1/27/14
//  plots one frame of some scalar already calculated.  linmaskwaves must be peviously made
// Renamed and modified 5/30/17


function PlotOneFrameHeat02(PlotNameStr, linmaskstackx, linmaskstacky, zwave,scaleflag,ticflag,bkgndflag,bx,by,markersize,xmin,ymin,xmax,ymax)
    string PlotNameStr   // Name of window for this plot 
    wave  linmaskstackx  // one column wave of x  coordinates of points to be plotted
    wave  linmaskstacky  // one column wave of y  coordinates of points to be plotted
    wave zwave        //  wave containing scalar value to be plotted
    variable scaleflag  // = 0 to autoscale this frame; = 1 to scale from zero to max(zwave)
    variable ticflag  // = 0 for no tics, = 1  no tics
    variable bkgndflag // = 0 for no bkgnd posts, = 1 for bkgndposts
    wave bx    // x pos of all posts including bkgnd posts
    wave by    // y pos of all posts including bkgnd posts
    variable markersize // size of hexagons
    variable xmin //  these 4 set the plotting range (set xmin =-1 for default)
    variable ymin
    variable xmax
    variable ymax
    
    variable imageXsize = 1024 // default size of image  (1024x1024 for fastcamera, 1350x1000 for Roper)
    variable imageYsize = 1024
    
    variable zmin1 = 0.5     // these seem to be good values for Cells 1-3  
    variable zmax1 = 1.85
    variable i,j
    string TBstr,backgroundstr
    string nameCellw = NameOfWave(linmaskstacky)   // for use with modifygraph
    string nameBkdw = NameOfWave(by)

 
//  make a dummy image
    sprintf backgroundstr, "%sbkd",PlotNameStr
    if (xmin >=0) // use user-given image dimensions
        make/O/N=(xmax-xmin+1,ymax-ymin+1)  $backgroundstr
        SetScale/I x xmin,xmax,"", $backgroundstr
        SetScale/I y ymin,ymax,"", $backgroundstr
     else // default
        make/O/N=(imageXsize,ImageYsize)  $backgroundstr
    endif
    WAVE bkdwave = $backgroundstr

    dowindow/k $PlotNameStr

    newimage/k=1/N=$PlotNameStr bkdwave
    if (xmin >=0)
        SetAxis top xmin - 0.5, xmax + 0.5
        SetAxis/R left ymax+0.5,ymin -0.5
     else
         SetAxis top  -0.5, imageXsize + 0.5
        SetAxis/R left imageYsize+0.5, -0.5
     endif
     
// append the wave we care about
  appendtograph/t linmaskstacky vs linmaskstackx
      
    if (scaleflag == 0) // autoscale this frame
        ModifyGraph zcolor($nameCellw) = {zwave,*,*,rainbow,0}
         
        
    else //  scale from zero to maximum value of zwave
//         ModifyGraph zcolor = {colorwave,wavemin(zwave),wavemax(zwave),rainbow,0}
         ModifyGraph zcolor($nameCellw) = {zwave,zmin1,zmax1,rainbow,0}
   		
    endif   

    ModifyGraph mode=3,marker=55,msize=markersize

    if (ticflag == 0)
        ModifyGraph tick=3,noLabel=2
    endif
    
if (bkgndflag != 0)
     appendtograph/t by vs bx
        ModifyGraph mode=3 
//      ModifyGraph marker($nameBkdw)=55,msize($nameBkdw)=markersize,rgb($nameBkdw)=(56797,56797,56797)  // gray hexagons
//       ModifyGraph marker($nameBkdw)=55,msize($nameBkdw)=markersize,rgb($nameBkdw)=(49151,53155,65535) // blue hexagons
      ModifyGraph marker($nameBkdw)=19,msize($nameBkdw)=8,rgb($nameBkdw)=(56797,56797,56797)  // gray circles
//       ModifyGraph marker($nameBkdw)=19,msize($nameBkdw)=8,rgb($nameBkdw)=(49151,53155,65535)  // blue circles
  endif   
  ModifyImage $backgroundstr explicit=1
//  ModifyImage $backgroundstr eval={0,56797,56797,56797} //gray background fill
  ModifyImage $backgroundstr eval={0,49151,53155,65535} // light blue background fill

    doupdate
    
   //  TextBox/C/N=text0/F=0/A=LT/X=5/Y= 0/B=1 TBstr
//	SavePICT/P=datapath/E=-6/B=72 as "EDiffPlots"+StringFromList(n,CellListStr)+num2str(iframe)

  end




function seperatewave(dwave,dwave_index)
wave dwave , dwave_index
duplicate/o dwave dwave_sep

variable datalength = dimsize(dwave, 0)

variable i,j
j=0

for(i=0;i<datalength; i+=1)
    if(dwave_index[i] > 0 )
      dwave_sep[j][0] = dwave[i][0]
      dwave_sep[j][1] = dwave[i][1]

//      if(dwave[i] == 7)
//       dwave_sep[j] = -0.2
//      endif 
       j+=1
            
    endif
    
endfor

redimension/N=(j,2) dwave_sep

end      




Window MSDvsforce() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(209.25,292.25,630.75,514.25) CFc6msdc vs CFc6forcec[*][0]
	ModifyGraph width={Aspect,2},height=144
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph zColor(CFc6msdc)={CFc6maskc,*,*,cindexRGB,0,directcolorwave}
	ModifyGraph log(left)=1
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph fSize=18
	ModifyGraph lblMargin(left)=13,lblMargin(bottom)=6
	ModifyGraph standoff=0
	Label left "MSD(\\F'Symbol't\\F'Arial' = 10 s) (\\F'Symbol'm\\F'Arial'm\\S2\\M)"
	Label bottom "Traction force (nN)"
	SetAxis left 1.1766452e-05,0.005
	SetAxis bottom 0,40
EndMacro

Window MSDvsSlope() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(355.5,155,727.5,338.75) U2OSc3msdc vs U2OSc3ac
	ModifyGraph width={Aspect,2.5},height=108
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph zColor(U2OSc3msdc)={U2OSc3maskc,*,*,cindexRGB,0,directcolorwave}
	ModifyGraph log(left)=1
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph fSize=18
	ModifyGraph lblMargin(bottom)=4
	Label bottom "MSD exponent"
EndMacro
