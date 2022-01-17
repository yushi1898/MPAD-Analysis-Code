#pragma rtGlobals=1		// Use modern global access method and strict wave access.



//3/23/2016    changed frequency dependence procedure so that now it 

//magnetic field to force conversion: 0.009V ->0.3354nN


//double frequency measurment
function DLIAnalysis()
Dowindow/K W_DLI_Analysis
NewPanel/K=1/N=W_DLI_Analysis/W=(100,100,385,270)
SetDrawEnv linefgc= (0,0, 53520),fillpat= 0, linethick=2;	DrawRect 5,5,265,65
SetDrawEnv linefgc= (0,53520, 0),fillpat= 0, linethick=2;	DrawRect 5,75,265,160
SetDrawEnv linefgc= (64512,14848,14848),arrow= 1,linethick= 4.00,arrowfat= 1.00; DrawLine 275,10,275,150
SetDrawLayer UserBack

Button  READWAVE, pos={10,10},size={250,20}, proc = b_readwaves, title="Load posts DLI"
Button  RELATIVESIG, pos = {10,35}, size = {250,20}, proc = b_relativeplot, title = "Correct time variance"
Button  READHALLSEN pos={10,80},size={250,20}, proc = b_readhall, title="Load Hall sensor readout"
Button  CALCDLIMAG, pos = {10,105}, size = {250,20}, proc = b_calcmagDLI, title = "Calc Hall sensor DLI"
//Button  Vanhoffcurve, pos={10,80},size={250,20}, proc=genvanhoff, title = "Generate Vanhove"
Button  CALCMODU, pos={10,130},size={250,20}, proc=b_calcmodulus, title = "Calc Modulus"

end


function b_readwaves(ctrlName): ButtonControl
string ctrlName

readwaves()

end

function b_relativeplot(ctrlName): ButtonControl
string ctrlName

wave DLI_sum,DLIref_sum,freq

relativeplot(DLI_sum,DLIref_sum,freq)

end

function b_readhall(ctrlName): ButtonControl
string ctrlName

readhall()

end

function  b_calcmagDLI(ctrlName): ButtonControl
string ctrlName

wave hall_sum
wave freq
calcmagDLI(hall_sum,freq)

end


function  b_calcmodulus(ctrlName): ButtonControl
string ctrlName

wave DLI_sum
wave DLIref_sum
wave relative_sig
wave hallDLI_ed
wave freq

calcmodulus(DLI_sum,DLIref_sum,relative_sig,hallDLI_ed,freq)

end




function readwaves()
variable wl = 17
variable npost
variable i, ipost, j 
wave heatmap, heatmaperror, heatmapref, heatmaperrorref
String heatmapS , heatmapSE, heatmapREFSE, heatmapREFS



newpath/O DFpath 
//make/N=(wl)/o Freq= {0.1, 0.2, 0.5, 0.8, 1 ,2 ,4 ,5, 8, 10, 20 ,35,55,  80, 95, 115, 135}

wave freq
//4_6ms
make/N = 17/O corcoefficient = {1, 1.00265, 1.00695, 1.00839, 1.00952, 1.0106, 1.01071, 1.01021, 1.00921, 1.00834, 0.99706, 0.973722, 0.914896, 0.807937, 0.730442, 0.614863, 0.494301}
//3ms
//make/N = 17/O corcoefficient = {1, 1.00313, 1.00657, 1.00756, 1.00834, 1.00826, 1.009, 1.00974, 1.00969, 1.00957, 1.00351, 0.995937, 0.970181, 0.919107, 0.886404, 0.828752, 0.762414}

for (i = 0; i<wl; i+=1)
    if (freq[i]>=1) 
      sprintf heatmapS "DFDLIheatmap_%d", freq[i]
      sprintf heatmapSE "DFDLIheatmaper_%d", freq[i]
      sprintf heatmapREFSE "RefDFDLIheatmaper_%d", freq[i]
      sprintf heatmapREFS "RefDFDLIheatmap_%d", freq[i]  
    else 
      sprintf heatmapS "DFDLIheatmap0_%d",(freq[i]*10)
      sprintf heatmapSE "DFDLIheatmaper0_%d", (freq[i]*10)
      sprintf heatmapREFSE "RefDFDLIheatmaper0_%d", (freq[i]*10)
      sprintf heatmapREFS "RefDFDLIheatmap0_%d", (freq[i]*10)  
    endif
   
    loadwave/A/O/P=DFpath heatmapS+".ibw"
    if (i ==0)
    npost = dimsize(heatmap, 0)
    make/N=(npost, 4, wl)/O DLI_sum, DLIREF_sum
    make/N=(npost, 2, wl)/O DLI_err, DLIREF_err
    endif
    for(ipost = 0; ipost<npost; ipost+=1)
        for (j = 0; j<4; j+=1) 
              DLI_sum[ipost][j][i] = heatmap[ipost][j]
           if (j ==0 || j ==2)
              DLI_sum[ipost][j][i] = heatmap[ipost][j] / corcoefficient[i]
          endif
          
           if ( (j == 1 || j ==3) && (freq[i] >50 && freq[i] <100))
                 DLI_sum[ipost][j][i] = - heatmap[ipost][j]
              endif
              

       

        endfor
    endfor
    loadwave/A/O/P=DFpath heatmapSE+".ibw"
    for(ipost = 0; ipost<npost; ipost+=1)
        for (j = 0; j<2; j+=1) 
       DLI_err[ipost][j][i] = heatmaperror[ipost][j]
        endfor
    endfor
    loadwave/A/O/P=DFpath heatmapREFSE+".ibw"
    for(ipost = 0; ipost<npost; ipost+=1)
        for (j = 0; j<2; j+=1) 
       DLIREF_err[ipost][j][i] = heatmaperrorref[ipost][j]
        endfor
    endfor
    loadwave/A/O/P=DFpath heatmapREFS+".ibw"
    for(ipost = 0; ipost<npost; ipost+=1)
        for (j = 0; j<4; j+=1) 
              DLIREF_sum[ipost][j][i] = heatmapref[ipost][j]
                  if (j ==0 || j ==2)
              DLIREF_sum[ipost][j][i] = heatmapref[ipost][j] / corcoefficient[i]
          endif
        endfor
    endfor

endfor

end

function sigVSrefplot(sigwave, refwave, freq)
wave sigwave
wave refwave
wave freq

variable wl = dimsize(freq,0)
variable i
variable cor = 0 

for (i = 0; i<wl; i+=1)
     if (i==0)
        display sigwave[][cor][i] vs refwave[][cor][i]
     else 
         appendtograph sigwave[][cor][i] vs refwave[][cor][i]
     endif
     
endfor

end

function relativeplot(sigwave, refwave,freq)
wave sigwave
wave refwave
wave freq

variable wl =dimsize(freq,0)
variable i 
variable cor = 0
variable magthresh = 0.002
variable npost = dimsize(sigwave, 0)
variable ipost
variable mpcounter=0

wave DLI_err, DLIREF_err

string mptrace
make/N=(npost,2,wl)/O relativesig, relativeerr

make/N = (wl)/O refeditedx

make/N = 17/O corcoefficient = {1, 1.00265, 1.00695, 1.00839, 1.00952, 1.0106, 1.01071, 1.01021, 1.00921, 1.00834, 0.99706, 0.973722, 0.914896, 0.807937, 0.730442, 0.614863, 0.494301}

make/n=17/o noseg = {17,17,17,17,5,5,5,5,5,5,2,2,2,2,2,2,2}
dowindow/K relagraph
dowindow/K siggraph
dowindow/K refgraph






for (ipost = 0 ; ipost <npost; ipost+=1)
       if ((sigwave[ipost][0][0]>magthresh) && (refwave[ipost][0][0] > magthresh))
               mpcounter+=1
               refeditedx = refwave[ipost][0][p]

               wavestats/Q refeditedx

 //              refeditedx = refwave[ipost][0][p]/V_avg
                refeditedx = V_avg/refwave[ipost][0][p]
                wavestats/Q refeditedx

               
               
               for (i = 0; i<wl; i+=1)
                     relativesig[ipost][0][i] = sigwave[ipost][0][i]*refeditedx[i]/V_avg
                     
                     relativesig[ipost][1][i] = sigwave[ipost][1][i] - refwave[ipost][1][i]
                     
                     relativeerr[ipost][0][i] = relativesig[ipost][0][i]/V_avg*sqrt((DLI_err[ipost][0][i]/sigwave[ipost][0][i])^2+(DLIREF_err[ipost][0][i]/refwave[ipost][0][i])^2)/sqrt(noseg[i])
                     relativeerr[ipost][1][i] = relativesig[ipost][1][i]*sqrt((DLI_err[ipost][1][i]/sigwave[ipost][1][i])^2+(DLIREF_err[ipost][1][i]/refwave[ipost][1][i])^2)
                                          
               endfor
               
               sprintf mptrace "p_%d" ipost
               
             if(mpcounter==1)
                  display/N=relagraph relativesig[ipost][0][]/TN=$mptrace vs freq 
                  modifygraph/W=relagraph mode=4
                  errorbars/W=relagraph $mptrace Y, wave=(relativeerr[ipost][0][], relativeerr[ipost][0][] )   
                  display/N=siggraph sigwave[ipost][0][]/TN=$mptrace vs freq 
                  modifygraph/W=siggraph mode=4         
                  errorbars/W=siggraph $mptrace Y, wave=(DLI_err[ipost][0][], DLI_err[ipost][0][] )                            
                  display/N=refgraph refwave[ipost][0][]/TN=$mptrace vs freq      
                  modifygraph/W=refgraph mode=4          
                 errorbars/W=refgraph $mptrace Y, wave=(DLIREF_err[ipost][0][], DLIREF_err[ipost][0][] )                                                          
             else 
                  appendtograph/W=relagraph relativesig[ipost][0][]/TN=$mptrace vs freq
                  modifygraph/W=relagraph mode=4
                  errorbars/W=relagraph $mptrace Y, wave=(relativeerr[ipost][0][], relativeerr[ipost][0][] )             
                  appendtograph/W=siggraph sigwave[ipost][0][]/TN=$mptrace vs freq
                  modifygraph/W=siggraph mode=4                  
                  errorbars/W=siggraph $mptrace Y, wave=(DLI_err[ipost][0][], DLI_err[ipost][0][] )                        
                  appendtograph/W=refgraph refwave[ipost][0][]/TN=$mptrace vs freq
                  modifygraph/W=refgraph mode=4      
                  errorbars/W=refgraph $mptrace Y, wave=(DLIREF_err[ipost][0][], DLIREF_err[ipost][0][] )                                    
                                                      
             endif
        else
               for (i = 0; i<wl; i+=1)
                     relativesig[ipost][0][i] = 0
                     relativesig[ipost][1][i] = 0     
               endfor
        endif
        

        
endfor

modifygraph/W=relagraph log(bottom)=1
modifygraph/W=siggraph log(bottom)=1
modifygraph/W=refgraph log(bottom)=1

end





#pragma rtGlobals=3		// Use modern global access method and strict wave access.
function readwavesbk()
variable wl = 10
variable npost
variable i, ipost, j 
wave heatmap_bk, heatmaperror_bk, heatmapref_bk, heatmaperrorref_bk
String heatmapS , heatmapSE, heatmapREFSE, heatmapREFS, heatmap_R



newpath/O DFpath 
make/N=(wl)/o Freq= {0.1, 0.2, 0.5, 1 ,2 ,5, 10, 20 ,35, 80}

for (i = 0; i<wl; i+=1)
    if (freq[i]>=1) 
      sprintf heatmapS "DFDLIheatmap_%d", freq[i]
      sprintf heatmapSE "DFDLIheatmaper_%d", freq[i]
      sprintf heatmapREFSE "RefDFDLIheatmaper_%d", freq[i]
      sprintf heatmapREFS "RefDFDLIheatmap_%d", freq[i]  
      sprintf heatmap_R "Rela_heatmap_%d", freq[i]
    else 
      sprintf heatmapS "DFDLIheatmap0_%d",(freq[i]*10)
      sprintf heatmapSE "DFDLIheatmaper0_%d", (freq[i]*10)
      sprintf heatmapREFSE "RefDFDLIheatmaper0_%d", (freq[i]*10)
      sprintf heatmapREFS "RefDFDLIheatmap0_%d", (freq[i]*10)  
      sprintf heatmap_R "Rela_heatmap0_%d", (freq[i]*10)
    endif
   
    loadwave/A/O/P=DFpath heatmapS+"_bk.ibw"
    duplicate/O heatmap_bk $heatmapS
    loadwave/A/O/P=DFpath heatmapSE+"_bk.ibw"
    duplicate/O heatmaperror_bk $heatmapSE
    loadwave/A/O/P=DFpath heatmapREFSE+"_bk.ibw"  
    duplicate/O heatmaperrorref_bk $heatmapREFSE
    loadwave/A/O/P=DFpath heatmapREFS+"_bk.ibw"     
    duplicate/O heatmapref_bk $heatmapREFS 
    duplicate/O heatmap_bk heatmap_relative
    heatmap_relative = heatmap_bk/heatmapref_bk
    duplicate/O heatmap_relative $heatmap_R
endfor

end

function showpost(ipost)
variable ipost

variable i
variable wl = 17
make/N=17/O freq={0.1,0.2,0.5,0.8,1.0,2.0,4.0,5.0,8.0,10.0,20.0,35.0,55.0,80.0,155.0,255.0,455.0}
String heatmap_R

for (i = 0; i<wl; i+=1)
    if (freq[i]>=1) 
      sprintf heatmap_R "Rela_heatmap_%d", freq[i]
    else 
      sprintf heatmap_R "Rela_heatmap0_%d", (freq[i]*10)
    endif
    
    if (i==0)
       display $heatmap_R[ipost][0][]
       modifygraph mode = 4 
    else 
       appendtograph $heatmap_R[ipost][0][]
       modifygraph mode = 4
    endif
endfor


end
                

function calcmagDLI(Hall_sum, freq)
wave Hall_sum
wave freq

variable fl = dimsize(freq, 0)
variable i 
variable fr = 100
variable freqthresh = fr/2
variable nof
make/N = 18000/O halltempL, halltempR

make/N=(fl, 4)/O  hallDLI, hallDLIref
make/N=(fl, 2)/O hallDLI_ed, hallDLIrelative
make/N=(fl)/O  hallrefmag, hallrefphase


wave DLIdata

variable reffreq = 7

variable exptime = 0.002

variable ampL, ampR, phaseL, phaseR
variable refampL, refampR, refphaseL, refphaseR


//4_6ms
//make/N = 17/O phasecorrection = {-3.14001, -3.12622, -3.08992, -3.05327, - 3.0362, -2.92371, -2.69938, -2.58692, -2.25123, -2.03717, -0.962077, 0.681201,-3.42915,-0.72637, 0.879343,3.08348, -0.971743}
//1ms
//make/N = 17/O phasecorrection = {-3.14057, -3.12996, -3.009603, -3.06718, -3.0443, -2.93862, -2.74025, -2.64367, -2.36274, -2.18631, -1.19882, 0.281164, 2.24278, 4.65587,6.16947, 8.11409, 10.0749}
//4_5ms
//make/N = 17/O phasecorrection = {-3.14074, -3.12607, -3.08916, -3.05674, -3.03304, -2.92271, -2.70025, -2.58973, -2.27524, -1.98821, -0.976614, 0.61697, 1.29518, -2.95496, -3.77504, -3.16129, -1.08963}

make/N = 17/O   phasecorrection //= { 0.00141372, 0.00921646, 0.0109753, 0.0142714, 0.0166573, 0.029963, 0.057912, 0.0720243, 0.114687, 0.142973, 0.283737, 0.496986, 0.77686, 1.12883, 1.33932, 1.6164, 1.89704}

phasecorrection = pi*exptime*freq[p]



for (i = 0; i<fl; i+=1)
        halltempL = Hall_sum[p][0][i]
        halltempR = Hall_sum[p][1][i]
        if (freq[i]< 1)  
           nof = 17000
         elseif (freq[i] >=1 &&freq[i] < 20)
           nof = 5000
          else
           nof = 2000
        endif
  
  
   
            DLI( halltempL, reffreq , nof)
            hallDLIref[i][0] = DLIdata[0]
            hallDLIref[i][1] = DLIdata[1]
            DLI( halltempR, reffreq, nof)
            hallDLIref[i][2] = DLIdata[0]
            hallDLIref[i][3] = DLIdata[1]        
        
        
        refampL = hallDLIref[i][0]
        refampR = hallDLIref[i][2]
        refphaseL = hallDLIref[i][1]
        refphaseR = hallDLIref[i][3]
        
        hallrefmag = sqrt(refampL^2 + refampR^2+2* refampL*refampR*cos(refphaseL - refphaseR))
        hallrefphase = atan((refampL*sin(refphaseL)+refampR*sin(refphaseR))/(refampL*cos(refphaseL)+refampR*cos(refphaseR)))-pi*0.5*(sign(refampL*cos(refphaseL)+refampR*cos(refphaseR))-1)       
       
        
        
        if(freq[i] < freqthresh)
            DLI( halltempL, freq[i] , nof)
            hallDLI[i][0] = DLIdata[0]
            hallDLI[i][1] = DLIdata[1]
            DLI( halltempR, freq[i], nof)
            hallDLI[i][2] = DLIdata[0]
            hallDLI[i][3] = DLIdata[1]
         elseif(freq[i] >freqthresh && freq[i] < fr)
            DLI(halltempL, abs(fr - freq[i]), nof)
            hallDLI[i][0] = DLIdata[0]
            hallDLI[i][1] = -DLIdata[1]
            DLI( halltempR,abs(fr - freq[i]), nof)
            hallDLI[i][2] = DLIdata[0]
            hallDLI[i][3] = -DLIdata[1]
          elseif(freq[i] >fr)
            DLI(halltempL, abs(fr - freq[i]), nof)
            hallDLI[i][0] = DLIdata[0]
            hallDLI[i][1] = DLIdata[1]
            DLI( halltempR,abs(fr - freq[i]), nof)
            hallDLI[i][2] = DLIdata[0]
            hallDLI[i][3] = DLIdata[1]          
        endif
        ampL = hallDLI[i][0]
        ampR = hallDLI[i][2]
        phaseL = hallDLI[i][1]
        phaseR = hallDLI[i][3]
 

 
        
        
        hallDLI_ed[i][0] = sqrt(ampL^2 + ampR^2+2* ampL*ampR*cos(phaseL - phaseR))
        hallDLI_ed[i][1] = atan((ampL*sin(phaseL)+ampR*sin(phaseR))/(ampL*cos(phaseL)+ampR*cos(phaseR)))-pi*0.5*(sign(ampL*cos(phaseL)+ampR*cos(phaseR))-1)       
        
        
        hallDLI_ed[i][1] = hallDLI_ed[i][1] - phasecorrection[i]// + pi

                
endfor

wavestats/Q hallrefmag

    for( i = 0; i<fl; i+=1)

       hallDLIrelative[i][0] = hallDLI_ed[i][0]/V_avg
       hallDLIrelative[i][1] = hallDLI_ed[i][1]
       
    endfor
       
       



end




function readhall()
variable wl =17
variable i, j  
string Halls

wave wave0, wave1

newpath/O DFpath 
make/N=(wl)/o Freq= {0.1, 0.2, 0.5, 0.8, 1 ,2 ,4 ,5, 8, 10, 20 ,35,55,  80, 95, 115, 135}

make/N=(25000, 2, wl)/O Hall_sum

for (i = 0 ; i<wl; i+=1)
//    if (freq[i]>=1) 
//      sprintf HallS "magread_%d.txt", freq[i]
//    else 
//      sprintf HallS "magread_0_%d.txt",(freq[i]*10) 
//    endif
    sprintf HallS "magread_%1f.txt", freq[i]
    
    

    loadwave/N/O/P=DFpath/G HallS
    
    
    
    for ( j=0; j < dimsize(wave0,0); j+=1)
    Hall_sum[j][0][i] = wave0[j]
    Hall_sum[j][1][i] = wave1[j]
    endfor
    
    
    endfor
    
end





function calcmodulus(sigwave, refwave, relawave, HallDLI, freq)
wave sigwave, refwave, relawave
wave HallDLI
wave freq


variable wl =dimsize(freq,0)
variable i 
variable cor = 0
variable magthresh = 0.001
variable npost = dimsize(sigwave, 0)
variable ipost
variable mpcounter=0
variable modumag
variable moduphase

//wave DLI_err, DLIREF_err

string mptrace

wave relativeerr, relativeerrsingle
make/N=(npost,4,wl)/O modulusw, modulusw_ed
make/o/n=(npost,wl) modmag_err


dowindow/K modustoregraph
dowindow/K modulossgraph
dowindow/K modumaggraph
dowindow/K moduphasegraph

for (ipost = 0 ; ipost <npost; ipost+=1)
       if ((sigwave[ipost][0][0]>magthresh)  && (refwave[ipost][0][0] > magthresh))
               mpcounter+=1
               for (i = 0; i<wl; i+=1)
                     modumag = HallDLI[i][0]/relawave[ipost][0][i]
                     moduphase = mod((HallDLI[i][1] - sigwave[ipost][1][i]),(2*pi))
                     if (moduphase>5)
                        moduphase = moduphase- 2*pi
                     endif
               
                     modulusw[ipost][0][i] = modumag
                     modulusw[ipost][1][i] = moduphase
                     modulusw[ipost][2][i] = abs(modumag*cos(moduphase))                  
                     modulusw[ipost][3][i] = abs(modumag*sin(moduphase))
                     
                     modulusw_ed[ipost][0][i] = modulusw[ipost][0][i]/modulusw[ipost][0][0]
                     modulusw_ed[ipost][1][i] = modulusw[ipost][1][i]
                     modulusw_ed[ipost][2][i] = modulusw[ipost][2][i]/modulusw[ipost][2][0]
                     modulusw_ed[ipost][3][i] = modulusw[ipost][3][i]/modulusw[ipost][3][0]   
                     
                     modmag_err[ipost][i] = modumag*relativeerrsingle[ipost][0][i]/relawave[ipost][0][i]                  
                     
//                     relativeerr[ipost][0][i] = relativesig[ipost][0][i]*sqrt((DLI_err[ipost][0][i]/sigwave[ipost][0][i])^2+(DLIREF_err[ipost][0][i]/refwave[ipost][0][i])^2)
//                     relativeerr[ipost][1][i] = relativesig[ipost][1][i]*sqrt((DLI_err[ipost][1][i]/sigwave[ipost][1][i])^2+(DLIREF_err[ipost][1][i]/refwave[ipost][1][i])^2)
                                          
               endfor
               
               sprintf mptrace "p_%d" ipost

               
             if(mpcounter==1)
                  display/N=modustoregraph modulusw[ipost][2][]/TN=$mptrace vs freq 
                  modifygraph/W=modustoregraph mode=4
                  display/N = modulossgraph modulusw[ipost][3][]/TN = $mptrace vs freq
                  modifygraph/W=modulossgraph mode=4           
                  display/N = modumaggraph modulusw[ipost][0][]/TN = $mptrace vs freq
                  modifygraph/W=modumaggraph mode=4     
                  display/N = moduphasegraph modulusw[ipost][1][]/TN = $mptrace vs freq
                  modifygraph/W=moduphasegraph mode=4                                                
//                  errorbars/W=relagraph $mptrace Y, wave=(relativeerr[ipost][0][], relativeerr[ipost][0][] )   
//                  display/N=siggraph sigwave[ipost][0][]/TN=$mptrace vs freq 
 //                 modifygraph/W=siggraph mode=4         
 //                 errorbars/W=siggraph $mptrace Y, wave=(DLI_err[ipost][0][], DLI_err[ipost][0][] )                            
//                  display/N=refgraph refwave[ipost][0][]/TN=$mptrace vs freq      
//                  modifygraph/W=refgraph mode=4          
//                 errorbars/W=refgraph $mptrace Y, wave=(DLIREF_err[ipost][0][], DLIREF_err[ipost][0][] )                                                          
             else 
                  appendtograph/W=modustoregraph modulusw[ipost][2][]/TN=$mptrace vs freq
                  modifygraph/W=modustoregraph mode=4
                  appendtograph/W=modulossgraph modulusw[ipost][3][]/TN=$mptrace vs freq
                  modifygraph/W=modulossgraph mode=4  
                  appendtograph/W=modumaggraph modulusw[ipost][0][]/TN=$mptrace vs freq
                  modifygraph/W=modumaggraph mode=4 
                  appendtograph/W=moduphasegraph modulusw[ipost][1][]/TN=$mptrace vs freq
                  modifygraph/W=moduphasegraph mode=4 
                                                                      
//                  errorbars/W=relagraph $mptrace Y, wave=(relativeerr[ipost][0][], relativeerr[ipost][0][] )             
//                  appendtograph/W=siggraph sigwave[ipost][0][]/TN=$mptrace vs freq
//                  modifygraph/W=siggraph mode=4                  
//                  errorbars/W=siggraph $mptrace Y, wave=(DLI_err[ipost][0][], DLI_err[ipost][0][] )                        
//                  appendtograph/W=refgraph refwave[ipost][0][]/TN=$mptrace vs freq
//                  modifygraph/W=refgraph mode=4      
//                  errorbars/W=refgraph $mptrace Y, wave=(DLIREF_err[ipost][0][], DLIREF_err[ipost][0][] )                                    
                                                      
             endif
        else
               for (i = 0; i<wl; i+=1)
                     modulusw[ipost][0][i] = 0
                     modulusw[ipost][1][i] = 0     
                     modulusw[ipost][2][i] = 0
                     modulusw[ipost][3][i] = 0   
                                            
               endfor
        endif
        

        
endfor

modifygraph/W=modustoregraph log(bottom)=1
modifygraph/W=modulossgraph log(bottom)=1
modifygraph/W=modumaggraph log(bottom)=1
modifygraph/W=moduphasegraph log(bottom)=1

end


function calcfitpara(cwave,pwave,ratio)
wave cwave, pwave
variable ratio

wave freqlog, noseg


duplicate/o cwave difwave
difwave[][0] = cwave[p][0]-pwave[p][0]*ratio
difwave[][1] = sqrt(cwave[p][1]^2+pwave[p][1]^2*ratio^2)

duplicate/o difwave difwave_log
difwave_log[][0] = log(difwave[p][0])
difwave_log[][1] = difwave[p][1]/difwave[p][0]//*sqrt(noseg[p])

variable delta, alpha, bbeta, alpha_err, bbeta_err
make/o/n=11 invsigmasq, xwsq,yw, xwyw, xw
invsigmasq = (1/difwave_log[p][1])^2
xwsq =  freqlog^2*invsigmasq
yw = difwave_log[p][0]*invsigmasq
xw = freqlog*invsigmasq
xwyw = difwave_log[p][0]*freqlog*invsigmasq

delta = sum(invsigmasq)*sum(xwsq) - (sum(xw))^2
alpha = 1/delta * (sum(xwsq)*sum(yw) - sum(xw)*sum(xwyw))
bbeta = 1/delta* (sum(invsigmasq)*sum(xwyw) - sum(xw)*sum(yw))

alpha_err = 1/delta*(sum(xwsq))
bbeta_err = 1/delta*(sum(invsigmasq))

print(bbeta)
print(bbeta_err)

end



function dispy2xratio(datawave, astr, postnum)
wave datawave, astr
variable postnum
variable i
variable radtheta = abs(astr[2])/180*pi
variable ymag, xmag


make/O/N=17  y2xratio

for (i = 0; i<17;i+=1)
     ymag = datawave[postnum][2][i]
     xmag = datawave[postnum][0][i]
//     y2xratio[i] = datawave[postnum][2][i]*/datawave[postnum][0][i]
       y2xratio[i] = (ymag*cos(radtheta)-xmag*sin(radtheta))/(xmag*cos(radtheta)+ymag*sin(radtheta))
endfor

wavestats/q y2xratio

print(V_avg)
//print(atan(V_avg)/pi*180-astr[2])



end


function calcfactor(moduwave, changedwave,ipost)
wave moduwave
wave changedwave
variable ipost

variable i
variable nof = dimsize(changedwave, 0)

make/o/n=(nof)  tempmodu
tempmodu = moduwave[ipost][0][p]

variable cfactor=0

for( i = 0; i<nof; i+=1)
     cfactor = cfactor + tempmodu[i]/changedwave[i]
//     cfactor = cfactor + 1
endfor


cfactor = cfactor/nof


print(cfactor)


end
     
