#pragma rtGlobals=1		// Use modern global access method and strict wave access.

// 12/28/12 Yu changed DFcalcDLI function to move the heatmap plotting part into the main procedure file
//changed initial frame shift to be fitted by the number of data readed    9/9/2015

// maskwave constants

static constant kEmpty = 0
static constant kCell = 1
static constant kGuide = 2
static constant kIgnore = 3
static constant kWirepost = 4
static constant kFitAll = -1
static constant kInterp = 5     // 5/24/06  new constant interp base fit
 
  
 function DLI(defwave, freq, nof)
 
// variable s_mag, phase
 
 wave defwave
 
 
 variable  nof   //=1000
 variable  freq   //=1
// prompt freq  " enter frquency"
// prompt nof    " enter frame numbers"
// doprompt  " prompt name"  freq, nof
//string  S_prefix= "H_"

duplicate/O defwave Sig_X



 variable  nop
 Variable  SumRx=0, SumRy=0, SumSx=0, SumSy=0
 variable  i=0
 variable  S_mag, S_phase, R_mag, R_phase, phase
 variable fr = 100


make/O/N=2 DLIdata

// duplicate/O $Sigstr Sig_X 
// duplicate/O $Refstr  Ref_V

make/N=(nof)/O Ref_V
//Ref_V=sin(2*pi*freq*x)      //changded for testing
Ref_V=sin(2*pi*freq*x/fr)
 
nop=floor(nof*freq/fr)

 
 

 
 for( i=0 ; i<nof; i+=1)
     SumSx=SumSx+Sig_X[i]*cos(2*Pi*i*nop/nof)
     SumSy=SumSy+Sig_X[i]*Sin(2*Pi*i*nop/nof)
     SumRx=SumRx+Ref_V[i]*cos(2*Pi*i*nop/nof)
     SumRy=SumRy+Ref_V[i]*sin(2*Pi*i*nop/nof)
  endfor
  
  SumSx = SumSx/ nof
  SumSy = SumSy/ nof
  SumRx /= nof
  SumRy /= nof  
  
  S_mag=cabs(cmplx(SumSx,SumSy))*2
  R_mag=cabs(cmplx(SumRx,SumRy))*2
  S_phase=atan(SumSy/SumSx)-pi*0.5*(sign(sumsx)-1)
  R_phase=atan(SumRy/SumRx)-pi*0.5*(sign(sumrx)-1)
  phase=S_phase-R_phase
//  return S_mag
  
  
  DLIdata[0]= S_mag
  DLIdata[1]=phase
  

//  print   S_mag, phase,R_mag
  
  
   
 
end  


function calcallDLI(ASTRW,datawave, bgshiftwave, maskwave,framerate, freq, heatflag)
wave astrw
wave datawave
wave bgshiftwave
wave maskwave
variable framerate
variable freq
variable heatflag




SVAR FitRes_String

String heatmapS , heatmapSE

variable winsize

if (freq>=1) 
  sprintf heatmapS "DLIheatmap_%d", freq
  sprintf heatmapSE "DLIheatmaper_%d", freq
else 
  sprintf heatmapS "DLIheatmap0_%d",(freq*10)
  sprintf heatmapSE "DLIheatmaper0_%d", (freq*10)
endif

variable npostx, nposty, npost
variable i, j 
variable ipost
variable len, nof

NVAR MicronsPerPixel  //12/10/21 DHR
variable pixelratio=MicronsPerPixel

wave linmaskx, linmasky, DLIdata

//duplicate/O datawave fitresults



len= dimsize(datawave,2)



npostx=astrw[3]
nposty=astrw[4]
npost=npostx*nposty


make/N=(len)/O tempdwx, tempdwy
make/N=(npost,4)/O heatmap
make/N= (npost,2)/O heatmaperror

if (freq>=1) 
   winsize = 10
   else 
    winsize = floor(len*freq/(framerate*8))
endif


nof = floor(floor(len / (framerate/freq)) * framerate/freq)           //temporary change
//nof = 8968
for ( ipost =0; ipost<npost; ipost+=1)
      i = floor(ipost/nposty)
      j = mod(ipost, nposty)
      if(maskwave[i][j]!=Kignore)
      tempdwx= datawave[ipost][0][p] - datawave[ipost][0][0]
      tempdwy= datawave[ipost][1][p] - datawave[ipost][1][0]
      tempdwx = tempdwx - bgshiftwave[p][0]
      tempdwy = tempdwy - bgshiftwave[p][1]
      tempdwx = tempdwx * pixelratio
      tempdwy = tempdwy * pixelratio

      DLI(tempdwx, freq, nof)
      heatmap[ipost][0] = DLIdata[0]
      heatmap[ipost][1] = DLIdata[1]
      TimeVtest(tempdwx, freq, winsize)
      wavestats/Q DLItimeV
      heatmaperror[ipost][0] = V_sdev
      
      DLI(tempdwy, freq, nof)
      heatmap[ipost][2] = DLIdata[0]
      heatmap[ipost][3] = DLIdata[1]
      TimeVtest(tempdwy, freq, winsize)
      wavestats/Q DLItimeV
      heatmaperror[ipost][1] = V_sdev      
      else 
      
      heatmap[ipost][0]=0
      heatmap[ipost][1]=0
      heatmap[ipost][2]=0
      heatmap[ipost][3]=0
      
      endif
      
      
endfor
      
dowindow/K HeatMapdisp
newimage/N=HeatMapdisp M_rotatedimage
appendtograph/T/W=Heatmapdisp linmasky vs linmaskx
modifygraph mode=3, marker=0
modifygraph zcolor(linmasky)={heatmap[][heatflag],*,*, Rainbow, 0}


Save/O/P=ToPractice heatmap as heatmapS+".ibw"
Save/O/P=ToPractice heatmaperror as heatmapSE+".ibw"


end


function TimeVtest(dataw, freq, winsize)
wave dataw
variable freq   // frequency of wave
variable winsize   // window size for DLI ( number of periods)
variable i
variable wavelen = dimsize(dataw,0)
variable fr = 100
variable winsizef
variable DLIcounter = 0
//winsizef = winsize*fr/freq

if (freq<1)
winsizef = 3000
else
winsizef = 1000
endif

make/O/N=(winsizef)  DLItempw
make/O/N=(wavelen-winsizef) DLItimeV
wave DLIdata

for ( i = 0; i<(wavelen-winsizef); i+=winsizef)
      duplicate/O/R=((i), (i+winsizef)) dataw DLItempw
//      duplicate/O/R=(i,i+winsizef-1) dataw DLItempw
//      duplicate/O/R=((i/112),(i+winsizef-1)/112) dataw DLItempw
      DLI(DLItempw, freq, winsizef)
      DLItimeV[DLIcounter] = DLIdata[0] 
      DLIcounter+=1
 endfor

redimension/N=(DLIcounter) DLItimeV

end
     



function notchfilter(dataw, freq)
wave dataw
variable freq

wave filtertemp, filterfinal
//duplicate/O dataw filtertemp, filterfinal

variable fr = 112
variable freqr

freqr =freq/fr

filterFIR/NMF={freqr, 0.02} dataw
filterFIR/NMF={freqr*2, 0.02} dataw
filterFIR/NMF={freqr*3, 0.02} dataw
filterFIR/NMF={freqr*4, 0.02} dataw
filterFIR/NMF={freqr*5, 0.02} dataw
filterFIR/NMF={freqr*6, 0.02} dataw


//filterFIR/NMF={freqr, 0.02} filtertemp

//filterfinal = dataw - filtertemp

end



function calcfreqdepend()

SVAR Maskstring

loadwave/A/O/P=ToPractice "DLIheatmap0_1.ibw"
duplicate/O heatmap DLIheatmap0_1
loadwave/A/O/P=ToPractice "DLIheatmap0_2.ibw"
duplicate/O heatmap DLIheatmap0_2
loadwave/A/O/P=ToPractice "DLIheatmap0_5.ibw"
duplicate/O heatmap DLIheatmap0_5
loadwave/A/O/P=ToPractice "DLIheatmap_1.ibw"
duplicate/O heatmap DLIheatmap1
loadwave/A/O/P=ToPractice "DLIheatmap_2.ibw"
duplicate/O heatmap DLIheatmap2
loadwave/A/O/P=ToPractice "DLIheatmap_5.ibw"
duplicate/O heatmap DLIheatmap5
loadwave/A/O/P=ToPractice "DLIheatmaper0_1.ibw"
duplicate/O heatmaperror DLIheatmaper0_1
loadwave/A/O/P=ToPractice "DLIheatmaper0_2.ibw"
duplicate/O heatmaperror DLIheatmaper0_2
loadwave/A/O/P=ToPractice "DLIheatmaper0_5.ibw"
duplicate/O heatmaperror DLIheatmaper0_5
loadwave/A/O/P=ToPractice "DLIheatmaper_1.ibw"
duplicate/O heatmaperror DLIheatmaper1
loadwave/A/O/P=ToPractice "DLIheatmaper_2.ibw"
duplicate/O heatmaperror DLIheatmaper2
loadwave/A/O/P=ToPractice "DLIheatmaper_5.ibw"
duplicate/O heatmaperror DLIheatmaper5


runcalcfreqdepend(DLIheatmap0_1, DLIheatmap0_2, DLIheatmap0_5, DLIheatmap1, DLIheatmap2, DLIheatmap5, DLIheatmaper0_1, DLIheatmaper0_2, DLIheatmaper0_5, DLIheatmaper1, DLIheatmaper2, DLIheatmaper5,  astr, $maskstring)

end

function runcalcfreqdepend(dwave0_1, dwave0_2, dwave0_5, dwave1, dwave2, dwave5, ewave0_1, ewave0_2, ewave0_5, ewave1, ewave2, ewave5, astr, maskwave)
wave dwave0_1,dwave0_2, dwave0_5, dwave1, dwave2, dwave5
wave ewave0_1, ewave0_2, ewave0_5, ewave1, ewave2, ewave5
wave astr
wave maskwave

variable nx = astr[3]
variable ny = astr[4]
variable npost = nx*ny
variable ipost =0 
variable i = 0, j =0

make/N=6/O freq
freq[0] = 0.1
freq[1] = 0.2
freq[2] = 0.5
freq[3] = 1
freq[4] = 2
freq[5] = 5

make/N=(npost, 4, 6)/O heatmapall=0
make/N=(npost, 2,6)/O heatmapaller
for(ipost = 0 ; ipost <npost; ipost +=1)
       for (i = 0; i<4 ; i+=1)
          heatmapall[ipost][i][0] = dwave0_1[ipost][i]
          heatmapall[ipost][i][1] = dwave0_2[ipost][i]
          heatmapall[ipost][i][2] = dwave0_5[ipost][i]
          heatmapall[ipost][i][3] = dwave1[ipost][i]
          heatmapall[ipost][i][4] = dwave2[ipost][i]
          heatmapall[ipost][i][5] = dwave5[ipost][i]
      endfor
       for ( i = 0; i<2 ; i+=1)
          heatmapaller[ipost][i][0] = ewave0_1[ipost][i]
          heatmapaller[ipost][i][1] = ewave0_2[ipost][i]
          heatmapaller[ipost][i][2] = ewave0_5[ipost][i]
          heatmapaller[ipost][i][3] = ewave1[ipost][i]
          heatmapaller[ipost][i][4] = ewave2[ipost][i]
          heatmapaller[ipost][i][5] = ewave5[ipost][i]
       endfor
endfor


make/N=(6)/O heatmapcell 

display heatmapall[0][0][]/TN=hm_0 vs freq

string hmtrace
sprintf hmtrace "hm_%d" ipost

for (ipost = 1; ipost<npost; ipost+=1)
     i = floor(ipost/ny)
     j = mod(ipost, ny)
     sprintf hmtrace "hm_%d" ipost

     if (maskwave[i][j] != 3 )

    if (maskwave[i][j] == 0)
        appendtograph/C=(0,0, 62976) heatmapall[ipost][0][]/TN = $hmtrace  vs freq
        errorbars $hmtrace Y, wave=(heatmapaller[ipost][0][], heatmapaller[ipost][0][] )       
        modifygraph mode=4
//        modifygraph rgb(heatmapall) =(0,0, 62976)
     elseif (maskwave[i][j] == 1)
         appendtograph/C=(652800, 0 , 0 ) heatmapall[ipost][0][]/TN = $hmtrace vs freq
         errorbars $hmtrace Y, wave=(heatmapaller[ipost][0][], heatmapaller[ipost][0][] )             
//         appendtograph/C=(652800, 0 , 0 ) heatmapcell vs freq
         modifygraph mode=4         
//         modifygraph rgb(heatmapcell) =(652800, 0 , 0 )
    endif
    
    endif
endfor


modifygraph log(bottom) =1

end



//calculate DLI for double frequency measurement
function DFcalcDLI(ASTRW,datawave, bgshiftwave, maskwave,framerate, freq, reffreq, heatflag, initialshift)
wave astrw
wave datawave
wave bgshiftwave
wave maskwave
variable framerate
variable freq, reffreq
variable heatflag
variable initialshift 

SVAR FitRes_String

//String heatmapS , heatmapSE, heatmapREFSE, heatmapREFS

variable winsize, winsizeref

//if (freq>=1) 
//  sprintf heatmapS "DFDLIheatmap_%d", freq
//  sprintf heatmapSE "DFDLIheatmaper_%d", freq
//  sprintf heatmapREFSE "RefDFDLIheatmaper_%d", freq
//  sprintf heatmapREFS "RefDFDLIheatmap_%d", freq  
//else 
//  sprintf heatmapS "DFDLIheatmap0_%d",(freq*10)
//  sprintf heatmapSE "DFDLIheatmaper0_%d", (freq*10)
//  sprintf heatmapREFSE "RefDFDLIheatmaper0_%d", (freq*10)
//  sprintf heatmapREFS "RefDFDLIheatmap0_%d", (freq*10)  
//endif

variable npostx, nposty, npost
variable i, j 
variable ipost
variable len, nof, nofref

NVAR MicronsPerPixel  //12/10/21 DHR
variable pixelratio = MicronsPerPixel

wave linmaskx, linmasky, DLIdata

wave wave2

string HallS

//duplicate/O datawave fitresults


//temporary change   8/26/2015
len = dimsize(datawave,2) - 1000



//len= dimsize(datawave,2)





npostx=astrw[3]
nposty=astrw[4]
npost=npostx*nposty


make/N=(len)/O tempdwx, tempdwy
make/N=(npost,4)/O heatmap, heatmapref
make/N= (npost,2)/O heatmaperror, heatmaperrorref, heatmaprelative


winsizeref = 10

if (freq>=1) 
   winsize = 10

   else 
    winsize = floor(len*freq/(framerate*8))
endif

//change nof according to frequency
///nof = floor(floor(len / (framerate/freq)) * framerate/freq)           //temporary change
//nofref = floor(floor(len / (framerate/reffreq)) * framerate/reffreq)
//nof = 8968

nof = len - 1
nofref = len - 1

for ( ipost =0; ipost<npost; ipost+=1)
      i = floor(ipost/nposty)
      j = mod(ipost, nposty)
      if(maskwave[i][j]!=Kignore)
      tempdwx= datawave[ipost][0][p+initialshift] - datawave[ipost][0][initialshift]
      tempdwy= datawave[ipost][1][p+initialshift] - datawave[ipost][1][initialshift]
      tempdwx = tempdwx - bgshiftwave[p+initialshift][0]
      tempdwy = tempdwy - bgshiftwave[p+initialshift][1]
      tempdwx = tempdwx * pixelratio
      tempdwy = tempdwy * pixelratio

      DLI(tempdwx, freq, nof)
      heatmap[ipost][0] = DLIdata[0]
      heatmap[ipost][1] = DLIdata[1]
      DLI(tempdwx, reffreq, nofref)
      heatmapref[ipost][0] = DLIdata[0]
      heatmapref[ipost][1] =  DLIdata[1]
      TimeVtest(tempdwx, freq, winsize)
      wavestats/Q DLItimeV
      heatmaperror[ipost][0] = V_sdev
      TimeVtest(tempdwx, reffreq, winsizeref)
      wavestats/Q DLItimeV
      heatmaperrorref[ipost][0] = V_sdev
 
      
      DLI(tempdwy, freq, nof)
      heatmap[ipost][2] = DLIdata[0]
      heatmap[ipost][3] = DLIdata[1]
      DLI(tempdwy, reffreq, nofref)
      heatmapref[ipost][2] = DLIdata[0]
      heatmapref[ipost][3] = DLIdata[1]    
      TimeVtest(tempdwy, freq, winsize)
      wavestats/Q DLItimeV
      heatmaperror[ipost][1] = V_sdev        
      TimeVtest(tempdwy, reffreq, winsizeref)
      wavestats/Q DLItimeV
      heatmaperrorref[ipost][1] = V_sdev
      
      heatmaprelative[ipost][0] = heatmap[ipost][0]/heatmapref[ipost][0]
      heatmaprelative[ipost][1] = heatmap[ipost][1]/heatmapref[ipost][1]
            
      else 
      
      heatmap[ipost][0]=0
      heatmap[ipost][1]=0
      heatmap[ipost][2]=0
      heatmap[ipost][3]=0
      heatmapref[ipost][0]=0
      heatmapref[ipost][1]=0
      heatmapref[ipost][2]=0
      heatmapref[ipost][3]=0
      
      endif
      
      
endfor


//Save/O/P=ToPractice heatmap as heatmapS+".ibw"
//Save/O/P=ToPractice heatmaperror as heatmapSE+".ibw"
//Save/O/P=ToPractice heatmapref as heatmapREFS+".ibw"
//Save/O/P=ToPractice heatmaperrorref as heatmapREFSE+".ibw"

end



function DFcalcDLI_break(ASTRW,datawave, bgshiftwave, maskwave,framerate, freq, reffreq, heatflag,nfreq)
wave astrw
wave datawave
wave bgshiftwave
wave maskwave
variable framerate
variable freq, reffreq
variable heatflag
variable nfreq


SVAR FitRes_String

//String heatmapS , heatmapSE, heatmapREFSE, heatmapREFS

variable winsize, winsizeref

//if (freq>=1) 
//  sprintf heatmapS "DFDLIheatmap_%d", freq
//  sprintf heatmapSE "DFDLIheatmaper_%d", freq
//  sprintf heatmapREFSE "RefDFDLIheatmaper_%d", freq
//  sprintf heatmapREFS "RefDFDLIheatmap_%d", freq  
//else 
//  sprintf heatmapS "DFDLIheatmap0_%d",(freq*10)
//  sprintf heatmapSE "DFDLIheatmaper0_%d", (freq*10)
//  sprintf heatmapREFSE "RefDFDLIheatmaper0_%d", (freq*10)
//  sprintf heatmapREFS "RefDFDLIheatmap0_%d", (freq*10)  
//endif

variable npostx, nposty, npost
variable i, j, l
variable ipost
variable len, nof, nofref

NVAR MicronsPerPixel  //12/10/21 DHR
variable pixelratio = MicronsPerPixel

wave linmaskx, linmasky, DLIdata

variable bkpiece 
//duplicate/O datawave fitresults

wave DLItimeV


//make/o/n=17 noseg = {18,18,18,18,6,6,6,6,6,6,3,3,3,3,3,3,3}
make/o/n=17 noseg = {6,6,6,6,6,6,6,6,6,6,3,3,3,3,3,3,3}
bkpiece = noseg[nfreq]


len= dimsize(datawave,2)



npostx=astrw[3]
nposty=astrw[4]
npost=npostx*nposty


make/N=(len)/O tempdwx, tempdwy
make/N=(npost,2,bkpiece)/O heatmap_bk, heatmapref_bk
make/N= (npost,2)/O heatmaperror_bk, heatmaperrorref_bk
make/o/n=(npost, 2,17) relativeerrsingle
wave heatmaprelative1

    winsize = floor(len*freq/(framerate*bkpiece))
    winsizeref =  floor(len*reffreq/(framerate*bkpiece))



nof = floor(floor(len / (framerate/freq)) * framerate/freq)           //temporary change
nofref = floor(floor(len / (framerate/reffreq)) * framerate/reffreq)
//nof = 8968
for ( ipost =0; ipost<npost; ipost+=1)
      i = floor(ipost/nposty)
      j = mod(ipost, nposty)
      if(maskwave[i][j]!=Kignore)
      tempdwx= datawave[ipost][0][p] - datawave[ipost][0][0]
      tempdwy= datawave[ipost][1][p] - datawave[ipost][1][0]
      tempdwx = tempdwx - bgshiftwave[p][0]
      tempdwy = tempdwy - bgshiftwave[p][1]
      tempdwx = tempdwx * pixelratio
      tempdwy = tempdwy * pixelratio

//      DLI(tempdwx, freq, nof)
//      heatmap_bk[ipost][0] = DLIdata[0]
//      heatmap[ipost][1] = DLIdata[1]
//      DLI(tempdwx, reffreq, nofref)
//      heatmapref[ipost][0] = DLIdata[0]
//      heatmapref[ipost][1] =  DLIdata[1]
      TimeVtest(tempdwx, freq, winsize)
      for ( l =0 ; l<bkpiece; l+=1)
           heatmap_bk[ipost][0][l] = DLItimeV[l]
      endfor
      wavestats/Q DLItimeV
      duplicate/o DLItimeV DLItimeV_sigx
      heatmaperror_bk[ipost][0] = V_sdev
      TimeVtest(tempdwx, reffreq, winsizeref)
      for (l = 0; l<bkpiece; l+=1)
            heatmapref_bk[ipost][0][l] = DLItimeV[l]
      endfor
      wavestats/Q DLItimeV
      duplicate/o DLItimeV DLItimeV_refx
      DLItimeV_refx = 1/DLItimeV
      heatmaperrorref_bk[ipost][0] = V_sdev
 
      
//      DLI(tempdwy, freq, nof)
//      heatmap[ipost][2] = DLIdata[0]
//      heatmap[ipost][3] = DLIdata[1]
//      DLI(tempdwy, reffreq, nofref)
//      heatmapref[ipost][2] = DLIdata[0]
//      heatmapref[ipost][3] = DLIdata[1]    
      TimeVtest(tempdwy, freq, winsize)
      for ( l =0 ; l<bkpiece; l+=1)
           heatmap_bk[ipost][1][l] = DLItimeV[l]
      endfor      
      wavestats/Q DLItimeV
      duplicate/o DLItimeV DLItimeV_sigy
      heatmaperror_bk[ipost][1] = V_sdev        
      TimeVtest(tempdwy, reffreq, winsizeref)
      for ( l =0 ; l<bkpiece; l+=1)
           heatmapref_bk[ipost][1][l] = DLItimeV[l]
      endfor      
      wavestats/Q DLItimeV
      duplicate/o DLItimeV DLItimeV_refy
      DLItimeV_refy = 1/DLItimeV
      heatmaperrorref_bk[ipost][1] = V_sdev
      
//      heatmaprelative[ipost][0] = heatmap[ipost][0]/heatmapref[ipost][0]
//      heatmaprelative[ipost][1] = heatmap[ipost][1]/heatmapref[ipost][1]
      
      
      duplicate/o DLItimeV DLItimeV_relax
      duplicate/o DLItimeV DLItimeV_relay
      
      
      DLItimeV_relax = DLItimeV_sigx*DLItimeV_refx/mean(DLItimeV_refx)
      DLItimeV_relay = DLItimeV_sigy*DLItimeV_refy/mean(DLItimeV_refy)
      
      for (l=0;l<bkpiece;l+=1)
          heatmaprelative1[ipost][l][nfreq] = DLItimeV_relax[l]
      endfor    
      
      wavestats/q DLItimeV_relax
      
      relativeerrsingle[ipost][0][nfreq] = V_sdev
      
            
      else 
      
      heatmap_bk[ipost][0]=0
      heatmap_bk[ipost][1]=0
//      heatmap[ipost][2]=0
//      heatmap[ipost][3]=0
      heatmapref_bk[ipost][0]=0
      heatmapref_bk[ipost][1]=0
//      heatmapref[ipost][2]=0
//      heatmapref[ipost][3]=0

      
      endif
      
      
endfor

//dowindow/K HeatMapdisp
//newimage/N=HeatMapdisp M_rotatedimage
//appendtograph/T/W=Heatmapdisp linmasky vs linmaskx
//modifygraph mode=3, marker=0
//modifygraph zcolor(linmasky)={heatmap[][heatflag],*,*, Rainbow, 0}
//DoUpdate

//Save/O/P=ToPractice heatmap as heatmapS+".ibw"
//Save/O/P=ToPractice heatmaperror as heatmapSE+".ibw"
//Save/O/P=ToPractice heatmapref as heatmapREFS+".ibw"
//Save/O/P=ToPractice heatmaperrorref as heatmapREFSE+".ibw"

end




function loadscatterplot(postnum)

variable postnum

wave freq

wave heatmaprelative1

string heatmapS, heatmaprefS
variable i
variable nofreq = dimsize(freq,0)
variable j = 0, l=0

//make/o/n=129 freqscat, sigscat, refscat, relascat
make/o/n=81 freqscat, sigscat, refscat, relascat

wave heatmap_bk, heatmapref_bk

for (i = 0; i<nofreq; i+=1)



          if (freq[i] <1)
             sprintf heatmapS "DFDLIheatmap0_%d_bk",(freq[i]*10)
             sprintf heatmapREFS "RefDFDLIheatmap0_%d_bk", (freq[i]*10)
          else 
             sprintf heatmapS "DFDLIheatmap_%d_bk", freq[i]
             sprintf heatmapREFS "RefDFDLIheatmap_%d_bk", freq[i] 
          endif


loadwave/q/o/p=topractice heatmapS+".ibw"
loadwave/q/o/p=topractice heatmaprefS+".ibw"

variable scatterpiece = dimsize(heatmap_bk,2)

for(l=0;l<scatterpiece; l+=1)

    freqscat[j] = freq[i]
    sigscat[j] = heatmap_bk[postnum][0][l]
    refscat[j] = heatmapref_bk[postnum][0][l]
    j+=1
//    relascat[j] = heatmaprelative1[postnum][l][i]
endfor


endfor


relascat = sigscat/refscat*mean(refscat)


end    