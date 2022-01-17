pl#pragma rtGlobals=1		// Use modern global access method and strict wave access.

//version 1_31 : fit strain fluctuation curve from minum to 1, added substracted posts slope
//version 1_2  : enable plotallfft to calculate SD and error

static constant kEmpty = 0
static constant kCell = 1
static constant kGuide = 2
static constant kIgnore = 3
static constant kWirepost = 4
static constant kFitAll = -1
static constant kInterp = 5     // 5/24/06  new constant interp base fit

function indpassivenoise(datawave, bintime, postid, freq)
wave datawave
variable bintime, postid, freq

variable framerate = 100
variable movielength = dimsize(datawave, 2)
variable iframe

make/N=(movielength)/O tempdwave
tempdwave = datawave[postid][0][p]
if (freq < 2)
notchfilter2(tempdwave, freq)
else
notchfilter4(tempdwave,freq)
endif
notchfilter3(tempdwave, 7)


variable binsize = bintime*framerate

variable i , j 
j=0
make/O/n=(binsize) stemp
variable tempwsize = floor((movielength-1)/binsize)

make/O/N= (tempwsize) noisetemp

for (iframe = 0 ; iframe < (movielength -1) ; iframe+=1)
     i = mod(iframe, binsize)
     j = floor(iframe/binsize)
     stemp[i] = tempdwave[iframe]         
     if ( i ==(binsize-1) )
        wavestats/Q stemp
        noisetemp[j] = V_sdev
     endif


endfor

for( i = 0; i<5; i+=1) 
   noisetemp[i] = noisetemp[5]
   noisetemp[j-i] = noisetemp[j-5]
endfor


     
      
end         
          


function calcpassivenoise(freq, postid)
wave freq
variable postid

variable bintime = 1
variable framerate = 100
variable avgtime = 5

SVAR FitRes_string

variable nofreq = dimsize(freq, 0)
variable i 
variable iframestart , iframe, jframe




string fitaffix, fitname

wave  noisetemp

make/O/N=(nofreq)  videolength
videolength = {180,180,180,180,60,60,60,60,60,60,30,30,30,30,30,30,30}

//make/O/N=(nofreq)  magnoise   first test , one video one data point
make/O/N=(sum(videolength)/avgtime) magnoise


make/O/N=(avgtime)  noisetemp_s

for (i=0; i<nofreq; i+=1)
  if (freq[i] <1)
         sprintf fitaffix "_0_%dhz"  (freq[i]*10)
  else
         sprintf fitaffix "_%dhz"  freq[i]
  endif
  
  fitname = FitRes_String+"_ed"+fitaffix
  
  if (freq[i] < framerate/2)
     indpassivenoise($fitname, bintime, postid, freq[i])
  else
      indpassivenoise($fitname, bintime, postid, abs(freq[i]-100) )
  endif
  

  
  if (i == 0)
      iframestart = 0
  else 
      iframestart = sum(videolength, 0, (i-1) )/avgtime
  endif
  
  
  for ( iframe = 0; iframe< videolength[i] ; iframe+=1)
       jframe = mod(iframe, avgtime)
       noisetemp_s[jframe] = noisetemp[iframe]
    if ( jframe == (avgtime-1))
          wavestats/Q noisetemp_s
        magnoise[floor(iframe/avgtime)+iframestart] = V_avg
    endif
  endfor

endfor


end



function notchfilter2(dataw, freq)
wave dataw
variable freq

wave filtertemp, filterfinal
//duplicate/O dataw filtertemp, filterfinal

variable fr = 100
variable freqr

freqr =freq/fr

filterFIR/NMF={freqr,0.001} dataw
filterFIR/NMF={freqr*2, 0.001} dataw
//filterFIR/NMF={freqr*3, 0.02} dataw
//filterFIR/NMF={freqr*4, 0.02} dataw
//filterFIR/NMF={freqr*5, 0.02} dataw
///filterFIR/NMF={freqr*6, 0.02} dataw


//filterFIR/NMF={freqr, 0.02} filtertemp

//filterfinal = dataw - filtertemp

end

function notchfilter4(dataw, freq)
wave dataw
variable freq

wave filtertemp, filterfinal
//duplicate/O dataw filtertemp, filterfinal

variable fr = 100
variable freqr

freqr =freq/fr

filterFIR/NMF={freqr,0.0001} dataw
///filterFIR/NMF={freqr*2, 0.02} dataw
//filterFIR/NMF={freqr*3, 0.02} dataw
//filterFIR/NMF={freqr*4, 0.02} dataw
//filterFIR/NMF={freqr*5, 0.02} dataw
///filterFIR/NMF={freqr*6, 0.02} dataw


//filterFIR/NMF={freqr, 0.02} filtertemp

//filterfinal = dataw - filtertemp

end



function notchfilter3(dataw, freq)
wave dataw
variable freq

wave filtertemp, filterfinal
//duplicate/O dataw filtertemp, filterfinal

variable fr = 100
variable freqr

freqr =freq/fr

filterFIR/NMF={freqr,0.005} dataw
filterFIR/NMF={freqr*2, 0.005} dataw
filterFIR/NMF={freqr*3, 0.005} dataw
//filterFIR/NMF={freqr*4, 0.02} dataw
//filterFIR/NMF={freqr*5, 0.02} dataw
///filterFIR/NMF={freqr*6, 0.02} dataw


//filterFIR/NMF={freqr, 0.02} filtertemp

//filterfinal = dataw - filtertemp

end

function TimeVtest2(dataw, freq, winsize)
wave dataw
variable freq   // frequency of wave
variable winsize   // window size for DLI ( number of periods)
variable i
variable wavelen = dimsize(dataw,0)
variable fr = 100
variable winsizef
variable DLIcounter = 0
//winsizef = winsize*fr/freq

winsizef = fr*winsize
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



function calcrefDLI(freq, postid)
wave freq
variable postid

variable reffreq = 7
variable framerate = 100
variable avgtime = 5

SVAR FitRes_string

variable nofreq = dimsize(freq, 0)
variable i 
variable iframestart , iframe, jframe, gapstart




string fitaffix, fitname
wave DLItimeV
wave  noisetemp

make/O/N=(nofreq)  videolength
videolength = {90,60,60,60,60,30,30,30,30,30}
//make/O/N=(nofreq)  magnoise   first test , one video one data point
make/O/N=(sum(videolength)/avgtime) refDLItimedep , timewave

make/O/N=(nofreq-1) gaplength
gaplength = {35,10,10,10,10,10,10,10,10}


for (i=0; i<nofreq; i+=1)
  if (freq[i] <1)
         sprintf fitaffix "_0_%dhz"  (freq[i]*10)
  else
         sprintf fitaffix "_%dhz"  freq[i]
  endif
  
  fitname = FitRes_String+"_ed"+fitaffix
  
  make/O/N=(videolength[i]*framerate)   DLItempwave
  
  duplicatesinglepost($fitname, DLItempwave, postid)
  
  TimeVtest2(DLItempwave, reffreq, avgtime)
  
  if (i == 0)
      iframestart = 0
      gapstart = avgtime/2
  else 
      iframestart = sum(videolength, 0, (i-1) )/avgtime
      gapstart = avgtime/2 + iframestart*avgtime + sum(gaplength, 0, (i-1))
  endif
  
  
  for ( iframe = 0; iframe< (videolength[i]/avgtime) ; iframe+=1)
        refDLItimedep[iframe+iframestart] = DLItimeV[iframe]
        timewave[iframe+iframestart] = gapstart+avgtime*iframe
  endfor
  
endfor  
    
  
  
  end
  
  
  
  
  
  
  
  
  
  
  
  
function duplicatesinglepost(datawave,destwave, postid)
wave datawave, destwave
variable postid

destwave = datawave[postid][0][p]

end

  


function calcpowerspectrum(freq, postid, magflag)
variable freq
variable postid
variable magflag



SVAR FitRes_string

variable framerate=100
variable reffreq = 7
string fitaffix, fitname
variable movielength = 18000
variable timelength = 180
variable cf = 0.1


if(freq ==0)
  movielength = dimsize(fitresults,2)-1
  timelength = movielength/framerate
endif
make/O/N=(movielength) pstemp

variable nlowf = (movielength*cf/framerate)-1

wave fitres_ed, fitresults

  if (freq <1)
         sprintf fitaffix "_0_%dhz"  (freq*10)
  else
         sprintf fitaffix "_%dhz"  freq
  endif

if (freq!=0)
fitname = FitRes_String+"_ed"+fitaffix
else 
fitname = "fitres_ed"

endif

duplicatesinglepost($fitname, pstemp, postid)

if(magflag == 1)
     notchfilter2(pstemp, freq)
     notchfilter3(pstemp, reffreq)
endif



setscale/P x  0 ,(1/framerate), "", pstemp

FFT/MAGS/DEST=PS_FFT   pstemp


variable nopfreq = 30
variable fftcounter=0
variable i=0
variable j=0


make/O/N= (nlowf+nopfreq+1) freqlogx, fft_ed
fft_ed = 0
make/O/N=(floor(movielength/2)+1) fft_x


for ( i = 0; i< nlowf; i +=1)
     freqlogx[i] = (i+1)/timelength
endfor

for (i = nlowf; i < (nlowf+nopfreq + 1); i+=1)
     freqlogx[i] = 10^((i-nlowf)*(log(50)-log(cf))/nopfreq+log(cf))
endfor

//freqlogx = 10^(x*(2+log(50))/nopfreq-2)
fft_x = (x)/timelength

for(i = 1; i<(floor(movielength/2)+1); i+=1)

    if ( (fft_x[i] >freqlogx[j]) || ( (i+1)==floor(movielength/2)) )
       fft_ed[j] = fft_ed[j]/fftcounter
       j=j+1
   //    print(fftcounter)
       fftcounter = 0
    endif            
       fft_ed[j] = fft_ed[j]+PS_FFT[i]
       fftcounter = fftcounter+1
    

endfor
    
    

    



fft_ed = fft_ed/((movielength/2)^2)/8^2     





end 

 
  

function plotallfft(maskwave, astr, freq)
wave maskwave
wave astr
variable freq

variable npostx = astr[3]
variable nposty = astr[4] 
variable npost = npostx*nposty

variable ipost, xpost, ypost

variable Vred, Vgreen, Vblue
variable magflag



variable nfreqlog = dimsize(freqlogx, 0)

variable i
string graphname , ffttracename
variable ncell = 0, nbg=0

wave fft_ed

sprintf graphname, "fft_%d", (freq*10)

dowindow/k $graphname

make/O/N=(npost, nfreqlog) fft_ed_all

make/O/N=(nfreqlog) fft_ed_bg, fft_ed_cell, fft_ed_bg_sd, fft_ed_cell_sd
fft_ed_bg = 0
fft_ed_cell = 0
fft_ed_bg_sd = 0
fft_ed_cell_sd = 0

for (ipost = 0 ; ipost<npost; ipost+=1)
     xpost = floor(ipost/nposty)
     ypost = mod(ipost, nposty)
     
     sprintf ffttracename, "fft_ed_%d", ipost
     if(maskwave[xpost][ypost] ==4)
        magflag=1
        else
        magflag = 0        
        endif
     
     if (maskwave[xpost][ypost] != 3)
        calcpowerspectrum(freq, ipost, magflag)
     
     
     for( i=0; i<nfreqlog; i+=1)
         fft_ed_all[ipost][i] = fft_ed[i]
     endfor
     
    matrixop/O fft_ed = ln(fft_ed)
     
     if (maskwave[xpost][ypost] == 0)
         fft_ed_bg = fft_ed_bg+fft_ed
         Vred =0
         Vgreen = 0
         Vblue = 50000
         nbg+=1
     elseif(maskwave[xpost][ypost] == 1)
         fft_ed_cell = fft_ed_cell+fft_ed
         Vred = 50000
         Vgreen = 0
         Vblue = 0
         ncell+=1
     elseif(maskwave[xpost][ypost] == 4)
         Vred = 0
         Vgreen = 50000
         Vblue = 0
     endif
     
     if (ipost==0)
        display/N=$graphname fft_ed vs freqlogx
        modifygraph rgb(fft_ed) = (Vred, Vgreen, Vblue)
//        modifygraph log(bottom)=1
//        modifygraph log(left)=1
      else
         appendtograph/W= $graphname fft_ed_all[ipost][]/TN=$ffttracename vs freqlogx
         modifygraph rgb($ffttracename) = (Vred, Vgreen, Vblue)
      endif
     
     
     endif

endfor

fft_ed_bg = fft_ed_bg/nbg
fft_ed_cell = fft_ed_cell/ncell

for ( ipost = 0; ipost<npost; ipost+=1)
      xpost = floor(ipost/nposty)
      ypost = mod(ipost, nposty)

      fft_ed = fft_ed_all[ipost][p]
      matrixop/O fft_ed = ln(fft_ed)
      
      if (maskwave[xpost][ypost] == 0)      
      fft_ed_bg_sd = fft_ed_bg_sd + (fft_ed - fft_ed_bg)^2
      elseif(maskwave[xpost][ypost] == 1)      
      fft_ed_cell_sd = fft_ed_cell_sd + (fft_ed - fft_ed_cell)^2
      endif
endfor

fft_ed_bg_sd = sqrt(fft_ed_bg_sd/(nbg^2))
fft_ed_cell_sd = sqrt(fft_ed_cell_sd/(ncell^2))            


duplicate/o fft_ed_bg_sd fft_ed_bg_sdu
duplicate/o fft_ed_bg_sd fft_ed_bg_sdl


duplicate/o fft_ed_cell_sd fft_ed_cell_sdu
duplicate/o fft_ed_cell_sd fft_ed_cell_sdl

   

matrixop/o fft_ed_bg = exp(fft_ed_bg)
matrixop/o fft_ed_cell = exp(fft_ed_cell)

fft_ed_bg_sdu = (exp(fft_ed_bg_sd)-1)*fft_ed_bg
fft_ed_cell_sdu = (exp(fft_ed_cell_sd)-1)*fft_ed_cell
fft_ed_bg_sdl =  (1 - exp(-fft_ed_bg_sd))*fft_ed_bg
fft_ed_cell_sdl = (1 - exp(-fft_ed_cell_sd))*fft_ed_cell   


end
        
        
        


function calcccfunction(postid)
variable postid

wave freq


calcrefDLI(freq, postid)
calcpassivenoise(freq, postid)


wave magnoise
wave refDLItimedep
wave timewave

//setscale/P x, 0, 5,"" magnoise
//setscale/P x, 0, 5,"" refDLItimedep

variable timelength = dimsize(magnoise, 0)

variable i , ccfunction_V

variable magmean, magsd, refmean, refsd

variable cclength = floor(timelength/2)

make/O/N=(cclength-1) ccfunction, automag, autoref

     wavestats/Q magnoise
     magmean = V_avg
     magsd = V_sdev
     wavestats/Q refDLItimedep
     refmean = V_avg
     refsd = V_sdev

for ( i=0; i<cclength-1; i+=1)
     duplicate/o/R=(0, timelength-i-1) magnoise magtemp
     duplicate/o/R=(0, timelength-i-1) refDLItimedep reftemp
     duplicate/o/R= (i+1, timelength-1) refDLItimedep reftemp_s
     duplicate/o/R= (i+1, timelength-1)  magnoise magtemp_s
     make/O/N=(timelength - i  ) ccfunctionV
     ccfunctionV = (magtemp-magmean) * (reftemp_s-refmean)
     ccfunction_V = mean(ccfunctionV)
     ccfunction[i] = ccfunction_V/(magsd*refsd)
     ccfunctionV=(magtemp-magmean)*(magtemp_s-magmean) 
     ccfunction_V = mean(ccfunctionV)
     automag[i] = ccfunction_V/(magsd^2)
     ccfunctionV=(reftemp-refmean)*(reftemp_s-refmean) 
     ccfunction_V = mean(ccfunctionV)
     autoref[i] = ccfunction_V/(refsd^2)     
   
     
     
endfor

FFT/out=3/dest=mag_c_fft automag
FFT/out=3/dest=ref_c_fft autoref
FFT/out=4/dest=cc_fft ccfunction

variable cohsize = dimsize(cc_fft,0)

make/o/N=(cohsize) cohfunction

cohfunction = cc_fft/(mag_c_fft*ref_c_fft)




end



function calcccoherence(postid, freq)
variable postid

variable freq

string fitaffix, fitname, resname, histname

//do all frequency at ones
//wave freq
//calcrefDLI(freq, postid)
//calcpassivenoise(freq, postid)

SVAR Fitres_string

  if (freq <1)
         sprintf fitaffix "_0_%dhz"  (freq*10)
  else
         sprintf fitaffix "_%dhz"  freq
  endif


fitname = FitRes_String+"_ed"+fitaffix
sprintf resname "coh_%d_%d" (postid), (freq*10)
histname = resname + "_hist"

//fitname="fitres_ed"

make/o/n=18000 dataw

duplicatesinglepost($fitname,dataw, postid)

indpassivenoise($fitname, 1, postid,freq)
TimeVtest2(dataw, 7, 1)


//wave magnoise
//wave refDLItimedep
//wave timewave

//variable timelength = dimsize(timewave, 0)

//setscale/P x, 0, 5,"" magnoise
//setscale/P x, 0, 5,"" refDLItimedep

//FFT/OUT=1/DEST=magnoise_fft magnoise
//FFT/OUT=1/DEST= ref_fft refDLItimedep
wave noisetemp
wave DLItimeV


variable timelength=dimsize(noisetemp, 0)

FFT/RP=[5,timelength-6]/OUT=1/Dest=magnoise_fft noisetemp
FFT/RP=[5,timelength-6]/OUT=1/Dest=ref_fft DLItimeV
FFT/RP=[5,timelength-6]/OUT=3/Dest=magnoise_auto noisetemp
FFT/RP=[5,timelength-6]/OUT=3/Dest=ref_auto DLItimeV



variable fftlength = dimsize(magnoise_fft, 0)
variable i

Matrixop/O magnoise_FFt_c = conj(magnoise_fft)
matrixop/O ref_fft_c = conj(ref_fft)

make/o/c/n=(fftlength) cohfunc_c, cohfunc
//make/o/n=(fftlength) cohfunc

//duplicate/o magnoise magnoise_c
//duplicate/o refDLItimedep ref_c
//duplicate/o refDLItimedep ref_c2

//correlate magniose_c, ref_c2

//correlate/auto magnoise_c, magnoise_c
//correlate/auto  ref_c, ref_c

//variable corlength=dimsize(magnoise_c,0)
//make/o/n=(fftlength) cohfunc

//FFT/RP=[(corlength-1)/2,corlength-1]/out=3/dest=magnoise_c_fft  magnoise_c
//FFT/RP=[(corlength-1)/2,corlength-1]/out=3/dest=ref_c_fft  ref_c
//FFT/RP=[(corlength-1)/2,corlength-1]/out=4/dest=cor_fft  ref_c2

//FFT/RP=[0,(corlength-1)/2]/out=3/dest=magnoise_c_fft  magnoise_c
//FFT/RP=[0,(corlength-1)/2]/out=3/dest=ref_c_fft  ref_c

variable magnoise_sum, ref_sum

magnoise_sum = sum(magnoise_auto)
ref_sum = sum(ref_auto)


for ( i =0; i<fftlength; i+=1)
      cohfunc_c[i] = (magnoise_fft[i]*ref_fft_c[i])/(magnoise_auto[i]*ref_auto[i])
      //cohfunc_c[i] = (magnoise_fft[i]*ref_fft_c[i])^2/(magnoise_sum*ref_sum)
      //cohfunc[i] = cor_fft[i]/(magnoise_c_fft[i]*ref_c_fft[i])
      cohfunc[i] = r2polar(cohfunc_c[i])
endfor


make/o/n=(fftlength) restemp

restemp = imag(cohfunc)
duplicate/o restemp $resname

make/o/n=10 $histname

histogram/B={-pi, 2*pi/10, 10} $resname, $histname


end





function calcccfunction2(postid)
variable postid

wave freq


wave fitres_TLMRc1_ed_0_1hz

make/o/n=18000 dataw

duplicatesinglepost(fitres_TLMRc1_ed_0_1hz,dataw, postid)

indpassivenoise(fitres_TLMRc1_ed_0_1hz, 1, postid, 0.1)
TimeVtest2(dataw, 0.1, 1)


//setscale/P t, 0, 5,"" DLItimeV
//setscale/P t, 0, 5,"" noisetemp


//calcrefDLI(freq, postid)
//calcpassivenoise(freq, postid)


wave noisetemp
wave DLItimeV
wave timewave

//setscale/P x, 0, 5,"" magnoise
//setscale/P x, 0, 5,"" refDLItimedep

variable timelength = dimsize(noisetemp, 0)

variable i , ccfunction_V

variable magmean, magsd, refmean, refsd

variable cclength = floor(timelength/2)

make/O/N=(cclength) ccfunction, automag, autoref

     wavestats/Q noisetemp
     magmean = V_avg
     magsd = V_sdev
     wavestats/Q DLItimeV
     refmean = V_avg
     refsd = V_sdev

for ( i=0; i<cclength; i+=1)
     duplicate/o/R=(0, timelength-i-1) noisetemp magtemp
     duplicate/o/R=(0, timelength-i-1) DLItimeV reftemp
     duplicate/o/R= (i+1, timelength-1) DLItimeV reftemp_s
     duplicate/o/R= (i+1, timelength-1)  noisetemp magtemp_s
     make/O/N=(timelength - i  ) ccfunctionV
     ccfunctionV = (magtemp-magmean) * (reftemp_s-refmean)
     ccfunction_V = mean(ccfunctionV)
     ccfunction[i] = ccfunction_V/(magsd*refsd)
     ccfunctionV=(magtemp-magmean)*(magtemp_s-magmean) 
     ccfunction_V = mean(ccfunctionV)
     automag[i] = ccfunction_V/(magsd^2)
     ccfunctionV=(reftemp-refmean)*(reftemp_s-refmean) 
     ccfunction_V = mean(ccfunctionV)
     autoref[i] = ccfunction_V/(refsd^2)     
   
     
     
endfor

FFT/out=3/dest=mag_c_fft automag
FFT/out=3/dest=ref_c_fft autoref
FFT/out=4/dest=cc_fft ccfunction

variable cohsize = dimsize(cc_fft,0)

make/o/N=(cohsize) cohfunction

cohfunction = cc_fft/(mag_c_fft*ref_c_fft)

end




   
     
//calculate coherence function between two different waves     
function calccoherence3(w1,w2)
wave w1, w2     

variable timelength=dimsize(w2, 0)
variable fftstart

if (mod(timelength,2)==1)
   fftstart=1
   else 
   fftstart=0
endif


FFT/RP=[fftstart,timelength-1]/OUT=1/Dest=w1_fft w1
FFT/RP=[fftstart,timelength-1]/OUT=1/Dest=w2_fft w2
FFT/RP=[fftstart,timelength-1]/OUT=3/Dest=w1_auto w1
FFT/RP=[fftstart,timelength-1]/OUT=3/Dest=w2_auto w2



variable fftlength = dimsize(w1_fft, 0)
variable i

Matrixop/O w1_FFt_c = conj(w1_fft)
matrixop/O w2_fft_c = conj(w2_fft)

make/o/c/n=(fftlength) cohfunc_c, cohfunc
//make/o/n=(fftlength) cohfunc


variable w1_sum,w2_sum

w1_sum = sum(w1_fft)
w2_sum = sum(w2_fft)


for ( i =0; i<fftlength; i+=1)
      cohfunc_c[i] = (w1_fft[i]*w2_fft_c[i])/(w1_auto[i]*w2_auto[i])
      //cohfunc_c[i] = (w1_fft[i]*w2_fft_c[i])/(w1_sum*w2_sum)
      //cohfunc[i] = cor_fft[i]/(magnoise_c_fft[i]*ref_c_fft[i])
      cohfunc[i] = r2polar(cohfunc_c[i])
endfor



end       


function calcbgdisp(astr, dataw)
wave astr, dataw
variable postnum = astr[3]*astr[4]
variable i 
variable avgrange=300
variable movielength = dimsize(dataw, 2)
make/o/n=(postnum, 2)   disp
make/o/n=(movielength) disptemp

for(i=0; i<postnum; i+=1)
    disptemp = dataw[i][0][p]
    disp[i][0] = mean(disptemp,movielength-avgrange-1, movielength-1) - mean(disptemp, 0,avgrange) 
    disptemp = dataw[i][1][p]
    disp[i][1] = mean(disptemp,movielength-avgrange-1, movielength-1) - mean(disptemp, 0,avgrange) 
    endfor
    
end
    
    



function  dispbgvector(astr, disp)

wave astr, disp

variable   i
wave     limaskval
variable dx, dy, length
variable totalnumpost = astr[3]*astr[4]
variable Mfact = 500
duplicate/o linmaskval cellMatrix
make/o/n=(totalnumpost,2) lengthandangle, bluelengthandangle
make/o/n=1 origx, origy
make/o/n=(1,2) refvector
origx = 0
origy = 0

lengthandangle = 0
bluelengthandangle = 0
 
      DoWindow/K bgvector
	NewImage/K=1/N=bgvector/S=1 M_rotatedimage         
	
	wave linmaskx, linmasky
	
	
duplicate/o linmasky cellmasky
duplicate/o linmaskx cellmaskx

cellmaskx = 0
cellmasky = 0	


for (i=0; i<(totalnumpost-1); i+=1 )	
     
      dx =disp[i][0]
      dy =disp[i][1]
			        
       length = sqrt(dx^2+dy^2)
	 
	 if( cellmatrix[i] == kcell)
	     cellmaskx[i] = linmaskx[i]
	     cellmasky[i] = linmasky[i]
	     lengthandangle[i][0] = 0 //Mfact*length
	     lengthandangle[i][1] = 0 //-angle(dx,dy)
	  elseif( cellmatrix[i] == kempty)
	     bluelengthandangle[i][0] = Mfact*length
	     bluelengthandangle[i][1] = -angle(dx,dy)
	  endif
	         		
                  
endfor			        

refvector[0][0] = Mfact*0.1
refvector[0][1] = 0



appendtograph/T linmasky vs linmaskx	
appendtograph/T cellmasky vs cellmaskx
appendtograph/T origy vs origx

TextBox/C/N=text1/A=LT/B=1/G=(65200,65200,65200)/X=1.20/Y=1.50/T=20 num2str(0.1)+"pixel"

modifygraph arrowmarker(linmasky)	= {bluelengthandangle,1,3,2,0}
modifygraph rgb(linmasky) = (0,0,65535)

modifygraph arrowmarker(cellmasky) = {lengthandangle,1,3,2,0}
modifygraph rgb(cellmasky) = (65535,0,0)

modifygraph arrowmarker(origy)={refvector, 3,3,2,0}
modifygraph rgb(origy) = (0,0,65535)
modifygraph mode=3, marker=19, msize=1

//if(cellflag==1) 
//appendtograph/T
end        			



function calcbgrotate(zeropos, dataw)
wave zeropos, dataw

variable postnum = dimsize(dataw,0)
variable movielength = dimsize(dataw,2)

variable i, iframe
variable centx, centy
variable dx, dy

wave cellmatrix, fitres_ed

make/o/n=(postnum, movielength) bgrotate
make/o/n=(movielength) bgangle
make/o/n=(postnum) bgrotate0, bgrotate1
bgrotate=0

make/o/n=(postnum,2,movielength) fitres_final, bgdist

for(iframe = 0; iframe<movielength; iframe+=1)

centx = 0
centy = 0

      for(i=0; i<postnum; i+=1)
          centx+=zeropos[i][0][iframe]
          centy+=zeropos[i][1][iframe]
          
          endfor
          
          centx = centx/postnum
          centy = centy/postnum
          
      for(i=0; i<postnum; i+=1)
         if (cellmatrix[i] ==kempty)
            dx = centx - dataw[i][0][iframe] 
            dy = centy - dataw[i][1][iframe]
            bgdist[i][iframe][0] = dx
            bgdist[i][iframe][1] = dy
            bgrotate[i][iframe] = -angle(dx, dy)  
          endif  
      endfor    
      
      
      if (iframe==0)
          bgrotate0 = bgrotate[p][iframe]     
      endif     
          bgrotate1 = bgrotate[p][iframe]     
                 
                 
        if(iframe==0)
           bgangle[iframe] =0 
        else   
           bgangle[iframe] = mean(bgrotate1, 0, postnum/2) - mean(bgrotate0,0, postnum/2)
        endif
          
                   
 endfor
 
 
 for ( iframe =0; iframe<movielength; iframe+=1)
 
       for( i =0; i<postnum; i+=1) 
            fitres_final[i][0][iframe] = fitres_ed[i][0][iframe] + bgdist[i][iframe][1]*bgangle[iframe]
            fitres_final[i][1][iframe] = fitres_ed[i][1][iframe] - bgdist[i][iframe][0]*bgangle[iframe]
        endfor
         
endfor

 
               
 

end          


function findroicenter(astr, dataw)
wave astr, dataw


wave zeropos
variable iframe
variable movielength = dimsize(dataw, 2)
variable postnum = astr[3]*astr[4]
variable centx, centy, dx, dy
variable i

make/o/n=(postnum, movielength,2)  postdist

for(iframe = 0; iframe<movielength; iframe+=1)

centx = 0
centy = 0

      for(i=0; i<postnum; i+=1)
          centx+=zeropos[i][0][iframe]
          centy+=zeropos[i][1][iframe]
          
          endfor
          
          centx = centx/postnum
          centy = centy/postnum
          
       for(i=0; i<postnum; i+=1)
            dx = centx - dataw[i][0][iframe] 
            dy = centy - dataw[i][1][iframe]
            postdist[i][iframe][0] = dx
            postdist[i][iframe][1] = dy 
      endfor  
      
      
      endfor
      
      
      end
      
      


function    calcdialation(astr, dataw)
wave astr, dataw
wave postdist
wave fitresults

variable movielength = dimsize(dataw, 2)
variable postnum = dimsize(dataw, 0)

variable bgcount

wave cellmatrix

findroicenter(astr, fitresults)

make/o/n=(postnum, movielength)  dialfactor          //dialation factor
make/o/n=(movielength)  dialfactoravg=0     // dialation factor averaged over post numbers

variable iframe
variable i



for (iframe = 0; iframe<movielength; iframe+=1)
     bgcount = 0 
     for (i=0; i<postnum; i+=1)
          
          dialfactor[i][iframe] = sqrt(postdist[i][iframe][0]^2+ postdist[i][iframe][1]^2)/sqrt(postdist[i][0][0]^2+ postdist[i][0][1]^2)
          
          if ( cellmatrix[i] == kempty )
              bgcount+=1
              dialfactoravg[iframe] = dialfactoravg[iframe] + dialfactor[i][iframe]
          endif
      
      endfor
      
      dialfactoravg[iframe] = dialfactoravg[iframe]/bgcount 

endfor


end
             
               

function   removebgdial(astr, dataw)
wave astr, dataw

duplicate/o dataw fitres_rd

wave postdist
wave dialfactoravg

variable movielength = dimsize(dataw, 2)
variable postnum = dimsize(dataw, 0)

variable i , iframe

calcdialation(astr, dataw)

for ( iframe = 0 ; iframe<movielength; iframe+=1)
     for( i =0; i<postnum; i+=1)
          fitres_rd[i][0][iframe] = dataw[i][0][iframe] + (dialfactoravg[iframe]-1)*postdist[i][0][0]
          fitres_rd[i][1][iframe] = dataw[i][1][iframe] + (dialfactoravg[iframe]-1)*postdist[i][0][1]
     endfor
endfor


end     
          
                    


function calcavgfft(maskwave, astr)
wave maskwave, astr
wave fft_ed_cell, fft_ed_bg, fft_ed_cell_sd,fft_ed_bg_sd, fft_ed_all

variable fftlength = dimsize(fft_ed_cell,0)

duplicate/o fft_ed_all fft_ed_all_avg

make/o/n=(fftlength)  fft_cell_avg=0, fft_bg_avg=0, fft_diff_avg=0, fft_diff_avg_u=0, fft_diff_avg_l=0, fft_bg_avg_u=0, fft_bg_avg_l=0, fft_cell_avg_u=0, fft_cell_avg_l=0
make/o/n=(fftlength)  fft_cell_avg_sd=0, fft_bg_avg_sd=0


plotallfft(maskwave, astr, 0.1)
fft_cell_avg+=fft_ed_cell
fft_bg_avg +=fft_ed_bg
fft_ed_all_avg += fft_ed_all
fft_cell_avg_sd += (fft_ed_cell_sd)^2
fft_bg_avg_sd += (fft_ed_bg_sd)^2

plotallfft(maskwave, astr, 0.2)
fft_cell_avg+=fft_ed_cell
fft_bg_avg +=fft_ed_bg
fft_ed_all_avg += fft_ed_all
fft_cell_avg_sd += (fft_ed_cell_sd)^2
fft_bg_avg_sd += (fft_ed_bg_sd)^2

plotallfft(maskwave, astr, 0.5)
fft_cell_avg+=fft_ed_cell
fft_bg_avg +=fft_ed_bg
fft_ed_all_avg += fft_ed_all
fft_cell_avg_sd += (fft_ed_cell_sd)^2
fft_bg_avg_sd += (fft_ed_bg_sd)^2
          
plotallfft(maskwave, astr, 0.8)
fft_cell_avg+=fft_ed_cell
fft_bg_avg +=fft_ed_bg
fft_ed_all_avg += fft_ed_all
fft_cell_avg_sd += (fft_ed_cell_sd)^2
fft_bg_avg_sd += (fft_ed_bg_sd)^2

fft_cell_avg= fft_cell_avg/4
fft_bg_avg= fft_bg_avg/4
fft_ed_all_avg = fft_ed_all_avg/4
fft_cell_avg_sd = sqrt(fft_cell_avg_sd)/4            
fft_bg_avg_sd = sqrt(fft_bg_avg_sd)/4               

fft_cell_avg_u = (exp(fft_cell_avg_sd)-1)*fft_cell_avg
fft_cell_avg_l = (1-exp(-fft_cell_avg_sd))*fft_cell_avg
fft_bg_avg_u = (exp(fft_bg_avg_sd)-1)*fft_bg_avg
fft_bg_avg_l = (1-exp(-fft_bg_avg_sd))*fft_bg_avg


  
fft_diff_avg = fft_cell_avg - fft_bg_avg
fft_diff_avg_u = sqrt((fft_cell_avg_u)^2 + (fft_bg_avg_u)^2)
fft_diff_avg_l = sqrt((fft_cell_avg_l)^2 + (fft_bg_avg_l)^2)
end


function editfreqlogx()
variable i
variable nopfreq = 30
wave freqlogx
duplicate/o freqlogx freqlogx_ed
for (i = 18; i < (18+nopfreq + 1); i+=1)
     freqlogx_ed[i] = 10^((log(freqlogx[i])+log(freqlogx[i-1]))/1.9)
endfor
end

//mag_all use fitted value at 0.5 Hz

function calcstrainstats(maskwave,astr, straindata)


wave maskwave,astr, straindata

variable i, j 
variable nxpost = astr[3]
variable nypost = astr[4]
variable xindex, yindex

variable npost, wavelength
wave freqlogx

wave fft_ed_bg

npost = dimsize(straindata, 0)
wavelength = dimsize(straindata, 1)

make/o/n=(npost) slope_all=0, mag_all=0, slope_diff_all=0, mag_diff_all=0


make/o/n=(wavelength) singlestrainlog, singlestrainlog_diff
make/o/n=(wavelength)  freqlogxlog
freqlogxlog = log(freqlogx)

wave W_coef
//fit through minum to 1hz
for ( i = 0; i<npost; i+=1)
      xindex = floor(i/nypost)
      yindex = mod(i,nypost)
      if (maskwave[xindex][yindex] !=3 )
      singlestrainlog = log(straindata[i][p])  
      singlestrainlog_diff = log(straindata[i][p] - fft_ed_bg[p])
      curvefit/NTHR=0/Q line singlestrainlog[5,18] /X=freqlogxlog/D
     slope_all[i] = W_coef[1]
// //     mag_all[i] = straindata[i][0]
      mag_all[i] = 10^(W_coef[0] + W_coef[1]*(log(0.5)))
      curvefit/NTHR=0/Q line singlestrainlog_diff[0,20] /X=freqlogxlog/D
      slope_diff_all[i] = W_coef[1]
      mag_diff_all[i] = 10^(W_coef[0] + W_coef[1]*(log(0.5)))


//       singlestrainlog = straindata[i][p]
//       Make/O/T/N=2 T_Constraints
//       T_Constraints[0] = {"K0 > 0","K1 > 0"}
//       curvefit/NTHR=0/Q/W=2 power singlestrainlog[0,29] /X=freqlogx/D/C=T_constraints
//       slope_all[i] = W_coef[2]
//       mag_all[i] = W_coef[0] +W_coef[1]*0.5^W_coef[2]

      endif

endfor

duplicate M_rotatedimage M_rotatedimage_strain

dowindow/K strain_heatmap
newimage/N=strain_heatmap M_rotatedimage_strain
appendtograph/T/W=strain_heatmap linmasky vs linmaskx
modifygraph mode =3, marker =0
modifygraph zcolor(linmaskY)={mag_all,*,*,Rainbow, 0}

end
      
      


function sortasbycontra(maskwave, astr,contradata, straindata, bithresh) 
wave maskwave, contradata, astr,straindata
variable bithresh

wave slope_all, mag_all

variable i,j 
variable nob =10   //number of bins
variable npost = dimsize(slope_all,0)
make/o/n=(nob) strainslopebins=0, strainmagbins=0, contrabins=0, countbins = 0
make/o/n=(npost) contrawave
contrawave = contradata[p][0]
variable V_max = 0
variable wavelength = dimsize(straindata,1)
variable nxpost = astr[3]
variable nypost = astr[4]
variable xindex ,yindex

for (i = 1; i<npost; i+=1)
             xindex = floor(i/nypost)
             yindex = mod(i,nypost)
      if (maskwave[xindex][yindex] ==1)       
          if (contrawave[i] > V_max)
              V_max = contrawave[i]
          endif
      endif
endfor

//temporary test
V_max = 11

variable contrabingap = V_max*1.01/nob
//variable bithresh = V_max*0.5


variable counthb=0, countlb=0
make/o/n=(wavelength)  fft_hb_ed=0, fft_lb_ed=0
make/o/n=(wavelength) fft_hb_ed_sd = 0, fft_lb_ed_sd=0
contrabins = x*contrabingap

for (i=1; i<npost; i+=1)   
             xindex = floor(i/nypost)
             yindex = mod(i,nypost)
      if (maskwave[xindex][yindex] ==1 ) 
            j = floor(contrawave[i]/contrabingap)
           strainslopebins[j]+=slope_all[i]   
           strainmagbins[j]+=log(mag_all[i])
           countbins[j]+=1
           if ( contrawave[i] >bithresh)
                fft_hb_ed += log(straindata[i][p])
                counthb +=1
           else
                fft_lb_ed += log(straindata[i][p])
                countlb +=1
            endif    
       endif
 endfor
 
 duplicate/o countbins countbin
 
 for(j=0; j<nob; j+=1)
      if (countbins[j]==0)
           countbins[j]=1
       endif
endfor           
 
 strainslopebins = strainslopebins/countbins        
 strainmagbins = strainmagbins/countbins
 strainmagbins = 10^strainmagbins
 
 fft_hb_ed = 10^(fft_hb_ed/counthb)
 fft_lb_ed = 10^(fft_lb_ed/countlb)
 
  //calculating the error associated with highly bent and lowly bent posts
for (i=1; i<npost; i+=1)   
             xindex = floor(i/nypost)
             yindex = mod(i,nypost)
      if (maskwave[xindex][yindex] ==1 ) 
           if ( contrawave[i] >bithresh)
                fft_hb_ed_sd += (log(straindata[i][p]) - log(fft_hb_ed))^2
                
           else
                fft_lb_ed_sd += (log(straindata[i][p]) - log(fft_lb_ed))^2
                
            endif    
       endif
 endfor 
 

fft_hb_ed_sd =  sqrt(fft_hb_ed_sd/ counthb^2)
fft_lb_ed_sd = sqrt(fft_lb_ed_sd/countlb^2)

duplicate/o fft_hb_ed_sd fft_hb_ed_sdu, fft_hb_ed_sdl
duplicate/o fft_lb_ed_sd fft_lb_ed_sdu, fft_lb_ed_sdl

fft_hb_ed_sdu = (10^(fft_hb_ed_sd)-1) * fft_hb_ed
fft_hb_ed_sdl = fft_hb_ed *( 1 - 10^(-1*fft_hb_ed_sd))

fft_lb_ed_sdu = (10^(fft_lb_ed_sd)-1) *fft_lb_ed
fft_lb_ed_sdl = fft_lb_ed *( 1 - 10^(-1*fft_lb_ed_sd))

 

 
 
 
 wave fft_ed_bg, freqlogx, fft_ed_bg_sdu, fft_ed_bg_sdl
 
 duplicate/o fft_hb_ed fft_hb_ed_diff, fft_hb_ed_sdu_diff, fft_hb_ed_sdl_diff
 fft_hb_ed_diff = fft_hb_ed - fft_ed_bg
 fft_hb_ed_sdu_diff = sqrt(fft_hb_ed_sdu^2 +  fft_ed_bg_sdu^2)
 fft_hb_ed_sdl_diff = sqrt(fft_hb_ed_sdl^2 +  fft_ed_bg_sdl^2)
 duplicate/o fft_lb_ed  fft_lb_ed_diff, fft_lb_ed_sdu_diff, fft_lb_ed_sdl_diff
 fft_lb_ed_diff = fft_lb_ed - fft_ed_bg
 fft_lb_ed_sdu_diff = sqrt(fft_lb_ed_sdu^2 +  fft_ed_bg_sdu^2) 
 fft_lb_ed_sdl_diff = sqrt(fft_lb_ed_sdl^2 +  fft_ed_bg_sdl^2)
 
 duplicate/o fft_hb_ed_diff fft_hb_ed_diff2
 duplicate/o fft_lb_ed_diff fft_lb_ed_diff2
 
 fft_hb_ed_diff2 = fft_hb_ed_diff - mean(fft_hb_ed_diff, wavelength-6, wavelength-2)
 fft_lb_ed_diff2 = fft_lb_ed_diff - mean(fft_lb_ed_diff, wavelength-6, wavelength-2)
 
 
 
 CurveFit/NTHR=0 Power  fft_hb_ed_diff2[0,20] /X=freqlogx /W=fft_lb_ed_sdu_diff /I=1 /D 
 CurveFit/NTHR=0 Power  fft_lb_ed_diff2[0,20] /X=freqlogx /W=fft_lb_ed_sdu_diff /I=1 /D
 
 
 
 duplicate/o fft_hb_ed modu
 modu = 49.7*freqlogx^0.159
 
 duplicate/o fft_hb_ed st_hb_ed
 duplicate/o fft_lb_ed st_lb_ed
 st_hb_ed = fft_hb_ed * modu^2
st_lb_ed = fft_lb_ed * modu^2
 
duplicate/o fft_hb_ed st_hb_diff
duplicate/o fft_lb_ed st_lb_diff
st_hb_diff = fft_hb_ed_diff * modu^2
st_lb_diff = fft_lb_ed_diff * modu^2
 
 wave W_coef, freqlogxlog
 
make/o/n=2  hbvlb_slope
duplicate/o st_hb_diff st_hb_ed_log
st_hb_ed_log = log(st_hb_diff)
CurveFit/Q/NTHR=0 line  st_hb_ed_log[0,20] /X=freqlogxlog /D 
hbvlb_slope[0] = W_coef[1]
duplicate/o st_lb_diff st_lb_ed_log
st_lb_ed_log = log(st_lb_diff)
CurveFit/Q/NTHR=0 line  st_lb_ed_log[0,20] /X=freqlogxlog /D 
hbvlb_slope[1] = W_coef[1]

duplicate/o slope_all slope_low, slope_bg
duplicate/o slope_diff_all slope_diff_low
duplicate/o mag_all mag_low, mag_bg
duplicate/o mag_diff_all mag_diff_low
for ( i = 0; i<npost; i+=1)
      if (contrawave[i] == 0)
          mag_low[i] =-1
          slope_low[i] = 1
          mag_diff_low[i] = -1
          slope_diff_low[i] = 1
       elseif (contrawave[i] > 0 && contrawave[i] <5)
          mag_bg[i] = -1
          slope_bg[i] = 1
       else
          mag_low[i] = -1
          slope_low[i] = 1
          mag_diff_low[i] = -1
          slope_diff_low[i] = 1          
          mag_bg[i] = -1
          slope_bg[i] = 1
       endif         
endfor




  
 end            
 
 

function fitpowerlaw(xwave, ywave)
wave xwave, ywave
duplicate/o ywave ywavelog
duplicate/o xwave xwavelog
ywavelog = log(ywave)
xwavelog = log(xwave)

curvefit/nthr=0 line ywavelog[17,33] /x=xwavelog /D

end 




function calcmagstats(magstats, magstats_t, forcesum, forcesum_t)

wave magstats, magstats_t, forcesum, forcesum_t

variable npost = dimsize(magstats, 0)
variable npost_t = dimsize(magstats_t, 0)
variable ipost
variable nlowpost =0
variable nlowpost_t = 0
variable statforce
variable mean_N, mean_t, var_N, var_t
variable tvalue


make/o/n=(npost) mag_low_n
make/o/n=(npost_t) mag_low_t_n

for (ipost = 0 ; ipost<npost; ipost+=1)
       statforce = forcesum[ipost][0]
       if ( (statforce >0) && (statforce<2))
           mag_low_n[nlowpost] = magstats[ipost]
           nlowpost +=1
       endif
        
endfor

redimension/N=(nlowpost) mag_low_n

duplicate/o mag_low_n mag_low_log_n


for (ipost = 0; ipost<npost_t; ipost+=1)
      statforce = forcesum_t[ipost][0]
      if ( (statforce>0) && (statforce<2))
           mag_low_t_n[nlowpost_t] = magstats_t[ipost]
           nlowpost_t +=1
      endif
      
endfor

redimension/N=(nlowpost_t) mag_low_t_n

duplicate/o mag_low_t_n mag_low_log_t_n

mag_low_log_n = log(mag_low_n)
mag_low_log_t_n = log(mag_low_t_n)

wavestats/q mag_low_log_n
mean_N = V_avg
var_N = V_sdev^2
print mean_N, var_N


wavestats/q mag_low_log_t_n
mean_t = V_avg
var_t = V_sdev^2
print mean_t, var_t

tvalue = (mean_N - mean_t)/sqrt(var_N/nlowpost + var_t/nlowpost_t)

print tvalue, nlowpost, nlowpost_t
print studentA(tvalue, nlowpost+nlowpost_t-2)


end
           
           
//calculate time trace based on deflection           
function calcdeftrace(datawave, zerowave)
wave datawave, zerowave

duplicate/o datawave fitres_deflected
fitres_deflected = datawave - zerowave

end


function plottracemap(datawave,zerowave, maskwave, astr)

//datawave should be the deflected trace of each posts
wave datawave, zerowave
wave maskwave, astr


string mapWN = "tracemap"
string deftracename
string bondstring

duplicate/o datawave fitres_deflected
fitres_deflected = datawave - zerowave[p][q][0]
 
variable pixelratio=125    //pixel to nm
 
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost


wave bondid
variable ibond, nbond, bondx, bondy,p1,p2
nbond = dimsize(bondid, 0)

variable avgflag=1           //1 for applying averaging/smoothing
variable bgflag =1           //1 for plotting background posts
variable matrixconst = 1              //matrix distance between posts in pixels
variable boxcarsize = 10
variable shrinkratio = matrixconst/astr[5]
make/o/n=(npost,2) matrixbase 
variable ixpost, iypost, ipost


killwindow $mapWN

for (ixpost = 0; ixpost<nxpost; ixpost+=1)
       for (iypost =0; iypost<nypost; iypost+=1)
              ipost = ixpost*nypost+iypost
              matrixbase[ipost][0] = ixpost*matrixconst+mod(iypost,2)*0.5*matrixconst
              matrixbase[ipost][1] = iypost*matrixconst*0.5*sqrt(3)
        endfor
endfor
 
duplicate/o datawave fitres_deflected_ed
fitres_deflected_ed = fitres_deflected - (1-shrinkratio)*(fitres_deflected[p][q][0])

variable i
variable iflag=0

duplicate/o fitres_deflected_ed def_trace
def_trace = matrixbase[p][q]+fitres_deflected_ed


if (avgflag==1)
    variable movielength = dimsize(fitres_deflected,2)
    variable avglength =floor(movielength/boxcarsize)
    make/o/n=(npost, 2, avglength) def_trace_ed=0
    variable avgcount = 0, trace_ed_count=0
    make/o/n=(boxcarsize)  trace_segx, trace_segy
    for(ipost=0; ipost<npost; ipost+=1)
         trace_ed_count=0
         avgcount=0
         for(i=0; i<movielength; i+=1)
            trace_segx[avgcount] = def_trace[ipost][0][i]
            trace_segy[avgcount] = def_trace[ipost][1][i]
            avgcount +=1
            if(mod(avgcount,boxcarsize) == 0)
               avgcount=0
               def_trace_ed[ipost][0][trace_ed_count]=mean(trace_segx)
               def_trace_ed[ipost][1][trace_ed_count]=mean(trace_segy)
               trace_ed_count+=1
            endif
         endfor
      endfor      


duplicate/o def_trace_ed def_trace
endif

variable tracelength = dimsize(def_trace,2)
make/o/n=(tracelength) timetrace
timetrace = x

def_trace = def_trace*pixelratio

if ( bgflag==0)

for (i=0; i<npost; i+=1)
       ixpost = floor(i/nypost)
       iypost = mod(i,nypost)
   if(maskwave[ixpost][iypost] == 1)
      iflag+=1
      sprintf deftracename "p_%d" i
      if(iflag == 1 )
         display/N=$mapWN def_trace[i][1][] vs def_trace[i][0][]
         modifygraph zcolor(def_trace) = {timetrace,*,*,Rainbow,0}
       else 
          appendtograph/W=$mapWN def_trace[i][1][]/TN=$deftracename  vs def_trace[i][0][]
          modifygraph zcolor($deftracename) = {timetrace,*,*,Rainbow,0}
       endif
   endif  
endfor  

elseif (bgflag==1)

for (i=0; i<npost; i+=1)
       ixpost = floor(i/nypost)
       iypost = mod(i,nypost)
   if(maskwave[ixpost][iypost] != 3)
      iflag+=1
      sprintf deftracename "p_%d" i
      if(iflag == 1 )
         display/N=$mapWN def_trace[i][1][] vs def_trace[i][0][]
         modifygraph zcolor(def_trace) = {timetrace,*,*,Rainbow,0}
       else 
          appendtograph/W=$mapWN def_trace[i][1][]/TN=$deftracename  vs def_trace[i][0][]
          if (maskwave[ixpost][iypost]==0)
          modifygraph zcolor($deftracename) = {timetrace,*,*,Rainbow,0}
          else
          modifygraph zcolor($deftracename) = {timetrace,*,*,Rainbow,0}
          endif
       endif
   endif  
endfor         

endif



 
// modifygraph/W=$mapWN width = {aspect,1}
ModifyGraph/W=$mapWN width={perUnit,0.06,bottom},height={perUnit,0.06,left}
SetAxis/W=$mapWN/A/R left


for ( ibond=0; ibond<nbond; ibond+=1)
       sprintf bondstring, "%d", ibond
       p1 = bondid[ibond][0]
       p2 = bondid[ibond][1]
       bondx = 0.5*(matrixbase[p1][0] + matrixbase[p2][0])*pixelratio
       bondy = 0.5*(matrixbase[p1][1] + matrixbase[p2][1])*pixelratio
       SetDrawEnv/W=$mapWN xcoord= bottom,ycoord= left, fsize=6;DelayUpdate
       DrawText/W=$mapWN bondx, bondy, bondstring
endfor


end



function plottracemap_seg(datawave,zerowave, maskwave, astr, startframe)

//datawave should be the deflected trace of each posts
wave datawave, zerowave
wave maskwave, astr



variable startframe
variable seglength = 70

wave linmaskval_intercept_int_sc

string mapWN = "tracemap"
string deftracename
string bondstring

duplicate/o datawave fitres_deflected
fitres_deflected = datawave - zerowave[p][q][0]
 
variable pixelratio=125    //pixel to nm
 
variable nxpost = astr[3]
variable nypost = astr[4]
variable npost = nxpost*nypost


wave bondid
variable ibond, nbond, bondx, bondy,p1,p2
nbond = dimsize(bondid, 0)

variable avgflag=1           //1 for applying averaging/smoothing
variable bgflag =0           //1 for plotting background posts
variable matrixconst = 0.6             //matrix distance between posts in pixels
variable boxcarsize = 2
variable shrinkratio = matrixconst/astr[5]
make/o/n=(npost,2) matrixbase 
variable ixpost, iypost, ipost


killwindow $mapWN

for (ixpost = 0; ixpost<nxpost; ixpost+=1)
       for (iypost =0; iypost<nypost; iypost+=1)
              ipost = ixpost*nypost+iypost
              matrixbase[ipost][0] = ixpost*matrixconst+mod(iypost,2)*0.5*matrixconst
              matrixbase[ipost][1] = iypost*matrixconst*0.5*sqrt(3)
        endfor
endfor


 
duplicate/o datawave fitres_deflected_ed
fitres_deflected_ed = fitres_deflected - (1-shrinkratio)*(fitres_deflected[p][q][startframe])

variable i
variable iflag=0

duplicate/o fitres_deflected_ed def_trace
def_trace = matrixbase[p][q]+fitres_deflected_ed


if (avgflag==1)
    variable movielength = seglength
    variable avglength =floor(movielength/boxcarsize)
    make/o/n=(npost, 2, avglength) def_trace_ed=0
    variable avgcount = 0, trace_ed_count=0
    make/o/n=(boxcarsize)  trace_segx, trace_segy
    for(ipost=0; ipost<npost; ipost+=1)
         trace_ed_count=0
         avgcount=0
         for(i=startframe; i<startframe+seglength; i+=1)
            trace_segx[avgcount] = def_trace[ipost][0][i]
            trace_segy[avgcount] = def_trace[ipost][1][i]
            avgcount +=1
            if(mod(avgcount,boxcarsize) == 0)
               avgcount=0
               def_trace_ed[ipost][0][trace_ed_count]=mean(trace_segx)
               def_trace_ed[ipost][1][trace_ed_count]=mean(trace_segy)
               trace_ed_count+=1
            endif
         endfor
      endfor      


duplicate/o def_trace_ed def_trace
endif

variable tracelength = dimsize(def_trace,2)
make/o/n=(tracelength) timetrace
timetrace = x

def_trace = def_trace*pixelratio

if ( bgflag==0)

for (i=0; i<npost; i+=1)
       ixpost = floor(i/nypost)
       iypost = mod(i,nypost)
//   if(maskwave[ixpost][iypost] == 1)
   if(linmaskval_intercept_int_sc[i] == 1)
      iflag+=1
      sprintf deftracename "p_%d" i
      if(iflag == 1 )
         display/N=$mapWN def_trace[i][1][] vs def_trace[i][0][]
         modifygraph zcolor(def_trace) = {timetrace,*,*,Rainbow,0}
       else 
          appendtograph/W=$mapWN def_trace[i][1][]/TN=$deftracename  vs def_trace[i][0][]
          modifygraph zcolor($deftracename) = {timetrace,*,*,Rainbow,0}
       endif
    elseif(linmaskval_intercept_int_sc[i] ==  4)
      iflag+=1
      sprintf deftracename "p_%d" i     
      if(iflag == 1 )
         display/N=$mapWN def_trace[i][1][] vs def_trace[i][0][]
         modifygraph zcolor(def_trace) = {timetrace,*,*,Spectrum,0}
       else 
          appendtograph/W=$mapWN def_trace[i][1][]/TN=$deftracename  vs def_trace[i][0][]
          modifygraph zcolor($deftracename) = {timetrace,*,*,Spectrum,0}
       endif         
   endif  
endfor  

elseif (bgflag==1)

for (i=0; i<npost; i+=1)
       ixpost = floor(i/nypost)
       iypost = mod(i,nypost)
   if(maskwave[ixpost][iypost] != 3)
      iflag+=1
      sprintf deftracename "p_%d" i
      if(iflag == 1 )
         display/N=$mapWN def_trace[i][1][] vs def_trace[i][0][]
         modifygraph zcolor(def_trace) = {timetrace,*,*,Rainbow,0}
       else 
          appendtograph/W=$mapWN def_trace[i][1][]/TN=$deftracename  vs def_trace[i][0][]
          if (maskwave[ixpost][iypost]==0)
          modifygraph zcolor($deftracename) = {timetrace,*,*,Rainbow,0}
          else
          modifygraph zcolor($deftracename) = {timetrace,*,*,Rainbow,0}
          endif
       endif
   endif  
endfor         

endif



 
// modifygraph/W=$mapWN width = {aspect,1}
ModifyGraph/W=$mapWN width={perUnit,0.06,bottom},height={perUnit,0.06,left}
SetAxis/W=$mapWN/A/R left


//for ( ibond=0; ibond<nbond; ibond+=1)
//       sprintf bondstring, "%d", ibond
//       p1 = bondid[ibond][0]
//       p2 = bondid[ibond][1]
//       bondx = 0.5*(matrixbase[p1][0] + matrixbase[p2][0])*pixelratio
//       bondy = 0.5*(matrixbase[p1][1] + matrixbase[p2][1])*pixelratio
//       SetDrawEnv/W=$mapWN xcoord= bottom,ycoord= left, fsize=6;DelayUpdate
//       DrawText/W=$mapWN bondx, bondy, bondstring
//endfor


end


              

