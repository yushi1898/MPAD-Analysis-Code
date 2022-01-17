#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//new function is used to calculate the correlation between magnitude of force vector map and 
//difference of distribution btw rad and tangental direction

//1.2:   enable option to calculate histogram based on local radial and tangental motion
//use this to do individual post or overall statistic
function   drawhisto(astrw, datawave, bgshiftwave, maskwave, zerow, lagtime, bins,posttype)
wave astrw, datawave, bgshiftwave, maskwave
wave zerow   // stroe initial positions of posts
variable lagtime, bins, posttype

wave linmaskval_intercept_int, fitres_ed

//SVAR  rootstring
//SVAR  bgshiftstring, fitres_string
//variable  lagtime=1
variable npost, ipost
variable numb, binw

variable nxpost, nypost

variable framerate = 10

variable laginterval = lagtime * framerate

variable fixlength = 5000

variable nof = dimsize(datawave, 2)

variable xdisp, ydisp

//wave xdisptemp, ydistemp, disp

make/O/N=(nof) xdisptemp, ydisptemp

variable radsdev, tansdev
variable posti, postj


numb = 50
binw = 2*abs(bins)/numb

npost = dimsize(datawave, 0)
nxpost = astrw[3]
nypost = astrw[4]

variable CX, CY
make/O/N=(npost, 2)   directw



//make/O/(N=numb)     temphistorad, temphistotan   //used to instore histogram distribution of a single post

make/O/N=(npost)      forcevecmag = 0,  disdif = 0, disrad=0, distan=0  //magnitude of force vector and distribution difference waves

variable ncell=0, nbg=0, ncenter=0, nedge=0


for(ipost=0; ipost<npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)
         if(maskwave[posti][postj]==1) 	
            CX+=fitres_ed[ipost][0][0]
            CY+=fitres_ed[ipost][1][0]
            ncell+=1   
         endif
endfor       


CX=CX/ncell
CY=CY/ncell    

ncell = 0


for(ipost=0; ipost<npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)
       
   if (linmaskval_intercept_int[ipost] !=3 )    
//     xdisp = datawave[ipost][0][0] - zerow[ipost][0][0]
//     ydisp = datawave[ipost][1][0] - zerow[ipost][1][0]


     xdisptemp = fitres_ed[ipost][0][p]
     ydisptemp =fitres_ed[ipost][1][p]
     
     xdisp = CX - mean(xdisptemp)
     ydisp = CY - mean(ydisptemp) 

     directw[ipost][0] = xdisp/sqrt(xdisp^2+ydisp^2)
     directw[ipost][1] = ydisp/sqrt(xdisp^2+ydisp^2)
     

     //added  to calculate the magnitude of force vector
     forcevecmag[ipost] = sqrt(xdisp^2+ydisp^2)
   endif
   
//   if(linmaskval_intercept_int[ipost] ==1)
   if(linmaskval_intercept_int[ipost] ==posttype)
      ncell+=1
   endif   
    
endfor     

//loadwave/P=topractice/N=datawave  
//duplicate/O  $fitres_string datawave
//duplicate/O  $bgshiftstring  bgshiftwave

variable i, j,DX, DY
variable dataR

variable icell=0




DataR=dimsize(datawave,2)

variable intervalnum, intervalnum_unwind                        //number of points in histogram bins

intervalnum_unwind =floor( DataR/laginterval )
intervalnum = DataR- laginterval
//intervalnum = 5000


//make/O/N=(DataR-laginterval) Histo
//make/O/N= fixlength  Histo

variable npavg = 1           //number of posts put in histosinglepost 
variable lseg = 18000
//variable nseg = floor(datar/lseg)
variable nseg = 1 
variable k 

make/O/N = (intervalnum) Histo_x, Histo_Y, histo_t, histo_r, Histo_d,Histo_x_sq
make/O/N = (numb, npost,4) singlehisto
make/O/N = (numb)  bghisto = 0, cellhisto=0, edgehisto=0, centerhisto=0
make/o/n=(intervalnum_unwind*ncell) histosingle_x, histosingle_y, histosingle_r, histosingle_t
make/o/n=(intervalnum_unwind*npavg)  histosinglepost_x, histosinglepost_y
make/o/n=(ncell)   cagesize,Lmax_msd, alpha_msd
make/o/n=(nseg*ncell) lmax_cell
wave msd_slope_max


ncell = 0
 
wave maskflag, msd

variable x4m, x2m, x4m_avg, x2m_avg

make/o/n=(floor(npost/npavg)) ngpara

j=0

for (ipost = 0 ; ipost <  npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)
//    if (linmaskval_intercept_int[ipost] == 1) 
   if(linmaskval_intercept_int[ipost] ==posttype)
//       For(i=0; i < intervalnum; i+=laginterval)


       for(i=0; i<intervalnum; i+=1)
       
       
       //was used for removing large deflected fit
//            if ((maskflag[ipost][i] == 1) || (maskflag[ipost][i+laginterval] ==1))
//               histo[i] = 100
//            else   
       
//       for(i=0; i<fixlength; i+=1)
////            Histo_x[i]=(datawave[ipost][0][i+laginterval]-bgshiftwave[i+laginterval][0]) - (datawave[ipost][0][i]-bgshiftwave[i][0])
////            Histo_y[i]=(datawave[ipost][1][i+laginterval]-bgshiftwave[i+laginterval][1]) - (datawave[ipost][1][i]-bgshiftwave[i][1])
    //      DY=(datawave[npost][4][i+lagtime]-bgshiftwave[i+lagtime][1]) - (datawave[npost][4][i]-bgshiftwave[i][1])
//           endif

//           histosingle_x[i+intervalnum*icell] = fitres_ed[ipost][0][i+laginterval] - fitres_ed[ipost][0][i]
//           histosingle_y[i+intervalnum*icell] = fitres_ed[ipost][1][i+laginterval] - fitres_ed[ipost][1][i]
//           histosinglepost_x[i] = fitres_ed[ipost][0][i+laginterval] - fitres_ed[ipost][0][i]
//           histosinglepost_y[i] =  fitres_ed[ipost][1][i+laginterval] - fitres_ed[ipost][1][i]
           
            Histo_x[i]=fitres_ed[ipost][0][i+laginterval] - fitres_ed[ipost][0][i]
            Histo_y[i]=fitres_ed[ipost][1][i+laginterval] - fitres_ed[ipost][1][i]      
           Histo_r[i] = histo_x[i]*directw[ipost][0] + histo_y[i]*directw[ipost][1]
           histo_t[i] = histo_x[i]*directw[ipost][1] - histo_y[i] *directw[ipost][0]            
            Histo_d[i] = sqrt(Histo_x[i]^2+Histo_y[i]^2)
            Histo_x_sq[i] = (Histo_x[i])^2 + Histo_y[i]^2

        Endfor
        

        for(k=0; k<nseg;k+=1)
//                wavestats/q/r=(lseg*k, lseg*k+lseg) histo_x
             wavestats/q histo_d
//               wavestats/q histo_y
            Lmax_cell[icell*nseg+k] =125* max(v_max, abs(v_min))
        endfor    
        alpha_msd[icell] = msd_slope_max[ipost]
        Lmax_msd[icell] = 125*sqrt(mean(Histo_x_sq))
//        Lmax_msd[icell] = sqrt(10^6*msd[floor((laginterval-100)/10)+99][ipost])             //changed on 11/3/2018
//          Lmax_msd[icell] = msd[laginterval][ipost]
//          wavestats/q histo_x_sq
//          Lmax_msd[icell] = sqrt(v_avg)
//        Lmax_cell[icell] = v_max
        
        
        
        for (i=0; i<intervalnum_unwind; i+=1)
           histosingle_x[i+intervalnum_unwind*icell] = fitres_ed[ipost][0][i*laginterval+laginterval] - fitres_ed[ipost][0][i*laginterval]
           histosingle_y[i+intervalnum_unwind*icell] = fitres_ed[ipost][1][i*laginterval+laginterval] - fitres_ed[ipost][1][i*laginterval]
           histosingle_r[i+intervalnum_unwind*icell] =  histosingle_x[i+intervalnum_unwind*icell]*directw[ipost][0] + histosingle_y[i+intervalnum_unwind*icell]*directw[ipost][1]
           histosingle_t[i+intervalnum_unwind*icell] = histosingle_x[i+intervalnum_unwind*icell]*directw[ipost][1] - histosingle_y[i+intervalnum_unwind*icell]*directw[ipost][0]               
           histosinglepost_x[i+j*intervalnum_unwind] = fitres_ed[ipost][0][i*laginterval+laginterval] - fitres_ed[ipost][0][i*laginterval]
           histosinglepost_y[i+(j)*intervalnum_unwind] =  fitres_ed[ipost][1][i*laginterval+laginterval] - fitres_ed[ipost][1][i*laginterval]                        
      endfor           
      
        

        if(j==(npavg-1))
        
         wavestats/q  histosinglepost_x
       
         x4m_avg=0
         x2m_avg=0

        cagesize[icell] = max(abs(v_min),abs(v_max))


        for(i=0; i<(intervalnum_unwind*npavg); i+=1)
//          for(i=0; i<intervalnum; i+=1)
               x4m = (histosinglepost_x[i]-V_avg)^4
               x2m = (histosinglepost_x[i]-V_avg)^2
               x4m_avg+= x4m
               x2m_avg+= x2m
        endfor

        x4m_avg = x4m_avg/(intervalnum_unwind*npavg)
        x2m_avg = x2m_avg/(intervalnum_unwind*npavg)


//        x4m_avg = x4m_avg/intervalnum
//        x2m_avg = x2m_avg/intervalnum

        ngpara[icell] = x4m_avg/(3*(x2m_avg^2)) - 1        
        
        j=-1
 
        icell +=1
                
        endif
        
        j+=1
        

      
      wavestats/Q   histo_r
      radsdev = V_sdev
      wavestats/Q  histo_t
      tansdev = V_sdev
      
      disdif[ipost] = radsdev - tansdev
      disrad[ipost] = radsdev
      distan[ipost] = tansdev
      
      
       Make/N=100/O Histox_Hist, Histoy_Hist, histor_HIst, Histot_hist;DelayUpdate
       histogram/B={bins, binw, numb} Histo_x,Histox_Hist
       histogram/B={bins, binw, numb} Histo_y,Histoy_Hist
       histogram/B={bins, binw, numb} Histo_r,Histor_Hist
       histogram/B={bins, binw, numb} Histo_t,Histot_Hist                    
       make/N=(numb)/O Histox
       
       histox=bins+x*binw
       
       singlehisto[][ipost][2] = histox_hist[p]
       singlehisto[][ipost][3] = histoy_hist[p]
       singlehisto[][ipost][0] = histor_hist[p]
       singlehisto[][ipost][1] = histot_hist[p]
       
endif       
        
       if (maskwave[posti][postj] == 0)
           bghisto = bghisto + histox_hist
           nbg += 1
       elseif (maskwave[posti][postj] == 1)    
           cellhisto = cellhisto + histox_hist
           centerhisto = centerhisto + histox_hist
           ncell+=1
           ncenter+=1
       elseif (maskwave[posti][postj] == 6)
           cellhisto = cellhisto + histox_hist
           edgehisto = edgehisto + histox_hist
           ncell+=1
           nedge+=1
       endif
       

endfor      
 
//bghisto = bghisto/(nbg*intervalnum)
//cellhisto = cellhisto/(ncell*intervalnum)
//centerhisto = centerhisto/(ncenter * intervalnum)
//edgehisto = edgehisto/(nedge * intervalnum)

bghisto = bghisto/nbg
cellhisto = cellhisto/ncell
centerhisto = centerhisto/ncenter 
edgehisto = edgehisto/nedge
      
      
redimension/n=(icell)   ngpara      


end


function   drawhisto_shortV(astrw, datawave,maskwave, zerow, lagtime, bins,posttype)
wave astrw, datawave,maskwave
wave zerow   // stroe initial positions of posts
variable lagtime, bins,posttype

wave linmaskval_intercept_int, fitres_ed_s

//SVAR  rootstring
//SVAR  bgshiftstring, fitres_string
//variable  lagtime=1
variable npost, ipost
variable numb, binw

variable nxpost, nypost

//variable framerate = 10           //100 fps
variable framerate = 1                // 10fps

variable laginterval = lagtime * framerate

variable fixlength = 5000

variable nof = dimsize(datawave, 2)

variable xdisp, ydisp

//wave xdisptemp, ydistemp, disp

make/O/N=(nof) xdisptemp, ydisptemp

variable radsdev, tansdev
variable posti, postj


numb = 50
binw = 2*abs(bins)/numb

npost = dimsize(datawave, 0)
nxpost = astrw[3]
nypost = astrw[4]

make/O/N=(npost, 2)   directw


//make/O/(N=numb)     temphistorad, temphistotan   //used to instore histogram distribution of a single post

make/O/N=(npost)      forcevecmag = 0,  disdif = 0, disrad=0, distan=0  //magnitude of force vector and distribution difference waves

variable ncell=0, nbg=0, ncenter=0, nedge=0

for(ipost=0; ipost<npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)
       
   if (linmaskval_intercept_int[ipost] !=3 )    
//     xdisp = datawave[ipost][0][0] - zerow[ipost][0][0]
//     ydisp = datawave[ipost][1][0] - zerow[ipost][1][0]
     
     xdisptemp = datawave[ipost][0][p] 
     ydisptemp = datawave[ipost][1][p] 
     
     xdisp = mean(xdisptemp) - zerow[ipost][0][0]
     ydisp = mean(ydisptemp) - zerow[ipost][1][0]

     directw[ipost][0] = xdisp/sqrt(xdisp^2+ydisp^2)
     directw[ipost][1] = ydisp/sqrt(xdisp^2+ydisp^2)
     //added  to calculate the magnitude of force vector
     forcevecmag[ipost] = sqrt(xdisp^2+ydisp^2)
   endif
   
//   if(linmaskval_intercept_int[ipost] ==1)
   if(linmaskval_intercept_int[ipost] ==1)
      ncell+=1
   endif   
    
endfor     

//loadwave/P=topractice/N=datawave  
//duplicate/O  $fitres_string datawave
//duplicate/O  $bgshiftstring  bgshiftwave

variable i, j,DX, DY
variable dataR

variable icell=0




DataR=dimsize(datawave,2)


variable intervalnum, intervalnum_unwind                        //number of points in histogram bins

intervalnum_unwind =floor( DataR/laginterval )
intervalnum = DataR- laginterval
//intervalnum = 5000


//make/O/N=(DataR-laginterval) Histo
//make/O/N= fixlength  Histo

variable npavg = 1           //number of posts put in histosinglepost 
//variable lseg = 18000
//variable nseg = floor(datar/lseg)
variable nseg = 1 
variable k 

make/O/N = (intervalnum) Histo_x, Histo_Y, histo_t, histo_r, Histo_d,Histo_x_sq
make/O/N = (numb, npost,4) singlehisto
make/O/N = (numb)  bghisto = 0, cellhisto=0, edgehisto=0, centerhisto=0
make/o/n=(intervalnum_unwind*ncell) histosingle_x, histosingle_y
make/o/n=(intervalnum_unwind*npavg)  histosinglepost_x, histosinglepost_y
make/o/n=(ncell)   cagesize,Lmax_msd
make/o/n=(nseg*ncell) lmax_cell



ncell = 0
 
wave maskflag, msd

variable x4m, x2m, x4m_avg, x2m_avg

make/o/n=(floor(npost/npavg)) ngpara

j=0


for (ipost = 0 ; ipost <  npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)
//    if (linmaskval_intercept_int[ipost] == 1) 
   if(linmaskval_intercept_int[ipost] ==1)
//       For(i=0; i < intervalnum; i+=laginterval)


       for(i=0; i<intervalnum; i+=1)


            Histo_x[i]=datawave[ipost][0][i+laginterval] - datawave[ipost][0][i]
            Histo_y[i]=datawave[ipost][1][i+laginterval] - datawave[ipost][1][i]      
            Histo_d[i] = sqrt(Histo_x[i]^2+Histo_y[i]^2)
            Histo_x_sq[i] = (Histo_y[i]*125)^2

            
           Histo_r[i] = histo_x[i]*directw[ipost][0] + histo_y[i]*directw[ipost][1]
           histo_t[i] = histo_x[i]*directw[ipost][1] - histo_y[i] *directw[ipost][0]            
            

        Endfor
        

        for(k=0; k<nseg;k+=1)
//                wavestats/q/r=(lseg*k, lseg*k+lseg) histo_x
             wavestats/q histo_d
//               wavestats/q histo_y
            Lmax_cell[icell*nseg+k] =125* max(v_max, abs(v_min))
        endfor    
        Lmax_msd[icell] = sqrt(10^6*msd[floor((laginterval-100)/10)+99][ipost])

        
      
         
        for (i=0; i<intervalnum_unwind; i+=1)
           histosingle_x[i+intervalnum_unwind*icell] =datawave[ipost][0][i*laginterval+laginterval] - datawave[ipost][0][i*laginterval]
           histosingle_y[i+intervalnum_unwind*icell] =datawave[ipost][1][i*laginterval+laginterval] - datawave[ipost][1][i*laginterval]
           histosinglepost_x[i+j*intervalnum_unwind] = datawave[ipost][0][i*laginterval+laginterval] - datawave[ipost][0][i*laginterval]
           histosinglepost_y[i+(j)*intervalnum_unwind] =  datawave[ipost][1][i*laginterval+laginterval] -datawave[ipost][1][i*laginterval]                        
      endfor        
        

        if(j==(npavg-1))
        
         wavestats/q  histosinglepost_x
       
         x4m_avg=0
         x2m_avg=0

        cagesize[icell] = max(abs(v_min),abs(v_max))


        for(i=0; i<(intervalnum_unwind*npavg); i+=1)
//          for(i=0; i<intervalnum; i+=1)
               x4m = (histosinglepost_x[i]-V_avg)^4
               x2m = (histosinglepost_x[i]-V_avg)^2
               x4m_avg+= x4m
               x2m_avg+= x2m
        endfor

        x4m_avg = x4m_avg/(intervalnum_unwind*npavg)
        x2m_avg = x2m_avg/(intervalnum_unwind*npavg)


//        x4m_avg = x4m_avg/intervalnum
//        x2m_avg = x2m_avg/intervalnum

        ngpara[icell] = x4m_avg/(3*(x2m_avg^2)) - 1        
        
        j=-1
 
        icell +=1
                
        endif
        
        j+=1
        

      
      wavestats/Q   histo_r
      radsdev = V_sdev
      wavestats/Q  histo_t
      tansdev = V_sdev
      
      disdif[ipost] = radsdev - tansdev
      disrad[ipost] = radsdev
      distan[ipost] = tansdev
      
      
       Make/N=100/O Histox_Hist, Histoy_Hist, histor_HIst, Histot_hist;DelayUpdate
       histogram/B={bins, binw, numb} Histo_x,Histox_Hist
       histogram/B={bins, binw, numb} Histo_y,Histoy_Hist
       histogram/B={bins, binw, numb} Histo_r,Histor_Hist
       histogram/B={bins, binw, numb} Histo_t,Histot_Hist                    
       make/N=(numb)/O Histox
       
       histox=bins+x*binw
       
       singlehisto[][ipost][2] = histox_hist[p]
       singlehisto[][ipost][3] = histoy_hist[p]
       singlehisto[][ipost][0] = histor_hist[p]
       singlehisto[][ipost][1] = histot_hist[p]
       
endif       
        
       if (maskwave[posti][postj] == 0)
           bghisto = bghisto + histox_hist
           nbg += 1
       elseif (maskwave[posti][postj] == 1)    
           cellhisto = cellhisto + histox_hist
           centerhisto = centerhisto + histox_hist
           ncell+=1
           ncenter+=1
       elseif (maskwave[posti][postj] == 6)
           cellhisto = cellhisto + histox_hist
           edgehisto = edgehisto + histox_hist
           ncell+=1
           nedge+=1
       endif
       

endfor      
 
//bghisto = bghisto/(nbg*intervalnum)
//cellhisto = cellhisto/(ncell*intervalnum)
//centerhisto = centerhisto/(ncenter * intervalnum)
//edgehisto = edgehisto/(nedge * intervalnum)

bghisto = bghisto/nbg
cellhisto = cellhisto/ncell
centerhisto = centerhisto/ncenter 
edgehisto = edgehisto/nedge
      
      
redimension/n=(icell)   ngpara      


end





function disphisto(datawave, xwave, maskwave, astrw)
wave datawave, maskwave, astrw, xwave
variable nxpost, nypost


nxpost = astrw[3]
nypost = astrw[4]

variable npost = nxpost*nypost

variable ipost, posti , postj

variable rval = 0, gval = 0, bval = 0

string mptrace

dowindow/K  histographrad
dowindow/K  histographtan

for (ipost = 0 ; ipost <  npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)
       rval = 0
       gval = 0
       bval = 0
       sprintf mptrace "p_%d" ipost     
 
       if (maskwave[posti][postj] == 0)
           bval = 50000
       elseif (maskwave[posti][postj] == 1)    
           rval = 50000
       endif
       
      if (maskwave[posti][postj] != 3)
       if(ipost==0)
        display/N=histographrad datawave[][ipost][0]/TN=$mptrace  vs xwave
        modifygraph tick=2, mirror=1
        modifygraph rgb($mptrace) = (rval, gval, bval)
        display/N=histographtan datawave[][ipost][1]/TN=$mptrace  vs xwave
        modifygraph tick=2, mirror=1
        modifygraph rgb($mptrace) = (rval, gval, bval)       
      else
         appendtograph/W=histographrad datawave[][ipost][0]/TN=$mptrace  vs xwave
         modifygraph/W=histographrad rgb($mptrace) = (rval, gval, bval)
         appendtograph/W=histographtan datawave[][ipost][1]/TN=$mptrace  vs xwave
         modifygraph/W=histographtan rgb($mptrace) = (rval, gval, bval)
      endif
    endif
    
endfor

end


//This function is used to calculate the correlation between magnitude of force vector map and 
//difference of distribution btw rad and tangental direction


//function histocorr(datawave, forcevect, astrw, maskwave)
//wave datawave, forcevect, astrw, maskwave
//variable nxpost, nypost


//nxpost = astrw[3]
//nypost = astrw[4]

//variable  npost = nxpost*nypost

//variable  ipost, posti , postj

//make/O/(N=50)     temphistorad, temphistotan   //used to instore histogram distribution of a single post

//make/O/(N=npost)      forcevecmag

//for( ipost = 0 ; ipost<npost; ipost+=1 )
 //    forcevecmag[ipost] = sqrt( forcevect[ipost][0]^2 + forcevect[ipost][1]^2 )
     
//     temphistorad = datawave[p][ipost][2]
 //    temphistotan = datawave[p][ipost][3]
     
//     curvefit
     
     
     
     
     
function   drawhisto_cellref(astrw, datawave, bgshiftwave, maskwave, zerow, lagtime, bins)
wave astrw, datawave, bgshiftwave, maskwave
wave zerow   // stroe initial positions of posts
variable lagtime, bins

//SVAR  rootstring
//SVAR  bgshiftstring, fitres_string
//variable  lagtime=1
variable npost, ipost
variable numb, binw

variable nxpost, nypost

variable framerate = 100

variable laginterval = lagtime * framerate

variable fixlength = 5000

variable nof = dimsize(datawave, 2)

variable xdisp, ydisp
     
          
make/O/N=(nof) xdisptemp, ydisptemp

variable radsdev, tansdev
variable posti, postj

variable CX, CY

variable ncell=0





numb = 50
binw = 2*abs(bins)/numb

npost = dimsize(datawave, 0)
nxpost = astrw[3]
nypost = astrw[4]

make/O/N=(npost, 2)   directw


for(ipost=0; ipost<npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)
         if(maskwave[posti][postj]==1) 	
            CX+=datawave[ipost][0][0]
            CY+=datawave[ipost][1][0]
            ncell+=1   
         endif
endfor         

CX=CX/ncell
CY=CY/ncell             
            
          

for(ipost=0; ipost<npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)
       
   if (maskwave[posti][postj] !=3 )    
//     xdisp = datawave[ipost][0][0] - zerow[ipost][0][0]
//     ydisp = datawave[ipost][1][0] - zerow[ipost][1][0]
     
     xdisptemp = datawave[ipost][0][p]
     ydisptemp = datawave[ipost][1][p]
     
     xdisp = CX - mean(xdisptemp)
     ydisp = CY - mean(ydisptemp) 

     directw[ipost][0] = xdisp/sqrt(xdisp^2+ydisp^2)
     directw[ipost][1] = ydisp/sqrt(xdisp^2+ydisp^2)
     //added  to calculate the magnitude of force vector
//     forcevecmag[ipost] = sqrt(xdisp^2+ydisp^2)
   endif
   

   
endfor

   

variable i, j
variable dataR


//variable ncell=0, nbg=0, ncenter=0, nedge=0


DataR=dimsize(datawave,2)

variable intervalnum

//intervalnum =  =floor( DataR/laginterval )
intervalnum = DataR- laginterval
//intervalnum = 5000


//make/O/N=(DataR-laginterval) Histo
//make/O/N= fixlength  Histo

make/O/N = (intervalnum) histo_x, histo_y, histo_ct, histo_cr
make/O/N = (numb, npost,2) singlehistoc
//make/O/N = (numb)  bghisto = 0, cellhisto=0, edgehisto=0, centerhisto=0



 
wave maskflag


for (ipost = 0 ; ipost <  npost; ipost+=1)
       posti = floor( ipost/nypost)
       postj = mod(ipost, nypost)

       for(i=0; i<intervalnum; i+=1)

            Histo_x[i]=(datawave[ipost][0][i+laginterval]-bgshiftwave[i+laginterval][0]) - (datawave[ipost][0][i]-bgshiftwave[i][0])
            Histo_y[i]=(datawave[ipost][1][i+laginterval]-bgshiftwave[i+laginterval][1]) - (datawave[ipost][1][i]-bgshiftwave[i][1])
           Histo_cr[i] = histo_x[i]*directw[ipost][0] + histo_y[i]*directw[ipost][1]
           histo_ct[i] = histo_x[i]*directw[ipost][1] - histo_y[i] *directw[ipost][0]

       endfor

       Make/N=100/O histocr_HIst, Histoct_hist;DelayUpdate
       histogram/B={bins, binw, numb} Histo_cr,Histocr_Hist
       histogram/B={bins, binw, numb} Histo_ct,Histoct_Hist                    
       make/N=(numb)/O Histox
       
       histox=bins+x*binw

       singlehistoc[][ipost][0] = histocr_hist[p]
       singlehistoc[][ipost][1] = histoct_hist[p]
   
     
endfor               


end