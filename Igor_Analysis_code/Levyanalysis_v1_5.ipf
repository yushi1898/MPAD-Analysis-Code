#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//latest version of large posts , 10 fps videos

//incoportate the function to calculate averaged cross correlation function
//integrated functions used in MSD analysis and Levy analysis

function  Glevyanlysis()



Dowindow/K  Gfit_levyanlysis
NewPanel/K=1/N=Gfit_levyanlysis/W=(100,100,385, 500)
SetDrawEnv linefgc= (0,0, 53520),fillpat= 0, linethick=2;	DrawRect 5,5,265,65
SetDrawEnv linefgc= (0,53520, 0),fillpat= 0, linethick=2;	DrawRect 5,75,265,160
SetDrawEnv linefgc= (53520,53520, 0),fillpat= 0, linethick=2;	DrawRect 5,170,265,260
SetDrawEnv linefgc= (64512,14848,14848),arrow= 1,linethick= 4.00,arrowfat= 1.00; DrawLine 275,10,275,370
SetDrawLayer UserBack

Button  CELLMASK, pos={15,10},size={110,20}, proc = gencellmask, title="Make Cell Mask"
Button  FORCEBIFURCATION, pos = {135,10}, size = {110,20}, proc = bifcell, title = "Force Bicurcation"
Button  CELLMASK_seg, pos={10,35},size={130,20}, proc = gencellmask_seg, title="Make Cell Mask Seg"
Button  FORCEBIFURCATION_seg, pos = {135,35}, size = {130,20}, proc = bifcell_seg, title = "Force Bicurcation Seg"
Button  Vanhoffcurve, pos={10,80},size={250,20}, proc=genvanhoff, title = "Generate Vanhove"
Button Lmaxdist, pos={10,105}, size = {250,20}, proc=genlmaxdist, title = "Genertae lmax"
Button anisotropy, pos={10,130}, size={250,20}, proc=genanisotropy, title = "Calculate anisotropy"
Button  twopmicrorheology, pos={10,180}, size={250,20}, proc=calc2pointMR, title = "2 point microrhoelogy"
Button disp2pmr, pos={10,210}, size= {250,20}, proc= dispind2pmr, title = "show individual 2P MR trace"
Button crosscorravg, pos={10,240}, size={250,20}, proc = calcCCcontrol, title = "calculate averaged crosscorrelation"


end


function gencellmask(ctrlName): ButtonControl
string ctrlName

wave fitres_ed
wave linmaskval
//input
variable framenum = 6000
variable slopethresh =0.5
variable bg_thresh = 0.2
variable recalc_msd_flag = 0

//famenum   : number of frames to calculate threshold
Prompt framenum, "Enter number of frames to calculate MSD slope"
Prompt slopethresh, "Enter threshold for slope"
Prompt bg_thresh,"Enter threshold for bg posts"
Prompt recalc_msd_flag, "type 0 for recalculate MSD"

DoPrompt "parameters" framenum, slopethresh,bg_thresh, recalc_msd_flag
cellmotdetect(fitres_ed, framenum, slopethresh,bg_thresh,recalc_msd_flag)
//cellmotdetect_seg(fitres_ed, framenum, slopethresh)

end



function gencellmask_seg(ctrlName): ButtonControl
string ctrlName

wave fitres_ed
wave linmaskval
//input
variable framenum = 900
variable slopethresh =0.5

//famenum   : number of frames to calculate threshold
Prompt framenum, "Enter number of frames to calculate MSD slope"
Prompt slopethresh, "Enter threshold for slope"

DoPrompt "parameters" framenum, slopethresh
//cellmotdetect(fitres_ed_s, framenum, slopethresh)
cellmotdetect_seg(fitres_ed, framenum, slopethresh)

end

function bifcell(ctrlName):ButtonControl
string ctrlName
wave linmaskval_intercept
variable posttype
prompt posttype, "enter micropost type"
doprompt "parameter", posttype
wave fitres_ed


//findinteriorposts_seg(fitres_ed, linmaskval_seg,posttype)
findinteriorposts_dis(fitres_ed, linmaskval_intercept,posttype)


End

function bifcell_seg(ctrlName):ButtonControl
string ctrlName
wave linmaskval_seg
variable posttype
prompt posttype, "enter micropost type"
doprompt "parameter", posttype
wave fitres_ed


findinteriorposts_seg(fitres_ed, linmaskval_seg,posttype)
//findinteriorposts_dis(fitres_ed, linmaskval_intercept,posttype)


End




function genvanhoff(ctrlName): Buttoncontrol
string ctrlName
newpath/c/o levywavefolder 

wave zeropos, astr, fitres_ed_s

SVAR FitRes_String, AttributeString,MaskString, bgshiftstring
variable startbin = 0.2

make/o/n=16 lagtime = {0.1,0.2,0.3,0.5, 1,2,3,5,10,20,30,50,100,200,300,500}

variable nlagtime = dimsize(lagtime,0)
variable posttype =1 

prompt posttype, "type of posts"
doprompt "parameter", posttype

string lagtimeaffix, singlehistostring_x, singlehistostring_y, singlehistostring_r, singlehistostring_t

variable i

wave histosingle_y, histosingle_x, histosingle_r, histosingle_t

for(i=0; i<nlagtime;  i+=1)
    if(lagtime[i] <1)      
        sprintf lagtimeaffix "_0%ds", lagtime[i]*10
      else
         sprintf lagtimeaffix "_%ds", lagtime[i]
      endif
singlehistostring_x = "histosingle_x"+lagtimeaffix+".txt"
singlehistostring_y = "histosingle_y"+lagtimeaffix+".txt"
singlehistostring_r = "histosingle_r"+lagtimeaffix+".txt"
singlehistostring_t = "histosingle_t"+lagtimeaffix+".txt"

//drawhisto(astr, $fitres_string, $bgshiftstring, $maskstring, zeropos, lagtime[i], -startbin,posttype)
drawhisto_shortV(astr, fitres_ed_s,$maskstring, zeropos, lagtime[i], -startbin,posttype)
if(posttype ==1 || posttype ==0)
save/G/o/A=2/p=levywavefolder histosingle_y as singlehistostring_y
save/G/o/A=2/p=levywavefolder histosingle_x as singlehistostring_x
elseif (posttype ==6)
save/G/o/A=2/p=levywavefolder histosingle_r as singlehistostring_r
save/G/o/A=2/p=levywavefolder histosingle_t as singlehistostring_t
endif


endfor


end


function genlmaxdist(ctrlName): Buttoncontrol
string ctrlName
newpath/c/o lmaxfolder
wave astr, zeropos
variable startbin = 0.2


SVAR FitRes_String, AttributeString,MaskString, bgshiftstring

//variable lagtime = 100
make/o/n=(2) lagtime = {20,100}
variable posttype

wave lmax_msd, lmax_cell

string lmaxstr, lstarstr, expstr, ngstr
variable i

for(i=0;i<dimsize(lagtime,0);i+=1)

//generage lmax and lstar for low traction post   //posttype = 1
lmaxstr = "Lmax_lb.txt"
lstarstr = "Lmsd_lb.txt"
expstr = "alpha_lb.txt"
sprintf ngstr  "ngpara_lb_%d.txt" ,lagtime[i]
posttype = 1

drawhisto(astr, $fitres_string, $bgshiftstring, $maskstring, zeropos, lagtime[i], -startbin,posttype)
wave alpha_msd, ngpara

//save/G/o/A=2/p=lmaxfolder lmax_cell as lmaxstr
save/G/o/A=2/p=lmaxfolder lmax_msd as lstarstr
//save/G/o/A=2/p=lmaxfolder alpha_msd as expstr
//save/G/o/A=2/p=lmaxfolder ngpara as ngstr

//generate lmax and lstar for high traction post     posttype = 6
lmaxstr = "Lmax_hb.txt"
lstarstr = "Lmsd_hb.txt"
expstr = "alpha_hb.txt"
sprintf ngstr  "ngpara_hb_%d.txt" ,lagtime[i]
posttype = 6
drawhisto(astr, $fitres_string, $bgshiftstring, $maskstring, zeropos, lagtime[i], -startbin,posttype)
//save/G/o/A=2/p=lmaxfolder lmax_cell as lmaxstr
save/G/o/A=2/p=lmaxfolder lmax_msd as lstarstr
//save/G/o/A=2/p=lmaxfolder alpha_msd as expstr
//save/G/o/A=2/p=lmaxfolder ngpara as ngstr

endfor

end

function genanisotropy(ctrlName): Buttoncontrol
string ctrlName
wave astr, fitres_ed_sub
newpath/c/o lmaxfolder
variable seglength = 900
variable posttype 
string anisotropystr,  angledifstr, anisotropynum_str
 
 

posttype = 1

findanisotropy_seg(fitres_ed_sub, astr, seglength, posttype)
findanisotropy_seg(fitres_ed_sub, astr, seglength, posttype)
wave cposanisotropy_raw, angledif, cposanisotropy_num

anisotropystr = "Anisotropy_lb.txt"
angledifstr = "angledif_lb.txt"
anisotropynum_str = "anisonum_lb.txt"
save/G/o/A=2/p=lmaxfolder cposanisotropy_raw as anisotropystr
save/G/o/A=2/p=lmaxfolder angledif as angledifstr
save/G/o/A=2/p=lmaxfolder cposanisotropy_num as anisotropynum_str


posttype = 6

findanisotropy_seg(fitres_ed_sub, astr, seglength, posttype)
findanisotropy_seg(fitres_ed_sub, astr, seglength, posttype)

anisotropystr = "Anisotropy_hb.txt"
angledifstr = "angledif_hb.txt"
anisotropynum_str = "anisonum_hb.txt"
save/G/o/A=2/p=lmaxfolder cposanisotropy_raw as anisotropystr
save/G/o/A=2/p=lmaxfolder angledif as angledifstr
save/G/o/A=2/p=lmaxfolder cposanisotropy_num as anisotropynum_str



end



//binsize in space fixed for this version
function calc2pointMR(ctrlName): Buttoncontrol
string ctrlName

//wave datawave, maskwave, origpos, astr
wave fitres_ed, linmaskval_intercept_int, zeropos,  astr
wave fitres_edrd
//datawave = fitres_ed
//maskwave = linmaskval_intercept_int
//origpos = zeropos

variable nbinslag = 200
variable framegap = 0.1
variable posttype = 1 
variable pixelratio = 0.125

variable particlerange = 12
variable NNP_r = 4

prompt  particlerange, "cutoff distance between particles"
prompt framegap, "frame rate"
prompt posttype, "post type"
doprompt "parameter", particlerange, framegap, posttype

variable maxgapx = ceil(particlerange/NNP_r)
variable maxgapy = 2*maxgapx




calcTPMR(astr,fitres_edrd, linmaskval_intercept_int,zeropos, particlerange, NNP_r, maxgapx, maxgapy, nbinslag, framegap,posttype,pixelratio)



end


function dispind2pmr(ctrlName):ButtonControl
string ctrlName
wave linmaskval_intercept_int, twoPM_corr, twopm_lagtime
variable posttype = 1
prompt posttype, "enter micropost type"
doprompt "parameter", posttype



displayind2pMcorr(twoPM_corr, linmaskval_intercept_int, posttype,twopm_lagtime)    


End



function calcCCcontrol(ctrlName): Buttoncontrol
string ctrlName

//wave datawave, maskwave, origpos, astr
wave fitres_ed, linmaskval_intercept_int, zeropos,  astr
//wave fitres_edrd
//datawave = fitres_ed
//maskwave = linmaskval_intercept_int
//origpos = zeropos

//variable nbinslag = 200
//variable framegap = 0.1
variable posttype = 1 
variable pixelratio = 0.125

variable particlerange = 50
variable NNP_r = 4

prompt  particlerange, "cutoff distance between particles"
//prompt framegap, "frame rate"
prompt posttype, "post type"
doprompt "parameter", particlerange, posttype

variable maxgapx = ceil(particlerange/NNP_r)
variable maxgapy = 2*maxgapx




calcCCfunc(astr,fitres_ed, linmaskval_intercept_int,zeropos, particlerange, NNP_r, maxgapx, maxgapy,posttype,pixelratio)



end



function calcTPMR(astr,datawave,maskwave, origpos, particlerange, NNP_r, maxgapx, maxgapy, nbinslag,framegap, posttype,pixelratio)
wave datawave, maskwave, origpos,astr
variable particlerange, maxgapx, maxgapy, nbinslag,framegap, NNP_r, posttype,pixelratio

variable npost = dimsize(datawave,0)
variable videolength = dimsize(datawave,2)
variable i,j, k
i=0
j=1

//generating log-space lag time wave
variable lagtimecutoff = videolength/3


make/o/n=(nbinslag) twoPM_lagtime_raw, twoPM_lagtime
variable lagtimelogmin, lagtimelogmax, lagtimebingap
lagtimelogmin = log(framegap)
lagtimelogmax = log(lagtimecutoff*framegap)
lagtimebingap = (lagtimelogmax - lagtimelogmin)/(nbinslag-1)

for(i=0;i<nbinslag;i+=1)
     twoPM_lagtime_raw[i] = floor(10^(lagtimelogmin+i*lagtimebingap)/framegap)
endfor

twoPM_lagtime[0] = twoPM_lagtime_raw[0]
for(i=1;i<nbinslag;i+=1)
     if(twoPM_lagtime_raw[i] != twoPM_lagtime_raw[i-1])
          twoPM_lagtime[j] = twoPM_lagtime_raw[i]
         j+=1            
     endif
endfor     


nbinslag = j+1
redimension/N=(j) twoPM_lagtime

        


//generating distance bin wave
variable nbinsdis = (particlerange-NNP_r)+2
make/o/n=(nbinsdis) twoPM_dis, dis_count, dis_count_all
twoPM_dis = x+NNP_r-0.5

variable posti, postj, idx, idy
variable dr, costheta, sintheta
variable dt, disindex
variable dx, dy
variable ixmin, ixmax, iymin, iymax
variable npostx= astr[3]
variable nposty = astr[4]

wave drr, dtt

//data strucutre of two point microrheology wave : row ---> post number, column--->lag time, layer--->distance, 
//chunk--->0 for longtitude 1 for tangential
make/o/n=(npost, nbinslag, nbinsdis, 2)  TwoPM_Corr=0 
make/o/n=(nbinslag, nbinsdis, 2) TwoPM_Corr_ap = 0


dis_count_all = 0
//changed for testing
for(posti = 0; posti<npost; posti+=1)
      if(maskwave[posti] == posttype)
      idx = floor(posti/nposty)
      idy = mod(posti, nposty)
      ixmin = max(0,idx-maxgapx)
      ixmax = min(idx+maxgapx, npostx-1)
      iymin = max(0,idy-maxgapy)
      iymax = min(idy+maxgapy, nposty-1)
      dis_count = 0
      for(i=ixmin; i<ixmax; i+=1)
          for(j=iymin; j<iymax; j+=1) 
                postj = i*nposty+j
                if ((postj!=posti) && (maskwave[postj] == posttype)) 
                    dx = (origpos[posti][0][0] - origpos[postj][0][0])*pixelratio
                    dy = (origpos[posti][1][0] - origpos[postj][1][0])*pixelratio
                    dr = sqrt(dx^2+dy^2)
                    costheta = dx/dr
                    sintheta = dy/dr
                    if(dr < (particlerange+0.5) )
                              disindex = round(dr) - 4
//                              dt =  twoPM_lagtime[k]
                              make/o/n=(videolength) dwtempx1, dwtempy1, dwtempx2, dwtempy2
                              dwtempx1 = datawave[posti][0][p]*pixelratio
                              dwtempy1 = datawave[posti][1][p]*pixelratio
                              dwtempx2 =  datawave[postj][0][p]*pixelratio
                              dwtempy2 =  datawave[postj][1][p]*pixelratio
                              calcindvid2pmr(dwtempx1, dwtempx2, dwtempy1, dwtempy2,twoPM_lagtime, costheta, sintheta)
                              TwoPM_corr[posti][][disindex][0] += drr[q]
                              TwoPM_corr[posti][][disindex][1] += dtt[q]
                              TwoPM_Corr_ap[][disindex][0] += drr[p]
                              TwoPM_Corr_ap[][disindex][1] += dtt[p]
                              dis_count[disindex] +=1
                              dis_count_all[disindex] +=1
                    endif
               
                endif
             
             endfor
       
       endfor         
    
    
    if(dis_count[0] !=0)      
       TwoPM_corr[posti][][][0] /= dis_count[r]
       TwoPM_corr[posti][][][1] /= dis_count[r]
    endif
     

      
      endif
      
      
endfor      
      
     TwoPM_corr_ap[][][0] /= dis_count_all[q]
     TwoPM_corr_ap[][][1] /= dis_count_all[q]      
      
      
end


function calcindvid2pmr(dwtempx1, dwtempx2, dwtempy1, dwtempy2, lagtimew, costheta, sintheta)
wave   dwtempx1, dwtempx2, dwtempy1, dwtempy2, lagtimew
variable costheta, sintheta

variable nbinslag = dimsize(lagtimew,0)
variable len = dimsize(dwtempx1,0)

make/o/n=(nbinslag) drr, dtt

variable i
variable dt


for(i=0;i<nbinslag;i+=1)
     dt = lagtimew[i]
     duplicate/o/r = [0, len-dt-1] dwtempx1, w2ptmrx1
     duplicate/o/r = [0,len-dt-1] dwtempy1, w2ptmry1
     duplicate/o/r = [dt, len-1] dwtempx1,w2ptmrx1_s
     duplicate/o/r = [dt, len-1] dwtempy1, w2ptmry1_s   
     duplicate/o/r = [0, len-dt-1] dwtempx2, w2ptmrx2
     duplicate/o/r = [0,len-dt-1] dwtempy2, w2ptmry2
     duplicate/o/r = [dt, len-1] dwtempx2,w2ptmrx2_s
     duplicate/o/r = [dt, len-1] dwtempy2, w2ptmry2_s  
     make/o/n=(len-dt) dxtemp1, dytemp1, dxtemp2, dytemp2, drrtemp, dtttemp      
     dxtemp1 = w2ptmrx1 - w2ptmrx1_s
     dxtemp2 = w2ptmrx2 - w2ptmrx2_s
     dytemp1 = w2ptmry1 - w2ptmry1_s
     dytemp2 = w2ptmry2 - w2ptmry2_s
     drrtemp = (dxtemp1*costheta + dytemp1*sintheta) * ( dxtemp2*costheta+dytemp2*sintheta)
     dtttemp = (dytemp1*costheta - dxtemp1*sintheta) * (dytemp2*costheta - dxtemp2*sintheta)
     drr[i] = mean(drrtemp)
     dtt[i] = mean(dtttemp)
endfor


end

    
    
function displayind2pMcorr(corrwave,maskwave, posttype,lagtimew)    
wave corrwave, maskwave, lagtimew
variable posttype


variable npost = dimsize(corrwave,0)
variable count = 0
variable i
make/o/n=(npost) twopm_mag_rr, twopm_mag_tt

for(i=0; i<npost;i+=1)
     if(maskwave[i] == posttype)
         twopm_mag_rr[i] = corrwave[i][108][0][0]
         twopm_mag_tt[i] = corrwave[i][108][0][1]
         if(count==0)
           display corrwave[i][][0][1] vs lagtimew
        else
           appendtograph corrwave[i][][0][1] vs lagtimew
        endif
        count+=1
     endif
endfor

end







function genrdwk(stepnum)
variable stepnum
variable i

make/o/n=(stepnum)  rdwk_x=0, rdwk_y=0

variable steprange = 1
//variable stepthresh = steprange/2

for(i=1;i<stepnum; i+=1)
      if(enoise(steprange)< 0 )
          rdwk_x[i] = rdwk_x[i-1] + 1
      else 
          rdwk_x[i] = rdwk_x[i-1] - 1
      endif
      
      if(enoise(steprange)< 0 )
         rdwk_y[i] = rdwk_y[i-1] + 1
      else 
         rdwk_y[i] =rdwk_y[i-1] - 1
      endif      

endfor



end  

function rdwkmatrix(postnum,stepnum)
variable postnum,stepnum

wave rdwk_x, rdwk_y

variable i

make/o/n=(postnum,2,stepnum) rdwkmatri
for(i=0;i<postnum;i+=1)
     genrdwk(stepnum)
     rdwkmatri[i][0][] = rdwk_x[r]
     rdwkmatri[i][1][] = rdwk_y[r]
endfor

end