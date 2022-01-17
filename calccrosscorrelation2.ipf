#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// 9/7/21 DHR
// calcCCfunc() edited by DHR to make compile in Igor 9
//   Correct operation not guaranteed!
function calcspearmancorr(dwave1, dwave2)
wave dwave1, dwave2
duplicate/o dwave1 dwave1O
duplicate/o dwave2 dwave2O
dwave1O = x
dwave2O = x
sort dwave1 dwave1O
sort dwave2 dwave2O

variable i

return statscorrelation(dwave1O, dwave2O)


end




//this function used to calculate averaging cross correaltion function
function calcCCfunc(astr,datawave,maskwave, origpos, particlerange, NNP_r, maxgapx, maxgapy, posttype,pixelratio)
wave datawave, maskwave, origpos,astr
variable particlerange, maxgapx, maxgapy, NNP_r, posttype,pixelratio

print "Edits by DHR 9/7/21 at end to make compile in Igor9.  Correct operation not guaranteed!"

variable npost = dimsize(datawave,0)
variable videolength = dimsize(datawave,2)
variable i,j, k
variable p_thresh = 0.5
i=0
j=1

        
//generating distance bin wave based on hex grid 10/9/18
variable nbinsdis = floor(particlerange/NNP_r)+2
make/o/n=(nbinsdis) twoPM_dis, dis_count, dis_count_all,dis_count_all_r_plus, dis_count_all_r_minus,dis_count_all_t_plus, dis_count_all_t_minus
//twoPM_dis = x+NNP_r-0.5

twoPM_dis[0] = NNP_r
twoPM_dis[1] = NNP_r*sqrt(3)
twoPM_dis[2] = NNP_r*2
twoPM_dis[3] = NNP_r*sqrt(7)
twoPM_dis[4] = NNP_r*3

variable dis_cutoff =4

for(i=dis_cutoff;i<nbinsdis;i+=1)
     twoPM_dis[i] = NNP_r*(i-1) + NNP_r*0.5
     endfor

variable posti, postj, idx, idy
variable dr, costheta, sintheta
variable dt, disindex
variable dx, dy
variable ixmin, ixmax, iymin, iymax
variable npostx= astr[3]
variable nposty = astr[4]
variable cor_coef_r, cor_coef_t 

wave drr, dtt

//data strucutre of two point microrheology wave : row ---> post number, column--->lag time, layer--->distance, 
//chunk--->0 for longtitude 1 for tangential
make/o/n=(npost, nbinsdis, 2)  CC_Corr=0 
make/o/n=(nbinsdis, 2) CC_Corr_ap = 0, CC_cor_ap_plus=0, CC_cor_ap_minus=0
make/o/n=(videolength) dwtempx1, dwtempy1, dwtempx2, dwtempy2, dwtempr1, dwtempr2, dwtempt1, dwtempt2

dis_count_all = 0
dis_count_all_r_plus =0
dis_count_all_t_plus =0
dis_count_all_r_minus = 0
dis_count_all_t_minus = 0
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
                    if(dr < (particlerange+p_thresh) )
                        if (dr<twoPM_dis[0])
                             disindex = 0
                         elseif(dr>twoPM_dis[nbinsdis-1]+p_thresh)
                             disindex = nbinsdis - 1
                         else    
                       for(k=0;k<(nbinsdis-1);k+=1)
                           if ((k<dis_cutoff) && (dr>twoPM_dis[k]-p_thresh) && (dr<twoPM_dis[k]+p_thresh) )
                              disindex = k
                              break
                           elseif ( (k>=dis_cutoff) && (dr>twoPM_dis[k]-NNP_r/2+p_thresh) && (dr< twoPM_dis[k] +NNP_r/2 - p_thresh ))
                              disindex = k
                              break
                          endif
                       endfor 
                       endif                         
//                             disindex = round(dr) - 4
//                              dt =  twoPM_lagtime[k]

                              dwtempx1 = datawave[posti][0][p]*pixelratio
                              dwtempy1 = datawave[posti][1][p]*pixelratio
                              dwtempx2 =  datawave[postj][0][p]*pixelratio
                              dwtempy2 =  datawave[postj][1][p]*pixelratio
                              dwtempr1 =  dwtempx1* costheta + dwtempy1*sintheta
                              dwtempt1 =  dwtempy1*costheta - dwtempx1*sintheta
                              dwtempr2 =   dwtempx2* costheta + dwtempy2*sintheta
                              dwtempt2 =  dwtempy2*costheta - dwtempx2*sintheta
                              cor_coef_r = calcspearmancorr(dwtempr1, dwtempr2)
                              cor_coef_t = calcspearmancorr(dwtempt1, dwtempt2)
                              
                              CC_corr[posti][disindex][0] += cor_coef_r
                              CC_corr[posti][disindex][1] += cor_coef_t
                              CC_corr_ap[disindex][0] += cor_coef_r
                              CC_corr_ap[disindex][1] += cor_coef_t
                              
                              if(cor_coef_r >=0)
                                 CC_cor_ap_plus[disindex][0]+=cor_coef_r
                                 dis_count_all_r_plus[disindex]+=1
                              else
                                 CC_cor_ap_minus[disindex][0] += cor_coef_r
                                 dis_count_all_r_minus[disindex]+=1
                              endif
                              
                              if(cor_coef_t>=0)
                                 CC_cor_ap_plus[disindex][1] +=cor_coef_t
                                 dis_count_all_t_plus[disindex]+=1
                              else
                                 CC_cor_ap_minus[disindex][1] += cor_coef_t   
                                 dis_count_all_t_minus[disindex]+=1
                              endif
                              
                              dis_count[disindex] +=1
                              dis_count_all[disindex] +=1
                    endif
               
                endif
             
             endfor
       
       endfor         
    
    
    if(dis_count[0] !=0)      
       CC_corr[posti][][0] /= dis_count[q]
       CC_corr[posti][][1] /= dis_count[q]
    endif
     

      
      endif
      
      
endfor      
      //average all posts
     CC_corr_ap[][0] /= dis_count_all[p]
     CC_corr_ap[][1] /= dis_count_all[p]      

//9/7/21  Edited by DHR to make compile     
     CC_cor_ap_plus[][0]/= dis_count_all_r_plus[p]
     CC_cor_ap_plus[][1]/=dis_count_all_t_plus[p]
     CC_cor_ap_minus[][0]/=dis_count_all_r_minus[p]
     CC_cor_ap_minus[][1]/=dis_count_all_t_minus[p]
// OLD CODE replaced with above - was missing [] for 2d waves on LHS
//     CC_cor_ap_plus[0]/= dis_count_all_r_plus[p]
//     CC_cor_ap_plus[1]/=dis_count_all_t_plus[p]
//     CC_cor_ap_minus[0]/=dis_count_all_r_minus[p]
//     CC_cor_ap_minus[1]/=dis_count_all_t_minus[p]
      
      
end



