#pragma rtGlobals=1		// Use modern global access method.

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

function ViewErrors(f1,f2,posti,postf,param1,FEwave)
    variable f1  // first frame
    variable f2  // last frame
    variable posti  // first post
    variable postf   // last post
    variable param1   // index of fitparam to be viewed
    wave FEwave    // fit error wave

    variable i,j
    
    make/O/n=(postf-posti+1,f2-f1+1) errorimage // always created. calling routine must copy it as needed
    SetScale/I x posti,postf,"", errorimage                   
    SetScale/I y f1,f2,"", errorimage                   
    errorimage = 0
    for (i = posti; i <= postf; i+=1)
        for (j = f1; j <= f2; j+=1)
            errorimage[i-posti][j-f1] = FEwave[i][param1][j]
        endfor
    endfor
    
    DoWindow/K FitErrImage
    NewImage/N=FitErrImage errorimage
    
    WaveStats errorimage
    Make/O/N=100/O W_Hist;DelayUpdate
    Histogram/B={V_min,(V_max-V_min)/99,100} errorimage,W_Hist
    DoWindow/K ErrorHist
    Display/N=ErrorHist w_hist
    ModifyGraph log(left)=1
    ModifyGraph mode=4


end

    
