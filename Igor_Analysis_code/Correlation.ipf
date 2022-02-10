#pragma rtGlobals=1		// Use modern global access method.

//this correlation function reads in the energies of all posts for all frames, stored in corrwave
//during makeitquiver ... makeitquiver must be run first!

function correlation(corrwave,f1,f2)

wave corrwave
variable f1, f2
variable i, j, ifile, k
variable horizposts, vertposts, totalnumposts, numframes = f2-f1 //should change this to be read in!
wave astr
wave cellmatrix
variable kcell=1

HorizPosts=ASTr[3]
VertPosts=ASTR[4]
TotalNumPosts=HorizPosts*VertPosts


make/o/n=(totalnumposts,numframes) corrwave2=0
make/o/n=(horizposts,vertposts) corrmap
//make/o/n=(horizposts,vertposts,numframes) corrmapmovie

make/o/n=6 dEPost0 =0
make/o/n=6 dEPostnn=0//, dEPostnn2 =0,dEPostnn3 =0,dEPostnn4=0,dEPostnn5 =0,dEPostnn6 =0//have to compare the changes in energy, start with between f20 and 79 (entire time field is on) adjust this later to take input
//6 nearest neighbors

variable post0, nn1, nn2, nn3, nn4, nn5, nn6 //used to indicate the post number of central post and 6 nearest neighbors
for(k = f1;k<=f2;k+=1) 

for(i = 0;i<=(horizposts);i+=1)
	for(j=0;j<=(vertposts);j+=1)

		post0 = vertposts*j+i
		if(cellmatrix[post0]==kcell)
		
		if(mod(post0,2)!=0) //define nearest neighbors relative to central post if odd number //these might not be necessary since in sum, they encompass the same posts
			nn1= post0 -1
			nn2 = post0-vertposts
			nn3 = post0+1
			nn4 = post0+vertposts+1
			nn5= post0 + vertposts
			nn6 = post0+vertposts-1
		else
			nn1 = post0-vertposts-1
			nn2 = post0-vertposts
			nn3 = post0-vertposts+1
			nn4 = post0+1
			nn5 = post0+vertposts
			nn6 = post0-1
			
		endif
		
		
		
		dEpost0= corrwave[post0][0][k] - corrwave[post0][0][f1] //difference with file 0...could be changed //right now fills 6 lines with thissame number
		dEpostnn[0] = corrwave[nn1][0][k] -  corrwave[nn1][0][f1]
		dEpostnn[1] = corrwave[nn2][0][k] -  corrwave[nn2][0][f1]
		dEpostnn[2] = corrwave[nn3][0][k] -  corrwave[nn3][0][f1]
		dEpostnn[3] = corrwave[nn4][0][k] -  corrwave[nn4][0][f1]
		dEpostnn[4] = corrwave[nn5][0][k] -  corrwave[nn5][0][f1]
		dEpostnn[5] = corrwave[nn6][0][k] -  corrwave[nn6][0][f1]

//		WaveStats/Q dEpost1 //used for normalizing to 1
//		Variable srcRMS= V_rms
//		Variable srcLen= numpnts(dEpost1)	
//
//		WaveStats/Q dEpost2
//		Variable destRMS= V_rms
//		variable destLen= numpnts(dEpost2)	
//			
		correlate dEpost0, dEpostnn
//			
//		dEpost2 /= (srcRMS * sqrt(srcLen) * destRMS * sqrt(destLen))
//	
		corrwave2 [post0][k] = dEpostnn[0]
		//corrwave2 [post0][k]/=sqrt(corrwave[vertposts*j+i][0][k+1]*corrwave[vertposts*j+i][0][k+1]) //attempting to divide by the average energy of the post in both frames, to scale properly
		//corrmapmovie[i][j][k] = dEpost2
		endif
	endfor
endfor

endfor

makelinearmaskwaves() //just want to make sure we have "perfect array" this goes back to the original frame, tc.
NewMovie/O/P=ToPractice/F=(1)/L as "CorrMov.mov" 
make/o/n=(totalnumposts) colormap =0
for(k = f1;k<f2;k+=1)

	NewImage/K=1/N=RotatedImage RotatedOrig1
	ModifyImage RotatedOrig1 ctab= {AsTR[6],AStr[7],Grays,0}
	      		
	AppendToGraph/T  makeshiftY vs makeshiftX
	 
	colormap= corrwave2[p][k]  		
	ModifyGraph zcolor = {colormap,*,*,rainbow,1}
	ModifyGraph mode=3,marker=55,msize=10.5
	
//	AppendToGraph/T  blueshiftY vs blueshiftX
	
	//ModifyGraph zcolor = {corrwave2,*,*,rainbow,1}
	//ModifyGraph mode=3,marker=19,msize=3     		
	  
	
	AddMovieFrame
	
	dowindow rotatedimage
 
endfor
closemovie	

end