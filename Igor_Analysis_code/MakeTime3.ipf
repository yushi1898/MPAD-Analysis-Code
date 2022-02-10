#pragma rtGlobals=1		// Use modern global access method.
function MakeTime(f1,f2)
//SVAR MaskString, FitRes_String, FitErr_String, RootString, ShiftString
variable f1,f2
variable V_zero, TotalCount, fstart, j



TotalCount = f2-f1+1
fstart = f1+1
make/O/N=(TotalCount) vTime

string s1
sprintf s1, "RHOD%02d.tif",f1
GetFileFolderInfo /Z/P=MovieFolder s1
//GetFileFolderInfo  s1
V_zero=V_modificationdate
vTime[0] = V_modificationdate - V_zero

for (j=fstart; j<=f2; j+=1)
	sprintf s1, "RHOD%02d.tif",j
	GetFileFolderInfo /Q/Z/P=MovieFolder s1
//	if (V_flag==0 && V_isFile)
	vTime[j] = (V_modificationdate - V_zero)/60
//	endif
endfor
//save/O/P=ToPractice vTime as "Time.ibw"
	
end


//window ForceTime() : GraphStyle
//ModifyGraph mode=4,marker=19
//Label left "Avg Potential Energy (aJ)"
//Label bottom "Time (min)"
end

