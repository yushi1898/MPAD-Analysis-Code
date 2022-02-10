#include "stdafx.h"
#include "DemoModuleModule.h"
#include "NIDAQmx.h"
#include "stdio.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define PI	3.1415926535
#define nop 10000

//edited on 3/9/2015  differet synchronized method

//Global Variable
FILE    *fp;
////FILE    *frp;
float64  data[nop*10]={0.0};
//float64  data2[10]={0.0};
float64  dataread[42000];
//float64    dataread[2];
int32    num=0;
int32    read;
int32    iframe=0;
int32    nof, now;
TaskHandle  taskHandle=0;
TaskHandle  taskHandleread=0;
char   filename[50];




//int32    stopindex=-1;
//int32    freq=1;

CDemoModuleModule::CDemoModuleModule()
{
    Name = _T("Demo Module");
    DLLName = _T("DemoModule.dll");

    DemoModuleDlg = new CDemoModuleDlg(NULL, this);
    CommandIncrement = NULL;
	CommandStart=NULL;
	CommandStop=NULL;
    EventOnButtonClick = NULL;
}
//-----------------------------------------------------------------------------
CDemoModuleModule::~CDemoModuleModule()
{
    //Remove the modules dialog
    GUIRequest->RemoveDialog(DemoModuleDlg);
    DemoModuleDlg->DestroyWindow();
    delete DemoModuleDlg;

    CommandInfo.clear();

    if(CommandIncrement)
        delete CommandIncrement;

	if(CommandStart)
		delete CommandStart;

	if(CommandStop)
		delete CommandStop;



    Events.clear(); //Clears the event list

    if(EventOnButtonClick)
        delete EventOnButtonClick;
}
//-----------------------------------------------------------------------------
void CDemoModuleModule::CreateCommands()
{
    //COMMAND_INCREMENT
    CommandIncrement = new CScriptCommandInfo();
    CommandIncrement->SetName(CString(_T("Increment")));
    CommandIncrement->SetDescription(CString(_T("Increment the number of frames captured and display that number.")));
    CommandIncrement->SetLine(CString(_T("Increment frame count.")));
    CommandIncrement->SetModule(this);
    CommandIncrement->SetCommandID(COMMAND_INCREMENT);
    CommandInfo.push_back(CommandIncrement);

	//COMMAND_STARTER
	CommandStart = new CScriptCommandInfo();
	CommandStart->SetName(CString(_T("Start")));
	CommandStart->SetDescription(CString(_T("reset the number of frames when start recording.")));
	CommandStart->SetLine(CString(_T("reset frame cound.")));
	CommandStart->SetModule(this);
	CommandStart->SetCommandID(COMMAND_START);
	CommandInfo.push_back(CommandStart);

	//COMMAND_STOP
	CommandStop = new CScriptCommandInfo();
	CommandStop->SetName(CString(_T("Stop")));
	CommandStop->SetDescription(CString(_T("reset the number of frames when Stop recording.")));
	CommandStop->SetLine(CString(_T("reset frame cound.")));
	CommandStop->SetModule(this);
	CommandStop->SetCommandID(COMMAND_STOP);
	CommandInfo.push_back(CommandStop);

}
//-----------------------------------------------------------------------------
void CDemoModuleModule::CreateVariables()
{
    IsComputerOnFire.Initialize(_T("IsComputerOnFire"),
        _T("True if the computer is on fire."));

    Variables.push_back(&IsComputerOnFire);

    IsComputerOnFire = false;
}
//-----------------------------------------------------------------------------
void CDemoModuleModule::CreateEvents()
{
    EventOnButtonClick = new CScriptEvent(_T("On Button Click"), 
        _T("This event is called when the user click on the button."), 
        EVENT_ON_BUTTON_CLICK, this);
    
    //Declare local parameters
    EventOnButtonClick->AddEventParameter(CParamInfo(_T("CVarUINT"), _T("frameCount"), _T("The frame count.")));

    RegisterEvent(EventOnButtonClick);
}
//-----------------------------------------------------------------------------
void CDemoModuleModule::Initialize()
{
    DemoModuleDlg->SetParentSetup(GetParentSetup());
    GUIRequest->AskForDialog(DemoModuleDlg);

    CModule::Initialize();
}
//-----------------------------------------------------------------------------
void CDemoModuleModule::PreEventScript(int eventID, PARAMETERS& localParams)
{

}
//-----------------------------------------------------------------------------
void CDemoModuleModule::PostEventScript(int eventID, PARAMETERS& localParams)
{

}
//-----------------------------------------------------------------------------
bool CDemoModuleModule::ExecuteCommand(int commandID, PARAMETERS& params, PARAMETERS& localParams)
{
	//subject to 2 frame shift   8/24/15

    int32       error=0;

	int32       i=0, j=0;

	float64     freq[17]={0.1,0.2,0.5,0.8,1.0,2.0,4.0,5.0,8.0,10.0,20.0,35.0,55.0,80.0,95.0,115.0,135.0};

	int32       datalength[17]={18000,18000,18000,18000,6000,6000,6000,6000,6000,6000,3000,3000,3000,3000,3000,3000,3000};

	int32       waitlength[17]={3000,3000,3000,3000, 500,500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500};

    switch(commandID)
    {
	    case COMMAND_START :
	{

			DemoModuleDlg->FrameCount=0;
			DemoModuleDlg->PostMessage(MSG_REFRESH);
				
		    for(i=0; i<nop*10; i++){
	             data[i] = sin(i*2*PI*freq[num]/nop)+sin(i*2*PI*7/nop);
//				  data[i] =4.5+ sin(i*2*PI*freq[num]/nop);           //diotest
               }



			if (num!=0){

//				DAQmxReadAnalogF64(taskHandleread,20000,0,DAQmx_Val_GroupByScanNumber,dataread,40000,&read,NULL);
//				DAQmxStopTask(taskHandleread);
//				for(j = 0; j<iframe; j+=2){
//				    fprintf(fp,"%d  %f  %d\n",j, dataread[j], dataread[j+1]);
//			        }


			    DAQmxStopTask(taskHandle);
		        DAQmxClearTask(taskHandle);

                DAQmxStopTask(taskHandleread);
				DAQmxClearTask(taskHandleread);
				fclose(fp);
			}

			nof = datalength[num];
			now = waitlength[num];

			sprintf( filename, "I:\\Yu\\magread_%1f.txt", freq[num]);

			iframe = 0;

                /*********************************************/
	              // DAQmx Configure Code
	            /*********************************************/
			fp = fopen(filename ,"w");
			DAQmxCreateTask("",&taskHandle);
			DAQmxCreateTask("",&taskHandleread);
			DAQmxCreateAOVoltageChan(taskHandle,"Dev1/ao0","",-10.0,10.0,DAQmx_Val_Volts,NULL);
			DAQmxCreateAIVoltageChan(taskHandleread,"Dev1/ai0","",DAQmx_Val_Cfg_Default,-0.10,0.10,DAQmx_Val_Volts,NULL);
			DAQmxCreateAIVoltageChan(taskHandleread,"Dev1/ai1","",DAQmx_Val_Cfg_Default,-0.10,0.10,DAQmx_Val_Volts,NULL);
//			DAQmxCreateAIVoltageChan(taskHandleread,"Dev1/ai0","",DAQmx_Val_Cfg_Default,-10.0,10.0,DAQmx_Val_Volts,NULL);
			DAQmxCfgSampClkTiming(taskHandle,"",10000.0,DAQmx_Val_Rising,DAQmx_Val_ContSamps,10000);
//			DAQmxCfgSampClkTiming(taskHandleread,"/Dev1/PFI0",100.0,DAQmx_Val_Rising,DAQmx_Val_ContSamps,2);
            DAQmxCfgSampClkTiming(taskHandleread,"/Dev1/PFI0",100.0,DAQmx_Val_Rising,DAQmx_Val_ContSamps,42000);
//			DAQmxCfgDigEdgeStartTrig(taskHandleread,"/Dev1/PFI0",DAQmx_Val_Rising);
			DAQmxWriteAnalogF64(taskHandle,100000,0,10.0,DAQmx_Val_GroupByChannel,data,NULL,NULL);
			
	            /*********************************************/
	            // DAQmx Start Code
	            /*********************************************/
			DAQmxStartTask(taskHandle);
//			DAQmxStartTask(taskHandleread);


			num+=1;

			return true;
		}
        case COMMAND_INCREMENT :
        {
		  if (num!=0){
            
			
			if (iframe == 0){
				DAQmxStartTask(taskHandleread);
			}
		    
//			DAQmxReadAnalogF64(taskHandleread,1,0.01,DAQmx_Val_GroupByScanNumber,dataread,2,&read,NULL);
           
			iframe = iframe + 1;

//			fprintf(fp,"%f  %f  %d\n", dataread[0], dataread[1], read );



            DemoModuleDlg->FrameCount++;
            DemoModuleDlg->PostMessage(MSG_REFRESH);


//			if (iframe == 0 ){
//				DAQmxStartTask(taskHandleread);
//			}
		


			if (iframe == nof){
                 DAQmxReadAnalogF64(taskHandleread,(nof+now),0,DAQmx_Val_GroupByScanNumber,dataread,2*(nof+now),&read,NULL);
			     for(j = 0; j<(2*(nof)); j+=2){
				    fprintf(fp,"%f  %f  %d\n", dataread[j], dataread[j+1], read );
			        }
				 DAQmxStopTask(taskHandleread);
				 DAQmxClearTask(taskHandleread);
				 fclose(fp);
			}
			}
            return true;
		  
        }
		case COMMAND_STOP :
		{
			num = 0;    
			

//		   DAQmxReadAnalogF64(taskHandleread,-1,-1,DAQmx_Val_GroupByChannel,dataread,40000,&read,NULL);
           DAQmxStopTask(taskHandle);
		   DAQmxStopTask(taskHandleread);
//				for(j = 0; j<iframe; j+=2){
//				    fprintf(fp,"%d  %f  %f\n",j, dataread[j], dataread[j+1]);
//			        }
		   DAQmxClearTask(taskHandle);
		   DAQmxClearTask(taskHandleread);

		   fclose(fp);



			return true;
		}
    }

return false;
}

//-----------------------------------------------------------------------------
int CDemoModuleModule::GetAutoLoadedCommandCount()
{
    return 1;
}
//-----------------------------------------------------------------------------
CAutoLoadedCommand CDemoModuleModule::GetAutoLoadedCommand(int index)
{
    if(index == 0)
    {
        CString command,command1,command2;
        command.Format(_T("%s:Increment"), Name);
		command1.Format(_T("%s:Start"), Name);
		command2.Format(_T("%s:Stop"), Name);
        return CAutoLoadedCommand(command1, _T("StreampixCore"), _T("On PostCreate AVI"));


    }

    return CAutoLoadedCommand();
}
//-----------------------------------------------------------------------------
void CDemoModuleModule::TriggerEvent()
{
//  int32   i=0,j=0;
 // int32   nop = 5000;

           num = 0;
           DAQmxStopTask(taskHandle);
		   DAQmxClearTask(taskHandle);
  



}
