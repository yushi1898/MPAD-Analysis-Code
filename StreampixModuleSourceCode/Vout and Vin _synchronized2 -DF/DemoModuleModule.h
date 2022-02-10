#if !defined(AFX_DEMOMODULEMODULE_H__FF178F02_A86C_492F_AF44_4937B99964AB__INCLUDED_)
#define AFX_DEMOMODULEMODULE_H__FF178F02_A86C_492F_AF44_4937B99964AB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Module.h"
#include "VarUINT.h"
#include "VarBool.h"
#include "DemoModuleDlg.h"

class CDemoModuleModule : public CModule
{
    enum COMMANDID
    {
        COMMAND_INCREMENT=0,
		COMMAND_START=1,
		COMMAND_STOP=2
    };


    enum EVENT_ID
    {
        EVENT_ON_BUTTON_CLICK=0
    };

private:
    CDemoModuleDlg* DemoModuleDlg;

    CVarBool IsComputerOnFire;

protected:
    //Sample command
    CScriptCommandInfo* CommandIncrement;
	CScriptCommandInfo* CommandStart;
	CScriptCommandInfo* CommandStop;
    
    //Sample event
    CScriptEvent* EventOnButtonClick; 

    virtual void CreateCommands();
    virtual void CreateVariables();
    virtual void CreateEvents();


public:
	CDemoModuleModule();
	virtual ~CDemoModuleModule();

    virtual void Initialize();

    virtual void PreEventScript(int eventID, PARAMETERS& localParams);
    virtual void PostEventScript(int eventID, PARAMETERS& localParams);
    virtual bool ExecuteCommand(int commandID, PARAMETERS& params, PARAMETERS& localParams);

    virtual int GetAutoLoadedCommandCount();
    virtual CAutoLoadedCommand GetAutoLoadedCommand(int index);

    void TriggerEvent();
};

#endif // !defined(AFX_DEMOMODULEMODULE_H__FF178F02_A86C_492F_AF44_4937B99964AB__INCLUDED_)
