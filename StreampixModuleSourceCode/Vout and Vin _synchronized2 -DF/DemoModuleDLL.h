#if !defined(AFX_DEMOMODULEDLL_H__77D507B3_D47A_47F4_A1A9_36761878C49F__INCLUDED_)
#define AFX_DEMOMODULEDLL_H__77D507B3_D47A_47F4_A1A9_36761878C49F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "ModuleDLL.h"
#include "DemoModuleModule.h"

class CDemoModuleDLL : public CModuleDLL
{
public:
	CDemoModuleDLL();
	virtual ~CDemoModuleDLL();

    virtual LPCTSTR GetDLLName();
    virtual LPCTSTR GetName();
    virtual LPCTSTR GetInfo();
    virtual CModule* Initialize();
};

#endif // !defined(AFX_DEMOMODULEDLL_H__77D507B3_D47A_47F4_A1A9_36761878C49F__INCLUDED_)
