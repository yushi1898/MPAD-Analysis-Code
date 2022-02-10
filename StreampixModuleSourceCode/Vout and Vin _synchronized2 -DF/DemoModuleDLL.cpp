#include "stdafx.h"
#include "DemoModuleDLL.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

CDemoModuleDLL::CDemoModuleDLL()
{

}
//-----------------------------------------------------------------------------
CDemoModuleDLL::~CDemoModuleDLL()
{

}
//-----------------------------------------------------------------------------
LPCTSTR CDemoModuleDLL::GetDLLName()
{
    return _T("DemoModule.dll");
}
//-----------------------------------------------------------------------------
LPCTSTR CDemoModuleDLL::GetName()
{
    return _T("VoltageoutNin_sync_DF");
}
//-----------------------------------------------------------------------------
LPCTSTR CDemoModuleDLL::GetInfo()
{
    return _T("Demo Module.\r\nCopyright Norpix Inc. 2007 \
\r\nA demonstration module.");
}
//-----------------------------------------------------------------------------
CModule* CDemoModuleDLL::Initialize()
{
    return new CDemoModuleModule;
}
