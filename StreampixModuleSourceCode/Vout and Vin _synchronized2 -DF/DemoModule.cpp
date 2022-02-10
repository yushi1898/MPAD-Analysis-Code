#include "stdafx.h"
#include <afxdllx.h>
#include "DemoModuleDLL.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

static AFX_EXTENSION_MODULE DemoModuleDLL = { NULL, NULL };

extern "C" int APIENTRY
DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID lpReserved)
{
	UNREFERENCED_PARAMETER(lpReserved);

	if (dwReason == DLL_PROCESS_ATTACH)
	{
		TRACE0("DEMOMODULE.DLL Initializing!\n");
		
		if (!AfxInitExtensionModule(DemoModuleDLL, hInstance))
			return 0;

		new CDynLinkLibrary(DemoModuleDLL);
	}
	else if (dwReason == DLL_PROCESS_DETACH)
	{
		TRACE0("DEMOMODULE.DLL Terminating!\n");
		AfxTermExtensionModule(DemoModuleDLL);
	}
	return 1;
}
//-----------------------------------------------------------------------------
CDemoModuleDLL MainClass; //The object single instance