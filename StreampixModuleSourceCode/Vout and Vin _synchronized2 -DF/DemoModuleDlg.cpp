#include "stdafx.h"
#include "DemoModuleDlg.h"
#include "DemoModuleModule.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CDemoModuleDlg::CDemoModuleDlg(CWnd* pParent, CDemoModuleModule* module)
	: CDockedDlg(CDemoModuleDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDemoModuleDlg)
	FrameCount = 0;
	//}}AFX_DATA_INIT

    Create(CDemoModuleDlg::IDD, pParent);
    
    SetGroupName(_T("Demo Module"));
    SetName(_T("Demo Module"));
    DockPriority = 1000;
    GroupDockPriority = 0;
    
    Module = module;
}
//-----------------------------------------------------------------------------
void CDemoModuleDlg::DoDataExchange(CDataExchange* pDX)
{
	CDockedDlg::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDemoModuleDlg)
	DDX_Text(pDX, IDC_FRAME_COUNT, FrameCount);
	//}}AFX_DATA_MAP
}
//-----------------------------------------------------------------------------
BEGIN_MESSAGE_MAP(CDemoModuleDlg, CDockedDlg)
	//{{AFX_MSG_MAP(CDemoModuleDlg)
	ON_BN_CLICKED(IDC_TRIGGER_EVENT_BUTTON, OnTriggerEventButton)
	//}}AFX_MSG_MAP
    ON_REGISTERED_MESSAGE(MSG_REFRESH, OnRefresh)
	ON_BN_CLICKED(IDC_BUTTON1, &CDemoModuleDlg::OnBnClickedButton1)
END_MESSAGE_MAP()
//-----------------------------------------------------------------------------
LRESULT CDemoModuleDlg::OnRefresh(WPARAM wParam, LPARAM lParam)
{
    //Refresh the display
    UpdateData(FALSE);
    return 0;
}
//-----------------------------------------------------------------------------
void CDemoModuleDlg::OnTriggerEventButton() 
{
	Module->TriggerEvent();	
}


void CDemoModuleDlg::OnBnClickedButton1()
{
	// TODO: Add your control notification handler code here
}
