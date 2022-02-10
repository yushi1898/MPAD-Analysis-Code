#if !defined(AFX_DEMOMODULEDLG_H__65149BB4_5EB2_45AB_9B39_BA282D702E44__INCLUDED_)
#define AFX_DEMOMODULEDLG_H__65149BB4_5EB2_45AB_9B39_BA282D702E44__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "DockedDlg.h"
#include "resource.h"

const UINT MSG_REFRESH = RegisterWindowMessage(_T("DEMO_REFRESH"));

class CDemoModuleModule;

class CDemoModuleDlg : public CDockedDlg
{
private:
    CDemoModuleModule* Module;

public:
	CDemoModuleDlg(CWnd* pParent, CDemoModuleModule* module);   // standard constructor

//MFC STUFF
	//{{AFX_DATA(CDemoModuleDlg)
	enum { IDD = IDD_DEMOMODULE_DLG };
	UINT	FrameCount;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDemoModuleDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

protected:
	//{{AFX_MSG(CDemoModuleDlg)
    afx_msg LRESULT OnRefresh(WPARAM wParam, LPARAM lParam);
	afx_msg void OnTriggerEventButton();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButton1();
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DEMOMODULEDLG_H__65149BB4_5EB2_45AB_9B39_BA282D702E44__INCLUDED_)
