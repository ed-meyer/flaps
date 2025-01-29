//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


#include <wx/wx.h>

#include "amvz.h"
#include "prefdialog.h"

using namespace std;

wxString PrinterStrings[] = {
	wxT("To File"),
	wxT("Other...")
};

wxString PrintStrings[] = {
	wxT("To File"),
	wxT("To Printer")
};

/*!
 * PrefDialog type definition
 */

IMPLEMENT_CLASS( PrefDialog, wxDialog )

/*!
 * PrefDialog event table definition
 */

//    EVT_BUTTON( wxID_OK, PrefDialog::OnOkClick )
BEGIN_EVENT_TABLE( PrefDialog, wxDialog )
    EVT_UPDATE_UI( ID_BG, PrefDialog::OnBgUpdate )
    EVT_BUTTON( wxID_RESET, PrefDialog::OnResetClick )
    EVT_BUTTON( wxID_HELP, PrefDialog::OnHelpClick )
END_EVENT_TABLE()

/*!
 * PrefDialog constructors
 */

PrefDialog::PrefDialog( )
{
    Init();
}

PrefDialog::PrefDialog( wxWindow* parent,
  wxWindowID id, const wxString& caption,
  const wxPoint& pos, const wxSize& size, long style )
{
    Init();

    Create(parent, id, caption, pos, size, style);
}

/// Initialisation
void
PrefDialog::
Init( ) {
    m_printer = PrinterStrings[0];
    m_print = PrintStrings[0];
    m_bg = false;
}

/*!
 * PrefDialog creator
 */

bool
PrefDialog::
Create( wxWindow* parent,
  wxWindowID id, const wxString& caption,
  const wxPoint& pos, const wxSize& size, long style )
{
    // We have to set extra styles before creating the
    // dialog
    
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);

    if (!wxDialog::Create( parent, id, caption, pos, size, style ))
        return false;

    CreateControls();
    SetDialogHelp();
	 // XXX setting validators instead of providing TransferDataToWindow()
	 // doesn't work for reasons I don't understand
    // SetDialogValidators();

    // This fits the dialog to the minimum size dictated by
    // the sizers

    GetSizer()->Fit(this);
    
    // This ensures that the dialog cannot be sized smaller
    // than the minimum size

    GetSizer()->SetSizeHints(this);

    // Centre the dialog on the parent or (if none) screen
 
    Centre();

    return true;
}

/*!
 * Control creation for PrefDialog
 */

void
PrefDialog::
CreateControls() {    
    // A top-level sizer

    wxBoxSizer* topSizer = new wxBoxSizer(wxVERTICAL);
    this->SetSizer(topSizer);
    
    // A second box sizer to give more space around the controls

    wxBoxSizer* boxSizer = new wxBoxSizer(wxVERTICAL);
    topSizer->Add(boxSizer, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    // A friendly message

    wxStaticText* descr = new wxStaticText( this, wxID_STATIC,
        wxT("Please choose a printer"),
		  wxDefaultPosition, wxDefaultSize, 0 );
    boxSizer->Add(descr, 0, wxALIGN_LEFT|wxALL, 5);

    // Spacer

    boxSizer->Add(5, 5, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    // Label for the printer text control

    wxStaticText* printerLabel = new wxStaticText ( this, wxID_STATIC,
        wxT("&Printer:"), wxDefaultPosition, wxDefaultSize, 0 );
    boxSizer->Add(printerLabel, 0, wxALIGN_LEFT|wxALL, 5);

    // A text control for the user's default printer

    //wxTextCtrl* printerCtrl = new wxTextCtrl ( this, ID_PRINTER,
	//		 wxT("Color H2"), wxDefaultPosition, wxDefaultSize, 0 );

    m_printerChoice = new wxChoice ( this, ID_PRINTER,
        wxDefaultPosition, wxSize(80, -1), WXSIZEOF(PrinterStrings),
            PrinterStrings, 0 );
    m_printerChoice->SetStringSelection(PrinterStrings[0]);
    boxSizer->Add(m_printerChoice, 0, wxGROW|wxALL, 5);

    // A horizontal box sizer to contain print, background, ?

    wxBoxSizer* junkBox = new wxBoxSizer(wxHORIZONTAL);
    boxSizer->Add(junkBox, 0, wxGROW|wxALL, 5);


    // Label for the print control

    wxStaticText* printLabel = new wxStaticText ( this, wxID_STATIC,
        wxT("&Print:"), wxDefaultPosition, wxDefaultSize, 0 );
    junkBox->Add(printLabel, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    // Create the print choice control

    wxChoice* printChoice = new wxChoice ( this, ID_PRINT,
        wxDefaultPosition, wxSize(80, -1), WXSIZEOF(PrintStrings),
            PrintStrings, 0 );
    printChoice->SetStringSelection(PrintStrings[0]);
    junkBox->Add(printChoice, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    // Add a spacer that stretches to push the bg control
    // to the right

    junkBox->Add(5, 5, 1, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    wxCheckBox* bgCheckBox = new wxCheckBox( this, ID_BG,
       wxT("&Dark Background"), wxDefaultPosition, wxDefaultSize, 0 );
    bgCheckBox ->SetValue(true);
    junkBox->Add(bgCheckBox, 0,
        wxALIGN_CENTER_VERTICAL|wxALL, 5);

    // A dividing line before the OK and Cancel buttons

    wxStaticLine* line = new wxStaticLine ( this, wxID_STATIC,
        wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
    boxSizer->Add(line, 0, wxGROW|wxALL, 5);

    // A horizontal box sizer to contain Reset, OK, Cancel and Help

    wxBoxSizer* okCancelBox = new wxBoxSizer(wxHORIZONTAL);
    boxSizer->Add(okCancelBox, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

    // The Reset button

    wxButton* reset = new wxButton( this, wxID_RESET, wxT("&Reset"),
        wxDefaultPosition, wxDefaultSize, 0 );
    okCancelBox->Add(reset, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    // The OK button

    wxButton* ok = new wxButton ( this, wxID_OK, wxT("&OK"),
        wxDefaultPosition, wxDefaultSize, 0 );
    okCancelBox->Add(ok, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    // The Cancel button

    wxButton* cancel = new wxButton ( this, wxID_CANCEL,
        wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
    okCancelBox->Add(cancel, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);

    // The Help button

    wxButton* help = new wxButton( this, wxID_HELP, wxT("&Help"),
        wxDefaultPosition, wxDefaultSize, 0 );
    okCancelBox->Add(help, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5);
}

// Set the validators for the dialog controls
bool
PrefDialog::
TransferDataToWindow() {
	wxChoice* printerChoice = (wxChoice*)FindWindow(ID_PRINTER);
	wxChoice* printChoice = (wxChoice*)FindWindow(ID_PRINT);
	wxCheckBox* bgCheckBox = (wxCheckBox*)FindWindow(ID_BG);

		printerChoice->SetSelection(0);
		printChoice->SetSelection(0);
		bgCheckBox->SetValue(m_bg);
		return true;
}
bool
PrefDialog::
TransferDataFromWindow() {
	wxChoice* printerChoice = (wxChoice*)FindWindow(ID_PRINTER);
	wxChoice* printChoice = (wxChoice*)FindWindow(ID_PRINT);
	wxCheckBox* bgCheckBox = (wxCheckBox*)FindWindow(ID_BG);

		int index = printerChoice->GetCurrentSelection();
		m_printer = PrinterStrings[index];
		index = printChoice->GetCurrentSelection();
		m_print = PrintStrings[index];
		m_bg = bgCheckBox->GetValue();
		if (m_bg) {
#ifdef NEVER // do not allow changing bg to black?
			BGRed = 0.0;
			BGGreen = 0.0;
			BGBlue = 0.0;
			BGAlpha = 0.0;
#endif // NEVER : do not allow changing bg to black?
		}
		//cout << "printer = " << m_printer << ", print = " << m_print
		//	<< ", bg = " << m_bg << endl;
		return true;
}

// Sets the help text for the dialog controls
void
PrefDialog::
SetDialogHelp() {
    wxString printerHelp = wxT("Pick a printer");
    wxString printHelp = wxT("Enter your full name.");
    wxString bgHelp = wxT("Check this if you want a dark background");

    FindWindow(ID_PRINTER)->SetHelpText(printerHelp);
    FindWindow(ID_PRINTER)->SetToolTip(printerHelp);

    FindWindow(ID_PRINTER)->SetHelpText(printerHelp);
    FindWindow(ID_PRINTER)->SetToolTip(printerHelp);

    FindWindow(ID_BG)->SetHelpText(bgHelp);
    FindWindow(ID_BG)->SetToolTip(bgHelp);
}

/*!
 * wxEVT_UPDATE_UI event handler for ID_CHECKBOX
 */

void
PrefDialog::
OnBgUpdate( wxUpdateUIEvent& event ) {
	//T("bg updated: now %s\n",m_bg?"true":"false");
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_RESET
 */

void
PrefDialog::
OnResetClick( wxCommandEvent& event ) {
    Init();
    TransferDataToWindow();
}

void
PrefDialog::
OnOkClick( wxCommandEvent& event ) {
	//wxString printer = m_printerChoice->GetValue();
	int index = m_printerChoice->GetCurrentSelection();
	cout << "got printer = " << index;
}

/*!
 * wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_HELP
 */

void
PrefDialog::
OnHelpClick( wxCommandEvent& event ) {
    // Normally we would wish to display proper online help.
    // For this example, we're just using a message box.
    /*
    wxGetApp().GetHelpController().DisplaySection(wxT("Personal record dialog"));
     */

    wxString helpText =
      wxT("Please enter your default printer.\n")
      wxT("Also indicate if you prefer a dark background.\n\n");

    wxMessageBox(helpText,
      wxT("Flaps Mode Visualizer Preferences Help"),
      wxOK|wxICON_INFORMATION, this);
}

