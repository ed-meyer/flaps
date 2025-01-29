//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
/////////////////////////////////////////////////////////////////////////////
// $Id$
/////////////////////////////////////////////////////////////////////////////

#ifndef _PERSONALRECORD_H_
#define _PERSONALRECORD_H_

/*!
 * Includes
 */

#include "wx/spinctrl.h"
#include "wx/statline.h"

/*!
 * Control identifiers
 */

enum {
    ID_PREFERENCES = 11000,
    ID_PRINTER = 11001,
    ID_PRINT = 11003,
    ID_BG = 11006
};

/*!
 * PrefDialog class declaration
 */

class PrefDialog: public wxDialog {
    DECLARE_CLASS( PrefDialog )
    DECLARE_EVENT_TABLE()

public:

//// Constructors

    PrefDialog( );

    PrefDialog( wxWindow* parent,
      wxWindowID id = ID_PREFERENCES,
      const wxString& caption = wxT("Preferences"),
      const wxPoint& pos = wxDefaultPosition,
      const wxSize& size = wxDefaultSize,
      long style = wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );

    /// Member initialisation
    void Init();

    /// Creation
    bool Create( wxWindow* parent,
      wxWindowID id = ID_PREFERENCES,
      const wxString& caption = wxT("Preferences"),
      const wxPoint& pos = wxDefaultPosition,
      const wxSize& size = wxDefaultSize,
      long style = wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );

    /// Creates the controls and sizers
    void CreateControls();

    /// Sets the validators for the dialog controls
    void SetDialogValidators();
	 bool TransferDataToWindow();
	 bool TransferDataFromWindow();

    /// Sets the help text for the dialog controls
    void SetDialogHelp();

    /// Printer accessors
    void SetPrinter(const wxString& printer) { m_printer = printer; }
    wxString GetPrinter() const { return m_printer; }

    /// PRINT accessors (male = false, female = true)
    void SetPrint(wxString const& print) { m_print = print; }
    wxString GetPrint() const { return m_print; }

    /// Does the person want light bg?
    void SetBg(bool bg) { m_bg = bg; }
    bool GetBg() const { return m_bg; }

//// PrefDialog event handler declarations

    /// wxEVT_UPDATE_UI event handler for ID_BG
    void OnBgUpdate( wxUpdateUIEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_RESET
    void OnResetClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_OK
    void OnOkClick( wxCommandEvent& event );

    /// wxEVT_COMMAND_BUTTON_CLICKED event handler for wxID_HELP
    void OnHelpClick( wxCommandEvent& event );

//// PrefDialog member variables

    /// Data members
    wxString    m_printer;
    wxString    m_print;
    bool        m_bg;

	 wxChoice* m_printerChoice;
};

#endif
    // _PERSONALRECORD_H_

