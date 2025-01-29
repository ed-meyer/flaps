//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef UTIL_H
#define UTIL_H 1

#include <iostream>
#include <string>
#include <time.h>
#include <wx/wx.h>
#include <wx/dir.h>
#include <wx/display.h>
#include <wx/filedlg.h>
#include <wx/filedlgcustomize.h>

#include "Regex.h"
#include "text.h"

void errormsg(const std::string& msg, bool quit);

std::ostream&
operator<<(std::ostream& s, const wxPoint& t);

std::ostream&
operator<<(std::ostream& s, const wxSize& t);

std::ostream&
operator<<(std::ostream& s, const wxBoxSizer& t);

std::ostream&
operator<<(std::ostream& s, const wxString& t);

// get one of the standard wxWidgets colors
//------------------------------------------------------------------
// AQUAMARINE         FIREBRICK            MEDIUM FOREST GREEN   RED
// BLACK              FOREST GREEN         MEDIUM GOLDENROD      SALMON
// BLUE               GOLD                 MEDIUM ORCHID         SEA GREEN
// BLUE VIOLET        GOLDENROD            MEDIUM SEA GREEN      SIENNA
// BROWN              GREY                 MEDIUM SLATE BLUE     SKY BLUE
// CADET BLUE         GREEN                MEDIUM SPRING GREEN   SLATE BLUE
// CORAL              GREEN YELLOW         MEDIUM TURQUOISE      SPRING GREEN
// CORNFLOWER BLUE    INDIAN RED           MEDIUM VIOLET RED     STEEL BLUE
// CYAN               KHAKI                MIDNIGHT BLUE         TAN
// DARK GREY          LIGHT BLUE           NAVY                  THISTLE
// DARK GREEN         LIGHT GREY           ORANGE                TURQUOISE
// DARK OLIVE GREEN   LIGHT STEEL BLUE     ORANGE RED            VIOLET
// DARK ORCHID        LIME GREEN           ORCHID                VIOLET RED
// DARK SLATE BLUE    MAGENTA              PALE GREEN            WHEAT
// DARK SLATE GREY    MAROON               PINK                  WHITE
// DARK TURQUOISE     MEDIUM AQUAMARINE    PLUM                  YELLOW
// DIM GREY           MEDIUM BLUE          PURPLE                YELLOW GREEN

wxColour getcolor(const std::string& name);

std::string buildinfo();
std::string OGLVersion();

bool isWayland();
bool isWSL();
wxSize getDisplaySize();

class FilterDialog : public wxDialog {
	wxTextCtrl* filterctrl{nullptr};
	//!! wxString filterstr;
public:
	FilterDialog(wxWindow* parent);
	wxString getfilters() { return filterctrl->GetValue(); }
};

class FileDialog : public wxDialog {
	std::vector<std::string> allfiles;
	wxTextCtrl* dirctrl{nullptr};
	wxTextCtrl* filterctrl{nullptr};
	wxListBox* filedlg{nullptr};
	wxString filterstr;
	wxString dirstr;
public:
	FileDialog(wxWindow* parent);
	std::vector<std::string> getFileNames();
};

// convert a string with "wild cards" to a regular expression
std::regex wild2rx(const std::string& wild);
// test a string against one or more regular expressions
bool containsrx(const std::vector<std::regex>& v, const std::string& s);

int uniqueId();

void snooze(double nanosec);

template<typename EventTag, typename Class,
	typename EventArg, typename EventHandler>
int
myBind(const EventTag& eventtype, void(Class::*method)(EventArg&),
		EventHandler* handler) {
// Replacement for wxWidgets Bind - the difference is this one returns
// a Window or menu item ID so you don't have to deal with creating ids.
// Bind should do this for you, returning an id that is guaranteed to be
// unique but since it doesn't all we can do is return one that is "unlikely"
// to be used. Note that in the wxWidgets docs
//    docs.wxwidgets.org/3.2/overview_windowids.html
// it says "To avoid all these restrictions, it is best to avoid using hard-coded
// IDs at all, they are not needed when using wxEvtHandler::Bind() for event
// handling (unlike with the previously used event tables)". No explanation is
// given for how to do this with Bind and IDs are included in the arguments.
// There are 2 overloads for Bind(), this being one of them. I don't use the
// other which is for non-member functions.
	int newid = uniqueId();
	handler->Bind(eventtype, method, handler, newid);
	return newid;
}

// #ifdef NEVER // nope
template<typename EventTag, typename Functor>
int
funBind(const EventTag& eventtype, Functor functor, wxEvtHandler* handler) {
// Replacement for wxWidgets Bind - the difference is this one returns
// a Window or menu item ID so you don't have to deal with creating ids.
// Bind should do this for you, returning an id that is guaranteed to be
// unique but since it doesn't all we can do is return one that is "unlikely"
// to be used. Note that in the wxWidgets docs
//    docs.wxwidgets.org/3.2/overview_windowids.html
// it says "To avoid all these restrictions, it is best to avoid using hard-coded
// IDs at all, they are not needed when using wxEvtHandler::Bind() for event
// handling (unlike with the previously used event tables)". No explanation is
// given for how to do this with Bind and IDs are included in the arguments.
// There are 2 overloads for Bind(), this being one of them. I don't use the
// other which is for non-member functions.
	int newid = uniqueId();
	handler->Bind(eventtype, functor, newid);
	return newid;
}
// #endif // NEVER // nope

std::string
sizer2string(const wxBoxSizer& sizer);

void
printsizes(const std::string& title, wxWindow* t);

void
pngwrite(const std::string& file, float* pixels, int width, int height);
//!! pngwrite(const std::string& file, unsigned char* pixels, int width, int height);

#endif // UTIL_H 1
