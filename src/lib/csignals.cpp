//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Utilities to make dealing with C signal-handling functions
// easier and more C++-like

#include <cassert>
#include <csignal>

#include "csignals.h"
#include "trace.h"
#include "message.h"

using namespace std;

string
Signal::
desc (int sig) {
// Returns a string description of a signal, given a
// signal number
	T_(Trace trc(1,"Signal::desc(",sig,")");)
	string rval;
	
	switch (sig) {
								// Std C signals:
		case SIGABRT:
			rval = "abort";
			break;
		case SIGFPE:
			rval = "floating-point exception";
			break;
		case SIGILL:
			rval = "illegal instruction";
			break;
		case SIGINT:
			rval = "interrupt";
			break;
		case SIGSEGV:
			rval = "segmentation violation";
			break;
		case SIGTERM:
			rval = "external termination request";
			break;
							// POSIX signals:
		case SIGALRM:
			rval = "alarm";
			break;
		case SIGCHLD:
			rval = "death of a child";
			break;
		case SIGCONT:
			rval = "continue if currently stopped";
			break;
		case SIGHUP:
			rval = "hangup";
			break;
		case SIGKILL:
			rval = "killed (cannot be caught or ignored)";
			break;
		case SIGPIPE:
			rval = "write on pipe with no one to read it";
			break;
		case SIGQUIT:
			rval = "quit (generated from terminal special character)";
			break;
		case SIGSTOP:
			rval = "stop (cannot be caught or ignored)";
			break;
		case SIGTSTP:
			rval = "interactive stop";
			break;
		case SIGTTIN:
			rval = "background read attempted from control terminal";
			break;
		case SIGTTOU:
			rval = "background write attempted to control terminal";
			break;
		case SIGUSR1:
			rval ="user-defined signal 1";
			break;
		case SIGUSR2:
			rval = "user-defined signal 2";
			break;
#ifdef SIGCPULIM
		case SIGCPULIM:						/* cpu time limit */
			rval = "cpu time limit";
			break;
#endif
		default:
			rval = Signal::macro(sig);
			if (rval.empty())
				rval = vastr("signal number ",sig);
	}
	T_(trc.dprint("returning ", rval);)
	return rval;
}

string
Signal::
macro(int sig) {
// Returns the name of the macro associated with signal "sig"
	string rval;
	if (sig == 1) rval = "SIGHUP";
	if (sig == 2) rval = "SIGINT";
	if (sig == 3) rval = "SIGQUIT";
	if (sig == 4) rval = "SIGILL";
	if (sig == 5) rval = "SIGTRAP";
	if (sig == 6) rval = "SIGABRT";
	if (sig == 7) rval = "SIGERR";
	if (sig == 8) rval = "SIGFPE";
	if (sig == 9) rval = "SIGKILL";
	if (sig == 10) rval = "SIGPRE";
	if (sig == 11) rval = "SIGORE";
	if (sig == 12) rval = "SIGSYS";
	if (sig == 13) rval = "SIGPIPE";
	if (sig == 14) rval = "SIGALRM";
	if (sig == 15) rval = "SIGTERM";
	if (sig == 16) rval = "SIGIO";
	if (sig == 17) rval = "SIGURG";
	if (sig == 18) rval = "SIGCLD";
	if (sig == 19) rval = "SIGPWR";
	if (sig == 20) rval = "SIGMT";
	if (sig == 21) rval = "SIGMTKILL";
	if (sig == 22) rval = "SIGBUFIO";
	if (sig == 23) rval = "SIGRECOVERY";
	if (sig == 24) rval = "SIGUME";
	if (sig == 25) rval = "SIGDLK";
	if (sig == 26) rval = "SIGCPULIM";
	if (sig == 27) rval = "SIGSHUTDN";
	if (sig == 28) rval = "SIGSTOP";
	if (sig == 29) rval = "SIGTSTP";
	if (sig == 30) rval = "SIGCONT";
	if (sig == 31) rval = "SIGTTIN";
	if (sig == 32) rval = "SIGTTOU";
	if (sig == 33) rval = "SIGWINCH";
	if (sig == 34) rval = "SIGRPE";
	if (sig == 35) rval = "SIGWRBKPT";
	if (sig == 36) rval = "SIGNOBDM";
	if (sig == 37) rval = "SIGAMI";
	if (sig == 38) rval = "SIGSMCE";
	if (sig == 39) rval = "SIGAPTEAM";
	if (sig == 40) rval = "SIGCANCEL";
	if (sig == 41) rval = "SIGPEFAILURE";
	if (sig == 42) rval = "SIGMSP";
	if (sig == 43) rval = "SIGCRAY10";
	if (sig == 44) rval = "SIGCRAY11";
	if (sig == 45) rval = "SIGCRAY12";
	if (sig == 46) rval = "SIGCRAY13";
	if (sig == 47) rval = "SIGCRAY14";
	if (sig == 48) rval = "SIGINFO";
	if (sig == 49) rval = "SIGUSR1";
	if (sig == 50) rval = "SIGUSR2";
	return rval;
}

int
Signal::
number(const string& macro) {
	if (macro == "SIGHUP") return 1;	/* hangup				*/
	if (macro == "SIGINT") return 2;	/* interrupt				*/
	if (macro == "SIGQUIT") return 3;	/* quit					*/
	if (macro == "SIGILL") return 4;	/* illegal instruction			*/
	if (macro == "SIGTRAP") return 5;	/* trace trap				*/
	if (macro == "SIGABRT") return 6;	/* abort				*/
	if (macro == "SIGERR") return 7;	/* error exit				*/
	if (macro == "SIGFPE") return 8;	/* floating point exception		*/
	if (macro == "SIGKILL") return 9;	/* kill (cannot be caught or ignored)	*/
	if (macro == "SIGPRE") return 10;	/* program range error			*/
	if (macro == "SIGORE") return 11;	/* operand range error			*/
	if (macro == "SIGSYS") return 12;	/* bad argument to system call		*/
	if (macro == "SIGPIPE") return 13;	/* write on a pipe & no one to read it	*/
	if (macro == "SIGALRM") return 14;	/* alarm clock				*/
	if (macro == "SIGTERM") return 15;	/* software termination signal from kill*/
	if (macro == "SIGIO") return 16;	/* input/output possible signal		*/
	if (macro == "SIGURG") return 17;	/* urgent condition on I/O channel	*/
	if (macro == "SIGCLD") return 18;	/* death of a child			*/
	if (macro == "SIGPWR") return 19;	/* power-fail				*/
	if (macro == "SIGMT") return 20;	/* Cray-2 multi-tasking wakeup signal	*/
	if (macro == "SIGMTKILL") return 21;	/* Cray-2 multi-tasking kill signal	*/
	if (macro == "SIGBUFIO") return 22;	/* Reserved for CRI library use on MPP */
	if (macro == "SIGRECOVERY") return 23;	/* Recovery signal			*/
	if (macro == "SIGUME") return 24;      /* uncorrectable memory error		*/
	if (macro == "SIGDLK") return 25;	/* X-MP true deadlock detected		*/
	if (macro == "SIGCPULIM") return 26;	/* cpu limit exceeded			*/
	if (macro == "SIGSHUTDN") return 27;	/* imminent system shutdown (advisory)	*/
	if (macro == "SIGSTOP") return 28;	/* sendable stop signal not from tty	*/
	if (macro == "SIGTSTP") return 29;	/* stop signal from tty			*/
	if (macro == "SIGCONT") return 30;	/* continue a stopped process		*/
	if (macro == "SIGTTIN") return 31;	/* to readers pgrp upon background tty read */
	if (macro == "SIGTTOU") return 32;	/* like TTIN for output if selected	*/
	if (macro == "SIGWINCH") return 33;	/* window size changes			*/
	if (macro == "SIGRPE") return 34;	/* Y-MP register parity error		*/
	if (macro == "SIGWRBKPT") return 35;	/* Y-MP C90 write breakpoint trap	*/
	if (macro == "SIGNOBDM") return 36;	/* Y-MP mode binary on C90 enabled BDM  */
	if (macro == "SIGAMI") return 37;	/* T90 address multiply interrupt	*/
	if (macro == "SIGSMCE") return 38;	/* Shared memory caching error		*/
	if (macro == "SIGAPTEAM") return 39;	/* MPP Application Team abnormal exit	*/
	if (macro == "SIGCANCEL") return 40;	/* pthread cancellation			*/
	if (macro == "SIGPEFAILURE") return 41;	/* MPP Application Team killed due to dead pe */
	if (macro == "SIGMSP") return 42;	/* MSP process killed due to connect timeout */
	if (macro == "SIGCRAY10") return 43;	/* reserved for Cray Research, Inc.	*/
	if (macro == "SIGCRAY11") return 44;	/* reserved for Cray Research, Inc.	*/
	if (macro == "SIGCRAY12") return 45;	/* reserved for Cray Research, Inc.	*/
	if (macro == "SIGCRAY13") return 46;	/* reserved for Cray Research, Inc.	*/
	if (macro == "SIGCRAY14") return 47;	/* reserved for Cray Research, Inc.	*/
	if (macro == "SIGINFO") return 48;	/* quota warning or limit reached	*/
	if (macro == "SIGUSR1") return 49;	/* user defined signal 1		*/
	if (macro == "SIGUSR2") return 50;	/* user defined signal 2		*/
	return 0;
}

void
Signal::
set (int sig, void (*handler)(int)) {
/*------------------------------------------------------------------
 * set a signal-handling function for a signal
 * All signals are blocked when inside the handler
 *
 * Input:
 *  sig			signal number (see signal.h)
 *  handler		pointer to the signal-catching function
 *					or SIG_DFL or SIG_IGN
 * Throws runtime_error if error
 *------------------------------------------------------------------*/
	struct sigaction act, oact;
	string exc;

	act.sa_handler = handler;

	// set signals to be blocked inside the handler to "all signals"
	if (sigfillset(&act.sa_mask) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("sigfillset failed setting ",sig," (",desc(sig),")"));

	// the sa_flags member is used by IBM for non-POSIX stuff
	// so it is very important to set it
#ifdef SA_RESTART
	act.sa_flags = SA_RESTART;
#else
	act.sa_flags = 0;
#endif
	if (sigaction (sig, &act, &oact) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("sigaction failed setting signal ",sig," (",desc(sig),")"));
}

void
Signal::
ignore (int sig) {
// set the action for receipt of a signal to "ignore"
// Input:
// sig			signal number (see signal.h)
// Throws runtime_error exception if unsuccessful
	struct sigaction act, oact;
	string exc;

	act.sa_handler = SIG_IGN;
	act.sa_flags = 0;
	if (sigaction (sig, &act, &oact) == -1)
		throw std::system_error(errno, std::generic_category(),
			 vastr("setting signal ",sig," (",desc(sig),") to \"ignore"));
}

void
Signal::
block (int sig) {
// Block "sig" from being delivered
// Input:
// sig			signal number (see signal.h)
// Throws runtime_error exception if unsuccessful
	sigset_t set;
	string exc;

	if (sigemptyset (&set) == -1)
		throw std::system_error(errno, std::generic_category(),
			"emptying signal set");

	if (sigaddset (&set, sig) == -1)
		throw std::system_error(errno, std::generic_category(),
			 vastr("adding signal ",sig," (", desc(sig),") to set"));

	if (sigprocmask(SIG_BLOCK, &set, 0) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("trying to block signal ",sig," (", desc(sig),")"));
}

void
Signal::
unblock (int sig) {
// Unblock "sig" from being delivered
//
// Input:
// sig			signal number (see signal.h)
//
// Throws system_error exception if unsuccessful
	sigset_t set;
	string exc;

	if (sigemptyset (&set) == -1) {
		throw std::system_error(errno, std::generic_category(),
			 "emptying signal set");
	}

	if (sigaddset (&set, sig) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("adding signal ",sig," (",desc(sig),") to set"));

	if (sigprocmask(SIG_UNBLOCK, &set, 0) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("trying to block signal ",sig," (",desc(sig),")"));
}

void
Signal::
def (int sig) {
// set the action for receipt of a signal to "default"
// 
// Input:
// sig			signal number (see signal.h)
// Throws runtime_error exception if unsuccessful
	struct sigaction act, oact;
	string exc;

	act.sa_handler = SIG_DFL;
	act.sa_flags = 0;
	if (sigaction (sig, &act, &oact) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("setting signal ",sig, " (",desc(sig),") to \"ignore\""));
}

bool
Signal::
istrapped (int sig) {
// Returns TRUE if the signal has not either been given
// a handler or ignored
	struct sigaction oact;
	bool rval{false};
	string exc;

	if (sigaction (sig, 0, &oact) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("checking handler for signal ",sig," (", desc(sig),")"));

	if (oact.sa_handler == SIG_DFL)
		rval = false;
	else
		rval = true;
	return rval;
}

bool
Signal::
ispending (int sig) {
//  Returns true if the signal is pending, false otherwise
	sigset_t set;
	bool rval{false};
	string exc;

	if (sigpending (&set) < 0)
		throw std::system_error(errno, std::generic_category(),
			vastr("checking pending set for signal ", sig," (",desc(sig),")"));

	if (sigismember (&set, sig))
		rval = true;
	else
		rval = false;
	return rval;
}

void
Signal::
catchall (void(*handler)(int)) {
// Set a signal handler for most standard signals.
	set (SIGINT, handler);
	set (SIGQUIT, handler);
	set (SIGILL, handler);
	set (SIGABRT, handler);
#ifdef SIGERR
	set (SIGERR, handler);
#endif /* SIGERR */
#ifdef SIGEMT
	set (SIGEMT, handler);
#endif
	set (SIGFPE, handler);
#ifdef SIGPRE
	set (SIGPRE, handler);
#endif /* SIGPRE */
#ifdef SIGBUS
	set (SIGBUS, handler);
#endif /* SIGBUS */
#ifdef SIGORE
	set (SIGORE, handler);
#endif /* SIGORE */
#ifdef SIGSEGV
	set (SIGSEGV, handler);
#endif /* SIGSEGV */
#ifdef SIGSYS
	set (SIGSYS, handler);
#endif /* SIGSYS */
	set (SIGPIPE, handler);
	set (SIGTERM, handler);
#ifdef SIGUME
	set (SIGUME, handler);
#endif /* SIGUME */
#ifdef SIGCPULIM
	set (SIGCPULIM, handler);
#endif
#ifdef SIGINFO
	set (SIGINFO, handler);
#endif /* SIGINFO */
}

#ifdef MAIN
#undef MAIN

#include <iostream>

void
handler (int sig) {
	T_(Trace trc(1,"handler");)
	T_(trc.dprint("handler got ",sig," (",Signal::desc(sig),")");)
	Signal::def(sig);
	raise(sig);
}

int
main(int argc, char* argv[]) {
	int sig{SIGINT}; // default signal to test

	if (argc > 1)
		str2int(argv[1], sig);

	cerr << "setting signal " << sig << " (" << Signal::desc(sig)
		<< ") to a handler\n";

	// use a lambda for the handler
	Signal::set(sig, [](int sig){
		cerr << "handler got signal #" << sig << ": " << Signal::desc(sig) << endl;
		Signal::def(sig);
		raise(sig);
	});
	if (Signal::istrapped(sig))
		cerr << "signal " << sig << " (" << Signal::desc(sig)
			<< ") is trapped (as it should be)\n";
	else
		cerr << "failed trapping signal " << sig << " ("
			<< Signal::desc(sig) << endl;
	cerr << "raising " << sig << ": should call the handler\n";
	raise (sig);
}
#endif /* MAIN */
