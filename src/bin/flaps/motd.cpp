//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <iostream>
#include <string>
#include <unistd.h>

#include "config.h"  // for VERSION
#include "text.h"

using namespace std;

std::string
timedate () {
// Returns string containing the current local time
// & date (without the newline added by asctime())
	time_t t = time(0);
	string rval(asctime(localtime(&t)));

	string::size_type idx = rval.find('\n');
	if (idx != string::npos)
		rval = rval.substr(0,idx);
	return rval;
}

std::string banner = "    -------------------------------------------------------------------\n\n"
"   FFFFFFFFFFFFFFFFFFFFFFlllllll                                                        \n"
"   F::::::::::::::::::::Fl:::::l                                                        \n"
"   F::::::::::::::::::::Fl:::::l                                                        \n"
"   FF::::::FFFFFFFFF::::Fl:::::l                                                        \n"
"     F:::::F       FFFFFF l::::l   aaaaaaaaaaaaa   ppppp   ppppppppp       ssssssssss   \n"
"     F:::::F              l::::l   a::::::::::::a  p::::ppp:::::::::p    ss::::::::::s  \n"
"     F::::::FFFFFFFFFF    l::::l   aaaaaaaaa:::::a p:::::::::::::::::p ss:::::::::::::s \n"
"     F:::::::::::::::F    l::::l            a::::a pp::::::ppppp::::::ps::::::ssss:::::s\n"
"     F:::::::::::::::F    l::::l     aaaaaaa:::::a  p:::::p     p:::::p s:::::s  ssssss \n"
"     F::::::FFFFFFFFFF    l::::l   aa::::::::::::a  p:::::p     p:::::p   s::::::s      \n"
"     F:::::F              l::::l  a::::aaaa::::::a  p:::::p     p:::::p      s::::::s   \n"
"     F:::::F              l::::l a::::a    a:::::a  p:::::p    p::::::pssssss   s:::::s \n"
"   FF:::::::FF           l::::::la::::a    a:::::a  p:::::ppppp:::::::ps:::::ssss::::::s\n"
"   F::::::::FF           l::::::la:::::aaaa::::::a  p::::::::::::::::p s::::::::::::::s \n"
"   F::::::::FF           l::::::l a::::::::::aa:::a p::::::::::::::pp   s:::::::::::ss  \n"
"   FFFFFFFFFFF           llllllll  aaaaaaaaaa  aaaa p::::::pppppppp      sssssssssss    \n"
"                                                    p:::::p                             \n"
"                                                    p:::::p                             \n"
"                                                   p:::::::p                            \n"
"                                                   p:::::::p                            \n"
"                                                   p:::::::p                            \n"
"                                                   ppppppppp                            \n"
"\n\n\n";

void
motd() {
	std::cout << banner;
	std::cout << "          Version " << VERSION << "  " << timedate()
			<< "   process " << getpid() << std::endl;
	cout << "          " << get_uname("nodename") << " " <<
		get_uname("release") << endl;
	std::cout << "     -------------------------------------------------------------------\n";
}
