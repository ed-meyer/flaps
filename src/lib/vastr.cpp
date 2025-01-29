//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "vastr.h"

std::string
vastr() {
// lowest-level vastr for template<> version
	return "";
}

#ifdef MAIN
#include <iostream>
using namespace std;
int
main(int argc, char** argv) {
	// test operator<<(vector<T>)
	vector<string> a{"string 1", "string 2", "string 3", "string 4",
		"string 5", "string 6", "string 7", "string 8", "string 9"};
	cout << "testing vector<string>: " << a << endl;
	return 0;
}
#endif // MAIN
