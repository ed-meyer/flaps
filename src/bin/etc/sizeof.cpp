//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "config.h"

#include <cstring>
#include <iostream>
#include <string>

#define IS(x) (arg == x)

static bool bigendian ();

using namespace std;

int
main (int argc, char **argv) {
	size_t len, bits;

	if (argc != 2) {
		cerr << "Usage: sizeof datatype (e.g. sizeof float)\n";
		return(1);
	}

	string arg{argv[1]};

	if (IS("long") || IS("long int"))
		len = sizeof(long int);
	else if (IS("long long int"))
		len = sizeof(long long int);
	else if (IS("unsigned"))
		len = sizeof(unsigned);
	else if (IS("unsigned int"))
		len = sizeof(unsigned int);
	else if (IS("unsigned long"))
		len = sizeof(unsigned long);
	else if (IS("unsigned long int"))
		len = sizeof(unsigned long int);
	else if (IS("long long")) {
#if defined(HAVE_LONG_LONG_INT)
		len = sizeof(long long int);
#else
		fprintf (stderr, "long long int is not available on this machine\n");
		return 1;
#endif
	} else if (IS("unsigned long long")) {
#if defined(HAVE_UNSIGNED_LONG_LONG)
		len = sizeof(unsigned long long int);
#else
		fprintf (stderr, "unsigned long long is not available on this machine\n");
		return 1;
#endif
	} else if (IS("unsigned short"))
		len = sizeof(unsigned short);
	else if (IS("unsigned short int"))
		len = sizeof(unsigned short int);
	else if (IS("char"))
		len = sizeof(char);
	else if (IS("unsigned char"))
		len = sizeof(unsigned char);
	else if (IS("int"))
		len = sizeof(int);
	else if (IS("short int") || IS("short"))
		len = sizeof(short int);
	else if (IS("size_t"))
		len = sizeof(size_t);
	else if (IS("size_type"))
		len = sizeof(std::string::size_type);
	else if (IS("double"))
		len = sizeof(double);
	else if (IS("long double"))
		len = sizeof(long double);
	else if (IS("float"))
		len = sizeof(float);

	// Pointers
	else if (IS("void*"))
		len = sizeof(void*);
	else if (IS("long*") || IS("long int"))
		len = sizeof(long int*);
	else if (IS("unsigned*"))
		len = sizeof(unsigned*);
	else if (IS("unsigned int*"))
		len = sizeof(unsigned int*);
	else if (IS("unsigned long*"))
		len = sizeof(unsigned long*);
	else if (IS("unsigned long int*"))
		len = sizeof(unsigned long int*);
	else if (IS("unsigned short*"))
		len = sizeof(unsigned short*);
	else if (IS("unsigned short int*"))
		len = sizeof(unsigned short int*);
	else if (IS("char*"))
		len = sizeof(char*);
	else if (IS("unsigned char*"))
		len = sizeof(unsigned char*);
	else if (IS("int*"))
		len = sizeof(int*);
	else if (IS("short int*") || IS("short*"))
		len = sizeof(short int*);
	else if (IS("double*"))
		len = sizeof(double*);
	else if (IS("long double*"))
		len = sizeof(long double*);
	else if (IS("float*"))
		len = sizeof(float*);
	else {
		cerr << "unrecognized datatype: " << arg << endl;
		return(1);
	}
	bits = len*8;

	string endian{"little"};
	if (bigendian())
		endian = "big";
	cout << "sizeof(" << arg << ") = " << len << " bytes ("
		<< bits << " bits) " << endian << "-endian\n";
}


bool
bigendian () {
  /* Are we little or big endian?  From Harbison&Steele.  */
  union {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}
