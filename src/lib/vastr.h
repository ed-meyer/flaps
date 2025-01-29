//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef VASTR_H
#define VASTR_H 1

// variadic template function to write an arbitrary
// number of arguments of any type provided each argument
// has operator<<

#include <deque>
#include <sstream>
#include <string>
#include <vector>

// operator<< for vector<T> for any type T that has operator<<
// short vectors get printed on one line, long ones with one
// element per line
template<typename Type>
std::ostream&
operator<<(std::ostream& s, const std::vector<Type>& x) {
	// put elements on a line until some max line length reached
	size_t maxlen{80};
	std::ostringstream os;
	// convert all items to strings
	if (x.empty())
		return s;
	std::vector<std::string> items;
	for (auto& xi : x) {
		os.str("");
		os << xi;
		items.push_back(os.str());
	}
	// create lines of length < maxlen
	std::vector<std::string> lines;
	size_t i{0};
	while (i < items.size()) {
		 os.str("");
		 std::string sep;
		 for (size_t j=i; j<items.size(); j++) {
			 os << sep << items[j];
			 i++;
			 sep = ", ";
			 if (os.str().size() > maxlen)
				 break;
		}
		lines.push_back(os.str());
	}
	if (lines.size() > 1) {
		s << " {\n";
		for (auto& li : lines)
			s << li << std::endl;
		s << "}";
	} else
		s << lines[0];
	return s;
}
// operator<< for vector<T*> for any type T that has operator<<
// short vectors get printed on one line, long ones with one
// element per line
template<typename Type>
std::ostream&
operator<<(std::ostream& s, const std::vector<Type*>& x) {
	if (x.size() < 6) {
		std::string sep;
		for (auto xi : x) {
			s << *xi << ", ";
			sep = ", ";
		}
	} else {
		s << " {\n";
		for (auto xi : x)
			s << *xi << std::endl;
		s << "}\n";
	}
	return s;
}

// operator<< for deque<T> for any type T that has operator<<
// short deques get printed on one line, long ones with one
// element per line
template<typename Type>
std::ostream&
operator<<(std::ostream& s, std::deque<Type> x) {
	// short vector? put on one line - XXX should roff it
	if (x.size() < 6) {
		std::string sep;
		for (auto& xi : x) {
			s << sep << xi ;
			sep = ", ";
		}
	} else {
		s << " {\n";
		for (auto& xi : x)
			s << xi << std::endl;
		s << "}\n";
	}
	return s;
}

std::string
vastr();  // final vastr (no arguments), implemented in vastr.c

template<typename T, typename... Types>
std::string
vastr(const T& firstArg, const Types&... args) {
	std::ostringstream os;
	os.str("");
	os << firstArg;
	os << vastr(args...);
	return os.str();
}

#endif // VASTR_H
