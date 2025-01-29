//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef CUSTOM_H
#define CUSTOM_H

#include "config.h"	// for BOOST
#include <any>
#ifdef HAVE_BOOST_CORE_DEMANGLE_HPP
#include <boost/core/demangle.hpp>
#endif
#include <exception>
#include <map>
#include <string>
#include <stdexcept>
#include <typeinfo>

#include "message.h"
#include "trace.h"

// namespace Custom contains functions for maintaining pointers to
// custom functions
namespace Custom {

template<typename T>
std::string
show(T& t) {
#ifdef HAVE_BOOST_CORE_DEMANGLE_HPP
	return boost::core::demangle(typeid(t).name());
#else
	return typeid(t).name();
#endif
}

bool is_defined(const std::string& name);
std::string datatype(const std::string& name);

// note that std::any requires C++17 (or Boost)
std::map<std::string, std::any>&
map();

template<typename T>
static bool get(const std::string& name, T& t) {
	auto p = Custom::map().find(name);
	if (p == Custom::map().end())
		throw std::runtime_error(vastr("function \"",name,"\" has not been defined"));
	try {
		t = std::any_cast<T>(p->second);
	} catch(std::bad_any_cast& s) {
		throw std::runtime_error(vastr("attempt to get custom function \"",name,
					"\": wrong function type (",show(t),"), should be ",
					datatype(name)));
	}
	return true;
}

// add a function to custom.so
std::string
add (std::string const& path);

// add multiple functions to custom.so
std::vector<std::string>
add (std::vector<std::string> const& names);

std::ifstream
open();

std::ifstream
prototype();

std::vector<std::string>
make_custom(const std::vector<std::string>& paths);

} // namespace Custom

#endif // CUSTOM_H
