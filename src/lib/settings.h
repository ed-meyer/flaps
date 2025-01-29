//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef SETTINGS_H
#define SETTINGS_H
/*
 * Settings handling classes that control the execution of Flaps
 * programs and certain functions. Settings are created from user-
 * supplied preferences and from within functions and programs at
 * run time.
 */

#include <any>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "fptype.h"
#include "message.h"
#include "Regex.h"

class Settings {

private:
	std::map<std::string, std::any> the_map;

public:
	// constructor
	Settings();
	~Settings();

	// default settings are accessed like this:
	//   Settings::defaults.get(...)
	static Settings defaults;

	// return (string) datatype of setting "name"
	std::string datatype(const std::string name);

	// is setting "name" defined?
	bool is_defined(const std::string& name) const;
	// remove a setting
	bool remove(const std::string& name);

	// write all settings to a string
	std::string toString() const;

	// to get a setting value call this with the return
	// value as the second argument
	// Returns: true if the setting is defined; t is set to the value
	//          false if it is not defined; t is unchanged
	// Throws runtimme_error if called with T the wrong datatype
	template<typename T> bool
	get(const std::string& name, T& t) const {
		if (!this->is_defined(name))
			return false;

		try {
			t = std::any_cast<T>(the_map.at(name));
			// t = boost::any_cast<T>(Settings::the_map().at(name));
		} catch(std::bad_any_cast& s) {
			throw std::runtime_error(vastr("attempting to get the value of setting (",
					name,"): wrong datatype"));
		}
		return true;
	}

	// create or change the value of setting "name" to "t"
	template<typename T> T
	set(const std::string& name, const T& t) {
		T add{t};
		// get existing value: might throw if the type changes
		T rval{};
		try {
			get(name,rval);
		} catch(std::runtime_error& s) {
			flaps::warning(vastr("changing the datatype of ",name));
		}

		Settings::the_map[name] = add;
		return rval;
	}

	// global settings are put in the environment so are accessible
	// to child processes
	static bool global (const std::string& preferences);

	friend std::ostream& operator<<(std::ostream& s, const Settings& t);

}; // class Settings

std::ostream&
operator<<(std::ostream& s, const std::any& t);

#endif  // SETTINGS_H
