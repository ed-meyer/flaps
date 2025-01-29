//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <stdexcept>

/*
 * Tracking exception classes
 */
class Finished : public std::runtime_error {
public:
	Finished(const std::string& m) : std::runtime_error(m) {}
};
class TrackingFinished : public std::runtime_error {
public:
	TrackingFinished(const std::string& m) : std::runtime_error(m) {}
};
class NoProgress : public std::runtime_error {
public:
	NoProgress(const std::string& m) : std::runtime_error(m) {}
};
class TrackingError : public std::runtime_error {
public:
	TrackingError(const std::string& m) : std::runtime_error(m) {}
};

class NullTangent : public std::runtime_error {
public:
	NullTangent (std::string const& m) : std::runtime_error(m) {}
};

class UncertainTangent : public std::runtime_error {
public:
	UncertainTangent (std::string const& m) : std::runtime_error(m) {}
};

class RejectPrev : public std::runtime_error {
public:
	RejectPrev (std::string const& m) : std::runtime_error(m) {}
};

class IterateOutOfRange : public std::runtime_error {
public:
	IterateOutOfRange (std::string const& m) : std::runtime_error(m) {}
};

// exceptions thown by functions searching for starting
// values from source runs:
class StartSearch : public std::runtime_error {
public:
	StartSearch (const std::string& m) : std::runtime_error(m) {}
};

class DiffValue : public StartSearch {
public:
	DiffValue (std::string const& m) : StartSearch(m) {}
};

class OutOfRange : public StartSearch {
public:
	OutOfRange (std::string const& m) : StartSearch(m) {}
};

#endif // EXCEPTION_H
