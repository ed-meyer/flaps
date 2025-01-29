//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef OBJECT_COUNTER
#define OBJECT_COUNTER 1

// object counting
// from:
// @book{vandevoorde2002c++,
//   title={C++ Templates},
//   author={Vandevoorde, David and Josuttis, Nicolai M},
//   year={2002},
//   publisher={Addison-Wesley Longman Publishing Co., Inc.}
// }

template<typename CountedType>
class ObjectCounter {
private:
protected:
	static size_t count;
	static size_t npeak;
	static size_t ntotal;
	// default constructor
	ObjectCounter() {
		++ObjectCounter<CountedType>::count;
		++ObjectCounter<CountedType>::ntotal;
		ObjectCounter<CountedType>::npeak = std::max(
				ObjectCounter<CountedType>::npeak, ObjectCounter<CountedType>::count);
	}
	// copy constructor
	ObjectCounter (ObjectCounter<CountedType> const&) {
		++ObjectCounter<CountedType>::count;
		ObjectCounter<CountedType>::npeak = std::max(
				ObjectCounter<CountedType>::npeak, ObjectCounter<CountedType>::count);
	}
	// destructor
	~ObjectCounter() {
		--ObjectCounter<CountedType>::count;
	}
public:
	// return number of existing objects:
	static size_t live() {
		return ObjectCounter<CountedType>::count;
	}
	static size_t peak() {
		return ObjectCounter<CountedType>::npeak;
	}
	static size_t total() {
		return ObjectCounter<CountedType>::ntotal;
	}
};

template<typename CountedType>
size_t ObjectCounter<CountedType>::count = 0;

template<typename CountedType>
size_t ObjectCounter<CountedType>::npeak = 0;

template<typename CountedType>
size_t ObjectCounter<CountedType>::ntotal = 0;

#endif // OBJECT_COUNTER
