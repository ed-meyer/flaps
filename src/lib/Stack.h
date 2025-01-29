//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef STACK_H
#define STACK_H
#include <deque>
#include <exception>
#include <iostream>

// Taken from Josuttis, Nicolai M., "The C++ Standard Library",
// section 12.1.3 "A User-Defined Stack Class"
template<typename T>
class Stack {
	protected:
		std::deque<T> c;  // container for the elements
	public:
		class empty_stack : public std::exception {
			public:
				virtual const char* what() const throw() {
					return "read empty stack";
				}
		};

		// number of elements
		typename std::deque<T>::size_type size() const {
			return c.size();
		}
		// is stack empty?
		bool empty() const {
			return c.empty();
		}
		// push an element onto the stack
		void push (const T& elem) {
			c.push_back(elem);
		}
		// pop an element from the stack and return its value
		T pop() {
			if (c.empty())
				throw empty_stack();
			T elem(c.back());
			c.pop_back();
			return elem;
		}
		// return value of next element
		T& top() {
			if (c.empty())
				throw empty_stack();
			return c.back();
		}
		// does the stack contain "elem"
		bool contains(const T& elem) {
			if (c.empty())
				return false;
			return (find(c.begin(),c.end(),elem) != c.end());
		}
		friend std::ostream& operator<<(std::ostream& s, const Stack<std::string>& t);
};

std::ostream&
operator<<(std::ostream& s, const Stack<std::string>& t);
			
#endif // STACK_H
