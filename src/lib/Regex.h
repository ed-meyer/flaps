//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// utilities to make std::regex easier to use

#ifndef Regex_h
#define Regex_h

#include <regex>
#include <string>

#include "message.h"

// regex_msg gives a more meaningful message regarding an exception
// from regex (regex_error). Taken from Josuttis, Standard C++ Library
// chapter 14.
// Usage:
//    try {
//    	regex pat(...
//    } catch (const regex_error& e) {
//    	cerr << "regex error: " << regex_msg(e.code()) << endl;
// Or use make_regex, e.g.
//    regex re = make_regex(pattern, regex_constants::icase);
template <typename T>
std::string regex_msg (T code) {
    switch (code) {
       case std::regex_constants::error_collate:
         return "error_collate: "
               "regex has invalid collating element name";
      case std::regex_constants::error_ctype:
        return "error_ctype: "
                 "regex has invalid character class name";
        case std::regex_constants::error_escape:
          return "error_escape: "
                "regex has invalid escaped char. or trailing escape";
       case std::regex_constants::error_backref:
         return "error_backref: "
               "regex has invalid back reference";
      case std::regex_constants::error_brack:
        return "error_brack: "
                 "regex has mismatched '[' and ']'";
        case std::regex_constants::error_paren:
          return "error_paren: "
                "regex has mismatched '(' and ')'";
       case std::regex_constants::error_brace:
         return "error_brace: "
               "regex has mismatched '{' and '}'";
      case std::regex_constants::error_badbrace:
        return "error_badbrace: "
                 "regex has invalid range in {} expression";
        case std::regex_constants::error_range:
          return "error_range: "
                "regex has invalid character range, such as '[b-a]'";
       case std::regex_constants::error_space:
         return "error_space: "
               "insufficient memory to convert regex into finite state";
      case std::regex_constants::error_badrepeat:
        return "error_badrepeat: "
                 "one of *?+{ not preceded by valid regex";
        case std::regex_constants::error_complexity:
          return "error_complexity: "
                "complexity of match against regex over pre-set level";
       case std::regex_constants::error_stack:
         return "error_stack: "
               "insufficient memory to determine regex match";
    }
     return "unknown/non-standard regex error code";
}

// create a regex...
std::regex
make_regex(const std::string& pat);
// ... ignoring case
std::regex
make_regex_ic(const std::string& pat);

#endif // Regex_h
