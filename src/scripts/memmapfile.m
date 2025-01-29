## Copyright (C) 2018 Guillaume Flandin
##
## This file is part of Octave.
##
## Octave is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {} {@var{m} =} memmapfile (@var{name})
## @deftypefnx {} {@var{m} =} memmapfile (@var{name}, @qcode{"Property"}, @var{value}, @dots{})
##
## Create an object of the memmapfile class allowing to access the contents of a
## file through memory mapping.
##
## @var{name} is the name of the file to be memory mapped.  It must exist at the
## time the object is created.
## 
## The memory mapping actually takes place the first time the contents of the
## file is accessed, and lasts until the object is destroyed or its settings
## are modified.
##
## Options can be provided by a list of property-value pairs:
## @table @asis
##
## @item @qcode{"Writable"}
## Specifies whether the contents of the mapped file can be modified.  Default
## is false.
##
## @item @qcode{"Offset"}
## Specifies the number of bytes to skip from the beginning of the file.  Any
## non-negative integer value is acceptable, and does not have to be a multiple
## of the page size.  Defaults is 0.
##
## @item @qcode{"Format"}
## Specifies the format of the data elements to be mapped to the file contents.
## If @qcode{"Format"} is a char array, it has to be one of @qcode{""uint8"},
## @qcode{"int8"}, @qcode{"uint16"}, @qcode{"int16"}, @qcode{"uint32"},
## @qcode{"int32"}, @qcode{"uint64"}, @qcode{"int64"}, @qcode{"single"} or
## @qcode{"double"}.  Default is @qcode{""uint8"}.
##
## Alternatively @qcode{"Format"} can be defined by a Nx3 cell array where each
## row has the form @{@var{prec}, @var{dim}, @var{name}@}.
## @var{prec} is one of the data types listed above.  @var{dim} is a dimension
## array.  @var{name} is the name of the field by which this data element will
## be accessed.
##
## For example, @qcode{"Format"} can be:
## @example
##         @{"uint16", [4, 3], "x";
##          "single",     16, "y"@}
## @end example
##
## @item @qcode{"Repeat"}
## Specifies the number of times the data (as defined by @qcode{"Format"}) have
## expected to be present on disk.  Default is Inf, meaning there will be as
## many data elements as available according to the size of the file.  If more
## data elements are specified than what the size of the file allows, none will
## be available.
##
## @end table
## @end deftypefn

## -*- texinfo -*-
## @deftypefn  {} {@var{name} =} m.Filename
## @deftypefnx {} {m.Filename =} @var{name}
## Return or specify the name of the file.
## @end deftypefn
## 
## @deftypefn {} {@var{tf} =} m.Writable
## @deftypefnx {} {m.Writable =} @var{tf}
## Return or specify whether the data can be modified.
## @end deftypefn
## 
## @deftypefnx {} {@var{offset} =} m.Offset
## @deftypefnx {} {m.Offset =} @var{offset}
## Return or specify the number of bytes to use as an offset in the file.
## @end deftypefn
##
## @deftypefnx {} {@var{fmt} =} m.Format
## @deftypefnx {} {m.Format =} @var{fmt}
## Return or specify the format of the data to be memory mapped.
## @end deftypefn
##
## @deftypefnx {} {@var{repeat} =} m.Repeat
## @deftypefnx {} {m.Repeat =} @var{repeat}
## Return or specify how many data elements should be mapped.
## @end deftypefn

## -*- texinfo -*-
## The mapped data will be accessible through the m.Data interface.
##
## If @qcode{"Format"} is a char array, m.Data will behave as if it were a
## @var{N}x1 numeric array of type @qcode{"Format"} where @var{N} is either
## @qcode{"Repeat"} or the maximal number of data elements the mapped area can
## store if @qcode{"Repeat"} is Inf or 0 if @qcode{"Repeat"} is larger than what
## the mapped area can store.
##
## If @qcode{"Format"} is a cell array, m.Data will behave as an @var{N}x1
## struct array whose field names are specified by the third column of
## @qcode{"Format"} and @var{N} obtained through the same mechanism than above.
## The data type and dimension of the numeric array of each field of the
## structure are defined by the first and second column of the @qcode{"Format"}
## cell array.
##
## For example:
## @example
## m = memmapfile ("bigdata.dat");
## m.Data(42)
## 
## m = memmapfile ("bigdata.dat", "Format", {"single", [4, 3], "speed"}, "Writable", true);
## m.Data(7).speed(1:2,:) = zeros (2, 3); 
## @end example

classdef memmapfile
  
  properties
    Filename = "";
    Writable = false;
    Offset   = 0;
    Format   = "uint8";
    Repeat   = Inf;
  endproperties
  
  properties (Dependent)
    Data;
  endproperties
  
  properties (Private)
    MemoryMap;
    Filesize = 0;
  endproperties
  
  methods (Access = public)

    ## Constructor
    function this = memmapfile (name, varargin)
      if (nargin == 0)
        error ("memmapfile: Invalid call to memmapfile.");
        # print_usage ();
      endif
      
      this.Filename = name;
      
      if (mod (numel (varargin), 2))
        error ("memmapfile: Property without a value.");
      endif
      for i = 1:2:numel (varargin)
        if (! ismember (varargin{i}, {"Writable", "Offset",  "Format", "Repeat"}))
          error ("memmapfile: Unknown property.");
        endif
        this.(varargin{i}) = varargin{i+1};
      endfor
      
      this.MemoryMap = __memmapfile_handle__ ();
    endfunction
    
    ## Property Setter and Getter Methods
    function name = get.Filename (this)
      name = this.Filename;
    endfunction
    
    function this = set.Filename (this, name)
      if (! ischar (name) || size (name, 1) != 1)
        error ("memmapfile: Property Filename is a file name.");
      endif
      # name = file_in_loadpath (name);
      name = make_absolute_filename (name);
      [info, err, msg] = stat (name);
      if (err != 0)
        error ("memmapfile: %s.", msg);
      endif
      if (! S_ISREG (info.mode)) # isfile
        error ("memmapfile: No such file.");
      endif
      this.Filesize = info.size;
      this.Filename = name;
    endfunction
    
    function writable = get.Writable (this)
      writable = this.Writable;
    endfunction
    
    function this = set.Writable (this, writable)
      if (! isscalar (writable) || ! islogical (writable))
        error ("memmapfile: Property Writable is a logical scalar.");
      endif
      this.Writable = writable;
    endfunction
    
    function offset = get.Offset (this)
      offset = this.Offset;
    endfunction
    
    function this = set.Offset (this, offset)
      if (! isscalar (offset) || ! isnumeric (offset) || ! isfinite (offset) || ! isreal (offset) || offset < 0 || offset != floor (offset))
        error ("memmapfile: Property Offset is a non-negative scalar integer.");
      endif
      this.Offset = offset;
    endfunction
    
    function fmt = get.Format (this)
      fmt = this.Format;
    endfunction
    
    function this = set.Format (this, fmt)
      t = available_types ();
      sts = true;
      if (ischar (fmt))
        if (! t.isKey (fmt))
          sts = false;
        endif
      elseif (iscell (fmt))
        if (size (fmt, 2) != 3 || size (fmt,1) < 1)
          sts = false;
        else
          for i = 1:size (fmt, 1)
            if (! ischar (fmt{i,1}) || ! t.isKey (fmt{i,1}))
              sts = false;
            endif
            if (! isnumeric (fmt{i,2}) || ! isrow (fmt{i,2}) || any (fmt{i,2} <= 0) || any (! isfinite (fmt{i,2})) || any (fmt{i,2} != floor (fmt{i,2})))
              sts = false;
            endif
            if (! ischar (fmt{i,3}))
              sts = false;
            endif
          endfor
          if (numel (unique (fmt(:,3))) != size (fmt,1))
            sts = false;
          endif
        endif
      else
        sts = false;
      endif
      if (! sts)
        error ("memmapfile: Property Format is invalid.");
      endif
      this.Format = fmt;
    endfunction
    
    function repeat = get.Repeat (this)
      repeat = this.Repeat;
    endfunction
    
    function this = set.Repeat (this, repeat)
      if (! isscalar (repeat) || ! isnumeric (repeat) || isnan (repeat) || ! isreal (repeat) || repeat <= 0 || repeat != floor (repeat))
        error ("memmapfile: Property Repeat is a positive scalar integer.");
      endif
      this.Repeat = repeat;
    endfunction
    
    ## Indexed Reference
    function v = subsref (this, s)
      if (s(1).type(1) != ".")
        error  ("memmapfile: Not a valid indexing expression.");
      endif
      switch (s(1).subs)
        
        case {"Filename", "Writable", "Offset",  "Format", "Repeat"}
          v = builtin ("subsref", this.(s(1).subs), s(2:end));
          
        case "Data"
          mmap (this);
          [n, b, c] = get_dimensions (this.Filesize - this.Offset, this.Format, this.Repeat);
          if (ischar (this.Format))
            # m.Data or m.Data(idx)
            if (numel (s) == 1)
              s(2) = substruct ("()",  {":"});
            endif
            if (! strcmp (s(2).type, "()"))
              error ("memmapfile: Not a valid indexing expression.");
            endif
            v = get_data (this.MemoryMap, this.Format, 0, [n, 1], s(2).subs);
            v = builtin ("subsref", v, s(3:end));
          else
            # m.Data or m.Data.field or m.Data.field(idx) or
            # m.Data(idx) or m.Data(idx).field or m.Data(idx).field(jdx)
            M = (1:n)'; # For the first (idx)(...)
            while (numel (s) > 1 && strcmp (s(2).type, "()"))
              M = subsref (M, s(2));
              s(2) = [];
            endwhile
            if (numel (s) == 1)
              # m.Data or m.Data(idx)
              for i = 1:numel(M)
                for j = 1:size (this.Format, 1)
                  dim = this.Format{j,2};
                  if numel (dim) == 1
                    dim = [dim, 1];
                  endif
                  offset = (M(i) - 1) * b + c(j);
                  v(i,1).(this.Format{j,3}) = get_data (this.MemoryMap, this.Format{j,1}, offset, dim);
                endfor
              endfor
            else
              # m.Data.field or m.Data.field(idx) or
              # m.Data(idx).field or m.Data(idx).field(jdx)
              if (numel (s) > 1 && strcmp (s(2).type, "."))
                if (numel (M) > 1)
                  error ("memmapfile: Not a valid indexing expression.");
                endif
                idx = ismember (this.Format(:,3), s(2).subs);
                if (! any (idx))
                  error ("memmapfile: Not a valid field '%s'.", s(2).subs);
                endif
                dim = this.Format{idx,2};
                if numel (dim) == 1
                  dim = [dim, 1];
                endif
                offset = (M - 1) * b + c(idx);
                if (numel (s) > 2)
                  if (! strcmp (s(3).type, "()"))
                    error ("memmapfile: Not a valid indexing expression.");
                  endif
                  subs = s(3).subs;
                else
                  subs = [];
                endif
                v = get_data (this.MemoryMap, this.Format{idx,1}, offset, dim, subs);
                v = builtin ("subsref", v, s(4:end));
              else
                error ("memmapfile: Not a valid indexing expression.");
              endif
            endif
          endif
          
        otherwise
          error  ("memmapfile: No such property.");
          
      endswitch
    endfunction
    
    ## Indexed Assignment
    function this = subsasgn (this, s, v)
      if (s(1).type(1) != ".")
        error  ("memmapfile: Not a valid indexing expression.");
      endif
      switch (s(1).subs)
        case {"Filename", "Writable", "Offset",  "Format", "Repeat"}
          old_value = this.(s(1).subs);
          this.(s(1).subs) = builtin ("subsasgn", this.(s(1).subs), s(2:end), v);
          if (! isequal (old_value, this.(s(1).subs)))
            munmap (this);
          endif
        case "Data"
          if (! this.Writable)
            error ("memmapfile: File access is set to read-only. Change the Writable property if needed.");
          endif
          mmap (this);
          [n, b, c] = get_dimensions (this.Filesize - this.Offset, this.Format, this.Repeat);
          if (ischar (this.Format))
            # m.Data = v; or m.Data(idx) = v;
            if (numel (s) == 1)
              s(2) = substruct ("()",  {":"});
            endif
            if (! strcmp (s(2).type, "()") || numel (s) > 2)
              error ("memmapfile: Not a valid indexing expression.");
            endif
            set_data (this.MemoryMap, this.Format, 0, [n, 1], s(2).subs, v);
          else
            # m.Data.field = v; or m.Data.field(idx) = v; or
            # m.Data(idx).field = v; or m.Data(idx).field(jdx) = v;
            # Matlab does not support: m.Data = v; or m.Data(idx) = v;
            # where v would be a struct
            if (numel (s) < 2)
              error ("memmapfile: Not a valid indexing expression.");
            endif
            M = (1:n)'; # For the first (idx)
            if (strcmp (s(2).type, "()"))
              M = subsref (M, s(2));
              s(2) = [];
            endif
            if (numel (M) != 1 || numel (s) < 2 || s(2).type(1) != ".")
              error ("memmapfile: Not a valid indexing expression.");
            endif
            idx = ismember (this.Format(:,3), s(2).subs);
            if (!any (idx))
              error ("memmapfile: Not a valid field '%s'.", s(2).subs);
            endif
            dim = this.Format{idx, 2};
            if numel (dim) == 1
              dim = [dim, 1];
            endif
            offset = (M - 1) * b + c(idx);
            if (numel (s) == 2)
              subs = {":"};
            elseif (numel (s) == 3 && strcmp (s(3).type, "()"))
              subs = s(3).subs;
            else
              error ("memmapfile: Not a valid indexing expression.");
            endif
            set_data (this.MemoryMap, this.Format{idx, 1}, offset, dim, subs, v);
          endif
        otherwise
          error  ("memmapfile: No such property.");
      endswitch
    endfunction
    
    ## Disp Method
    function disp (this)
      b2s = @(x) ifelse (any (x), "true", "false");
      if (ischar (this.Format))
        fmt = sprintf ("\"%s\"", this.Format);
      else
        fmt = "";
        for i = 1:size (this.Format, 1)
          fmt = [fmt, blanks(16), sprintf("\"%s\" [%s] \"%s\"\n", ...
            this.Format{i,1}, num2str (this.Format{i,2}), this.Format{i,3})];
        endfor
        fmt = ["{" fmt(17:end-1) "}"];
      endif
      
      n = get_dimensions (this.Filesize - this.Offset, this.Format, this.Repeat);
      if (ischar (this.Format))
        dat = sprintf ("%dx1 %s array", n, this.Format);
      else
        dat = sprintf("%dx1 struct array with fields:\n", n);
        for i = 1:size (this.Format, 1)
          dat = [dat, blanks(16), this.Format{i,3}, sprintf("\n")];
        endfor
        dat(end) = [];
      endif
      
      printf ("  memmapfile object with properties:\n\n");
      printf (["    Filename : \"%s\"\n" ...
               "    Writable : %s\n" ...
               "    Offset   : %d\n" ...
               "    Format   : %s\n" ...
               "    Repeat   : %d\n" ...
               "    Data     : %s\n\n"],
               this.Filename, b2s (this.Writable), this.Offset, ...
               fmt, this.Repeat, dat);
    endfunction
    
    ## Save process for memmapfile object
    function s = saveobj (this)
      s = this;
      s.MemoryMap = 0;
      s.Filesize = 0;
    endfunction
    
    ## Disabled Methods
    function newobj = horzcat (varargin)
      error ("memmapfile: horizontal concatenation is not allowed.");
    endfunction
    
    function newobj = vertcat (varargin)
      error ("memmapfile: vertical concatenation is not allowed.");
    endfunction
    
    function s = struct (this)
      error ("memmapfile: conversion into a struct is not allowed.");
    endfunction
    
  endmethods
  
  methods (Static)
    ## Load process for memmapfile object
    function this = loadobj (s)
      if isstruct(s)
        this = memmapfile (s.Filename);
        this.Writable = s.Writable;
        this.Offset = s.Offset;
        this.Format = s.Format;
        this.Repeat = s.Repeat;
      else
        this = s;
      endif
     endfunction
  endmethods
  
  methods (Access = private)
    function mmap (this)
      if (! is_mmapped (this.MemoryMap))
        [n, b] = get_dimensions (this.Filesize - this.Offset, this.Format, this.Repeat);
        if (! n)
          error ("memmapfile: File is too small for the requested memory mapping.");
        endif
        mmap (this.MemoryMap, this.Filename, n * b, this.Offset, this.Writable);
      endif
    endfunction
    
    function munmap (this)
      munmap (this.MemoryMap);
    endfunction
  endmethods
    
endclassdef

## List of available data types and their associated number of bytes
function m = available_types ()
  m = containers.Map(...
    {"uint8", "int8", "uint16", "int16", "uint32", "int32", "uint64", "int64", "single", "double"},...
    [1 1 2 2 4 4 8 8 4 8]);
endfunction

## Size and dimension of memory mapped variables
# n: Number of elements
# b: Size of one element in bytes
# c: Offsets within one element in bytes
function [n, b, c] = get_dimensions (bytes, fmt, repeat)
  t = available_types ();
  if (ischar (fmt))
    b = t(fmt);
    c = 0;
  else
    c = zeros (1, size(fmt, 1));
    for i = 1:size(fmt, 1)
      c(i) = t(fmt{i,1}) * prod(fmt{i,2});
    endfor
    b = sum (c);
    c = cumsum ([0 c]);
  endif
  n = floor (bytes / b);
  if (n < 0)
    n = 0;
  elseif (! isinf (repeat))
    if (n >= repeat)
      n = repeat;
    else
      n = 0;
    endif
  endif
endfunction

## Test memmapfile with char array Format
%!test
%! unwind_protect
%!   mmapfile = tempname ();
%!   D = uint8 (1:16)';
%!   fid = fopen (mmapfile, "w");
%!   fwrite (fid, D, "uint8");
%!   fclose (fid);
%!   m = memmapfile (mmapfile);
%!   assert (m.Filename, mmapfile);
%!   assert (m.Writable, false);
%!   assert (m.Offset, 0);
%!   assert (m.Format, "uint8");
%!   assert (m.Repeat, Inf);
%!   assert (m.Data, D);
%!   assert (m.Data(1), D(1));
%!   assert (m.Data(3:5), D(3:5));
%!   assert (m.Data(3:5)(2), D(4));
%!   assert (m.Data(6,1), D(6));
%!   assert (m.Data(end), D(end));
%!   m.Writable = true;
%!   assert (m.Writable, true);
%!   m.Offset = 10;
%!   assert (m.Offset, 10);
%!   m.Format = "single";
%!   assert (m.Format, "single");
%!   m.Format = {"uint64", 5, "X"};
%!   assert (m.Format, {"uint64", 5, "X"});
%!   m.Format = {"uint64", 5, "X"; "double", 6, "Y"};
%!   assert (m.Format, {"uint64", 5, "X"; "double", 6, "Y"});
%!   m.Repeat = 16;
%!   assert (m.Repeat, 16);
%!   m = memmapfile (mmapfile, "Offset", 2, "Format", "uint32", "Repeat", 1);
%!   assert (m.Data, typecast (D(3:6), "uint32"));
%!   assert (m.Data(1), typecast (D(3:6), "uint32"));
%!   m = memmapfile (mmapfile, "Format", "uint16", "Repeat", 8, "Writable", true);
%!   m.Data(4:6) = 1:3;
%!   assert (m.Data(4:6), uint16 (1:3)');
%!   d = typecast (D, "uint16"); d(4:6) = 1:3;
%!   assert (m.Data, d);
%!   m = memmapfile (mmapfile, "Format", "uint16", "Offset", 8, "Repeat", 4);
%!   m.Writable = true;
%!   m.Data(3) = 5;
%!   assert (m.Data(3), uint16 (5));
%!   d(7) = 5;
%!   assert (m.Data, d(5:end));
%! unwind_protect_cleanup
%!   %try, fclose (fid); end_try_catch
%!   delete (mmapfile);
%! end_unwind_protect

## Test memmapfile with cell array Format
%! unwind_protect
%!   mmapfile = tempname ();
%!   D = (1:5000)';
%!   fid = fopen (mmapfile, "w");
%!   fwrite (fid, D, "double");
%!   fclose (fid);
%!   m = memmapfile (mmapfile, "Format", {"uint16", [5, 8], "X";  "double", [4, 5], "Y" });
%!   assert (size (m.Data), [166 1]);
%!   assert (isstruct (m.Data));
%!   assert (fieldnames (m.Data)(:)', {"X","Y"});
%!   assert (size (m.Data(1).X), [5, 8]);
%!   assert (class (m.Data(2).X), "uint16");
%!   assert (size (m.Data(1).X(4:6)), [1, 3]);
%!   assert (size (m.Data(1).X(3,5)), [1, 1]);
%!   assert (size (m.Data(1).X(3:4,5:2:8)), [2, 2]);
%!   assert (size (m.Data(3).Y), [4, 5]);
%!   assert (class (m.Data(4).Y), "double");
%!   assert (m.Data(24).X(7), uint16 (40960));
%!   assert (m.Data(42).Y(9), 1249);
%!   m.Offset = 1200;
%!   assert (size (m.Data), [161 1]);
%!   assert (m.Data(24).X(7), uint16 (20480));
%!   m.Format{2,1} = "uint32";
%!   assert (size (m.Data), [242 1]);
%!   assert (m.Data(42).Y(12), uint32 (1083101184))
%!   m.Writable = true;
%!   d = reshape (uint16 (1:5*8), 5, 8);
%!   d(3, 5) = 100;
%!   d(38) = 200;
%!   m.Data(24).X = d(:);
%!   m.Data(24).X(3,5) = 100;
%!   m.Data(24).X(38) = 200;
%!   assert (m.Data(24).X, d);
%! unwind_protect_cleanup
%!   %try, fclose (fid); end_try_catch
%!   delete (mmapfile);
%! end_unwind_protect

## Test input validation
%!test
%! mmapfile = canonicalize_file_name (which ("figure"));
%!error memmapfile ();
%!error memmapfile (tempname ());
%!error memmapfile (1);
%!error memmapfile (fileparts (mmapfile ()));
%!error memmapfile (mmapfile (), "Writable");
%!error memmapfile (mmapfile (), "NonExisting");
%!error memmapfile (mmapfile (), "Filename");
%!error memmapfile (mmapfile (), "Filename", "file.txt");
%!error memmapfile (mmapfile (), "Writable", "Unknown");
%!error memmapfile (mmapfile (), "Writable", 2);
%!error memmapfile (mmapfile (), "Offset", "Unknown");
%!error memmapfile (mmapfile (), "Offset", pi ());
%!error memmapfile (mmapfile (), "Offset", -2);
%!error memmapfile (mmapfile (), "Offset", [0 42]);
%!error memmapfile (mmapfile (), "Format", 0);
%!error memmapfile (mmapfile (), "Format", "Unknown");
%!error memmapfile (mmapfile (), "Format", {"uint16"});
%!error memmapfile (mmapfile (), "Format", {"uint42",2,"X"});
%!error memmapfile (mmapfile (), "Format", {"uint64","a","X"});
%!error memmapfile (mmapfile (), "Format", {"uint64",[4,3],5});
%!error memmapfile (mmapfile (), "Format", {"uint16",[2,3],"X"}');
%!error memmapfile (mmapfile (), "Format", {"uint8",1,"X";"uint16",2,"X"});
%!error memmapfile (mmapfile (), "Repeat", "Unknown");
%!error memmapfile (mmapfile (), "Repeat", pi ());
%!error memmapfile (mmapfile (), "Repeat", -2);
%!error memmapfile (mmapfile (), "Repeat", 0);
%!error memmapfile (mmapfile (), "Repeat", [0 42]);
