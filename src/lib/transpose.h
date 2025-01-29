//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef TRANSPOSE_H
#define TRANSPOSE_H

namespace flaps {
template<typename T>
void
transpose (size_t nr, size_t nc, T const* a, T* at) {
	size_t i, j;
	for (i=0; i<nr; i++)
		for (j=0; j<nc; j++)
			at[(j+i*nc)] = a[(i+j*nr)];
}
} // namespace flaps

#endif
