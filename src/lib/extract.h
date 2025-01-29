//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef EXTRACT_H
#define EXTRACT_H

#include <algorithm>
#include <vector>

#include "message.h"

namespace flaps {
template<typename T>
void
extract (size_t nra, size_t nca, T const* a, size_t nrows, int const* rows,
			size_t ncols, int const* cols, size_t ldb, T* b) {
/*------------------------------------------------------------------
 * Copy elements of matrix "a" specified by int arrays "rows" &
 * "cols" to matrix "b". rows & cols are 1-based row and column numbers
 *
 * It is ok for the row and column numbers to be out of range of the
 * matrix we are extracting from - we simply ignore them
 *------------------------------------------------------------------*/
	size_t i, j;

	if (rows && cols) {
		for (i=0; i<nrows; i++) {
			if (rows[i] > 0 && rows[i] <= (int)nra) {
				for (j=0; j<ncols; j++) {
					if (cols[j] > 0 && cols[j] <= (int)nca)
						b[IJ(i,j,ldb)] = a[IJ(rows[i]-1,cols[j]-1,nra)];
				}
			}
		}
	} else if (rows) {
		for (i=0; i<nrows; i++)
			if (rows[i] > 0 && rows[i] <= (int)nra) {
				for (j=0; j<nca; j++)
					b[IJ(i,j,ldb)] = a[IJ(rows[i]-1,j,nra)];
			}
	} else if (cols) {
		for (j=0; j<ncols; j++)
			if (cols[j] > 0 && cols[j] <= (int)nca) {
				for (i=0; i<nra; i++)
					b[IJ(i,j,ldb)] = a[IJ(i,cols[j]-1,nra)];
			}
	} else {
		for (i=0; i<nra; i++) {
			for (j=0; j<nca; j++) {
				b[IJ(i,j,ldb)] = a[IJ(i,j,nra)];
			}
		}
	}
}

}  // namespace flaps


#endif
