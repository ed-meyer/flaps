//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Import/Export of ASCII files of various formats:
//   Nastran output4 (.op4)
//   Universal file (.uf)
//   Matrix Market (.mm)
//   xfig (.fig)

#ifndef EXIM_H
#define EXIM_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "Curve.h"
#include "lti.h"
#include "matrix.h"
#include "text.h"

// NASTRAN Output4
namespace Output4 {

std::vector<Matrix*> importer(std::string const& path, const std::string& output="");
bool exporter(const std::string& path, const std::vector<Matrix*>& ml, bool append=false);
} // namespace output4

// Universal file
namespace UF {
// XXX deprecated:
std::vector<Matrix*> importer (std::string const& path, const std::string& output="");
bool importer (const std::string& path, std::vector<int>& nodenumbers,
	std::vector<double>& coords, std::vector<int>& segments,
	std::vector<std::complex<double>>& gct);
bool exporter(const std::string& path,
		const std::vector<Matrix*>& ml, bool append=false);
bool exporter (const std::string& path, const std::vector<int> nodes,
		const std::vector<double>& coords, const std::vector<int>& conn,
		const std::vector<double>& gct);
bool exporter (const std::string& path, const std::vector<int> nodes,
		const std::vector<double>& coords, const std::vector<int>& conn,
		const std::vector<std::complex<double>>& gct);
// convert from UF tracelines (aka "penlift") format to OpenGL segments...
std::vector<int> penlift2segments(const std::vector<int>& penlift);
// ... or to 0b indexes into the nodes array
std::vector<int> penlift2segidx(const std::vector<int>& penlift, const std::vector<int>& nodes); 
} // namespace UF

// Matlab
namespace Matlab {
std::vector<Matrix*> importer(std::string const& path, const std::string& output="");
bool exporter(const std::string& path, const std::vector<Matrix*>& ml, bool append=false);
} // namespace Matlab


// Matrix Market: namespace MM
namespace MM {

Matrix*
importer (std::string const& path, const std::string& output="");
std::vector<double>
importer (const std::string& path, int& nr, int& nc, bool& is_complex);
// Note MM has 3 exporters, one for Matrix, one for double arrays,
// and one for complex<double> arrays. This is because matview is
// used to visualize matrices and it takes .mm files
bool
exporter (const std::string& path, const Matrix* mat);
bool
exporter (const std::string& path, const std::string& title,
		const double* cp, int nr, int nc);
bool
exporter (const std::string& path, const std::string& title,
		const std::complex<double>* cp, int nr, int nc);

} // namespace MM

// xfig (.fig) files: namespace fig
// only an exporter is provided (no importer) and it takes a vector<double>
// of nodal coordinates: (x,y,z) for each node
namespace Fig {

bool exporter (const std::string& path, std::vector<double>& coord);

} // namespace Fig

// ASCII plot files (.apf): namespace Apf
//
namespace Apf {
std::vector<Curve*> importer(std::string const& path, const std::string& output="");
// 3 diff ways to plot
bool exporter(const std::vector<Curve*>& curves,
		const std::vector<std::string>& toplot,
		const std::string& filename, bool append);
bool exporter (std::string const& cid,
		std::vector<std::pair<std::string,std::string>> const& params,
		std::vector<double> const& data,
		const std::string& path, bool append);
bool plot (const std::string& path, const std::string& cid, 
		const std::vector<std::vector<double>>& xs,
		const std::vector<std::string>& names, bool append);

} // namespace Apf

namespace LTI {
bool
importer (std::string const& path, Matrix& A, Matrix& B, Matrix& C, Matrix& D,
	std::vector<Itd>& Atd, std::vector<double>& itd, std::vector<double>& otd,
	std::vector<int>& Sidx, double rho);
}	// namespace LTI

#endif // EXIM_H
