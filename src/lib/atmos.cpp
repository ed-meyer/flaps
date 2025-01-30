//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

#include "Ad.h"
#include "atmos.h"
#include "trace.h"

using namespace std;

constexpr double Radius{6356.766};    // polar radius of the Earth (km)
constexpr double Gmr{34.153195};    // gas constant
// ==========================================================================
Ad
geopotential(const Ad& alt) {
// given the geometric altitude in km, returns the geopotential altitude in km
	Ad rval(alt*Radius/(alt + Radius));
	return rval;
}

Ad
theta (const Ad& z) {
// theta as a function of geometric altitude z in meters
	Ad h = geopotential(z/1000.0);  // h in km
	if (h.value() < 11.0) {         // Troposphere
      return (288.15-6.5*h)/288.15;
	}
	return 216.65/288.15;
}

Ad
atmos::
delta (const Ad& z) {
// delta as a function of geometric altitude in meters
	Ad h = geopotential(z/1000.0);  // h in km
	if (h.value() < 11.0) {        // Troposphere
		Ad th = theta(z);
      return pow(th, Gmr/6.5);
	}
	return 0.2233611*exp(-Gmr*(h-11.0)/216.65);
}

Ad
sigma (const Ad& z) {
// sigma as a function of geometric altitude in meters
	Ad th = theta(z);
	Ad del = atmos::delta(z);
	return del/th;
}

Ad
atmos::
temp (Ad const& z) {
// temperature in degrees Kelvin (K)
	double tzero{288.15};      // sea level temperature, kelvins
	return tzero*theta(z);
}

Ad
atmos::
press (Ad const& z) {
	double pzero{101325.0};        // sea-level pressure, N/sq.m
	return pzero*atmos::delta(z);
}

Ad
atmos::
rho (Ad const& z) {
	double rhozero{1.225};           // sea level density, kg/cu.m
	return rhozero*sigma(z);
}

Ad
atmos::
vsonic (Ad const& z) {
	double azero{340.294};    // sea-level speed of sound, m/sec
	return azero*sqrt(theta(z));
}

double
atmos::
altmin() {
// returns the min allowable altitude in meters
	return -2000.0;
}

double
atmos::
altmax() { return 86000.0; }

double
atmos::
rhoref() {
	double rhozero{1.225};           // sea level density, kg/cu.m
	return rhozero;
}

double
atmos::
spressref() {
	double pzero{101325.0};        // sea-level pressure, N/sq.m
	return pzero;
}

Ad
atmos::
vcas(const Ad& vtas, const Ad& alt) {
// vcas = veas*(1 + 1/8(1-\delta)M^2 + 3/640(1-10\delta+9\delta^2)M^4)
// but since veas = sqrt(rho/rhoref)*vtas, rho = atmos::rho(alt),
//           M = vtas/vsound, vsound = atmos::vsound(alt)
// all we need is vtas & alt
// from wikipedia Equivalent_airspeed
//
	T_(Trace trc(3,"atmos::vcas");)

	Ad mach = vtas/atmos::vsonic(alt);
	Ad mach2{mach*mach};
	Ad rho{atmos::rho(alt)};
	Ad rhoref{atmos::rhoref()};
	Ad veas = sqrt(rho/rhoref)*vtas;
	// delta is the pressure ratio
	//!! Ad delta = atmos::press(alt)/atmos::spressref();
	Ad delt = atmos::delta(alt);

	Ad t = 1.0 + 0.125*(1.0 - delt)*mach2 +
		0.0046875*(1.0 - 10.0*delt + 9.0*delt*delt)*mach2*mach2;
	Ad rval = veas*t;
	T_(T_(trc.dprint("returning ",rval);))
	return rval;
}

Ad
atmos::
cdpress(const Ad& vtas, const Ad& alt) {
// calibrated dynamic pressure (impact pressure) as a
// function of Mach and static pressure in the standard atmosphere
// qc = spress*((1+0.2*mach^2)^3.5 - 1)
// Reference:
//   Clancy, L.J. "Aerodynamics", Wiley, 1975
//   wikipedia.org/wiki/Impact_pressure
	T_(Trace trc(3,"atmos::cdpress");)

	Ad spress = atmos::press(alt);
	Ad mach{vtas/atmos::vsonic(alt)};
	T_(T_(trc.dprint("spress ",spress,", mach ",mach);))
	if (mach.value() < 0.0) {
		T_(T_(trc.dprint("mach < 0: ",mach.value());))
		return spress;
	}
	Ad t(1.0);
	Ad rval(0.0);

	constexpr double expon{3.5};
	t = pow(1.0 + mach*mach*0.2, expon) - 1.0;
	rval = spress*t;
	T_(T_(trc.dprint("returning ",rval);))
	return rval;
}

#ifdef MAIN
#undef MAIN

#include "main.h"

using namespace std;

int
main() {
	T_(Trace trc(1,"atmos");)

	const int nstep = 1001;
	int i, j;
	const double zmin{atmos::altmin()};
	const double zmax{atmos::altmax()};
	const double del = (zmax - zmin)/((double)(nstep-1));
	FILE *fp;
	struct _plot {
		char const* filename;
		char const* title;
		char const* var;
		char const* xtitle;
		Ad (*eval)(Ad const&);
		double conv;
		char const* uscs;
	} plot[] = {
		{"temp.apf", "Temperature", "temp", "Kelvin", atmos::temp, 1.0, "Kelvin" },
		{"press.apf", "Pressure", "p", "Pascal", atmos::press, 1.0, "psi"},
		{"rho.apf", "Density", "rho", "Kg/m^3", atmos::rho, 1.94032e-3, "slug/ft^3"},
		{"vsonic.apf", "Sonic Velocity", "sv", "m/s", atmos::vsonic, 3.281, "ft/s" },
		{0, 0, 0, 0, (Ad (*)(Ad const&))NULL, 0.0 }
	};

	for (i=0; plot[i].filename; i++) {
		fp = fopen(plot[i].filename, "w");
		assert (fp != (FILE*)NULL);
		fprintf (fp, "$ %s\n", plot[i].title);
		fprintf (fp, "+ 1 Altitude (ft)\n");
		fprintf (fp, "+ 2 %s\n", plot[i].xtitle);
		fprintf (fp, "+ 3 %s %s\n", plot[i].xtitle,plot[i].uscs);
		fprintf (fp, "curve\n");
		fprintf (fp, "alt %s uscs\n", plot[i].var);
		for (j=0; j<nstep; j++) {
			Ad z(zmin + (double)j*del);
			Ad t = (*plot[i].eval)(z);
			double tuscs = t.value()*plot[i].conv;
			fprintf (fp, "%g %g %g\n", z.value(), t.value(), tuscs);
		}
		fprintf(fp, "*EOF\n");
		fclose (fp);
		cout << "Created \"" << plot[i].filename << "\" with parameters: alt, "
				<< plot[i].var << ", uscs\n";
	}
}
#endif
