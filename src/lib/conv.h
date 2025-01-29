//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// conversion factors, mostly from SI to USCS: multiply a number
// in equation units (SI) to get presentation units (USCS)
// All are in namespace flaps
//   m2ft            ft/meter
//   kg2slug         slug/kg
//   kg2lbm          lbm/kg
//   radps2Hz        Hz/(radian/second)
//   radps2rpm       (rev/min)/(radian/second)
//   rad2deg         degree/radian 57.295779
//   mps2knot        knot/(meter/second)
//   N2lbf           lbf/Newton
//   Pa2psf          (lbf/ft^2)/Pascal
//   kgpm32slugpft3  (slug/ft^3)/(kg/m^3)
//   kgpm32lbmpin3   (lbm/in^3)/(kg/m^3)

#ifndef conv_h
#define conv_h

// #ifdef HAVE_NUMBERS
#if __cplusplus > 201703L
#include <numbers>   // requires -std=c++20
namespace flaps {
constexpr double pi{std::numbers::pi};
}
#else
namespace flaps {
constexpr double pi{3.141592653589793238462643383279502884L};
}
#endif

namespace flaps {
// Length: meters to ft
constexpr double m2ft{1.0/(12.0*0.0254)};
// Mass: kg to slug
constexpr double kg2slug{1.0/14.593904};
// Mass: kg to lbm
constexpr double kg2lbm{1.0/0.45359237};
// Frequency: rad/s to Hz  0.15915494309
constexpr double radps2Hz{1.0/(2.0*flaps::pi)};
// Frequency: rad/s to Hz  0.15915494309
constexpr double radps2rpm{60.0/(2.0*flaps::pi)};
// Angle: rad to degrees
constexpr double rad2deg{180.0/flaps::pi};
// Velocity: m/s to knot
constexpr double mps2knot{3600.0/1852.0};
// Force: Newton to lb_f
constexpr double N2lbf{1.0/4.4482216152605};
// Pressure: Pa (N/m^2) to lb_f/ft^2
constexpr double Pa2psf{N2lbf/(m2ft*m2ft)};
// Density: kg/m^3 to slug/ft^3
constexpr double kgpm32slugpft3{0.3048*0.3048*0.3048*kg2slug};
// Density: kg/m^3 to lbm/in^3
constexpr double kgpm32lbmpin3{0.0254*0.0254*0.0254*kg2lbm};
} // namespace flaps

#endif // conv_h
