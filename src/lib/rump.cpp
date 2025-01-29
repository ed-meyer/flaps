
#include <iostream>
#include <cmath>

using namespace std;

int
main(int argc, char** argv) {
	long double c1 = 333.75;
	long double c2 = 11.0;
	long double c3 = 121.0;
	long double c4 = 5.5;
	long double a{77617.0};
	long double b{33096.0};
	// long double f = c1*b*b*b*b*b*b  + a*a*(c2*a*a*b*b - b*b*b*b*b*b - c3*b*b*b*b - 2) + c4*pow(b,8.0) + a/(2.0*b);
	long double f = c1*pow(b, 6.0) + a*a*(c2*a*a*b*b - pow(b,6.0) - c3*pow(b, 4.0) - 2.0) + c4*pow(b,8.0) + a/(2.0*b);
	cout << "long double f = " << f << endl;
}
