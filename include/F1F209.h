// F1F209 class
// Notes: - Variables: beam energy (Es) and scattered energy (Ep) are in GeV, 
//          scattering angle (th) is in degrees.  A is in g/mol. 
//        - Input to the GetXS method: Z, A, Es, Ep, th 
//        - The resulting inelastic cross section, calculated from the fortran subroutine 
//          F1F209, is in units of pb/MeV/sr. 

#ifndef F1F209_H
#define F1F209_H

#include <cstdlib> 
#include <iostream>
#include <iomanip> 
#include <cmath>
#include "eInclusiveCrossSection.h"

extern "C"{
   void f1f209(double *,double *,double *,double *,double *,double *,double *,double *);
}

using namespace std; 

class F1F209: public eInclusiveCrossSection{

	private: 

		double GetQ2(double,double,double);
		double GetW(double,double,double);

	public: 
		F1F209();
		~F1F209();

		// concrete implementation of GetBornXS()
		double GetBornXS();   

};

#endif
