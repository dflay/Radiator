

#ifndef RADIATOR_H
#define RADIATOR_H

#include <cstdlib> 
#include <iostream>
#include <iomanip> 
#include <cmath>
#include "eInclusiveCrossSection.h"
#include "F1F209.h"

#define PI 3.14159265359
#define ALPHA 1./137. 
#define DEG_TO_RAD PI/180.
#define HBAR_C 624.4197  // in GeV*nb^(1/2) 

using namespace std; 

class F1F209; 

class Radiator{

	private: 
		double fZ,fA;
                // for GetRadiatedXS method, we'll need to 
                // make some private data members to make 
                // the integration go smoothly  
                double fDeltaE;
                double fMT;
                double fEs,fEp,fThDeg;
                double fR,fCFACT; 
                double fTa,fTb,fT,fEta,fXi,fb; 

                void CalculateB(); 
                void CalculateXi();
		void CalculateEta(); 
		void CalculateR();                          
                void CalculateCFACT();                       
                double CalculateEsIntegral();
                double CalculateEpIntegral();

		double GetQ2(double,double,double);
		double GetW(double,double,double);
	
                // For the GetRadiatedXS method
		double GetTr(double);                        // Tr(Q2): changes PER INTEGRATION BIN 
		double GetFTilde(double);                    // FTilde(Q2): changes PER INTEGRATION BIN 
		double GetPhi(double);                       // phi(v), v = arbitrary value  
		double GetEsMin(double);                     // EsMin(Ep) 
		double GetEpMax(double);                     // EpMax(Es)                        
		double GetQ2(double,double);                 // Q2(Es,Ep) 
		double GetSpence(double);                    // Spence(x), x = arbitrary value 
               
		double EsIntegrand(const double);
                double EpIntegrand(const double); 
                double Integrate(double (Radiator::*)(const double),double,double,double,int);
                double AdaptiveSimpsonAux(double (Radiator::*)(const double),double,double,double,double,double,double,double,int);

		eInclusiveCrossSection *fInclXS;

	public: 
		Radiator();
		~Radiator();

		void Init();
		void Print();

                void SetA(double A){fA = A;} 
                void SetZ(double Z){fZ = Z;} 
                void SetTh(double th){fThDeg = th;} 
                void SetTa(double ta){fTa = ta;}
                void SetTb(double tb){fTb = tb;}
        
		void SetCrossSection(eInclusiveCrossSection *XS){fInclXS = XS;}
 
                double Radiate(); 

};

#endif
