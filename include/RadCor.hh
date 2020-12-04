#ifndef RADCOR_HH
#define RADCOR_HH 

#include <cstdlib>
#include <iostream>
#include <iomanip> 
#include <cmath> 
#include <vector> 
#include "Target.h"
#include "Spectrum.h"
#include "Interpolation.h"
#include "Parameters.h"
#include "Kinematics.h"

using namespace::std;
using std::vector; 

class RadCor; 

class RadCor{

        // see Stein et. al., Appendix A [Phys. Rev. D 12, 1884 (1975)]
        // for formulas FTilde, xi, eta, phi(v), EsMin(Ep), EpMax, R(Es,Ep),Tr(Q2), CFACT(Es,Ep) 
        // EsIntegrand(Es'), EpIntegrand(Ep'), etc.  
 
	private: 
		int fUnits;                                 // Units of the XS data
		int fConstOfInterp;                         // Constant of interpolation (0 = y, 1 = yScale, 2 = x) [from input file] 
		bool fIsElastic,fIsQuasiElastic;            // boolean quantities to turn on elastic tail and quasi-elastic (and inelastic) tail subtraction 
                double fb,fXi,fEta,fFTilde,fR,fTr;          // Constants used in the evaluation of the tails
                double fTa,fTb,fT,fThRad,fThDeg,fA,fZ;      // Values obtained from the Target class
                double fDeltaE;                             // Parameter introduced to avoid infrared divergence [MeV] 
                double fAlpha;                              // Fine structure constant [-] 
                double fEs,fEp,fQ2;                         // (Es,Ep,Q2) from data 
                double fMe,fMp,fMT;                         // Electron mass (MeV), proton mass (MeV), target mass (MeV)  
                double fAnsEs,fAnsEp;                       // The integrated result of the Es and Ep integrals 
                double fCFACT;                              // Multiplicative factor in the radiative tail calculation 
                double PI; 
                double fMott;
		double fRCSystErrEst;                       // Estimate of the systematic error from radiative corrections [from input file] 
                vector<double> fNormFactor;                 // Normalization factor [from input file] 
                vector<Spectrum *> fS,fSNew;                // Vector of cross section spectra: fSNew holds the unfolded results
                                                            // for a given iteration  
		vector<double> fXSNew,fXSStatErrNew;
                vector<double> fXSSystErrNew;
                Interpolation *fInterp; 
                Kinematics *fK; 

        public: 
                RadCor();
                virtual ~RadCor(); 

                void Init();
                void Print(); 
		void Clear();

                void Run(int);
                void Iterate(int); 

                void SetA(double A){fA = A;} 
                void SetZ(double Z){fZ = Z;} 
                void SetThDeg(double th){fThDeg = th;} 
                void SetThRad(double thr){fThRad = thr;} 
                void SetTa(double ta){fTa = ta;}
                void SetTb(double tb){fTb = tb;}

                void SetData(vector<Spectrum *>); 
                void SetTargetParameters(Target *); 
                void SetParameters(Parameters *); 
                void CalculateVariables(int,int);           // Set up to calculate the radiative corrections for a given Es and Ep [vector]               
                void Unfold(int,int);                       // Subtract off Es and Ep integrals from XSf (Born XS model)       
                void Fold(int,int);                         // Add on the Es and Ep integrals to a Born XS        
		void CalculateError(int,int);               // Calculate the stat. and syst. error from input values and the RC's 
		void UpdateModel();                         // Update the fS vector with the values saved in the fSNew vector 
		void SaveNewModel(int);                     // Save the newly calculated unfolded XS data to a Spectrum object and 
                                                            // push_back the fSNew vector

                void CalculateB(); 
                void CalculateXi();
		void CalculateEta(); 
		void CalculateR();                          // R(Es,Ep) 
                void CalculateCFACT();                      // CFACT(Es,Ep) 
                void CalculateEsIntegral();
                void CalculateEpIntegral();

                double EsIntegrand(const double);
                double EpIntegrand(const double); 
                double Integrate(double (RadCor::*)(const double),double,double,double,int);
                double AdaptiveSimpsonAux(double (RadCor::*)(const double),double,double,double,double,double,double,double,int);

		double GetTr(double);                        // Tr(Q2): changes PER INTEGRATION BIN 
		double GetFTilde(double);                    // FTilde(Q2): changes PER INTEGRATION BIN 
		double GetPhi(double);                       // phi(v), v = arbitrary value  
		double GetEsMin(double);                     // EsMin(Ep) 
		double GetEpMax(double);                     // EpMax(Es)                        
		double GetQ2(double,double);                 // Q2(Es,Ep) 
		double GetSpence(double);                    // Spence(x), x = arbitrary value  

                vector<Spectrum *> GetSpectra(){return fS;}  // Return vector of cross section spectra after RC's 
                 

};

#endif  
