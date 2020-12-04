#include "../include/Radiator.h"
//________________________________________________________________________
Radiator::Radiator(){
	Init();
}
//________________________________________________________________________
Radiator::~Radiator(){

}
//________________________________________________________________________
void Radiator::Init(){
        fDeltaE   = 0.01;           // in GeV
        fMT       = 0; 
        fZ        = 0; 
        fA        = 0; 
        fA        = 0;
        fZ        = 0;
	fb        = 0;
	fXi       = 0;
	fEta      = 0;
	fTa       = 0;
	fTb       = 0;
	fT        = 0;
        fThDeg    = 0; 
	fEs       = 0;
	fEp       = 0;
        fR        = 0; 
        fCFACT    = 0;
	fMT       = 0;
}
//_____________________________________________________________________________________________
double Radiator::Radiate(){

        // set important variables 
	fZ     = fInclXS->GetZ();
	fA     = fInclXS->GetA();
        fEs    = fInclXS->GetEs(); 
        fEp    = fInclXS->GetEp();
        fThDeg = fInclXS->GetTh(); 
        fMT    = fA*PROTON_MASS;      // set the target mass 
        fT     = fTa + fTb; 

	if( (fTa==0)||(fTb==0) ){
		cout << "[GetRadiatedXS]: Radiation lengths are zero! Check your input... " << endl;
		exit(1);
        }

	CalculateEta();
	CalculateB();
	CalculateXi();
	CalculateR();
        CalculateCFACT(); 

        double BornXS = fInclXS->GetBornXS(); 
        double AnsEs  = CalculateEsIntegral(); 
        double AnsEp  = CalculateEpIntegral(); 
        double RadXS  = fCFACT*BornXS + AnsEs + AnsEp;  

        return RadXS;

}
//________________________________________________________________________
double Radiator::GetQ2(double Es,double Ep,double th){
	double thr  = th*DEG_TO_RAD; 
	double SIN  = sin(thr/2.); 
	double SIN2 = SIN*SIN;
        double Q2   = 4.*Es*Ep*SIN2; 
        return Q2;  
}
//________________________________________________________________________
double Radiator::GetW(double Es,double Ep,double th){
	double Nu = Es-Ep;
	double Q2 = GetQ2(Es,Ep,th); 
        double W2 = PROTON_MASS*PROTON_MASS + 2.*PROTON_MASS*Nu - Q2;
	return sqrt(W2);
}
//_____________________________________________________________________________________________
double Radiator::GetPhi(double v){

	double phi = 1.0 - v + (3.0/4.0)*v*v; 
	return phi; 

}
//_____________________________________________________________________________________________
double Radiator::GetTr(double Q2){

        // General terms
	double M2 = ELECTRON_MASS*ELECTRON_MASS; 
        // Individual terms
	double T1 = (1.0/fb)*(ALPHA/PI); 
	double T2 = log(Q2/M2) - 1.0; 
        // Put it all together 
	double Tr = T1*T2;
        return Tr; 

}
//_____________________________________________________________________________________________
double Radiator::GetFTilde(double Q2){

        // General terms
	double M2     = ELECTRON_MASS*ELECTRON_MASS;
	double PI2    = PI*PI; 
        double thr    = fThDeg*DEG_TO_RAD;  
        double COS    = cos(thr/2.0);
        double COS2   = COS*COS; 
        double SPENCE = GetSpence(COS2); 
        // Individual terms 
	double T1     = 1.0 + 0.5772*fb*fT; 
	double T2     = (2.0*ALPHA/PI)*( (-14.0/9.0) + (13.0/12.0)*log(Q2/M2) );
	double T3     = (-1.0)*(ALPHA/(2.0*PI))*log( pow(fEs/fEp,2.0) ); 
	double T4     = (ALPHA/PI)*( (PI2/6.0) - SPENCE );
        // Put it all together
        double FTilde = T1+T2+T3+T4;  
        return FTilde; 

}
//_____________________________________________________________________________________________
double Radiator::GetEsMin(double Ep){

        // elastic 
        // double num   = Ep;
        // double thr   = fThDeg*DEG_TO_RAD;  
        // double SIN   = sin(thr/2.0); 
        // double SIN2  = SIN*SIN;  
        // double denom = 1.0 - (2.0*Ep/fMT)*SIN2; 

        // pion production threshold 
        double num   = PION_MASS*PION_MASS + 2.*PROTON_MASS*PION_MASS + 2.*PROTON_MASS*Ep;
        double thr   = fThDeg*DEG_TO_RAD;  
        double SIN   = sin(thr/2.0); 
        double SIN2  = SIN*SIN;  
        double denom = 2.*PROTON_MASS - (4.0*Ep)*SIN2; 

        double EsMin = num/denom; 
        return EsMin; 

}
//_____________________________________________________________________________________________
double Radiator::GetEpMax(double Es){

        // elastic 
        // double num   = Es;
        // double thr   = fThDeg*DEG_TO_RAD;  
        // double SIN   = sin(thr/2.0); 
        // double SIN2  = SIN*SIN;  
        // double denom = 1.0 + (2.0*Es/fMT)*SIN2; 

        // pion production threshold 
        double num   = 2.*PROTON_MASS*Es - 2.*PROTON_MASS*PION_MASS - PION_MASS*PION_MASS;
        double thr   = fThDeg*DEG_TO_RAD;  
        double SIN   = sin(thr/2.0); 
        double SIN2  = SIN*SIN;  
        double denom = 2.*PROTON_MASS + (4.0*Es)*SIN2; 

        double EpMax = num/denom; 
        return EpMax; 

}
//_____________________________________________________________________________________________
double Radiator::GetSpence(double x){

        // converted from radcor.f: 
        double num=0,denom=0,Index=0; 
	double PI2 = PI*PI; 
	double ans = (PI2/6.0) - log(x)*log(1.-x); 

        for(int i=0;i<50;i++){
		Index   = (double)i + 1.0; 
		num     = pow(x,i+1.0);
		denom   = pow(Index,2.);
		if(denom>0){
			ans -= num/denom; 
		}else{
			ans -= 0;
		}
        }

        return ans; 

}
//_____________________________________________________________________________________________
void Radiator::CalculateEta(){

	double Z23   = pow(fZ,-2.0/3.0);
	double Z13   = pow(fZ,-1.0/3.0);
	double num   = log(1440.0*Z23); 
	double denom = log(183.0*Z13); 
	fEta         = num/denom; 

}
//_____________________________________________________________________________________________
void Radiator::CalculateB(){

	double Z13 = pow(fZ,-1.0/3.0);
        double T1  = 1.0;
        double T2  = (1.0/9.0)*( (fZ+1.0)/(fZ+fEta) ); 
        double T3  = 1.0/log(183.0*Z13); 
        fb         = (4.0/3.0)*(T1 + T2*T3); 

}
//_____________________________________________________________________________________________
void Radiator::CalculateXi(){

	double Z13 = pow(fZ,-1.0/3.0);
        double T1  = PI*ELECTRON_MASS/(2.0*ALPHA); 
        double T2  = fT/( (fZ+fEta)*log(183.0*Z13) );
        fXi        = T1*T2; 

}
//_____________________________________________________________________________________________
void Radiator::CalculateR(){
 
        double thr   = fThDeg*DEG_TO_RAD; 
	double SIN   = sin(thr/2.0); 
	double SIN2  = SIN*SIN; 
	double num   = fMT + 2.0*fEs*SIN2;
	double denom = fMT - 2.0*fEp*SIN2; 
        fR           = num/denom; 

}
//_____________________________________________________________________________________________
void Radiator::CalculateCFACT(){

        // General terms 
        double Q2     = GetQ2(fEs,fEp,fThDeg); 
        double Tr     = GetTr(Q2); 
        double FTilde = GetFTilde(Q2); 
        // First term
	double Term1  = fR*fDeltaE/fEs;
	double Exp1   = fb*(fTb+Tr);  
	double T1     = pow(Term1,Exp1); 
        // Second term 
	double Term2  = fDeltaE/fEp; 
	double Exp2   = fb*(fTa+Tr); 
	double T2     = pow(Term2,Exp2);    
        // Third term
        double num    = fXi/fDeltaE; 
        double denom  = 1.0 - fb*(fTa+fTb+2.0*Tr); 
        double T3     = 1.0 - num/denom; 
        // Put it all together 
        fCFACT       = FTilde*T1*T2*T3; 

}
//_____________________________________________________________________________________________
double Radiator::EsIntegrand(const double EsPrime){

        // general terms
        double Q2       = GetQ2(EsPrime,fEp,fThDeg); 
        double FTilde   = GetFTilde(Q2);  
        double Tr       = GetTr(Q2); 
	double dEs      = fEs-EsPrime; 
	double v        = dEs/fEs; 
	double phi      = GetPhi(v); 
	fInclXS->SetEs(EsPrime);
	fInclXS->SetEp(fEp);
	fInclXS->SetTh(fThDeg);
	double Sig      = fInclXS->GetBornXS();  
	if(Sig!=Sig){
		cout << "[Radiator::EsIntegrand]: Invalid cross section! " << endl;
		exit(1);
        }
        double SigTilde = FTilde*Sig;  
	// first term 
	double Term1    = dEs/(fEp*fR);
	double Exp1     = fb*(fTa+Tr); 
	double T1       = pow(Term1,Exp1); 
	// second term
	double Term2    = dEs/fEs; 
	double Exp2     = fb*(fTb+Tr); 
	double T2       = pow(Term2,Exp2); 
	// third term 
	double T3       = fb*( ((fTb+Tr)/dEs)*phi + fXi/(2.0*pow(dEs,2.0)) ); 
	// put it all together
	double FES      = T1*T2*T3*SigTilde; 

	return FES; 

}
//_____________________________________________________________________________________________
double Radiator::EpIntegrand(const double EpPrime){

        // general terms 
        double Q2       = GetQ2(fEs,EpPrime,fThDeg); 
        double Tr       = GetTr(Q2);
        double FTilde   = GetFTilde(Q2); 
	double dEp      = EpPrime-fEp; 
	double v        = dEp/EpPrime; 
	double phi      = GetPhi(v);  
	fInclXS->SetEs(fEs);
	fInclXS->SetEp(EpPrime);
	fInclXS->SetTh(fThDeg);
	double Sig      = fInclXS->GetBornXS();  
	if(Sig!=Sig){
		cout << "[Radiator::EpIntegrand]: Invalid cross section! " << endl;
		exit(1);
        }
        double SigTilde = FTilde*Sig;  
	// first term 
	double Term1    = dEp/(EpPrime);
	double Exp1     = fb*(fTa+Tr); 
	double T1       = pow(Term1,Exp1); 
	// second term
	double Term2    = (dEp*fR)/fEs; 
	double Exp2     = fb*(fTb+Tr); 
	double T2       = pow(Term2,Exp2); 
	// third term 
	double T3       = fb*( ((fTa+Tr)/dEp)*phi + fXi/(2.0*pow(dEp,2.0)) ); 
	// put it all together
	double FEP      = T1*T2*T3*SigTilde; 

	return FEP; 

}
//_____________________________________________________________________________________________
double Radiator::CalculateEsIntegral(){

	int depth      = 10; 
        double epsilon = 1e-10; 
	double min     = GetEsMin(fEp); 
	double max     = fEs - fR*fDeltaE; 
	double AnsEs   = Integrate(&Radiator::EsIntegrand,min,max,epsilon,depth);
	return AnsEs;
}
//_____________________________________________________________________________________________
double Radiator::CalculateEpIntegral(){

	int depth      = 10; 
        double epsilon = 1e-10; 
	double min     = fEp + fDeltaE; 
	double max     = GetEpMax(fEs); 
	double AnsEp   = Integrate(&Radiator::EpIntegrand,min,max,epsilon,depth);
	return AnsEp;
}
//_____________________________________________________________________________________________
double Radiator::Integrate(double (Radiator::*f)(const double),double A,double B,double epsilon,int Depth){
        // Adaptive Simpson's Rule
        double C   = (A + B)/2.0;
        double H   = B - A;
	double fa  = (this->*f)(A);
        double fb  = (this->*f)(B);
        double fc  = (this->*f)(C);
        double S   = (H/6.0)*(fa + 4.0*fc + fb);
        double ans = AdaptiveSimpsonAux(f,A,B,epsilon,S,fa,fb,fc,Depth);

        return ans; 
}
//_____________________________________________________________________________________________
double Radiator::AdaptiveSimpsonAux(double (Radiator::*f)(const double),double A,double B,double epsilon,
                                  double S,double fa,double fb,double fc,int bottom){
        // Recursive auxiliary function for AdaptiveSimpson() function
        double C      = (A + B)/2.0;
        double H      = B - A;
        double D      = (A + C)/2.0;
        double E      = (C + B)/2.0;
        double fd     = (this->*f)(D);
        double fe     = (this->*f)(E);
        double Sleft  = (H/12.0)*(fa + 4.0*fd + fc);
        double Sright = (H/12.0)*(fc + 4.0*fe + fb);
        double S2     = Sleft + Sright;
        if( (bottom <= 0) || (fabs(S2 - S) <= 15.0*epsilon) ){
                return S2 + (S2 - S)/15;
        }
        double arg = AdaptiveSimpsonAux(f,A,C,epsilon/2.0,Sleft, fa,fc,fd,bottom-1) +
                     AdaptiveSimpsonAux(f,C,B,epsilon/2.0,Sright,fc,fb,fe,bottom-1);
        return arg;
}
//_____________________________________________________________________________________________
void Radiator::Print(){

        cout << "------------------------------------"              << endl;
        cout << "Radiative correction quantities: "                 << endl;
        cout << "DeltaE = " << fixed      << setprecision(4) << fDeltaE << " [GeV]"  << endl;
        cout << "Constants for given thicknesses: "                 << endl;
        cout << "Tb     = " << scientific << setprecision(4) << fTb     << " [#X0]"  << endl;
        cout << "Ta     = " << scientific << setprecision(4) << fTa     << " [#X0]"  << endl;
        cout << "eta    = " << scientific << setprecision(4) << fEta    << " [-]"    << endl;
        cout << "b      = " << fixed      << setprecision(4) << fb      << " [-]"    << endl;
        cout << "xi     = " << scientific << setprecision(4) << fXi     << " [GeV]"  << endl;
        cout << "Values that change for each (Es,Ep): "   << endl;
        cout << "Es         = " << fixed      << setprecision(4) << fEs     << " [GeV]"   << endl;
        cout << "Ep         = " << fixed      << setprecision(4) << fEp     << " [GeV]"   << endl;
        cout << "R          = " << fixed      << setprecision(4) << fR      << " [-]"     << endl;
        cout << "CFACT      = " << scientific << setprecision(4) << fCFACT  << " [-]"     << endl;


}
