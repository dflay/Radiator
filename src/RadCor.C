#include "../include/RadCor.h"
//_____________________________________________________________________________________________
RadCor::RadCor(){
	Init();
	Clear();
}
//_____________________________________________________________________________________________
RadCor::~RadCor(){
	Clear();
        delete fInterp;
        delete fK; 
}
//_____________________________________________________________________________________________
void RadCor::Init(){

	fIsElastic=false;
	fIsQuasiElastic=false;
        PI=4.*atan(1); 
	fAlpha=1.0/137.03604; 
	fDeltaE=0.0;            // in MeV
	fMe=0.511;              // in MeV
        fMp=938.0;              // in MeV
        // fMp=931.494;            // in MeV
        fA=0;
        fZ=0;
	fb=0;
	fXi=0;
	fEta=0;
	fTa=0;
	fTb=0;
	fT=0;
        fTr=0;
        fThDeg=0; 
	fThRad=0;
	fEs=0;
	fEp=0;
        fQ2=0; 
        fR=0; 
        fFTilde=0; 
        fAnsEs=0;
        fAnsEp=0; 
        fCFACT=0;
	fMT=0;
	fUnits=-1;
        fConstOfInterp=0;
	fRCSystErrEst=0;
        fInterp = new Interpolation(); 
        fK = new Kinematics();
}
//_____________________________________________________________________________________________
void RadCor::Clear(){
	fS.clear();
	fSNew.clear();
	fXSNew.clear();
	fXSStatErrNew.clear();
	fXSSystErrNew.clear();
	fNormFactor.clear();
}
//_____________________________________________________________________________________________
void RadCor::SetData(vector<Spectrum *> S){

	const int N = S.size();
	for(int i=0;i<N;i++){
		fS.push_back(S[i]); 
	}
        fInterp->SetData(fS);
	fInterp->SetA(fA);
        fInterp->SetScale(fUnits);
        fInterp->SetConstOfInterp(fConstOfInterp);  

} 
//_____________________________________________________________________________________________
void RadCor::SetTargetParameters(Target *T){

	fTa    = T->GetTa();
	fTb    = T->GetTb();
	fA     = T->GetA();
	fZ     = T->GetZ(); 
	fThDeg = T->GetThDeg(); 
        fThRad = T->GetThRad(); 
	fT     = fTa + fTb;
        fMT    = fMp*fA;        // Set the target mass  

        CalculateEta(); 
        CalculateXi(); 
        CalculateB(); 

}
//_____________________________________________________________________________________________
void RadCor::SetParameters(Parameters *P){

        fConstOfInterp = P->GetConstOfInterp(); 
	fRCSystErrEst  = P->GetRCSystErrEst(); 
        fUnits         = P->GetUnits(); 
        fDeltaE        = P->GetDeltaE(); 

        const int N = fS.size(); 
        for(int i=0;i<N;i++){
		fNormFactor.push_back( P->GetNormFactor(i) ); 
        }
        
}
//_____________________________________________________________________________________________
void RadCor::Run(int NumOfIter){

	for(int i=0;i<NumOfIter;i++){
		Iterate(i); 
	}

}
//_____________________________________________________________________________________________
void RadCor::Iterate(int Iter){

	int M=0;
        const int N = fS.size(); 
	cout << "----- Iteration " << Iter+1 << " -----" << endl;
	for(int j=0;j<N;j++){                                 // 2nd-level: loop over Es 
		M = fS[j]->GetSize();
		cout << "Es = " << fixed << setprecision(1)
		     << fS[j]->GetEs() << endl;
		for(int k=0;k<M;k++){                         // 3rd level: loop over Ep 
			// cout <<"bin = " << k << endl;
			CalculateVariables(j,k);              // calculate important quantities
			CalculateEsIntegral();
			CalculateEpIntegral();
			Unfold(j,k);                          // Subtract off Es and Ep integrals from newest unfolded model 
		}
		SaveNewModel(j);                              // For each Es, we save the unfolded XS result to a new vector of spectra
	}
	UpdateModel();

}
//_____________________________________________________________________________________________
void RadCor::CalculateVariables(int i,int j){

        // For a given Ep bin, we grab the appropriate value of Es and Ep 
	fEs = fS[i]->GetEs(); 
	fEp = fS[i]->GetEp(j); 
	fQ2 = fK->GetQ2(fEs,fEp,fThDeg); 

        CalculateR(); 
	CalculateCFACT();

}
//_____________________________________________________________________________________________
void RadCor::CalculateEta(){

	double Z23   = pow(fZ,-2.0/3.0);
	double Z13   = pow(fZ,-1.0/3.0);
	double num   = log(1440.0*Z23); 
	double denom = log(183.0*Z13); 
	fEta        = num/denom; 

}
//_____________________________________________________________________________________________
void RadCor::CalculateB(){

	double Z13 = pow(fZ,-1.0/3.0);
        double T1  = 1.0;
        double T2  = (1.0/9.0)*( (fZ+1.0)/(fZ+fEta) ); 
        double T3  = 1.0/log(183.0*Z13); 
        fb         = (4.0/3.0)*(T1 + T2*T3); 

}
//_____________________________________________________________________________________________
void RadCor::CalculateXi(){

	double Z13 = pow(fZ,-1.0/3.0);
        double T1  = PI*fMe/(2.0*fAlpha); 
        double T2  = fT/( (fZ+fEta)*log(183.0*Z13) );
        fXi       = T1*T2; 

}
//_____________________________________________________________________________________________
void RadCor::CalculateR(){

	double SIN   = sin(fThRad/2.0); 
	double SIN2  = SIN*SIN; 
	double num   = fMT + 2.0*fEs*SIN2;
	double denom = fMT - 2.0*fEp*SIN2; 
        fR          = num/denom; 

}
//_____________________________________________________________________________________________
void RadCor::CalculateCFACT(){

        // General terms 
        double Q2     = fK->GetQ2(fEs,fEp,fThDeg); 
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
void RadCor::CalculateEsIntegral(){

	int depth      = 10; 
        double epsilon = 1e-10; 
	double min     = GetEsMin(fEp); 
	double max     = fEs - fR*fDeltaE; 
	fAnsEs         = Integrate(&RadCor::EsIntegrand,min,max,epsilon,depth);
}
//_____________________________________________________________________________________________
void RadCor::CalculateEpIntegral(){

	int depth      = 10; 
        double epsilon = 1e-10; 
	double min     = fEp + fDeltaE; 
	double max     = GetEpMax(fEs); 
	fAnsEp         = Integrate(&RadCor::EpIntegrand,min,max,epsilon,depth);
}
//____________________________________________________________________________________________
void RadCor::Unfold(int i,int j){

        // NOTE: XSi is NOT scaled by the MottXS!  We must do it here
        //       so that the final answer makes sense. 
        double MottXS = fS[i]->GetMottXS(); 
	double XS     = fS[i]->GetXSi(j)/MottXS;  
	double arg    = ( XS - (fAnsEs + fAnsEp)/MottXS )/fCFACT; 

        // store result to a vector.  This will be propagated to
        // the vector of spectra at the end of the iteration loop
        fXSNew.push_back(arg); 

        CalculateError(i,j); 
}
//____________________________________________________________________________________________
void RadCor::Fold(int i,int j){

        // FIXME: Here, we tack on the radiation 

        // NOTE: XSi is NOT scaled by the MottXS!  We must do it here
        //       so that the final answer makes sense. 
        double MottXS = fS[i]->GetMottXS(); 
	double XS     = fS[i]->GetXSi(j)/MottXS;  
	double arg    = ( XS - (fAnsEs + fAnsEp)/MottXS )/fCFACT; 

        // store result to a vector.  This will be propagated to
        // the vector of spectra at the end of the iteration loop
        fXSNew.push_back(arg); 

        CalculateError(i,j); 
}
//____________________________________________________________________________________________
void RadCor::CalculateError(int i,int j){

       // FIXME: This calculation is still WRONG! The errors are way too big!  

       double MottXS     = fS[i]->GetMottXS(); 
       double XSf        = fXSNew[j]; 
       double XSi        = fS[i]->GetXSi(j)/MottXS;
       double XSiStatErr = fS[i]->GetXSiStatErr(j)/MottXS; 
       double XSiSystErr = fS[i]->GetXSiSystErr(j)/MottXS; 
       double RAT        = XSf/XSi; 

       // initial stat. and syst. error
       // NOTE: these are in units/MeV/sr.  We divide by XSi to convert to %   
       double iStat = XSiStatErr/XSi; 
       double iSyst = XSiSystErr/XSi;

       // final stat. error (we multiply by XSi to convert it back into units/MeV/sr) 
       double FinStat = XSi*RAT*iStat; 
       // final syst. error from radiative corrections 
       double RCSyst  = abs(XSf-XSi)*fRCSystErrEst; 
       RCSyst        *= (1.0/XSf); 
       // final syst. error = in-quadrature sum of intial and RC syst. errors 
       // multiplying by XSf converts it back into units/MeV/sr.  Remember, 
       // XSf is unitless... so we have to multiply by the Mott XS. 
       double A       = iSyst*RAT;
       double B       = RCSyst; 
       double A2      = A*A;
       double B2      = B*B; 
       double FinSyst = XSf*sqrt(A2 + B2); 
 
       fXSStatErrNew.push_back(FinStat); 
       fXSSystErrNew.push_back(FinSyst); 

}
//____________________________________________________________________________________________
void RadCor::SaveNewModel(int i){

        // NOTE: index i is for the ith spectra  
 
	// cout << "----------- Saving New Model ----------- " << endl;
	Spectrum *uS = new Spectrum(); 
	uS->SetEs( fS[i]->GetEs() ); 
	uS->SetThDeg( fS[i]->GetThDeg() ); 
	const int N = fS[i]->GetSize();  // Need size of Ep vector 
	for(int j=0;j<N;j++){
		uS->SetEp( fS[i]->GetEp(j) );
		uS->SetXSi(0.0); 
		uS->SetXSiStatErr(0.0); 
		uS->SetXSiSystErr(0.0);  
		uS->SetXSf( fXSNew[j] ); 
		uS->SetXSfStatErr( fXSStatErrNew[j] ); 
		uS->SetXSfSystErr( fXSSystErrNew[j] ); 
                // cout << fixed      << fS[i]->GetEp(j)   << "\t XS(i) = " 
                //      << scientific << fS[i]->GetXSf(j)  << "\t XS(f) = " << fXSNew[j] << endl;
	}
	// uS->Print();
        fSNew.push_back(uS); 

        // set up for next spectrum: 
        fXSNew.clear();
        fXSStatErrNew.clear();
        fXSSystErrNew.clear();

}
//____________________________________________________________________________________________
void RadCor::UpdateModel(){

	// Update the XS spectra at the conclusion of iteration i. 
        int M=0; 
        double NewXS=0,NewStat=0,NewSyst=0; 

        // Clear XSf, XSfErr vectors for all spectra: 
        const int N = fSNew.size(); 

        for(int j=0;j<N;j++){
		fS[j]->ClearXSf(); 
		fS[j]->ClearXSfErr(); 
	}

        // Fill XSf, XSfErr vectors with updated values from the SNew vector of spectra:   
	for(int j=0;j<N;j++){                                       // Es loop
		M = fSNew[j]->GetSize();
		for(int k=0;k<M;k++){                               // Ep loop 
			NewXS   = fSNew[j]->GetXSf(k);
			NewStat = fSNew[j]->GetXSfStatErr(k);
			NewSyst = fSNew[j]->GetXSfSystErr(k);
			fS[j]->SetXSf(NewXS);
			fS[j]->SetXSfStatErr(NewStat);
			fS[j]->SetXSfSystErr(NewSyst);
		}
	}

        // Update the Interpolation object: 
        fInterp->SetData(fS); 

        // now that the model is updated, we DELETE the Spectrum objects from fSNew
        // to set up for the next iteration 
        for(int j=0;j<N;j++){
		delete fSNew[j]; 
        }

        fSNew.clear();


}
//_____________________________________________________________________________________________
double RadCor::GetPhi(double v){

	double phi = 1.0 - v + (3.0/4.0)*v*v; 
	return phi; 

}
//_____________________________________________________________________________________________
double RadCor::GetTr(double Q2){

        // General terms
	double M2 = fMe*fMe; 
        // Individual terms
	double T1 = (1.0/fb)*(fAlpha/PI); 
	double T2 = log(Q2/M2) - 1.0; 
        // Put it all together 
	double Tr = T1*T2;
        return Tr; 

}
//_____________________________________________________________________________________________
double RadCor::GetFTilde(double Q2){

        // General terms
	double M2     = fMe*fMe;
	double PI2    = PI*PI;  
        double COS    = cos(fThRad/2.0);
        double COS2   = COS*COS; 
        double SPENCE = GetSpence(COS2); 
        // Individual terms 
	double T1     = 1.0 + 0.5772*fb*fT; 
	double T2     = (2.0*fAlpha/PI)*( (-14.0/9.0) + (13.0/12.0)*log(Q2/M2) );
	double T3     = (-1.0)*(fAlpha/(2.0*PI))*log( pow(fEs/fEp,2.0) ); 
	double T4     = (fAlpha/PI)*( (PI2/6.0) - SPENCE );
        // Put it all together
        double FTilde = T1+T2+T3+T4;  
        return FTilde; 

}
//_____________________________________________________________________________________________
double RadCor::GetEsMin(double Ep){

        double num   = Ep;
        double SIN   = sin(fThRad/2.0); 
        double SIN2  = SIN*SIN;  
        double denom = 1.0 - (2.0*Ep/fMT)*SIN2; 
        double EsMin = num/denom; 
        return EsMin; 

}
//_____________________________________________________________________________________________
double RadCor::GetEpMax(double Es){

        double num   = Es;
        double SIN   = sin(fThRad/2.0); 
        double SIN2  = SIN*SIN;  
        double denom = 1.0 + (2.0*Es/fMT)*SIN2; 
        double EpMax = num/denom; 
        return EpMax; 

}
//_____________________________________________________________________________________________
double RadCor::GetSpence(double x){

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
double RadCor::EsIntegrand(const double EsPrime){

        // general terms
        double Q2       = fK->GetQ2(EsPrime,fEp,fThDeg); 
        double FTilde   = GetFTilde(Q2);  
        double Tr       = GetTr(Q2); 
	double dEs      = fEs-EsPrime; 
	double v        = dEs/fEs; 
	double phi      = GetPhi(v); 
	double Sig      = fInterp->FTCS(EsPrime,fEp); 
	// double Sig      = fCS->GetXS(EsPrime,fEp,fThDeg);  
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
double RadCor::EpIntegrand(const double EpPrime){

        // general terms 
        double Q2       = fK->GetQ2(fEs,EpPrime,fThDeg); 
        double Tr       = GetTr(Q2);
        double FTilde   = GetFTilde(Q2); 
	double dEp      = EpPrime-fEp; 
	double v        = dEp/EpPrime; 
	double phi      = GetPhi(v);  
	double Sig      = fInterp->FTCS(fEs,EpPrime);  
	// double Sig      = fCS->GetXS(fEs,EpPrime,fThDeg);  
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
double RadCor::Integrate(double (RadCor::*f)(const double),double A,double B,double epsilon,int Depth){
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
double RadCor::AdaptiveSimpsonAux(double (RadCor::*f)(const double),double A,double B,double epsilon,
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
void RadCor::Print(){

	cout << "------------------------------------"              << endl;
	cout << "Radiative correction quantities: "                 << endl;
        cout << "RC systematic error estimate = " 
             << fixed << setprecision(2) << fRCSystErrEst  << endl;
        cout << "DeltaE = " << fixed      << setprecision(1) << fDeltaE << " [MeV]"  << endl;
        cout << "Constants for given thicknesses: "                 << endl; 
	cout << "Tb     = " << scientific << setprecision(4) << fTb     << " [#X0]"  << endl;
	cout << "Ta     = " << scientific << setprecision(4) << fTa     << " [#X0]"  << endl;
	cout << "eta    = " << fixed      << setprecision(4) << fEta    << " [-]"    << endl;
	cout << "b      = " << fixed      << setprecision(4) << fb      << " [-]"    << endl;
	cout << "xi     = " << fixed      << setprecision(4) << fXi     << " [-]"    << endl;
        cout << "Values that change for each (Es,Ep): "   << endl; 
        cout << "Es         = " << fixed      << setprecision(1) << fEs     << " [MeV]"   << endl;  
        cout << "Ep         = " << fixed      << setprecision(1) << fEp     << " [MeV]"   << endl;  
        cout << "Q2         = " << scientific << setprecision(4) << fQ2     << " [MeV^2]" << endl;  
        cout << "R          = " << fixed      << setprecision(4) << fR      << " [-]"     << endl;
        cout << "CFACT      = " << scientific << setprecision(4) << fCFACT  << " [-]"     << endl;
        cout << "EsInt      = " << scientific << setprecision(4) << fAnsEs  << " [-]"     << endl;
        cout << "EpInt      = " << scientific << setprecision(4) << fAnsEp  << " [-]"     << endl;


}
