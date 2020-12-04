#include "../include/F1F209.hh"
//________________________________________________________________________
F1F209::F1F209(){

}
//________________________________________________________________________
F1F209::~F1F209(){

}
//________________________________________________________________________
double F1F209::GetBornXS(){
   double F1,F2; 
   double Q2 = GetQ2(fEs,fEp,fTh);
   double W  = GetW(fEs,fEp,fTh);
   double W2 = W*W; 

   f1f209(&fZ,&fA,&fEs,&fTh,&Q2,&W2,&F1,&F2);

   double thr    = fTh*DEG_TO_RAD;
   double TAN    = tan(thr/2.);
   double TAN2   = TAN*TAN;  
   double Nu     = fEs-fEp;
   double MottXS = GetMottXS(fEs,fTh); 
   double XS     = MottXS*( (1./Nu)*F2 + 2.*TAN2*(1./PROTON_MASS)*F1 );
   return XS; 
}
//________________________________________________________________________
double F1F209::GetQ2(double Es,double Ep,double th){
   double thr  = th*DEG_TO_RAD; 
   double SIN  = sin(thr/2.); 
   double SIN2 = SIN*SIN;
   double Q2   = 4.*Es*Ep*SIN2; 
   return Q2;  
}
//________________________________________________________________________
double F1F209::GetW(double Es,double Ep,double th){
   double Nu = Es-Ep;
   double Q2 = GetQ2(Es,Ep,th); 
   double W2 = PROTON_MASS*PROTON_MASS + 2.*PROTON_MASS*Nu - Q2;
   return sqrt(W2);
}
