#include <cstdlib> 
#include <cmath> 
#include <iostream> 
#include <iomanip> 
#include <fstream> 
#include <vector> 
#include "./include/eInclusiveCrossSection.h"
#include "./include/F1F209.h"
#include "./include/Radiator.h"

using namespace std; 

double gStep=0; 

void SetEp(double,double,vector<double> &);
void PrintToFile(double,double,vector<double>,vector<double>,vector<double>);

int main(){

	// NOTE: Es, Ep are in GeV 

        // double Es,EpMin,EpMax,th,Z,A,Ta,Tb; 
        double Z,A,Es,Ep,Th; 
	double Tb,Ta;
 
	cout << "Enter Z: ";
	cin  >> Z;
	cout << "Enter A (g/mol): ";
	cin  >> A; 
        cout << "Enter beam energy (GeV): ";
	cin  >> Es;
	cout << "Enter momentum (GeV): ";
	cin  >> Ep; 
	cout << "Enter scattering angle (deg): ";
	cin  >> Th;
        cout << "Enter radiation length before scattering (Tb, in #X0): "; 
        cin  >> Tb; 
        cout << "Enter radiation length after scattering (Ta, in #X0): ";
        cin  >> Ta; 

	F1F209 F1F2; 
        eInclusiveCrossSection *CS = &F1F2; 
	CS->SetZ(Z);
	CS->SetA(A); 
	CS->SetEs(Es);
	CS->SetEp(Ep);
	CS->SetTh(Th);

        Radiator *Rad = new Radiator(); 
        Rad->SetCrossSection(CS); 
        Rad->SetTb(Tb);
        Rad->SetTa(Ta);

        double MottXS = CS->GetMottXS(Es,Th); 
	double BornXS = CS->GetBornXS();
	double RadXS  = Rad->Radiate();

        cout << "Mott: " << MottXS << " nb"        << endl;
        cout << "Born: " << BornXS << " nb/GeV/sr" << endl;
        cout << "Rad:  " << RadXS  << " nb/GeV/sr" << endl;

        delete Rad; 

	return 0;

}
//_____________________________________________________________
void SetEp(double EpMin,double EpMax,vector<double> &Ep){

	double v = EpMin; 
	while(v < EpMax){
		Ep.push_back(v);
		v += gStep;
	}

}
//_____________________________________________________________
void PrintToFile(double Es,double th,vector<double> Ep,vector<double> XSBorn,vector<double> XSRad){

        int N = Ep.size();
	int Pass=0; 
	if(Es==4.73) Pass = 4; 
	if(Es==5.89) Pass = 5; 

        char fn[2148]; 
	sprintf(fn,"./output/F1F209/xs-%d-%.0f.dat",Pass,th); 

	ofstream outfile;
	outfile.open(fn);
	if(outfile.fail()){
		cout << "Cannot open the file: " << fn << endl;
		exit(1);
	}else{
		for(int i=0;i<N;i++){
			outfile << scientific << setprecision(4) << Ep[i]     << "\t" 
                                << scientific << setprecision(4) << XSBorn[i] << "\t" 
                                << scientific << setprecision(4) << XSRad[i]  << endl;
		}
		outfile.close();
		cout << "The data has been written to the file: " << fn << endl;
	}


}
