#include <cassert>
#include <iostream>
#include <fstream>
#include <stdlib.h> //rand()
#include <sstream> // std::istringstream
#include <time.h> // random seed
#include <vector> //array.size
#include <math.h>  // floor
#include <stdio.h> // printf
#include <unistd.h>

int nthermstate = 0;
int ntemp = 0;
int nwalkers = 0;

std::vector<double> Hlambda;
std::vector<double> Temper;

struct energypair {
	double benergy;
	double penergy;
};

std::vector<std::vector<std::vector<energypair> > > energypairlist;

class RandomWalker {
private:
	int indexT, indexL, eindex;
	energypair energy;

public:
	void init_walker(int tempid, int stateid);
	void change_index(int newindexT, int newindexL, int neweindex);
	int get_indexT();
	int get_indexL();
	int get_eindex();
	double get_penergy();
	double get_benergy();	
	void walking();
};

int RandomWalker::get_indexL() {
	return indexL;
}

int RandomWalker::get_indexT() {
	return indexT;
}

int RandomWalker::get_eindex() {
	return eindex;
}
																			  
double RandomWalker::get_penergy() {
	return energy.penergy;
}

double RandomWalker::get_benergy() {
	return energy.benergy;
}

void RandomWalker::change_index(int newindexT, int newindexL, int neweindex) {
	indexT = newindexT;
	indexL = newindexL;
	eindex = neweindex;
}

void RandomWalker::init_walker(int tempid, int stateid) {
	double randnum;
	
	indexT = tempid;
	indexL = stateid;
	randnum = (double)rand()/(double)(RAND_MAX + 1.0);
	eindex = floor(randnum*energypairlist[indexT][indexL].size());
	assert(eindex < energypairlist[indexT][indexL].size());
	energy.penergy = energypairlist[indexT][indexL][eindex].penergy;
	energy.benergy = energypairlist[indexT][indexL][eindex].benergy;
}

void RandomWalker::walking() {
	double randnum;

	randnum = (double)rand()/(double)(RAND_MAX + 1.0);
	eindex = floor(randnum*energypairlist[indexT][indexL].size());
	assert(eindex < energypairlist[indexT][indexL].size());
	energy.penergy = energypairlist[indexT][indexL][eindex].penergy;
	energy.benergy = energypairlist[indexT][indexL][eindex].benergy;	
}

class RepExSystem {
private:
	// RandomWalker walker[nwalkers];
	std::vector<RandomWalker> walker;
	std::vector<std::vector<double> > ratioPFT;	
	std::vector<std::vector<double> > ratioPFL;	
public:
	void init_RepExSystem(int horien);
	void printout(float temper, int nlambda);
	void printFE(float temper, int icycle);
	void print4ST(int icycle);	
	void calPF();
	void update_database(int mwalker, int nwalker);
	void swap(int m, int n);
	void exchange_attempts(int stype, int nattempts, int ncycle);
	void MDupdate(int nstep);
};

void RepExSystem::init_RepExSystem(int horien) {
	for (int i=0; i<ntemp; i++) {
		for (int j=0; j<nthermstate; j++) {
			RandomWalker tmpwalker;
			tmpwalker.init_walker(i, j);
			walker.push_back(tmpwalker);
			nwalkers++;
		}
	}
	assert(nwalkers == (nthermstate*ntemp));
	for (int i=0; i<ntemp; i++) {
		std::vector<double> oneDarray(nthermstate, 0.0);
		ratioPFT.push_back(oneDarray);
	}
	for (int i=0; i<ntemp; i++) {
		std::vector<double> oneDarray(nthermstate, 0.0);
		ratioPFL.push_back(oneDarray);
	}	
}


void RepExSystem::MDupdate(int nstep) {
	for (int i=0; i<nstep; i++) {
		for (int j=0; j<nwalkers; j++) {
			walker[j].walking();
		}
	}
}

void RepExSystem::update_database(int mwalker, int nwalker) {
	int mindexT, mindexL, meindex, nindexT, nindexL, neindex;
	double mpenergy, mbenergy, npenergy, nbenergy;
	
	mindexT = walker[mwalker].get_indexT();
	mindexL = walker[mwalker].get_indexL();
	meindex = walker[mwalker].get_eindex();
	mpenergy = walker[mwalker].get_penergy();
	mbenergy = walker[mwalker].get_benergy();
	nindexT = walker[nwalker].get_indexT();
	nindexL = walker[nwalker].get_indexL();
	neindex = walker[nwalker].get_eindex();
	npenergy = walker[nwalker].get_penergy();
	nbenergy = walker[nwalker].get_benergy();

	energypairlist[mindexT][mindexL][meindex].penergy = npenergy;
	energypairlist[mindexT][mindexL][meindex].benergy = nbenergy;
	energypairlist[nindexT][nindexL][neindex].penergy = mpenergy;
	energypairlist[nindexT][nindexL][neindex].benergy = mbenergy;
}

void RepExSystem::exchange_attempts(int stype, int nattempts, int ncycle) {
	double factor;
	int m, n, swappair, mindexT, mindexL, nindexT, nindexL;
	double mpenergy, mbenergy, npenergy, nbenergy;
	double randnum;
	
	if (stype == 0) {
		if (nattempts < 0) {
			nattempts *= -nthermstate;
		}
		for (int i=0; i<nattempts; i++) {
			randnum = (double)rand()/(double)(RAND_MAX + 1.0);
			m = floor(randnum*nwalkers);
			assert(m<nwalkers);
			randnum = (double)rand()/(double)(RAND_MAX + 1.0);
			n = floor(randnum*nwalkers);
			assert(n<nwalkers);
			while (n == m) {
				randnum = (double)rand()/(double)(RAND_MAX + 1.0);
				n = floor(randnum*nwalkers);
				assert(n<nwalkers);
			}
			swappair = 1;

			mindexT = walker[m].get_indexT();
			mindexL = walker[m].get_indexL();
			mpenergy = walker[m].get_penergy();
			mbenergy = walker[m].get_benergy();
			nindexT = walker[n].get_indexT();
			nindexL = walker[n].get_indexL();
			npenergy = walker[n].get_penergy();
			nbenergy = walker[n].get_benergy();

			factor = exp(-(1.0/0.0019872041/Temper[mindexT] -
						   1.0/0.0019872041/Temper[nindexT])*(npenergy - mpenergy)
						 - (Hlambda[mindexL]/0.0019872041/Temper[mindexT]
							- Hlambda[nindexL]/0.0019872041/Temper[nindexT])
						 *(nbenergy - mbenergy));
			if (factor < 1) {
				randnum = (double)rand()/(double)(RAND_MAX + 1.0);
				if (randnum > factor) {
					swappair = 0;
				}
			}
			if (swappair == 1) {
				update_database(m, n);
				swap(m, n);
			}
		}
	}
}


void RepExSystem::calPF() {
	std::vector<std::vector<energypair> > calenergy(ntemp, std::vector<energypair>(nthermstate));	
	
	for (int i=0; i<nwalkers; i++) {
		int pindex = walker[i].get_indexT()*nthermstate + walker[i].get_indexL();
		int iindexT = walker[i].get_indexT();
		int iindexL = walker[i].get_indexL();
		calenergy[iindexT][iindexL].benergy = walker[i].get_benergy();
		calenergy[iindexT][iindexL].penergy = walker[i].get_penergy();
	}
	for (int i=0; i<ntemp; i++) {
		for (int j=0; j<nthermstate; j++) {
			if (j==0) {
				ratioPFT[i][j] += 1.0;
			}
			else {
				ratioPFT[i][j] += exp(-(1.0/0.0019872041/Temper[i])*(Hlambda[j] - Hlambda[j-1])*calenergy[i][j-1].benergy);
			}
		}
	}
	for (int j=0; j<nthermstate; j++) {
		for (int i=0; i<ntemp; i++) {
			if (i==0) {
				ratioPFL[i][j] += 1.0;
			}
			else {
				ratioPFL[i][j] += exp(-(1.0/0.0019872041/Temper[i] - 1.0/0.0019872041/Temper[i-1])*
									  (calenergy[i-1][j].penergy + Hlambda[j]*calenergy[i-1][j].benergy));
			}
		}
	}
}

void RepExSystem::printFE(float temper, int icycle) {
	for (int i=0; i<ntemp; i++) {
		if ((abs(temper-0.0)<1.0e-10) || (abs(temper-Temper[i])<1.0e-10)) {		
			double deltaG = 0.0;
			for (int j=0; j<nthermstate; j++) {
				deltaG += -0.0019872041*Temper[i]*log(ratioPFT[i][j]/icycle);
				printf("%12g ", deltaG);			
			}
			printf("\n");
		}
	}
	printf("\n");
}

void RepExSystem::print4ST(int icycle) {
	std::vector<std::vector<double> > deltaGL(ntemp, std::vector<double>(nthermstate));	
	for (int i=0; i<ntemp; i++) {
		deltaGL[i][0] = log(ratioPFT[i][0]/icycle);
		for (int j=1; j<nthermstate; j++) {
			deltaGL[i][j] = deltaGL[i][j-1] + log(ratioPFT[i][j]/icycle);
		}
	}

	std::vector<std::vector<double> > deltaGT(ntemp, std::vector<double>(nthermstate));
	for (int j=0; j<nthermstate; j++) {
		deltaGT[0][j] = log(ratioPFL[0][j]/icycle);
		for (int i=1; i<ntemp; i++) {
			deltaGT[i][j] = deltaGT[i-1][j] + log(ratioPFL[i][j]/icycle);
		}
	}

	std::vector<std::vector<double> > deltaG(ntemp, std::vector<double>(nthermstate, 0.0));
	for (int i=0; i<ntemp; i++) {
		for (int j=0; j<nthermstate; j++) {
			deltaG[i][j] += (deltaGL[0][j] - deltaGL[0][0]);
			deltaG[i][j] += (deltaGT[i][j] - deltaGT[0][j]);
			deltaG[i][j] += (deltaGT[i][0] - deltaGT[0][0]);
			deltaG[i][j] += (deltaGL[i][j] - deltaGL[i][0]);			
			deltaG[i][j] /= 2.0;
		}
	}


	for (int i=0; i<ntemp; i++) {
		for (int j=0; j<nthermstate; j++) {
			printf("%12g ", deltaG[i][j]);			
		}
		printf("\n");
	}
	printf("\n");
}

void RepExSystem::printout(float temper, int nlambda) {
	std::vector<energypair> printenergy(nwalkers);
	// energypair printenergy[nwalkers]; // this is not good, not safe
	int nprint;
	
	for (int i=0; i<nwalkers; i++) {
		int pindex = walker[i].get_indexT()*nthermstate + walker[i].get_indexL();
		printenergy[pindex].benergy = walker[i].get_benergy();
		printenergy[pindex].penergy = walker[i].get_penergy();
	}
	nprint = 0;
	for (int i=0; i<ntemp; i++) {
		for (int j=0; j<nthermstate; j++) {
			if ((abs(temper-0.0)<1.0e-10) || (abs(temper-Temper[i])<1.0e-10)) {
				if ((nthermstate-j)<=nlambda) {
					printf("%10g %10g ", printenergy[nprint].benergy, printenergy[nprint].penergy);
				}
			}
			nprint++;
		}
	}
	printf("\n");
}


void RepExSystem::swap(int m, int n) {
	int newmindexT, newmindexL, newnindexT, newnindexL, newmeindex, newneindex;

	newmindexT = walker[n].get_indexT();
	newmindexL = walker[n].get_indexL();
	newmeindex = walker[n].get_eindex();
	newnindexT = walker[m].get_indexT();
	newnindexL = walker[m].get_indexL();
	newneindex = walker[m].get_eindex();
	walker[m].change_index(newmindexT, newmindexL, newmeindex);
	walker[n].change_index(newnindexT, newnindexL, newneindex);
}

/*****************************************************************************************************/
int main(int argc, char* argv[]) {
	// std::cout << "RepExchange.o temper nattempts totalcycle" << "\n";
	time_t start,end;
	time (&start);
	
	double dummyx, dummyy, dummyz;
	char* lambdafile;
	char filename[99];
	int icycle, totalcycle, nattempts, nlambda, equilibrium;
	int printfreq=1;
	int stype=0;
	int info = 0;	
	float temper;
	std::ifstream read_file;
	std::string line;		
	RepExSystem mysim;
	int option_char;

	opterr = 0;	
	while ((option_char = getopt(argc, argv, "hf:t:l:n:a:e:q:i:")) != -1) {
		switch (option_char) {
		case 'h':
			printf("%s -f lambdalist -t temperture -l nlambdastate -n totalcycle -a nattempts -e equilibrium -q [printfreq] -i [info]\n", argv[0]);
			return 0;
		case 'f': lambdafile = optarg; break;
		case 't': temper = atof(optarg); break;
		case 'l': nlambda = atoi(optarg); break;
		case 'n': totalcycle = atoi(optarg) ; break;
		case 'a': nattempts = atoi(optarg); break;
		case 'e': equilibrium = atoi(optarg); break;
		case 'q': printfreq = atoi(optarg); break;
		case 'i': info = atoi(optarg); break;
		case '?': 
			if ((optopt == 'f') || (optopt == 't') || (optopt == 'l') || (optopt == 'n') || (optopt == 'a') || (optopt == 'e') || (optopt == 'q') || (optopt == 'i'))  
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
				fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
			return 1;
		default:
			abort ();
		}
	}
	
	srand (time(NULL)); //random seed

	// read lambda and temperature
	read_file.open(lambdafile);	
	assert(read_file.is_open());
	std::getline(read_file, line);
	{
		std::istringstream vstring(line);
		while ( vstring >> dummyx) {
			Temper.push_back(dummyx);
			ntemp++;
		}
	}
	std::getline(read_file, line);
	{
		std::istringstream vstring(line);
		while ( vstring >> dummyx) {
			Hlambda.push_back(dummyx);
			nthermstate++;
		}
	}
	read_file.close();

	// read energy list
	for (int i=0; i<ntemp; i++) {
		std::vector<std::vector<energypair> >  twoDarray;					
		for (int j=0; j<nthermstate; j++) {
			std::vector<energypair>  oneDarray;
			int inttemp = floor(Temper[i]);			
			sprintf(filename, "input/T_%d_lambda_%02d.cluster", inttemp, j);
			read_file.open(filename);
			assert(read_file.is_open());
			while (read_file >> dummyx >> dummyy >> dummyz) {
				energypair tmpenergy;
				tmpenergy.benergy = dummyx;
				tmpenergy.penergy = dummyz;
				oneDarray.push_back(tmpenergy);
			}
			read_file.close();
			twoDarray.push_back(oneDarray);
		}
		energypairlist.push_back(twoDarray);
	}
	
	mysim.init_RepExSystem(0);
	// mysim.printout(temper, nlambda);

	icycle = 0;
	while(icycle < (totalcycle+equilibrium)) {
		mysim.MDupdate(1);
		mysim.exchange_attempts(stype, nattempts, icycle);
		if ((info == 1) || (info == 2)) {
			mysim.calPF();
		}
		icycle++;		
		if ((icycle>equilibrium) && ((icycle-equilibrium)%printfreq==0)) {
			switch (info) {
			case 0: mysim.printout(temper, nlambda); break;
			case 1: mysim.printFE(temper, icycle); break;
			case 2: mysim.print4ST(icycle); break;
			}
		}
	}

	time (&end);
	double dif = difftime (end,start);
	printf ("# Elasped time is %.2lf seconds.\n", dif );
	
	return 0;
}
		
