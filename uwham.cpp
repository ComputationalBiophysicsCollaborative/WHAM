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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

// g++ wham.cpp -lm -lgsl -lgslcblas 

const double kB = 0.0019872041;

int nstates = 0; // how many walks or thermodynamic states
// int npara = 0; // how many parameters
int ndata = 0; // how mnay data for one observation
int totalndata = 0; // total number of observations
std::vector<int> nobs; // how many observations at each state
std::vector<double> density_of_states; // this is the density of states, totalndata
std::vector<double> partition_function;
std::vector<std::vector<double> > bias_matrix; // this is the bias matrix totalndata*nstates

double potential_energy(std::vector<double>& parameter, std::vector<double>& variable) {
	double energy;

	energy = (variable[2] + parameter[1]*variable[0])/(parameter[0]*kB); // potential energy in unit of k_B*T
	return energy;
}

void print_out_vector(gsl_vector *pvector) {
	int ncolumn;

	ncolumn = pvector->size;
	printf("The vector is\n");
	for (int i=0; i<ncolumn; i++) {
		printf("%12g ", gsl_vector_get(pvector, i));
	}
	printf("\n");		
}


void print_out_matrix(gsl_matrix *pmatrix) {
	int nrow, ncolumn;

	nrow = pmatrix->size1;
	ncolumn = pmatrix->size2;
	printf("The matrix is\n");
	for (int i=0; i<nrow; i++)
		for (int j=0; j<ncolumn; j++)
			printf(j==(ncolumn-1)?"%12g\n":"%12g ", gsl_matrix_get(pmatrix,i,j));
}


void convex_newton_raphson_solver(int niteration, double max_error, double factor) {
	int iteration=0;
	double error=1.0e10;
	std::vector<double> oldfe;
	std::vector<double> newfe;

	void print_out_matrix(gsl_matrix *pmatrix);
	double cal_partition_function(int index);
	double cal_density_of_states(int index);
	void normalize_partition_functions();	
	double cal_free_energy(int index);

	double alpha = -1.0/totalndata;			
	int nstatesm1 = nstates - 1;
	int signum;	
	gsl_vector * gradient = gsl_vector_alloc(nstates);	
	gsl_matrix * hessian = gsl_matrix_alloc(nstates, nstates);
	gsl_matrix * weight_matrix = gsl_matrix_alloc(totalndata, nstates);
	gsl_matrix * pimatrix = gsl_matrix_alloc(nstates, nstates);
	gsl_vector * onem = gsl_vector_alloc(nstates);
	gsl_vector * onen = gsl_vector_alloc(totalndata);
	gsl_matrix * tmpm1 = gsl_matrix_alloc(nstates, totalndata);
	gsl_matrix * tmpm2 = gsl_matrix_alloc(nstates, nstates);			
	gsl_vector * tmpv1 = gsl_vector_alloc(nstates);
	gsl_matrix * tmphessian = gsl_matrix_alloc(nstatesm1, nstatesm1);
	gsl_matrix * invhessian = gsl_matrix_alloc(nstatesm1, nstatesm1);		
	gsl_permutation * perm = gsl_permutation_alloc(nstatesm1);
	gsl_vector * tmpgradient = gsl_vector_alloc(nstatesm1);
	gsl_vector * product = gsl_vector_alloc(nstatesm1);		

	// initialize vectors and matrices
	gsl_vector_set_all(onem, 1.0);
	gsl_vector_set_all(onen, 1.0);
	gsl_matrix_set_all(pimatrix, 0.0);
	for (int i=0; i<nstates; i++) {
		double tmpelem = 1.0*nobs[i]/totalndata;
		gsl_matrix_set(pimatrix, i, i, tmpelem);
	}

	// printf("iteration=%-12d", iteration);
	printf("\033[1;32miteration=%-12d\033[0m", iteration);			
	for (int i=0; i<nstates; i++) {
		double fe = cal_free_energy(i);
		printf("%12g ", fe);
		oldfe.push_back(fe);
		newfe.push_back(fe);
	}
	printf("\n");	

	// fix \zeta[0] = 0.0	
	normalize_partition_functions(); 
	// update density of states
	for (int i=0; i<totalndata; i++) {
		density_of_states[i] = cal_density_of_states(i);
	}

	while ((iteration < niteration) && (error > max_error)) {
		iteration++;
		// printf("iteration=%-12d", iteration);
		printf("\033[1;32miteration=%-12d\033[0m", iteration);				

		// calculate the weight matrix
		for (int i=0; i<totalndata; i++) {
			for (int j=0; j<nstates; j++) {
				double tmpelem = density_of_states[i]*bias_matrix[j][i]/partition_function[j]*totalndata;
				gsl_matrix_set(weight_matrix, i, j, tmpelem);
			}
		}
		// calculate the gradient matrix			
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, pimatrix, weight_matrix, 0.0, tmpm1);
		gsl_blas_dgemv(CblasNoTrans, alpha, tmpm1, onen, 0.0, gradient);
		gsl_matrix_set_all(hessian, 0.0);
		for (int i=0; i< nstates; i++) {
			double tmpelem = -1.0*gsl_vector_get(gradient, i);
			gsl_matrix_set(hessian, i, i, tmpelem);
		}
		gsl_blas_dgemv(CblasNoTrans, 1.0, pimatrix, onem, 0.0, tmpv1);
		gsl_vector_add(gradient, tmpv1);
		// calculate the hessian matrix
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmpm1, weight_matrix, 0.0, tmpm2);			
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, alpha, tmpm2, pimatrix, 1.0, hessian);			
		// calculate partition function changes
		for (int i=1; i<nstates; i++) {
			for (int j=1; j<nstates; j++) {
				double tmpelem = gsl_matrix_get(hessian, i, j);
				gsl_matrix_set(tmphessian, (i-1), (j-1), tmpelem);
			}
		}
		gsl_linalg_LU_decomp(tmphessian, perm, &signum);
		gsl_linalg_LU_invert(tmphessian, perm, invhessian);
		// print_out_matrix(invhessian);		
		for (int i=1; i<nstates; i++) {
			double tmpelem = gsl_vector_get(gradient, i);
			gsl_vector_set(tmpgradient, (i-1), tmpelem);
		}
		gsl_blas_dgemv(CblasNoTrans, 1.0, invhessian, tmpgradient, 0.0, product);

		// // update partition function
		factor = ((factor*iteration)<1.0)?(factor*iteration):1.0;
		partition_function[0] = 1.0;
		for (int i=1; i<nstates; i++) {
			partition_function[i] *= exp(-factor*gsl_vector_get(product, (i-1)));
		}
		// // update density of states
		for (int i=0; i<totalndata; i++) {
			density_of_states[i] = cal_density_of_states(i);
		}		

		error = 0.0;
		for (int i=0; i<nstates; i++) {
			double fe = cal_free_energy(i);
			printf("%12g ", fe);
			newfe[i] = fe;
			double tmperror = fabs((newfe[i] - oldfe[i])/oldfe[i]);
			if (tmperror > error) {
				error = tmperror;
			}
			oldfe[i] = newfe[i];			
		}
		printf("\tmax error=%12g", error);
		printf("\n");
	}				
}

void newton_raphson_solver(int niteration, double max_error, double factor) {
	int iteration=0;
	double error=1.0e10;
	std::vector<double> oldfe;
	std::vector<double> newfe;
	
	std::vector<double> nrfunc(nstates, 0.0);
	std::vector<std::vector<double> > weight_matrix(nstates, std::vector<double>(totalndata, 0.0));
	gsl_matrix * jacobian = gsl_matrix_alloc(nstates, nstates);

	void print_out_matrix(gsl_matrix *pmatrix);
	double cal_partition_function(int index);
	double cal_density_of_states(int index);
	void normalize_partition_functions();	
	double cal_free_energy(int index);
		
	// printf("iteration=%-12d", iteration);
	printf("\033[1;32miteration=%-12d\033[0m", iteration);			
	for (int i=0; i<nstates; i++) {
		double fe = cal_free_energy(i);
		printf("%12g ", fe);
		oldfe.push_back(fe);
		newfe.push_back(fe);
	}
	printf("\n");	

	// update density of states
	for (int i=0; i<totalndata; i++) {
		density_of_states[i] = cal_density_of_states(i);
	}

	while ((iteration < niteration) && (error > max_error)) {
		iteration++;
		// printf("iteration=%-12d", iteration);
		printf("\033[1;32miteration=%-12d\033[0m", iteration);				
		// weight matrix
		for (int i=0; i<nstates; i++) {
			for (int j=0; j<totalndata; j++) {
				weight_matrix[i][j] = density_of_states[j]*bias_matrix[i][j]/partition_function[i];
			}
		}
		// g_i(\theta)
		for (int i=0; i<nstates; i++) {
			double sum = 0.0;
			for (int j=0; j<totalndata; j++) {
				sum += weight_matrix[i][j];
			}
			nrfunc[i] = nobs[i] - nobs[i]*sum;
			// printf("%12g\n", nrfunc[i]);
		}
		for (int i=0; i<nstates; i++) {
			for (int j=0; j<nstates; j++) {
				double sum = 0.0;		
				for (int k=0; k<totalndata; k++) {
					if (i == j) {
						sum += -nobs[i]*weight_matrix[i][k]*(1 - nobs[i]*weight_matrix[i][k]);
					}
					else {
						sum += nobs[i]*weight_matrix[i][k]*nobs[j]*weight_matrix[j][k];
					}
				}
				gsl_matrix_set(jacobian, i, j, sum);
			}
		}
		// print_out_matrix(jacobian);

		gsl_matrix * inv_jacobian = gsl_matrix_alloc(nstates, nstates);		
		gsl_permutation * perm = gsl_permutation_alloc(nstates);
		int signum;
		// calculate the matrix inverse
		gsl_linalg_LU_decomp(jacobian, perm, &signum);
		// print_out_matrix(jacobian);
		gsl_linalg_LU_invert(jacobian, perm, inv_jacobian);
		// print_out_matrix(inv_jacobian);

		// // update partition function
		factor = ((factor*iteration)<1.0)?(factor*iteration):1.0; 
		// factor = (iteration==1)?0.1:1.0;		
		for (int i=0; i<nstates; i++) {
			partition_function[i] = log(partition_function[i]);
			for (int j=0; j<nstates; j++) {
				partition_function[i] += factor*gsl_matrix_get(inv_jacobian, i, j)*nrfunc[j];
			}
			partition_function[i] = exp(partition_function[i]);
		}

		normalize_partition_functions();
		// // update density of states
		for (int i=0; i<totalndata; i++) {
			density_of_states[i] = cal_density_of_states(i);
		}		

		error = 0.0;
		for (int i=0; i<nstates; i++) {
			double fe = cal_free_energy(i);
			printf("%12g ", fe);
			newfe[i] = fe;
			double tmperror = fabs((newfe[i] - oldfe[i])/oldfe[i]);
			if (tmperror > error) {
				error = tmperror;
			}
			oldfe[i] = newfe[i];			
		}
		printf("\tmax error=%12g", error);			
		printf("\n");		
	}
}


void cal_bias_matrix(std::vector<std::vector<double> >& datalist, std::vector<std::vector<double> >& paralist, int printtag)  {
	std::vector<std::vector<double> > energy_matrix; // this is the potential matrix totalndata*nstates	
	double potential_energy(std::vector<double>& parameter, std::vector<double>& variable);
	FILE * write_file;		
	
	for (int i=0; i<nstates; i++) {
		double tmpe;
		double tmpb;
		std::vector<double> tmp_potentials;
		std::vector<double> tmp_biases;		
		for (int j=0; j<totalndata; j++) {
			tmpe = potential_energy(paralist[i], datalist[j]);
			tmp_potentials.push_back(tmpe);
			tmpb = exp(-tmpe);
			tmp_biases.push_back(tmpb);
		}
		energy_matrix.push_back(tmp_potentials);
		bias_matrix.push_back(tmp_biases);
	}
	// print out potential energy matrix]
	if (printtag == 1) {
		write_file = fopen("MatrixPE", "w");		
		for (int j=0; j<totalndata; j++) {
			for (int i=0; i<nstates; i++) {
				fprintf(write_file, "%12g\t", energy_matrix[i][j]);
			}
			fprintf(write_file, "\n");
		}
		fclose(write_file);
	}
}

void print_out_alldata(std::vector<std::vector<double> >& datalist) {
	FILE * write_file;	
	
	write_file = fopen("all.data", "w");
	for (int i=0; i<totalndata; i++) {
		for (int j=0; j<ndata; j++) {
			fprintf(write_file, "%12g\t", datalist[i][j]);
		}
		fprintf(write_file, "\n");
	}
	fclose(write_file);	
}


double cal_partition_function(int index) {
	double pf = 0.0;
	for (int j=0; j<totalndata; j++) {
		pf += density_of_states[j]*bias_matrix[index][j];
	}
	return pf;
}

void normalize_partition_functions() {
	for (int i=1; i<nstates; i++) {
		partition_function[i] /= partition_function[0];
	}
	partition_function[0] = 1.0;
}

double cal_density_of_states(int index) {
	double ds = 0.0;

	for (int i=0; i<nstates; i++) {
		ds += nobs[i]*bias_matrix[i][index]/partition_function[i];
	}
	ds = 1.0/ds;
	return ds;
}

double cal_free_energy(int index) {
	double fe=0.0;
	if (index != 0) 
		fe = -log(partition_function[index]/partition_function[0]);
	return fe;
}


void print_out_weights(std::vector<int>& printtag) {
	FILE * write_file;
	std::vector<double> sum(nstates, 0.0);
	double tmpw;

	write_file = fopen("weights.data", "w");
	fprintf(write_file, "# ");
	for (int i=0; i<nstates; i++) {
		if (printtag[i] != 0) {
			int j = i + 1;
			fprintf(write_file, "%12d ", j);
		}
	}
	fprintf(write_file, "\n");	
	for (int j=0; j<totalndata; j++) {
		for (int i=0; i<nstates; i++) {		
			tmpw = density_of_states[j]*bias_matrix[i][j];
			sum[i] += tmpw;
		}
	}
	for (int j=0; j<totalndata; j++) {
		for (int i=0; i<nstates; i++) {
			if (printtag[i] != 0) {			
				tmpw = density_of_states[j]*bias_matrix[i][j]/sum[i];
				fprintf(write_file, "%12g ", tmpw);
			}
		}
		fprintf(write_file, "\n");
	}
	fclose(write_file);	
}

void self_consistent_iteration(int niteration, double max_error) {
	int iteration=0;
	double error=1.0e10;
	std::vector<double> oldfe;
	std::vector<double> newfe;

	double cal_partition_function(int index);
	double cal_density_of_states(int index);
	double cal_free_energy(int index);	
	
	// printf("iteration=%-12d", iteration);
	// printf("\033[1;32miteration=%-12d\033[0m", iteration);		
	for (int i=0; i<nstates; i++) {
		double fe = cal_free_energy(i);
		// printf("%12g ", fe);
		oldfe.push_back(fe);
		newfe.push_back(fe);
	}
	// printf("\n");	

	while ((iteration < niteration) && (error > max_error)) {
		iteration++;		
		// printf("iteration=%-12d", iteration);
		printf("\033[1;32miteration=%-12d\033[0m", iteration);				
		// update density of states
		for (int i=0; i<totalndata; i++) {
			density_of_states[i] = cal_density_of_states(i);
		}
		// update partition function
		for (int i=0; i<nstates; i++) {
			partition_function[i] = cal_partition_function(i);
		}

		error = 0.0;
		for (int i=0; i<nstates; i++) {
			double fe = cal_free_energy(i);
			printf("%12g ", fe);
			newfe[i] = fe;
			double tmperror = fabs((newfe[i] - oldfe[i])/oldfe[i]);
			if (tmperror > error) {
				error = tmperror;
			}
			oldfe[i] = newfe[i];			
		}
		printf("\tmax error=%12g", error);
		printf("\n");		
	}
}

/*****************************************************************************************************/
int main(int argc, char* argv[]) {
	// std::cout << "RepExchange.o temper nattempts totalcycle" << "\n";
	time_t start,end;
	time (&start);

	double dummyx;
	int dummym, dummyn;
	std::string dummys;
	char* parafile;
	const char* potfile="";
	std::ifstream read_file;
	std::string line;
	int niteration = RAND_MAX;
	int nequilibrium = 1;
	int method;
	double tolerance = 1.0e-25;
	double gamma = 1.0;
	int printpe = 0; // print out the potential matrix or not
	int option_char;

	std::vector<int> printtag;

	void print_out_alldata(std::vector<std::vector<double> >& datalist);	
	void cal_bias_matrix(std::vector<std::vector<double> >& datalist, std::vector<std::vector<double> >& paralist, int printtag);
	double cal_partition_function(int index);	
	void self_consistent_iteration(int niteration, double max_error);
	void newton_raphson_solver(int niteration, double max_error, double factor);
	void convex_newton_raphson_solver(int niteration, double max_error, double factor);	
	void print_out_weights(std::vector<int>& printtag);
	
	opterr = 0;	
	while ((option_char = getopt(argc, argv, "hpx:f:m:g:q:n:e:")) != -1) {
		switch (option_char) {
		case 'h':
			printf("%s [-x potential_energy_file] -f parameters_file -m iteration_method [-g converge_factor -q equilibrium_number] -n iteration_number -e tolerance\n", argv[0]);
			return 0;
		case 'p': printpe = 1; break; // print out potential matrix, only for testing
		case 'x': potfile = optarg; break;
		case 'f': parafile = optarg; break;
		case 'm': method = atoi(optarg); break;
		case 'q': nequilibrium = atoi(optarg); break;
		case 'g': gamma = atof(optarg); break;
		case 'n': niteration = atoi(optarg); break;
		case 'e': tolerance = atof(optarg); break;			
		case '?': 
			if ((optopt == 'x') || (optopt == 'f') || (optopt == 'm') || (optopt == 'g') || (optopt == 'q') || (optopt == 'n') || (optopt == 'e'))
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

	// read parameters
	if (potfile[0] == '\0') { 
		std::vector<std::vector<double> > paralist;
		std::vector<std::vector<double> > datalist;
		std::vector<std::string> filename;
		
		read_file.open(parafile);	
		assert(read_file.is_open());
		while (std::getline(read_file, line)) {
			std::vector<double>  paraarray;		
			std::istringstream vstring(line);
			vstring >> dummys;
			filename.push_back(dummys);
			vstring >> dummyn;
			printtag.push_back(dummyn);
			while (vstring >> dummyx) {
				paraarray.push_back(dummyx);
			}
			paralist.push_back(paraarray);
		}
		read_file.close();

		nstates = paralist.size();
		// npara = paralist[0].size();
		// printf("%12d\t%12d\n", nstates, npara);
	
		// read data list
		for (int i=0; i<nstates; i++) {
			nobs.push_back(0);
			const char* tmpname = filename[i].c_str();			
			read_file.open(tmpname);
			assert(read_file.is_open());
			while (std::getline(read_file, line)) {
				std::vector<double>  tmparray;		
				std::istringstream vstring(line);
				while (vstring >> dummyx) {
					tmparray.push_back(dummyx);
				}
				datalist.push_back(tmparray);
				nobs[i]++;
				totalndata++;
			}
			read_file.close();		
		}
		ndata = datalist[0].size();

		// write all data into one file
		print_out_alldata(datalist);

		// calculate the potential matrix
		cal_bias_matrix(datalist, paralist, printpe);		
	}
	// read energy matrix
	else {
		read_file.open(parafile);
		assert(read_file.is_open());
		while (read_file >> dummym >> dummyn) {
			nobs.push_back(dummym);
			printtag.push_back(dummyn);
		}
		read_file.close();
		nstates = nobs.size();
		
		read_file.open(potfile);	
		assert(read_file.is_open());
		std::getline(read_file, line);
		std::istringstream firststring(line);
		while (firststring >> dummyx) {
			std::vector<double> biasarray;							
			biasarray.push_back(exp(-dummyx));
			bias_matrix.push_back(biasarray);
		}
		while (std::getline(read_file, line)) {
			std::istringstream vstring(line);
			int k = 0;
			while (vstring >> dummyx) {
				bias_matrix[k].push_back(exp(-dummyx));
				k++;
			}
		}
		read_file.close();		
		totalndata = bias_matrix[0].size();
	}
	
	
    // initialize density of states
	for (int i=0; i<totalndata; i++) {
		density_of_states.push_back(1.0/totalndata);
	}
	// initialize partition functions 
	for (int i=0; i<nstates; i++) {
		partition_function.push_back(1.0);			
	}
	// run self_consistent_iteration once
	// printf("Equilibration:\n");
	std::cout << "\033[1;31mEquilibration:\033[0m\n";
	self_consistent_iteration(nequilibrium, tolerance);
	
	// iteration
	// printf("Iterations:\n");
	std::cout << "\033[1;31mIterations:\033[0m\n";	
	switch (method) {
	case 0: self_consistent_iteration(niteration, tolerance); break;
	case 1: convex_newton_raphson_solver(niteration, tolerance, gamma); break;						
	case 2: newton_raphson_solver(niteration, tolerance, gamma); break;
	}

	print_out_weights(printtag);
	
	time (&end);
	double dif = difftime (end,start);
	printf ("# Elasped time is %.2lf seconds.\n", dif );
	
	return 0;
}
