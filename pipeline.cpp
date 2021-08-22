#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstring>
#include <omp.h>

#include <sstream>
//газ аргон
//double sigma = 3.41*pow(10.0, -10);
// ener = 119.8* 1.38*10^23
//double rho = 840.0 * pow(sigma, 3); //не уверен что это разумное значение

double t = 1;
double dt = 0.01;

void progress(int total, int ready) { //декоративная штука
	system("clear");
	int totaldotz = 40;
	double progr = double(ready) / double(total);
	int dots = int(progr * totaldotz);
	int ii = 0;
	std::cout << int(progr * 100) << "% [";
	for (int i = 0; i < dots; i++) std::cout << "#";
	for (int i = 0; i < totaldotz - dots - 1 ; i++) std::cout << " ";
	std::cout << "]" << std::endl;
}



double sigma(double * arr, int n) {
	double aver = 0, mean = 0;
	for (int i = 0; i < n; i++) aver += arr[i] / double(n);
	for (int i = 0; i < n; i++) mean += (arr[i] * arr[i] - aver * aver) / double(n);
	return mean;
}




class model {

public:
	double rho;
	double *rx;
	double *ry;
	double *rz;
	int *nx;
	int *ny;
	int *nz;
	double *vx;
	double *vy;
	double *vz;
	double *fx;
	double *fy;
	double *fz;
	int numofpart;
	double V;
	double L;
	double dt;
	double *EN;
	double *Pp;
	double en = 0.0, pot_en = 0.0, kin_en = 0.0, T0;
	double rr3, ecut;
	double vir = 0.0, ecor, pcor, rc2;


	double P = 0, P1 = 0, div = 0;
	double dt2 = dt * dt;

	std::string name;



	model(int numofpart0, double rho0 = 0.8, double dt0 = 0.01, double T = 2.0, double rc = 1) {
		dt = dt0;
		numofpart = numofpart0;
		rho = rho0;
		T0 = T;
		rc2 = rc;
	}

	// int getN() {return numofpart;}

	void nullv() {
		for (int i = 0; i < numofpart; i++) {
			vx[i] = 0;
			vy[i] = 0;
			vz[i] = 0;
		}
	}

	void init(bool use_e_corr, int timesteps) {
		rx = new double [numofpart];
		ry = new double [numofpart];
		rz = new double [numofpart];
		nx = new int [numofpart];
		ny = new int [numofpart];
		nz = new int [numofpart];
		vx = new double [numofpart];
		vy = new double [numofpart];
		vz = new double [numofpart];
		fx = new double [numofpart];
		fy = new double [numofpart];
		fz = new double [numofpart];
		EN = new double [timesteps];
		Pp = new double [timesteps];
		rr3 = 1.0 / (rc2 * rc2 * rc2);
		ecut = 4 * (rr3 * rr3 * rr3 * rr3 - rr3 * rr3);
		ecor = use_e_corr ? 8 * M_PI * rho * (rr3 * rr3 * rr3 / 9.0 - rr3 / 3.0) : 0.0;
		pcor = use_e_corr ? 16.0 / 3.0 * M_PI * rho * rho * (2. / 3.*rr3 * rr3 * rr3 - rr3) : 0.0;
		V = numofpart / rho;
		L = pow(V, 0.3333333);
		std::cout << L;
	}


	void place() { // размещаем частицы в сетку, обрезаем количество до полного куба
		std::cout << L;

		double lenx = L;
		double leny = L;
		double lenz = L;
		double dist = pow(double(L * L * L) / double(numofpart), 1.0 / 3.0);
		int tmpnum = 0;

		for (int itx = 0; itx < int(lenx / dist); itx++)
			for (int ity = 0; ity < int(leny / dist); ity++)
				for (int itz = 0; itz < int(lenz / dist); itz++) {
					rx[tmpnum] = (itx + 0.5) * dist;
					ry[tmpnum] = (ity + 0.5) * dist;
					rz[tmpnum] = (itz + 0.5) * dist;
					tmpnum++;
				}

		gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
		unsigned long int Seed = 23410981;
		gsl_rng_set(r, Seed);
		for (int i = 0; i < tmpnum; i++) {
			vx[i] = gsl_ran_exponential(r, 3.0);
			vy[i] = gsl_ran_exponential(r, 3.0);
			vz[i] = gsl_ran_exponential(r, 3.0);
			nx[i] = 0;
			ny[i] = 0;
			nz[i] = 0;
		}

		double cmvx = 0, cmvy = 0, cmvz = 0;
		for (int i = 0; i < tmpnum; i++) {
			cmvx += vx[i];
			cmvy += vy[i];
			cmvz += vz[i];
		}
		for (int i = 0; i < tmpnum; i++) {
			vx[i] -= cmvx / tmpnum;
			vy[i] -= cmvy / tmpnum;
			vz[i] -= cmvz / tmpnum;
			kin_en += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
		}
		double T = kin_en / tmpnum * 2. / 3.;
		double fac = sqrt(T0 / T);

		if (T0 != 0.0) {
			kin_en = 0.0;
			for (int i = 0; i < tmpnum; i++) {
				vx[i] *= fac;
				vy[i] *= fac;
				vz[i] *= fac;
				kin_en += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
			}
		}
		numofpart = tmpnum;

	}

	void xyz_in() {
		std::ifstream fin;
		double tmp;
		fin.open("liq.xyz");
		if (!fin.is_open()) std::cout << "ошибка открытия";
		std::getline(fin, name);
		std::getline(fin, name);
		int i = 0;
		while (i < numofpart) {
			fin >> tmp >> rx[i] >> ry[i] >> rz[i] >> vx[i] >> vy[i] >> vz[i];
			i++;
			//std::cout<<"yes";
		}
		fin.close();
	}



	void total_e () {
		int i, j;
		double dx, dy, dz, r2, r6i;
		double e = 0.0, hL = L / 2.0, f;


		for (i = 0; i < numofpart; i++) fx[i] = fy[i] = fz[i] = 0.0;
		vir = 0.0;
		pot_en = 0.0;
		#pragma omp parallel for private(dx, dy, dz, r2, r6i, f, j) reduction(+:vir, e)
		for (i = 0; i < (numofpart - 1); i++) {
			for ( j = i + 1; j < numofpart; j++) {
				dx  = (rx[i] - rx[j]);
				dy  = (ry[i] - ry[j]);
				dz  = (rz[i] - rz[j]);
				if (dx > hL)       dx -= L;
				else if (dx < -hL) dx += L;
				if (dy > hL)       dy -= L;
				else if (dy < -hL) dy += L;
				if (dz > hL)       dz -= L;
				else if (dz < -hL) dz += L;
				r2 = dx * dx + dy * dy + dz * dz;
				if (r2 < rc2) {
					r6i   = 1.0 / (r2 * r2 * r2);
					e    += 4 * (r6i * r6i - r6i) - ecut;
					f     = 48 * (r6i * r6i - 0.5 * r6i);
					vir += f;
					#pragma omp atomic
					fx[i] += dx * f / r2;
					#pragma omp atomic
					fx[j] -= dx * f / r2;
					#pragma omp atomic
					fy[i] += dy * f / r2;
					#pragma omp atomic
					fy[j] -= dy * f / r2;
					#pragma omp atomic
					fz[i] += dz * f / r2;
					#pragma omp atomic
					fz[j] -= dz * f / r2;


				}
			}
		}

		pot_en = e + numofpart * ecor;
	}


	void move() {
		double lenx = L;
		double leny = L;
		double lenz = L;
		//двигаем все
		#pragma omp parallel for
		for (int i = 0; i < numofpart; i++) {
			rx[i] += vx[i] * dt + 0.5 * dt * dt * fx[i];
			ry[i] += vy[i] * dt + 0.5 * dt * dt * fy[i];
			rz[i] += vz[i] * dt + 0.5 * dt * dt * fz[i];
			if (rx[i] > lenx) {rx[i] -= lenx; nx[i]++;}
			if (ry[i] > leny) {ry[i] -= leny; ny[i]++;}
			if (rz[i] > lenz) {rz[i] -= lenz; nz[i]++;}
			if (rx[i] < 0.0) {rx[i] += lenx; nx[i]--;}
			if (ry[i] < 0.0) {ry[i] += leny; ny[i]--;}
			if (rz[i] < 0.0) {rz[i] += lenz; nz[i]--;}
		}
	}


	void addvel() { //меняем скорости
		kin_en = 0;
		#pragma omp parallel for reduction(+:kin_en)
		for (int i = 0; i < numofpart; i++) {
			kin_en += (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]) / 2.0;
			//меняем скорости
			vx[i] += 0.5 * dt * fx[i];
			vy[i] += 0.5 * dt * fy[i];
			vz[i] += 0.5 * dt * fz[i];
		}
	}


	void print_ev(int timestep, bool expanded) { //создает файлы .xyz
		std::ofstream fout;
		std::string str, def = ".xyz", name;
		std::string s = std::to_string(timestep);
		name = s + def;
		fout.open(name);
		if (!fout.is_open()) std::cout << "ошибка открытия";
		fout << numofpart << std::endl;
		fout << "Lattice=\" " << L << " 0.0 0.0 0.0 " << L << " 0.0 0.0 0.0 " << L << " \"\n";
		//fout << "num" << '\t' << "x" << '\t' << "y" << '\t' << "z" << '\t' << "vx" << '\t' << "vy" << '\t' << "vz" << std::endl;
		if(expanded == 0)
			for (int i = 0; i < numofpart; i++) {
				fout << 16 << '\t' << rx[i] << '\t' << ry[i] << '\t' << rz[i] << '\t' << vx[i] << '\t' << vy[i] << '\t' << vz[i] << std::endl;
			}
		else
			for (int i = 0; i < numofpart; i++) {
				fout << 16 << '\t' << rx[i] + nx[i]*L << '\t' << ry[i] + ny[i]*L << '\t' << rz[i] + nz[i]*L<< '\t' << vx[i] << '\t' << vy[i] << '\t' << vz[i] << std::endl;
			}
		fout.close();
	}


	void print_scal(int t) {

		std::ofstream fout;
		fout.open("scalars.txt", std::ios::app);
		if (!fout.is_open()) std::cout << "ошибка открытия scalars";
		fout << '\t' << t << '\t' << pot_en + kin_en << '\t' << pot_en << '\t' << kin_en << "\t" << kin_en * (2.0 / 3.0) / numofpart << "\t" << rho*kin_en * 2. / 3. / numofpart + vir / 3.0 / V << std::endl; //температура последняя
		fout.close();


		P = rho * kin_en * 2. / 3. / numofpart + vir / 3.0 / V;
		P1 = (rho * kin_en * 2. / 3. / numofpart) / P;
		div = kin_en / (pot_en + kin_en);

		EN[t] = pot_en + kin_en;
		Pp[t] = P;

		fout.open("sost.txt", std::ios::app);
		if (!fout.is_open()) std::cout << "ошибка открытия sost";
		fout << P << '\t' << div << '\t' << P1 << std::endl; //температура последняя
		fout.close();
	}

	void esterren( int n) {
		//for(int i = 0 ; i < n; i ++) std::cout << EN[i] << " " << std::endl;
		double * errs = new double[14];
		double * errsn = new double[14];

		std::ofstream fout;
		fout.open("errs.txt");
		errs[0] = sigma(EN, n);
		errsn[0] = errs[0] / (double(n) - 1.0);
		fout << errsn[0] << std::endl;
		std::cout << "ошибка " << errs[0] << "ошибка исп " << errsn[0] << " " << "n " << n  << std::endl ;
		for (int j = 1; j < 14; j++) {

			//double * newarr = new double[n];
			if (n > 10) {
				n = n / 2;
				for (int i = 0; i < n; i++) EN[i] = (EN[2 * i] + EN[2 * i + 1]) / 2.0;

				errs[j] = sigma(EN, n);
				errsn[j] = errs[j] / (double(n) - 1.0);
				std::cout << "ошибка " << errs[j] << "ошибка исп " << errsn[j] << " " << "n " << n  << std::endl ;

				//delete [] newarr;
				fout << errsn[j] << '\t' << sqrt(2 * pow(errs[j], 2) / pow((double(n) - 1.0), 3)) << std::endl;
			}
		}
		fout.close();
		delete [] errs, errsn;
	}

	void printp() {
		double paver = 0;
		#pragma parallel for
		for (int i = 800; i < 1000; i++) paver += Pp[i] / 200.0;
		std::ofstream fout;
		fout.open("pressur.txt", std::ios::app);
		if (!fout.is_open()) std::cout << "ошибка открытия pressure";
		fout << rc2 << '\t' << paver << std::endl; //температура последняя
		fout.close();
	}


	void therm() {
		double lambda = sqrt(1 + 0.1 * (T0 / (2.0 * kin_en / 3.0 / numofpart) - 1.0));
		kin_en = 0.0;
		#pragma omp parallel for reduction(+:kin_en)
		for (int i = 0; i < numofpart; i++) {
			vx[i] *= lambda;
			vy[i] *= lambda;
			vz[i] *= lambda;
			kin_en += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
		}
		kin_en *= 0.5;
	}

};

int get_numofpart(double &L) {
	int numofpart;
	std::string n;
	std::ifstream fin;
	fin.open("liq.xyz");
	if (!fin.is_open()) std::cout << "ошибка открытия";
	std::string::size_type sz;

	fin >> numofpart;
	fin >> n;
	fin >> n;
	L = std::stod (n, &sz);
	fin.close();
	return numofpart;
}



int main(  int argc, char * argv[]) {
	omp_set_num_threads(4);
	bool nishere = 0;
	bool tishere = 0;
	bool dtishere = 0;
	bool use_e_corr = 1;
	bool wannaddt = 0;
	bool expanded = 1;
	double rho = 0.8;
	double use_preset = 0;
	setlocale(LC_ALL, "ru");
	int numofpart, timesteps;

	/*
	for (int i=1;i<argc;i++) {
	    if (!strcmp(argv[i],"-N")) {
	        numofpart=atoi(argv[++i]);
	        nishere = 1;
	    }
	    else if (!strcmp(argv[i],"-dt")) {
	        dt=atof(argv[++i]);
	        dtishere =1;
	    }
	    else if (!strcmp(argv[i],"-t")) {
	        t=atof(argv[++i]);
	        tishere =1;
	    }
	}
	if(!nishere){
	std::cout<<"введите желаемое колво частиц"<< std::endl;
	std::cin>>numofpart;
	}
	 if(!dtishere){
	    std::cout<<"введите количество шагов по времени"<< std::endl;
	    std::cin>>timesteps;
	}
	else{
	    timesteps = t/dt;
	}
	*/
	double rc2 = 100;
	timesteps = 1000;

	if (use_preset) {
		double L;
		numofpart = get_numofpart(L);
		double V = pow(L, 3);
		rho = numofpart / V;
	}
	else numofpart = 625;


	model *mod;
	mod = new model(numofpart, rho , 0.0001,  2.0, rc2);
	mod->init(use_e_corr, timesteps);

	if (use_preset) mod->xyz_in();
	else mod->place();
	//mod->nullv();

	mod->total_e();


	for (int t = 0; t < timesteps; t++) {
		mod->move();
		mod->addvel();
		mod->total_e();
		mod->addvel();
		mod->therm();
		if (t % 10 == 0) mod->print_ev(t, expanded);
		mod->print_scal(t);

		progress(timesteps, t);
	}
	//mod->printp();
	//mod->esterren(timesteps);
	delete mod;

	//fluct(timesteps, dt);
	system("python visualise.py"); //график энергий (тут пока без уточнения насчет обрезки потенциала)
	return 0;
}
// запускать
//g++ clapcontr.cpp -lm -lgsl -lgslcblas -fopenmp; ./a.out ; ovito ; rm  *.xyz
//rm  *.xyz scalars.txt sost.txt errs.txt  ddt.txt