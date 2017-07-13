#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <vector>
#include <string>


#pragma warning(disable:4996)

using namespace std;

const int N = 21; //Size of the grid
const int FLIP_NO = 1000; //number of spins to flip a tick
const int TICKS = 1000; //number of ticks

int grid[N][N]; //Grid of spins

				//The following are the only options for delta-E of flipped spins
double exp4; //Exponent of e^(4*-1/T)
double exp8; //Exponent of e^(8*-1/T)

void printgrid() {
	// Prints out the grid to the terminal
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			cout << grid[i][j] << ' ';
		}
		cout << endl;
	}
}

void initialise() {
	// Initialise the grid with -1 or 1 in each place
	int random;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {

			random = rand() % 2;

			if (random == 0) {
				grid[i][j] = 1;
			}
			else {
				grid[i][j] = -1;
			}
		}
	}
}

int borders(int a) {
	// Ensures periodic boundary conditions
	if (a<0) {
		a += N;
	}
	else {
		if (a >(N - 1)) {
			a -= N;
		}
	}
	return a;
}

int energy(int x, int y) {
	// Returns the potential energy increase of flipping a spin at x,y
	int result = grid[x][y] * (grid[borders(x + 1)][y] + grid[borders(x - 1)][y] + grid[x][borders(y + 1)] + grid[x][borders(y - 1)]);
	return result;
}

void tick(bool equil) {
	// Flip at most FLIP_NO spins, or 1 if equilibrating is on i.e. if equil is true
	int random;
	int en;
	double ran;
	int flips;

	// Set the number of spins to attempt to flip
	// Use 1 for equilibriation so it tries to flip exactly N^3 times
	if (equil) {
		flips = 1;
	}
	else {
		flips = FLIP_NO;
	}

	/* This is the main physics, attempts to flip a random spin
	and accepts it only if the energy change is beneficial or with
	a probability given by the Boltzmann weight factor */

	for (int i = 0; i < flips; i++) {

		int x = rand() % N;
		int y = rand() % N;

		en = 2 * energy(x, y);          // Calculate delta-E for flipping (x,y)
		ran = rand() / (1.0 * RAND_MAX);

		if (en <= 0) {

			grid[x][y] *= -1;
		}
		else {

			// Test each possible value of delta-E, this is only valid for 1, -1 spins
			if (en == 4 && exp4 > ran) {
				grid[x][y] *= -1;
			}
			if (en == 8 && exp8 > ran) {
				grid[x][y] *= -1;
			}
		}
	}
}

int mag() {
	// Sums all the spins to give magnetization
	int sum = 0;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			sum += grid[i][j];
		}
	}
	return sum;
}

int tot_energy() {
	int en = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			en += -1 * energy(i, j);
		}
	}
	return en;
}

int main() {

	//Seed random 
	srand(time(NULL));


	// Set up the temperatures

	int temp_steps;     //Number of steps of temperatures

	cout << "Number of T steps: ";
	cin >> temp_steps;



	//cout << temp_steps << '\n';

	vector<double> T(temp_steps); // Temperature, set by input
	for (int i = 0; i < temp_steps; i++) {
		T[i] = 1.0 + (3.0 * (1.0*(i)) / (1.0*temp_steps));
	}


	// Open data file
	ofstream file("ising1.dat");

	//Loop over all of the temperatures, and run each one
	for (int t = 0; t < temp_steps; t++) {

		initialise();

		double magn = 0;       //Sum of magnetizations to be averaged later
		double en = 0;         //Sum of energies to be averaged later
		double spec = 0;       //Sum to produce specific heat capacity
		double susc = 0;       //Sum to produce susceptibility

		double ene1 = 0;       //Temporary energy variable
		double mag1 = 0;       //Temporary magnetisation variable

		exp4 = exp((-1.0 * 4) / (T[t]));
		exp8 = exp((-1.0 * 8) / (T[t]));

		/* Tick the simulation N^3 times in order to equilibrate, or twice that if under T=1.75
		there seems to be a lot of noise in the lower regime, not sure why - but this helps a bit */

		int equi_ticks = 0;
		if (T[t] <1.75) {
			equi_ticks = (N*N*N);
		}
		else {
			equi_ticks = (N*N*N);
		}

		for (int i = 0; i < equi_ticks; i++) {
			tick(true);
		}

		// Now tick TICKS times and sum averages as it goes along 

		for (int i = 0; i < TICKS; i++) {
			tick(false);
			ene1 = (tot_energy() / (4.0));
			mag1 = (mag() / (1.0));
			magn += mag1;
			en += ene1;
			spec += ene1*ene1;
			susc += mag1*mag1;
		}

		// Put all the data in the file 

		file << T[t] << ' ' << abs((magn / (TICKS*N*N))) << ' ' << (en / (TICKS*N*N)) << ' ';
		file << ((spec / TICKS - en*en / (TICKS*TICKS)) / (N*T[t] * T[t])) << ' ' << ((susc / TICKS - magn*magn / (TICKS*TICKS)) / (N*T[t] * T[t])) << '\n';

	}
	//cout << run_time123 << '\n';

	file.close();

	//Archive the data

	time_t t = time(0);
	struct tm * now = localtime(&t);
	char buffer[200];
	strftime(buffer, 200, "%Y-%m-%d_at_%H-%M_", now);        //17 characters here
	string str(buffer);
	string newstr = str + "with_" + to_string(temp_steps) + "_steps.dat";
	ifstream src("ising1.dat");
	ofstream dst("Archive/" + newstr);
	dst << src.rdbuf();


	//Plot the result with gnuplot
	system("gnuplot -p -e \" load 'plotscript' \" ");

}



