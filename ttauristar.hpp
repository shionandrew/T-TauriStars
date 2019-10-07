/**
 * \file ttauristar.hpp
 *
 * \authors Cynthia Yan with additions by Shion Andrew
 *
 * \brief Declares the TTauriStar class.
 */

#ifndef TTAURISTAR_HPP_INCLUDED
#define TTAURISTAR_HPP_INCLUDED

#include <cmath>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <unistd.h>

using namespace std;
/**
 * \class TTauriStar
 *
 * \brief Represents a T Tauri Star throughout its evolution.
 *
 * \details The first parameter is a pointer matrix to the mass-age-radius table,
 *          mass and age are the final values for a particular star,
 *          massdotfactor is the randomization Mdot parameter,
 *          bfieldstrength is the randomization magnetic field strength.
 */
class TTauriStar {
public:

	/**
	 * \brief Parameterized constructor for T Tauri Star
	 *
	 * \param cmktable: a
	 * \param mass: final mass of star
	 * \param age: final age of star
	 * \param massdotfactor:
	 * \param bfieldstrength:
	 */
	TTauriStar(vector<vector<double>> cmktable, double mass, double age,
		double massdotfactor, double bfieldstrength);

	/**
	 * \brief update runs the simulation for each star until convergence.
	 *
	 * \returns period
	 */
	double update();

	/**
	 * \brief getvector fetches various evolutionary vectors for each star,
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns the vector you want
	 */
	vector<double> getvector(int n);

	/**
	 * \brief fetch names of these vectors for use in plotting
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns a string telling the name
	 */
	string getname(int n);

	/**
	 * \brief fetch units of these vectors for use in plotting
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns a string telling the unit
	 */
	string getunit(int n);

	/**
	 * \brief makes plots through gnuplot vector m vs. vector n is plotted.
	 *
	 * \param integer codes of the vectors you want
	 *
	 * \returns
	 */
	void plot(int m, int n);


	/**
	 * \brief finds slope of two vectors, assuming linear relationships
	 *
	 * \param integer codes of the vectors you want
	 *
	 * \returns slope
	 */
	double slope( vector<double>& x, vector<double>& y);


private:
	double calculatemassdot();
	double calculateradius();
	double calculatebfield();
	double calculaterm();
	double calculatediskdensity();
	void calculatemasses();
	void calculateperiods();

	// TTauriStar data members (Matrix)
	vector<vector<double>> cmktable_;
	double mass_;								/// (solar masses)
	double mass0_;              /// (solar mass) mass that we are trying to converge to, equal to the original input mass
	double mass2_;              /// (solar mass) mass calculated going forward to test for convergence against mass0
	double age_;								/// Myears
	double massdotfactor_;      /// equal to the original input parameter
	double bfieldstrength_;			/// (KGauss)
	double massdot_;            /// M_sun/yr
	double period_;							/// days
	double radius_;							/// R_sun
	double bfield_;							/// seems like kG from rm calculation
	double rm_;									/// R_sun
	double periodrm_;						/// days
	double diskdensity_;        /// surface density at R_M, g/cm^2
	double propendtime_;        /// end time for the propeller effect
	double acceff_;             /// fraction of deposition by accretion
	bool valid_;  /// if the mass is > 3, there are not values in the table, so we drop stars with those masses.
	//This should be changed in the future with better data tables.

	vector<double> ages_;
	vector<double> masses_;         /// mass of the protostar changed backward
	vector<double> massdots_;       /// mass accretion rate of the protostar
	vector<double> periods_;        /// period of the protostar
	vector<double> acceffs_;
	vector<double> radii_;
	vector<double> rms_;
	vector<double> diskdensities_;
	vector<double> phase_;					// vector storing "phases" the star. Solely for the purpose of understanding the stages of the star's evolution in our graph
	vector<double> rmPeriods_;			// stores keplerian period calculated at magnetospheric radius (days)
};

#endif
