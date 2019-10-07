/**
 * \file simulation.cpp
 *
 * \authors Cynthia Yan
 *
 * \brief Provides the main() function for simulating a ttauri star cluster or one star; defines other helper functions not associated with TTauriStar class.
 */

#include "ttauristar.hpp"

using namespace std;
/**
 * \brief read cmk data file
 *
 * \param fname      string representing the cmk data filename
 * \returns          a vector of vectors of doubles
 */
vector<vector<double>> readcmk(string fname)
{
	ifstream inputFile(fname);

	if (!inputFile.good()) {
        throw invalid_argument( "Couldn't open cmk file for reading" );
    }

    // throw the first line
    string line;
    getline(inputFile, line);

    // allocate space on the heap to store the cmk data as a vector of vectors of doubles
    vector<vector<double>> cmktable;

    vector<double> cmkmasses;
    vector<double> cmkages;
    vector<double> cmkradii;

    const double PI = 3.141592653589793;
    const double SIGMAB = 5.67e-5;   // stefan-boltzmann constatn in erg/K^4 cm^2 s
    const double SOLARRADIUS = 7e10; // Solar radius in cm
    const double SOLARLUMINOSITY = 3.86e33;  // Solar luminosity in ergs per second
    // Solar tempertaure in Kelvins
    const double SOLARTEMPERATURE = pow(SOLARLUMINOSITY/(4*PI*pow(SOLARRADIUS,2)*SIGMAB),0.25);


    double mass; // Mass in M_sun
    double logage; //???
    double age; // age in Myrs
    double logL; // log(Luminosity)
    double logT; // log(Temperature)
    double radius;

    // read one line at a time
    while (getline(inputFile, line)){
    	stringstream lineStream(line);
        // read the four doubles separated by spaces
        // and put them in the corresponding vector
        lineStream >> mass;
        lineStream >> logage;
        lineStream >> logL;
        lineStream >> logT;
        // calculate age
        age = pow(10.0,logage - 6.0);
        // calculate radius
        radius = pow(10,logL/2) * pow(pow(10,logT)/SOLARTEMPERATURE,-2);
        cmkmasses.push_back(mass);
        cmkages.push_back(age);
        cmkradii.push_back(radius);
    }
    // store mass, age, radius vectors in cmktable
    cmktable.push_back(cmkmasses);
    cmktable.push_back(cmkages);
    cmktable.push_back(cmkradii);

    // close file
    inputFile.close();
    return cmktable;
}

/**
 * \brief read cluster mass and age file
 *
 * \param fname      string representing the cluster filename
 * \returns          a vector of vectors of doubles
 */
vector<vector<double>> readcluster(string fname)
{
    ifstream inputFile(fname);

    if (!inputFile.good()) {
        throw invalid_argument( "Couldn't open cluster file for reading" );
    }

    // throw the first line
    string line;
    getline(inputFile, line);

    // allocate space on the heap to store the cmk data as a vector of vectors of doubles
    vector<vector<double>> startable;

    vector<double> logmasses;
    vector<double> logages;

    double mass;
    double age;

    // read one line at a time
    while (getline(inputFile, line)){
        stringstream lineStream(line);
        // read the four doubles separated by spaces
        // and put them in the corresponding vector
        lineStream >> mass;
        lineStream >> age;
        // push back log10 of masses and ages
        logmasses.push_back(log10(mass));
        if (log10(age) < -3) {
            cout << "mass " << mass << "age " << age << endl;
        }
        logages.push_back(log10(age));
    }
    // store mass, age, radius vectors in cmktable
    startable.push_back(logmasses);
    startable.push_back(logages);

    // close file
    inputFile.close();
    return startable;
}

/**
 * \brief generate the distribution corresponding to a cluster
 *
 * \param fname      vectors of logmasses and logages
 * \param n          size of the simulated star table
 * \returns          a distribution
 */
vector<vector<double>> generatedistribution(vector<vector<double>> startable, double n)
{
    // create mass bins
    vector<double> logmassbins{-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1};
    // number of mass bins
    size_t nummassbin = logmassbins.size() - 1;
    // initialize weights with zeros
    vector<double> logmassweights(nummassbin,0.0);
    // create age bins
    vector<double> logagebins{-1.2,-1,-0.8,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
    /*
    double logagemin = -2.0;
    double logagemax = 3.0;
    for (size_t i = 0; i <= 10; ++i) {
        logagebins.push_back(logagemin+(logagemax-logagemin)*i/10);
    }
    */
    // number of age bins
    size_t numagebin = logagebins.size() - 1;
    // initilize ages with zeros
    vector<vector<double>> logageweights;
    for (size_t i = 0; i < nummassbin; ++i) {
        logageweights.push_back(vector<double>(numagebin,0.0));
    }
    // iterate through startable to get the weight
    // get logmass and logage from the table
    vector<double> logmasses = startable[0];
    vector<double> logages = startable[1];
    // the number of stars in the table
    size_t numstar = logmasses.size();
    // iterate over the stars in the cluster
    for (size_t i = 0; i < numstar; ++i) {
        // find the interval that the logmass belongs to
        for (size_t j = 0; j < nummassbin; ++j) {
            if (logmassbins[j] <= logmasses[i] && logmasses[i] < logmassbins[j+1]) {
                // increase the mass count in the appropriate interval
                logmassweights[j] += 1.0;
                // look at its age
                for (size_t k = 0; k < numagebin; ++k) {
                    if (logagebins[k] <= logages[i] && logages[i] < logagebins[k+1]) {
                        // increase the age weight
                        logageweights[j][k] += 1.0;
                    }
                }
            }
        }
    }
    /*
    for (size_t i = 0; i < nummassbin; ++i) {
        cout << "m" << logmassweights[i] << endl;
        for (size_t j = 0; j < numagebin; ++j) {
            cout << logageweights[i][j] << endl;
        }
    }
    */
    // create distributions
    piecewise_linear_distribution<> logmassdist(logmassbins.begin(),
        logmassbins.end(),logmassweights.begin());
    vector<piecewise_linear_distribution<>> logagedists;
    for (size_t i = 0; i < nummassbin; ++i) {
        logagedists.push_back(piecewise_linear_distribution<>(logagebins.begin(),
        logagebins.end(),logageweights[i].begin()));
    }
    // vectors to store simulated mass
    vector<double> simlogmasses;
    vector<double> simlogages;
    // random number generator
    random_device rd;
    mt19937 gen(rd());

    // choose n stars from the distribution
    for (size_t i = 0; i < n; ++i){
        double logmass = logmassdist(gen);
        double logage = 0.0;
        for (size_t i = 0; i < nummassbin; ++i) {
            if (logmassbins[i] <= logmass && logmass < logmassbins[i+1]) {
                logage = logagedists[i](gen);
            }
        }
        simlogmasses.push_back(logmass);
        simlogages.push_back(logage);
    }
    // pack the two vectors as a table
    vector<vector<double>> simstartable;
    simstartable.push_back(simlogmasses);
    simstartable.push_back(simlogages);
    return simstartable;
}

/**
 * \brief plot startable, a 2-D vector of star values
 *
 * \param startable  vector of <logmass,logage>
 * \returns
 */
void plotstartable(vector<vector<double>> startable, size_t cluster, bool simulated)
{
    FILE * temp = fopen("star.temp", "w");
    FILE* gp=popen("gnuplot -persistent","w");
    for(size_t k=0;k<startable[0].size();k++) {
        fprintf(temp,"%f %f \n",startable[0][k],startable[1][k]);
    }
    fprintf(gp, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    if (cluster == 1) {
        if (simulated == false) {
            fprintf(gp, "%s \n", "set output 'ONC.eps'");
            fprintf(gp, "%s%s %s \n", "set title \"","ONC","\"");
        } else {
            fprintf(gp, "%s \n", "set output 'simulatedONC.eps'");
            fprintf(gp, "%s%s %s \n", "set title \"","simulatedONC","\"");
        }
        fprintf(gp, "%s \n", "set yrange [-5:3]");
    } else {
        if (simulated == false) {
            fprintf(gp, "%s \n", "set output 'NGC2264.eps'");
            fprintf(gp, "%s%s %s \n", "set title \"","NGC2264","\"");
        } else {
            fprintf(gp, "%s \n", "set output 'simulatedNGC2264.eps'");
            fprintf(gp, "%s%s %s \n", "set title \"","simulatedNGC2264","\"");
        }
    }

    fprintf(gp, "%s \n", "set xlabel \"LogMass\"");
    fprintf(gp, "%s \n", "set ylabel \"LogAge\"");
    fprintf(gp, "%s \n", "plot 'star.temp'");

    fclose(temp);
}

/**
 * \brief plot period histogram,
 *
 * \param periods1 for stars with mass < 0.25 and periods2 for stars with mass > 0.25
 * \returns
 */
void plothistogram(vector<double> periods1, vector<double> periods2)
{
    FILE * temp1 = fopen("periods1.temp", "w");
    FILE * temp2 = fopen("periods2.temp", "w");
    FILE* gp3=popen("gnuplot -persistent","w");
    for(size_t k=0;k<periods1.size();k++) {
        fprintf(temp1,"%f \n",periods1[k]);
    }
    for(size_t k=0;k<periods2.size();k++) {
        fprintf(temp2,"%f \n",periods2[k]);
    }

    fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    fprintf(gp3, "%s \n", "set output 'distribution.eps'");
    fprintf(gp3, "%s\n", "binwidth=1");
    fprintf(gp3, "%s\n", "set boxwidth binwidth");
    fprintf(gp3, "%s\n", "bin(x,width)=width*floor(x/width) + binwidth/2.0");
    fprintf(gp3, "%s%s %s \n", "set multiplot layout 2,1 title \"","Period Distribution","\"");
    fprintf(gp3, "%s \n", "set ylabel \"Number of stars\"");
    fprintf(gp3, "%s \n", "unset xlabel");
    fprintf(gp3, "%s \n", "set xrange [0:14]");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m < 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s \n", "plot 'periods1.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes notitle");
    fprintf(gp3, "%s \n", "set xlabel");
    fprintf(gp3, "%s \n", "set xlabel \"Period (days)\"");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m > 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s \n", "plot 'periods2.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes notitle");
    fprintf(gp3, "%s \n", "unset multiplot");

    fclose(temp1);
    fclose(temp2);
}

/**
 * \brief plot observed ONC period distribution
 *
 * \returns
 */

void plotdistribution()
{
    // herbst (2002) data
    ifstream infile("herbst.txt");
    string line;
    // skip the first line
    getline(infile, line);
    double period, mass, disgard;
    vector<double> periods1, periods2;
    while (infile >> disgard >> period >> disgard >> disgard >> mass >> disgard) {
        if (period > 0.001) {
            if (mass < 0.25) {
                periods1.push_back(period);
            } else {
                periods2.push_back(period);
            }
        }
    }
    // plot
    plothistogram(periods1, periods2);
}

/**
 * \brief simulation of n stars
 *
 * \param cluster 1=ONC, 2=NGC
 * \returns
 */
void simulation(vector<vector<double>> simstartable)
{
    // compute cmk table
    vector<vector<double>> cmktable = readcmk("cmkdata.txt");

    // generator
    default_random_engine genm;
    default_random_engine genb;
    // normal distribution of log(massfactor)
    normal_distribution<double> logmdotfactordist(0.0,0.32);
    // normal distribution of bfieldstrength
    normal_distribution<double> bfielddist(1.67,0.1);

    // loop through the stars
    std::clock_t start;
    double duration;
    start = std::clock();
    FILE * datafile;
    datafile = fopen("simulationONC.txt", "w");
    // write first line
    fprintf(datafile,"%s \n","mass(solarmass) age(Myr)        period(days)    log(mdotfactor) bfield(kG)");
    vector<double> periods1, periods2;
    for (size_t i = 0; i < simstartable[0].size(); ++i) {
        double mass = pow(10,simstartable[0][i]);
        double age = pow(10,simstartable[1][i]);
        double logmdotfactor = 0; //logmdotfactordist(genm); //pick fixed value
        double bfieldstrength = 1.67; //bfielddist(genb); //pick fixed value
        TTauriStar star = TTauriStar(cmktable, mass, age, pow(10,logmdotfactor), bfieldstrength);
        double period = star.update();
        // write to file
        fprintf(datafile,"%f        %f        %f        %f        %f \n",mass,age,period,logmdotfactor,bfieldstrength);
        if (period > 0.001) {
            if (mass < 0.25) {
                periods1.push_back(period);
            } else {
                periods2.push_back(period);
            }
            //star.plot(1,2);
        }
    }
    fclose(datafile);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"the program takes "<< duration << " s" << endl;
    // plot
    plothistogram(periods1, periods2);
}




/*The main program
 */
int main()
{
		// choose cluster simulation of single-star simulation
    size_t typeofsimulation;
    cout << "Type 1 for single-star simulation, type 2 for cluster simulation" << endl;
    cin >> typeofsimulation;
    cout << "The value you entered is " << typeofsimulation << endl;

    if (typeofsimulation == 1) {
        // mass
        double mass;
        cout << "mass: ";
        cin >> mass;
        cout << "The mass is " << mass << endl;

        // age
        double age; //Myrs
        cout << "age: ";
        cin >> age;
        cout << "The age is " << age << endl;

        // compute cmk table
        vector<vector<double>> cmktable = readcmk("cmkdata.txt");
				
        // create a star
        TTauriStar star = TTauriStar(cmktable, mass, age, 1, 1.67);
        star.update();

<<<<<<< HEAD
        int TIME_INDEX = 1;
        int N_Y_INDEX = 10;
=======
>>>>>>> 32ab899eeb927028c033fa342b2328a902b90dfe

				// Saving plots to folder
        stringstream folderNameStream;
        folderNameStream << std::fixed << setprecision(2);
        folderNameStream << "_Mass_" << mass;
        //folderNameStream << "Star: Mass_" << mass << " Age_ " << age;
        string folderName = folderNameStream.str();

				// make new folder and change directories into new folder
        mkdir(folderName.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
        chdir(folderName.c_str());

			  // Plot all of the star's parameters as a function of time
				int TIME_INDEX = 1;
        int N_Y_INDEX = 9;
        for(int y_index = 2; y_index <= N_Y_INDEX; y_index += 1) {
            star.plot(TIME_INDEX, y_index);
        }


    } else if (typeofsimulation == 2) {
        // Simuate a cluster: 1 = ONC, 2 = NGC 2264
        size_t cluster;
        cout << "Which cluster do you want to simulate?" << endl;
        cout << "Type 1 for ONC, type 2 for NGC2264" << endl;
        cin >> cluster;
        cout << "The value you entered is " << cluster << endl;
        // sample from a given cluster
        vector<vector<double>> startable;
        if (cluster == 1) {
            startable = readcluster("hillenbrand.txt");
        } else if (cluster == 2) {
            startable = readcluster("dahm.txt");
        } else {
            throw invalid_argument( "Invalid input" );
        }
        // plotstartable(startable, cluster, false);
        // simulate a startable of size n
        vector<vector<double>> simstartable = generatedistribution(startable,1000);
        // plotstartable(simstartable, cluster, true);
        simulation(simstartable);
    } else {
        throw invalid_argument( "Invalid input" );
    }

    // Plot the observed period distribution for ONC
    // plotdistribution();

    return 0;
}
