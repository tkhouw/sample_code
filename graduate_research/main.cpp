#include "exchange.h"
#include <sstream>

//Conducts simulations using exchange.cpp

int main() {

    //First we read in variables
    //The variables are described in read.in
    std::ifstream inFile;
    std::string line;
    double floaty;

    inFile.open("read.in");
    if (!inFile) {
	std::cerr << "Unable to open read.in";
	exit(1);
    }

    double p;
    std::getline(inFile,line);
    std::istringstream(line) >> p;
    std::cout << p << std::endl;

    double lambda;
    std::getline(inFile,line);
    std::istringstream(line) >> lambda;
    std::cout << lambda << std::endl;

    int T;
    std::getline(inFile,line);
    std::istringstream(line) >> T;
    std::cout << T << std::endl;

    bool sf;
    std::getline(inFile,line);
    std::istringstream(line) >> sf;
    std::cout << sf << std::endl;

    double a;
    std::getline(inFile,line);
    std::istringstream(line) >> a;
    std::cout << a << std::endl;

    bool er;
    std::getline(inFile,line);
    std::istringstream(line) >> er;
    std::cout << er << std::endl;

    double r;
    std::getline(inFile,line);
    std::istringstream(line) >> r;
    std::cout << r << std::endl;

    int N;
    std::getline(inFile,line);
    std::istringstream(line) >> N;
    std::cout << N << std::endl;

    double mu;
    std::getline(inFile,line);
    std::istringstream(line) >> mu;
    std::cout << mu << std::endl;

    double breakRate;
    std::getline(inFile,line);
    std::istringstream(line) >> breakRate;
    std::cout << breakRate << std::endl;

    double minDegree;
    std::getline(inFile,line);
    std::istringstream(line) >> minDegree;
    std::cout << minDegree << std::endl;

    double intensity;
    std::getline(inFile,line);
    std::istringstream(line) >> intensity;
    std::cout << intensity << std::endl;

    double noise;
    std::getline(inFile,line);
    std::istringstream(line) >> noise;
    std::cout << noise << std::endl;

    double additiveNoise;
    std::getline(inFile,line);
    std::istringstream(line) >> additiveNoise;
    std::cout << additiveNoise << std::endl;

    bool doFitness;
    std::getline(inFile,line);
    std::istringstream(line) >> doFitness;
    std::cout << doFitness << std::endl;

    bool doExchangeGrowth;
    std::getline(inFile,line);
    std::istringstream(line) >> doExchangeGrowth;
    std::cout << doExchangeGrowth << std::endl;

    bool hubify;
    std::getline(inFile,line);
    std::istringstream(line) >> hubify;
    std::cout << hubify << std::endl;

    int correlationMaxLength;
    std::getline(inFile,line);
    std::istringstream(line) >> correlationMaxLength;
    std::cout << correlationMaxLength << std::endl;

    bool sq;
    std::getline(inFile,line);
    std::istringstream(line) >> sq;
    std::cout << sq << std::endl;

    double rsq;
    std::getline(inFile,line);
    std::istringstream(line) >> rsq;
    std::cout << rsq << std::endl;

    double f;
    std::getline(inFile,line);
    std::istringstream(line) >> f;
    std::cout << f << std::endl;

    bool nonzeroOnly;
    std::getline(inFile,line);
    std::istringstream(line) >> nonzeroOnly;
    std::cout << nonzeroOnly << std::endl;

    bool ring;
    std::getline(inFile,line);
    std::istringstream(line) >> ring;
    std::cout << ring << std::endl;

    double rring;
    std::getline(inFile,line);
    std::istringstream(line) >> rring;
    std::cout << rring << std::endl;

    bool sfc;
    std::getline(inFile,line);
    std::istringstream(line) >> sfc;
    std::cout << sfc << std::endl;

    bool tde;
    std::getline(inFile,line);
    std::istringstream(line) >> tde;
    std::cout << tde << std::endl;

    int halfk;
    std::getline(inFile,line);
    std::istringstream(line) >> halfk;
    std::cout << halfk << std::endl;

    bool doLearn;
    std::getline(inFile,line);
    std::istringstream(line) >> doLearn;
    std::cout << doLearn << std::endl;

    double learnStep;
    std::getline(inFile,line);
    std::istringstream(line) >> learnStep;
    std::cout << learnStep << std::endl;

    double nukeLimit;
    std::getline(inFile,line);
    std::istringstream(line) >> nukeLimit;
    std::cout << nukeLimit << std::endl;

    std::string loadWealths;
    std::getline(inFile,line);
    std::istringstream(line) >> loadWealths;
    std::cout << loadWealths << std::endl;

    int seed;
    std::getline(inFile,line);
    std::istringstream(line) >> seed;
    std::cout << seed << std::endl;

    int nData;
    std::getline(inFile,line);
    std::istringstream(line) >> nData;
    std::cout << nData << std::endl;

    bool moreExchanges;
    std::getline(inFile,line);
    std::istringstream(line) >> moreExchanges;
    std::cout << moreExchanges << std::endl;
    int exchangeSteps;
    if (moreExchanges) {
	exchangeSteps = N;
    }
    else {
	exchangeSteps = 1;
    }

    int corrDim;
    if (sq) {
	corrDim = 2;
    }
    else {
	corrDim = 1;
    }

    if ((correlationMaxLength > 0) && (!(sq) && !(ring))) {
	std::cout << "Warning: You shouldn't use the correlation function for non-lattice topologies" << std::endl;
    }
    if (!seed) {    
	seed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	std::cout << "In here seed is " << seed << std::endl;
        //std::time_t t = std::time(0);  // t is an integer type
	//seed = t;
    }

    std::string loadMatrix = "";///mnt/d/runs/conff2.2min2N10/0/data/matrix.txt";
    std::string dataDir = "data/";

    Network net(N, seed);
    net.lambda = lambda;
    net.p = p;
    net.mu = mu;
    net.f = f;
    net.growthNoise = noise;
    net.additiveGrowthNoise = additiveNoise;
    net.doFitness = doFitness;
    net.doExchangeGrowth = doExchangeGrowth;
    net.doLearn = doLearn;
    net.learnStep = learnStep;
    net.nukeLimit = nukeLimit;
    std::cout << "Seed: " << net.seed << std::flush;
    std::string seedFilename = dataDir + "seed.txt";
    std::ofstream seedFile(seedFilename);
    seedFile << net.seed;
    seedFile.close();
    net.fullyConnectMatrix();
    //net.makeAnalogLattice(1, 2000);
    if (loadMatrix.empty()) {
        if (sf) {
            std::cout << "Generating scale-free network (preferential attachment)..." << std::flush;
            net.BarabasiAlbertMatrix(halfk);

            std::cout << "removing bad nodes..." << std::flush;
	    net.prune();
    
        }
        else if (er) {
            std::cout << "Generating Erdos-Renyi network..." << std::flush;
            net.erdosRenyiMatrix(r);

            std::cout << "removing bad nodes..." << std::flush;
	    net.prune();

        }
	else if (sq) {
            std::cout << "Generating 2D lattice..." << std::flush;
            net.make2DLattice(rsq);
	}
	else if (ring) {
	    std::cout << "Generating 1D ring..." << std::flush;
	    net.makeRing(rring);
        }
        if (sfc) {
            std::cout << "Generating scale-free network (configuration model)..." << std::flush;
            net.scaleFreeConfigMatrix(a, 2);

            std::cout << "removing bad nodes..." << std::flush;
	    net.prune();
    
        }

    }
    else { //if the user wants to load a previously-made matrix
        std::cout << "Loading matrix from file " << loadMatrix << "... " << std::endl;
	net.loadMatrixFromFile(loadMatrix);

        std::cout << "removing bad nodes..." << std::flush;
        net.prune();
    }
    if (!(loadWealths=="none")) {
	std::cout << "Loading wealths from file " << loadWealths << std::endl;
	net.loadWealthsFromFile(loadWealths);
    }

    std::cout << net.N << std::endl;

    if (hubify) {
        std::cout << "Hubifying matrix..." << std::flush;
        net.hubifyMatrix(0.01, 0.9);
    }

    if (true) {
	//Writing degree distribution to file
	std::string degreeFilename = dataDir + "degrees.txt";
        net.writeDegreeDistribution(degreeFilename);
	
	//Writing connection matrix to file
	std::string matrixFilename = dataDir + "matrix.txt";
        //net.writeMatrixToFile(matrixFilename);

	std::cout << " Done!" << std::endl;
    }

    //We run the system for T time steps before taking data
    std::cout << "Running transient..." << std::flush;
    int transientTime = T;
    if (!(loadWealths=="none")) {
	transientTime = 0;
        std::cout << "Never mind! Skipping transient!" << std::flush;
    }
    for (int t=0; t<transientTime; t++) {
	if (!nonzeroOnly) {
            net.exchangeStep(exchangeSteps);
	}
	else {
	    net.exchangeStepNonzeroOnly(exchangeSteps);
	}
	if (!tde) {
            net.growthSpurt();
	}
	else {
            net.trickleDownEnergyGrowthSpurt();
	}
        //net.inequalityWeakenedGrowthSpurt(intensity);
	net.normalizeWealth();
        if (doFitness) {
            net.normalizeFitnesses();
        }
	net.updateWealthAverages();
        //Special breaking and rebinding step:
        if ((!(t%1000)) && (breakRate > 0.0)) {
            net.breakAndBind(breakRate, minDegree);
        }
	if (doLearn) {
	    if (!(t % (T/5000))) {
                int block = (t / (T/5000)) - 5000;
                std::string b = std::to_string(block);
	        std::string filename;
                filename = dataDir + "guesses_" + b + ".txt";
                std::cout << "Saving guess distribution to file " << filename << std::endl;
                net.writeGuessesToFile(filename);
		std::string stepSizesFilename = dataDir + "stepSizes_" + b + ".txt";
                net.writeStepSizesToFile(stepSizesFilename);
	    }
	}
    }
    std::cout << " Done!" << std::endl;
    std::string wealthMetricFilename = dataDir + "wealthMetric.txt";
    std::cout << "Opening output file " << wealthMetricFilename << std::endl;
    std::ofstream wealthMetricFile(wealthMetricFilename);
    std::string susceptibilityFilename = dataDir + "susceptibility.txt";
    std::cout << "Opening output file " << susceptibilityFilename << std::endl;
    std::ofstream susceptibilityFile(susceptibilityFilename);
    std::string energyFilename = dataDir + "energy.txt";
    std::cout << "Opening output file " << energyFilename << std::endl;
    std::ofstream energyFile(energyFilename);
    std::string phiFilename = dataDir + "phi.txt";
    std::cout << "Opening output file " << phiFilename << std::endl;
    std::ofstream phiFile(phiFilename);
    std::string clusterSizesFilename = dataDir + "clusterSizes.txt";
    std::cout << "Taking data now" << std::endl;
    net.writeClusterSizes(clusterSizesFilename);
    net.time = 0; //resetting time so as to not mess up wealth metric
    net.resetWealthAverages(); //we also need to reset wealth averages
    for (int t=0; t<T+1; t++) {

        //simulation part
	if (!nonzeroOnly) {
            net.exchangeStep(exchangeSteps);
	}
	else {
	    net.exchangeStepNonzeroOnly(exchangeSteps);
	}
	if (!tde) {
            net.growthSpurt();
	}
	else {
            net.trickleDownEnergyGrowthSpurt();
	}
        //net.inequalityWeakenedGrowthSpurt(intensity);
	net.normalizeWealth();
        if (doFitness) {
            net.normalizeFitnesses();
        }
	net.updateWealthAverages();

        //recording output part
        wealthMetricFile << net.wealthMetric() << std::endl;
        susceptibilityFile << net.wealthVariance(false) << std::endl;
        energyFile << net.energy() << std::endl;
        phiFile << net.phi() << std::endl;

        //recording wealth distribution part
	if (!(t % (T/nData))) {
	    int block = t / (T/nData);
	    std::string b = std::to_string(block);
	    std::string filename;
            filename = dataDir + "wealthDistribution-" + b + ".txt";
            std::cout << "Saving wealth distribution to file " << filename << std::endl;
            net.writeAssetsToFile(filename);
            if (correlationMaxLength) {
                std::string corrFilename = dataDir + "correlations-" + b + ".txt";
                net.writeCorrelations(corrFilename, correlationMaxLength, corrDim);
            }
            if (doLearn) {
		std::string guessFilename = dataDir + "guesses_" + b + ".txt";
                net.writeGuessesToFile(guessFilename);
		std::string stepSizesFilename = dataDir + "stepSizes_" + b + ".txt";
                net.writeStepSizesToFile(stepSizesFilename);
	    }

	}
        //Special breaking and rebinding step:
        if ((!(t%10000)) && (breakRate > 0.0)) {
            net.breakAndBind(breakRate, minDegree);
            if (!(t%10000)) {            
	        int bock = t / (T/10000);
                std::string B = std::to_string(bock);
                std::string degreeFilename = dataDir + "degrees-" + B + ".txt";
                net.writeDegreeDistribution(degreeFilename);
            }
        }

    }
    wealthMetricFile.close();
    susceptibilityFile.close();
    energyFile.close();
    phiFile.close();
    std::cout << "Finished!" << std::endl;

    return 0;
}

