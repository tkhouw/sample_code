#include <iostream>
#include <array>
#include <vector>
#include <cstdlib> //for random numbers
#include <cmath> //for pow function
#include <fstream> // for output files
#include <string>
#include <numeric> // for accumulate
#include <algorithm> //for max_element finder function
#include <ctime> //for seeding RNG
#include <math.h> //for exp and log functions
#include <chrono> //for ms timestamp for seed

struct Agent {

    double wealth;
    Agent(double w): wealth(w), wealthAverage(w) {
        fitness = 1.0;
        stepSizeM = 0;
    }
    Agent(): wealth(0.0), wealthAverage(0.0) {
        fitness = 1.0;
        stepSizeM = 0;
    }
    double wealthAverage; //historical average of wealth value for each agent
    double wealthSquaredAverage; //historical average of squared wealth value for each agent
    double fitness;
    double guess;
    int stepSizeM;

};

struct Network {
    
    //Parameters
    double mu; //growth parameter
    double lambda; //distribution parameter delta w \propto w^lambda
    double f; //the fraction of wealth transferred in an exchange
    double p; //one half of the net benefit accrued to both agents after a trade expressed as a fraction of the transferred wealth
    int N; //the total number of agents in the system
    int time; //the number of time steps that this network has experienced so far
    double totalAbsoluteWealth; //all of the wealth ever added to a system from a growth spurt / N
    double growthNoise; //multiplicative noise to Kang rule:  wi += (1 + X)*(wi^lambda / Sum wi), X in [-growthNoise, growthNoise]
    double additiveGrowthNoise; //noise to add to Kang rule:  wi += X + (wi^lambda / Sum wi), X in [-addedGrowthNoise, addedGrowthNoise]

    bool doFitness; //whether we use agent fitnesses in the simulation
    bool doExchangeGrowth; //whether to use pairwise lambda ghrowth based on exchanges

    double truth = 0.5;
    bool doLearn;
    double learnStep; // > 1.0 for convergence

    double nukeLimit; // maxWealth before you are immune to lambda
    int seed;

    std::vector<std::vector<double>> matrix;
    std::vector<std::vector<int>> nonZeros; //will contain indices of nonzero elements in matrix
    std::vector<Agent> agents;
    Network(int n, int sead) { //N is number of agents in network; p is the net benefit parameter

	N = n;
	seed = sead;
	//The first thing we need to do is seed the RNG...
	srand(seed);
	
	//Default values of parameters (can be changed in main)
	mu = 0.1;
	lambda = 1.0;
	f = 0.1;
	p = 0.0;
	time = 0;
        totalAbsoluteWealth = (double)N;
        growthNoise = 0.0;
        additiveGrowthNoise = 0.0;
        doFitness = false;
        doExchangeGrowth = false;
        learnStep = 2.0;
        nukeLimit = (double)N;	

	fullyConnectMatrix();
        checkNonzeroConnections();

	Agent initialAgent(1.0);
	agents.assign(N, initialAgent);

	randomizeWealths(); //randomizing wealths while keeping a total of N
	randomizeGuesses();
        //equalizeWealths();

	for (int a=0; a<N; a++) {//for each agent
	    agents.at(a).wealthAverage = agents.at(a).wealth;
	    agents.at(a).wealthSquaredAverage = pow(agents.at(a).wealth, 2);
	}

    }

    void growthSpurt(); //uses parameters mu and lambda
    void exchange(int agentIndex1, int agentIndex2);
    void writeAssetsToFile(std::string filename);
    void normalizeWealth(); //makes the total wealth sum to N
    void exchangeStep(int n); //evolves the sytem n exchange steps.  1 exchange step is N exchanges
    void fullyConnectMatrix(); //Makes the connection matrix all 1+p except no self-connections
    void BarabasiAlbertMatrix(int halfk); //Makes the connection matrix randomly scale-free with avg degree k
    void writeDegreeDistribution(std::string filename); //Writes degree distribution of network to file
    void writeMatrixToFile(std::string filename);
    void loadMatrixFromFile(std::string filename);
    void erdosRenyiMatrix(double r); //r is the probability that any two nodes are connected
    double wealthMetric(bool ignoreRichest=false);
    void resetWealthAverages(); //fills wealth averages vector with the current wealth of each agent
    void updateWealthAverages();  //uses recursive formula to update the historical wealth average of each agent
    int richestAgent(); //returns the index of the wealthiest agent
    double wealthVariance(bool ignoreRichest); //returns the susceptibility chi (X) from the historical wealth averages of each agent if ignoreRichest=True.  Otherwise it's just the variance.
    double energy(); //gives the current energy (-1 + (1/N)*S)
    void randomizeWealths(); //randomizing wealths while keeping a total of N
    double checkDistance2D(int i, int j);
    int checkX(int i);
    int checkY(int i);
    void make2DLattice(double R);
    void destroyMatrix();  //sets all matrix elements to zero
    double phi();  //returns order parameter (wealth of all but richest agent)
    std::vector<int> getClusterSizes(); //Returns a list of the size of each cluster in the network
    void recursiveNeighborFunction(int a, int& clusterSize, bool* visited); //trippy recursive func. for getClusterSize
    void writeClusterSizes(std::string filename); //outputs cluster sizes to file
    std::vector<int> getLargestCluster(); //Returns a list of each agent in the largest connected cluster in the network
    void removeAgents(std::vector<int> agentsToRemove); //It will be as if these agents never existed 
    void scaleFreeMatrix(double exponent); //Uses Chung-Lu method to construct a general power-law degreely distributed matrix: P(k) ~ k^-exponent
    void prune(); //removes agents that are not part of the largest connected cluster
    void preferentialAttachmentMatrix(double lambda); //Uses B-A algorithm to construct scale-free matrices if lambda=1
    void scaleFreeConfigMatrix(double exponent, int minDegree); //uses configuration model to generate power-law graph, and then cleans it up
    void randomizeFitnessesGaussian(); //Assigns random fitness to each agent in the network
    void breakAndBind(double breakRate, int minDegree); //evolves matrix depending on wealths of each agent
    void growthSpurt(double mew); //uses parameter lambda and custom mew!
    void trickleDownGrowthSpurt(double intensity);
    void equalityDrivenGrowthSpurt();
    void inequalityDrivenGrowthSpurt();
    void inequalityWeakenedGrowthSpurt(double intensity); //growth becomes very negative if inequality is too high
    int poorestAgent(); //returns index of agent with least wealth
    double totalWealth(); //returns the sum of all agents' wealth
    double getDegree(int agentIndex); //returns degree of specified node.  Should be an int if p=0
    double averageDegree(); //returns the average degree of all nodes in the graph
    void normalizeFitnesses();  //rescales all fitnesses to N.  This does not change the dynamics of the system.
    void makeHub(int index, double probability); //turns a single node into a hub
    void hubifyMatrix(double hubProbability, double hubStrength); //Adds hubs to a matrix using makeHub
    int getXNeighbor(int homeIndex, int distance);
    int getYNeighbor(int homeIndex, int distance);
    int getIndex(int x, int y); //returns index of agent at (x,y)
    std::vector<double> correlations2D(int maxLength); //returns average correlation function at this moment in time 
    void writeCorrelations(std::string filename, int maxLength, int dim);
    void makeRing(double R); //Creates 1D PBC topology
    int checkDistance1D(int i, int j); //gives distance between agents i and j on ring
    std::vector<double> correlations1D(int maxLength); //returns 1D correlation as function of separation
    int checkX3D(int i);
    int checkY3D(int i);
    int checkZ3D(int i);
    int getXNeighbor3D(int homeIndex, int distance);
    int getYNeighbor3D(int homeIndex, int distance);
    int getZNeighbor3D(int homeIndex, int distance);
    int getIndex(int x, int y, int z); //returns index of agent at (x,y,z)
    std::vector<double> correlations3D(int maxLength); //returns average correlation function at this moment in time 
    void makeLattice(int dim, double R); //constructs lattice of either 2 or 3 dimensions
    double checkDistance3D(int i, int j); //3D distance between two agents
    void makeAnalogLattice(int dim, double correlationLength); //constructs lattice of either 2 or 3 dimensions
    void equalizeWealths(); //sets all wealths equal to 1.0
    void checkNonzeroConnections(); //updates and overwrites the list of connected agents
    void exchangeStepNonzeroOnly(int n); //perfors n exchange steps but it avoids all exchanges with matrix element zero
    void trickleDownEnergyGrowthSpurt(); //uses energy as mu value
    void updateGuess(int agentIndex, bool justWon); //changes an agent's guess randomly depending on whether it just won
    void writeGuessesToFile(std::string filename);
    void writeStepSizesToFile(std::string filename);
    void randomizeGuesses(); //randomizing eah agent's guess
    void loadWealthsFromFile(std::string filename);  //Sets agent wealths according to file.  Must have right N
    double dependentLambda(double wealth); //returns lambda value for an individual agent based on their wealth

};


				    
