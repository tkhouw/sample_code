#include <iostream>
#include "exchange.h"


void Network::growthSpurt() {
//introduces new wealth into the system and distributes among agents

    double denominator = 0.0; //running sum
    double totalWealth = 0.0;
    for (int a = 0; a<N; a++) {
	double wealth = agents.at(a).wealth;
	denominator += pow(wealth, dependentLambda(wealth));
	totalWealth += wealth;
    }

    double deltaW = mu*totalWealth;
    double X = 0.0; //multiplicative noise term
    double X2 = 0.0; //additive noise term
    if (growthNoise != 0.0) {
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
        X = 2.0*growthNoise*(x - 0.5);
    }   
    if (additiveGrowthNoise != 0.0) {
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
        X2 = 2.0*additiveGrowthNoise*(x - 0.5);
    }   
    for (int a=0; a<N; a++) {
	double numerator;
        numerator = pow(agents.at(a).wealth, dependentLambda(agents.at(a).wealth));
	agents.at(a).wealth += X2 + (1.0 + X)*(deltaW*numerator / denominator);
    }

}

void Network::normalizeWealth() {
//Rescales so total system wealth becomes N
    double totalWealth = 0.0;
    for (int a=0; a<N; a++) { 
	totalWealth = totalWealth + agents.at(a).wealth;
    }
    for (int a=0; a<N; a++) { 
	agents.at(a).wealth = agents.at(a).wealth*(N / totalWealth);
    }

}


void Network::exchange(int agentIndex1, int agentIndex2) {

    //double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
    double exchangeAmount = f*std::min(agents.at(agentIndex1).wealth, agents.at(agentIndex2).wealth);

    int winner;
    int loser;
    
    if (!doLearn) {
        //Now we roll the die to see who wins the money
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
    
        //Fitness time! We will determine the probability that agent1 wins
        double probability = 0.5;
        if (doFitness) {
            probability = agents.at(agentIndex1).fitness / (agents.at(agentIndex1).fitness + agents.at(agentIndex2).fitness);
        }
        
        if (x > probability) {
    	winner = agentIndex1;
    	loser = agentIndex2;
        }
        else {
    	winner = agentIndex2;
    	loser = agentIndex1;
        }
    }
    else {
	double guessDistance1 = std::abs(agents.at(agentIndex1).guess - truth);
	double guessDistance2 = std::abs(agents.at(agentIndex2).guess - truth);
	if (guessDistance1 < guessDistance2) {
	    winner = agentIndex1;
	    loser = agentIndex2;
	}
	else {
	    winner = agentIndex2;
	    loser = agentIndex1;
	}

	//Now we need to update the guesses of the loser and winner
	updateGuess(winner, 1);
	updateGuess(loser, 0);
    }

    double originalLoserWealth = agents.at(loser).wealth;

    //the matrix element is 1+p or 0
    double bonus = std::max(exchangeAmount*matrix.at(loser).at(winner) - exchangeAmount, 0.0);

    //Wealth exchange
    agents.at(winner).wealth += std::min(matrix.at(loser).at(winner), 1.0)*exchangeAmount + bonus;
    agents.at(loser).wealth -= std::min(matrix.at(loser).at(winner), 1.0)*exchangeAmount - bonus;

    if (doFitness) {
    //Now we find the relative difference between their original wealths and their new wealths.  This will determine how the loser's fitness changes.
        double difference = std::max(originalLoserWealth - agents.at(winner).wealth, 0.0);
        double lossFactor = difference / originalLoserWealth;
        agents.at(loser).fitness *= (1.0 + lossFactor); //makes them more fit!
    }

    if (doExchangeGrowth) {
        double growthAmount = (mu/2.0)*(agents.at(winner).wealth + agents.at(loser).wealth);
        double winnerToLambda = pow(agents.at(winner).wealth, lambda);
        double loserToLambda = pow(agents.at(loser).wealth, lambda);
        agents.at(winner).wealth += growthAmount*winnerToLambda / (winnerToLambda + loserToLambda);
        agents.at(loser).wealth += growthAmount*loserToLambda / (winnerToLambda + loserToLambda);
        //Now we need to make sure we keep total wealth constant:
        agents.at(winner).wealth = agents.at(winner).wealth / (1.0 + (mu/2.0));
        agents.at(loser).wealth = agents.at(loser).wealth / (1.0 + (mu/2.0));
    }

}


void Network::exchangeStep(int n) {
    
    for (int i=0; i<n; i++) {
	//For each step we will do N exchanges
	for (int x=0; x<N; x++) {
	    //Choosing a random pair to exchange
	    int index1 = rand() % N;
	    int index2 = N; //initializing to bad value
	    while ((index2 == index1) || (index2 == N)) {
		index2 = rand() % N;
	    }
	    
            exchange(index1, index2);

	}
	time += 1;
    }
    
}

void Network::writeAssetsToFile(std::string filename) {

    std::ofstream assetFile(filename);
    for (int a=0; a<N; a++) {
	assetFile << agents.at(a).wealth << std::endl;
    }
    assetFile.close();

}

void Network::fullyConnectMatrix() {

    matrix.clear();
    //Default weight matrix is fully connected with weight one
    for (int i=0; i<N; i++) {
        std::vector<double> lilVector(N, 1.0 + p);
    	lilVector.at(i) = 0.0; //no self connections
    	matrix.push_back(lilVector);
    }

    checkNonzeroConnections();

}

void Network::BarabasiAlbertMatrix(int halfk) { //a is a paramter that will determine the exponent of the degree distribution

    matrix.clear();
    //First we will make a little string as the adjacency matrix.  Then we will fill the rest by preferential attachment
    std::vector<double> firstRow(N, 0.0);
    std::vector<double> secondRow(N, 0.0);
    //We will connect these two nodes:
    firstRow.at(1) = 1.0 + p;
    secondRow.at(0) = 1.0 + p;
    matrix.push_back(firstRow);
    matrix.push_back(secondRow);
    for (int i=2; i<halfk; i++) {
        std::vector<double> initRow(N, 0.0);
        initRow.at(i-1) = 1.0 + p;
	matrix.at(i-1).at(i) = initRow.at(i-1);
        matrix.push_back(initRow);
    }

    //Now we will do preferential attachment
    for (int i=std::max(2, halfk); i<N; i++) {
	std::vector<double> row(N, 0.0);

        double totalDegree = 0.0;
	for (int j=0; j<(int)matrix.size(); j++) {
	    totalDegree += std::accumulate(matrix.at(j).begin(), matrix.at(j).end(), 0);
	}

	while (std::accumulate(row.begin(), row.end(), 0) < halfk) {
        //while this row is still too disconnected from all other nodes
            double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
	    double runningSumTrick = 0.0;
	    double degreesSoFar = 0.0;
	    for (int j=0; j<(int)matrix.size(); j++) {
		double thisDegree = std::accumulate(matrix.at(j).begin(), matrix.at(j).end(), 0);
		double P = (thisDegree/(totalDegree - degreesSoFar));
		if (row.at(j) == 0.0) { //if we haven't decided to connect to this node yet
		    runningSumTrick += P;
		    if (x < runningSumTrick) {
		        row.at(j) = 1.0 + p;
			degreesSoFar += thisDegree;
		        break;
		    }
		}
	    }
	}

        //Now we make sure that the matrix elements of the nodes we already made reflect the new connections that we just made
	for (int l=0; l<(int)matrix.size(); l++) {
	    matrix.at(l).at(i) = row.at(l);
	}
	
	//Finally we add this new row to the matrix!
	matrix.push_back(row);
    }

    checkNonzeroConnections();

}//Finished??


void Network::writeDegreeDistribution(std::string filename) {
//Assumes connections in matrix are either 1 or 0

    std::ofstream File(filename);
    for (int j; j < N; j++) { //for each row in the adjacency matrix
	File << std::accumulate(matrix.at(j).begin(), matrix.at(j).end(), 0) << std::endl;
    }
    File.close();

}

void Network::writeMatrixToFile(std::string filename) {

    std::ofstream File(filename);
    for (int i=0; i<N; i++) {
	for (int j=0; j<N; j++) {
	    File << matrix.at(i).at(j) << std::endl;
	}
    }
    File.close();

}

void Network::loadMatrixFromFile(std::string filename) {

    //First we obliterate the matrix we already have
    matrix.clear();

    //Reading in file now
    std::ifstream File(filename);
    std::vector<double> elements; //will hold some matrix elements from file
    double c; //c is the connection between two nodes (an element in the adjacency matrix)
    int n = 0; //The number of matrix elements in the file you read
    while (File >> c) {
	elements.push_back(c);
	if ((int)elements.size() == N) {//if we now have an entire row of the matrix stored in elements
	    matrix.push_back(elements);
	    elements.clear(); //Making room for the next row now
	}
	n += 1;
    }

    //Checking that you read an appropriate file
    if (n != (N*N)) {
	std::cout << "Warning: wrong number of matrix elements in " << filename << ".  Expected " << (N*N) << " for a network of size " << N << "; found " << n << "." << std::endl;
    }
    //Checking the shape of the matrix
    bool goodShape = ((int)matrix.size() == N);
    bool zeroDiagonal = true; //assume at first that the diagonal is all zero
    bool squareShape = true; //assume true at first
    for (int i=0; i<(int)matrix.size(); i++) {
	//Making sure each row has N elements
	goodShape = (goodShape && ((int)matrix.at(i).size() == N));
	//Making sure the diagonal elements are zero
	zeroDiagonal = (zeroDiagonal && (matrix.at(i).at(i) == 0.0));
        //Making sure each row has as many elements as the nu,ber of rows
	squareShape = (squareShape && ((int)matrix.at(i).size() == (int)matrix.size()));
    }
    if (!goodShape) {
	std::cout << "Warning: the shape of the matrix loaded from " << filename << " is messed up.  Expected a square " << N << " by " << N << " matrix." << std::endl;
        if (squareShape) {
            std::cout << "Found square " << (int)matrix.size() << "x" << (int)matrix.size() << " matrix instead." << std::endl;
        }
        else {
            std::cout << "The loaded matrix isn't even square!" << std::endl;
        }
    }
    if (!zeroDiagonal) {
	std::cout << "Warning: non-zero diagonal elements are present in the matrix loaded from " << filename << ".  This means an agent will be able to trade with himself." << std::endl;
    }

    if ((n == (N*N)) && (goodShape) && (zeroDiagonal)) {
	std::cout << "Successfully loaded " << N << "x" << N << " adjacency matrix from " << filename << "!" << std::endl;
    }

    checkNonzeroConnections();

}

void Network::erdosRenyiMatrix(double r) {

    matrix.clear();
    for (int i=0; i<N; i++) {
	std::vector<double> lilVector(N, 0.0);
	matrix.push_back(lilVector);
    }
    for (int i=0; i<N; i++) {
	for (int j=0; j<i; j++) {
            double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
	    if (x < r) {
	        matrix.at(i).at(j) = 1.0 + p;
	        matrix.at(j).at(i) = 1.0 + p;
	    }
	}
    }

    checkNonzeroConnections();

}

void Network::resetWealthAverages() {

    for (int a=0; a<N; a++) {//for each agent
        agents.at(a).wealthAverage = agents.at(a).wealth;
        agents.at(a).wealthSquaredAverage = pow(agents.at(a).wealth, 2);
    }

}

void Network::updateWealthAverages() {

    for (int w=0; w<N; w++) {
	double newWealthAvg = (agents.at(w).wealthAverage * (time - 1) / (double)time) + (agents.at(w).wealth / (double)time);
	double newWealthSquaredAvg = (agents.at(w).wealthSquaredAverage * (time - 1) / (double)time) + (pow(agents.at(w).wealth, 2) / (double)time);
	//double newWealthAvg = 0.0001*agents.at(w).wealth + 0.9999*agents.at(w).wealthAverage; //EMA recursive 
	//double newWealthSquaredAvg = 0.0001*pow(agents.at(w).wealth, 2) + 0.9999*agents.at(w).wealthSquaredAverage; //EMA recursive 

	agents.at(w).wealthAverage = newWealthAvg;
	agents.at(w).wealthSquaredAverage = newWealthSquaredAvg;
    }

}

double Network::wealthMetric(bool ignoreRichest/*=false*/) {

    int richestAgentIndex = -1; //initialize to bad value
    if (ignoreRichest) {
	richestAgentIndex = richestAgent();
    }

    //Finding average historical wealth average:
    std::vector<double> wealthAverages;
    for (int a=0; a<N; a++) {
	if ((a != richestAgentIndex) || (!ignoreRichest)) {
	    wealthAverages.push_back(agents.at(a).wealthAverage);
	}
    }
    double totalWealthAverage = std::accumulate(wealthAverages.begin(), wealthAverages.end(), 0.0);
    double averageWealthAverage = totalWealthAverage / (double)wealthAverages.size();

    double wealthMetric = 0.0; //running sum at first
    for (int a=0; a<(int)wealthAverages.size(); a++) {
	wealthMetric += pow(wealthAverages.at(a) - averageWealthAverage, 2);
    }
    wealthMetric = wealthMetric / (double)wealthAverages.size();
    return wealthMetric;

}

int Network::richestAgent() {

    std::vector<double> wealths;
    for (int a=0; a<N; a++) {
	wealths.push_back(agents.at(a).wealth);
    }
    int richestAgent = std::max_element(wealths.begin(), wealths.end()) - wealths.begin(); 

    return richestAgent;

}

double Network::wealthVariance(bool ignoreRichest) {

    int richestAgentIndex = -1; //start with bad value
    if (ignoreRichest) {
        richestAgentIndex = richestAgent();
    }

    double runningSum = 0.0;
    int nTerms = 0;
    for (int a=0; a<N; a++) {
		if ((a != richestAgentIndex) || (!ignoreRichest)) {
	    double difference = agents.at(a).wealthSquaredAverage - pow(agents.at(a).wealthAverage, 2);
	    runningSum += difference;
            nTerms += 1;
        }
    }
    runningSum = runningSum / (double)nTerms;

    return runningSum;

}

double Network::energy() {

    //Assumes total wealth N

    //E = -1 + (1/N)*Sum[ w^2 ]
    double sum = 0;
    for (int a=0; a<N; a++) {
	sum += pow(agents.at(a).wealth, 2);
    }

    return -1.0 + sum / (double)N;

}

void Network::randomizeWealths() {

    double totalWealth = 0.0;
    for (int a=0; a<N; a++) {
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
	agents.at(a).wealth = x;
	totalWealth += x;
    }

    //Now we will normalize the total wealth to be N:
    double scaleFactor = N / totalWealth;
    for (int a=0; a<N; a++) {
	agents.at(a).wealth = agents.at(a).wealth * scaleFactor;
    }

}
    
void Network::make2DLattice(double R) { //R is range of interaction

    if ((N - ((int)sqrt(N) * (int)sqrt(N))) != 0) {
        std::cout << "make2DLattice warning! N = " << N << " is not a square!" << std::endl;
    }

    for (int i=0; i<N; i++) {
        for (int j=i+1; j<N; j++) {

            double interaction;

            double distance = checkDistance2D(i, j);
            if (distance <= R) {
                interaction = 1.0 + p;
            }
            else {
                interaction = 0.0;
            }

            matrix.at(i).at(j) = interaction;
            matrix.at(j).at(i) = interaction;

        }
    }

    checkNonzeroConnections();

}


double Network::checkDistance2D(int i, int j) {

    int rawDeltaX = abs(checkX(i) - checkX(j));
    int rawDeltaY = abs(checkY(i) - checkY(j)); 

    int deltaX = std::min((int)(sqrt(N) - rawDeltaX), rawDeltaX);
    int deltaY = std::min((int)(sqrt(N) - rawDeltaY), rawDeltaY);

    double distance = sqrt(deltaX*deltaX + deltaY*deltaY);
    return distance;

}

int Network::checkX(int i) {

    int x = i % (int)sqrt(N);
    return x;

}

int Network::checkY(int i) {

    int y = (int)(i/sqrt(N)) % (int)sqrt(N);
    return y;

}

void Network::destroyMatrix() {

    matrix.clear();

    for (int a=0; a<N; a++) {
        std::vector<double> lilEmpty(N, 0.0);
        matrix.push_back(lilEmpty);
    }

    checkNonzeroConnections();

}

double Network::phi() {

    //No longer assumes total wealth in the network is N
    int richestAgentIndex = richestAgent();

    double allWealth = totalWealth();
    double phi = allWealth - agents.at(richestAgentIndex).wealth;
    phi = phi / allWealth;

    return phi;

}

std::vector<int> Network::getClusterSizes() {
//Returns vector containing the size of each cluster in the matrix, not ordered

    //We need to create a pointer-array to pass as a reference to recursiveNeighborFunction
    bool* visited = new bool[N];
    for (int a=0; a<N; a++) {
        visited[a] = false;
    }

    //Creating output vector
    std::vector<int> clusterSizes;
    for (int a=0; a<N; a++) {
        if (!visited[a]) {//if we haven't counted the cluster that a belongs to yet

            int clusterSize = 0;
            recursiveNeighborFunction(a, clusterSize, visited); //modifies visited and clusterSize
            clusterSizes.push_back(clusterSize);

        }
    }

    delete [] visited;

    return clusterSizes;

}
                
                

void Network::recursiveNeighborFunction(int a, int& clusterSize, bool* visited) {
//used only for getClusterSize function

    visited[a] = true;
    clusterSize += 1;

    for (int j=0; j<N; j++) {
        if ((!visited[j]) && (matrix.at(a).at(j) != 0.0)) {
        //if j is an adjacent node that we haven't visited yet
            recursiveNeighborFunction(j, clusterSize, visited);
        }
    }

}

void Network::writeClusterSizes(std::string filename) {

    std::ofstream clusterSizesFile(filename);

    std::vector<int> clusterSizes = getClusterSizes();
    for (int c=0; c<(int)clusterSizes.size(); c++) {
        clusterSizesFile << clusterSizes.at(c) << std::endl;
    }

    clusterSizesFile.close();

}

    
std::vector<int> Network::getLargestCluster() {
//Returns vector containing the index of each agent in the largest cluster in the matrix, not ordered

    //Creating all-zero vector to store whether each agent has been visited
    std::vector<bool> previouslyVisited(N, false);

    //Creating output vector
    std::vector<int> clusterSizes;
    std::vector<std::vector<int>> clusters; //holds the index of agents that belong to each connected cluster
    for (int a=0; a<N; a++) {
        if (!previouslyVisited.at(a)) {//if we haven't counted the cluster that a belongs to yet

            //We need to create a pointer-array to pass as a reference to recursiveNeighborFunction
            bool* visited = new bool[N];
            for (int a=0; a<N; a++) {
                visited[a] = previouslyVisited.at(a);
            }

            int clusterSize = 0;
            recursiveNeighborFunction(a, clusterSize, visited); //modifies visited and clusterSize

            //Now we will find thw difference between previouslyVisited and visited to find the members of this cluster
            std::vector<int> cluster; //will hold the agents in this one cluster
            for (int a=0; a<N; a++) {
                if ((visited[a]) && (!previouslyVisited.at(a))) { 
                    cluster.push_back(a);
                }

                //updating list of nodes belonging to a completely discovered cluster
                previouslyVisited.at(a) = visited[a]; 
            }

            clusterSizes.push_back(clusterSize);
            clusters.push_back(cluster);

            delete[] visited;   
            
        }
    }

    int indexOfMaxCluster = -1; //initialize to bad value
    int sizeOfMaxCluster = -1; //initialize to bad value
    for (int c=0; c<(int)clusterSizes.size(); c++) { //for each cluster
        if ((int)clusters.at(c).size() != clusterSizes.at(c)) {
            std::cout << "getLargestCluster warning!  Apparently the cluster-finding algorithm isn't working!\n\n";
        }
        if (clusterSizes.at(c) > sizeOfMaxCluster) {
            sizeOfMaxCluster = clusterSizes.at(c);
            indexOfMaxCluster = c;
        }
    }


    return clusters.at(indexOfMaxCluster);

}


void Network::removeAgents(std::vector<int> agentsToRemove) {

    //Update matrix
    //update agent list

    std::vector<Agent> newAgents; //initializing empty
    std::vector<std::vector<double>> newMatrix; //also empty
    int nNewMatrixElements = 0; //used to make sure we have the right number of matrix elements at the end

    for (int a=0; a<N; a++) {
        if (std::find(agentsToRemove.begin(), agentsToRemove.end(), a) == agentsToRemove.end()) {
        //if agentsToRemove does not contain a:
            newAgents.push_back(agents.at(a));

            //Removing elements from the proper row of the matrix
            std::vector<double> newRow;
            for (int j=0; j<(int)matrix.at(a).size(); j++) {//for each agent in this row of the matrix
                if (std::find(agentsToRemove.begin(), agentsToRemove.end(), j) == agentsToRemove.end()) {
                //if agentsToRemove does not contain j:
                    newRow.push_back(matrix.at(a).at(j));
                }
            }
            newMatrix.push_back(newRow);
            nNewMatrixElements += (int)newRow.size();

        }
    }

    //Now that we have constructed our new agent list and adjacency matrix, we will replace the old ones
    agents = newAgents;
    matrix = newMatrix;
    N = (int)newAgents.size();

    //Checking the size of the new matrix
    if (nNewMatrixElements != (N*N)) {
        std::cout << "removeAgents Warning! The new matrix and/or agent list was not constructed properly\n";
    }
 
    checkNonzeroConnections();

}


void Network::scaleFreeMatrix(double exponent) {

    //Constructing normalization factor S
    double S = 0.0;
    for (int k=1; k<(N-1); k++) {
	S += pow(k, -1.0*exponent);
    }

    //These weights will be used to construct the connection probability for each pair
    std::vector<double> weights;
    for (int a=1; a<(N+1); a++) {
        //weights.push_back(kappa*pow(a, -1.0 / (double)(exponent - 1.0)));
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
        //double expectedDegree = 1.0 / pow(S*(1.0-x), -1.0 / exponent);
        double expectedDegree = 0.0;
        double baseline = 0.0;
        while (x >= baseline) {
            expectedDegree += 1.0;
            baseline += pow(expectedDegree, -1.0*exponent) / S;
        }
        weights.push_back(expectedDegree);
    }
    double totalWeight = std::accumulate(weights.begin(), weights.end(), 0.0);

    //Now it's time to make the connections!
    for (int i=0; i<N; i++) {

	for (int j=i+1; j<N; j++) {

	    double matrixElement;
            double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
	    double probability = weights.at(i)*weights.at(j) / totalWeight;


	    if (x < probability) {
		matrixElement = 1.0 + p;
	    }
	    else {
		matrixElement = 0.0;
	    }
	    
            matrix.at(i).at(j) = matrixElement;
            matrix.at(j).at(i) = matrixElement;

	}
    }

    checkNonzeroConnections();

}

void Network::prune() {//removes nodes that don't belong to the main connected cluster

    std::vector<int> largestCluster = getLargestCluster();
    std::vector<int> willRemove; //will hold the nodes to be removed from the network
    for (int a=0; a<N; a++) {
        if (std::find(largestCluster.begin(), largestCluster.end(), a) == largestCluster.end()) {
        //if a is not in the largest cluster
            //add to removal list
            willRemove.push_back(a);
        }
    }
    removeAgents(willRemove); //removing all but largest connected cluster from network

}

void Network::preferentialAttachmentMatrix(double lambda) {
//lambda is a parameter for forming bonds analogous to the other lambda in the growthSpurt function

    matrix.clear();
    //First we will add two nodes to the adjacency matrix.  Then we will fill the rest by preferential attachment
    std::vector<double> firstRow(N, 0.0);
    std::vector<double> secondRow(N, 0.0);
    //We will connect these two nodes:
    firstRow.at(1) = 1.0 + p;
    secondRow.at(0) = 1.0 + p;
    matrix.push_back(firstRow);
    matrix.push_back(secondRow);

    //Now we will do preferential attachment
    for (int i=2; i<N; i++) {
	std::vector<double> row(N, 0.0);

        std::vector<int> degrees; 
        double totalDegree = 0.0;
	for (int j=0; j<(int)matrix.size(); j++) {

            int degree = std::accumulate(matrix.at(j).begin(), matrix.at(j).end(), 0.0);
            degrees.push_back(degree);
            totalDegree += pow(degree, lambda);

	}

        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
        double baseline = 0.0;
        int jay = -1; //the index of the agent that we will attach to preferentially
        while ((x >= baseline) && (jay < ((int)matrix.size() - 1))) {
            jay += 1;
            baseline += pow(degrees.at(jay), lambda) / totalDegree;
        }
        row.at(jay) = 1.0 + p;
        matrix.at(jay).at(i) = 1.0 + p;
	
	//Finally we add this new row to the matrix!
	matrix.push_back(row);
    }

    checkNonzeroConnections();

}//Finished!

void Network::scaleFreeConfigMatrix(double exponent, int minDegree) {

    //Uses configuration model to generate a scale free network
 
    //First step is to create an empty matrix:
    matrix.clear();
    for (int a=0; a<N; a++) {
        std::vector<double> emptyRow(N, 0.0);
        matrix.push_back(emptyRow);
    }

    //Constructing normalization factor S
    double S = 0.0;
    for (int k=minDegree; k<(N-1); k++) {
	S += pow(k, -1.0*exponent);
    }

    //Constructing degree of each node
    std::vector<int> degrees;
    for (int a=0; a<N; a++) {
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
        int expectedDegree = minDegree - 1;
        double baseline = 0.0;
        while (x >= baseline) {
            expectedDegree += 1;
            baseline += pow(expectedDegree, -1.0*exponent) / S;
        }
        degrees.push_back(expectedDegree);
    }
    int totalDegree = std::accumulate(degrees.begin(), degrees.end(), 0.0);
    if ((totalDegree % 2) != 0) { //if total degree is odd
        //we cannot accept this
        //choosing random agent to give an extra degree to
	int agentIndex = rand() % N;
        degrees.at(agentIndex) += 1;
        totalDegree += 1;
    }

    //Now it's time to conenct the stubs we made!
    std::vector<int> stillHaveStubs; //list of agent indices that are yet to be fully connected
    for (int a=0; a<N; a++) {
        stillHaveStubs.push_back(a);
    }
    while ((int)stillHaveStubs.size() > 0) {

        //Choose random pair of agents
	int index1 = rand() % (int)stillHaveStubs.size();
        int agentIndex1 = stillHaveStubs.at(index1);
	int index2 = rand() % (int)stillHaveStubs.size();
        int agentIndex2 = stillHaveStubs.at(index2);
        //We're necessarily allowing self-loops and multi-links..  We'll manually take these out later

        //Connect these guys!
        matrix.at(agentIndex1).at(agentIndex2) = 1.0 + p;
        matrix.at(agentIndex2).at(agentIndex1) = 1.0 + p;

        degrees.at(agentIndex1) -= 1;
        degrees.at(agentIndex2) -= 1;
        if (degrees.at(agentIndex1) <= 0) {
            auto foundIter = std::find(stillHaveStubs.begin(), stillHaveStubs.end(), agentIndex1);
            stillHaveStubs.erase(foundIter);
        }
        if (degrees.at(agentIndex2) <= 0) {
            auto foundIter = std::find(stillHaveStubs.begin(), stillHaveStubs.end(), agentIndex2);
            if (foundIter != stillHaveStubs.end()) { //incase agentIndex1 equaled agentIndex2
                stillHaveStubs.erase(foundIter);
            }
        }

    }


    //Now it's time to clean up this network
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if (i == j) { //getting rid of self-loops
                matrix.at(i).at(j) = 0.0;
            }
            else { //getting rid of multi-links
                matrix.at(i).at(j) = std::min(matrix.at(i).at(j), (1.0 + p));
            }
        }
    }
    //Good grief, man

    checkNonzeroConnections();

}

void Network::breakAndBind(double breakRate, int minDegree) {
//Destroys bonds between nodes with probability breakRate.  Then attaches new bonds between possibly different nodes
//proportionally to w^lambda.
//Bonds will not break if it would leave a node with degree less than minDegree.

    std::vector<int> degrees;
    for (int a=0; a<N; a++) {
        double rawDegree = std::accumulate(matrix.at(a).begin(), matrix.at(a).end(), 0);
        int degree = rawDegree / (1.0 + p); //generally, a connected matrix element is 1 + p
        degrees.push_back(degree);
    }

    int nBroken = 0; //The number of bonds we've broken so far
    //The update is not done in parallel to avoid degree dropping below minDegree
    for (int a=0; a<N; a++) {
        if (degrees.at(a) > minDegree) {
        //if this node is well connected enough to break some bonds
            for (int j=a+1; j<N; j++) {
                if ((degrees.at(a) > minDegree) && (degrees.at(j) > minDegree)) {
                //if both of these nodes still have high enough degree
                    if (matrix.at(a).at(j) != 0.0) {//if they're connected agents
                        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
                        if (x < breakRate) { //if we've decided to break this bond
                            matrix.at(a).at(j) = 0.0;
                            matrix.at(j).at(a) = 0.0;
                            degrees.at(a) -= 1;
                            degrees.at(j) -= 1;
                            nBroken += 1;
                        }
                    }
                }
            }
        }
    }

    //Now it's time to make new bonds!
    //Calculating sum constant S:
    double S = 0.0;
    for (int a=0; a<N; a++) {
        S += pow(agents.at(a).wealth, lambda);
    }
    int nMade = 0; //the number of new bonds we have made so far
    while (nMade < nBroken) {//for each new bond we want to make

        //Choosing first agent to make the new bond at                        
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
        int agentIndex1 = -1; //initialize to bad value
        double baseline = 0.0;
        while ((x > baseline) || (agentIndex1 < 0)) {
            agentIndex1 += 1;
            baseline += pow(agents.at(agentIndex1).wealth, lambda) / S;
        }
        //Choosing second agent:
        int agentIndex2 = agentIndex1; //initializing to bad value
        while (agentIndex2 == agentIndex1) {
            double y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
            agentIndex2 = -1; //initialize to bad value
            baseline = 0.0;
            while ((y > baseline) || (agentIndex2 < 0)) {
                agentIndex2 += 1;
                baseline += pow(agents.at(agentIndex2).wealth, lambda) / S;
            }
        }

        if (matrix.at(agentIndex1).at(agentIndex2) == 0.0) {
        //if the nodes aren't already connected
            matrix.at(agentIndex1).at(agentIndex2) = 1.0 + p;
            matrix.at(agentIndex2).at(agentIndex1) = 1.0 + p;
            nMade += 1;
        }

    }

    checkNonzeroConnections();

}
                

void Network::trickleDownGrowthSpurt(double intensity) {

    double fi = phi();

    //We say that the normal richestWealth is 1
    double mew = mu*(1.0 - intensity*log(fi)); //should be positive

    growthSpurt(mew);

}



void Network::equalityDrivenGrowthSpurt() {

    double fi = phi();
    double mew = mu*fi;
    growthSpurt(mew);

}

void Network::inequalityDrivenGrowthSpurt() {

    double fi = phi();
    double mew = mu*(1.0 - fi);
    growthSpurt(mew);

}

void Network::growthSpurt(double mew) {

    //Check out this cool trick!
    double originalMu = mu;
    mu = mew;
    growthSpurt();
    mu = originalMu;

}


void Network::inequalityWeakenedGrowthSpurt(double intensity) { 

    double fi = phi();
    double smallestWealth = agents.at(poorestAgent()).wealth;


    double antimu = intensity*log(fi)*mu; //should be negative
    growthSpurt(antimu + mu);

    //Now we make sure that no one's wealth becomes negative
    for (int a=0; a<N; a++) {
        //totalAbsoluteWealth += std::max(totalAbsoluteWealth*(smallestWealth - agents.at(a).wealth), 0.0) / (double)N;
        agents.at(a).wealth = std::max(agents.at(a).wealth, smallestWealth);
    }

}


int Network::poorestAgent() {

    std::vector<double> wealths;
    for (int a=0; a<N; a++) {
	wealths.push_back(agents.at(a).wealth);
    }
    int poorestAgent = std::min_element(wealths.begin(), wealths.end()) - wealths.begin(); 

    return poorestAgent;

}

double Network::totalWealth() {
    
    double total = 0.0;
    for (int a=0; a<N; a++) {
        total += agents.at(a).wealth;
    }
    return total;

}

double Network::getDegree(int agentIndex) {
//writes the degree of node given by agentIndex.  Assumes connections are either one or zero.

    return std::accumulate(matrix.at(agentIndex).begin(), matrix.at(agentIndex).end(), 0);

}

double Network::averageDegree() {

    double totalDegree = 0.0;
    for (int a=0; a<N; a++) {
        totalDegree += getDegree(a);
    }

    return (totalDegree / (double)N);

}

void Network::normalizeFitnesses() {

    double totalFitness = 0.0;
    for (int a=0; a<N; a++) {
        totalFitness += agents.at(a).fitness;
    }

    double scaleFactor = (double)N / totalFitness;
    for (int a=0; a<N; a++) {
        agents.at(a).fitness = scaleFactor*agents.at(a).fitness;
    }

}

void Network::makeHub(int index, double probability) {
//Connects the node given by index to each of its neighbors with probability probability

    for (int a=0; a<N; a++) {//for each agent
        if (a != index) { //except for the new hub itself
            double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
            if (x < probability) {
                //making bidirectional connection
                matrix.at(index).at(a) = 1.0 + p;
                matrix.at(a).at(index) = 1.0 + p;
            }
        }
    }

}


void Network::hubifyMatrix(double hubProbability, double hubStrength) {
//Chooses random nodes in the network to become hubs.  New hubs are created each with probability hubProbability.
//The average number of neighbors that a new hub has is determined by hubStrength

    for (int a=0; a<N; a++) {
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
        if (x < hubProbability) {
            makeHub(a, hubStrength);
        }
    }
                   
    checkNonzeroConnections();

}


int Network::getXNeighbor(int homeIndex, int distance) {
//Returns the index of the agent distance away from homeIndex on the x axis for a 2D lattice with PBC
    int L = sqrt(N);
    int rawX = checkX(homeIndex) + distance;
    while (rawX < 0) {//we need to make it positive because of surprising behavior
	rawX += L;
    }
    int X = rawX % L;
    int Y = checkY(homeIndex);
    return getIndex(X, Y);

}

int Network::getYNeighbor(int homeIndex, int distance) {
//Returns the index of the agent distance away from homeIndex on the x axis for a 2D lattice with PBC
    int L = sqrt(N);
    int rawY = checkY(homeIndex) + distance;
    while (rawY < 0) {//we need to make it positive because of surprising behavior
	rawY += L;
    }
    int Y = rawY % L;
    int X = checkX(homeIndex);
    return getIndex(X, Y);

}

int Network::getIndex(int x, int y) {
//Returns agent index given agent's x and y coordinates
//Assumes 2D lattice
//x and y must each be in [0,L)
    
    int L = sqrt(N);
    return x + (L*y);
}

std::vector<double> Network::correlations2D(int maxLength) {
//Returns a vector of length maxLength+1.  The first element is self correlation which should be one.
//In general output[i] is the correlation between agents separated a distance i.
//C[i] = Avg[ (w0 - wavg) * (wi - wavg) / ( stddev(w0 -wavg) * stddev(wi - wavg) ) ]
//C[i] = Avg[ (w0 - wavg) * (wi - wavg) / var(w - wavg) ]

    std::vector<double> correlations(maxLength+1, 0.0); //initializing to zero

    for (int a=0; a<N; a++) {
	double innovationa = agents.at(a).wealth - agents.at(a).wealthAverage;
	double stdDeva = pow((agents.at(a).wealthSquaredAverage - agents.at(a).wealthAverage*agents.at(a).wealthAverage), 0.5);
	for (int r=0; r < (maxLength+1); r++) {

	    int xNeighbor = getXNeighbor(a, r);
	    double innovationx = agents.at(xNeighbor).wealth - agents.at(xNeighbor).wealthAverage;
	    double stdDevx = pow((agents.at(xNeighbor).wealthSquaredAverage - agents.at(xNeighbor).wealthAverage*agents.at(xNeighbor).wealthAverage), 0.5);
            correlations.at(r) += innovationa * innovationx / (stdDeva * stdDevx);
	    
	    int yNeighbor = getYNeighbor(a, r);
	    double innovationy = agents.at(yNeighbor).wealth - agents.at(yNeighbor).wealthAverage;
	    double stdDevy = pow((agents.at(yNeighbor).wealthSquaredAverage - agents.at(yNeighbor).wealthAverage*agents.at(yNeighbor).wealthAverage), 0.5);
            correlations.at(r) += innovationa * innovationy / (stdDeva * stdDevy);

	}

        if (getIndex(checkX(a), checkY(a)) != a) {
            std::cout << "getIndex and/or checkX etc. functions are not working right (2D)" << std::endl;
        }

    }

    //Now we normalize by the appropriate term
    for (int r=0; r < (maxLength+1); r++) {
	correlations.at(r) = correlations.at(r) / (2 * N); //2 for 2 directions
    }

    return correlations;

}

void Network::writeCorrelations(std::string filename, int maxLength, int dim) {
//Computes correlations and writes to file

    std::ofstream correlationsFile(filename);
    std::vector<double> correlations;
    if (dim == 2) {
        correlations = correlations2D(maxLength);
    }
    else {
        if (dim == 1) {
            correlations = correlations1D(maxLength);
        }
        else {
            if (dim == 3) {
                correlations = correlations3D(maxLength);
            }
        }
    }

    for (int r=0; r<=maxLength; r++) {
	correlationsFile << correlations.at(r) << std::endl;
    }
    correlationsFile.close();

}

void Network::makeRing(double R) {
//Connects agents in a 1D chain with periodic boundary conditions
//R is interaction range between agents

    
    for (int i=0; i<N; i++) {
        for (int j=i+1; j<N; j++) {

            double interaction;

            double distance = checkDistance1D(i, j);
            if (distance <= R) {
                interaction = 1.0 + p;
            }
            else {
                interaction = 0.0;
            }

            matrix.at(i).at(j) = interaction;
            matrix.at(j).at(i) = interaction;

        }
    }

    checkNonzeroConnections();

}

int Network::checkDistance1D(int i, int j) {
//Checks distance between two agents on a 1D chain with periodic boundary conditions

    int rawSeparation = abs(j - i);
    int distance = std::min((N - rawSeparation), rawSeparation);
    return distance;

}

std::vector<double> Network::correlations1D(int maxLength) {
//Returns a vector of length maxLength+1.  The first element is self correlation which should be one.
//In general output[i] is the correlation between agents separated a distance i.
//C[i] = Avg[ (w0 - wavg) * (wi - wavg) / ( stddev(w0 -wavg) * stddev(wi - wavg) ) ]
//C[i] = Avg[ (w0 - wavg) * (wi - wavg) / var(w - wavg) ]

    std::vector<double> correlations(maxLength+1, 0.0); //initializing to zero

    for (int a=0; a<N; a++) {
	double innovationa = agents.at(a).wealth - agents.at(a).wealthAverage;
	double stdDeva = pow((agents.at(a).wealthSquaredAverage - agents.at(a).wealthAverage*agents.at(a).wealthAverage), 0.5);
	for (int r=0; r < (maxLength+1); r++) {

	    int neighbor = (a + r) % N;
	    double innovationr = agents.at(neighbor).wealth - agents.at(neighbor).wealthAverage;
	    double stdDevr = pow((agents.at(neighbor).wealthSquaredAverage - agents.at(neighbor).wealthAverage*agents.at(neighbor).wealthAverage), 0.5);
            correlations.at(r) += innovationa * innovationr / (stdDeva * stdDevr);

	}

        if (((a+0)%N) != a) {
            std::cout << "getIndex and/or checkX etc. functions are not working right (1D)" << std::endl;
        }

    }

    //Now we normalize by the number of agents
    for (int r=0; r < (maxLength+1); r++) {
	correlations.at(r) = correlations.at(r) / N;
    }

    return correlations;

}

int Network::getIndex(int x, int y, int z) { //overloaded function
//Returns agent index given agent's x, y, and z coordinates
//Assumes 3D lattice
//x, y, and z must each be in [0,L)
    
    int L = round(pow(N, (1.0 / 3.0)));
    return x + (L*y) + (L*L*z);
}


int Network::checkX3D(int i) {

    int x = i % (int)round(pow(N, (1.0/3.0)));
    return x;

}

int Network::checkY3D(int i) {

    int L = round(pow(N, (1.0/3.0)));

    int y = ((int)(i/(double)L)) % L;
    return y;

}

int Network::checkZ3D(int i) {

    int L = round(pow(N,(1.0/3.0)));
    int L2 = L*L;

    int z = ((int)(i/(double)L2)) % L;
    return z;

}

int Network::getXNeighbor3D(int homeIndex, int distance) {
//Returns the index of the agent distance away from homeIndex on the x axis for a 3D lattice with PBC
    int L = round(pow(N, (1.0/3.0)));
    int rawX = checkX3D(homeIndex) + distance;
    while (rawX < 0) {//we need to make it positive because of surprising behavior
	rawX += L;
    }
    int X = rawX % L;
    int Y = checkY3D(homeIndex);
    int Z = checkZ3D(homeIndex);
    return getIndex(X, Y, Z);

}

int Network::getYNeighbor3D(int homeIndex, int distance) {
//Returns the index of the agent distance away from homeIndex on the x axis for a 3D lattice with PBC
    int L = round(pow(N, (1.0/3.0)));
    int rawY = checkY3D(homeIndex) + distance;
    while (rawY < 0) {//we need to make it positive because of surprising behavior
	rawY += L;
    }
    int Y = rawY % L;
    int X = checkX3D(homeIndex);
    int Z = checkZ3D(homeIndex);
    return getIndex(X, Y, Z);

}

int Network::getZNeighbor3D(int homeIndex, int distance) {
//Returns the index of the agent distance away from homeIndex on the x axis for a 3D lattice with PBC
    int L = round(pow(N, (1.0/3.0)));
    int rawZ = checkZ3D(homeIndex) + distance;
    while (rawZ < 0) {//we need to make it positive because of surprising behavior
	rawZ += L;
    }
    int Z = rawZ % L;
    int X = checkX3D(homeIndex);
    int Y = checkY3D(homeIndex);
    return getIndex(X, Y, Z);

}

std::vector<double> Network::correlations3D(int maxLength) {
//Returns a vector of length maxLength+1.  The first element is self correlation which should be one.
//In general output[i] is the correlation between agents separated a distance i.
//C[i] = Avg[ (w0 - w0avg) * (wi - wiavg) / ( stddev(w0) * stddev(wi) ) ]

    std::vector<double> correlations(maxLength+1, 0.0); //initializing to zero

    for (int a=0; a<N; a++) {
	double innovationa = agents.at(a).wealth - agents.at(a).wealthAverage;
	double stdDeva = pow((agents.at(a).wealthSquaredAverage - agents.at(a).wealthAverage*agents.at(a).wealthAverage), 0.5);
	for (int r=0; r < (maxLength+1); r++) {
	    
	    int xNeighbor = getXNeighbor3D(a, r);
	    double innovationx = agents.at(xNeighbor).wealth - agents.at(xNeighbor).wealthAverage;
	    double stdDevx = pow((agents.at(xNeighbor).wealthSquaredAverage - agents.at(xNeighbor).wealthAverage*agents.at(xNeighbor).wealthAverage), 0.5);
            correlations.at(r) += innovationa * innovationx / (stdDeva * stdDevx);

	    int yNeighbor = getYNeighbor3D(a, r);
	    double innovationy = agents.at(yNeighbor).wealth - agents.at(yNeighbor).wealthAverage;
	    double stdDevy = pow((agents.at(yNeighbor).wealthSquaredAverage - agents.at(yNeighbor).wealthAverage*agents.at(yNeighbor).wealthAverage), 0.5);
            correlations.at(r) += innovationa * innovationy / (stdDeva * stdDevy);

	    int zNeighbor = getZNeighbor3D(a, r);
	    double innovationz = agents.at(zNeighbor).wealth - agents.at(zNeighbor).wealthAverage;
	    double stdDevz = pow((agents.at(zNeighbor).wealthSquaredAverage - agents.at(zNeighbor).wealthAverage*agents.at(zNeighbor).wealthAverage), 0.5);
            correlations.at(r) += innovationa * innovationz / (stdDeva * stdDevz);

	}

        if (getIndex(checkX3D(a), checkY3D(a), checkZ3D(a)) != a) {
            std::cout << "getIndex and/or checkX3D etc. functions are not working right (3D)" << std::endl;
        }

    }

    //Now we normalize by the number of terms we added
    for (int r=0; r < (maxLength+1); r++) {
	correlations.at(r) = correlations.at(r) / (double)(3 * N); //3 for 3 directions
    }

    return correlations;

}

void Network::makeLattice(int dim, double R) { //R is range of interaction

    int L = round(pow(N, 1.0/dim));

    if ((N - (pow(L, dim))) != 0) {
        std::cout << "makeLattice warning! N = " << N << " is not a square/cube!" << std::endl;
    }

    for (int i=0; i<N; i++) {
        for (int j=i+1; j<N; j++) {

            double interaction;

            double distance = -1.0; //initialize to bad value
            if (dim == 2) {
                distance = checkDistance2D(i, j);
            }
            if (dim == 3) {
                distance = checkDistance3D(i, j);
            }
            if (distance < 0.0) {
                std::cout << "ValueError!  dim must be 2 or 3 (makeLattice function) ";
            }
            if (distance <= R) {
                interaction = 1.0 + p;
            }
            else {
                interaction = 0.0;
            }

            matrix.at(i).at(j) = interaction;
            matrix.at(j).at(i) = interaction;

        }
    }

    checkNonzeroConnections();

}

double Network::checkDistance3D(int i, int j) {

    int L = round(pow(N, 1.0/3.0));

    int rawDeltaX = abs(checkX3D(i) - checkX3D(j));
    int rawDeltaY = abs(checkY3D(i) - checkY3D(j)); 
    int rawDeltaZ = abs(checkZ3D(i) - checkZ3D(j)); 

    int deltaX = std::min((L - rawDeltaX), rawDeltaX);
    int deltaY = std::min((L - rawDeltaY), rawDeltaY);
    int deltaZ = std::min((L - rawDeltaZ), rawDeltaZ);

    double distance = sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    return distance;

}

void Network::makeAnalogLattice(int dim, double correlationLength) { //R is range of interaction

    int L = round(pow(N, 1.0/dim));
    double c = 1.0 - (1.0 / correlationLength); //recursive diminishing factor

    int maxLength = 4*correlationLength;
    std::vector<double> diminisher(maxLength, 0.0);
    for (int r=0; r<maxLength; r++) {
        diminisher.at(r) = pow(c, r); //this is expensive, so we store values in an array
    }

    if ((N - (pow(L, dim))) != 0) {
        std::cout << "makeLattice warning! N = " << N << " is not a square/cube!" << std::endl;
    }

    for (int i=0; i<N; i++) {
        for (int j=i+1; j<N; j++) {


            double distance = -1.0; //initialize to bad value
            if (dim == 1) {
                distance = checkDistance1D(i, j);
            }
            if (dim == 2) {
                distance = checkDistance2D(i, j);
            }
            if (dim == 3) {
                distance = checkDistance3D(i, j);
            }
            if (distance < 0.0) {
                std::cout << "ValueError!  dim must be 1, 2, or 3 (makeLattice function) ";
            }

            double interaction;
            if (distance < (maxLength-1)) {
                interaction = (1.0 + p)*diminisher.at((int)round(distance));
            }
            else {
                interaction = 0.0;
            }

            matrix.at(i).at(j) = interaction;
            matrix.at(j).at(i) = interaction;

        }
    }

    checkNonzeroConnections();

}

void Network::equalizeWealths() {

    for (int a=0; a<N; a++) {
        agents.at(a).wealth = 1.0;
    }

}

void Network::exchangeStepNonzeroOnly(int n) {
    //i is number of steps you'd like to take

    for (int j=0; j<n; j++) {//for each step we want to take
	for (int k=0; k<N; k++) {//for each of N exchanges per timestep
	    //Choose a pair that has a nonzero connection
            double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
	    int randomIndex = std::min(x*(int)nonZeros.size(), (double)((int)nonZeros.size() - 1));
	    int index1 = nonZeros.at(randomIndex).at(0);
	    int index2 = nonZeros.at(randomIndex).at(1);

	    exchange(index1, index2);
	}

	time += 1;

    }

}

void Network::checkNonzeroConnections() {
    //creates vector of 2-element vectors corresponding to the indices of nonzero elements of the connectivity matrix

    nonZeros.clear(); //emptying our collection before we update
    for (int i=0; i<N; i++) {
	for (int j=0; j<N; j++) {
	    if (matrix.at(i).at(j) != 0.0) {
		std::vector<int> nonZero({i, j});
		nonZeros.push_back(nonZero);
            }
	}
    }

}

void Network::trickleDownEnergyGrowthSpurt() {

    double mew = energy(); //should be positive
    mew = std::pow(mew, mew); //extreme!

    growthSpurt(mew);

}

void Network::updateGuess(int agentIndex, bool justWon) {
	int delta = 2*justWon - 1;
	agents.at(agentIndex).stepSizeM = std::max(agents.at(agentIndex).stepSizeM + delta, 0);

        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
        agents.at(agentIndex).guess = agents.at(agentIndex).guess + x / std::pow(learnStep, agents.at(agentIndex).stepSizeM);
	//Making sure it's in the range of [0,1]:
	agents.at(agentIndex).guess -= std::floor(agents.at(agentIndex).guess);
}


void Network::writeGuessesToFile(std::string filename) {

    std::ofstream assetFile(filename);
    for (int a=0; a<N; a++) {
	assetFile << agents.at(a).guess << std::endl;
    }
    assetFile.close();

}

void Network::writeStepSizesToFile(std::string filename) {

    std::ofstream assetFile(filename);
    for (int a=0; a<N; a++) {
	assetFile << agents.at(a).stepSizeM << std::endl;
    }
    assetFile.close();

}

void Network::randomizeGuesses() {

    for (int a=0; a<N; a++) {
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX); //random number in [0, 1]
	agents.at(a).guess = x;
    }

}

void Network::loadWealthsFromFile(std::string filename) {

    //Reading in file now
    std::ifstream File(filename);
    double c; //c is the connection between two nodes (an element in the adjacency matrix)
    int n = 0; //The number of matrix elements in the file you read
    while (File >> c) {
	agents.at(n).wealth = c;
	n += 1;
    }

    //Checking that you read an appropriate file
    if (n != N) {
	std::cout << "Warning: wrong number of wealths in " << filename << ".  Expected " << N << " for a network of size " << N << "; found " << n << "." << std::endl;
    }
    else {
	//std::cout << "Successfully loaded " << n << " wealths from " << filename << "!" << std::endl;
    }

}

double Network::dependentLambda(double wealth) {
    //gives an sgent's individual lambda as a function of its wealth
    double individualLambda;
    if (wealth < nukeLimit) {
	individualLambda = lambda;
    }
    else {
	individualLambda = wealth / nukeLimit;
    }
    return individualLambda;
}
