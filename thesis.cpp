#include <bits/stdc++.h>
using namespace std;

int nodes,links,zones;
vector<vector<int>>node_costs;

void Sioux_Falls_Network_Assign(){

    // The total number of nodes, links and zones present in the Sioux Falls Network
    nodes=24;
    links=76;
    zones=24;
}

void Anaheim_Assign(){

    // The total number of nodes, links and zones present in the Anaheim Network
    nodes=416;
    links=914;
    zones=38;
}

void Sioux_Falls_Network(vector<double>&linkFlows, vector<double>&performanceFunctions,vector<vector<double>>&demandMatrix){

    freopen("Sioux_Falls_Costs.txt","r",stdin);

    vector<int>ran;
    double val;
    char ch;
    for(int i=0;i<links;i++){
        for(int j=0;j<10;j++){
            cin>>val;
            if(j==0 || j==1 || j==3 || j==4){
                ran.push_back(val);
            }
        }
        cin>>ch;
        node_costs.push_back(ran);
        ran.clear();
    }

    for(int i=0;i<node_costs.size();i++){
        demandMatrix[node_costs[i][0]][node_costs[i][1]]=node_costs[i][2];
        demandMatrix[node_costs[i][1]][node_costs[i][0]]=node_costs[i][2];
        performanceFunctions[i]=node_costs[i][3];
    }

    for(int i=0;i<links;i++){
        for(int j=0;j<4;j++){
            int val;
            cin>>val;
            if(j==2){
                linkFlows[i]=val;
            }
        }
    }

}

void Anaheim_Network(vector<double>&linkFlows, vector<double>&performanceFunctions,vector<vector<double>>&demandMatrix){

    freopen("Anaheim_Costs.txt","r",stdin);

    vector<int>ran;
    double val;
    char ch;
    for(int i=0;i<links;i++){
        for(int j=0;j<10;j++){
            cin>>val;
            if(j==0 || j==1 || j==3 || j==4){
                ran.push_back(val);
            }
        }
        cin>>ch;
        node_costs.push_back(ran);
        ran.clear();
    }

    for(int i=0;i<node_costs.size();i++){
        demandMatrix[node_costs[i][0]][node_costs[i][1]]=node_costs[i][2];
        demandMatrix[node_costs[i][1]][node_costs[i][0]]=node_costs[i][2];
        performanceFunctions[i]=node_costs[i][3];
    }

    for(int i=0;i<links;i++){
        for(int j=0;j<4;j++){
            int val;
            cin>>val;
            if(j==2){
                linkFlows[i]=val;
            }
        }
    }

}

vector<double> frankWolfeAlgorithm(vector<double>& linkFlows, const vector<double>& performanceFunctions, const vector<std::vector<double>>& demandMatrix,
 int numIterations) {
    int numLinks = links;
    int numNodes = nodes;

    vector<double> newLinkFlows(numLinks);

    for (int iter = 0; iter < numIterations; ++iter) {
        vector<double> travelCosts(numLinks);
        for (int i = 0; i < numLinks; ++i) {
            for (int j = 0; j < numLinks; ++j) {
                travelCosts[i] += performanceFunctions[j] * linkFlows[j];
            }
        }

        vector<double> routeChoiceProportions(numLinks);
        for (int i = 0; i < numLinks; ++i) {
            routeChoiceProportions[i] = exp(-travelCosts[i]);
        }

        double normalizationFactor = 0.0;
        for (int i = 0; i < numLinks; ++i) {
            normalizationFactor += routeChoiceProportions[i];
        }

        for (int i = 0; i < numLinks; ++i) {
            newLinkFlows[i] = 0.0;
            for (int j = 0; j < numNodes; ++j) {
                newLinkFlows[i] += demandMatrix[j][i] * routeChoiceProportions[i] / normalizationFactor;
            }
        }

        double stepSize = 2.0 / (2 + iter);
        for (int i = 0; i < numLinks; ++i) {
            linkFlows[i] = (1 - stepSize) * linkFlows[i] + stepSize * newLinkFlows[i];
        }
    }

    return linkFlows;
}

int main() {

    freopen("Sioux_Falls_Network_Output.txt","w",stdout);
    Sioux_Falls_Network_Assign();
    int numIterations1 = 100;
    vector<double> linkFlows(links,0.0);
    vector<double> performanceFunctions(links,0.0);
    vector<vector<double>>demandMatrix(nodes+1,vector<double>(nodes+1,0.0));

    Sioux_Falls_Network(linkFlows,performanceFunctions,demandMatrix);

    
    freopen("Anaheim_Output.txt","w",stdout);
    Anaheim_Assign();
    int numIterations2 = 1000;
    vector<double> linkFlows(links,0.0);
    vector<double> performanceFunctions(links,0.0);
    vector<vector<double>>demandMatrix(nodes+1,vector<double>(nodes+1,0.0));

    Anaheim_Network(linkFlows,performanceFunctions,demandMatrix);

    // Iterative solution algorithm (Frank-Wolfe)
    linkFlows = frankWolfeAlgorithm(linkFlows, performanceFunctions, demandMatrix, numIterations1);
    linkFlows = frankWolfeAlgorithm(linkFlows, performanceFunctions, demandMatrix, numIterations2);

    // Calculate LCI scores for each link (assuming uniform path weights)
    vector<double> pathWeights(linkFlows.size(), 1.0);
    vector<double> LCIScores(links);
    for (int i = 0; i < linkFlows.size(); ++i) {
        LCIScores[i]=linkFlows[i]*pathWeights[i];
    }

    // Obtain criticality rankings of links based on LCI scores
    vector<int> criticalityRankings(linkFlows.size());
    iota(criticalityRankings.begin(), criticalityRankings.end(), 0);
    sort(criticalityRankings.begin(), criticalityRankings.end(), [&LCIScores](int i, int j) { return LCIScores[i] > LCIScores[j]; });

    // Visualize results (not implemented in this example)
    cout << "Criticality rankings are given as:"<<endl<<endl;

    for (int i = 0; i < links; ++i) {
        cout << "Link " << node_costs[i][0]<<" -> "<<node_costs[i][1]<< " has Criticality rank of "<<criticalityRankings[i]<<endl;
    }

    return 0;
}