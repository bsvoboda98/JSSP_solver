#include "AlgorithmManager.h"

#include <string>
#include <iostream>

int main(int argc, char* argv[]) {
    std::vector<std::string> config;
    for(int i = 1; i < argc; i++) {

        std::string arg = argv[i];
        config.push_back(arg);
    }


    if(!config.empty()) {
        if(config[0] == "-co") {
            /**
             * @param config -co
             *               -g popSize maxGen crossRate mutateRate mutateImpact
             *               -o orCount orImpact
             *               -i instances
             *               -f crossover mutate select
             */
            AlgorithmManager::runCombined(config);
        }
        else if(config[0] == "-bo") {
            /**
             * @param config -bo
             *               -g popSize maxGen crossRate mutateRate mutateImpact
             *               -o orCount orImpact
             *               -i instances
             *               -f crossover mutate select
             */
            AlgorithmManager::runBottleneck(config);
        }
        else if(config[0] == "-mi") {
            /**
             * @param config -mi
             *               -g popSize maxGen crossRate mutRate mutImpact
             *               -i instances
             *               -f crossover mutate select
             */
            AlgorithmManager::runMultipleCombined(config);
        }
        else if(config[0] == "-jo") {
            /**
             * @param config -jo
             *               -g popSize maxGen crossRate mutateRate mutateImpact
             *               -o orCount orImpact
             *               -i instances
             *               -f crossover mutate select
             */
            AlgorithmManager::runJob(config);
        }
        else if(config[0] == "-ga") {
            /**
             * @param config -ga
             *               -g popSize maxGen crossRate mutateRate mutateImpact
             *               -i instances
             *               -f crossover mutate select
             */
            AlgorithmManager::runGA(config);
        }
        else if(config[0] == "-or") {
            /**
             * @param config -or
             *               -i instances
             */
            AlgorithmManager::runOR(config);
        }
        else if(config[0] == "-m") {
            /**
             * @param config -m
             *               -g popSize maxGen crossRate mutateRate mutateImpact
             *               -i instances
             *               -f crossover mutate select
             */
            AlgorithmManager::runMemetic(config);
        }
        else if(config[0] == "-t") {
            AlgorithmManager::runTest(config);
        }
        else if(config[0] == "-help") {
            //print help
        }
        return 0;
    }
    else {
        std::string userInput;

        std::cout << "WHICH SOLVER SHOULD BE USED? \n combined (c) ~standard\n bottleneck (b)\n multiple intervals (mi)\n jobs (j) \n GA \n OR" << std::endl;
        std::string solverDefinition;
        while(true) {
            std::getline(std::cin, solverDefinition);
            if(solverDefinition.empty()) break;
            if(solverDefinition == "c") break;
            if(solverDefinition == "b") break;
            if(solverDefinition == "mi") break;
            if(solverDefinition == "j") break;
            if(solverDefinition == "GA") break;
            if(solverDefinition == "OR") break;
            if(solverDefinition == "test") break;
            std::cout << "INVALID INPUT" << std::endl;
        }

        if(solverDefinition == "combined" || solverDefinition == "c" ||solverDefinition.empty()) {
            AlgorithmManager::runCombined({});
        }

        if(solverDefinition == "bottleneck" || solverDefinition == "b") {
            AlgorithmManager::runBottleneck({});
        }

        if(solverDefinition == "multiple intervals" || solverDefinition == "mi") {
            AlgorithmManager::runMultipleCombined({});
        }

        if(solverDefinition == "j") {
            AlgorithmManager::runJob({});
        }

        if(solverDefinition == "OR") {
            AlgorithmManager::runOR({});
        }

        if(solverDefinition == "AR") {
            AlgorithmManager::runGA({});
        }

        if(solverDefinition == "test") {
            AlgorithmManager::runTest({});
        }

        return 1;

    }

}




