//
// Created by svobee on 21.04.24.
//
#include "AlgorithmManager.h"
#include "Structs.h"
#include "GeneticAlgorithm.h"
#include "OrtoolsAlgorithm.h"
#include "FunctionFactory.h"
#include "include/ts.h"
#include "include/jssp.h"


#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <random>
#include <filesystem>
#include <ctime>

#include "Logger.h"


/**
 * @param config -co
 *               -g popSize maxGen crossRate mutateRate mutateImpact
 *               -o orCount orImpact
 *               -i instances
 *               -f crossover mutate select
 */
void AlgorithmManager::runCombined(const std::vector<std::string> config) {
    bool timeConstraint = false;
    int seconds;
    bool bm = false;
    int populationSize;
    int maxGenerations;
    int numOR;
    float ortoolsImpact;
    float crossoverRate;
    float mutationRate;
    float mutationImpact;
    FunctionFactory functions;
    std::vector<JobShop> instances = {};

    if(!config.empty()) {
        for(int i = 1; i < config.size(); i++) {
            if(config[i] == "-g") {
                populationSize = std::stoi(config[i+1]);
                maxGenerations = std::stoi(config[i+2]);
                crossoverRate = std::stof(config[i+3]);
                mutationRate = std::stof(config[i+4]);
                mutationImpact = std::stof(config[i+5]);
            }
            if(config[i] == "-o") {
                numOR = std::stoi(config[i+1]);
                ortoolsImpact = std::stof(config[i+2]);
            }
            if(config[i] == "-i") {
                instances = InstanceHandler::getJobShopsByString(config[i+1]);
            }
            if(config[i] == "-f") {
                functions = FunctionFactory(config[i+1], config[i+2], config[i+3]);
            }
            if(config[i] == "-bm") {
                bm = true;
            }
            if(config[i] == "-t") {
                timeConstraint = true;
                seconds = std::stoi(config[i+1]);
            }
        }
    }
    else {
        populationSize = getPopulationSize();
        maxGenerations = getGenerationCount();
        numOR = getORCount();
        ortoolsImpact = getORImpact();
        crossoverRate = getCrossoverRate();
        mutationRate = getMutationRate();
        mutationImpact = getMutationImpact();
        functions = FunctionFactory(1);
        instances = InstanceHandler::getJobShops();
    }
    std::vector<std::vector<Schedule>> population = getPopulation(instances, populationSize);

    for(int i = 0; i < instances.size(); i++) {
        GeneticAlgorithm::evaluate(population[i], instances[i]);
        const int firstSolution = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
        Logger logger(firstSolution, instances[i], config);
        if(timeConstraint) {
            std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::high_resolution_clock::now();
            while(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - startTime).count() < seconds) {
                population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
                population[i] = OrtoolsAlgorithm::solveInterval(population[i], instances[i], ortoolsImpact);
                Algorithm::evaluate(population[i], instances[i]);
                int curOpt = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
                logger.log(curOpt);
                if(curOpt == instances[i].bestKnownSolution) {
                    break;
                }
            }
        }
        else {
            for(int o = 0; o < numOR; o++) {
                population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
                population[i] = OrtoolsAlgorithm::solveInterval(population[i], instances[i], ortoolsImpact);
                Algorithm::evaluate(population[i], instances[i]);
                int curOpt = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
                logger.log(curOpt);
                if(curOpt == instances[i].bestKnownSolution) {
                    break;
                }
            }
        }
        GeneticAlgorithm::printSolution(instances[i], population[i], firstSolution);
        if(bm) {
            logger.save();
        }
    }
}





/**
 * @param config -bo
 *               -g popSize maxGen crossRate mutateRate mutateImpact
 *               -o orCount orImpact
 *               -i instances
 *               -f crossover mutate select
 *               -bm
 *               -t secondsToRun
 */
void AlgorithmManager::runBottleneck(std::vector<std::string> config) {
    bool bm = false;
    bool timeConstraint = false;
    int seconds;
    int populationSize;
    int maxGenerations;
    int numOR;
    float ortoolsImpact;
    float crossoverRate;
    float mutationRate;
    float mutationImpact;
    FunctionFactory functions;
    std::vector<JobShop> instances = {};

    if(!config.empty()) {
        for(int i = 1; i < config.size(); i++) {
            if(config[i] == "-g") {
                populationSize = std::stoi(config[i+1]);
                maxGenerations = std::stoi(config[i+2]);
                crossoverRate = std::stof(config[i+3]);
                mutationRate = std::stof(config[i+4]);
                mutationImpact = std::stof(config[i+5]);
            }
            if(config[i] == "-o") {
                numOR = std::stoi(config[i+1]);
                ortoolsImpact = std::stof(config[i+2]);
            }
            if(config[i] == "-i") {
                instances = InstanceHandler::getJobShopsByString(config[i+1]);
            }
            if(config[i] == "-f") {
                functions = FunctionFactory(config[i+1], config[i+2], config[i+3]);
            }
            if(config[i] == "-bm") {
                bm = true;
            }
            if(config[i] == "-t") {
                timeConstraint = true;
                seconds = std::stoi(config[i+1]);
            }

        }
    }
    else {
        populationSize = getPopulationSize();
        maxGenerations = getGenerationCount();
        numOR = getORCount();
        ortoolsImpact = getORImpact();
        crossoverRate = getCrossoverRate();
        mutationRate = getMutationRate();
        mutationImpact = getMutationImpact();
        functions = FunctionFactory(1);
        instances = InstanceHandler::getJobShops();
    }
    std::vector<std::vector<Schedule>> population = getPopulation(instances, populationSize);



    for(int i = 0; i < instances.size(); i++) {
        GeneticAlgorithm::evaluate(population[i], instances[i]);
        const int firstSolution = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
        Logger logger(firstSolution, instances[i], config);

        if(timeConstraint) {
            std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::high_resolution_clock::now();
            while(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - startTime).count() < seconds) {
                population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
                population[i] = OrtoolsAlgorithm::solveBottleneck(population[i], instances[i], ortoolsImpact);
                Algorithm::evaluate(population[i], instances[i]);
                int curOpt = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
                logger.log(curOpt);
                if(curOpt == instances[i].bestKnownSolution) {
                    break;
                }
            }
        }else {
            for(int o = 0; o < numOR; o++) {
                population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
                population[i] = OrtoolsAlgorithm::solveBottleneck(population[i], instances[i], ortoolsImpact);
                Algorithm::evaluate(population[i], instances[i]);
                int curOpt = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
                logger.log(curOpt);
                if(curOpt == instances[i].bestKnownSolution) {
                    break;
                }
            }
        }

        GeneticAlgorithm::printSolution(instances[i], population[i], firstSolution);
        if(bm) {
            logger.save();
        }
    }
}


/**
 * @param config -mi
 *               -g popSize maxGen crossRate mutRate mutImpact
 *               -o orCount maxIntervalCount orImapct
 *               -i instances
 *               -f crossover mutate select
 *               -bm
 *               -t secondsToRun
 */
void AlgorithmManager::runMultipleCombined(std::vector<std::string> config) {
    bool timeConstraint = false;
    int seconds;
    bool bm = false;
    int populationSize;
    int maxGenerations;
    int numOR;
    int maxIntervalCount;
    float ortoolsImpact;
    float crossoverRate;
    float mutationRate;
    float mutationImpact;
    FunctionFactory functions;
    std::vector<JobShop> instances = {};

    if(!config.empty()) {
        for(int i = 1; i < config.size(); i++) {
            if(config[i] == "-g") {
                populationSize = std::stoi(config[i+1]);
                maxGenerations = std::stoi(config[i+2]);
                crossoverRate = std::stof(config[i+3]);
                mutationRate = std::stof(config[i+4]);
                mutationImpact = std::stof(config[i+5]);
            }
            if(config[i] == "-o") {
                numOR = std::stoi(config[i+1]);
                maxIntervalCount = std::stoi(config[i+2]);
                ortoolsImpact = std::stof(config[i+3]);

            }
            if(config[i] == "-i") {
                instances = InstanceHandler::getJobShopsByString(config[i+1]);
            }
            if(config[i] == "-f") {
                functions = FunctionFactory(config[i+1], config[i+2], config[i+3]);
            }
            if(config[i] == "-bm") {
                bm = true;
            }
            if(config[i] == "-t") {
                timeConstraint = true;
                seconds = std::stoi(config[i+1]);
            }
        }
    }
    else {
        populationSize = getPopulationSize();
        maxGenerations = getGenerationCount();
        numOR = getORCount();
        maxIntervalCount = getMaxIntervalCount();
        ortoolsImpact = getORImpact();
        crossoverRate = getCrossoverRate();
        mutationRate = getMutationRate();
        mutationImpact = getMutationImpact();
        functions = FunctionFactory(1);
        instances = InstanceHandler::getJobShops();
    }
    std::vector<std::vector<Schedule>> population = getPopulation(instances, populationSize);

    for(int i = 0; i < instances.size(); i++) {
        GeneticAlgorithm::evaluate(population[i], instances[i]);
        const int firstSolution = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
        Logger logger(firstSolution, instances[i], config);

        if(timeConstraint) {
            std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::high_resolution_clock::now();
            while(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - startTime).count() < seconds) {
                population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
                population[i] = OrtoolsAlgorithm::solveMultipleIntervals(population[i], instances[i], ortoolsImpact, maxIntervalCount);
                Algorithm::evaluate(population[i], instances[i]);
                int curOpt = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
                logger.log(curOpt);
                if(curOpt == instances[i].bestKnownSolution) {
                    break;
                }
            }
        }
        else {
            for(int o = 0; o < numOR; o++) {
                population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
                population[i] = OrtoolsAlgorithm::solveMultipleIntervals(population[i], instances[i], ortoolsImpact, maxIntervalCount);
                Algorithm::evaluate(population[i], instances[i]);
                int curOpt = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
                logger.log(curOpt);
                if(curOpt == instances[i].bestKnownSolution) {
                    break;
                }
            }
        }
            GeneticAlgorithm::printSolution(instances[i], population[i], firstSolution);
        if(bm) {
            logger.save();
        }
    }
}





/**
 * @param config -jo
 *               -g popSize maxGen crossRate mutateRate mutateImpact
 *               -o orCount orImpact
 *               -i instances
 *               -f crossover mutate select
 */
void AlgorithmManager::runJob(std::vector<std::string> config) {
    bool timeConstraint = false;
    int seconds;
    bool bm = false;
    int populationSize;
    int maxGenerations;
    int numOR;
    float ortoolsImpact;
    float crossoverRate;
    float mutationRate;
    float mutationImpact;
    FunctionFactory functions;
    std::vector<JobShop> instances = {};

    if(!config.empty()) {
        for(int i = 1; i < config.size(); i++) {
            if(config[i] == "-g") {
                populationSize = std::stoi(config[i+1]);
                maxGenerations = std::stoi(config[i+2]);
                crossoverRate = std::stof(config[i+3]);
                mutationRate = std::stof(config[i+4]);
                mutationImpact = std::stof(config[i+5]);
            }
            if(config[i] == "-o") {
                numOR = std::stoi(config[i+1]);
                ortoolsImpact = std::stof(config[i+2]);

            }
            if(config[i] == "-i") {
                instances = InstanceHandler::getJobShopsByString(config[i+1]);
            }
            if(config[i] == "-f") {
                functions = FunctionFactory(config[i+1], config[i+2], config[i+3]);
            }
            if(config[i] == "-bm") {
                bm = true;
            }
            if(config[i] == "-t") {
                timeConstraint = true;
                seconds = std::stoi(config[i+1]);
            }
        }
    }
    else {
        populationSize = getPopulationSize();
        maxGenerations = getGenerationCount();
        numOR = getORCount();
        ortoolsImpact = getORImpact();
        crossoverRate = getCrossoverRate();
        mutationRate = getMutationRate();
        mutationImpact = getMutationImpact();
        functions = FunctionFactory(1);
        instances = InstanceHandler::getJobShops();
    }
    std::vector<std::vector<Schedule>> population = getPopulation(instances, populationSize);



    for(int i = 0; i < instances.size(); i++) {
        GeneticAlgorithm::evaluate(population[i], instances[i]);
        const int firstSolution = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
        Logger logger(firstSolution, instances[i], config);
        if(timeConstraint) {
            std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::high_resolution_clock::now();
            while(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - startTime).count() < seconds) {
                population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
                population[i] = OrtoolsAlgorithm::solveJobs(population[i], instances[i], ortoolsImpact);
                Algorithm::evaluate(population[i], instances[i]);
                int curOpt = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
                logger.log(curOpt);
                if(curOpt == instances[i].bestKnownSolution) {
                    break;
                }
            }
        }else {
            for(int o = 0; o < numOR; o++) {
                population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
                population[i] = OrtoolsAlgorithm::solveJobs(population[i], instances[i], ortoolsImpact);
                Algorithm::evaluate(population[i], instances[i]);
                int curOpt = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
                logger.log(curOpt);
                if(curOpt == instances[i].bestKnownSolution) {
                    break;
                }
            }
        }
        GeneticAlgorithm::printSolution(instances[i], population[i], firstSolution);
        if(bm) {
            logger.save();
        }
    }



}

void AlgorithmManager::runTest(std::vector<std::string> config) {
    std::vector<JobShop> instances = InstanceHandler::getJobShopsByString("swv17");

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> keyDis(0.0, 1.0);
    std::vector<float> key;
    std::vector<Schedule> population = {{Schedule({ 3,
        13,
        30,
        31,
        36,
        40,
        42,
        44,
        48,
        49,
        29,
        7,
        25,
        37,
        22,
        38,
        49,
        23,
        24,
        26,
        9,
        36,
        10,
        41,
        6,
        30,
        48,
        18,
        15,
        28,
        36,
        18,
        1,
        31,
        20,
        22,
        44,
        4,
        41,
        7,
        23,
        3,
        29,
        23,
        22,
        24,
        34,
        23,
        17,
        4,
        49,
        11,
        46,
        22,
        47,
        33,
        3,
        30,
        32,
        34,
        6,
        49,
        45,
        32,
        40,
        18,
        4,
        16,
        21,
        12,
        22,
        49,
        43,
        4,
        16,
        47,
        17,
        19,
        0,
        21,
        24,
        6,
        0,
        14,
        27,
        12,
        41,
        18,
        34,
        46,
        13,
        49,
        36,
        0,
        13,
        25,
        14,
        7,
        41,
        3,
        46,
        14,
        12,
        13,
        40,
        0,
        19,
        32,
        36,
        23,
        34,
        6,
        12,
        14,
        15,
        32,
        35,
        39,
        46,
        42,
        23,
        49,
        40,
        43,
        22,
        18,
        46,
        47,
        0,
        2,
        8,
        2,
        37,
        32,
        39,
        49,
        22,
        5,
        41,
        6,
        37,
        48,
        35,
        45,
        7,
        2,
        15,
        35,
        47,
        28,
        3,
        18,
        43,
        24,
        45,
        47,
        3,
        9,
        28,
        10,
        28,
        29,
        18,
        30,
        35,
        37,
        35,
        13,
        9,
        34,
        9,
        1,
        42,
        47,
        35,
        9,
        43,
        13,
        23,
        48,
        15,
        31,
        42,
        37,
        12,
        19,
        9,
        48,
        17,
        46,
        45,
        27,
        42,
        24,
        15,
        30,
        7,
        48,
        36,
        1,
        39,
        3,
        13,
        45,
        1,
        12,
        42,
        14,
        48,
        29,
        25,
        17,
        44,
        33,
        44,
        5,
        10,
        38,
        46,
        3,
        9,
        41,
        16,
        39,
        40,
        0,
        17,
        44,
        3,
        9,
        8,
        31,
        4,
        29,
        40,
        43,
        5,
        10,
        49,
        38,
        9,
        19,
        33,
        14,
        5,
        43,
        21,
        6,
        28,
        38,
        36,
        5,
        47,
        10,
        25,
        23,
        32,
        43,
        30,
        35,
        32,
        10,
        21,
        33,
        19,
        38,
        30,
        34,
        44,
        14,
        28,
        43,
        25,
        41,
        15,
        12,
        4,
        27,
        38,
        30,
        25,
        29,
        32,
        5,
        44,
        4,
        20,
        24,
        7,
        17,
        39,
        33,
        43,
        31,
        5,
        40,
        45,
        40,
        34,
        10,
        20,
        25,
        31,
        24,
        29,
        2,
        44,
        24,
        13,
        21,
        47,
        41,
        37,
        20,
        45,
        1,
        15,
        44,
        26,
        24,
        29,
        32,
        46,
        48,
        14,
        42,
        2,
        31,
        11,
        20,
        37,
        44,
        42,
        29,
        1,
        26,
        24,
        37,
        26,
        15,
        20,
        45,
        5,
        31,
        11,
        2,
        17,
        18,
        36,
        49,
        46,
        17,
        25,
        19,
        20,
        16,
        11,
        12,
        11,
        26,
        31,
        11,
        18,
        2,
        1,
        22,
        16,
        22,
        8,
        27,
        25,
        46,
        18,
        20,
        10,
        20,
        19,
        12,
        16,
        39,
        21,
        38,
        48,
        6,
        16,
        7,
        8,
        42,
        36,
        39,
        17,
        39,
        8,
        19,
        21,
        27,
        8,
        16,
        26,
        8,
        23,
        7,
        19,
        28,
        7,
        6,
        0,
        26,
        40,
        34,
        27,
        4,
        26,
        36,
        20,
        15,
        27,
        8,
        11,
        26,
        23,
        27,
        7,
        35,
        10,
        34,
        45,
        26,
        4,
        14,
        35,
        37,
        16,
        33,
        22,
        8,
        11,
        6,
        45,
        37,
        17,
        5,
        28,
        28,
        34,
        27,
        31,
        21,
        8,
        39,
        0,
        38,
        11,
        12,
        41,
        47,
        33,
        40,
        0,
        15,
        0,
        35,
        16,
        6,
        38,
        14,
        2,
        10,
        9,
        47,
        33,
        5,
        27,
        30,
        3,
        1,
        19,
        43,
        2,
        33,
        25,
        29,
        4,
        41,
        1,
        38,
        13,
        42,
        21,
        11,
        28,
        32,
        13,
        39,
        48,
        2,
        33,
        30,
        1,
        21})}};
    GeneticAlgorithm::evaluate(population, instances[0]);
    std::cout << population[0].fitness << std::endl;
}

void AlgorithmManager::runMemetic(std::vector<std::string> config) {
    int populationSize = 20;
    bool timeConstraint = false;
    int seconds;
    int tabuSearchIteration = 100;
    bool bm = false;
    int maxGenerations = 100;
    float crossoverRate = 0.7;
    float mutationRate = 0.7;
    float mutationImpact = 0.5;
    FunctionFactory functions = FunctionFactory("or", "or", "rs");
    std::vector<JobShop> instances = InstanceHandler::getJobShopsByString("abz5");

    if(!config.empty()) {
        for(int i = 1; i < config.size(); i++) {
            if(config[i] == "-g") {
                populationSize = std::stoi(config[i+1]);
                maxGenerations = std::stoi(config[i+2]);
                crossoverRate = std::stof(config[i+3]);
                mutationRate = std::stof(config[i+4]);
                mutationImpact = std::stof(config[i+5]);
            }
            if(config[i] == "-i") {
                instances = InstanceHandler::getJobShopsByString(config[i+1]);
            }
            if(config[i] == "-f") {
                functions = FunctionFactory(config[i+1], config[i+2], config[i+3]);
            }
            if(config[i] == "-bm") {
                bm = true;
            }
            if(config[i] == "-ts") {
                tabuSearchIteration = std::stoi(config[i+1]);
            }
            if(config[i] == "-t") {
                timeConstraint = true;
                seconds = std::stoi(config[i+1]);
            }
        }
    }
    else {
        populationSize = getPopulationSize();
        maxGenerations = getGenerationCount();
        crossoverRate = getCrossoverRate();
        mutationRate = getMutationRate();
        mutationImpact = getMutationImpact();
        functions = FunctionFactory(1);
        instances = InstanceHandler::getJobShops();
    }
    std::vector<std::vector<Schedule>> population = getPopulation(instances, populationSize);



    for(int i = 0; i < instances.size(); i++) {
        GeneticAlgorithm::evaluate(population[i], instances[i]);
        int firstSolution = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
        Logger logger(firstSolution, instances[i], config);
        JSSPInstance instance = JobShopToJSSPInstance(instances[i]);
        std::cout << instances[i].name << std::endl;

        if(timeConstraint) {
            std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::high_resolution_clock::now();
            for(int j = 0; j < populationSize; j++) {
                TabuSearch ts = TabuSearch(instance);
                Solution sol = scheduleToMachineSchedule(population[i][j], instances[i]);
                sol = ts.optimize_it(sol,tabuSearchIteration);
                population[i][j] = machineScheduleToSchedule(sol, instances[i]);
            }
            Algorithm::evaluate(population[i], instances[i]);
            Schedule bestSchedule = GeneticAlgorithm::findBestSchedule(population[i]);
            logger.log(bestSchedule.fitness);
            logger.logSchedule(bestSchedule);
            population[i] = GeneticAlgorithm::solveByTimeAndMemetic(population[i], instances[i], seconds, crossoverRate, mutationRate, mutationImpact, functions, logger, startTime, tabuSearchIteration);
            //GeneticAlgorithm::printSolutionWithTime(instances[i], population[i], firstSolution, startTime);
        }
        else {
            for(int j = 0; j < populationSize; j++) {
                TabuSearch ts = TabuSearch(instance);
                Solution sol = scheduleToMachineSchedule(population[i][j], instances[i]);
                sol = ts.optimize_it(sol,tabuSearchIteration);
                population[i][j] = machineScheduleToSchedule(sol, instances[i]);
            }
            Algorithm::evaluate(population[i], instances[i]);
            Schedule bestSchedule = GeneticAlgorithm::findBestSchedule(population[i]);
            logger.log(bestSchedule.fitness);
            logger.logSchedule(bestSchedule);
            population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
            GeneticAlgorithm::printSolution(instances[i], population[i], firstSolution);
        }
        if(bm) {
            logger.save();
        }
    }
}


/**
 *
 * @param config -or
 *               -i instances
 */
void AlgorithmManager::runOR(std::vector<std::string> config) {
    bool bm = false;
    int seconds = 0;
    std::vector<JobShop> instances;
    if(!config.empty()) {
        for(size_t i = 1; i < config.size(); i++) {
            if(config[i] == "-i") {
                instances = InstanceHandler::getJobShopsByString(config[i+1]);
            }
            else if(config[i] == "-bm") {
                bm = true;
            }
            else if(config[i] == "-t") {
                seconds = std::stoi(config[i+1]);
            }
        }
    }
    else {
        instances = InstanceHandler::getJobShops();
    }
    for(auto i : instances) {
        Logger logger(i, config);
        Schedule finalSchedule = OrtoolsAlgorithm::solveAll(i, logger, seconds);
        if(bm) {
            if(finalSchedule.fitness == -1) continue;
            logger.logSchedule(finalSchedule);
            logger.save();
        }
    }
}




/**
 * @param config -ga
 *               -g popSize maxGen crossRate mutateRate mutateImpact
 *               -i instances
 *               -f crossover mutate select
 */
void AlgorithmManager::runGA(std::vector<std::string> config) {
    bool bm = false;
    bool timeConstraint = false;
    int seconds;
    int populationSize;
    int maxGenerations;
    float crossoverRate;
    float mutationRate;
    float mutationImpact;
    FunctionFactory functions;
    std::vector<JobShop> instances = {};

    if(!config.empty()) {
        for(int i = 1; i < config.size(); i++) {
            if(config[i] == "-g") {
                populationSize = std::stoi(config[i+1]);
                maxGenerations = std::stoi(config[i+2]);
                crossoverRate = std::stof(config[i+3]);
                mutationRate = std::stof(config[i+4]);
                mutationImpact = std::stof(config[i+5]);
            }
            if(config[i] == "-i") {
                instances = InstanceHandler::getJobShopsByString(config[i+1]);
            }
            if(config[i] == "-f") {
                functions = FunctionFactory(config[i+1], config[i+2], config[i+3]);
            }
            if(config[i] == "-bm") {
                bm = true;
            }
            if(config[i] == "-t") {
                timeConstraint = true;
                seconds = std::stoi(config[i+1]);
            }
        }
    }
    else {
        populationSize = getPopulationSize();
        maxGenerations = getGenerationCount();
        crossoverRate = getCrossoverRate();
        mutationRate = getMutationRate();
        mutationImpact = getMutationImpact();
        functions = FunctionFactory(1);
        instances = InstanceHandler::getJobShops();
    }
    std::vector<std::vector<Schedule>> population = getPopulation(instances, populationSize);



    for(int i = 0; i < instances.size(); i++) {
        GeneticAlgorithm::evaluate(population[i], instances[i]);
        int firstSolution = GeneticAlgorithm::findBestSchedule(population[i]).fitness;
        Logger logger(firstSolution, instances[i], config);
        std::cout << instances[i].name << std::endl;

        if(timeConstraint) {
            std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::high_resolution_clock::now();
            population[i] = GeneticAlgorithm::solveByTime(population[i], instances[i], seconds, crossoverRate, mutationRate, mutationImpact, functions, logger, startTime);
            Schedule bestSchedule = GeneticAlgorithm::findBestSchedule(population[i]);
            logger.log(bestSchedule.fitness);
            logger.logSchedule(bestSchedule);
        }
        else {

            population[i] = GeneticAlgorithm::solve(population[i], instances[i], maxGenerations, crossoverRate, mutationRate, mutationImpact, functions, logger);
            GeneticAlgorithm::printSolution(instances[i], population[i], firstSolution);
            Schedule bestSchedule = GeneticAlgorithm::findBestSchedule(population[i]);
            logger.log(bestSchedule.fitness);
            logger.logSchedule(bestSchedule);
        }
        if(bm) {
            logger.save();
        }
    }
}


/**
 * const std::vector<float>& key    : random key encoded schedule
 * int jobs                         : count of jobs in the JobShop
 *
 * converting the random key encoded schedule into a working procedure encoded schedule
 */
std::vector<int> AlgorithmManager::decode(const std::vector<float>& key, int jobs){

    //copy of the original schedule (random key encoded) in order to sort it in the next step
    std::vector<float> sortedKey = key;
    std::sort(sortedKey.begin(), sortedKey.end());

    //map will be filled with floats and their respective position in the sorted list
    std::map<float, int> map;
    for(int i = 0; i < sortedKey.size(); i++){
        map[sortedKey[i]] = i;
    }

    //combine the original float list (and its order) with the positions of the floats in the sorted list
    //calculating the modulo of those values
    std::vector<int> decoded;
    for(const float& val : key){
        decoded.push_back(map[val] % jobs); // NOLINT(*-inefficient-vector-operation)
    }

    return decoded;
}

int AlgorithmManager::getPopulationSize() {
    std::string userInput;
    std::cout << "WHICH POPULATION SIZE SHOULD BE USED? \n100 ~standard" << std::endl;
    int populationSize;
    while(true){
        std::getline(std::cin, userInput);
        if(userInput.empty()){
            populationSize = 100;
            break;
        }
        try {
            int num = std::stoi(userInput);
            populationSize = num;
            break;
        } catch (std::invalid_argument& e) {
            std::cout << "INVALID INPUT" << std::endl;
        } catch (std::out_of_range& e) {
            std::cout << "INVALID INPUT" << std::endl;
        }
    }
    return populationSize;
}

int AlgorithmManager::getGenerationCount() {
    std::string userInput;
    std::cout << "HOW MANY GENERATION SHOULD BE BETWEEN EACH OR-TOOLS PROCESS?\n100 ~standard" << std::endl;
    int maxGenerations;
    while(true){
        std::getline(std::cin, userInput);
        if(userInput.empty()){
            maxGenerations = 100;
            break;
        }
        try {
            int num = std::stoi(userInput);
            maxGenerations = num;
            break;
        } catch (std::invalid_argument&) {
            std::cout << "INVALID INPUT" << std::endl;
        } catch (std::out_of_range&) {
            std::cout << "INVALID INPUT" << std::endl;
        }
    }
    return maxGenerations;
}

int AlgorithmManager::getORCount() {
    std::cout << "HOW MANY TIMES SHOULD THE OR-TOOLS SOLVER BE SCHEDULED? \n 20 ~standard" << std::endl;

    std::string userInput;
    int numOR;
    while(true){
        std::getline(std::cin, userInput);
        if(userInput.empty()){
            numOR = 20;
            break;
        }
        try {
            int num = std::stoi(userInput);
            numOR = num;
            break;
        } catch (std::invalid_argument&) {
            std::cout << "INVALID INPUT" << std::endl;
        } catch (std::out_of_range&) {
            std::cout << "INVALID INPUT" << std::endl;
        }
    }
    return numOR;
}

float AlgorithmManager::getORImpact() {
    std::cout << "HOW BIG SHOULD THE OR-TOOLS CHANGE BE?\n0.2 ~standard" << std::endl;

    std::string userInput;
    float ortoolsImpact;
    while(true){
        std::getline(std::cin, userInput);
        if(userInput.empty()){
            ortoolsImpact = 0.2;
            break;
        }
        try {
            const float num = std::stof(userInput);
            ortoolsImpact = num;
            break;
        } catch (std::invalid_argument&) {
            std::cout << "INVALID INPUT" << std::endl;
        } catch (std::out_of_range&) {
            std::cout << "INVALID INPUT" << std::endl;
        }
    }
    return ortoolsImpact;
}

float AlgorithmManager::getCrossoverRate() {
    std::cout << "WHICH CROSSOVER RATE SHOULD BE USED?\n0.7 ~standard" << std::endl;

    std::string userInput;
    float crossoverRate;

    while(true){
        std::getline(std::cin, userInput);
        if(userInput.empty()){
            crossoverRate = 0.7;
            break;
        }
        try {
            const float num = std::stof(userInput);
            crossoverRate = num;
            break;
        } catch (std::invalid_argument&) {
            std::cout << "INVALID INPUT" << std::endl;
        } catch (std::out_of_range&) {
            std::cout << "INVALID INPUT" << std::endl;
        }
    }
    return crossoverRate;
}

float AlgorithmManager::getMutationRate() {
    std::cout << "WHICH MUTATION RATE SHOULD BE USED?\n1.0 ~standard" << std::endl;

    std::string userInput;
    float mutationRate;

    while(true){
        std::getline(std::cin, userInput);
        if(userInput.empty()){
            mutationRate = 1.0;
            break;
        }
        try {
            const float num = std::stof(userInput);
            if(num < 0.0 || num > 1.0){
                std::cout << "INVALID INPUT" << std::endl;
            }
            else{
                mutationRate = num;
                break;
            }
        } catch (std::invalid_argument&) {
            std::cout << "INVALID INPUT" << std::endl;
        } catch (std::out_of_range&) {
            std::cout << "INVALID INPUT" << std::endl;
        }
    }
    return mutationRate;
}

float AlgorithmManager::getMutationImpact() {

    std::cout << "WHICH MUTATION IMPACT SHOULD BE USED?\n0.5 ~standard" << std::endl;

    std::string userInput;
    float mutationImpact;

    while(true){
        std::getline(std::cin, userInput);
        if(userInput.empty()){
            mutationImpact = 0.5;
            break;
        }
        try {
            const float num = std::stof(userInput);
            mutationImpact = num;
            break;
        } catch (std::invalid_argument&) {
            std::cout << "INVALID INPUT" << std::endl;
        } catch (std::out_of_range&) {
            std::cout << "INVALID INPUT" << std::endl;
        }
    }
    return mutationImpact;
}

int AlgorithmManager::getMaxIntervalCount() {
    std::cout << "HOW MANY INTERVALS SHOULD BE POSSIBLE?\n3 ~standard" << std::endl;

    std::string userInput;
    int maxIntervalCount;
    while(true){
        std::getline(std::cin, userInput);
        if(userInput.empty()){
            maxIntervalCount = 3;
            break;
        }
        try {
            const int num = std::stoi(userInput);
            maxIntervalCount = num;
            break;
        } catch (std::invalid_argument&) {
            std::cout << "INVALID INPUT" << std::endl;
        } catch (std::out_of_range&) {
            std::cout << "INVALID INPUT" << std::endl;
        }
    }
    return maxIntervalCount;
}

std::vector<std::vector<Schedule>> AlgorithmManager::getPopulation(const std::vector<JobShop> &instances, const int populationSize) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> keyDis(0.0, 1.0);
    std::vector<float> key;
    std::vector<std::vector<Schedule>> population = {};

    for(int i = 0; i < instances.size(); i++) {
        population.push_back({});
        for(int j = 0; j < populationSize; j++){
            while (key.size() != instances[i].jobCount * instances[i].machineCount) {
                float randomNum = keyDis(gen);
                if(std::find(key.begin(), key.end(), randomNum) == key.end()){
                    key.push_back(randomNum);
                }
            }
            population.back().emplace_back(Schedule(decode(key, instances[i].jobCount)));
            key.clear();
        }
    }
    return population;
}


Solution AlgorithmManager::scheduleToMachineSchedule(Schedule schedule, JobShop &js) {
    std::vector<std::vector<int>> machineSchedule(js.machineCount);
    std::vector<int> currentTaskIndex(js.jobCount, 0);
    for(int i = 0; i < schedule.permutation.size(); i++) {
        int jobIndex = schedule.permutation[i];
        int taskIndex = currentTaskIndex[jobIndex];
        int machine = js.jobs[jobIndex].tasks[taskIndex].machine;
        machineSchedule[machine].emplace_back(jobIndex);
        currentTaskIndex[jobIndex]++;
    }
    return {machineSchedule, schedule.fitness};
}


Schedule AlgorithmManager::machineScheduleToSchedule(
    Solution machineSchedule, JobShop js) {

    std::vector<int> currentJobIndex(js.jobCount, 0), jobTime(js.jobCount, 0), currentMachineIndex(js.machineCount, 0), machineTime(js.machineCount, 0) ;
    std::vector<int> permutation = {};

    int taskCount = js.machineCount * js.jobCount;

    int errorDetection = 0;
    while (taskCount > 0) {
        errorDetection++;
        for(int i = 0; i < js.machineCount; i++) {
            if(currentMachineIndex[i] >= js.jobCount) continue;
            int job = machineSchedule.solution[i][currentMachineIndex[i]];
            if(js.jobs[job].tasks[currentJobIndex[job]].machine == i) {
                const int duration = js.jobs[job].tasks[currentJobIndex[job]].duration;
                const int startTime = std::max(machineTime[i], jobTime[job]);
                jobTime[job] = startTime + duration;
                machineTime[i] = startTime + duration;
                permutation.push_back(job);
                currentMachineIndex[i]++;
                currentJobIndex[job]++;
                taskCount--;
                errorDetection = 0;
            }

        }
        if(errorDetection > taskCount) {
            std::cout << "Error detected in schedule (TabuSearch)..." << std::endl;
            exit(1);
        }
    }
    const int makespan = *std::max_element(machineTime.begin(), machineTime.end());
    auto schedule = Schedule(permutation, makespan);
    return schedule;


}

JSSPInstance AlgorithmManager::JobShopToJSSPInstance(JobShop &js) {
    std::vector<std::vector<Operation>> operations(js.jobCount);
    for(int i = 0; i < js.jobCount; i++) {
        for(auto task : js.jobs[i].tasks) {
            int machine = task.machine;
            int duration = task.duration;
            int job = task.job;
            operations[i].emplace_back(Operation({machine, duration, job}));
        }
    }
    return JSSPInstance(operations, js.jobCount, js.machineCount);
}
