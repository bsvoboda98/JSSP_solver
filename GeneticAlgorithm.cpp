//
// Created by Bennet on 24.01.2024.
//

#include "GeneticAlgorithm.h"
#include "Structs.h"
#include "include/ts.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>

/**
 * This function solves the instances that are transferred. To do this, it creates a population for each instance with initialize().
 * The sequence of functions consisting of crossover, mutate, select and evaluate ist then executed for each generation.
 * At the end, the solution is visualized in the console using printSolution().
 *
 * @param js   : vector of JobShops containing the instance for the evaluation steps
 * @param population    : vector of initial population
 * @param maxGenerations : number of Generations to perfome
 * @param crossoverRate : relative quantity of crossovers
 * @param mutationRate : relative quantity of mutations
 * @param mutationImpact : impact of mutations
 * @param functions : collection of functions (crossover, mutate and select) which will be used in this solver
 */
std::vector<Schedule> GeneticAlgorithm::solve(std::vector<Schedule> population, JobShop js, const int maxGenerations,
    const float crossoverRate, const float mutationRate, const float mutationImpact, const FunctionFactory functions, Logger &logger) {



    evaluate(population, js);
    const int populationSize = population.size();
    int curOpt = findBestSchedule(population).fitness;

    for(int i = 0; i < maxGenerations; i++){
        functions.crossover(population, crossoverRate, js, logger); //FunctionFactory
        functions.mutate(population, mutationRate, mutationImpact, js, logger);//FunctionFactory
        evaluate(population, js);
        population = functions.select(population, populationSize);//FunctionFactory
        curOpt = findBestSchedule(population).fitness;
        logger.log(curOpt);
        if(curOpt == js.bestKnownSolution) {
            break;
        }


    }
    return population;
}

std::vector<Schedule> GeneticAlgorithm::solveByTimeAndMemetic(std::vector<Schedule> population, JobShop js, int seconds,
    float crossoverRate, float mutationRate, float mutationImpact, FunctionFactory functions, Logger &logger, std::chrono::time_point<std::chrono::system_clock> &startTime, int tabuSearchIterations) {

    evaluate(population, js);
    const int populationSize = population.size();
    int curOpt = findBestSchedule(population).fitness;

    int noChangeCount = 0;
    int stuckAt = 0;
    float modifiedMutationImpact = mutationImpact;

    while(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - startTime).count() < seconds) {
        functions.crossover(population, crossoverRate, js, logger); //FunctionFactory
        if(population.size() > populationSize) {
            for(int i = populationSize; i < population.size(); i++) {
                population[i] = applyMemetic(population[i], js, tabuSearchIterations);
            }
        }
        else {
            for(int i = 0; i < population.size(); i++) {
                population[i] = applyMemetic(population[i], js, tabuSearchIterations);
            }
        }
        functions.mutate(population, mutationRate, modifiedMutationImpact, js, logger);//FunctionFactory
        evaluate(population, js);
        population = functions.select(population, populationSize);//FunctionFactory
        curOpt = findBestSchedule(population).fitness;
        logger.log(curOpt);
        if(curOpt == js.bestKnownSolution) {
            break;
        }
        if(stuckAt == curOpt) {
            noChangeCount++;
            if(noChangeCount > 3) {
                modifiedMutationImpact+=0.05;
                if(modifiedMutationImpact > 0.8) {
                    modifiedMutationImpact = 0.8;
                }
                noChangeCount = 0;
            }
        }
        else{
            noChangeCount = 0;
            modifiedMutationImpact = mutationImpact;
            stuckAt = curOpt;
        }
    }
    return population;
}

Schedule GeneticAlgorithm::applyMemetic(Schedule schedule, JobShop &js, int iterations) {
    JSSPInstance instance = JobShopToJSSPInstance(js);
    TabuSearch ts = TabuSearch(instance);
    Solution sol = scheduleToMachineSchedule(schedule, js);
    sol = ts.optimize_it(sol, iterations);
    return machineScheduleToSchedule(sol, js);
}


std::vector<Schedule> GeneticAlgorithm::solveByTime(std::vector<Schedule> population, JobShop js, int seconds,
                                                    float crossoverRate, float mutationRate, float mutationImpact, FunctionFactory functions, Logger &logger, std::chrono::time_point<std::chrono::system_clock> &startTime) {

    evaluate(population, js);
    const int populationSize = population.size();
    int curOpt = findBestSchedule(population).fitness;

    int noChangeCount = 0;
    int stuckAt = 0;
    float modifiedMutationImpact = mutationImpact;

    while(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - startTime).count() < seconds) {
        functions.crossover(population, crossoverRate, js, logger); //FunctionFactory
        functions.mutate(population, mutationRate, modifiedMutationImpact, js, logger);//FunctionFactory
        evaluate(population, js);
        population = functions.select(population, populationSize);//FunctionFactory
        curOpt = findBestSchedule(population).fitness;
        logger.log(curOpt);
        if(curOpt == js.bestKnownSolution) {
            break;
        }
        if(stuckAt == curOpt) {
            noChangeCount++;
            if(noChangeCount > 3) {
                modifiedMutationImpact+=0.05;
                if(modifiedMutationImpact > 0.8) {
                    modifiedMutationImpact = 0.8;
                }
                noChangeCount = 0;
            }
        }
        else{
            noChangeCount = 0;
            modifiedMutationImpact = mutationImpact;
            stuckAt = curOpt;
        }
    }
    return population;
}


/**
 * This function finds the Schedule with the best fitness of a generation
 *
 * @param population : population to find the best Schedule from
 * @return           : the best Schedule.
 */
Schedule GeneticAlgorithm::findBestSchedule(const std::vector<Schedule>& population) {

    auto bestS = std::min_element(population.begin(), population.end());
    Schedule bestSchedule = *bestS;

    return bestSchedule;
}

/**
 * In this function basic data will be collected for displaying a comparable benchmark.
 *
 * @param js            : JobShop to get data
 * @param s             : population to find the best schedule from
 * @param firstSolution : best fitness of the first generation (random generated)
 */
void GeneticAlgorithm::printSolution(JobShop& js, std::vector<Schedule>& s, int firstSolution) {

    Schedule bestSchedule = findBestSchedule(s);
    int bestFitness = bestSchedule.fitness;
    float initQuality = (float) js.bestKnownSolution / (float)firstSolution;
    float quality = (float) js.bestKnownSolution / (float) bestFitness;
    std::cout << std::defaultfloat << std::setprecision(2)
            << "\n\n\n\n"
            << "Name \t\t\t" << js.name << " \n"
            << "Bks/lb \t\t\t" << js.bestKnownSolution << ", " << js.lowerBound << "\n"
            << "#Jobs\t\t\t" << js.jobs.size() << "\n"
            << "#Machines\t\t" << js.machineCount << "\n"
            << "Init solution\t" << firstSolution << "\n"
            << "Solution\t\t" << bestFitness << "\n"
            << "Init Quality\t" << initQuality << "\n"
            << "Quality\t\t\t" << quality << "\n"
            << "Improvement\t\t" << quality - initQuality << std::defaultfloat << std::setprecision(6) << std::endl;
            std::cout << "Schedule" << std::endl;
            for(int i = 0; i < bestSchedule.permutation.size(); i++) {
                std::cout << bestSchedule.permutation[i] << ",";
            }
            std::cout << std::endl;


}

void GeneticAlgorithm::printSolutionWithTime(JobShop &js, std::vector<Schedule> &s, int firstSolution,
    std::chrono::time_point<std::chrono::system_clock> &startTime) {
    int duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - startTime).count();
    Schedule bestSchedule = findBestSchedule(s);
    int bestFitness = bestSchedule.fitness;
    float initQuality = (float) js.bestKnownSolution / (float)firstSolution;
    float quality = (float) js.bestKnownSolution / (float) bestFitness;
    std::cout << std::defaultfloat << std::setprecision(2)
            << "\n\n\n\n"
            << "Name \t\t\t" << js.name << " \n"
            << "Bks/lb \t\t\t" << js.bestKnownSolution << ", " << js.lowerBound << "\n"
            << "#Jobs\t\t\t" << js.jobs.size() << "\n"
            << "#Machines\t\t" << js.machineCount << "\n"
            << "Init solution\t" << firstSolution << "\n"
            << "Solution\t\t" << bestFitness << "\n"
            << "Init Quality\t" << initQuality << "\n"
            << "Quality\t\t\t" << quality << "\n"
            << "Improvement\t\t" << quality - initQuality << "\n"
            << "Time\t\t\t" << duration << std::defaultfloat << std::setprecision(6) << std::endl;
            std::cout << "Schedule" << std::endl;
            for(int i = 0; i < bestSchedule.permutation.size(); i++) {
                std::cout << bestSchedule.permutation[i] << ",";
            }
            std::cout << std::endl;
}

Solution GeneticAlgorithm::scheduleToMachineSchedule(Schedule schedule, JobShop &js) {
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


Schedule GeneticAlgorithm::machineScheduleToSchedule(
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

JSSPInstance GeneticAlgorithm::JobShopToJSSPInstance(JobShop &js) {
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
