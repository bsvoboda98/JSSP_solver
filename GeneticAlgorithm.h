//
// Created by Bennet on 24.01.2024.
//

#ifndef GENETICALGORIHTM_H
#define GENETICALGORIHTM_H

#include <utility>
#include <vector>
#include "Structs.h"
#include "Logger.h"
#include <functional>

#include "Algorithm.h"
#include "FunctionFactory.h"
#include "InstanceHandler.h"
#include "include/jssp.h"

class GeneticAlgorithm : public Algorithm{
    public:

        static std::vector<Schedule> solve(std::vector<Schedule> population, JobShop js,  int maxGenerations,  float crossoverRate, float mutationRate, float mutationImpact, FunctionFactory functions, Logger&);
        static std::vector<Schedule> solveByTimeAndMemetic(std::vector<Schedule> population, JobShop js, int seconds, float crossoverRate, float mutationRate, float mutationImpact, FunctionFactory functions, Logger &logger, std::chrono::time_point<std::chrono::system_clock> &startTime, int tabuSearchIterations);
        static Schedule applyMemetic(Schedule schedule, JobShop &js, int iterations);
        static std::vector<Schedule> solveByTime(std::vector<Schedule> population, JobShop js,  int seconds,  float crossoverRate, float mutationRate, float mutationImpact, FunctionFactory functions, Logger&, std::chrono::time_point<std::chrono::system_clock>& startTime);
        static Schedule findBestSchedule(const std::vector<Schedule>& population);
        static void printSolution(JobShop& js, std::vector<Schedule>& s, int firstSolution);
        static void printSolutionWithTime(JobShop& js, std::vector<Schedule>& s, int firstSolution, std::chrono::time_point<std::chrono::system_clock>& startTime);

        static Solution scheduleToMachineSchedule(Schedule schedule, JobShop &js);

        static Schedule machineScheduleToSchedule(Solution machineSchedule, JobShop js);

        static JSSPInstance JobShopToJSSPInstance(JobShop &js);
};

#endif //GENETICALGORIHTM_H















