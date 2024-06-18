//
// Created by svobee on 21.04.24.
//
#include <vector>

#include "Structs.h"
#include "include/jssp.h"

#ifndef ALGORITHMMANAGER_H
#define ALGORITHMMANAGER_H
class AlgorithmManager {
public:
    static void runCombined(std::vector<std::string> config);
    static void runGA(std::vector<std::string> config);
    static void runOR(std::vector<std::string> config);
    static void runBottleneck(std::vector<std::string> config);
    static void runMultipleCombined(std::vector<std::string> config);
    static void runJob(std::vector<std::string> config);
    static void runTest(std::vector<std::string> config);
    static void runMemetic(std::vector<std::string> config);
    static std::vector<int> decode(const std::vector<float>& key, int jobs);

private:
    static int getPopulationSize();
    static int getGenerationCount();
    static int getORCount();
    static float getORImpact();
    static float getCrossoverRate();
    static float getMutationRate();
    static float getMutationImpact();
    static int getMaxIntervalCount();
    static std::vector<std::vector<Schedule>> getPopulation(const std::vector<JobShop>& instances, int populationSize);
    static Solution scheduleToMachineSchedule(Schedule schedule, JobShop &js);
    static Schedule machineScheduleToSchedule(Solution machineSchedule, JobShop js);
    static JSSPInstance JobShopToJSSPInstance(JobShop &js);
};
#endif //ALGORITHMMANAGER_H
