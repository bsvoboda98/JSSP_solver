//
// Created by svobee on 19.04.24.
//

#ifndef ORTOOLSALGORITHM_H
#define ORTOOLSALGORITHM_H

#include "Algorithm.h"
#include "Structs.h"
#include "Logger.h"
#include <ortools/sat/cp_model.h>
#include <utility>
#include <functional>


class OrtoolsAlgorithm : public Algorithm {
public:
    OrtoolsAlgorithm();
    static std::vector<Schedule> solve(std::vector<Schedule> population, JobShop js);


    static Schedule solveAll(JobShop& js, Logger& logger, int seconds);
    static std::vector<Schedule> solveInterval(std::vector<Schedule>& schedules, JobShop& js, float impact);
    static std::vector<Schedule> solveBottleneck(std::vector<Schedule>& schedules, JobShop& js, float impact);
    static std::vector<Schedule> solveMultipleIntervals(std::vector<Schedule>& schedules, JobShop& js, float impact, int maxIntervals);
    //static std::vector<Schedule> solveOuterIntervals(std::vector<Schedule>& schedules, JobShop& js, float impact);
    static std::vector<Schedule> solveJobs(std::vector<Schedule>& schedules, JobShop& js, float impact);

    static std::function<void(std::vector<Schedule>&, float, float, JobShop&, Logger&)> getORMutation();
    static std::function<void(std::vector<Schedule>&, float, float, JobShop&, Logger&)> getORMutationInterval();
    static std::function<void(std::vector<Schedule>&, float, float, JobShop&, Logger&)> getORMutationBottleneck();
    static std::function<void(std::vector<Schedule>&, float, JobShop&, Logger&)> getORCrossover();




private:
    static void collectData(Schedule& schedule, JobShop& js);
    static std::vector<int> getPermutation(operations_research::sat::CpSolverResponse& response, std::vector<std::vector<operations_research::sat::IntervalVar>>& jobIntervals, std::vector<std::pair<int, int>> tasks);
    static int collectDataAndBottleneck(Schedule& schedule, JobShop& js, float impact);
    static std::vector<std::pair<int, bool>> findLongestCommonSequence(const std::vector<int> & vector, const std::vector<int> & parent2);

    static std::vector<std::pair<int, bool>> findLongestCommonMachineSequences(const std::vector<int> &parent1, const std::vector<int> &parent2, JobShop &js);
    static std::vector<std::vector<int>> scheduleToMachineSchedule(std::vector<int> schedule, JobShop &js);
    static std::vector<std::pair<int, bool>> machineScheduleToSchedule(std::vector<std::vector<std::pair<int,bool>>> machineSchedule, JobShop js);
};

#endif //ORTOOLSALGORITHM_H
