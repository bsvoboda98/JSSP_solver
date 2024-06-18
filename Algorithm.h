//
// Created by svobee on 19.04.24.
//

#ifndef ALGORITHM_H
#define ALGORITHM_H
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include "include/json.hpp"
#include "Structs.h"

class Algorithm {
public:
    /**
     * std::vector<Schedule>& population    : population which should be evaluated (update fitness)
     * JobShop& jobShop                     : jobShop which contains the needed information about jobs and tasks to reconstruct the schedule
     *
     * in this function, the schedule is applied to the JobShop and the fitness of the schedule is determined
     */
    static void evaluate(std::vector<Schedule>& population, JobShop& jobShop) {

        for (auto& schedule : population) {

            std::vector<int> machineTime(jobShop.machineCount, 0); // Tracks the finishing time of the current job for each machine
            std::unordered_map<int, int> jobTaskIndex; // Tracks the next task to be scheduled for each job
            std::unordered_map<int, std::unordered_map<int, int>> endTimes;

            for (int jobId : schedule.permutation) {

                if (jobTaskIndex[jobId] < jobShop.jobs[jobId].tasks.size()) {

                    Task& task = jobShop.jobs[jobId].tasks[jobTaskIndex[jobId]];

                    //calculates earliest startTime + Duration -> endTime of the task
                    endTimes[jobId][jobTaskIndex[jobId]] = std::max(machineTime[task.machine] + task.duration, endTimes[jobId][jobTaskIndex[jobId] - 1] + task.duration);

                    //sets the MachineTime of the machine used (represents until when the machine is occupied) to the endTime of the schedules Task
                    machineTime[task.machine] = endTimes[jobId][jobTaskIndex[jobId]];

                    //increases the task counter of the current job
                    jobTaskIndex[jobId]++;
                }
            }
            // Calculate the fitness as the maximum finishing time across all machines
            schedule.fitness = *std::max_element(machineTime.begin(), machineTime.end());

        }
    }
};

#endif //ALGORITHM_H
