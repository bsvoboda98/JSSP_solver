//
// Created by svobee on 19.04.24.
//

#include "OrtoolsAlgorithm.h"
#include "Logger.h"
#include "ortools/base/logging.h"
#include "ortools/util/time_limit.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "ortools/util/sorted_interval_list.h"
#include "include/json.hpp"
#include <unordered_map>
#include <algorithm>
#include <random>
#include <utility>
#include "include/ts.h"


OrtoolsAlgorithm::OrtoolsAlgorithm() {
}

std::vector<Schedule> OrtoolsAlgorithm::solve(std::vector<Schedule> population, JobShop js) {
    std::vector<Schedule> newPopulation;
    for (auto schedule: population) {
        collectData(schedule, js);
    }

    return population;
}


Schedule OrtoolsAlgorithm::solveAll(JobShop &js, Logger &logger, int seconds) {
    js.upperBound = js.bestKnownSolution * 2;

    operations_research::sat::CpModelBuilder cp_model; //cp_model with all tasks and constraints.


    std::vector<operations_research::sat::IntVar> lastJobs; // lastJobs to minimize the makespan (AddMaxEquality)

    //for each job one vector. These vectors contain every task in the correct order as IntervalVar.
    std::vector<std::vector<operations_research::sat::IntervalVar> > jobIntervals(js.jobCount);

    //add every task as IntervalVar in the vector of its job.
    for (int i = 0; i < js.jobCount; i++) {
        //iterate ofer every task of the curret job (i)
        for (int j = 0; j < js.machineCount; j++) {
            const int duration = js.jobs[i].tasks[j].duration; //duration of task
            operations_research::sat::IntVar start = cp_model.NewIntVar({0, js.upperBound}); //start of the interval
            operations_research::sat::IntVar end = cp_model.NewIntVar({0, js.upperBound}); //end of the interval
            //if its the last job (each job has as many jobs as machines) store the endTime as the last job. So the model can minimize the endtime of the last tasks.
            if (j == js.machineCount - 1) {
                lastJobs.push_back(end);
            }
            operations_research::sat::IntervalVar interval = cp_model.NewIntervalVar(start, duration, end);
            //create a interval and link start, end and duration to it.
            jobIntervals[i].push_back(interval); //adding the interval to its jobs vector.
        }
    }

    //add every interval to its tasks machine-vector
    std::vector<std::vector<operations_research::sat::IntervalVar> > machineIntervals(js.machineCount);
    for (int i = 0; i < js.jobCount; i++) {
        for (int j = 0; j < js.machineCount; j++) {
            machineIntervals[js.jobs[i].tasks[j].machine].emplace_back(jobIntervals[i][j]);
        }
    }

    //the intervals of each machine shall not overlap
    for (const auto &interval: machineIntervals) {
        cp_model.AddNoOverlap(interval);
    }

    //the start time of task j+1 should be before the end time of task j. For every job. (keep order of tasks)
    for (int i = 0; i < js.jobCount; i++) {
        for (int j = 0; j < js.machineCount - 1; j++) {
            cp_model.AddGreaterOrEqual(jobIntervals[i][j + 1].StartExpr(), jobIntervals[i][j].EndExpr());
        }
    }

    //define makespan as IntVar. Make it eaqual to the max endtime of the lastJobs.
    operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, js.upperBound});
    cp_model.AddMaxEquality(makespan, lastJobs);


    //let the model minimize the makespan.
    cp_model.Minimize(makespan);

    //build the model
    operations_research::sat::CpModelProto fixedModel = cp_model.Build();

    operations_research::sat::Model model; // Configure the solver
    operations_research::sat::SatParameters parameters;
    parameters.set_enumerate_all_solutions(true);
    model.Add(NewSatParameters(parameters));
    if (seconds > 0) {
        parameters.set_max_time_in_seconds(seconds);
        model.Add(NewSatParameters(parameters));
    }

    // Create an atomic Boolean to stop the solver after a certain condition is met
    std::atomic<bool> stopped(false);
    model.GetOrCreate<operations_research::TimeLimit>()->RegisterExternalBooleanAsLimit(&stopped);
    model.Add(operations_research::sat::NewFeasibleSolutionObserver(
        [&](const operations_research::sat::CpSolverResponse &response) {
            logger.log(operations_research::sat::SolutionIntegerValue(response, makespan));
        }));

    //solve the builded model.
    operations_research::sat::CpSolverResponse response = operations_research::sat::SolveCpModel(fixedModel, &model);


    //print solution.
    /**if(response.status() == operations_research::sat::OPTIMAL) {
        std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
        std::cout << "bks : " << js.bestKnownSolution << std::endl;
    }
    else if(response.status() != operations_research::sat::UNKNOWN) {
        std::cout << "optimal solution not found!" << std::endl;
        std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
        std::cout << "bks : " << js.bestKnownSolution << std::endl;
    }
    else {
        std::cout << "optimal solution not found!" << std::endl;
    }**/
    if(response.status() != operations_research::sat::OPTIMAL && response.status() != operations_research::sat::FEASIBLE) {
        Schedule error = Schedule({0});
        error.fitness = -1;
        return error;
    }
    std::vector<Schedule> finalSchedule = {Schedule(getPermutation(response, jobIntervals, {}))};
    evaluate(finalSchedule, js);
    return finalSchedule[0];
}


std::vector<Schedule> OrtoolsAlgorithm::solveInterval(std::vector<Schedule> &schedules, JobShop &js, float impact) {
    std::vector<Schedule> newSchedules;
    for (auto schedule: schedules) {
        collectData(schedule, js);
        std::random_device rd;
        std::mt19937 gen(rd());
        int intervalSize = schedule.fitness * impact;
        std::uniform_int_distribution<> intervalGenerator(0, schedule.fitness - intervalSize);
        int intervalStart = intervalGenerator(gen);


        //remove start of schedule
        std::vector<Task> interval, afterInterval;
        std::vector<int> currentTaskIndex(js.jobCount, 0), jobTime(js.jobCount, 0), machineTime(js.machineCount, 0);

        const int taskCount = js.jobCount * js.machineCount;

        operations_research::sat::CpModelBuilder cp_model;

        std::vector<operations_research::sat::IntVar> lastJobs;
        operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, js.upperBound});


        std::vector<std::vector<operations_research::sat::IntervalVar> > intervalJobVars(js.jobCount);
        std::vector<std::vector<operations_research::sat::IntervalVar> > intervalMachineVars(js.jobCount);

        std::vector<std::vector<operations_research::sat::IntervalVar> > afterJobVars(js.jobCount);
        std::vector<std::vector<operations_research::sat::IntervalVar> > afterMachineVars(js.jobCount);

        //store starting times of tasks before the interval in order to create the new schedule at the end.
        std::vector<std::pair<int, int> > beforeTasks;

        for (int i = 0; i < taskCount; i++) {
            const int jobIndex = schedule.permutation[i];
            const int nextTask = currentTaskIndex[jobIndex];
            Task task = js.jobs[jobIndex].tasks[nextTask];
            const int machine = task.machine;

            if (task.endTime < intervalStart) {
                machineTime[machine] = task.endTime;
                jobTime[jobIndex] = task.endTime;
                beforeTasks.push_back(std::make_pair(task.startTime, jobIndex));
            } else if (task.startTime < intervalStart + intervalSize) {
                const int duration = task.duration;
                int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
                operations_research::sat::IntVar start = cp_model.NewIntVar({taskLowerBound, js.upperBound});
                operations_research::sat::IntVar end = cp_model.NewIntVar({taskLowerBound, js.upperBound});
                operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);
                if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                    lastJobs.push_back(end);
                }
                intervalJobVars[task.job].push_back(taskInterval);
                intervalMachineVars[task.machine].push_back(taskInterval);
            } else {
                afterInterval.push_back(task);
                const int duration = task.duration;
                int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
                operations_research::sat::IntVar start = cp_model.NewIntVar({0, js.upperBound});
                operations_research::sat::IntVar end = cp_model.NewIntVar({0, js.upperBound});
                if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                    lastJobs.push_back(end);
                }
                operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);
                cp_model.AddGreaterThan(taskInterval.StartExpr(), taskLowerBound);
                afterJobVars[task.job].push_back(taskInterval);
                afterMachineVars[task.machine].push_back(taskInterval);
            }

            currentTaskIndex[jobIndex]++;
        }

        //keep order of interval tasks for every job and make it smaller than the first task of the same job after the interval.
        for (int i = 0; i < js.jobCount; i++) {
            if (!afterJobVars[i].empty() && !intervalJobVars[i].empty()) {
                cp_model.AddGreaterOrEqual(afterJobVars[i][0].StartExpr(),
                                           intervalJobVars[i][intervalJobVars[i].size() - 1].EndExpr());
            }
            for (int j = 0; j < (int) intervalJobVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(intervalJobVars[i][j + 1].StartExpr(), intervalJobVars[i][j].EndExpr());
            }
        }

        //allow no overlap of tasks on the same machine in the intervall. Every task has to be before the first task on this machine after the interval
        for (int i = 0; i < js.machineCount; i++) {
            cp_model.AddNoOverlap(intervalMachineVars[i]);
            if (afterMachineVars[i].empty() || intervalMachineVars[i].empty()) continue;
            for (int j = 0; j < intervalMachineVars[i].size(); j++) {
                cp_model.AddGreaterOrEqual(afterMachineVars[i][0].StartExpr(), intervalMachineVars[i][j].EndExpr());
            }
        }

        //keep the order of every task of each job after the interval
        for (int i = 0; i < js.jobCount; i++) {
            for (int j = 0; j < (int) afterJobVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(afterJobVars[i][j + 1].StartExpr(), afterJobVars[i][j].EndExpr());
            }
        }
        //keep the order of every task on every machine after the interval
        for (int i = 0; i < js.jobCount; i++) {
            for (int j = 0; j < (int) afterMachineVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(afterMachineVars[i][j + 1].StartExpr(), afterMachineVars[i][j].EndExpr());
            }
        }

        cp_model.AddMaxEquality(makespan, lastJobs);
        cp_model.Minimize(makespan);
        //build the model
        operations_research::sat::CpModelProto fixedModel = cp_model.Build();

        //solve the builded model.
        operations_research::sat::CpSolverResponse response = operations_research::sat::Solve(fixedModel);
        //print solution.
        /*if(response.status() == operations_research::sat::OPTIMAL) {
            std::cout << "initSolution" << js.upperBound << std::endl;
            std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
            std::cout << "bks : " << js.bestKnownSolution << std::endl;
        }
        else {
            std::cout << "optimal solution not found!" << std::endl;
        }*/


        //concat the interval and after jobVars in order to create the permutation of the schedule in the next step.
        std::vector<std::vector<operations_research::sat::IntervalVar> > jobInterval(js.jobCount);

        for (int i = 0; i < js.jobCount; i++) {
            intervalJobVars[i].insert(intervalJobVars[i].end(), afterJobVars[i].begin(), afterJobVars[i].end());
        }

        std::vector<int> permutation = getPermutation(response, intervalJobVars, beforeTasks);
        newSchedules.push_back(Schedule(permutation));
    }
    evaluate(newSchedules, js);
    return newSchedules;
}


std::vector<Schedule> OrtoolsAlgorithm::solveBottleneck(std::vector<Schedule> &schedules, JobShop &js, float impact) {
    std::vector<Schedule> newSchedules;
    for (auto schedule: schedules) {
        int bottleneck = collectDataAndBottleneck(schedule, js, impact);
        int intervalSize = schedule.fitness * impact;
        int intervalStart = bottleneck - intervalSize / 2;


        //remove start of schedule
        std::vector<int> currentTaskIndex(js.jobCount, 0), jobTime(js.jobCount, 0), machineTime(js.machineCount, 0);

        const int taskCount = js.jobCount * js.machineCount;

        operations_research::sat::CpModelBuilder cp_model;

        std::vector<operations_research::sat::IntVar> lastJobs;
        operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, js.upperBound});


        std::vector<std::vector<operations_research::sat::IntervalVar> > intervalJobVars(js.jobCount);
        std::vector<std::vector<operations_research::sat::IntervalVar> > intervalMachineVars(js.jobCount);

        std::vector<std::vector<operations_research::sat::IntervalVar> > afterJobVars(js.jobCount);
        std::vector<std::vector<operations_research::sat::IntervalVar> > afterMachineVars(js.jobCount);

        //store starting times of tasks before the interval in order to create the new schedule at the end.
        std::vector<std::pair<int, int> > beforeTasks;

        for (int i = 0; i < taskCount; i++) {
            const int jobIndex = schedule.permutation[i];
            const int nextTask = currentTaskIndex[jobIndex];
            Task task = js.jobs[jobIndex].tasks[nextTask];
            const int machine = task.machine;

            if (task.endTime < intervalStart) {
                machineTime[machine] = task.endTime;
                jobTime[jobIndex] = task.endTime;
                beforeTasks.push_back(std::make_pair(task.startTime, jobIndex));
            } else if (task.startTime < intervalStart + intervalSize) {
                const int duration = task.duration;
                int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
                operations_research::sat::IntVar start = cp_model.NewIntVar({taskLowerBound, js.upperBound});
                operations_research::sat::IntVar end = cp_model.NewIntVar({taskLowerBound, js.upperBound});
                operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);
                if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                    lastJobs.push_back(end);
                }
                intervalJobVars[task.job].push_back(taskInterval);
                intervalMachineVars[task.machine].push_back(taskInterval);
            } else {
                const int duration = task.duration;
                int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
                operations_research::sat::IntVar start = cp_model.NewIntVar({0, js.upperBound});
                operations_research::sat::IntVar end = cp_model.NewIntVar({0, js.upperBound});
                if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                    lastJobs.push_back(end);
                }
                operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);
                cp_model.AddGreaterThan(taskInterval.StartExpr(), taskLowerBound);
                afterJobVars[task.job].push_back(taskInterval);
                afterMachineVars[task.machine].push_back(taskInterval);
            }

            currentTaskIndex[jobIndex]++;
        }

        //keep order of interval tasks for every job and make it smaller than the first task of the same job after the interval.
        for (int i = 0; i < js.jobCount; i++) {
            if (!afterJobVars[i].empty() && !intervalJobVars[i].empty()) {
                cp_model.AddGreaterOrEqual(afterJobVars[i][0].StartExpr(),
                                           intervalJobVars[i][intervalJobVars[i].size() - 1].EndExpr());
            }
            for (int j = 0; j < (int) intervalJobVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(intervalJobVars[i][j + 1].StartExpr(), intervalJobVars[i][j].EndExpr());
            }
        }

        //allow no overlap of tasks on the same machine in the intervall. Every task has to be before the first task on this machine after the interval
        for (int i = 0; i < js.machineCount; i++) {
            cp_model.AddNoOverlap(intervalMachineVars[i]);
            if (afterMachineVars[i].empty() || intervalMachineVars[i].empty()) continue;
            for (int j = 0; j < intervalMachineVars[i].size(); j++) {
                cp_model.AddGreaterOrEqual(afterMachineVars[i][0].StartExpr(), intervalMachineVars[i][j].EndExpr());
            }
        }

        //keep the order of every task of each job after the interval
        for (int i = 0; i < js.jobCount; i++) {
            for (int j = 0; j < (int) afterJobVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(afterJobVars[i][j + 1].StartExpr(), afterJobVars[i][j].EndExpr());
            }
        }
        //keep the order of every task on every machine after the interval
        for (int i = 0; i < js.jobCount; i++) {
            for (int j = 0; j < (int) afterMachineVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(afterMachineVars[i][j + 1].StartExpr(), afterMachineVars[i][j].EndExpr());
            }
        }

        cp_model.AddMaxEquality(makespan, lastJobs);
        cp_model.Minimize(makespan);
        //build the model
        operations_research::sat::CpModelProto fixedModel = cp_model.Build();

        //solve the builded model.
        operations_research::sat::CpSolverResponse response = operations_research::sat::Solve(fixedModel);
        //print solution.


        /*if(response.status() == operations_research::sat::OPTIMAL) {
            std::cout << "initSolution" << js.upperBound << std::endl;
            std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
            std::cout << "bks : " << js.bestKnownSolution << std::endl;
        }
        else {
            std::cout << "optimal solution not found!" << std::endl;
        }*/


        //concat the interval and after jobVars in order to create the permutation of the schedule in the next step.
        std::vector<std::vector<operations_research::sat::IntervalVar> > jobInterval(js.jobCount);

        for (int i = 0; i < js.jobCount; i++) {
            intervalJobVars[i].insert(intervalJobVars[i].end(), afterJobVars[i].begin(), afterJobVars[i].end());
        }

        std::vector<int> permutation = getPermutation(response, intervalJobVars, beforeTasks);
        newSchedules.push_back(Schedule(permutation));
    }
    evaluate(newSchedules, js);
    return newSchedules;
}

std::vector<Schedule> OrtoolsAlgorithm::solveMultipleIntervals(std::vector<Schedule> &schedules, JobShop &js,
                                                               float impact, int maxIntervals) {
    std::vector<Schedule> newSchedules;
    for (auto schedule: schedules) {
        collectData(schedule, js);
        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_int_distribution<> intervalCountGenerator(1, maxIntervals);

        int intervalCount = intervalCountGenerator(gen);

        const int intervalSize = (schedule.fitness * impact) / intervalCount;
        std::uniform_int_distribution<> intervalGenerator(0, schedule.fitness - intervalSize);
        std::vector<int> intervalStarts;

        while (intervalStarts.size() < intervalCount) {
            int start = intervalGenerator(gen);

            bool overlap = false;
            for (const auto &interval: intervalStarts) {
                if (start + intervalSize > interval && start < interval + intervalSize) {
                    overlap = true;
                    break;
                }
            }
            if (!overlap) {
                intervalStarts.push_back(start);
            }
        }

        std::sort(intervalStarts.begin(), intervalStarts.end());


        //remove start of schedule

        std::vector<int> currentTaskIndex(js.jobCount, 0), jobTime(js.jobCount, 0), machineTime(js.machineCount, 0);

        const int taskCount = js.jobCount * js.machineCount;

        operations_research::sat::CpModelBuilder cp_model;

        std::vector<operations_research::sat::IntVar> lastJobs;
        operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, js.upperBound});


        std::vector intervalJobVars(intervalCount,
                                    std::vector<std::vector<operations_research::sat::IntervalVar> >(js.jobCount));
        std::vector intervalMachineVars(intervalCount,
                                        std::vector<std::vector<operations_research::sat::IntervalVar> >(js.jobCount));

        std::vector betweenIntervalJobVars(intervalCount - 1,
                                           std::vector<std::vector<
                                               operations_research::sat::IntervalVar> >(js.jobCount));
        std::vector betweenIntervalMachineVars(intervalCount - 1,
                                               std::vector<std::vector<operations_research::sat::IntervalVar> >(
                                                   js.jobCount));

        std::vector<std::vector<operations_research::sat::IntervalVar> > afterJobVars(js.jobCount);
        std::vector<std::vector<operations_research::sat::IntervalVar> > afterMachineVars(js.jobCount);

        //store starting times of tasks before the interval in order to create the new schedule at the end.
        std::vector<std::pair<int, int> > beforeTasks;

        for (int i = 0; i < taskCount; i++) {
            const int jobIndex = schedule.permutation[i];
            const int nextTask = currentTaskIndex[jobIndex];
            Task task = js.jobs[jobIndex].tasks[nextTask];
            const int machine = task.machine;


            if (task.endTime < intervalStarts[0]) {
                machineTime[machine] = task.endTime;
                jobTime[jobIndex] = task.endTime;
                beforeTasks.push_back(std::make_pair(task.startTime, jobIndex));
                currentTaskIndex[task.job]++;
                continue;
            }

            const int duration = task.duration;
            int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
            operations_research::sat::IntVar start = cp_model.NewIntVar({taskLowerBound, js.upperBound});
            operations_research::sat::IntVar end = cp_model.NewIntVar({taskLowerBound, js.upperBound});
            operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);

            for (int j = 0; j < intervalCount; j++) {
                if (j > 0) {
                    if (task.startTime < intervalStarts[j]) {
                        betweenIntervalJobVars[j - 1][task.job].push_back(taskInterval);
                        betweenIntervalMachineVars[j - 1][task.machine].push_back(taskInterval);
                        if (currentTaskIndex[task.job] == js.machineCount - 1) {
                            lastJobs.push_back(end);
                        }
                        break;
                    }
                }
                if (task.startTime < intervalStarts[j] + intervalSize) {
                    intervalJobVars[j][task.job].push_back(taskInterval);
                    intervalMachineVars[j][task.machine].push_back(taskInterval);
                    if (currentTaskIndex[task.job] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                    break;
                }
            }

            if (task.startTime >= intervalStarts[intervalCount - 1] + intervalSize) {
                if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                    lastJobs.push_back(end);
                }
                cp_model.AddGreaterThan(taskInterval.StartExpr(), taskLowerBound);
                afterJobVars[task.job].push_back(taskInterval);
                afterMachineVars[task.machine].push_back(taskInterval);
            }

            currentTaskIndex[task.job]++;
        }

        //keep order of interval tasks for every job and make it smaller than the first task of the same job after the interval.
        for (int k = 0; k < intervalCount; k++) {
            for (int i = 0; i < js.jobCount; i++) {
                if (k < intervalCount - 1) {
                    if (!betweenIntervalJobVars[k][i].empty() && !intervalJobVars[k][i].empty()) {
                        cp_model.AddGreaterOrEqual(betweenIntervalJobVars[k][i][0].StartExpr(),
                                                   intervalJobVars[k][i][intervalJobVars[k][i].size() - 1].EndExpr());
                    }
                    if (!intervalJobVars[k + 1][i].empty() && !intervalJobVars[k][i].empty()) {
                        cp_model.AddGreaterOrEqual(intervalJobVars[k + 1][i][0].StartExpr(),
                                                   intervalJobVars[k][i][intervalJobVars[k][i].size() - 1].EndExpr());
                    }
                }
                if (k > 0 && !intervalJobVars[k][i].empty()) {
                    int betweenJobSize = betweenIntervalJobVars[k - 1][i].size();
                    if (betweenJobSize > 0) {
                        cp_model.AddGreaterOrEqual(intervalJobVars[k][i][0].StartExpr(),
                                                   betweenIntervalJobVars[k - 1][i][betweenJobSize - 1].EndExpr());
                    } else {
                        int beforeJobSize = intervalJobVars[k - 1][i].size();
                        if (beforeJobSize > 0) {
                            cp_model.AddGreaterOrEqual(intervalJobVars[k][i][0].StartExpr(),
                                                       intervalJobVars[k - 1][i][beforeJobSize - 1].EndExpr());
                        }
                    }
                }
                if (!afterJobVars[i].empty() && !intervalJobVars[k][i].empty()) {
                    cp_model.AddGreaterOrEqual(afterJobVars[i][0].StartExpr(),
                                               intervalJobVars[k][i][intervalJobVars[k][i].size() - 1].EndExpr());
                }
                for (int j = 0; j < (int) intervalJobVars[k][i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(intervalJobVars[k][i][j + 1].StartExpr(),
                                               intervalJobVars[k][i][j].EndExpr());
                }
            }
        }

        //allow no overlap of tasks on the same machine in the intervall. Every task has to be before the first task on this machine after the interval
        for (int k = 0; k < intervalCount; k++) {
            for (int i = 0; i < js.machineCount; i++) {
                cp_model.AddNoOverlap(intervalMachineVars[k][i]);
                if (k < intervalCount - 1) {
                    if (!betweenIntervalMachineVars[k][i].empty()) {
                        for (int m = 0; m < (int) intervalMachineVars[k][i].size(); m++) {
                            cp_model.AddGreaterOrEqual(betweenIntervalMachineVars[k][i][0].StartExpr(),
                                                       intervalMachineVars[k][i][m].EndExpr());
                        }
                    }
                    if (!intervalMachineVars[k + 1][i].empty()) {
                        for (int m = 0; m < (int) intervalMachineVars[k][i].size(); m++) {
                            for (int b = 0; b < (int) betweenIntervalMachineVars[k][i].size(); b++) {
                                cp_model.AddGreaterOrEqual(betweenIntervalMachineVars[k][i][b].StartExpr(),
                                                           intervalMachineVars[k][i][m].EndExpr());
                            }
                        }
                    }
                }
                if (k > 0) {
                    int betweenMachineSize = betweenIntervalMachineVars[k - 1][i].size();
                    if (betweenMachineSize > 0) {
                        for (int j = 0; j < (int) intervalMachineVars[k][i].size(); j++) {
                            cp_model.AddGreaterOrEqual(intervalMachineVars[k][i][j].StartExpr(),
                                                       betweenIntervalMachineVars[k - 1][i][betweenMachineSize - 1].
                                                       EndExpr());
                        }
                    } else {
                        for (int j = 0; j < (int) intervalMachineVars[k][i].size(); j++) {
                            for (int b = 0; b < (int) intervalMachineVars[k - 1][i].size(); b++) {
                                cp_model.AddGreaterOrEqual(intervalMachineVars[k][i][j].StartExpr(),
                                                           intervalMachineVars[k - 1][i][b].EndExpr());
                            }
                        }
                    }
                }
                if (afterMachineVars[i].empty() || intervalMachineVars[k][i].empty()) continue;
                for (int j = 0; j < (int) intervalMachineVars[k][i].size(); j++) {
                    cp_model.AddGreaterOrEqual(afterMachineVars[i][0].StartExpr(),
                                               intervalMachineVars[k][i][j].EndExpr());
                }
            }
        }

        //keep the order of every task of each job between the intervals
        for (int k = 0; k < intervalCount - 1; k++) {
            for (int i = 0; i < js.jobCount; i++) {
                for (int j = 0; j < (int) betweenIntervalJobVars[k][i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(betweenIntervalJobVars[k][i][j + 1].StartExpr(),
                                               betweenIntervalJobVars[k][i][j].EndExpr());
                }
            }
        }

        //keep the order of every task of each machine between the intervals
        for (int k = 0; k < intervalCount - 1; k++) {
            for (int i = 0; i < js.machineCount; i++) {
                for (int j = 0; j < (int) betweenIntervalMachineVars[k][i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(betweenIntervalMachineVars[k][i][j + 1].StartExpr(),
                                               betweenIntervalMachineVars[k][i][j].EndExpr());
                }
            }
        }

        //keep the order of every task of each job after the interval
        for (int i = 0; i < js.jobCount; i++) {
            for (int j = 0; j < (int) afterJobVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(afterJobVars[i][j + 1].StartExpr(), afterJobVars[i][j].EndExpr());
            }
        }
        //keep the order of every task on every machine after the interval
        for (int i = 0; i < js.jobCount; i++) {
            for (int j = 0; j < (int) afterMachineVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(afterMachineVars[i][j + 1].StartExpr(), afterMachineVars[i][j].EndExpr());
            }
        }

        cp_model.AddMaxEquality(makespan, lastJobs);
        cp_model.Minimize(makespan);
        //build the model
        operations_research::sat::CpModelProto fixedModel = cp_model.Build();

        //solve the builded model.
        operations_research::sat::CpSolverResponse response = operations_research::sat::Solve(fixedModel);

        //print solution.
        if (response.status() == operations_research::sat::OPTIMAL) {
            std::cout << "initSolution" << js.upperBound << std::endl;
            std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) <<
                    std::endl;
            std::cout << "bks : " << js.bestKnownSolution << std::endl;
        } else {
            std::cout << "optimal solution not found!" << std::endl;
        }

        std::vector<std::vector<operations_research::sat::IntervalVar> > jobVars(js.jobCount);

        //concat the interval and after jobVars in order to create the permutation of the schedule in the next step.
        for (int k = 0; k < intervalCount; k++) {
            if (k > 0) {
                for (int i = 0; i < js.jobCount; i++) {
                    jobVars[i].insert(jobVars[i].end(), betweenIntervalJobVars[k - 1][i].begin(),
                                      betweenIntervalJobVars[k - 1][i].end());
                }
            }
            for (int i = 0; i < js.jobCount; i++) {
                jobVars[i].insert(jobVars[i].end(), intervalJobVars[k][i].begin(), intervalJobVars[k][i].end());
            }
        }
        for (int i = 0; i < js.jobCount; i++) {
            jobVars[i].insert(jobVars[i].end(), afterJobVars[i].begin(), afterJobVars[i].end());
        }

        std::vector<int> permutation = getPermutation(response, jobVars, beforeTasks);
        newSchedules.push_back(Schedule(permutation));
    }
    evaluate(newSchedules, js);
    return newSchedules;
}

std::vector<Schedule> OrtoolsAlgorithm::solveJobs(std::vector<Schedule> &schedules, JobShop &js, float impact) {
    std::vector<Schedule> newSchedules;
    for (auto schedule: schedules) {
        collectData(schedule, js);
        std::random_device rd;
        std::mt19937 gen(rd());
        int jobSelection = js.jobCount * impact;
        std::uniform_int_distribution<> intervalGenerator(0, js.jobCount);
        std::vector<int> jobsToSolve;
        while (jobsToSolve.size() < jobSelection) {
            int selectedJob = intervalGenerator(gen);
            if (std::find(jobsToSolve.begin(), jobsToSolve.end(), selectedJob) == jobsToSolve.end()) {
                jobsToSolve.push_back(selectedJob);
            }
        }


        //remove start of schedule
        std::vector<Task> interval, afterInterval;
        std::vector<int> currentTaskIndex(js.jobCount, 0), jobTime(js.jobCount, 0), machineTime(js.machineCount, 0);

        const int taskCount = js.jobCount * js.machineCount;

        operations_research::sat::CpModelBuilder cp_model;

        std::vector<operations_research::sat::IntVar> lastJobs;
        operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, js.upperBound});


        std::vector<std::vector<operations_research::sat::IntervalVar> > allJobVars(js.jobCount);
        std::vector<std::vector<operations_research::sat::IntervalVar> > allMachineVars(js.jobCount);

        //std::vector<std::vector<operations_research::sat::IntervalVar>> notSelectedJobVars(js.jobCount);
        std::vector<std::vector<operations_research::sat::IntervalVar> > notSelectedMachineVars(js.jobCount);

        //store starting times of tasks before the interval in order to create the new schedule at the end.
        std::vector<std::pair<int, int> > beforeTasks;

        for (int i = 0; i < taskCount; i++) {
            const int jobIndex = schedule.permutation[i];
            const int nextTask = currentTaskIndex[jobIndex];
            Task task = js.jobs[jobIndex].tasks[nextTask];
            const int machine = task.machine;

            const int duration = task.duration;
            operations_research::sat::IntVar start = cp_model.NewIntVar({0, js.upperBound});
            operations_research::sat::IntVar end = cp_model.NewIntVar({0, js.upperBound});
            operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);


            if (std::find(jobsToSolve.begin(), jobsToSolve.end(), jobIndex) != jobsToSolve.end()) {
                allJobVars[jobIndex].push_back(taskInterval);
                allMachineVars[machine].push_back(taskInterval);
                if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                    lastJobs.push_back(end);
                }
            } else {
                allJobVars[jobIndex].push_back(taskInterval);
                notSelectedMachineVars[machine].push_back(taskInterval);
                allMachineVars[machine].push_back(taskInterval);
                if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                    lastJobs.push_back(end);
                }
            }

            currentTaskIndex[jobIndex]++;
        }

        //keep order of interval tasks for every job and make it smaller than the first task of the same job after the interval.
        for (int i = 0; i < js.jobCount; i++) {
            for (int j = 0; j < js.machineCount - 1; j++) {
                cp_model.AddGreaterOrEqual(allJobVars[i][j + 1].StartExpr(), allJobVars[i][j].EndExpr());
            }
        }

        //allow no overlap of tasks on the same machine in the intervall. Every task has to be before the first task on this machine after the interval
        for (int i = 0; i < js.machineCount; i++) {
            cp_model.AddNoOverlap(allMachineVars[i]);
            for (int j = 0; j < notSelectedMachineVars[i].size() - 1; j++) {
                cp_model.AddGreaterOrEqual(notSelectedMachineVars[i][j + 1].StartExpr(),
                                           notSelectedMachineVars[i][j].EndExpr());
            }
        }

        //keep the order of every task of each job after the interval
        for (int i = 0; i < js.jobCount; i++) {
        }
        //keep the order of every task on every machine after the interval
        for (int i = 0; i < js.jobCount; i++) {
        }

        cp_model.AddMaxEquality(makespan, lastJobs);
        cp_model.Minimize(makespan);
        //build the model
        operations_research::sat::CpModelProto fixedModel = cp_model.Build();

        //solve the builded model.
        operations_research::sat::CpSolverResponse response = operations_research::sat::Solve(fixedModel);
        //print solution.
        if (response.status() == operations_research::sat::OPTIMAL) {
            std::cout << "initSolution" << js.upperBound << std::endl;
            std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) <<
                    std::endl;
            std::cout << "bks : " << js.bestKnownSolution << std::endl;
        } else {
            std::cout << "optimal solution not found!" << std::endl;
        }


        //concat the interval and after jobVars in order to create the permutation of the schedule in the next step.
        std::vector<std::vector<operations_research::sat::IntervalVar> > jobInterval(js.jobCount);


        std::vector<int> permutation = getPermutation(response, allJobVars, {});
        newSchedules.push_back(Schedule(permutation));
    }
    evaluate(newSchedules, js);
    return newSchedules;
}

std::function<void(std::vector<Schedule> &, float, float, JobShop &js, Logger &)> OrtoolsAlgorithm::getORMutation() {
    return [](std::vector<Schedule> &population, float mutationRate, float mutationImpact, JobShop &js,
              Logger &logger) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> mutationGenerator(0, 1);
        std::vector<Schedule> newPop;
        for (int p = 0; p < population.size(); p++) {
            if (mutationGenerator(gen) > mutationRate) continue;
            collectData(population[p], js);
            int jobSelection = js.jobCount * mutationImpact;
            std::uniform_int_distribution<> intervalGenerator(0, js.jobCount);
            std::vector<int> jobsToSolve;
            while (jobsToSolve.size() < jobSelection) {
                int selectedJob = intervalGenerator(gen);
                if (std::find(jobsToSolve.begin(), jobsToSolve.end(), selectedJob) == jobsToSolve.end()) {
                    jobsToSolve.push_back(selectedJob);
                }
            }


            //remove start of schedule
            std::vector<Task> interval, afterInterval;
            std::vector<int> currentTaskIndex(js.jobCount, 0), jobTime(js.jobCount, 0), machineTime(js.machineCount, 0);

            const int taskCount = js.jobCount * js.machineCount;

            operations_research::sat::CpModelBuilder cp_model;

            std::vector<operations_research::sat::IntVar> lastJobs;
            int64_t horizon = 0;
            for(auto job : js.jobs) {
                for(auto task : job.tasks) {
                    horizon += task.duration;
                }
            }
            operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, horizon});


            std::vector<std::vector<operations_research::sat::IntervalVar> > allJobVars(js.jobCount);
            std::vector<std::vector<operations_research::sat::IntervalVar> > allMachineVars(js.jobCount);

            //std::vector<std::vector<operations_research::sat::IntervalVar>> notSelectedJobVars(js.jobCount);
            std::vector<std::vector<operations_research::sat::IntervalVar> > notSelectedMachineVars(js.jobCount);

            //store starting times of tasks before the interval in order to create the new schedule at the end.
            std::vector<std::pair<int, int> > beforeTasks;

            for (int i = 0; i < taskCount; i++) {
                const int jobIndex = population[p].permutation[i];
                const int nextTask = currentTaskIndex[jobIndex];
                Task task = js.jobs[jobIndex].tasks[nextTask];
                const int machine = task.machine;

                const int duration = task.duration;
                operations_research::sat::IntVar start = cp_model.NewIntVar({0, horizon});
                operations_research::sat::IntVar end = cp_model.NewIntVar({0, horizon});
                operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);


                if (std::find(jobsToSolve.begin(), jobsToSolve.end(), jobIndex) != jobsToSolve.end()) {
                    allJobVars[jobIndex].push_back(taskInterval);
                    allMachineVars[machine].push_back(taskInterval);
                    if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                } else {
                    allJobVars[jobIndex].push_back(taskInterval);
                    notSelectedMachineVars[machine].push_back(taskInterval);
                    allMachineVars[machine].push_back(taskInterval);
                    if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                }

                currentTaskIndex[jobIndex]++;
            }

            //keep order of interval tasks for every job and make it smaller than the first task of the same job after the interval.
            for (int i = 0; i < js.jobCount; i++) {
                for (int j = 0; j < js.machineCount - 1; j++) {
                    cp_model.AddGreaterOrEqual(allJobVars[i][j + 1].StartExpr(), allJobVars[i][j].EndExpr());
                }
            }

            //allow no overlap of tasks on the same machine in the intervall. Every task has to be before the first task on this machine after the interval
            for (int i = 0; i < js.machineCount; i++) {
                cp_model.AddNoOverlap(allMachineVars[i]);
                for (int j = 0; j < notSelectedMachineVars[i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(notSelectedMachineVars[i][j + 1].StartExpr(),
                                               notSelectedMachineVars[i][j].EndExpr());
                }
            }


            cp_model.AddMaxEquality(makespan, lastJobs);
            cp_model.Minimize(makespan);
            //build the model
            operations_research::sat::CpModelProto fixedModel = cp_model.Build();

            operations_research::sat::Model model;
            operations_research::sat::SatParameters parameters;
            parameters.set_max_time_in_seconds(25.0);
            model.Add(NewSatParameters(parameters));


            //solve the builded model.
            operations_research::sat::CpSolverResponse response = operations_research::sat::SolveCpModel(
                fixedModel, &model);


            if (response.status() != operations_research::sat::CpSolverStatus::OPTIMAL && response.status() != operations_research::sat::CpSolverStatus::FEASIBLE) {
                continue;
            }
            int curSol = operations_research::sat::SolutionIntegerValue(response, makespan);
            logger.log(curSol);

            //print solution.
            /**if(response.status() == operations_research::sat::OPTIMAL) {
                std::cout << "initSolution" << js.upperBound << std::endl;
                std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
                std::cout << "bks : " << js.bestKnownSolution << std::endl;
            }
            else {
                std::cout << "optimal solution not found!" << std::endl;
            }**/


            //concat the interval and after jobVars in order to create the permutation of the schedule in the next step.
            std::vector<std::vector<operations_research::sat::IntervalVar> > jobInterval(js.jobCount);


            std::vector<int> permutation = getPermutation(response, allJobVars, {});
            std::vector<Schedule> newSchedule = {Schedule(permutation)};
            evaluate(newSchedule, js);

            if (population[p].fitness >= newSchedule[0].fitness) {
                population[p] = newSchedule[0];
            }
            //newPop.push_back(Schedule(permutation));
        }
        /**evaluate(newPop, js);
        for(int k = 0; k < newPop.size(); k++) {
            population.emplace_back(newPop[k]);
        }**/
        return population;
    };
}

std::function<void(std::vector<Schedule> &, float, float, JobShop &, Logger &)> OrtoolsAlgorithm::
getORMutationInterval() {
    return [](std::vector<Schedule> &population, float mutationRate, float mutationImpact, JobShop &js, Logger &logger) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> mutationGenerator(0, 1);
        std::vector<Schedule> newSchedules;
        for (int p = 0; p < population.size(); p++) {
            if(mutationGenerator(gen) > mutationRate) continue;
            collectData(population[p], js);
            int intervalSize = population[p].fitness * mutationImpact;
            std::uniform_int_distribution<> intervalGenerator(0, population[p].fitness - intervalSize);
            int intervalStart = intervalGenerator(gen);


            //remove start of schedule
            std::vector<Task> interval, afterInterval;
            std::vector<int> currentTaskIndex(js.jobCount, 0), jobTime(js.jobCount, 0), machineTime(js.machineCount, 0);

            const int taskCount = js.jobCount * js.machineCount;

            operations_research::sat::CpModelBuilder cp_model;

            std::vector<operations_research::sat::IntVar> lastJobs;
            operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, js.upperBound});


            std::vector<std::vector<operations_research::sat::IntervalVar> > intervalJobVars(js.jobCount);
            std::vector<std::vector<operations_research::sat::IntervalVar> > intervalMachineVars(js.jobCount);

            std::vector<std::vector<operations_research::sat::IntervalVar> > afterJobVars(js.jobCount);
            std::vector<std::vector<operations_research::sat::IntervalVar> > afterMachineVars(js.jobCount);

            //store starting times of tasks before the interval in order to create the new schedule at the end.
            std::vector<std::pair<int, int> > beforeTasks;

            for (int i = 0; i < taskCount; i++) {
                const int jobIndex = population[p].permutation[i];
                const int nextTask = currentTaskIndex[jobIndex];
                Task task = js.jobs[jobIndex].tasks[nextTask];
                const int machine = task.machine;

                if (task.endTime < intervalStart) {
                    machineTime[machine] = task.endTime;
                    jobTime[jobIndex] = task.endTime;
                    beforeTasks.push_back(std::make_pair(task.startTime, jobIndex));
                } else if (task.startTime < intervalStart + intervalSize) {
                    const int duration = task.duration;
                    int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
                    operations_research::sat::IntVar start = cp_model.NewIntVar({taskLowerBound, js.upperBound});
                    operations_research::sat::IntVar end = cp_model.NewIntVar({taskLowerBound, js.upperBound});
                    operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);
                    if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                    intervalJobVars[task.job].push_back(taskInterval);
                    intervalMachineVars[task.machine].push_back(taskInterval);
                } else {
                    afterInterval.push_back(task);
                    const int duration = task.duration;
                    int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
                    operations_research::sat::IntVar start = cp_model.NewIntVar({0, js.upperBound});
                    operations_research::sat::IntVar end = cp_model.NewIntVar({0, js.upperBound});
                    if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                    operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);
                    cp_model.AddGreaterThan(taskInterval.StartExpr(), taskLowerBound);
                    afterJobVars[task.job].push_back(taskInterval);
                    afterMachineVars[task.machine].push_back(taskInterval);
                }

                currentTaskIndex[jobIndex]++;
            }

            //keep order of interval tasks for every job and make it smaller than the first task of the same job after the interval.
            for (int i = 0; i < js.jobCount; i++) {
                if (!afterJobVars[i].empty() && !intervalJobVars[i].empty()) {
                    cp_model.AddGreaterOrEqual(afterJobVars[i][0].StartExpr(),
                                               intervalJobVars[i][intervalJobVars[i].size() - 1].EndExpr());
                }
                for (int j = 0; j < (int) intervalJobVars[i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(intervalJobVars[i][j + 1].StartExpr(), intervalJobVars[i][j].EndExpr());
                }
            }

            //allow no overlap of tasks on the same machine in the intervall. Every task has to be before the first task on this machine after the interval
            for (int i = 0; i < js.machineCount; i++) {
                cp_model.AddNoOverlap(intervalMachineVars[i]);
                if (afterMachineVars[i].empty() || intervalMachineVars[i].empty()) continue;
                for (int j = 0; j < intervalMachineVars[i].size(); j++) {
                    cp_model.AddGreaterOrEqual(afterMachineVars[i][0].StartExpr(), intervalMachineVars[i][j].EndExpr());
                }
            }

            //keep the order of every task of each job after the interval
            for (int i = 0; i < js.jobCount; i++) {
                for (int j = 0; j < (int) afterJobVars[i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(afterJobVars[i][j + 1].StartExpr(), afterJobVars[i][j].EndExpr());
                }
            }
            //keep the order of every task on every machine after the interval
            for (int i = 0; i < js.jobCount; i++) {
                for (int j = 0; j < (int) afterMachineVars[i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(afterMachineVars[i][j + 1].StartExpr(), afterMachineVars[i][j].EndExpr());
                }
            }

            cp_model.AddMaxEquality(makespan, lastJobs);
            cp_model.Minimize(makespan);
            //build the model
            operations_research::sat::CpModelProto fixedModel = cp_model.Build();

            operations_research::sat::Model model;
            operations_research::sat::SatParameters parameters;
            parameters.set_max_time_in_seconds(25.0);
            model.Add(NewSatParameters(parameters));


            //solve the builded model.
            operations_research::sat::CpSolverResponse response = operations_research::sat::SolveCpModel(
                fixedModel, &model);

            if (response.status() != operations_research::sat::CpSolverStatus::OPTIMAL && response.status() != operations_research::sat::CpSolverStatus::FEASIBLE) {
                continue;
            }
            int curSol = operations_research::sat::SolutionIntegerValue(response, makespan);
            logger.log(curSol);
            //print solution.
            /*if(response.status() == operations_research::sat::OPTIMAL) {
                std::cout << "initSolution" << js.upperBound << std::endl;
                std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
                std::cout << "bks : " << js.bestKnownSolution << std::endl;
            }
            else {
                std::cout << "optimal solution not found!" << std::endl;
            }*/


            //concat the interval and after jobVars in order to create the permutation of the schedule in the next step.
            std::vector<std::vector<operations_research::sat::IntervalVar> > jobInterval(js.jobCount);

            for (int i = 0; i < js.jobCount; i++) {
                intervalJobVars[i].insert(intervalJobVars[i].end(), afterJobVars[i].begin(), afterJobVars[i].end());
            }

            std::vector<int> permutation = getPermutation(response, intervalJobVars, beforeTasks);
            std::vector<Schedule> newSchedule = {Schedule(permutation)};
            evaluate(newSchedule, js);

            if (population[p].fitness >= newSchedule[0].fitness) {
                population[p] = newSchedule[0];
            }
        }
        return population;
    };
}

std::function<void(std::vector<Schedule> &, float, float, JobShop &, Logger &)> OrtoolsAlgorithm::
getORMutationBottleneck() {
    std::vector<Schedule> newSchedules;
    return [](std::vector<Schedule> &population, float mutationRate, float mutationImpact, JobShop &js, Logger &logger) {
        for (int p = 0; p < population.size(); p++) {
            int bottleneck = collectDataAndBottleneck(population[p], js, mutationImpact);
            int intervalSize = population[p].fitness * mutationImpact;
            int intervalStart = bottleneck - intervalSize / 2;


            //remove start of schedule
            std::vector<int> currentTaskIndex(js.jobCount, 0), jobTime(js.jobCount, 0), machineTime(js.machineCount, 0);

            const int taskCount = js.jobCount * js.machineCount;

            operations_research::sat::CpModelBuilder cp_model;

            std::vector<operations_research::sat::IntVar> lastJobs;
            operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, js.upperBound});


            std::vector<std::vector<operations_research::sat::IntervalVar> > intervalJobVars(js.jobCount);
            std::vector<std::vector<operations_research::sat::IntervalVar> > intervalMachineVars(js.jobCount);

            std::vector<std::vector<operations_research::sat::IntervalVar> > afterJobVars(js.jobCount);
            std::vector<std::vector<operations_research::sat::IntervalVar> > afterMachineVars(js.jobCount);

            //store starting times of tasks before the interval in order to create the new schedule at the end.
            std::vector<std::pair<int, int> > beforeTasks;

            for (int i = 0; i < taskCount; i++) {
                const int jobIndex = population[p].permutation[i];
                const int nextTask = currentTaskIndex[jobIndex];
                Task task = js.jobs[jobIndex].tasks[nextTask];
                const int machine = task.machine;

                if (task.endTime < intervalStart) {
                    machineTime[machine] = task.endTime;
                    jobTime[jobIndex] = task.endTime;
                    beforeTasks.push_back(std::make_pair(task.startTime, jobIndex));
                } else if (task.startTime < intervalStart + intervalSize) {
                    const int duration = task.duration;
                    int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
                    operations_research::sat::IntVar start = cp_model.NewIntVar({taskLowerBound, js.upperBound});
                    operations_research::sat::IntVar end = cp_model.NewIntVar({taskLowerBound, js.upperBound});
                    operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);
                    if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                    intervalJobVars[task.job].push_back(taskInterval);
                    intervalMachineVars[task.machine].push_back(taskInterval);
                } else {
                    const int duration = task.duration;
                    int taskLowerBound = std::max(jobTime[task.job], machineTime[task.machine]);
                    operations_research::sat::IntVar start = cp_model.NewIntVar({0, js.upperBound});
                    operations_research::sat::IntVar end = cp_model.NewIntVar({0, js.upperBound});
                    if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                    operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);
                    cp_model.AddGreaterThan(taskInterval.StartExpr(), taskLowerBound);
                    afterJobVars[task.job].push_back(taskInterval);
                    afterMachineVars[task.machine].push_back(taskInterval);
                }

                currentTaskIndex[jobIndex]++;
            }

            //keep order of interval tasks for every job and make it smaller than the first task of the same job after the interval.
            for (int i = 0; i < js.jobCount; i++) {
                if (!afterJobVars[i].empty() && !intervalJobVars[i].empty()) {
                    cp_model.AddGreaterOrEqual(afterJobVars[i][0].StartExpr(),
                                               intervalJobVars[i][intervalJobVars[i].size() - 1].EndExpr());
                }
                for (int j = 0; j < (int) intervalJobVars[i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(intervalJobVars[i][j + 1].StartExpr(), intervalJobVars[i][j].EndExpr());
                }
            }

            //allow no overlap of tasks on the same machine in the intervall. Every task has to be before the first task on this machine after the interval
            for (int i = 0; i < js.machineCount; i++) {
                cp_model.AddNoOverlap(intervalMachineVars[i]);
                if (afterMachineVars[i].empty() || intervalMachineVars[i].empty()) continue;
                for (int j = 0; j < intervalMachineVars[i].size(); j++) {
                    cp_model.AddGreaterOrEqual(afterMachineVars[i][0].StartExpr(), intervalMachineVars[i][j].EndExpr());
                }
            }

            //keep the order of every task of each job after the interval
            for (int i = 0; i < js.jobCount; i++) {
                for (int j = 0; j < (int) afterJobVars[i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(afterJobVars[i][j + 1].StartExpr(), afterJobVars[i][j].EndExpr());
                }
            }
            //keep the order of every task on every machine after the interval
            for (int i = 0; i < js.jobCount; i++) {
                for (int j = 0; j < (int) afterMachineVars[i].size() - 1; j++) {
                    cp_model.AddGreaterOrEqual(afterMachineVars[i][j + 1].StartExpr(), afterMachineVars[i][j].EndExpr());
                }
            }

            cp_model.AddMaxEquality(makespan, lastJobs);
            cp_model.Minimize(makespan);
            //build the model
            operations_research::sat::CpModelProto fixedModel = cp_model.Build();

            operations_research::sat::Model model;
            operations_research::sat::SatParameters parameters;
            parameters.set_max_time_in_seconds(25.0);
            model.Add(NewSatParameters(parameters));


            //solve the builded model.
            operations_research::sat::CpSolverResponse response = operations_research::sat::SolveCpModel(
                fixedModel, &model);

            if (response.status() != operations_research::sat::CpSolverStatus::OPTIMAL && response.status() != operations_research::sat::CpSolverStatus::FEASIBLE) {
                continue;
            }
            int curSol = operations_research::sat::SolutionIntegerValue(response, makespan);
            logger.log(curSol);
            //print solution.
            /*if(response.status() == operations_research::sat::OPTIMAL) {
                std::cout << "initSolution" << js.upperBound << std::endl;
                std::cout << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
                std::cout << "bks : " << js.bestKnownSolution << std::endl;
            }
            else {
                std::cout << "optimal solution not found!" << std::endl;
            }*/


            //concat the interval and after jobVars in order to create the permutation of the schedule in the next step.
            std::vector<std::vector<operations_research::sat::IntervalVar> > jobInterval(js.jobCount);

            for (int i = 0; i < js.jobCount; i++) {
                intervalJobVars[i].insert(intervalJobVars[i].end(), afterJobVars[i].begin(), afterJobVars[i].end());
            }

            std::vector<int> permutation = getPermutation(response, intervalJobVars, beforeTasks);
            std::vector<Schedule> newSchedule = {Schedule(permutation)};
            evaluate(newSchedule, js);

            if (population[p].fitness >= newSchedule[0].fitness) {
                population[p] = newSchedule[0];
            }
        }
        return population;
    };


}

std::function<void(std::vector<Schedule> &, float, JobShop &, Logger &)> OrtoolsAlgorithm::getORCrossover() {
    return [](std::vector<Schedule> &population, float crossoverRate, JobShop &js, Logger &logger) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> mutationGenerator(0, 1);
        std::shuffle(population.begin(), population.end(), gen);
        std::vector<Schedule> newPop;

        for (int p = 0; p < population.size() - 1; p += 2) {
            if (mutationGenerator(gen) > crossoverRate) continue;
            collectData(population[p], js);
            std::vector<int> parent1 = population[p].permutation;
            std::vector<int> parent2 = population[p + 1].permutation;
            std::vector<std::pair<int, bool> > lcs = findLongestCommonSequence(parent1, parent2);
            //std::vector<std::pair<int, bool>> lcs = findLongestCommonMachineSequences(parent1, parent2, js);


            const int taskCount = js.jobCount * js.machineCount;
            int64_t horizon = 0;
            for(auto job : js.jobs) {
                for(auto task : job.tasks) {
                    horizon += task.duration;
                }
            }

            operations_research::sat::CpModelBuilder cp_model;

            std::vector<operations_research::sat::IntVar> lastJobs;
            operations_research::sat::IntVar makespan = cp_model.NewIntVar({0, horizon});


            std::vector<std::vector<operations_research::sat::IntervalVar> > allJobVars(js.jobCount);
            std::vector<std::vector<operations_research::sat::IntervalVar> > allMachineVars(js.jobCount);

            //std::vector<std::vector<operations_research::sat::IntervalVar>> notSelectedJobVars(js.jobCount);
            std::vector<std::vector<operations_research::sat::IntervalVar> > notSelectedMachineVars(js.jobCount);

            std::vector<int> currentTaskIndex(js.jobCount, 0);


            for (int i = 0; i < taskCount; i++) {
                const int jobIndex = parent1[i];
                const int nextTask = currentTaskIndex[jobIndex];
                Task task = js.jobs[jobIndex].tasks[nextTask];
                const int machine = task.machine;

                const int duration = task.duration;
                operations_research::sat::IntVar start = cp_model.NewIntVar({0, horizon});
                operations_research::sat::IntVar end = cp_model.NewIntVar({0, horizon});
                operations_research::sat::IntervalVar taskInterval = cp_model.NewIntervalVar(start, duration, end);


                if (lcs[i].second == false) {
                    allJobVars[jobIndex].push_back(taskInterval);
                    allMachineVars[machine].push_back(taskInterval);
                    if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                } else {
                    allJobVars[jobIndex].push_back(taskInterval);
                    notSelectedMachineVars[machine].push_back(taskInterval);
                    allMachineVars[machine].push_back(taskInterval);
                    if (currentTaskIndex[jobIndex] == js.machineCount - 1) {
                        lastJobs.push_back(end);
                    }
                }
                currentTaskIndex[jobIndex]++;
            }

            //keep order of interval tasks for every job and make it smaller than the first task of the same job after the interval.
            for (int i = 0; i < js.jobCount; i++) {
                for (int j = 0; j < js.machineCount - 1; j++) {
                    cp_model.AddGreaterOrEqual(allJobVars[i][j + 1].StartExpr(), allJobVars[i][j].EndExpr());
                }
            }

            //allow no overlap of tasks on the same machine in the intervall. Every task has to be before the first task on this machine after the interval
            for (int i = 0; i < js.machineCount; i++) {
                cp_model.AddNoOverlap(allMachineVars[i]);
                if (notSelectedMachineVars[i].size() > 1) {
                    for (int j = 0; j < notSelectedMachineVars[i].size() - 1; j++) {
                        cp_model.AddGreaterOrEqual(notSelectedMachineVars[i][j + 1].StartExpr(),
                                                   notSelectedMachineVars[i][j].EndExpr());
                    }
                }
            }

            cp_model.AddMaxEquality(makespan, lastJobs);
            cp_model.Minimize(makespan);

            operations_research::sat::Model model;

            operations_research::sat::SatParameters parameters;
            parameters.set_max_time_in_seconds(25.0);
            model.Add(NewSatParameters(parameters));


            //build the model
            operations_research::sat::CpModelProto fixedModel = cp_model.Build();

            //solve the builded model.


            operations_research::sat::CpSolverResponse response = operations_research::sat::SolveCpModel(
                fixedModel, &model);



            if (response.status() != operations_research::sat::CpSolverStatus::OPTIMAL && response.status() != operations_research::sat::CpSolverStatus::FEASIBLE) {
                continue;
            }
            int curSol = operations_research::sat::SolutionIntegerValue(response, makespan);
            logger.log(curSol);

            //print solution
            /**if(response.status() == operations_research::sat::OPTIMAL) {
                 << "initSolution" << js.upperBound << std::endl;
                 << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
                 << "bks : " << js.bestKnownSolution << std::endl;
            }
            else {
                 << "optimal solution not found!" << std::endl;
                 << "initSolution" << js.upperBound << std::endl;
                 << "solution: " << operations_research::sat::SolutionIntegerValue(response, makespan) << std::endl;
                 << "bks : " << js.bestKnownSolution << std::endl;
            }**/


            int worseParent;
            if (population[p].fitness < population[p + 1].fitness) {
                worseParent = p + 1;
            } else {
                worseParent = p;
            }

            if (operations_research::sat::SolutionIntegerValue(response, makespan) <= population[worseParent].fitness) {
                std::vector<int> permutation = getPermutation(response, allJobVars, {});
                std::vector<Schedule> newSchedule = {Schedule(permutation)};
                evaluate(newSchedule, js);
                population[worseParent] = newSchedule[0];
            }
            //newPop.push_back(Schedule(permutation));
        }

        /**evaluate(newPop, js);
        for(int k = 0; k < newPop.size(); k++) {
            population.emplace_back(newPop[k]);
        }**/

        return population;
    };
}


std::vector<int> OrtoolsAlgorithm::getPermutation(operations_research::sat::CpSolverResponse &response,
                                                  std::vector<std::vector<operations_research::sat::IntervalVar> > &
                                                  jobIntervals, std::vector<std::pair<int, int> > tasks) {
    std::vector<std::pair<int, int> > jobTaskPairs;
    std::vector<int> schedule;


    for (int i = 0; i < jobIntervals.size(); i++) {
        for (int j = 0; j < jobIntervals[i].size(); j++) {
            tasks.push_back(std::make_pair(SolutionIntegerValue(response, jobIntervals[i][j].StartExpr()), i));
        }
    }

    std::sort(tasks.begin(), tasks.end());

    for (const auto &pair: tasks) {
        schedule.push_back(pair.second);
    }

    return schedule;
}

int OrtoolsAlgorithm::collectDataAndBottleneck(Schedule &schedule, JobShop &js, const float impact) {
    std::vector<std::pair<int, int> > idleTimes;
    std::vector<int> jobTaskIndex(js.jobCount, 0), jobTime(js.jobCount, 0);
    std::vector<int> machineTime(js.machineCount, 0);
    for (int jobId: schedule.permutation) {
        Task &task = js.jobs[jobId].tasks[jobTaskIndex[jobId]];

        if (jobTaskIndex[jobId] > js.jobs[jobId].tasks.size()) {
            std::cout << "WARNING: Wrong number of tasks for job " << jobId << "!" << std::endl;
            std::cout << "#tasks : " << js.jobs[jobId].tasks.size() << "; jobTaskIndex : " << jobTaskIndex[jobId] <<
                    std::endl;
            continue;
        }

        task.startTime = std::max(machineTime[task.machine], jobTime[jobId]);
        task.endTime = task.startTime + task.duration;

        if (task.startTime - machineTime[task.machine] > 0) {
            idleTimes.push_back(std::make_pair(machineTime[task.machine], task.startTime - machineTime[task.machine]));
        }

        jobTime[jobId] = task.endTime;
        machineTime[task.machine] = task.endTime;
        jobTaskIndex[jobId]++;
    }
    schedule.fitness = *std::max_element(machineTime.begin(), machineTime.end());
    js.upperBound = schedule.fitness;

    std::sort(idleTimes.begin(), idleTimes.end());

    int AccessPoint = 0;
    int impactSize = schedule.fitness * impact;
    int radius = impactSize / 4;
    int stepSize = radius / 20;
    int maxIdleTime = 0;
    int bottleneck = schedule.fitness / 2;

    for (int i = radius; i < schedule.fitness - radius; i += stepSize) {
        int currentIdleTime = 0;
        int j = AccessPoint;
        while (idleTimes[j].first < i - radius && j < idleTimes.size()) {
            j++;
        }
        AccessPoint = j;

        while (idleTimes[j].first < i + radius && j < idleTimes.size()) {
            currentIdleTime += idleTimes[j].second;
            j++;
        }

        if (currentIdleTime > maxIdleTime) {
            maxIdleTime = currentIdleTime;
            bottleneck = i;
        }
    }
    if (bottleneck < impactSize / 2) {
        bottleneck = impactSize / 2;
    }
    return bottleneck;
}

std::vector<std::pair<int, bool> > OrtoolsAlgorithm::findLongestCommonSequence(const std::vector<int> &parent1,
                                                                               const std::vector<int> &parent2) {
    int m = parent1.size();
    int n = parent2.size();

    std::vector<std::vector<int> > dp(m + 1, std::vector<int>(n + 1, 0));

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (parent1[i - 1] == parent2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
            } else {
                dp[i][j] = std::max(dp[i - 1][j], dp[i][j - 1]);
            }
        }
    }
    int counter = 0;
    std::vector<bool> includedLCS(m, false);
    int i = m, j = n;
    while (i > 0 && j > 0) {
        if (parent1[i - 1] == parent2[j - 1]) {
            includedLCS[i - 1] = true;
            counter++;
            i--;
            j--;
        } else if (dp[i - 1][j] > dp[i][j - 1]) {
            i--;
        } else {
            j--;
        }
    }

    std::vector<std::pair<int, bool> > result;
    for (int r = 0; r < m; r++) {
        result.push_back({parent1[r], includedLCS[r]});
    }
    return result;
}


std::vector<std::pair<int, bool> > OrtoolsAlgorithm::findLongestCommonMachineSequences(const std::vector<int> &parent1,
    const std::vector<int> &parent2, JobShop &js) {
    std::vector<std::vector<int> > parent1MachineSchedule = scheduleToMachineSchedule(parent1, js);
    std::vector<std::vector<int> > parent2MachineSchedule = scheduleToMachineSchedule(parent2, js);
    int counter = 0;

    std::vector<std::vector<std::pair<int, bool> > > machineLCS(js.machineCount);
    for (int k = 0; k < js.machineCount; k++) {
        int m = parent1MachineSchedule[k].size();
        int n = parent2MachineSchedule[k].size();
        std::vector<std::vector<int> > dp(m + 1, std::vector<int>(n + 1, 0));

        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                if (parent1MachineSchedule[k][i - 1] == parent2MachineSchedule[k][j - 1]) {
                    dp[i][j] = dp[i - 1][j - 1] + 1;
                } else {
                    dp[i][j] = std::max(dp[i - 1][j], dp[i][j - 1]);
                }
            }
        }
        std::vector<bool> includedLCS(m, false);
        int i = m, j = n;
        while (i > 0 && j > 0) {
            if (parent1MachineSchedule[k][i - 1] == parent2MachineSchedule[k][j - 1]) {
                includedLCS[i - 1] = true;
                counter++;
                i--;
                j--;
            } else if (dp[i - 1][j] > dp[i][j - 1]) {
                i--;
            } else {
                j--;
            }
        }

        for (int r = 0; r < m; r++) {
            machineLCS[k].push_back({parent1MachineSchedule[k][r], includedLCS[r]});
        }
    }
    return machineScheduleToSchedule(machineLCS, js);
}

std::vector<std::vector<int> > OrtoolsAlgorithm::scheduleToMachineSchedule(std::vector<int> schedule, JobShop &js) {
    std::vector<std::vector<int> > machineSchedule(js.machineCount);
    std::vector<int> currentTaskIndex(js.jobCount, 0);
    for (int i = 0; i < schedule.size(); i++) {
        int jobIndex = schedule[i];
        int taskIndex = currentTaskIndex[jobIndex];
        int machine = js.jobs[jobIndex].tasks[taskIndex].machine;
        machineSchedule[machine].emplace_back(jobIndex);
        currentTaskIndex[jobIndex]++;
    }
    return machineSchedule;
}

std::vector<std::pair<int, bool> > OrtoolsAlgorithm::machineScheduleToSchedule(
    std::vector<std::vector<std::pair<int, bool> > > machineSchedule, JobShop js) {
    std::vector<int> currentJobIndex(js.jobCount, 0), jobTime(js.jobCount, 0), currentMachineIndex(js.machineCount, 0),
            machineTime(js.machineCount, 0);
    std::vector<std::pair<int, bool> > permutation = {};

    int taskCount = js.machineCount * js.jobCount;

    int errorDetection = 0;
    while (taskCount > 0) {
        errorDetection++;
        for (int i = 0; i < js.machineCount; i++) {
            if (currentMachineIndex[i] >= js.jobCount) continue;
            std::pair<int, bool> job = machineSchedule[i][currentMachineIndex[i]];
            if (js.jobs[job.first].tasks[currentJobIndex[job.first]].machine == i) {
                const int duration = js.jobs[job.first].tasks[currentJobIndex[job.first]].duration;
                const int startTime = std::max(machineTime[i], jobTime[job.first]);
                jobTime[job.first] = startTime + duration;
                machineTime[i] = startTime + duration;
                permutation.push_back(job);
                currentMachineIndex[i]++;
                currentJobIndex[job.first]++;
                taskCount--;
                errorDetection = 0;
            }
        }
        if (errorDetection > taskCount) {
            std::cout << "Error detected in schedule (ortools)..." << std::endl;
            exit(1);
        }
    }
    return permutation;
}

void OrtoolsAlgorithm::collectData(Schedule &schedule, JobShop &js) {
    std::vector<int> machineTime(js.machineCount, 0); // Tracks the finishing time of the current job for each machine
    std::unordered_map<int, int> jobTaskIndex; // Tracks the next task to be scheduled for each job
    std::unordered_map<int, std::unordered_map<int, int> > endTimes;

    for (int jobId: schedule.permutation) {
        if (jobTaskIndex[jobId] > js.jobs[jobId].tasks.size()) {
            std::cout << "WARNING: Wrong number of tasks for job " << jobId << "!" << std::endl;
            continue;
        }

        Task &task = js.jobs[jobId].tasks[jobTaskIndex[jobId]];

        //calculates earliest startTime + Duration -> endTime of the task
        int currentTask = jobTaskIndex[jobId];
        int const calculatedEndTime = std::max(machineTime[task.machine] + task.duration,
                                               endTimes[jobId][jobTaskIndex[jobId] - 1] + task.duration);
        endTimes[jobId][currentTask] = calculatedEndTime;

        task.endTime = calculatedEndTime;
        task.startTime = calculatedEndTime - task.duration;

        //sets the MachineTime of the machine used (represents until when the machine is occupied) to the endTime of the schedules Task
        machineTime[task.machine] = endTimes[jobId][currentTask];

        //increases the task counter of the current job
        jobTaskIndex[jobId]++;
    }

    // Calculate the fitness as the maximum finishing time across all machines
    js.upperBound = *std::max_element(machineTime.begin(), machineTime.end());

    //TODO: only for testing rn
    schedule.fitness = js.upperBound;
}
