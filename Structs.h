//
// Created by Bennet on 30.01.2024.
//
#include <memory>
#include <vector>
#include <string>
#include <ctime>
#include "include/json.hpp"


#ifndef JSSP_SOLVER_STRUCTS_H
#define JSSP_SOLVER_STRUCTS_H

struct Schedule {
    std::vector<int> permutation;
    int fitness = 0;

    bool operator<(const Schedule& other) const{
        return fitness < other.fitness;
    }

    bool operator==(const Schedule& other) const {
        return permutation == other.permutation;
    }

    explicit Schedule(std::vector<int> p, int f) : permutation((std::move(p))) {fitness = f;}
    explicit Schedule(std::vector<int> p ) : permutation(std::move(p)) {}
    Schedule()= default;
};

struct Task {
    int machine = 0;
    int duration = 0;
    int startTime = 0;
    int endTime = 0;
    int job = 0;

    Task(const int m, const int d, const int j) : machine(m), duration(d), job(j) {}
    Task()= default;
};

struct Job {
    std::vector<Task> tasks;

    explicit Job(std::vector<Task> v) : tasks(std::move(v)){}
};

struct JobShop{

    JobShop(std::vector<Job> j, int mc, int jc, int bks, int lb, std::string n) :
            jobs(std::move(j)), machineCount(mc), jobCount(jc), bestKnownSolution(bks), lowerBound(lb), name(n)  {}


    std::vector<Job> jobs;
    int machineCount;
    int jobCount;
    int bestKnownSolution;
    int lowerBound;
    int upperBound;
    std::string name;
};
#endif //JSSP_SOLVER_STRUCTS_H
