//
// Created by svobee on 11.05.24.
//
#include <ctime>
#include <iostream>
#include "include/json.hpp"
#include "Structs.h"
#include <filesystem>
#include <fstream>
#include <limits.h>
#include <limits>


#ifndef LOGGER_H
#define LOGGER_H

class Logger {
public:
    Logger(const int makespan, JobShop &js, std::vector<std::string> config) {
        initlogg(js, config);
        log(makespan);
    }
    Logger(JobShop &js, std::vector<std::string>config) {
        initlogg(js, config);
    }

    void log(const int bestFitness) {
        const std::chrono::duration<double> elapsed_seconds = (std::chrono::system_clock::now() - startTime);
        const double seconds = elapsed_seconds.count();
        if(bestFitness < curOpt) {
            curOpt = bestFitness;
            results["benchmark"][std::to_string(seconds)] = curOpt;
        }
    }

    void save() {
        std::filesystem::path dirPath;
        std::string usedFunctions = "";
        if(results["algorithm"].get<std::string>() == "ga") {
            usedFunctions = results["functions"]["mutate"].get<std::string>() + "_" + results["functions"]["crossover"].get<std::string>() + "_";
            dirPath = std::filesystem::path("../benchmark/" + results["functions"]["mutate"].get<std::string>() + "_" + results["functions"]["crossover"].get<std::string>() + "_" + std::to_string(results["geneticAlgorithm"]["mutateImpact"].get<double>()) + "/");

        }
        else if(results["algorithm"].get<std::string>() == "m") {
            usedFunctions = "m" + results["functions"]["mutate"].get<std::string>() + "_" + results["functions"]["crossover"].get<std::string>() + "_";
            dirPath = std::filesystem::path("../benchmark/" + usedFunctions + std::to_string(results["geneticAlgorithm"]["mutateImpact"].get<double>()) + "/");

        }
        else {
            dirPath = std::filesystem::path("../benchmark/" + results["algorithm"].get<std::string>() + "/");
        }

        if(!std::filesystem::exists(dirPath)) {
            std::filesystem::create_directories(dirPath);
        }

        const std::string fileName = results["benchmarkClass"].get<std::string>() + results["benchmarkNumber"].get<std::string>() + ".json";


        std::filesystem::path filePath = dirPath / fileName;

        std::ofstream file(filePath);
        file << results.dump(4);
        file.close();
    }

    void logSchedule(Schedule schedule) {
        results["schedule"] = schedule.permutation;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::high_resolution_clock::now();
    int curOpt = INT_MAX;
    nlohmann::json results;

    void initlogg(JobShop &js, std::vector<std::string> &config) {
        auto now = std::chrono::system_clock::now();
        std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
        std::tm* localTime = std::localtime(&currentTime);

        bool g = false;
        bool f = false;
        bool t = false;


        std::string algorithm = config[0].substr(1);
        int maxGen;
        int popSize;
        int timeLimit;
        double crossRate;
        double mutateRate;
        double mutateImpact;

        std::string crossover;
        std::string mutate;
        std::string select;

        for(int i = 1; i < config.size(); i++) {
            if(config[i] == "-g") {
                g = true;
                popSize = std::stoi(config[i+1]);
                maxGen = std::stoi(config[i+2]);
                crossRate = std::round(std::stod(config[i+3]) * 100) / 100;
                mutateRate = std::round(std::stod(config[i+4]) * 100) / 100;
                mutateImpact = std::round(std::stod(config[i+5]) * 100) / 100;
                i += 5;
            }
            else if(config[i] == "-f") {
                f = true;
                crossover = config[i+1];
                mutate = config[i+2];
                select = config[i+3];
                i += 3;
            }
            else if(config[i] == "-t") {
                t = true;
                timeLimit = std::stoi(config[i+1]);
            }
        }

        std::stringstream benchmarkClassStream;
        std::stringstream benchmarkNumberStream;

        for(char c : js.name) {
            if(std::isdigit(c)) {
                benchmarkNumberStream << c;
            }
            else {
                benchmarkClassStream << c;
            }
        }

        std::string benchmarkClass = benchmarkClassStream.str();
        std::string benchmarkNumber = benchmarkNumberStream.str();

        results =
        {
            {"benchmarkClass", benchmarkClass},
            {"benchmarkNumber", benchmarkNumber},
            {"algorithm", algorithm},
            {"bks", js.bestKnownSolution},
            {"machineCount", js.machineCount},
            {"jobCount", js.jobCount},
            {"benchmark", {}}
        };

        results["timestamp"]["year"] = localTime->tm_year + 1900;
        results["timestamp"]["month"] = localTime->tm_mon + 1;
        results["timestamp"]["day"] = localTime->tm_mday;
        results["timestamp"]["hour"] = localTime->tm_hour;
        results["timestamp"]["minute"] = localTime->tm_min;

        if(g) {
            results["geneticAlgorithm"]["popSize"] = popSize;
            if(!t) {
                results["geneticAlgorithm"]["maxGen"] = maxGen;
            }
            results["geneticAlgorithm"]["crossRate"] = crossRate;
            results["geneticAlgorithm"]["mutateRate"] = mutateRate;
            results["geneticAlgorithm"]["mutateImpact"] = mutateImpact;
        }

        if(f) {
            results["functions"]["crossover"] = crossover;
            results["functions"]["mutate"] = mutate;
            results["functions"]["select"] = select;
        }

        if(t) {
            results["timeLimit"] = timeLimit;
        }
    }




};
#endif //LOGGER_H
