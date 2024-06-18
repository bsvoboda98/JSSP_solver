//
// Created by Bennet on 30.01.2024.
//
#include <random>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <chrono>
#include "Structs.h"
#include "Logger.h"
#include <functional>
#include <map>
#include <utility>
#include <numeric>

#include "OrtoolsAlgorithm.h"


#ifndef JSSP_SOLVER_FUNCTIONFACTORY_H
#define JSSP_SOLVER_FUNCTIONFACTORY_H


class FunctionFactory {
public:
    std::function<void(std::vector<Schedule>&, float, JobShop&, Logger&)> crossover;
    std::function<void(std::vector<Schedule>&, float, float, JobShop&, Logger&)> mutate;
    std::function<std::vector<Schedule>(std::vector<Schedule>&, int)> select;

    FunctionFactory()= default;

    explicit FunctionFactory(int i){
        this->crossover = getCrossover();
        this->mutate = getMutate();
        this->select = getSelect();
    }

    FunctionFactory(std::string crossover, std::string mutate, std::string select) {

            if(crossover == "spc"){
                //single point crossover function.
                this->crossover = getSinglePointCrossover();
            }
            else if (crossover == "tpc") {
                //two point crossover function.
                this->crossover = getTwoPointCrossover();
            }
            else if (crossover == "cc") {
                //cycle crossover function.
                this->crossover = getCycleCrossover();
            }
            else if (crossover == "or") {
                this->crossover = OrtoolsAlgorithm::getORCrossover();
            }


            if(mutate == "swm"){
                //swap mutation
                this->mutate = getSwapMutation();
            }
            else if(mutate == "im"){
                this->mutate = getInversionMutation();
            }
            else if(mutate == "rim"){
                this->mutate = getRangedInversionMutation();
            }
            else if(mutate == "scm"){
                this->mutate = getScrambleMutation();
            }
            else if(mutate == "or") {
                this->mutate = OrtoolsAlgorithm::getORMutation();
            }
            else if(mutate == "ori") {
                this->mutate = OrtoolsAlgorithm::getORMutationInterval();
            }
            else if(mutate == "orb") {
                this->mutate = OrtoolsAlgorithm::getORMutationBottleneck();
            }

            if(select == "dc"){
                //deterministic crowding
                this->select = getDeterministicCrowding();
            }
            else if(select == "rws"){
                this->select = getRouletteWheelSelection();
            }

            else if(select == "rs"){
                this->select = getRankSelection();
            }
            else if(select == "ts"){
                this->select = getTournamentSelection();
            }
            else if(select == "ra") {
                this -> select = getRandomSelection();
            }
            else if(select == "no") {
                this->select = getNoSelection();
            }

    }

    /**
     * Interaction with the user. Reading from input which crossover shall be used.
     *
     * @return the function that will be used in the GA as crossover as lambda function
     */
    static std::function<void(std::vector<Schedule>&, float, JobShop&, Logger&)> getCrossover(){
        std::cout
                << "WHICH CROSSOVER FUNCTION SHOULD BE USED?\n"
                << "AVAILABLE CROSSOVER VARIANTS:\n"
                << "single point crossover (spc)\n"
                << "two point crossover (tpc)\n"
                << "cycle crossover (cc) ~ standard" << std::endl;


        while(true){
            std::string userInput;
            std::getline(std::cin, userInput);

            if(userInput == "spc"){
                //single point crossover function.
                return getSinglePointCrossover();
            }

            if (userInput == "tpc") {
                //two point crossover function.
                return getTwoPointCrossover();
            }

            if (userInput == "cc" || userInput == "") {
                //cycle crossover function.
                return getCycleCrossover();
            }
            std::cout << "INVALID INPUT" << std::endl;
        }
    }


    /**
     * Interaction with the user. Reading from input which mutate function shall be used.
     *
     * @return the function that will be used in the GA as mutate as lambda function
     */
    static std::function<void(std::vector<Schedule>&, float, float, JobShop&, Logger&)> getMutate(){
        std::cout
                << "WHICH MUTATE FUNCTION SHOULD BE USED?\n"
                << "AVAILABLE MUTATION VARIANTS:\n"
                << "swap mutation (swm)\n"
                << "inversion mutation (im)\n"
                << "ranged inversion mutation (rim)\n"
                << "scramble mutation (scm) ~ standard" << std::endl;

        while(true){
            std::string userInput;
            std::getline(std::cin, userInput);

            if(userInput == "swm"){
                //swap mutation
                return getSwapMutation();
            }
            if(userInput == "im"){
                return getInversionMutation();
            }
            if(userInput == "rim"){
                return getRangedInversionMutation();
            }
            if(userInput == "scm" || userInput == ""){
                return getScrambleMutation();
            }
            std::cout << "INVALID INPUT" << std::endl;
        }
    }





    /**
     * Interaction with the user. Reading from input which select function shall be used.
     *
     * @return the function that will be used in the GA as mutate as lambda function
     */
     static std::function<std::vector<Schedule>(std::vector<Schedule>&, int)> getSelect(){

        std::cout
                << "WHICH SELECT FUNCTION SHOULD BE USED?\n"
                << "AVAILABLE MUTATION VARIANTS:\n"
                << "deterministic crowding (dc)\n"
                << "roulette wheel selection (rws) ~ standard\n"
                << "rank selection (rs)\n"
                << "tournament selection (ts)" << std::endl;

        while(true){
            std::string userInput;
            std::getline(std::cin, userInput);

            if(userInput == "dc"){
                //deterministic crowding
                return getDeterministicCrowding();
            }
            else if(userInput == "rws" || userInput == ""){
                return getRouletteWheelSelection();
            }

            else if(userInput == "rs"){
                return getRankSelection();
            }
            else if(userInput == "ts"){
                return getTournamentSelection();
            }
            else{
                std::cout << "INVALID INPUT" << std::endl;
            }
        }
     }





     static std::pair<std::vector<int>, std::vector<int>> cross(std::vector<int> offspring1, std::vector<int> offspring2, int jobCount){
         int offspringLength = offspring1.size();

         std::map<int, int> offspringCount1;
         std::map<int, int> offspringCount2;
         std::map<int, int> offspringShared;
         std::map<int, int> offspringUnique1;
         std::map<int, int> offspringUnique2;

         for (int j = 0; j <= jobCount; j++) {
             offspringCount1[j] = 0;
             offspringCount2[j] = 0;
         }

         for (int j = 0; j < offspringLength; j++) {
             offspringCount1[offspring1[j]]++;
             offspringCount2[offspring2[j]]++;
         }

         for (int j = 0; j <= jobCount; j++) {
             offspringShared[j] = std::min(offspringCount1[j], offspringCount2[j]);
             offspringUnique1[j] = offspringCount1[j] - offspringShared[j];
             offspringUnique2[j] = offspringCount2[j] - offspringShared[j];
         }

         std::vector<int> crossover1;
         for (int j = 0; j < offspring1.size(); j++) {
             if (offspringUnique1[offspring1[j]] > 0) {
                 crossover1.push_back(offspring1[j]);
                 offspringUnique1[offspring1[j]]--;
                 offspring1.erase(offspring1.begin() + j);
                 j--;
             }
         }

         std::vector<int> crossover2;
         for (int j = 0; j < offspring2.size(); j++) {
             if (offspringUnique2[offspring2[j]] > 0) {
                 crossover2.push_back(offspring2[j]);
                 offspringUnique2[offspring2[j]]--;
                 offspring2.erase(offspring2.begin() + j);
                 j--;
             }
         }

         offspringLength = offspring1.size();
         for (int j = 0; j < offspringLength; j++) {
             crossover1.push_back(offspring2[j]);
             crossover2.push_back(offspring1[j]);
         }

         return std::make_pair(crossover1, crossover2);

     }





    /**
     * std::vector<Schedule>& population    : population on which to perform the crossover operation
     *
     * performing single point crossover on the given population.
     * Meaning that 2 schedules are split at one point and concatenated in pairs.
     */
    static std::function<void(std::vector<Schedule>&, float, JobShop&, Logger&)> getSinglePointCrossover(){
         auto lambda = [](std::vector<Schedule>& population, float crossoverRate, JobShop &js, Logger &logger){
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    int size = (int) population[0].permutation.size();
                    std::uniform_int_distribution<> dis(0, size - 1);
                    std::uniform_real_distribution<> floatDis(0.0f, 1.0f);
                    int popCount = population.size();
                    int jobCount = *std::max_element(population[0].permutation.begin(), population[0].permutation.end());

                    for (size_t i = 0; i < popCount; i += 2) {
                        if(floatDis(gen) > crossoverRate){
                            continue;
                        }

                        int crossoverPosition = dis(gen);

                        std::vector<int> parent1 = population[i].permutation;
                        std::vector<int> parent2 = population[i + 1].permutation;


                        //copy the parents
                        std::vector<int> child1(parent1.begin(), parent1.begin() + crossoverPosition);
                        std::vector<int> child2(parent2.begin(), parent2.begin() + crossoverPosition);


                        std::vector<int> offspring1(parent1.begin() + crossoverPosition, parent1.end());
                        std::vector<int> offspring2(parent2.begin() + crossoverPosition, parent2.end());

                        auto [crossover1, crossover2] = cross(offspring1, offspring2, jobCount);

                        child1.insert(child1.end(), crossover1.begin(), crossover1.end());
                        child2.insert(child2.end(), crossover2.begin(), crossover2.end());


                        population.push_back(Schedule(child1));
                        population.push_back(Schedule(child2));


                    }
                };
                return lambda;
     }





    /**
     * std::vector<Schedule>& population    : population on which to perform the crossover operation
     *
     * performing two point crossover on the given population.
     * Meaning subsequences from 2 schedules are exchanged. Those subsequences are in the same position in their schedules.
     */
    static std::function<void(std::vector<Schedule>&, float, JobShop&, Logger&)> getTwoPointCrossover(){
        auto lambda = [](std::vector<Schedule> &population, float crossoverRate, JobShop &js, Logger &logger) {
            int crossoverImpact = 20;
            std::random_device rd;
            std::mt19937 gen(rd());
            int size = (int) population[0].permutation.size();
            std::uniform_int_distribution<> dis(0, size - 1 - crossoverImpact);
            std::uniform_real_distribution<> floatDis(0.0f, 1.0f);
            int popCount = population.size();
            int jobCount = *std::max_element(population[0].permutation.begin(), population[0].permutation.end());

            for (size_t i = 0; i < popCount; i += 2) {
                if(floatDis(gen) > crossoverRate){
                    continue;
                }

                std::vector<int> parent1 = population[i].permutation;
                std::vector<int> parent2 = population[i + 1].permutation;

                //where to cut the schedules
                int crossoverPosition1 = dis(gen);
                int crossoverPosition2 = dis(gen);

                std::vector<int> child1(parent1.begin(), parent1.begin() + crossoverPosition1);
                std::vector<int> child2(parent2.begin(), parent2.begin() + crossoverPosition2);

                std::vector<int> offspring1(parent1.begin() + crossoverPosition1, parent1.begin() + crossoverPosition1 + crossoverImpact);
                std::vector<int> offspring2(parent2.begin() + crossoverPosition2, parent2.begin() + crossoverPosition2 + crossoverImpact);

                auto [crossover1, crossover2] = cross(offspring1, offspring2, jobCount);


                child1.insert(child1.end(), crossover1.begin(), crossover1.end());
                child2.insert(child2.end(), crossover2.begin(), crossover2.end());

                child1.insert(child1.end(), parent1.begin() + crossoverPosition1 + crossoverImpact, parent1.end());
                child2.insert(child2.end(), parent2.begin() + crossoverPosition2 + crossoverImpact, parent2.end());

                population.push_back(Schedule(child1));
                population.push_back(Schedule(child2));
            }
        };
        return lambda;
    }





    /**
     * std::vector<Schedule>& population    : population on which to perform the crossover operation
     *
     * performing two point crossover on the given population.
     * Meaning subsequences from 2 schedules are exchanged. Those subsequences are in the same position in their schedules.
     */
    static std::function<void(std::vector<Schedule>&, float, JobShop&, Logger&)> getCycleCrossover() {
        auto lambda = [](std::vector<Schedule> &population, float crossoverRate, JobShop &js, Logger &logger) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> floatDis(0.0f, 1.0f);

            int popCount = population.size();
            int size = population[0].permutation.size();
            for (int i = 0; i < popCount; i += 2) {
                if(floatDis(gen) > crossoverRate){
                    continue;
                }

                std::vector<int> parent1 = population[i].permutation;
                std::vector<int> parent2 = population[i + 1].permutation;

                std::vector<bool> visited(parent1.size(), false);

                int currentCycle = 0;
                std::map<int, int> cycles;

                for (int j = 0; j < size; j++) {
                    cycles[j] = 0;
                }


                for (int j = 0; j < size; j++) {
                    if (!visited[j]) {
                        currentCycle++;
                        int cycleStart = j;
                        int cycleEnd = j;
                        while (true) {
                            for (int k = 0; k < visited.size(); k++) {
                                if (!visited[k] && parent1[cycleEnd] == parent2[k]) {
                                    cycleEnd = k;
                                    visited[k] = true;
                                    cycles[k] = currentCycle;
                                    break;
                                }
                            }
                            if (cycleEnd == cycleStart) { break; }
                        }
                    }
                }

                int cycleCount = 0;
                for (const auto entry: cycles) {
                    if (entry.second > cycleCount) {
                        cycleCount = entry.second;
                    }
                }

                std::vector<int> child1, child2;

                std::uniform_int_distribution<> cDistrib(1, cycleCount);

                for (int k = 0; k < size; k++) {
                    if (cycles[k] % 3 == 0) {
                        child1.push_back(parent1[k]);
                        child2.push_back(parent2[k]);
                    } else {
                        child1.push_back(parent2[k]);
                        child2.push_back(parent1[k]);
                    }
                }
                population.push_back(Schedule(child1));
                population.push_back(Schedule(child2));
                child1.clear();
                child2.clear();
            }
            return population;
        };
        return lambda;
    }





    /**
     * std::vector<Schedule>& population    : population on which to perform the mutate operation
     *
     * this function applies a random mutation to a part of the population (controlled by the mutationRate)
     * since we have to ensure unique floats and we dont want to have while(!unique) loops, the operation fails if the new value is not unique. TODO: find a better way to do keep the random key encoding
     * if the new Value (newVal) is smaller than 0/greater than 1 we will decrease/increase the value by (val * mutationImpact)/((1 - val) * mutationImpact)
     */
    static std::function<void(std::vector<Schedule> &, float, float, JobShop&, Logger&)> getSwapMutation() {

        return [](std::vector<Schedule>& population, float mutationRate, float mutationImpact, JobShop &js, Logger &logger) -> void {
            int impactRange = 20;
            float fImpact = 10.0f * mutationImpact - 5.0f;
            int impact = static_cast<int>(fImpact);

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);
            int size = population[0].permutation.size();
            std::uniform_int_distribution<> index_dis(0, size - impactRange - impact);

            int popCount = population.size();

            for (int i = 0; i < popCount; i++){

                Schedule newSchedule = population[i];
                int pos1 = index_dis(gen);
                int pos2 = index_dis(gen);

                if(pos1 > pos2){
                    int temp = pos2;
                    pos2 = pos1;
                    pos1 = temp;
                }

                if(pos2 - pos1 < impactRange + impact){
                    if(pos2 > size / 2){
                        pos1 -= impactRange + impact;
                    }
                    else{
                        pos2 += impactRange + impact;
                    }
                }

                std::swap_ranges(newSchedule.permutation.begin() + pos1, newSchedule.permutation.begin() + pos1 + impactRange + impact, newSchedule.permutation.begin() + pos2);
                population.push_back(newSchedule);
            }
        };
    }





    static std::function<void(std::vector<Schedule>&, float, float, JobShop&, Logger&)> getInversionMutation() {
        return [](std::vector<Schedule>& population, float mutationRate, float mutationImpact, JobShop& js, Logger &logger) -> void {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);

            int popCount = population.size();
            int jobCount = *std::max_element(population[0].permutation.begin(), population[0].permutation.end());
            int size = population[0].permutation.size();


            for(int i = 0; i < popCount; i++){
                float threshold = dis(gen);
                if(threshold < mutationRate){
                    Schedule newSchedule = population[i];
                    for(int j = 0; j < size; j++){
                        newSchedule.permutation[j] = jobCount - newSchedule.permutation[j];
                    }
                    population.push_back(newSchedule);
                }
            }
        };

    }




    static std::function<void(std::vector<Schedule>&, float, float, JobShop&, Logger&)> getRangedInversionMutation() {
        return [](std::vector<Schedule>& population, float mutationRate, float mutationImpact, JobShop &js, Logger &logger) -> void {
            float fImpact = 10.0f * mutationImpact - 5.0f;
            int impact = static_cast<int>(fImpact);

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);
            int size = population[0].permutation.size();
            std::uniform_int_distribution<> index_dis(0, size - 10 - impact);

            int popCount = population.size();
            int jobCount = *std::max_element(population[0].permutation.begin(), population[0].permutation.end());


            for(int i = 0; i < popCount; i++){
                float threshold = dis(gen);
                if(threshold < mutationRate){
                    int pos = index_dis(gen);
                    Schedule newSchedule = population[i];
                    for(int j = pos; j < pos - 5; j++){
                        newSchedule.permutation[j] = jobCount - newSchedule.permutation[j];
                    }
                    population.push_back(newSchedule);
                }
            }
        };

    }





    static std::function<void(std::vector<Schedule>&, float, float, JobShop&,Logger&)> getScrambleMutation() {
        return [](std::vector<Schedule>& population, float mutationRate, float mutationImpact, JobShop &js,Logger &logger) -> void {
            int impactRange = 25;
            float fImpact = 10.0f * mutationImpact - 5.0f;
            int impact = static_cast<int>(fImpact);
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);
            int size = population[0].permutation.size();
            std::uniform_int_distribution<> index_dis(0, size - impactRange - impact);

            int popCount = population.size();

            for(int i = 0; i < popCount; i++){
                float threshold = dis(gen);
                if(threshold < mutationRate){
                    int pos = index_dis(gen);
                    Schedule newSchedule = population[i];
                    std::shuffle(newSchedule.permutation.begin() + pos, newSchedule.permutation.begin() + pos + impactRange + impact, gen);
                    population.push_back(newSchedule);
                }
            }
        };
    }






    /**
     * std::vector<Schedule>& population    : population where to select the "populationSize" best schedules from
     *
     * compare new schedules with their parents and keep the better one in the population
     */
    static std::function<std::vector<Schedule>(std::vector<Schedule>&, int)> getDeterministicCrowding(){
        return [](std::vector<Schedule>& population, int populationSize) -> std::vector<Schedule> {
            std::vector<Schedule> bestPop;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine rng(seed);

            for (int i = 0; i < populationSize; i++) {
                if (i + populationSize < population.size()) {
                    if (population[i].fitness < population[i + populationSize].fitness) {
                        bestPop.push_back(population[i]);
                    } else {
                        bestPop.push_back(population[i + populationSize]);
                    }
                } else {
                    bestPop.push_back(population[i]);
                }
            }

            std::shuffle(bestPop.begin(), bestPop.end(), rng);
            return bestPop;
        };
    }




    static std::function<std::vector<Schedule>(std::vector<Schedule>&, int)> getRouletteWheelSelection(){
        return [](std::vector<Schedule>& population, int populationSize) -> std::vector<Schedule> {
            std::vector<Schedule> selectedPop;

            std::sort(population.begin(), population.end());
            int worstFit = population.back().fitness;

            selectedPop.insert(selectedPop.end(), population.begin(), population.begin() + 5);
            population.erase(population.begin(), population.begin() + 5);

            int sum = 0;
            std::vector<int> wheel;
            int probability;

            for(auto & i : population){
                probability = worstFit + 500 - i.fitness;
                wheel.push_back(probability);
                sum += probability;
            }

            std::random_device rd;
            std::mt19937 gen(rd());

            for(int i = 0; i < populationSize - 5; i++){
                std::uniform_int_distribution<> index_dis(0, sum - 1);
                int spin = index_dis(gen);
                for(int j = 0; j < population.size(); j++){
                    if(spin < wheel[j]){
                        selectedPop.push_back(population[j]);
                        sum -= wheel[j];
                        wheel[j] = 0;
                        break;
                    }
                    else{
                        spin -= wheel[j];
                    }
                }
            }
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine rng(seed);
            std::shuffle(selectedPop.begin(), selectedPop.end(), rng);
            return selectedPop;
        };
    }





    static std::function<std::vector<Schedule>(std::vector<Schedule>&, int)> getRankSelection() {
        return [](std::vector<Schedule> &population, int populationSize) -> std::vector<Schedule> {
            std::vector<Schedule> selectedPop;

            std::sort(population.begin(), population.end());

            selectedPop.insert(selectedPop.end(), population.begin(), population.begin() + 5);
            population.erase(population.begin(), population.begin() + 5);

            std::vector<int> rankProbability;

            int size = population.size();
            int sum = 0;
            for (int i = 0; i < size; i++) {
                rankProbability.push_back(size + 1 - i);
                sum += rankProbability[i];
            }

            std::random_device rd;
            std::mt19937 gen(rd());
            for(int i = 0; i < populationSize; i++){
                std::uniform_int_distribution<> index_dis(0, sum - 1);
                int spin = index_dis(gen);
                for(int j = 0; j < size; j++) {
                    if (spin < rankProbability[j]) {
                        selectedPop.push_back(population[j]);
                        rankProbability[j] = 0;
                        sum -= rankProbability[j];
                        break;
                    } else {
                        spin -= rankProbability[j];
                    }
                }
            }

            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine rng(seed);
            std::shuffle(selectedPop.begin(), selectedPop.end(), rng);
            return selectedPop;

        };
    }



    static std::function<std::vector<Schedule>(std::vector<Schedule>&, int)> getTournamentSelection() {
        return [](std::vector<Schedule> &population, int populationSize) -> std::vector<Schedule> {
            std::vector<Schedule> selectedPop;
            std::vector<int> side1, side2;

            for(int i = 0; i < populationSize; i++){
                side1.push_back(i);
                side2.push_back(i + populationSize);
            }

            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine rng(seed);

            std::shuffle(side1.begin(), side1.end(), rng);
            std::shuffle(side2.begin(), side2.end(), rng);

            for(int i = 0; i < populationSize; i++){
                if(side2[i] < population.size()){
                    if(population[side1[i]].fitness < population[side2[i]].fitness){
                        selectedPop.push_back(population[side1[i]]);
                    }
                    else{
                        selectedPop.push_back(population[side2[i]]);
                    }
                }

                else{
                    selectedPop.push_back(population[side1[i]]);
                }
            }
            std::shuffle(selectedPop.begin(), selectedPop.end(), rng);
            return selectedPop;
        };
    }




    static std::function<std::vector<Schedule>(std::vector<Schedule>&, int)> getRandomSelection() {

        return [](std::vector<Schedule> &population, int populationSize) -> std::vector<Schedule> {
            std::vector<Schedule> selectedPop;

            std::sort(population.begin(), population.end());

            selectedPop.insert(selectedPop.end(), population.begin(), population.begin() + 5);
            population.erase(population.begin(), population.begin() + 5);

            std::vector<int> indices(population.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::random_shuffle(indices.begin(), indices.end());

            for(int i = 0; i < populationSize - 5; i++) {
                selectedPop.push_back(population[indices[i]]);
            }

            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(selectedPop.begin(), selectedPop.end(), g);


            return selectedPop;
        };
    }
    static std::function<std::vector<Schedule>(std::vector<Schedule>&, int)> getNoSelection() {

        return [](std::vector<Schedule> &population, int populationSize) -> std::vector<Schedule> {
            return population;
        };
    }

};


#endif //JSSP_SOLVER_FUNCTIONFACTORY_H
