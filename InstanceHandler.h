//
// Created by Bennet on 30.01.2024.
//
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <tuple>
#include "Structs.h"


#ifndef JSSP_SOLVER_INSTANCEHANDLER_H
#define JSSP_SOLVER_INSTANCEHANDLER_H


class InstanceHandler {
public:
    static std::vector<JobShop> getJobShopsByString(std::string userInput) {
        std::vector<std::string> instances = parseInput(userInput);
        return readFile(instances);
    }

    /**
     * Interaction with user. Presents all available instances in the console.
     * Afterwards the user types in the instances he wants to use.
     *
     * The userinput will be used for calling the function "parseInput" which is transformating the user input into a list of instance names.
     *
     * @return the received list is used to call the "read file" function which is returning all JobShops with a name from the list.
     */
    static std::vector<JobShop> getJobShops() {
        std::string userInput;
        std::cout
                << "ENTER INSTANCES :\n"
                << "TYPE \"all\" FOR ALL INSTANCES.\n\n"
                << "AVAILABLE INSTANCES :\n"
                << "abz5-9\n"
                << "dmu01-80\n"
                << "ft06, ft10, ft20\n"
                << "la01-40\n"
                << "orb01-10\n"
                << "swv01-20\n"
                << "ta01-80\n"
                << "yn1-4" << std::endl;


        std::getline(std::cin, userInput);

        std::vector<std::string> instances = parseInput(userInput);

        return readFile(instances);
    }

    /**
     * Every two numbers represent a task (int machine, int duration).
     * This function reads all numbers and form tasks out of them.
     *
     * @param line  : string which contains data for a job.
     * @return Job created with a vector of all tasks the function has created.
     */
    static Job getJob(std::string &line, int job) {
        std::stringstream findNumbers;
        int number;
        int machine;
        int duration;

        std::vector<Task> tasks = {};
        findNumbers.str("");
        findNumbers.clear();
        findNumbers << line;
        while (findNumbers >> number) {
            machine = number;
            findNumbers >> number;
            duration = number;
            tasks.push_back(Task(machine, duration, job));
        }
        return Job(tasks);
    }

    /**
     * function reads the file every line each.
     *
     * step 1 : look for a line which starts with " instance". The rest of the line will be stored as name
     * step 2 : skip 2 lines since there is no information.
     *          in the 3rd line look for "lower Bound" and take the number before this.
     *          After this take the number after "best known solution".
     * step 3 : In the next line are 2 numbers. The first one is the machineCount, the second one is the jobCount.
     * step 4 : every line contains Job from here to the end of the Job Shop (line starting with " +")
     *          So call getJob(line) with every line from here.
     *
     * @param instances vector which contains the instance names which will be included
     * @return vector which contains all JobShop instances which will be solved by the solver
     */
    static std::vector<JobShop> readFile(std::vector<std::string> instances) {
        std::ifstream file;
        std::string line;

        std::vector<JobShop> js;
        int machineCount;
        int jobCount;
        int bks;
        int lb;
        std::string name;

        file.open("../instance_data.txt");

        while (std::getline(file, line)) {
            //STEP 1
            if (line.substr(0, 9) != " instance") continue;

            std::vector<Job> jobs = {};
            name = line.substr(10, line.size());

            auto it = std::find(instances.begin(), instances.end(), name);
            if (instances[0] != "all" && it == instances.end())continue;

            //STEP 2
            for (int i = 0; i < 3; i++) {
                std::getline(file, line);
            }
            int strpos = line.find("lower bound");
            int length = 0;
            while (line.substr(strpos + 12 + length, 1) != ";") {
                length++;
            }

            std::stringstream bound;
            bound << line.substr(strpos + 12, length);
            bound >> lb;
            std::stringstream best;
            best << line.substr(strpos + 12 + length + 23);
            best >> bks;

            //STEP 3
            std::getline(file, line);
            std::stringstream findNumbers(line);
            int number;
            findNumbers >> number;
            jobCount = number;

            findNumbers >> number;
            machineCount = number;

            //STEP 4
            std::getline(file, line);

            int job = 0;
            while (line.substr(0, 2) != " +") {
                jobs.push_back((getJob(line, job)));
                std::getline(file, line);
                job++;
            }
            js.push_back(JobShop(jobs, machineCount, jobCount, bks, lb, name));
        }
        file.close();

        return js;
    }

    /**
     * function converts a list of ranges of instances (e.g. abz5-7) to a list of single instances (e.g. abz5, abz6, abz7)
     *
     * step 1 : check if the input is "all" or " ", then all instances gets listed.
     * step 2 : separate line by ',' and look at those substrings one by one
     * step 3 : remove all ' ' at the beginning of the line and build the instance name with each alphabetic index.
     * step 4 : take the rest of the line and remove all ' ' at the beginning and end of the line.
     *          replace '-' by ' ' so we have to numbers separated by ' '. We now can read them as words and store them as int.
     * step 5 : Since there are instances with 2 digit numbers (e.g. 01) and some with 1 digit (e.g. 1) we have to store the width of the first number.
     *          The width is used to set a fixed width of the final instance name.
     *          The function builds a struct like "prefix00" and inserts the number of the instance.
     *
     *
     * @param input string that contains the userInput to define which instances should be used
     * @return
     */
    static std::vector<std::string> parseInput(std::string input) {
        std::map<std::string, std::pair<std::string, std::string> > maxRangeMap = {
            {"abz", std::make_pair("5", "9")},
            {"dmu", std::make_pair("01", "80")},
            {"ft", std::make_pair("06", "20")},
            {"la", std::make_pair("01", "40")},
            {"orb", std::make_pair("01", "10")},
            {"swv", std::make_pair("01", "20")},
            {"ta", std::make_pair("01", "80")},
            {"yn", std::make_pair("1", "4")}
        };
        //STEP 1
        if (input == "all" || input == "") {
            return {"all"};
        }

        std::vector<std::string> instances;
        std::istringstream iss(input);
        std::string token;

        //STEP 2
        while (std::getline(iss, token, ',')) {
            //STEP 3
            std::string prefix;
            std::string startStr, endStr;

            size_t i = 0;

            while (!token.empty() && token[0] == ' ') { token = token.substr((1)); }

            while (i < token.size() && std::isalpha(token[i])) {
                prefix += token[i];
                i++;
            }

            //STEP 4
            std::string numbers = token.substr(i);

            while (!numbers.empty() && numbers[0] == ' ') { numbers = numbers.substr((1)); }
            while (!numbers.empty() && numbers.back() == ' ') { numbers.pop_back(); }

            size_t dash = numbers.find('-');
            size_t space = numbers.find(' ');

            if (dash != std::string::npos || space != std::string::npos) {
                std::replace(numbers.begin(), numbers.end(), '-', ' ');
                std::istringstream range(numbers);

                range >> startStr >> endStr;
            } else {
                startStr = endStr = numbers;
            }
            if(startStr.empty()) {
                startStr = maxRangeMap[prefix].first;
                endStr = maxRangeMap[prefix].second;
            }
            int start = std::stoi(startStr);
            int end = std::stoi(endStr);
            int width = startStr.length();
            //STEP 5
            for (int j = start; j <= end; j++) {
                std::ostringstream name;
                name << prefix;
                name.fill('0');
                name.width(width);
                name << j;
                instances.push_back(name.str());
            }
        }
        return instances;
    }
};


#endif //JSSP_SOLVER_INSTANCEHANDLER_H
