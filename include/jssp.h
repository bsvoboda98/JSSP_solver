#ifndef JOB_SHOP_SCHEDULING_ILLNER_2023_JSSP_H
#define JOB_SHOP_SCHEDULING_ILLNER_2023_JSSP_H


#include <string>
#include <vector>
#include <tuple>
#include <random>

using std::vector;
using std::string;

/**
 * internal struct for managing the operations inside the instance class
 */
struct Operation {
    int machine;
    int duration;
    int job;
};
/**
 * struct for managing solutions
 */
struct Solution {
    vector<vector<int>> solution;
    int makespan;
};
/**
 * Solution struct containing solution history of the algorithm. Designated for benchmarking
 */
struct BMResult {
    vector<vector<int>> solution;
    int makespan;
    vector<std::tuple<double,int>> history;
};

class JSSPInstance {
public:
    const vector<vector<Operation>> instance;
    const int jobCount, machineCount;
    const string filename;

    explicit JSSPInstance(string &filename): instance(readInstance(filename)), jobCount(std::get<0>(readMetrics(filename))),
                                                machineCount(std::get<1>(readMetrics(filename))),
                                                filename(filename), seed(rd()) { rng.seed(seed); };
    explicit JSSPInstance(string &filename, int seed): instance(readInstance(filename)), jobCount(std::get<0>(readMetrics(filename))),
                                             machineCount(std::get<1>(readMetrics(filename))),
                                             filename(filename), seed(seed) { rng.seed(seed); };
    explicit JSSPInstance(std::vector<std::vector<Operation>> operations, int jc, int mc): instance(operations), jobCount(jc), machineCount(mc), seed(rd()){rng.seed(seed);};

    unsigned int getSeed() {
        std::uniform_int_distribution<std::mt19937::result_type> dist(0, INT32_MAX);
        return dist(rng);
    };

    // reads a solution from file
    static Solution readSolution(string const &filename);

    // reads a solution from file
    static void writeSolutionToFile(Solution const &solution, string const &filename);

    // calculate makespan of a solution of this instance. Will not terminate, if solution is invalid.
    int calcMakespan(vector<vector<int>> const &solution) const;

    // number of operations of this instance
    [[nodiscard]] int operationCount() const;

    // repairs any invalid solution based on a random metric and returns makespan of resulting solution
    int calcMakespanAndFixSolution(vector<vector<int>> &solution, unsigned int _seed=0) const;

private:
    std::random_device rd;
    unsigned int const seed;
    std::mt19937 rng;

    // reads an instance from file
    static vector<vector<Operation>> readInstance(string &filename);

    // read first line from file (#jobs, #machines)
    static std::tuple<int,int> readMetrics(string &filename);

    // random metric for calcMakespanAndFixSolution
    void recover_solution(vector<vector<int>> &solution, vector<int> &sol_ptr, vector<int> &job_ptr, std::mt19937 &local_rnd) const;

    inline static bool contains_op(int m_no, const vector<Operation> & job) {
        for (auto op: job) {if (m_no == op.machine) {return true;}} return false;
    }
};


#endif //JOB_SHOP_SCHEDULING_ILLNER_2023_JSSP_H
