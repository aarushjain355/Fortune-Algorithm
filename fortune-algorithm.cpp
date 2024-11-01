#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

struct Point {
    int x;
    int y;
};

struct Edge {
    Point point1;
    Point point2;
};

class FortunesAlgo {

    public:

        FortunesAlgo() : sites(std::vector<Point>()), edges(std::vector<Edge>()) {}

        int* generate_random_points() {

        }

    private:
        std::vector<Point> sites;
        std::vector<Edge> edges;

};

int main() {
    int i = 0;
    std::cout << "hELLO" << endl;
    std::cout << "y" << endl;
    return 0;
}