#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <thread>

using namespace std;

struct Point {
    int x;
    int y;
};

struct Edge {
    Point point1;
    Point point2;
};

// struct Parabola {

//     std::vector<Point> Parabola_points;
// };

struct Node {
    int i;
    Node* left;
    Node* right;
};

class FortunesAlgo {

    public:

        double x_min;
        double x_max;
        double y_min;
        double y_max;

        FortunesAlgo() : sites(std::vector<Point>()), edges(std::vector<Edge>()) {}

        void generate_random_points() {
            Point* arr = new Point[5];
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> distr(0.0, 10.0);
            for (int i = 0; i < 5; i++) {
                Point new_point;
                new_point.x = distr(gen);
                new_point.y = distr(gen);
                sites.push_back(new_point);
            }
            
        }

        void fortunes_algorithm() {
            
            auto max_x = std::max_element(sites.begin(), sites.end(),
                [](const Point& lhs, const Point& rhs) {
                    return lhs.x < rhs.x;
                });
            auto min_x = std::min_element(sites.begin(), sites.end(),
                [](const Point& lhs, const Point& rhs) {
                    return lhs.x > rhs.x;
                });
            auto max_y = std::max_element(sites.begin(), sites.end(),
                [](const Point& lhs, const Point& rhs) {
                    return lhs.y < rhs.y;
                });
            auto min_y = std::min_element(sites.begin(), sites.end(),
                [](const Point& lhs, const Point& rhs) {
                    return lhs.y > rhs.y;
                });

            Edge top, right, bottom, left;
            Point top_right, top_left, bottom_right, bottom_left;
            top_right.x = min_x->x - 0.2*(max_x->x - min_x->x);
            top_right.y = max_y->y + 0.2*(max_y->y - min_y->y);
            top_left.x = max_x->x + 0.2*(max_x->x - min_x->x);
            top_left.y = max_y->y + 0.2*(max_y->y - min_y->y);
            bottom_right.x = min_x->x - 0.2*(max_x->x - min_x->x);
            bottom_right.y = min_y->y - 0.2*(max_y->y - min_y->y);
            bottom_left.x = max_x->x + 0.2*(max_x->x - min_x->x);
            bottom_left.y = min_y->y - 0.2*(max_y->y - min_y->y);
            top.point1 = top_left;
            top.point2 = top_right;
            right.point1 = top_right;
            right.point2 = bottom_right;
            bottom.point1 = bottom_right;
            bottom.point2 = bottom_left;
            left.point1 = bottom_left;
            left.point2 = top_left;
            if (head != nullptr) {
                cout << "hello";
            }
            // for(size_t i = top.point1.y; i >= 0; i-=0.1) {
            //     auto target = std::find_if(sites.begin(), sites.end(), [i](const Point& d) {
            //         return d.y == i;
            //     });

            //     if (target != sites.end()) {
            //         double site_x = sites[target - sites.begin()].x;
            //         double site_y = sites[target - sites.begin()].y;
            //         Node* new_parabola;
            //         if (head == nullptr) {
            //             cout << "hELLO";
            //         }

            //         //std::thread greater_than_thread(parabola_generator, site_x, site_y, i, top_left.y, true);
            //         //std::thread less_than_thread(parabola_generator, site_x, site_y, i, top_left.y, false);
            //         //greater_than_thread.join();
            //         //less_than_thread.join();
            //         //parabola_generator(site_x, site_y, i, top_left.y, true);
            //     }
                
            // }
            

        }

        static void parabola_generator(double x_site, double y_site, double y_line, double y_boundary, bool greater_than) {

            double a = 1;
            double b = 2*x_site/(2*y_site - 2*y_line);
            double c = (x_site*x_site + y_site*y_site - y_line*y_line)/(2*y_site - 2*y_line) - y_boundary;

            if (greater_than) {
                double x_max = std::max(((-1*b + sqrt(b*b - 4*a*c))/2*a), (-1*b - sqrt(b*b - 4*a*c))/2*a);
            } else {
                double x_min = std::min(((-1*b + sqrt(b*b - 4*a*c))/2*a), (-1*b - sqrt(b*b - 4*a*c))/2*a);
            }
        }



    private:
        std::vector<Point> sites;
        std::vector<Edge> edges;
        Node* head;

};

int main(int argc, char *argv[]) {
    
   // cout << "Hello";
    FortunesAlgo object;
    object.generate_random_points();
    object.fortunes_algorithm();
    //cout << "Hello";
    return 0;
}