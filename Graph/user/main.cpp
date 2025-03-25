// ConsoleApplication1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <vector>
#include <iostream>

#include "UBox.hpp"
#include "colors.h"

#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <math.h>
#include <deque>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>

#include <boost/foreach.hpp>

using namespace std;


typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygonT;

vector<UBox> graphUBox;

double maxY_calc = 0, maxX_calc = 0;

template <typename Point>
class check_coordinates
{
private:
        inline void check(double x, double * max)
    {
        if (*max < x)*max = x;;
    }

public:
    check_coordinates()
    {}

    inline void operator()(Point& p)
    {
        using boost::geometry::get;
        check(get<0>(p), &maxX_calc);
        check(get<1>(p), &maxY_calc);

    }
};


int graphCreator();
int dijkstraRunner(/*pair<int, int> start, pair<int, int> finish*/int s, int f);
int devideBox(UBox element, polygonT pMain);
int getNearestGraphNode(int x, int y);
double getDistance(int p1, int p2);

int main(int, char* []) {
    std::pair<int, int> start(1, 1);
    std::pair<int, int> finish(10, 10);
    maxX_calc = max(start.first, finish.first);
    maxY_calc = max(start.second, finish.second);
    graphCreator();
    dijkstraRunner(getNearestGraphNode(start.first, start.second), getNearestGraphNode(finish.first, finish.second));
}

int graphCreator()
{
    typedef boost::geometry::model::d2::point_xy<double> point;

    polygonT green, blue;

    boost::geometry::read_wkt(
        "POLYGON((2 2,2 8,8 8,8 2,2 2))", green);

    boost::geometry::for_each_point(green, check_coordinates<point>());

    //cout << "X: " << maxX_calc << " Y: " << maxY_calc << endl;
    //создадим сетку из больших квадратов и проверим пересечение с полигоном
    for (int x = 0; x < maxX_calc / CELL_SIZE + 1; x++) {
        for (int y = 0; y < maxY_calc / CELL_SIZE + 1; y++) {
            UBox b1{ (x + 0.5) * CELL_SIZE, (y + 0.5) * CELL_SIZE};
            polygonT p1;
            boost::geometry::read_wkt(b1.PolygonAsString(), p1);

            std::deque<polygonT> output;
            boost::geometry::intersection(green, p1, output);

            int i = 0;
            //std::cout << "green && p1:" << " (x:" << x << ", y:" << y << ")" << std::endl;
            //если у квадрата есть пересечение с полигоном то его тип - 1, если внутри полигона - тип 2
            BOOST_FOREACH(polygonT const& p, output)
            {
                //std::cout << i++ << ": " << boost::geometry::area(p) << std::endl;
                if (boost::geometry::area(p) == CELL_SIZE * CELL_SIZE) {
                    b1.type = 2;
                }
                else {
                    b1.type = 1;
                }
            }
            //добавим в качестве соседей соседние квадраты
            b1.AddNeigbourByCenter(&graphUBox, b1.center_x, b1.center_y - CELL_SIZE);
            b1.AddNeigbourByCenter(&graphUBox, b1.center_x, b1.center_y + CELL_SIZE);
            b1.AddNeigbourByCenter(&graphUBox, b1.center_x + CELL_SIZE, b1.center_y);
            b1.AddNeigbourByCenter(&graphUBox, b1.center_x - CELL_SIZE, b1.center_y);
                                   
            b1.AddNeigbourByCenter(&graphUBox, b1.center_x - CELL_SIZE, b1.center_y - CELL_SIZE);
            b1.AddNeigbourByCenter(&graphUBox, b1.center_x + CELL_SIZE, b1.center_y - CELL_SIZE);
            b1.AddNeigbourByCenter(&graphUBox, b1.center_x - CELL_SIZE, b1.center_y + CELL_SIZE);
            b1.AddNeigbourByCenter(&graphUBox, b1.center_x + CELL_SIZE, b1.center_y + CELL_SIZE);

            graphUBox.push_back(b1);
        }
    }
    //разобьем большие квадраты с типом 1 на маленькие
    for (int i = 0; i < graphUBox.size(); i++)
    {
        if (graphUBox[i].type == 1) {
            devideBox(graphUBox[i], green);
        }
        graphUBox[i].Print();
    }
    return 0;
}

int devideBox(UBox element, polygonT pMain) {
    if (element.side >= CELL_MIN_SIZE * 2) {
        int start = graphUBox.size();
        graphUBox.push_back({ element.center_x - 0.25 * element.side, element.center_y - 0.25 * element.side, element.side / 2 });
        graphUBox.push_back({ element.center_x - 0.25 * element.side, element.center_y + 0.25 * element.side, element.side / 2 });
        graphUBox.push_back({ element.center_x + 0.25 * element.side, element.center_y - 0.25 * element.side, element.side / 2 });
        graphUBox.push_back({ element.center_x + 0.25 * element.side, element.center_y + 0.25 * element.side, element.side / 2 });

        for (int i = start; i < start + 4; i++)
        {
            for (int j = start; j < start + 4; j++)
            {
                if (j != i)graphUBox[i].neighbours.push_back(j);
            }
            UBox nB = graphUBox[i];
            polygonT p1, p2, p3;
            string str1 = nB.ExtPolygonAsString();

            boost::geometry::read_wkt(nB.ExtPolygonAsString(), p1);
            boost::geometry::read_wkt(nB.PolygonAsString(), p3);
            for (int neighbour : element.neighbours) {
                string str2 = graphUBox[neighbour].PolygonAsString();
                boost::geometry::read_wkt(graphUBox[neighbour].PolygonAsString(), p2);

                std::deque<polygonT> output;
                boost::geometry::intersection(p1, p2, output);
                int inter = 0;
                BOOST_FOREACH(polygonT const& p, output)
                {
                    inter++;
                }

                if (inter > 0) {
                    graphUBox[i].neighbours.push_back(neighbour);
                    graphUBox[neighbour].neighbours.push_back(i);
                }
            }

            std::deque<polygonT> output;
            boost::geometry::intersection(p3, pMain, output);
            BOOST_FOREACH(polygonT const& p, output)
            {
                //std::cout << nB.center_x <<", "<<nB.center_y << ": " << boost::geometry::area(p) << std::endl;
                if (boost::geometry::area(p) == nB.side * nB.side) {
                    graphUBox[i].type = 2;
                }
                else {
                    graphUBox[i].type = 1;
                }
            }
        }
    }



    return 0;
}

int dijkstraRunner(/*pair<int, int> start, pair<int, int> finish*/int start, int finish)
{
    //cout << "start " << start << "finish " << finish;
    using namespace boost;

    typedef adjacency_list < listS, vecS, directedS,
        no_property, property < edge_weight_t, double > > graph_t;

    typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
    typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;

    typedef std::pair<int, int> Edge;

    const int num_nodes = graphUBox.size();

    vector<Edge> edge_vector;
    vector<double> weight_vector;

    for (int i = 0; i < num_nodes; i++) {
        if (graphUBox[i].type == 0) {
            for (int j = 0; j < graphUBox[i].neighbours.size(); j++) {
                if (graphUBox[graphUBox[i].neighbours[j]].type == 0) {
                    edge_vector.push_back(make_pair(i, graphUBox[i].neighbours[j]));
                    weight_vector.push_back(getDistance(i, graphUBox[i].neighbours[j]));
                }
            }
        }
    }
    
    Edge edge_array[1000];
    double* weights = new double[weight_vector.size()];

    for (int i = 0; i < 1000; i++)
    {
        if (i < edge_vector.size()) {
            edge_array[i] = edge_vector[i];
            weights[i] = weight_vector[i];
        }
        else {
            edge_array[i] = Edge(0, 0);
        }
    }

    int num_arcs = sizeof(edge_array) / sizeof(Edge);

    graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
    property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);

    std::vector<vertex_descriptor> p(num_vertices(g));
    std::vector<double> d(num_vertices(g));
    vertex_descriptor s = vertex(start, g);//vertex_descriptor s = vertex(21, g);

    dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));

   /* std::cout << "distances and parents from: " << start << " (" << start * 4 / STD_CELL_CNT + 2 << "." << (start % STD_CELL_CNT) * 4 + 2 << ")"
        << " to " << finish << " (" << finish * 4 / STD_CELL_CNT + 2 << "." << (finish % STD_CELL_CNT) * 4 + 2 << ")" << std::endl;
        */
    graph_traits < graph_t >::vertex_iterator vi, vend;
    int parent = p[finish];
    vector<int>path;
    path.push_back(finish);
    while (parent != start) {
        path.push_back(parent);
        parent = p[parent];
    }
    path.push_back(parent);
    for (int i = path.size() - 1; i > 0; i--) {
        std::cout << "(" << graphUBox[path[i]].center_x << " " << graphUBox[path[i]].center_y << ")->";
    }
    std::cout << "(" << graphUBox[path[0]].center_x << " " << graphUBox[path[0]].center_y << ") - " << d[finish] << endl;

    for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
        int node = (*vi);
        /*if (node == finish) {
            std::cout << BLUE;
        }*/
        std::cout << "distance(" << (*vi) << ") = " << d[*vi] << ", ";
        std::cout << "parent(" << (*vi) << ") = " << p[*vi] << std::
            endl;
        std::cout << RESET;
    }
    std::cout << std::endl;

    return EXIT_SUCCESS;
}

int getNearestGraphNode(int x, int y) {
    double minDistance = 1000, distance, a, b;
    int vertex = 0;
    for (int i = 0; i < graphUBox.size(); i++) {
        a = x - graphUBox[i].center_x;
        b = y - graphUBox[i].center_y;
        distance = sqrt(a * a + b * b);
        if (distance < minDistance) {
            minDistance = distance;
            vertex = i;
        }
    }
    return vertex;
}

double getDistance(int p1, int p2) {
    double a = graphUBox[p2].center_x - graphUBox[p1].center_x;
    double b = graphUBox[p2].center_y - graphUBox[p1].center_y;
    return sqrt(a * a + b * b);
}

string pointsAsString(pair<double, double> p) {
    std::string xStr = std::to_string(p.first);
    std::string yStr = std::to_string(p.second);
    std::string ResultString{ xStr + " " + yStr };
    return ResultString;
}