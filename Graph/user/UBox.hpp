#include <vector>
#include <string>
#include <iostream>

#define CELL_SIZE  4
#define CELL_MIN_SIZE  2

std::string pointsAsString(std::pair<double, double> p);

class UBox
{
public:
    double center_x;
    double center_y;
    int side;
    int type;
    std::vector<int> neighbours;
    std::vector<std::pair<double, double>> points;  //нужны для упрощения формирования полигонов

    UBox(double x, double y, int s = CELL_SIZE, int t = 0)
    {
        center_x = x;
        center_y = y;
        side = s;
        type = t;
        points.push_back(std::make_pair(center_x - side / 2, center_y - side / 2));
        points.push_back(std::make_pair(center_x - side / 2, center_y + side / 2));
        points.push_back(std::make_pair(center_x + side / 2, center_y + side / 2));
        points.push_back(std::make_pair(center_x + side / 2, center_y - side / 2));
        //ext points
        points.push_back(std::make_pair(points[0].first - 1, points[0].second - 1));
        points.push_back(std::make_pair(points[1].first - 1, points[1].second + 1));
        points.push_back(std::make_pair(points[2].first + 1, points[2].second + 1));
        points.push_back(std::make_pair(points[3].first + 1, points[3].second - 1));
    }

    std::string PolygonAsString() {
        double x, y;
        std::string Polygon{ "POLYGON((" + pointsAsString(points[0]) + "," +
            pointsAsString(points[1]) + "," +
            pointsAsString(points[2]) + "," +
            pointsAsString(points[3]) + "," +
            pointsAsString(points[0]) + "))" };
        //std::cout << Polygon << std::endl;
        return Polygon;
    }

    std::string ExtPolygonAsString() {
        double x, y;
        std::string Polygon{ "POLYGON((" + pointsAsString(points[4]) + "," +
            pointsAsString(points[5]) + "," +
            pointsAsString(points[6]) + "," +
            pointsAsString(points[7]) + "," +
            pointsAsString(points[4]) + "))" };
        //std::cout << Polygon << std::endl;
        return Polygon;
    }

    void Print() {
        std::string Polygon{ "Points: " + pointsAsString(points[0]) + "," +
            pointsAsString(points[1]) + "," +
            pointsAsString(points[2]) + "," +
            pointsAsString(points[3]) };
        //std::cout << Polygon << " - type: " << type << ", side: " << side << std::endl;
    }

    void AddNeigbourByCenter(std::vector<UBox>* collection, double x, double y) {
        for (int i = 0; i < collection->size(); i++)
        {

            if (((*collection)[i].center_x == x) && ((*collection)[i].center_y == y)) {
                neighbours.push_back(i);
                (*collection)[i].neighbours.push_back(collection->size());
            }
        }
    }
};