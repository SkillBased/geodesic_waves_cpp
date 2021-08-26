#include <vector>
#include <math.h>
#include <iostream>

struct Point3 {
    double x, y, z;
};

struct Face {
    std::vector <int> vertice_ids;
    std::vector <std::vector <double> > distances;

    bool operator==(const Face& other) {
        int n = this->vertice_ids.size();
        if (n != other.vertice_ids.size()) return false;
        for (int i = 0; i < n; ++i) {
            if (this->vertice_ids[i] != other.vertice_ids[i]) return false;
        }
        return true;
    }

    // return a common edge for two faces, if there is one
    // ans.second == -1 guarrantees no common edge
    // ans.second != -1 guarrantees a common edge
    std::pair<int, int> operator%(const Face& other) {
        std::pair<int, int> ans(-1, -1);
        if (*this == other) return ans;
        for (int i = 0; i < this->vertice_ids.size(); ++i) {
            for (int j = 0; j < other.vertice_ids.size(); ++j) {
                if (this->vertice_ids[i] == other.vertice_ids[j]) {
                    if (ans.first == -1) {
                        ans.first = this->vertice_ids[i];
                    } else {
                        ans.second = this->vertice_ids[i];
                        return ans;
                    }
                }
            }
        }
        return ans;
    }
};

struct Polyhedra {
    std::vector<Point3> vertices;
    std::vector<Face> faces;
};

Point3 InitPoint(double a, double b, double c) {
    Point3 pt;
    pt.x = a;
    pt.y = b;
    pt.z = c;
    return pt;
}

double DistBetween(Point3 a, Point3 b) {
    double sqared = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z);
    return sqrt(sqared);
}

Face InitFace(Polyhedra poly) {
    Face face;
    int vid;
    while (std::cin >> vid) {
        if (vid == 0) break;
        face.vertice_ids.push_back(vid - 1);
    }
    std::vector<double> vdist;
    for (int i = 0; i < face.vertice_ids.size(); ++i) {
        vdist.clear();
        for (int j = 0; j < face.vertice_ids.size(); ++j) {
            vdist.push_back(DistBetween(poly.vertices[i], poly.vertices[j]));
        }
        face.distances.push_back(vdist);
    }
    return face;
}

Polyhedra InitPolyhedra() {
    Polyhedra poly;
    int count;
    std::cin >> count;
    double a, b, c;
    for (int i = 0; i < count; ++i) {
        std::cin >> a >> b >> c;
        poly.vertices.push_back(InitPoint(a, b, c));
    }
    std::cin >> count;
    for (int i = 0; i < count; ++i) {
        poly.faces.push_back(InitFace(poly));
    }
    std::cout << "Polyhedra constructed\n";
    return poly;
}
