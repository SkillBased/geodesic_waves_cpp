#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include "polyhedra.h"

const long double EPS = 1e-4;

struct Point2 {
    long double x, y, dist;
    int pid;
    bool visible;

    bool operator<(const Point2& other) {
        if (abs(this->dist - other.dist) < EPS) return this->pid < other.pid;
        return this->dist < other.dist;
    }
};

struct Line2 {
    long double a, b, c;

    long double DistToPoint(Point2 pt) {
        return (this->a * pt.x + this->b * pt.y + this->c);
    }
};


struct EdgeFold {
    bool set;
    Face first, second;

    Face unfold(Face current) {
        if (!set) return current;
        if (current == this->first) return this->second;
        return this->first;
    }
};


struct Mapper {
    std::vector <std::vector <EdgeFold> > unfolds;
    std::vector <std::set <long double> > data;

    void reset() {
        for (int i = 0; i < this->data.size(); ++i) {
            data[i].clear();
        }
    }
};

Mapper MapPolyhedra(Polyhedra poly) {
    Mapper umap;
    int n = poly.faces.size();
    int v = poly.vertices.size();
    std::vector <EdgeFold> folds;
    for (int i = 0; i < v; ++i) {
        folds.clear();
        for (int j = 0; j < v; ++j) {
            EdgeFold fold;
            fold.set = false;
            folds.push_back(fold);
        }
        umap.unfolds.push_back(folds);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            std::pair<int, int> edge = poly.faces[i] % poly.faces[j];
            if (edge.second != -1) {
                EdgeFold fold;
                fold.set = true;
                fold.first = poly.faces[i];
                fold.second = poly.faces[j];
                umap.unfolds[edge.first][edge.second] = fold;
                umap.unfolds[edge.second][edge.first] = fold;
//                printf("paired edges ");
//                for (int vid : fold.first.vertice_ids) printf("%d", vid);
//                printf(" <-> ");
//                for (int vid : fold.second.vertice_ids) printf("%d", vid);
//                printf(" over edge %d - %d\n", edge.first, edge.second);
            }
        }
    }
    std::set <long double> data_line;
    for (int i = 0; i < v; ++i) {
        data_line.clear();
        umap.data.push_back(data_line);
    }

    std::cout << "Polyhedra mapped\n";
//    for (int a = 0; a < 4; ++a) {
//        for (int b = 0; b < 4; ++b) {
//            printf("paired over edge %d - %d: ", a, b);
//            for (int vid : umap.unfolds[a][b].first.vertice_ids) printf("%d", vid);
//            printf(" <-> ");
//            for (int vid : umap.unfolds[a][b].second.vertice_ids) printf("%d", vid);
//            printf("\n");
//        }
//    }
    return umap;
}

Line2 LineFromPoints(Point2 a, Point2 b) {
    Line2 line;
    line.a = a.y - b.y;
    line.b = b.x - a.x;
    line.c = (b.x - a.x) * a.y + (a.y - b.y) * a.x;
    long double norm = sqrt(line.a * line.a + line.b * line.b);
    line.a = line.a / norm;
    line.b = line.b / norm;
    line.c = -1 * line.c / norm;
    return line;
}

std::vector<Point2> UnfoldFace(Face face, Point2 ref, Point2 a, Point2 b) {
    std::vector<Point2> constructed;
    for (int i = 0; i < face.vertice_ids.size(); ++i) {
        if (face.vertice_ids[i] == a.pid) {constructed.push_back(a); continue;}
        if (face.vertice_ids[i] == b.pid) {constructed.push_back(b); continue;}
        Point2 pt;
        pt.pid = face.vertice_ids[i];
        int apos, bpos, ppos;
        for (int i = 0; i < face.vertice_ids.size(); ++i) {
            if (face.vertice_ids[i] == a.pid) apos = i;
            if (face.vertice_ids[i] == b.pid) bpos = i;
            if (face.vertice_ids[i] == pt.pid) ppos = i;
        }
        Line2 ab = LineFromPoints(a, b);
        //align point to a and b
        long double d = face.distances[apos][bpos];
        long double da = face.distances[apos][ppos];
        long double db = face.distances[bpos][ppos];
        long double q = (da * da - db * db + d * d) / (2 * d);
        long double h = sqrt(da * da - q * q);
        long double initx = (b.x - a.x) * q / d + a.x;
        long double inity = (b.y - a.y) * q / d + a.y;
        // choose the correct intersection based on the previous face root
        pt.x = initx + h * (b.y - a.y) / d;
        pt.y = inity - h * (b.x - a.x) / d;
        if (ab.DistToPoint(ref) * ab.DistToPoint(pt) > 0) {
            pt.x = initx - h * (b.y - a.y) / d;
            pt.y = inity + h * (b.x - a.x) / d;
        }
//        printf("point %d at %.4g %.4g\n", pt.pid, pt.x, pt.y);
        pt.dist = sqrt(pt.x * pt.x + pt.y * pt.y);
        // check if the new point is visible through the unfold edge
        // and it has to be in AoB sector (as in between Ao and Bo lines)
        // it has to be on the other side of AB line from origin
        Point2 o;
        o.x = 0;
        o.y = 0;
        Line2 oa = LineFromPoints(o, a);
        Line2 ob = LineFromPoints(o, b);
        pt.visible = false;
        if (oa.DistToPoint(pt) * ob.DistToPoint(pt) < 0) {
            if (ab.DistToPoint(o) * ab.DistToPoint(pt) < 0) {
                pt.visible = true;
            }
        }
        constructed.push_back(pt);
    }
    return constructed;
}

void Unfold(Mapper& umap, Face face, std::vector <Point2> verts, int limit) {
    if (limit <= 0) return;
    Point2 ref, a, b;
    int n = verts.size();
    for (int i = 0; i < n; ++i) {
        ref = verts[i];
        a = verts[(i + 1) % n];
        b = verts[(i + 2) % n];
        if (a.visible || b.visible) {
            if (!umap.unfolds[a.pid][b.pid].set) continue;

           
            Face nface = umap.unfolds[a.pid][b.pid].unfold(face);
            std::vector <Point2> nverts = UnfoldFace(nface, ref, a, b);
            for (Point2 pt : verts) {
                if (!pt.visible || pt.pid == a.pid || pt.pid == b.pid) continue;
                umap.data[pt.pid].insert(pt.dist);
            }
            Unfold(umap, nface, nverts, limit - 1);
        }
    }
}

Mapper Run(Polyhedra poly, Mapper& umap, int vid, int unfolds=25) {
    Face init;
    bool defined = false;
    for (Face face : poly.faces) {
        for (int v : face.vertice_ids) {
            if (v == vid) {
                init = face;
                defined = true;
                break;
            }
        }
        if (defined) break;
    }
    Point2 o, a, ref;
    o.x = 0;
    o.y = 0;
    o.dist = 0;
    o.visible = true;
    o.pid = init.vertice_ids[0];
    a.x = init.distances[0][1];
    a.y = 0;
    a.dist = init.distances[0][1];
    a.visible = true;
    a.pid = init.vertice_ids[1];
    ref.x = -1;
    ref.y = -1;
    ref.dist = sqrt(2.);
    ref.visible = false;
    ref.pid = -1;
    umap.data[a.pid].insert(a.dist);
    std::vector <Point2> verts = UnfoldFace(init, ref, o, a);
    Unfold(umap, init, verts, unfolds);
    printf("Unfold for start point with id %d complete\n", vid);
    return umap;
}