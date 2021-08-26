#include "unfolder.h"

struct Entry {
    int vid;
    long double time;
    long double precision = 1e-4;

    bool operator<(const Entry& other) const {
        if (abs(this->time - other.time) < this->precision) return this->vid < other.vid;
        return this->time < other.time;
    }
};


struct ScanNet {
    std::vector<Mapper> map;
    int radius;
    std::vector <int> scanners;
    std::set <Entry> entries; 

    void Launch(bool debug=false) {
        int timestamp = 0;
        while (!entries.empty()) {
            Entry cur = *(entries.begin());
            entries.erase(entries.begin(), ++entries.begin());
            Mapper net = map[cur.vid];
            if (cur.time > timestamp) {
                printf("sim reached timestamp %d\n", timestamp);
                ++timestamp;
                if (debug) {
                    printf("scanner states at this point: ");
                    for (int i = 0; i < scanners.size(); ++i) {
                        printf("%5d ", scanners[i]);
                    }
                    printf("\n");
                }
            }
            scanners[cur.vid] += 1;
            std::set <long double> :: iterator callptr;
            for (int i = 0; i < net.data.size(); ++i) {
                std::set <long double> triggers = net.data[i];
                for (callptr = triggers.begin(); callptr != triggers.end(); ++callptr) {
                    Entry ping;
                    ping.vid = i;
                    ping.time = cur.time + *callptr;
                    if (ping.time <= radius) {
                        entries.insert(ping);
                    }
                }
            }
        }
        printf("simulation complete\n");
    }
};


ScanNet Prep(int start = 0, int distance = 25, long double precision = 1e-5) {
    ScanNet sim;
    sim.radius = distance;
    Polyhedra poly = InitPolyhedra();
    Mapper xmap = MapPolyhedra(poly);
    for (int i = 0; i < poly.vertices.size(); ++i) {
        sim.scanners.push_back(0);
    }
    for (int i = 0; i < poly.vertices.size(); ++i) {
        Mapper umap = Run(poly, xmap, i, distance);
        std::set <long double> stripped;
        for (int i = 0; i < umap.data.size(); ++i) {

            stripped.clear();
            long double last = 0;
            std::set <long double> :: iterator it = umap.data[i].begin();
            for (; it != umap.data[i].end(); ++it) {
                if (*it - last > precision) {
                    stripped.insert(*it);
                }
            }
            umap.data[i] = stripped;
        }
        sim.map.push_back(umap);
    }
    Entry init;
    init.time = 0;
    init.vid = start;
    sim.entries.insert(init);
    printf("simulation prepared\n");
    return sim;
}
