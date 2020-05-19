#ifndef EPINETCPP_EPIMAP_H
#define EPINETCPP_EPIMAP_H

#include "gdal/gdal_priv.h"
#include "gdal/cpl_conv.h"
#include <cstdint>
#include <vector>
#include <exception>
#include "Common.h"
#include <random>
#include <algorithm>
#include <iterator>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_set_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Point_set_2<K>::Vertex_handle  Vertex_handle;
typedef K::Point_2 Point_2;

CPLErr
rasterCopy(GDALRWFlag direction, GDALRasterBand *band, uint8_t *buffer, int xOffset, int yOffset, int xSize, int ySize);

class EpiMap {

private:

    static std::mt19937 &mt()
    {
        // initialize once per thread
        thread_local static std::random_device srd;
        thread_local static std::mt19937 smt(srd());
        return smt;
    }

    static double exp_variate(double v);

    void seir_expose(event event, double epsilon, double tMax);

    void seir_infect(event event, double beta, double gamma, double tMax);

    void seir_remove(event event);

    std::map<uint32_t, uint32_t> index = {}; // coordinate index: x*y -> index in nodes vector
    std::vector<Point_2> Lr; // all points vector (points to nodes correspondence is kept in the index)
    CGAL::Point_set_2<K> PSet; // 2D pointset representation of points vector for CGAL
    std::priority_queue<event, std::vector<event>, std::greater<event>> Q = {};
    uint32_t N = 0; // nodes count
    uint8_t* buffer; // pixel buffer for direct change

public:
    size_t width;
    size_t height;
    std::vector<node> nodes = {}; // all nodes vector (the source of truth)
    explicit EpiMap(GDALRasterBand* band);
    ~EpiMap();
    uint32_t getN() const;
    std::vector<point> radiusSearch(point origin, double radius);
    void simulate(double tMax, double beta, double epsilon, double gamma);
    std::vector<node> find_contacts(node nd);
    static std::vector<node> randomSample(size_t n, const std::vector<node>& _nodes);
    void buildSearchTree();
};


#endif //EPINETCPP_EPIMAP_H
