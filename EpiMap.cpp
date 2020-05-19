//
// Created by atemerev on 5/10/20.
//

#include "EpiMap.h"

#include <random>

std::vector<node> EpiMap::randomSample(size_t n, const std::vector<node>& _nodes) {
    size_t len = _nodes.size();
    if (n >= len || len < 100) {
        std::vector<node> result = _nodes;
        std::shuffle(result.begin(), result.end(), std::mt19937(std::random_device()()));
        result.resize(std::min(len, n));
        return result;
    } else {
        std::vector<node> sampled;
        std::unordered_set<size_t> indexes = {};
        std::uniform_int_distribution<> distr(0, len);
        while (indexes.size() < n) {
            int r = distr(mt());
            indexes.insert(r);
        }
        for (size_t idx : indexes) {
            node node = _nodes[idx];
            sampled.push_back(node);
        }
        return sampled;
    }
}

EpiMap::EpiMap(GDALRasterBand* band) {

    mt();
    this->width = band->GetXSize();
    this->height = band->GetYSize();
    // todo create direct in-mem copy of the band, instead of the buffer
    this->buffer = (uint8_t*) CPLMalloc(sizeof(uint8_t) * width * height);

    CPLErr err = rasterCopy(GF_Read, band, buffer, 0, 0, width, height);
    if (err) {
        throw std::runtime_error("Cannot allocate raster");
    }

    for (uint32_t i = 0; i < height; i++) {
        for (uint32_t j = 0; j < width; j++) {
            uint32_t idx = i * width + j;
            if (buffer[idx] > 0) {
                // append the node to the vector of all nodes
                node n = {N, 'S', static_cast<double>(j), static_cast<double>(i)};
                nodes.push_back(n);
                // append point to the vector of points
                Point_2 point(j, i);
                Lr.push_back(point);
                index.insert({idx, N});
                N++;
                buffer[idx] = 1;
            }
        }
    }
}

EpiMap::~EpiMap() {
    CPLFree(this->buffer);
}

void EpiMap::buildSearchTree() {
    PSet.insert(Lr.begin(),Lr.end());
}

uint32_t EpiMap::getN() const {
    return N;
}

std::vector<point> EpiMap::radiusSearch(point origin, double radius) {
    Point_2 p0(origin.x, origin.y);
    Point_2 pc(origin.x + radius, origin.y);
    CGAL::Circle_2<K> rc(p0, pc);
    std::vector<Vertex_handle> LV;
    PSet.range_search(rc, std::back_inserter(LV));

    std::vector<Vertex_handle>::const_iterator it;
    std::vector<point> result;
    for (it = LV.begin();it != LV.end(); it++) {
        Point_2 point2 = (*it)->point();
        point p = {point2.x(), point2.y()};
        result.push_back(p);
    }
    return result;
}

void EpiMap::simulate(double tMax, double beta, double epsilon, double gamma) {
    std::vector<node> infected = randomSample(2, nodes);
    // assigning initial infected
    for (node n : infected) {
        event event = {0, n.index, Infect};
        Q.push(event);
    }

    while (!Q.empty()) {
        event e = Q.top();
        std::cout << e.time << ", queue: " << e.node_index << ", action: " << e.action << std::endl;
        switch (e.action) {
            case Expose:
                seir_expose(e, epsilon, tMax);
                break;
            case Infect:
                seir_infect(e, beta, gamma, tMax);
                break;
            case Remove:
                seir_remove(e);
                break;
        }
        Q.pop();
    }
}

void EpiMap::seir_expose(event e, double epsilon, double tMax) {
//    std::cout << "Expose" << std::endl;
    node& node = nodes.at(e.node_index); // a reference
    if (node.status == 'S') { // only susceptible nodes are exposed
        node.status = 'E';
        // todo update the GDAL map, update the result tallies
        double inf_time = e.time + exp_variate(epsilon);
        if (inf_time < tMax) {
            event inf_event = {inf_time, node.index, Infect};
            Q.push(inf_event);
        }
    }
}

void EpiMap::seir_infect(event e, double beta, double gamma, double tMax) {
//    std::cout << "Infect" << std::endl;
    node& nd = nodes.at(e.node_index); // a reference
    if (nd.status == 'E' || nd.status == 'S') { // only exposed nodes can be infected
        nd.status = 'I';

        // todo update the GDAL map, update the result tallies

        // inserting removal event in the future
        double rem_time = e.time + exp_variate(gamma);
        if (rem_time < tMax) {
            event rem_event = {rem_time, nd.index, Remove};
            Q.push(rem_event);
        }
        // finding susceptible contacts
        std::vector<node> contacts = find_contacts(nd);
        // exposing contacts
//        std::cout << "contacts: " << contacts.size() << std::endl;
        for (node u : contacts) {
            if (u.status == 'S') {
                double exp_time = e.time + exp_variate(beta);
//            std::cout << exp_time << " | " << e.node_index << " -> Infect: " << u.index << std::endl;
                event exp_event = {exp_time, u.index, Expose};
                Q.push(exp_event);
            }
        }
    }
}

void EpiMap::seir_remove(event e) {
//    std::cout << "Remove" << std::endl;
    node nd = nodes.at(e.node_index); // a reference
    if (nd.status == 'I') { // only infected nodes can be removed
        nd.status = 'R';
        // todo update the GDAL map, update the result tallies
    } else {
        std::cout << "Unexpected status for node " << nd.index << ": " << nd.status << std::endl;
        throw std::exception();
    }
}

std::vector<node> EpiMap::find_contacts(node nd) {
    std::vector<node> all_contacts;
    double radius = exp_variate(0.005); // todo extract constant
    std::uniform_int_distribution<> distr(0, 10);
    int num_contacts = distr(mt());
    point origin = nd.loc;
    std::vector<point> contact_coords = radiusSearch(origin, radius);
//    std::cout << "contacts: " << contact_coords.size() << std::endl;

    for (point p : contact_coords) {
//        std::cout << "x: " << p.x << ", y: " << p.y << std::endl;
        uint32_t p_idx = std::lround(p.y) * this->width + std::lround(p.x);
        uint32_t n_idx = index.at(p_idx);
        node contact = nodes.at(n_idx);
        if (contact.index != nd.index) {
            all_contacts.push_back(contact);
        }
    }
    std::vector<node> result = randomSample(num_contacts, all_contacts);
    return result;
}

double EpiMap::exp_variate(double v) {
    std::uniform_real_distribution<double> distr(0, 1);
    double u = distr(mt());
    return -std::log(u) / v;
}

CPLErr
rasterCopy(GDALRWFlag direction, GDALRasterBand *band, uint8_t *buffer, int xOffset, int yOffset, int xSize, int ySize) {
    return band->RasterIO(direction, xOffset, yOffset, xSize, ySize, buffer, xSize, ySize, GDT_Byte, 0, 0);
}
