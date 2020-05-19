#ifndef EPINETCPP_COMMON_H
#define EPINETCPP_COMMON_H

enum Action {Expose, Infect, Remove};

struct point {
    double x;
    double y;
};

struct node {
    uint32_t index;
    char status;
    point loc;
};

struct event {
    double time;
    uint32_t node_index;
    Action action;

    bool operator >(const event& other) const {
        return (time > other.time);
    }
};

#endif //EPINETCPP_COMMON_H
