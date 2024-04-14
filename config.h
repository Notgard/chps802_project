#ifndef __CONFIG_H__
#define __CONFIG_H__

#define MAX(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define MIN(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

#define MIN_RAND_VAL -1E6
#define MAX_RAND_VAL 1E6

/// @brief linear system wrapper structure which contains the system matrix as a 1D array (data) 
///        and a pointer array (storage) to facilitate access to elements inside the linear system matrix
typedef struct linear_system_t {
    int nb_unknowns;
    double * data;
    double ** storage;
} linear_system_t;

#endif