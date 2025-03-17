#include <stdlib.h>

#include "et_data.h"

Point createPoint(int dim) {
    Point p;
    p.dim = dim;
    p.coords = (double *)malloc(dim * sizeof(double));
    p.isSP = 0; 
    p.neighbors.count = 0;
    p.neighbors.neighbors = NULL;
    return p;
}

void freePoint(Point *p) {
    free(p->coords);
    p->coords = NULL;

    free(p->neighbors.neighbors);
    p->neighbors.neighbors = NULL;
    p->neighbors.count = 0;

    p->dim = 0;
}
