#include <math.h>

#include "et_data.h"

double distance(Point a, Point b) {
    double sum = 0.0;
    for (int i = 0; i < a.dim; i++) {
        sum += pow(b.coords[i] - a.coords[i], 2);
    }
    return sqrt(sum);
}
