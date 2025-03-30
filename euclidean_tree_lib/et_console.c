#include <stdio.h>

#include "et_data.h"

void printPoint(Point p) {
    printf("Point(%dD, %s): [", p.dim, p.isSP ? "Steiner" : "Terminal");
    for (int i = 0; i < p.dim; i++) {
        printf("%lf", p.coords[i]);
        if (i < p.dim - 1) printf(", ");
    }
    printf("]\n");

    printf("Neighbors: %d\n", p.neighbors.count);
}

void printSteinerPoints(Point *points, int count) {
    printf("Steiner points:\n");
    for (int i = 0; i < count; i++) {
        printf("s%d = ", i + 1);
        printPoint(points[i]);
    }
}

void printTerminals(Point *points, int count) {
    printf("Initial terminal points:\n");
    for (int i = 0; i < count; i++) {
        printf("t%d = ", i + 1);
        printPoint(points[i]);
    }
}
