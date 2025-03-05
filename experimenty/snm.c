#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

typedef struct Point Point;

typedef struct {
    int count;        
    Point* neighbors; 
} SteinerNeighbors;

struct Point {
    double x; // X-coordinate
    double y; // Y-coordinate
    SteinerNeighbors neighbors; // Neighbors of this point
};

double distance(Point a, Point b) {
    return sqrt(pow(b.x - a.x, 2) + pow(b.y - a.y, 2));
}

void updateSteinerPoints(Point steinerPoints[], int numSteinerPoints) {
    int n = 2 * numSteinerPoints;
    gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_vector *b = gsl_vector_calloc(n);
    gsl_vector *x = gsl_vector_calloc(n);
    gsl_permutation *p = gsl_permutation_alloc(n);
    int signum;

    for (int i = 0; i < numSteinerPoints; i++) {
        double sumW = 0.0, sumX = 0.0, sumY = 0.0;

        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point neighbor = steinerPoints[i].neighbors.neighbors[j];
            double Dj = distance(steinerPoints[i], neighbor);
            if (Dj < 1e-6) continue; // Predchádza deleniu nulou

            double w = 1.0 / Dj;
            sumW += w;
            sumX += w * neighbor.x;
            sumY += w * neighbor.y;
        }

        gsl_matrix_set(A, 2 * i, 2 * i, sumW);
        gsl_matrix_set(A, 2 * i + 1, 2 * i + 1, sumW);
        gsl_vector_set(b, 2 * i, sumX);
        gsl_vector_set(b, 2 * i + 1, sumY);
    }

    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, b, x);

    for (int i = 0; i < numSteinerPoints; i++) {
        steinerPoints[i].x = gsl_vector_get(x, 2 * i);
        steinerPoints[i].y = gsl_vector_get(x, 2 * i + 1);
    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_permutation_free(p);
}

double calculateDistance(Point terminals[], Point steinerPoints[], int numSteinerPoints) {
    double sum = 0.0;

    for (int i = 0; i < numSteinerPoints; i++) {
        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point neighbor = steinerPoints[i].neighbors.neighbors[j];
            sum += distance(steinerPoints[i], neighbor);
            printf("Vzdialenosť medzi s%d a n%d: %.6f\n", i + 1, j + 1, distance(steinerPoints[i], neighbor));
        }
    }
    sum -= distance(steinerPoints[0],steinerPoints[1]);

    return sum;
}


void iterateSteinerAlgorithm(Point steinerPoints[], Point terminals[], int numSteinerPoints, int iterations) {
    for (int i = 0; i < iterations; i++) {
        updateSteinerPoints(steinerPoints, numSteinerPoints);
        printf("\nIterácia %d - Celková vzdialenosť: %.6f\n", i + 1, calculateDistance(terminals, steinerPoints, numSteinerPoints));

        printf("\nIterácia %d - Aktualizované Steinerove body:\n", i + 1);
        for (int j = 0; j < numSteinerPoints; j++) {
            printf("s%d = (%.6f, %.6f)\n", j + 1, steinerPoints[j].x, steinerPoints[j].y);
        }
    }
}

int main(void)
{
    Point terminals[] = {
        {1.0, 3.0}, // T1
        {2.0, 1.0}, // T2
        {5.0, 4.0}, // T3
        {4.0, 0.6}  // T4
    };

    Point steinerPoints[] = {
        {3.0, 2.0, {0, NULL}}, // S1
        {4.0, 2.0, {0, NULL}}  // S2
    };

    Point neighborsS1[] = {terminals[0], terminals[1], steinerPoints[1]}; // Susedia S1 (T1 a T2)
    Point neighborsS2[] = {terminals[2], terminals[3], steinerPoints[0]}; // Susedia S2

    SteinerNeighbors steinerConnections[] = {
        {3, neighborsS1},
        {3, neighborsS2}
    };

    steinerPoints[0].neighbors = steinerConnections[0];
    steinerPoints[1].neighbors = steinerConnections[1];

    printf("Pôvodné Steinerove body:\n");
    for (int i = 0; i < 2; i++) {
        printf("s%d = (%.6f, %.6f)\n", i + 1, steinerPoints[i].x, steinerPoints[i].y);
    }

    printf("terminal points:\n");
    for (int i = 0; i < 4; i++) {
        printf("s%d = (%.6f, %.6f)\n", i + 1, terminals[i].x, terminals[i].y);
    }

    iterateSteinerAlgorithm(steinerPoints, terminals, 2, 1);
};
