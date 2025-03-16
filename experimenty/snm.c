#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

typedef struct Point Point;

typedef struct {
    int count;        
    Point** neighbors; 
} SteinerNeighbors;

struct Point {
    int dim;      // Pocet dimenzii
    double *coords; // Pole suradnic
    int isSP; // 1 if Steiner point, 0 if terminal point
    SteinerNeighbors neighbors; // Neighbors of this point
};

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

double distance(Point a, Point b) {
    double sum = 0.0;
    for (int i = 0; i < a.dim; i++) {
        sum += pow(b.coords[i] - a.coords[i], 2);
    }
    return sqrt(sum);
}


void updateSteinerPoints(Point steinerPoints[], int numSteinerPoints, int dim) {
    int n = dim * numSteinerPoints;
    gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_vector *b = gsl_vector_calloc(n);
    gsl_vector *x = gsl_vector_calloc(n);
    gsl_permutation *p = gsl_permutation_alloc(n);
    int signum;

    for (int i = 0; i < numSteinerPoints; i++) {
        double sumW = 0.0;
        double beta[dim];
        for (int d = 0; d < dim; d++)
            beta[d] = 0.0;

        // vypocet "menovatel" sumW
        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point* neighbor = steinerPoints[i].neighbors.neighbors[j];
            double Dj = distance(steinerPoints[i], *neighbor);
            if (Dj < 1e-6) continue;
            sumW += 1.0 / Dj;
        }

        // Nastavenie koeficientov do matice A a vektora b
        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point* neighbor = steinerPoints[i].neighbors.neighbors[j];
            double Dj = distance(steinerPoints[i], *neighbor);
            if (Dj < 1e-6) continue;

            double alpha = (1.0 / Dj) / sumW;

            for (int d = 0; d < dim; d++) {
                int row = dim * i + d; // aktualny riadok
                if (neighbor->isSP == 1) { //zavisle od xi(steiner)
                    int col = dim * (neighbor - steinerPoints) + d; // stlpec v matici
                    gsl_matrix_set(A, row, col, alpha);
                } else { //nezavisle od xi
                    beta[d] += neighbor->coords[d] / Dj; 
                }
            }
        }

        // Nastavenie diagonaly (-1) a vektora b(coef(co neni pri xi))
        for (int d = 0; d < dim; d++) {
            gsl_vector_set(b, dim * i + d, beta[d] / sumW);
            gsl_matrix_set(A, dim * i + d, dim * i + d, -1);
        }
    }

    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, b, x);

    for (int i = 0; i < numSteinerPoints; i++) {
        for (int d = 0; d < dim; d++) {
            steinerPoints[i].coords[d] = -gsl_vector_get(x, dim * i + d);
        }
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
            Point* neighbor = steinerPoints[i].neighbors.neighbors[j];
            sum += distance(steinerPoints[i], *neighbor);
           // printf("Vzdialenosť medzi s%d a n%d: %.6f\n", i + 1, j + 1, distance(steinerPoints[i], *neighbor));
        }
    }
    sum -= distance(steinerPoints[0], steinerPoints[1]); // redundantnej hrana  S1 - S2

    return sum;
}
void iterateSteinerAlgorithm(Point steinerPoints[], Point terminals[], int numSteinerPoints, int dim, int iterations) {

    for (int i = 0; i < iterations; i++) {
        updateSteinerPoints(steinerPoints, numSteinerPoints, dim);

        double totalLength = calculateDistance(terminals, steinerPoints, numSteinerPoints);
        printf("\nIteration %d - Total length: %.6f\n", i + 1, totalLength);

        printf("\nIteration %d - Updated Steiner points:\n", i + 1);
        printSteinerPoints(steinerPoints, numSteinerPoints);
    }
}

void outputTikZ2D(Point steinerPoints[], Point terminals[], int numSteinerPoints, int numOfTerminals) {
    FILE *file = fopen("outputTikZ.txt", "w");

    // Vypis Steinerovych bodov
    for (int i = 0; i < numSteinerPoints; i++) {
        fprintf(file, "\\addvertex{S%d}{%.6f}{%.6f}\n", i + 1,
                steinerPoints[i].coords[0], steinerPoints[i].coords[1]);
    }

    // Vypis terminalnych bodov
    for (int i = 0; i < numOfTerminals; i++) {
        fprintf(file, "\\addvertex{T%d}{%.6f}{%.6f}\n", i + 1,
                terminals[i].coords[0], terminals[i].coords[1]);
    }

    // Vypis hran
    for (int i = 0; i < numSteinerPoints; i++) {
        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point *neighbor = steinerPoints[i].neighbors.neighbors[j];

            if (neighbor->isSP == 1) {
                int neighborIndex = (int)(neighbor - steinerPoints) + 1; 
                fprintf(file, "\\addedge{S%d}{S%d}\n", i + 1, neighborIndex);
            } else {
                int neighborIndex = (int)(neighbor - terminals) + 1; 
                fprintf(file, "\\addedge{S%d}{T%d}\n", i + 1, neighborIndex);
            }
        }
    }

    fclose(file);
    printf("'outputTikZ.txt' bolo vytvorene\n");
}



int main(void) {
    int terminalCount = 4;
    Point terminals[terminalCount];
    terminals[0] = createPoint(2); //T1
    terminals[0].coords[0] = 1.0;
    terminals[0].coords[1] = 3.0;

    terminals[1] = createPoint(2); //T2
    terminals[1].coords[0] = 2.0;
    terminals[1].coords[1] = 1.0;

    terminals[2] = createPoint(2); //T3
    terminals[2].coords[0] = 5.0;
    terminals[2].coords[1] = 4.0;

    terminals[3] = createPoint(2); //T4
    terminals[3].coords[0] = 4.0;
    terminals[3].coords[1] = 0.6;

    int steinerCount = 2;
    Point steinerPoints[steinerCount];
    steinerPoints[0] = createPoint(2); //S1
    steinerPoints[0].coords[0] = 3.0;
    steinerPoints[0].coords[1] = 2.0;
    steinerPoints[0].isSP = 1;

    steinerPoints[1] = createPoint(2); //S2
    steinerPoints[1].coords[0] = 4.0;
    steinerPoints[1].coords[1] = 2.0;
    steinerPoints[1].isSP = 1;

    Point* neighborsS1[] = {&terminals[0], &terminals[1], &steinerPoints[1]}; // S1 neighbors
    Point* neighborsS2[] = {&terminals[2], &terminals[3], &steinerPoints[0]}; // S2 neighbors

    steinerPoints[0].neighbors.count = 3;
    steinerPoints[0].neighbors.neighbors = neighborsS1;

    steinerPoints[1].neighbors.count = 3;
    steinerPoints[1].neighbors.neighbors = neighborsS2;

    printSteinerPoints(steinerPoints, steinerCount);
    printTerminals(terminals, terminalCount);

    iterateSteinerAlgorithm(steinerPoints, terminals, 2, 2, 3);

    outputTikZ2D(steinerPoints, terminals, 2, 4);

}
