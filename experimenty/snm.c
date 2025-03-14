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
    double x; // X-coordinate
    double y; // Y-coordinate
    int isSP; // 1 if Steiner point, 0 if terminal point
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
        double sumW = 0.0;
        double betaX = 0.0, betaY = 0.0;
        
        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point* neighbor = steinerPoints[i].neighbors.neighbors[j];
            double Dj = distance(steinerPoints[i], *neighbor);
            if (Dj < 1e-6) continue; 
            sumW += 1.0 / Dj;
        }

        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point* neighbor = steinerPoints[i].neighbors.neighbors[j]; 

            double Dj = distance(steinerPoints[i], *neighbor);
            double alpha = (1.0 / Dj) / sumW;

            int rowX = 2 * i;     // Riadok pre x
            int rowY = 2 * i + 1; // Riadok pre y
            
            if(neighbor->isSP == 1){
                int colX = 2 * (neighbor - steinerPoints);     // Index x v A
                int colY = 2 * (neighbor - steinerPoints) + 1; // Index y v A
                gsl_matrix_set(A, rowX, colX, alpha);
                gsl_matrix_set(A, rowY, colY, alpha);
            }else{
                betaX += neighbor->x / Dj;
                betaY += neighbor->y / Dj;
            }
        }
        gsl_vector_set(b, 2 * i, betaX / sumW);
        gsl_vector_set(b, 2 * i + 1, betaY / sumW);

        gsl_matrix_set(A, 2 * i, 2 * i, -1);
        gsl_matrix_set(A, 2 * i + 1, 2 * i + 1, -1);
    }
    
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, b, x);
    
    for (int i = 0; i < numSteinerPoints; i++) {
        steinerPoints[i].x = -gsl_vector_get(x, 2 * i);
        steinerPoints[i].y = -gsl_vector_get(x, 2 * i + 1);
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
            Point* neighbor = steinerPoints[i].neighbors.neighbors[j]; // Použitie ukazovateľa
            sum += distance(steinerPoints[i], *neighbor);
            printf("Vzdialenosť medzi s%d a n%d: %.6f\n", i + 1, j + 1, distance(steinerPoints[i], *neighbor));
        }
    }
    sum -= distance(steinerPoints[0], steinerPoints[1]); // Odčítanie redundantnej hrany medzi S1 a S2

    return sum;
}
void iterateSteinerAlgorithm(Point steinerPoints[], Point terminals[], int numSteinerPoints, int iterations) {
    for (int i = 0; i < iterations; i++) {
        updateSteinerPoints(steinerPoints, numSteinerPoints);
        
        double totalLength = calculateDistance(terminals, steinerPoints, numSteinerPoints);
        printf("\nIterácia %d - Celková vzdialenosť: %.6f\n", i + 1, totalLength);

        printf("\nIterácia %d - Aktualizované Steinerove body:\n", i + 1);
        for (int j = 0; j < numSteinerPoints; j++) {
            printf("s%d = (%.6f, %.6f)\n", j + 1, steinerPoints[j].x, steinerPoints[j].y);
        }
    }
}

void outputTikZ(Point steinerPoints[], Point terminals[], int numSteinerPoints, int numOfTerminals) {
    FILE *file = fopen("outputTikZ.txt", "w");

    // Vypis vrcholov
    for (int i = 0; i < numSteinerPoints; i++) {
        fprintf(file, "\\addvertex{S%d}{%.6f}{%.6f}\n", i + 1, steinerPoints[i].x, steinerPoints[i].y);
    }
    for (int i = 0; i < numOfTerminals; i++) {
        fprintf(file, "\\addvertex{T%d}{%.6f}{%.6f}\n", i + 1, terminals[i].x, terminals[i].y);
    }

    // Vypis hran
    for (int i = 0; i < numSteinerPoints; i++) {
        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point *neighbor = steinerPoints[i].neighbors.neighbors[j];

            if (neighbor->isSP == 1) {
                int neighborIndex = (neighbor - steinerPoints) + 1; 
                fprintf(file, "\\addedge{S%d}{S%d}\n", i + 1, neighborIndex);
            } else {
                int neighborIndex = (neighbor - terminals) + 1; 
                fprintf(file, "\\addedge{S%d}{T%d}\n", i + 1, neighborIndex);
            }
        }
    }
    
    fclose(file);
    printf("'outputTikZ.txt' bolo vytvorene\n");
}


int main(void) {
    Point terminals[] = {
        {1.0, 3.0, 0, {0, NULL}}, // T1
        {2.0, 1.0, 0, {0, NULL}}, // T2
        {5.0, 4.0, 0, {0, NULL}}, // T3
        {4.0, 0.6, 0, {0, NULL}}  // T4
    };

    Point steinerPoints[] = {
        {3.0, 2.0, 1, {0, NULL}}, // S1
        {4.0, 2.0, 1, {0, NULL}}  // S2
    };

    Point* neighborsS1[] = {&terminals[0], &terminals[1], &steinerPoints[1]}; // S1 neighbors
    Point* neighborsS2[] = {&terminals[2], &terminals[3], &steinerPoints[0]}; // S2 neighbors

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

    iterateSteinerAlgorithm(steinerPoints, terminals, 2, 3);

    outputTikZ(steinerPoints, terminals, 2, 4);

}
