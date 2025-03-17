#include <gsl/gsl_linalg.h>

#include "et_console.h"
#include "et_data.h"
#include "et_geometry.h"

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
            if(neighbor->isSP == 1)
                sum += distance(steinerPoints[i], *neighbor)/2; // delenie 2 pretoze vzdialenost je pocitana 2krat
            else
                sum += distance(steinerPoints[i], *neighbor);
        }
    }
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
