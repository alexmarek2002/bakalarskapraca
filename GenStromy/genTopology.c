#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>    
#include <gsl/gsl_linalg.h>

typedef struct Point Point;

typedef struct {
    int count;        
    Point** neighbors; 
} NeighborList;

struct Point {
    int dim;      // Pocet dimenzii
    double *coords; // Pole suradnic
    int degree;
    int isSP; // 1 if Steiner point, 0 if terminal point
    int indexInSteiner;
    NeighborList neighbors; // Neighbors of this point
};

typedef struct {         // struktura pre sekvenicie daneho poctu steinerovych bodov
    int S;               // počet Steinerových bodov
    int count;           // počet validných sekvencií
    int** sequences;     // pole sekvencií
    int seqLength;       // dĺžka každej sekvencie
} SequenceSet;

Point createPoint(int dim) {
    Point p;
    p.dim = dim;
    p.coords = (double *)malloc(dim * sizeof(double));
    p.isSP = 0; 
    p.neighbors.count = 0;
    p.indexInSteiner = -1; 
    p.neighbors.neighbors = NULL;
    return p;
}
Point* copyPoints(Point* src, int N) {
    Point* dst = malloc(N * sizeof(Point));
    for (int i = 0; i < N; i++) {
        dst[i].dim = src[i].dim;
        dst[i].coords = malloc(dst[i].dim * sizeof(double));
        for (int d = 0; d < dst[i].dim; d++) {
            dst[i].coords[d] = src[i].coords[d];
        }
        dst[i].degree = src[i].degree;
        dst[i].isSP = src[i].isSP;
        dst[i].indexInSteiner = src[i].indexInSteiner;

        // neighbors kópiu nepotrebujeme do detailu, stačí count=0, 
        dst[i].neighbors.count = 0;
        dst[i].neighbors.neighbors = NULL;
    }
    return dst;
}
void freeClonedPoints(Point* pts, int N) {
    if (!pts) return;
    for (int i = 0; i < N; i++) {
        free(pts[i].coords);
    }
    free(pts);
}
void initNeighborList(Point* p, int maxNeighbors) {
    p->neighbors.neighbors = malloc(maxNeighbors * sizeof(Point*));
    p->neighbors.count = 0;
}
void addNeighbor(Point* p, Point* n) {
    p->neighbors.neighbors[p->neighbors.count++] = n;
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
void assignSteinerIndices(Point points[], int N) {
    int idx = 0;
    for (int i = 0; i < N; i++) {
        if (points[i].isSP == 1) {
            points[i].indexInSteiner = idx++;
        } else {
            points[i].indexInSteiner = -1;
        }
    }
}
/*
void updateSteinerPoints(Point points[], int N, int dim) {

    int numSteinerPoints = 0;
    for (int i = 0; i < N; i++) {
        if (points[i].isSP == 1) numSteinerPoints++;
    }

    if (numSteinerPoints == 0) return; 
    printf("%d", numSteinerPoints);
    int n = dim * numSteinerPoints;

    gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_vector *b = gsl_vector_calloc(n);
    gsl_vector *x = gsl_vector_calloc(n);
    gsl_permutation *p = gsl_permutation_alloc(n);
    int signum;

    for (int i = 0; i < N; i++) {
        if (!points[i].isSP) continue; 
        int si = points[i].indexInSteiner;
        if (si < 0) continue;

        double sumW = 0.0;
        double beta[dim];
        for (int d = 0; d < dim; d++)
            beta[d] = 0.0;

        // vypocet "menovatel" sumW
        for (int j = 0; j < points[i].neighbors.count; j++) {
            Point* neighbor = points[i].neighbors.neighbors[j];
            double Dj = distance(points[i], *neighbor);
            if (Dj < 1e-6) continue;
            sumW += 1.0 / Dj;
        }

        // Nastavenie koeficientov do matice A a vektora b
        for (int j = 0; j < points[i].neighbors.count; j++) {
            Point* neighbor = points[i].neighbors.neighbors[j];
            double Dj = distance(points[i], *neighbor);
            if (Dj < 1e-6) continue;

            double alpha = (1.0 / Dj) / sumW;

            for (int d = 0; d < dim; d++) {
                int row = dim * i + d; // aktualny riadok
                if (neighbor->isSP == 1) { //zavisle od xi(steiner)
                    int colSteiner = neighbor->indexInSteiner;
                    if (colSteiner >= 0) {
                        int col = dim * colSteiner + d;
                        gsl_matrix_set(A, row, col, alpha);
                    }
                } else { //nezavisle od xi
                    beta[d] += neighbor->coords[d] / Dj; 
                }
            }
        }

        // Nastavenie diagonaly (-1) a vektora b(coef(co neni pri xi))
        for (int d = 0; d < dim; d++) {
            int row = dim * i + d;
            gsl_vector_set(b, row, beta[d] / sumW);
            gsl_matrix_set(A, row, row, -1.0);
        }
    }

    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, b, x);

    for (int i = 0; i < N; i++) {
        if (!points[i].isSP) continue;
        int i = points[i].indexInSteiner;
        if (i < 0) continue;
        for (int d = 0; d < dim; d++) {
            points[i].coords[d] = -gsl_vector_get(x, dim * i + d);
        }
    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_permutation_free(p);
}
    */
   void updateSteinerPoints(Point points[], int N, int dim) {

    int numSteinerPoints = 0;
    for (int i = 0; i < N; i++) {
        if (points[i].isSP == 1) {
            numSteinerPoints++;
        }
    }
    if (numSteinerPoints == 0) return; 

    int n = dim * numSteinerPoints;

    gsl_matrix *A = gsl_matrix_calloc(n, n);
    gsl_vector *b = gsl_vector_calloc(n);
    gsl_vector *x = gsl_vector_calloc(n);
    gsl_permutation *p = gsl_permutation_alloc(n);
    int signum;

    for (int i = 0; i < N; i++) {
        if (!points[i].isSP) continue; 
        int si = points[i].indexInSteiner;   
        if (si < 0) continue;

        double sumW = 0.0;
        double beta[dim];
        for (int d = 0; d < dim; d++) {
            beta[d] = 0.0;
        }

        for (int j = 0; j < points[i].neighbors.count; j++) {
            Point* neighbor = points[i].neighbors.neighbors[j];
            double Dj = distance(points[i], *neighbor);
            if (Dj < 1e-8) Dj = 1e-8; 
            sumW += 1.0 / Dj;
        }

        for (int j = 0; j < points[i].neighbors.count; j++) {
            Point* neighbor = points[i].neighbors.neighbors[j];
            double Dj = distance(points[i], *neighbor);
            if (Dj < 1e-8) Dj = 1e-8;

            double alpha = (1.0 / Dj) / sumW;

            for (int d = 0; d < dim; d++) {
                int row = dim * si + d;      

                if (neighbor->isSP == 1) {
                    int colSteiner = neighbor->indexInSteiner;
                    if (colSteiner >= 0) {
                        int col = dim * colSteiner + d;
                        gsl_matrix_set(A, row, col, alpha);
                    }
                } else {
                    beta[d] += neighbor->coords[d] / Dj; 
                }
            }
        }

        for (int d = 0; d < dim; d++) {
            int row = dim * si + d;          
            gsl_vector_set(b, row, beta[d] / sumW);
            gsl_matrix_set(A, row, row, -1.0);
        }
    }

    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, b, x);

    for (int i = 0; i < N; i++) {
        if (!points[i].isSP) continue;
        int si = points[i].indexInSteiner;   
        if (si < 0) continue;

        for (int d = 0; d < dim; d++) {
            // Použijeme index si (nie i!)
            points[i].coords[d] = -gsl_vector_get(x, dim * si + d);
        }
    }

    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_permutation_free(p);
}
 
double calculateDistance(Point points[], int N) {
    double sum = 0.0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < points[i].neighbors.count; j++) {
            Point* neighbor = points[i].neighbors.neighbors[j];
            double d = distance(points[i], *neighbor);
            if (points[i].isSP && neighbor->isSP)
                sum += d / 2.0;
            else
                sum += d;
        }
    }
    return sum;
}
double iterateSteinerAlgorithm(Point points[], int N, int dim, int iterations) {
    double totalLength = 0.0;
    for (int i = 0; i < iterations; i++) {

        assignSteinerIndices(points, N);
        updateSteinerPoints(points, N, dim);

        totalLength = calculateDistance(points, N);
        printf("\nIteration %d - Total length: %.6f\n", i + 1, totalLength);

        printf("Iteration %d - Updated Steiner points:\n", i + 1);
        for (int p = 0; p < N; p++) {
            if (points[p].isSP) {
                printf("S[%d]: ", p);
                printPoint(points[p]);
            }
        }
    }
    return totalLength;
}

void pruferToTree(int* prufer, int length, int edges[][2]) {
    int n = length + 2;
    int degree[n + 1];  

    for (int i = 1; i <= n; i++) {
        degree[i] = 1;
    }

    for (int i = 0; i < length; i++) {
        degree[prufer[i]]++;
    }

    for (int i = 0; i < length; i++) {
        int leaf = -1;

        for (int j = 1; j <= n; j++) {
            if (degree[j] == 1) {
                leaf = j;
                break;
            }
        }

        edges[i][0] = leaf;
        edges[i][1] = prufer[i];

        degree[leaf]--;
        degree[prufer[i]]--;
    }

    int a = -1, b = -1;
    for (int i = 1; i <= n; i++) {
        if (degree[i] == 1) {
            if (a == -1) a = i;
            else b = i;
        }
    }

    edges[length][0] = a;
    edges[length][1] = b;

    int degreeCheck[n + 1];
    for (int i = 0; i <= n; i++) degreeCheck[i] = 0;

    for (int i = 0; i < n - 1; i++) {
        degreeCheck[edges[i][0]]++;
        degreeCheck[edges[i][1]]++;
    }
}

void printPoints(Point* points, int N) {
    for (int i = 0; i < N; i++) {
        printf("Vrchol %d: ", i + 1);
        for (int j = 0; j < points[i].dim; j++) {
            printf("%lf ", points[i].coords[j]);
        }
        printf("\n");
    }
}
void addSteinerPoint(Point** pointsPtr, int* N, int dim, double* cords) {
    *pointsPtr = realloc(*pointsPtr, (*N + 1) * sizeof(Point));
    Point* points = *pointsPtr;

    points[*N].dim = dim;
    points[*N].coords = malloc(dim * sizeof(double));

    for (int i = 0; i < dim; i++) {
        points[*N].coords[i] = cords[i];
    }

    points[*N].degree = 0;
    points[*N].isSP = 1;
    points[*N].indexInSteiner = -1;
    points[*N].neighbors.count = 0;
    points[*N].neighbors.neighbors = NULL;

    (*N)++;
}

void freeAllPruferMemory(int** pruferSequencies, int numberOfTrees, int* sequence) {
    for (int i = 0; i < numberOfTrees; i++) {
        free(pruferSequencies[i]);
    }
    free(pruferSequencies);
    free(sequence);
}
// rekurzivna funkcia na generovane pruferovych sekvencii s omedzenim ze kazdy bod (cislo 1-N) sa moze nachadzat nanjvys 2 krat
void generateLimitedPrufer(int** pruferSequencies, int* sequence, int* counts, int position, int N, int sequenciesLength, int* indexPtr, int* limits, Point* points) {
    if (position == sequenciesLength) { // toto je vlastne hranica pokial mame rekurzivne volat funkciu ... potialto sa iteruje
        if (pruferSequencies != NULL) { // pre prve volanie sa nezapisuju seq iba sa zistuje pocet
            for (int i = 0; i < sequenciesLength; i++) {
                pruferSequencies[*indexPtr][i] = sequence[i];
            }
            // (***) tuto by som pridal funkcionalitu z filterValidSequences, tak, ze by sa presla sekvencia zratali sa jedntlive cisla a ak (isSP==1 pocet==2) tak sa zapisu len tie .
            // to by usetrilo pamat a vyrazne, (implementujem ked nebude robit treba dolezitejsie veci alebo sa stane implementacia dolezitejsou ako ostante veci)
        }
        (*indexPtr)++;
        return;
    }

    for (int i = 0; i < N; i++) {
        if (counts[i] < limits[i]) {
            sequence[position] = i + 1; // hodnoty v Prüfer kóde sú od 1
            counts[i]++;
            generateLimitedPrufer(pruferSequencies, sequence, counts, position + 1, N, sequenciesLength, indexPtr, limits, points);
            counts[i]--;
        }
    }
}
// FILTER pre vygenerovane pruferove sekvencie
// cize funkcia ktora prejde vsetky sekvenice a necha len tie kde isSP==1 su 2 krat v sekvencii, maju 3 hrany
// ak bude kod moc pomaly tak na miesto kde su hviezdicky (***)(generateLimitedPrufer()) pridam tuto funkcionalitu a nebudem to muset prechadzat cele (neviem preco som to takto vymyslel teraz mi to pride hlupe)
void filterValidSequences(int*** pruferSequencies, int* count, int seqLength, int N, Point* points) {
    int** input = *pruferSequencies;
    int** filtered = malloc(*count * sizeof(int*));
    int validCount = 0;

    for (int i = 0; i < *count; i++) {
        int tempCounts[N];
        memset(tempCounts, 0, sizeof(tempCounts));

        for (int j = 0; j < seqLength; j++) {
            int val = input[i][j];
            if (val >= 1 && val <= N) {
                tempCounts[val - 1]++;
            }
        }

        bool valid = true;
        for (int j = 0; j < N; j++) {
            if (points[j].isSP && tempCounts[j] != 2) {
                valid = false;
                break;
            }
        }

        if (valid) {
            filtered[validCount] = malloc(seqLength * sizeof(int));
            memcpy(filtered[validCount], input[i], seqLength * sizeof(int));
            validCount++;
        }
    }

    // Free original sequences
    for (int i = 0; i < *count; i++) {
        free(input[i]);
    }
    free(input);

    // Resize and assign back
    *pruferSequencies = realloc(filtered, validCount * sizeof(int*));
    *count = validCount;
}
// uklada do pamate (neprepisoval som funkciu lebo sa moze este zist)
void GenAllSteinerSequencies(Point* points, int* Nptr, int dim, SequenceSet** results, int* setCount) {
    int baseN = *Nptr;
    int maxS = baseN - 2;
    *results = malloc((maxS + 1) * sizeof(SequenceSet));
    *setCount = 0;

    for (int s = 1; s <= maxS; s++) { // takze postupne pridavame 1 sp (steiner point), 2 sp, 3 sp, ... n-2 sp a pre kazdy pocet sp generujeme vsetky sekvencie
                                                    // ak budeme neskor uvazovat len pripady kedy vsetky terminaly==listy tak nebude treba iterat, lebo pre taketo stromy je jedina pripustna
                                                    // topologia - full steiner topology - N-2 sp 
        double dummyCoords[2] = {1, 1};
        addSteinerPoint(&points, Nptr, dim, dummyCoords); // aktualizuje *Nptr
//prerob - seg fault-uje to...
        int N = *Nptr;
        int seqLength = N - 2;
        int* sequence = malloc(seqLength * sizeof(int));
        int* counts = calloc(N, sizeof(int));
        int* limits = malloc(N * sizeof(int));
        for (int i = 0; i < N; i++) {
            limits[i] = points[i].isSP ? 2 : 1;
        }

        int count = 0;
        generateLimitedPrufer(NULL, sequence, counts, 0, N, seqLength, &count, limits, points);

        int** pruferSequences = malloc(count * sizeof(int*));
        for (int i = 0; i < count; i++) {
            pruferSequences[i] = malloc(seqLength * sizeof(int));
        }
        memset(counts, 0, N * sizeof(int));
        int index = 0;
        generateLimitedPrufer(pruferSequences, sequence, counts, 0, N, seqLength, &index, limits, points);

        filterValidSequences(&pruferSequences, &count, seqLength, N, points);

        (*results)[s].S = s;
        (*results)[s].count = count;
        (*results)[s].seqLength = seqLength;
        (*results)[s].sequences = pruferSequences;
        (*setCount)++;

        printf("S = %d → %d platných sekvencií\n", s, count);

        free(sequence);
        free(counts);
        free(limits);
    }
}
// neuklda do pamate 
void GenAllSteinerSequenciesP(Point* points, int* Nptr, int dim)
{
    int baseN = *Nptr;
    int maxS = baseN - 2;
    double shortestLength = INFINITY;
    Point* bestPoints = NULL;  
    int (*bestEdges)[2] = NULL;
    int bestN = 0;  

    srand(time(NULL));
    for (int s = 1; s <= maxS; s++) {
        double* dummyCoords = malloc(dim * sizeof(double));
        for (int d = 0; d < dim; d++) {
            dummyCoords[d] = 0.0;
            for (int p = 0; p < baseN; p++) { 
                dummyCoords[d] += points[p].coords[d];
            }
            dummyCoords[d] /= baseN;
        
        }

        addSteinerPoint(&points, Nptr, dim, dummyCoords);
        free(dummyCoords);

        int N = *Nptr;
        int seqLength = N - 2;

        int* sequence = malloc(seqLength * sizeof(int));
        int* counts   = calloc(N, sizeof(int));
        int* limits   = malloc(N * sizeof(int));

        for (int i = 0; i < N; i++) {
            limits[i] = points[i].isSP ? 2 : 1;
        }

        int count = 0;
        generateLimitedPrufer(NULL, sequence, counts, 0, N, seqLength, &count, limits, points);

        int** pruferSequences = malloc(count * sizeof(int*));
        for (int i = 0; i < count; i++) {
            pruferSequences[i] = malloc(seqLength * sizeof(int));
        }

        memset(counts, 0, N * sizeof(int));
        int index = 0;
        generateLimitedPrufer(pruferSequences, sequence, counts, 0, N, seqLength, &index, limits, points);

        filterValidSequences(&pruferSequences, &count, seqLength, N, points);

        
        printf("\nS = %d → %d platných sekvencií. Vypisujem hrany:\n", s, count);

        for (int i = 0; i < count; i++) {
            printf("\nSekvencia %d \n", i + 1);
            int edges[N - 1][2];
            pruferToTree(pruferSequences[i], seqLength, edges);

            for (int p = 0; p < N; p++) {
                if (points[p].neighbors.neighbors != NULL) {
                    free(points[p].neighbors.neighbors); 
                }
                initNeighborList(&points[p], N - 1);
            }

            for (int e = 0; e < N-1; e++) {
                int k = edges[e][0]-1;
                int l = edges[e][1]-1;
                addNeighbor(&points[k], &points[l]);
                addNeighbor(&points[l], &points[k]);
            }
            int iterations = 30;
            double totalLength = iterateSteinerAlgorithm(points, N, dim, iterations);
            if (totalLength < shortestLength) {
                shortestLength = totalLength;
            
                if (bestPoints) {
                    freeClonedPoints(bestPoints, bestN);
                    bestPoints = NULL;
                }
                if (bestEdges) {
                    free(bestEdges);
                    bestEdges = NULL;
                }

                bestN = N;
                bestPoints = copyPoints(points, N);
                bestEdges = malloc( (N - 1) * sizeof(*bestEdges) );
            
                for (int e = 0; e < (N - 1); e++) {
                    bestEdges[e][0] = edges[e][0];
                    bestEdges[e][1] = edges[e][1];
                }
            }   
        }

        for (int i = 0; i < count; i++) {
            free(pruferSequences[i]);
        }
        free(pruferSequences);
        free(sequence);
        free(counts);
        free(limits);
    }
    printf("\n=========== Najlepšie riešenie ===========\n");
    printf("Najkratšia dĺžka = %.6f\n", shortestLength);

    if (bestPoints) {
        printf("Počet vrcholov = %d\n", bestN);
        for (int i = 0; i < bestN; i++) {
            printf("Point[%d]: (", i);
            for (int d = 0; d < bestPoints[i].dim; d++) {
                printf("%lf", bestPoints[i].coords[d]);
                if (d < bestPoints[i].dim - 1) printf(", ");
            }
            printf(")  %s\n", bestPoints[i].isSP ? "Steiner" : "Terminal");
        }
    
        if (bestEdges) {
            printf("Hrany:\n");
            for (int e = 0; e < bestN - 1; e++) {
                printf("  %d -- %d\n", bestEdges[e][0], bestEdges[e][1]);
            }
            free(bestEdges);
        }
    
        freeClonedPoints(bestPoints, bestN);
    }
}


int main() {
    gsl_set_error_handler_off();
/*
    int N = 4;
    int dim = 2;
    double pointsCoordinates[4][2] = {
        {1.0, 3.0},
        {2.0, 1.0},
        {5.0, 4.0},
        {4.0, 0.6},
    };
*/

    int N = 6;
    int dim = 3;
    double pointsCoordinates[6][3] = {
        { 1.0,  0.0,  0.0},
        {-1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        { 0.0, -1.0,  0.0},
        { 0.0,  0.0,  1.0},
        { 0.0,  0.0, -1.0}
    };

    Point* points = malloc(N * sizeof(Point));
    for (int i = 0; i < N; i++) {
        points[i].dim = dim;
        points[i].coords = malloc(dim * sizeof(double));
        for (int d = 0; d < dim; d++) {
            points[i].coords[d] = pointsCoordinates[i][d];
        }
        points[i].degree = 0;
        points[i].isSP = 0; 
        points[i].indexInSteiner = -1;
        points[i].neighbors.count = 0;
        points[i].neighbors.neighbors = NULL;
    }

    GenAllSteinerSequenciesP(points, &N, dim);

    for (int i = 0; i < N; i++) {
        free(points[i].coords);
        free(points[i].neighbors.neighbors);
    }
    free(points);


    /*
    // Zavolaj generátor sekvencií pre všetky počty Steiner bodov (1 až N-2)
    SequenceSet* results;
    int setCount;
    GenAllSteinerSequencies(points, &N, dim, &results, &setCount);

    for (int i = 1; i <= setCount; i++) {
        printf("\n==== Výpis sekvencií pre S = %d ====\n", results[i].S);
        printf("Počet platných sekvencií: %d\n", results[i].count);

        for (int j = 0; j < results[i].count; j++) {
            printf("  Sekvencia %d: ", j + 1);
            for (int k = 0; k < results[i].seqLength; k++) {
                printf("%d ", results[i].sequences[j][k]);
            }
            printf("\n");
        }
    }
 */       
/*
    // Uvoľnenie pamäte za všetky sekvencie a pole results
    for (int i = 1; i <= setCount; i++) {
        for (int j = 0; j < results[i].count; j++) {
            free(results[i].sequences[j]);
        }
        free(results[i].sequences);
    }
    free(results);

    // Uvoľnenie pamäte za všetky pôvodné (aj pridané) body
    for (int i = 0; i < N; i++) {
        free(points[i].cords);
    }
    free(points);
*/
    return 0;
}
/*
int main(){
    int N = 8;
    int dim = 2;
    int pointsCordinates[8][2] = {
        {0, 0},
        {1, 0},
        {2, 1},
        {2, 2},
        {1, 3},
        {0, 3},
        {-1,2},
        {-1,1}
    };

    Point* points = (Point*)malloc(N * sizeof(Point));
    for (int i = 0; i < N; i++) {
        points[i].dim = dim;
        points[i].cords = (double*)malloc(dim * sizeof(double));
        points[i].cords[0] = pointsCordinates[i][0];
        points[i].cords[1] = pointsCordinates[i][1];
        points[i].degree = 0;
        points[i].isSP = 0;
    }
    
    int newPoint1[2] = {1, 1};
    addSteinerPoint(&points, &N, dim, newPoint1);
    int newPoint2[2] = {1.2, 1.2};
    addSteinerPoint(&points, &N, dim, newPoint2);
    int newPoint3[2] = {1, 1.2};
    addSteinerPoint(&points, &N, dim, newPoint3);

    printPoints(points, N);
//////////////////////////////////////////////////////

int* limits = malloc((10) * sizeof(int));
for (int i = 0; i <= 9; i++) {
 //   limits[i] = points[i].isSP ? 2 : 1;
    limits[i] = 2;
}
printf("Limits: \n");
for (int i = 0; i <= 8; i++) {
    printf("%d.  ",i);
    printf("%d ", limits[i]);
    printf("isPS: %d\n", points[i].isSP);
}
//////////////////////////////////////////////////////
    int sequencesLength = N - 2;
    int* sequence = malloc(sequencesLength * sizeof(int));
    int* counts = calloc(N, sizeof(int)); // pomocna premena na pocitanie kolko krat sa dany prvok vyskytuje v sekvencii

    int pocet = 0; // index sekvencie ktoru generujeme 
    generateLimitedPrufer(NULL, sequence, counts, 0, N, sequencesLength, &pocet, limits, points); // 1. faza
    printf("Počet vygenerovaných sekvencií: %d\n", pocet);

    int** pruferSequencies = malloc(pocet * sizeof(int*));//2. faza
    for (int i = 0; i < pocet; i++) {
        pruferSequencies[i] = malloc(sequencesLength * sizeof(int));
    } 

    memset(counts, 0, N * sizeof(int));
    int index = 0;
    generateLimitedPrufer(pruferSequencies, sequence, counts, 0, N, sequencesLength, &index, limits, points);

    filterValidSequences(&pruferSequencies, &pocet, sequencesLength, N, points);
    printf("Počet platných sekvencií: %d\n", pocet);

    for (int i = 0; i < pocet; i++) {
        printf("Sekvencia %d: ", i + 1);
        for (int j = 0; j < sequencesLength; j++) {
            printf("%d ", pruferSequencies[i][j]);
        }
        printf("\n");
    }

    freeAllPruferMemory(pruferSequencies, pocet, sequence);
    for (int i = 0; i < N; i++) free(points[i].cords);
    free(points);
    free(counts);
    free(limits);
    
    return 0;   
}
*/