#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

typedef struct {
    int dim;
    double* cords;
    int degree;
    int isSP;
} Point;

typedef struct {         // struktura pre sekvenicie daneho poctu steinerovych bodov
    int S;               // počet Steinerových bodov
    int count;           // počet validných sekvencií
    int** sequences;     // pole sekvencií
    int seqLength;       // dĺžka každej sekvencie
} SequenceSet;

void printPoints(Point* points, int N) {
    for (int i = 0; i < N; i++) {
        printf("Vrchol %d: ", i + 1);
        for (int j = 0; j < points[i].dim; j++) {
            printf("%lf ", points[i].cords[j]);
        }
        printf("\n");
    }
}
void addSteinerPoint(Point** pointsPtr, int* N, int dim, int* cords) {
    *pointsPtr = realloc(*pointsPtr, (*N + 1) * sizeof(Point));
    Point* points = *pointsPtr;

    points[*N].dim = dim;
    points[*N].cords = malloc(dim * sizeof(double));

    for (int i = 0; i < dim; i++) {
        points[*N].cords[i] = cords[i];
    }

    points[*N].degree = 0;
    points[*N].isSP = 1;
    (*N)++;
}

void freeAllPruferMemory(int** pruferSequencies, int numberOfTrees, int* sequence) {
    for (int i = 0; i < numberOfTrees; i++) {
        free(pruferSequencies[i]);
    }
    free(pruferSequencies);
    free(sequence);
}
void generateLimitedPrufer(int** pruferSequencies, int* sequence, int* counts, int position, int N, int sequenciesLength, int* indexPtr, int* limits, Point* points) {
    if (position == sequenciesLength) {
        if (pruferSequencies != NULL) { // pre prve volanie sa nezapisuju seq iba sa zistuje pocet
            for (int i = 0; i < sequenciesLength; i++) {
                pruferSequencies[*indexPtr][i] = sequence[i];
            }
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

void GenAllSteinerSequencies(Point* points, int* Nptr, int dim, SequenceSet** results, int* setCount) {
    int baseN = *Nptr;
    int maxS = baseN - 2;
    *results = malloc((maxS + 1) * sizeof(SequenceSet));
    *setCount = 0;

    for (int s = 1; s <= maxS; s++) {

        int dummyCoords[2] = {1, 1};
        addSteinerPoint(&points, Nptr, dim, dummyCoords); // aktualizuje *Nptr

        int N = *Nptr;
        int seqLength = N - 2;
        int* sequence = malloc(seqLength * sizeof(int));
        int* counts = calloc(N, sizeof(int));
        int* limits = malloc(N * sizeof(int));
        for (int i = 0; i < N; i++) {
            limits[i] = points[i].isSP ? 3 : 2;
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
/*
int main() {
    int N = 8;
    int dim = 2;
    int pointsCoordinates[8][2] = {
        {0, 0},
        {1, 0},
        {2, 1},
        {2, 2},
        {1, 3},
        {0, 3},
        {-1, 2},
        {-1, 1}
    };

    // Inicializácia pôvodných bodov
    Point* points = malloc(N * sizeof(Point));
    for (int i = 0; i < N; i++) {
        points[i].dim = dim;
        points[i].cords = malloc(dim * sizeof(double));
        points[i].cords[0] = pointsCoordinates[i][0];
        points[i].cords[1] = pointsCoordinates[i][1];
        points[i].degree = 0;
        points[i].isSP = 0;
    }

    // Zavolaj generátor sekvencií pre všetky počty Steiner bodov (1 až N-2)
    SequenceSet* results;
    int setCount;
    GenAllSteinerSequencies(points, &N, dim, &results, &setCount);


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
    limits[i] = points[i].isSP ? 2 : 1;
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
/*
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
    */
    return 0;   
}