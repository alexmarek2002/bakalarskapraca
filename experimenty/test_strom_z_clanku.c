#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

#include "et_data.h"
#include "et_console.h"
#include "et_tikz.h"
#include "et_smith.h"

int main(void) {
// povodny elementarny strom (z clanku)     
/*
    int dim = 2;
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
*/

// o trosku komplexnejsi strom    
/*    
    int dim = 2;
    int terminalCount = 5;
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
    terminals[3].coords[0] = 6.0;
    terminals[3].coords[1] = 1.5;

    terminals[4] = createPoint(2); //T5
    terminals[4].coords[0] = 4.9;
    terminals[4].coords[1] = 0.0;

    int steinerCount = 3;
    Point steinerPoints[steinerCount];
    steinerPoints[0] = createPoint(2); //S1
    steinerPoints[0].coords[0] = 3.0;
    steinerPoints[0].coords[1] = 2.0;
    steinerPoints[0].isSP = 1;

    steinerPoints[1] = createPoint(2); //S2
    steinerPoints[1].coords[0] = 4.0;
    steinerPoints[1].coords[1] = 2.0;
    steinerPoints[1].isSP = 1;

    steinerPoints[2] = createPoint(2); //S3
    steinerPoints[2].coords[0] = 5.0;
    steinerPoints[2].coords[1] = 1.0;
    steinerPoints[2].isSP = 1;

    Point* neighborsS1[] = {&terminals[0], &terminals[1], &steinerPoints[1]}; // S1 neighbors
    Point* neighborsS2[] = {&terminals[2], &steinerPoints[2], &steinerPoints[0]}; // S2 neighbors
    Point* neighborsS3[] = {&terminals[3], &terminals[4], &steinerPoints[1]}; // S2 neighbors


    steinerPoints[0].neighbors.count = 3;
    steinerPoints[0].neighbors.neighbors = neighborsS1;

    steinerPoints[1].neighbors.count = 3;
    steinerPoints[1].neighbors.neighbors = neighborsS2;

    steinerPoints[2].neighbors.count = 3;
    steinerPoints[2].neighbors.neighbors = neighborsS3;
*/

// strom 3D
int dim = 3;
int terminalCount = 5;
Point terminals[terminalCount];
terminals[0] = createPoint(3); // T1
terminals[0].coords[0] = 1.0;
terminals[0].coords[1] = 1.0;
terminals[0].coords[2] = 2.0;

terminals[1] = createPoint(3); // T2
terminals[1].coords[0] = 2.0;
terminals[1].coords[1] = 0.0;
terminals[1].coords[2] = 1.0;

terminals[2] = createPoint(3); // T3
terminals[2].coords[0] = 2.0;
terminals[2].coords[1] = 4.0;
terminals[2].coords[2] = 0.0;

terminals[3] = createPoint(3); // T4
terminals[3].coords[0] = 6.0;
terminals[3].coords[1] = 3.0;
terminals[3].coords[2] = 0.0;

terminals[4] = createPoint(3); // T5
terminals[4].coords[0] = 5.5;
terminals[4].coords[1] = 2.0;
terminals[4].coords[2] = 3.0;

int steinerCount = 3;
Point steinerPoints[steinerCount];
steinerPoints[0] = createPoint(3); // S1
steinerPoints[0].coords[0] = 2.0;
steinerPoints[0].coords[1] = 2.0;
steinerPoints[0].coords[2] = 1.0;
steinerPoints[0].isSP = 1;

steinerPoints[1] = createPoint(3); // S2
steinerPoints[1].coords[0] = 3.0;
steinerPoints[1].coords[1] = 3.0;
steinerPoints[1].coords[2] = 1.0;
steinerPoints[1].isSP = 1;

steinerPoints[2] = createPoint(3); // S3
steinerPoints[2].coords[0] = 5.0;
steinerPoints[2].coords[1] = 3.0;
steinerPoints[2].coords[2] = 2.0;
steinerPoints[2].isSP = 1;

// neighbors
Point* neighborsS1[] = {&terminals[0], &terminals[1], &steinerPoints[1]};
Point* neighborsS2[] = {&terminals[2], &steinerPoints[2], &steinerPoints[0]};
Point* neighborsS3[] = {&terminals[3], &terminals[4], &steinerPoints[1]};

steinerPoints[0].neighbors.count = 3;
steinerPoints[0].neighbors.neighbors = neighborsS1;

steinerPoints[1].neighbors.count = 3;
steinerPoints[1].neighbors.neighbors = neighborsS2;

steinerPoints[2].neighbors.count = 3;
steinerPoints[2].neighbors.neighbors = neighborsS3;


    printSteinerPoints(steinerPoints, steinerCount);
    printTerminals(terminals, terminalCount);

    iterateSteinerAlgorithm(steinerPoints, terminals, steinerCount, dim, 100);

    if (dim == 2)
        outputTikZ2D(steinerPoints, terminals, steinerCount, terminalCount);
    else if (dim == 3)  
        outputTikZ3D(steinerPoints, terminals, steinerCount, terminalCount);

}
