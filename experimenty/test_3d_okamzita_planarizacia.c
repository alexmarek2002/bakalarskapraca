#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

#include "et_data.h"
#include "et_console.h"
#include "et_tikz.h"
#include "et_smith.h"

int main(void) {
    // strom 3D
    int dim = 3;
    int terminalCount = 4;
    Point terminals[terminalCount];
    int steinerCount = 2;
    Point steinerPoints[steinerCount];

    /* terminaly vo vrcholoch stvorca [-1,+1]x[-1,+1] a v rovine z==0 */
    terminals[0] = createPoint(3); // T1
    terminals[0].coords[0] = -1.0;
    terminals[0].coords[1] = +1.0;
    terminals[0].coords[2] =  0.0;
    steinerPoints[0].isSP = 0;

    terminals[1] = createPoint(3); // T2
    terminals[1].coords[0] = -1.0;
    terminals[1].coords[1] = -1.0;
    terminals[1].coords[2] =  0.0;
    steinerPoints[0].isSP = 0;

    terminals[2] = createPoint(3); // T3
    terminals[2].coords[0] = +1.0;
    terminals[2].coords[1] = +1.0;
    terminals[2].coords[2] =  0.0;
    steinerPoints[0].isSP = 0;

    terminals[3] = createPoint(3); // T4
    terminals[3].coords[0] = +1.0;
    terminals[3].coords[1] = -1.0;
    terminals[3].coords[2] =  0.0;
    steinerPoints[0].isSP = 0;

    /* terminaly na y-ovej osi toho stvorca, symetricky wrt x, ale z==1 */
    steinerPoints[0] = createPoint(3); // S1
    steinerPoints[0].coords[0] = -0.5;
    steinerPoints[0].coords[1] =  0.0;
    steinerPoints[0].coords[2] =  1.0;
    steinerPoints[0].isSP = 1;

    steinerPoints[1] = createPoint(3); // S2
    steinerPoints[1].coords[0] = +0.5;
    steinerPoints[1].coords[1] =  0.0;
    steinerPoints[1].coords[2] = -1.0;
    steinerPoints[1].isSP = 1;


    Point* neighbours_a[] = {&terminals[0], &terminals[1], &steinerPoints[1]};
    steinerPoints[0].neighbors.count = 3;
    steinerPoints[0].neighbors.neighbors = neighbours_a;

    Point* neighbours_b[] = {&terminals[2], &terminals[3], &steinerPoints[0]};
    steinerPoints[1].neighbors.count = 3;
    steinerPoints[1].neighbors.neighbors = neighbours_b;


    // testovaci vypis
    printSteinerPoints(steinerPoints, steinerCount);
    printTerminals(terminals, terminalCount);

    // ratajme
    iterateSteinerAlgorithm(steinerPoints, terminals, steinerCount, dim, 1);

    if (dim == 2)
        outputTikZ2D(steinerPoints, terminals, steinerCount, terminalCount);
    else if (dim == 3)  
        outputTikZ3D(steinerPoints, terminals, steinerCount, terminalCount);

}
