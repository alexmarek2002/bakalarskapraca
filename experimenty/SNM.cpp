#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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


main(){
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
Point neighborsS2[] = {terminals[2], terminals[3], steinerPoints[0]}; //Susedia S2

SteinerNeighbors steinerConnections[] = {
    {3, neighborsS1},
    {3, neighborsS2}
};

steinerPoints[0].neighbors = steinerConnections[0];
steinerPoints[1].neighbors = steinerConnections[1];

int numSteinerPoints = sizeof(steinerPoints) / sizeof(steinerPoints[0]);

for(i=0;i<=numSteinerPoints;i++){
	
	int numOfAdjacentPoints = steinerPoints[i].SteinerNeighbors.count;
	for(j=0;i<numOfAdjacentPoints;i++){
		
		Sum += (steinerPoints[i].neighbors.neighbors[j+1].x)/sqrt(pow(steinerPoints[i].neighbors.neighbors[j].x - steinerPoints[i].x)+pow(steinerPoints[i].neighbors.neighbors[j].y - steinerPoints[i].y))/(1.0)/sqrt(pow(steinerPoints[i].neighbors.neighbors[j].x - steinerPoints[i].x)+pow(steinerPoints[i].neighbors.neighbors[j].y - steinerPoints[i].y))
	};
	
};	

};
