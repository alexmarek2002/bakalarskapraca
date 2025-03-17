#ifndef ET_DATA
#define ET_DATA

/* ---------------------------------------------------------- data types --- */

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

/* --------------------------------------------------- symbolic constants ---*/
/*
#define ET_VERTEX_TERMINAL
#define ET_VERTEX_STEINER
#define ET_VERTEX_CENTER
*/
/* ----------------------------------------------------------- functions --- */

Point createPoint(int);
void  freePoint(Point*);

#endif
