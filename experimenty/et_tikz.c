#include <stdio.h>
#include <stdlib.h>

#include "et_data.h"

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

void outputTikZ3D(Point steinerPoints[], Point terminals[], int numSteinerPoints, int numOfTerminals) {
    FILE *file = fopen("outputTikZ3D.txt", "w");

    // terminal
    fprintf(file, "%% Terminály T1-T%d\n", numOfTerminals);
    fprintf(file, "\\addplot3[\n\tonly marks,\n\tmark=*,\n\tscatter,\n\tscatter src=explicit symbolic\n] coordinates {\n");
    for (int i = 0; i < numOfTerminals; i++) {
        fprintf(file, "    (%.6f, %.6f, %.6f) [terminal]\n",
                terminals[i].coords[0], terminals[i].coords[1], terminals[i].coords[2]);
    }
    fprintf(file, "};\n\n");

    //popis
    fprintf(file, "%% Popisy terminálov\n");
    for (int i = 0; i < numOfTerminals; i++) {
        fprintf(file, "\\node[anchor=south east, blue] at (axis cs:%.6f, %.6f, %.6f) {T%d};\n",
                terminals[i].coords[0], terminals[i].coords[1], terminals[i].coords[2], i + 1);
    }
    fprintf(file, "\n");

    // Steiner
    fprintf(file, "%% Steinerové body S1-S%d\n", numSteinerPoints);
    fprintf(file, "\\addplot3[\n\tonly marks,\n\tmark=square*,\n\tscatter,\n\tscatter src=explicit symbolic\n] coordinates {\n");
    for (int i = 0; i < numSteinerPoints; i++) {
        fprintf(file, "    (%.6f, %.6f, %.6f) [steiner]\n",
                steinerPoints[i].coords[0], steinerPoints[i].coords[1], steinerPoints[i].coords[2]);
    }
    fprintf(file, "};\n\n");

    //popisu 
    fprintf(file, "%% Popisy Steinerových bodov\n");
    for (int i = 0; i < numSteinerPoints; i++) {
        fprintf(file, "\\node[anchor=south west, red] at (axis cs:%.6f, %.6f, %.6f) {S%d};\n",
                steinerPoints[i].coords[0], steinerPoints[i].coords[1], steinerPoints[i].coords[2], i + 1);
    }
    fprintf(file, "\n");

    // hrany
    fprintf(file, "%% Hrany\n");
    for (int i = 0; i < numSteinerPoints; i++) {
        for (int j = 0; j < steinerPoints[i].neighbors.count; j++) {
            Point *neighbor = steinerPoints[i].neighbors.neighbors[j];
            fprintf(file, "\\addplot3[thick] coordinates {");
            fprintf(file, "(%.6f, %.6f, %.6f) ",
                    steinerPoints[i].coords[0], steinerPoints[i].coords[1], steinerPoints[i].coords[2]);
            fprintf(file, "(%.6f, %.6f, %.6f)};\n",
                    neighbor->coords[0], neighbor->coords[1], neighbor->coords[2]);
        }
    }

    fclose(file);
    printf("'outputTikZ3D.txt' bolo vytvorené\n");
}
