#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

/*
void generatePruferSequencies(int pruferSequencies[], int position, int N, int m) {
    if (position == m) {
        printf("Prüferov kód: ");
        for (int i = 0; i < m; i++)
            printf("%d ", pruferSequencies[i]);
        printf("\n");
        return;
    }
    for (int i = 1; i <= N; i++) {
        pruferSequencies[position] = i;
        generatePruferSequencies(pruferSequencies, position + 1, N, m);
    }
}
*/
//predchodca


void generatePruferSequencies(int** pruferSequencies, int* sequence, int position, int N, int sequenciesLenght, int* indexPtr) {
    if (position == sequenciesLenght) {
        for (int i = 0; i < sequenciesLenght; i++)
            pruferSequencies[*indexPtr][i] = sequence[i];
        (*indexPtr)++;
        return;
    }
    for (int i = 1; i <= N; i++) {
        sequence[position] = i;
        generatePruferSequencies(pruferSequencies, sequence, position + 1, N, sequenciesLenght, indexPtr);
    }
}

void printPruferSequencies(int** pruferSequencies, int numberOfTrees, int sequenciesLenght) {
    for (int i = 0; i < numberOfTrees; i++) {
        printf("Sekvencia %d: ", i + 1);
        for (int j = 0; j < sequenciesLenght; j++) {
            printf("%d ", pruferSequencies[i][j]);
        }
        printf("\n");
    }
}

//usporiadané sekvencie dĺžky k z množiny S, pričom každý prvok zo S sa môže použiť najviac r krát
//bruteforce
void generateLimitedPrufer(int** pruferSequencies, int* sequence, int* counts, int position, int N, int sequenciesLenght, int* indexPtr, int limit) {
    if (position == sequenciesLenght) {
        if (pruferSequencies != NULL) { // pre prve volanie sa nezapisuju seq iba sa zistuje pocet
            for (int i = 0; i < sequenciesLenght; i++) {
                pruferSequencies[*indexPtr][i] = sequence[i];
            }
        }
        (*indexPtr)++;
        return;
    }

    for (int i = 1; i <= N; i++) {
        if (counts[i] < limit) {
            sequence[position] = i;
            counts[i]++;  // aby sa cislo pouzilo len urcity pocet krat
            generateLimitedPrufer(pruferSequencies, sequence, counts, position + 1, N, sequenciesLenght, indexPtr, limit);
            counts[i]--; 
        }
    }
}

void freeAllMemory(int** pruferSequencies, int numberOfTrees, int* sequence) {
    for (int i = 0; i < numberOfTrees; i++) {
        free(pruferSequencies[i]);
    }
    free(pruferSequencies);
    free(sequence);
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

    /*printf("\nDEBUG — Stupne vrcholov:\n");
    for (int i = 1; i <= n; i++) {
        printf("Vrchol %d: stupeň %d\n", i, degreeCheck[i]);
        if (degreeCheck[i] > n - 1) {
            printf("CHYBA: Vrchol %d má príliš veľký stupeň!\n", i);
        }
    }*/
}


void treeToTikZ(int edges[][2], int N, int points[][2], FILE* file) {
    fprintf(file, "\\begin{tikzpicture}[scale=1.5]\n");
    for (int i = 0; i < N; i++) {
        int label = i + 1;
        fprintf(file,
            "  \\node[circle,fill=black,inner sep=1pt,label=above:{\\textcolor{red}{%d}}] (%d) at (%d,%d) {};\n",
            label, label, points[i][0], points[i][1]);
    }

    for (int i = 0; i < N - 1; i++) {
        int u = edges[i][0];
        int v = edges[i][1];
        if (u >= 1 && u <= N && v >= 1 && v <= N) {
            fprintf(file, "  \\draw (%d) -- (%d);\n", u, v);
        } else {
            fprintf(file, "  %% Neplatná hrana: (%d -- %d)\n", u, v);
        }
    }

    fprintf(file, "\\end{tikzpicture}\n\n");
}


void printSequenciesInTikZ(int** pruferSequencies, int numberOfTrees, int sequenciesLenght, int points[][2], int N, const char* subor, int toPrint ) {
    FILE* file = fopen(subor, "w");
    if (file == NULL) {
        perror("Chyba pri otváraní súboru");
        return;
    }

    fprintf(file, "\\documentclass{article}\n");
    fprintf(file, "\\usepackage{tikz}\n");
    fprintf(file, "\\begin{document}\n\n");

    int edges[sequenciesLenght + 1][2];

    if (toPrint > 0 && toPrint < numberOfTrees) {
        numberOfTrees = toPrint;
    }
    for (int i = 0; i < numberOfTrees; i++) {
        fprintf(file, "%% Prüferova sekvencia %d\n", i + 1);
        fprintf(file, "\\textbf{Strom č. %d}\\\\\n", i + 1);

        //printf("\n strom %d", i + 1);
        pruferToTree(pruferSequencies[i], sequenciesLenght, edges);

        treeToTikZ(edges, N, points, file);
        fprintf(file, "\\vspace{2cm}\n\n");
    }

    fprintf(file, "\\end{document}\n");
    fclose(file);
    printf("Vykreslenie stromov do súboru %s bolo úspešné.\n", subor);
}

int main() {

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
    int sequenciesLenght = N - 2; //dlzka jednej pruferovej sekvencie;

    //VSEOBECNE GENEROVANIE STROMOV NA N VRCHOLOCH
    /*
    long long numberOfTrees = pow((double)N, (double)(N - 2)); //celkovy pocet stromov (Cayley's formula)

    int** pruferSequencies = malloc(numberOfTrees * sizeof(int*));
    for (int i = 0; i < numberOfTrees; i++) {
        pruferSequencies[i] = malloc(sequenciesLenght * sizeof(int)); 
    }

    int* sequence = malloc(sequenciesLenght * sizeof(int)); // pomocna premenna (jedna sekvencia na ktorej budememe generovat a potom "zapisovat" do pruferSequencies)
    int index = 0; // pocitadlo na ktorom sme riadku v matici array v pruferSequencies

    generatePruferSequencies(pruferSequencies, sequence, 0, N, sequenciesLenght, &index);
 
    printPruferSequencies(pruferSequencies, numberOfTrees, sequenciesLenght);

    printSequenciesInTikZ(pruferSequencies, numberOfTrees, sequenciesLenght, pointsCordinates, N, "vykresliStromy.tex", NULL);

    freeAllMemory(pruferSequencies, numberOfTrees, sequence);
    return 0;
    */

    int limit = 2;
    int* sequence = malloc(sequenciesLenght * sizeof(int));
    int* counts = calloc(N + 1, sizeof(int));  // od 1 po N
    int pocet = 0;
    generateLimitedPrufer(NULL, sequence, counts, 0, N, sequenciesLenght, &pocet, limit);

    int** pruferSequencies = malloc(pocet * sizeof(int*));
    for (int i = 0; i < pocet; i++) {
        pruferSequencies[i] = malloc(sequenciesLenght * sizeof(int));
    }

    for (int i = 0; i <= N; i++) counts[i] = 0;
    int index = 0;
    generateLimitedPrufer(pruferSequencies, sequence, counts, 0, N, sequenciesLenght, &index, limit);

    //printPruferSequencies(pruferSequencies, pocet, sequenciesLenght);
    printf("Pocet stromov: %d\n", pocet);

    printSequenciesInTikZ(pruferSequencies, pocet, sequenciesLenght, pointsCordinates, N, "vykresliStromy.tex", 30);

    freeAllMemory(pruferSequencies, pocet, sequence);
    return 0;

}
