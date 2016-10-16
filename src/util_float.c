#ifndef UTIL_FLOAT_H
#define UTIL_FLOAT_H

#include <stdio.h>
#include <stdlib.h>

int allouer_vecteur(float ** vecteur, int size){
  *vecteur = malloc(sizeof(float) * size);
  return 1;
}

void desallouer_vecteur(float ** vecteur){
  free(vector);
}

int allouer_matrice(float ** matrice, int m, int n){
  *matrice = malloc(sizeof(float) * m * n);
  return m;
}

void desallouer_matrice(float * matrice) {
  free(matrice);
}

void affiche(int m, int n, float * a, int lda, FILE* flux){
  int i, j;
  for (i = 0; i < m; ++i){
    for (j = 0; j < n; ++j){
      fprintf(flux, "%f", a[j*lda+i]);
    }
    fprintf(flux, "\n\n");
  }
}

#endif
