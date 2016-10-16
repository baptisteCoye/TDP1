#include "util.h"
#include "cblas.h"

int allouer_vecteur(double ** vecteur, int size){
  *vecteur = malloc(sizeof(double) * size);
  return 1;
}

void desallouer_vecteur(double * vector){
  free(vector);
}

int allouer_matrice(double ** matrice, int m, int n){
  *matrice =  malloc(sizeof(double) * m * n);
  return m;
}

void desallouer_matrice(double * matrice){
  free(matrice);
}


void affiche(int m, int n, double * a, int lda, FILE* flux){
  int i, j;
  for (i = 0; i < m; ++i){
    for (j = 0; j < n; ++j){
      fprintf(flux, "%10.2f  ", a[j*lda+i]);
    }
    fprintf(flux, "\n\n");
  }
}

