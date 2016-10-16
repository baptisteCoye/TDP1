#ifndef UTIL_FLOAT_H
#define UTIL_FLOAT_H

#include <stdio.h>
#include <stdlib.h>

int allouer_vecteur(float ** vecteur, int size);

void desallouer_vecteur(float ** vecteur);

int allouer_matrice(float ** matrice, int m, int n);

void desallouer_matrice(float * matrice);

void affiche(int m, int n, float * a, int lda, FILE* flux);

#endif
