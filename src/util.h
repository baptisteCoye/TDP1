#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>

/*!
 * \file util.h
 *
 * \author Mathieu Deschamps
 * \date 7/10/2016
 *
 * \brief Ce fichier contient différentes fonctions utilitaires.
 */



int allouer_vecteur(double ** vecteur, int size);
void desallouer_vecteur(double * vecteur);

int allouer_matrice(double ** matrice, int m, int n);
void desallouer_matrice(double * matrice);

/*!
 * \brief affiche une matrice 
 * \details stockee par colonnes (Order == CBlasColMajor)
 *
 * \param m nombre de lignes
 * \param n nombre de colonnes
 * \param a la matrice a afficher
 * \param lda leading dimension, espace memoire qui separe deux elements consecutifs sur une ligne
 * \param flux flux sur lequel afficher la matrice
 */
void affiche(int m, int n, double * a, int lda, FILE* flux);

#endif /* UTIL_H */
