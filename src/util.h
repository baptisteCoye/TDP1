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


/*!
 * \brief alloue l'espace necessaire a la création d'un vecteur de taille size 
 * \details stockee par colonnes (Order == CBlasColMajor), retourne la Leading Dimension
 *
 * \param vecteur adresse du vecteur a creer
 * \param size taille du vecteur désiré
 */
int allouer_vecteur(double ** vecteur, int size);


/*!
 * \brief désalloue l'espace réservé a un vecteur déjà crée 
 * \param vecteur adresse du vecteur désallouer
 */
void desallouer_vecteur(double * vecteur);


/*!
 * \brief alloue l'espace necessaire a la création d'une matrice de taille (m,n)(m lignes, n colonnes) 
 * \details stockee par colonnes (Order == CBlasColMajor), retourn la Leading Dimension
 *
 * \param vecteur adresse de la matrice a creer
 * \param m nombre de lignes de la matrice
 * \param n nombre de colonnes de la matrice
 */
int allouer_matrice(double ** matrice, int m, int n);


/*!
 * \brief désalloue l'espace réservé a une matrice déjà crée 
 *
 * \param matrice adresse de la matrice à désallouer
 */
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
