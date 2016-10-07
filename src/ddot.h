#ifndef DDOT_H
#define DDOT_H

#include <stdio.h>
#include <stdlib.h>
/*!
 * \file ddot.h
 * \brief Ce fichier contient la fonction permettant de réaliser un produit scalaire.
 */


/*!
 * \brief définie le produit scalaire entre deux vecteurs de taille N 
 *
 * \param N taille des vecteurs dont le produit scalaire doit être réalisé
 * \param DX premier vecteur
 * \param INCX espacement entre 2 valeurs consécutives du vecteur DX dans la mémoire
 * \param DY deuxieme vecteur
 * \param INCY espacement entre 2 valeurs consécutives du vecteur DY dans la mémoire 
 */
double ddot(int N, double * DX, int INCX, double * DY, int INCY);

#endif
