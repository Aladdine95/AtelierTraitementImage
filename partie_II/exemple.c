#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

#define SEUIL 4
#define BLANC 255
#define NOIR 0
#define IMAGE_GREYSCALE 0
#define IMAGE_RGB 1

const int reponse_impulsionnelle[3][3] = {{1,1,1}, {1,1,1}, {1,1,1}};
const int Sobel_horizontal[3][3] = {{-1,0,1} , {-2,0,2}, {-1,0,1}}; //Gradient horizontal Y
const int Sobel_vertical[3][3] = {{-1,-2,-1} , {0,0,0}, {1,2,1}}; //Gradient vertical X

void ApplicationFiltres(byte **I, byte** Ix, byte** Iy, int i, int j);
void convert_rbg_to_greyscale(rgb8** image_rgb,byte** image_ndg,long nrl, long nrh, long ncl, long nch, double tauxRGB[3]);

int main(int argc, char** argv)
{
    byte **image_ndg; // Image en niveaux de gris
    byte **Ix;
	byte **Iy;
    byte **In;
    byte norme;
    
    long nrh,nrl,nch,ncl;


    image_ndg=LoadPGM_bmatrix("cubes.pgm",&nrl,&nrh,&ncl,&nch);

    Ix=bmatrix(nrl,nrh,ncl,nch);
    Iy=bmatrix(nrl,nrh,ncl,nch);
    In=bmatrix(nrl,nrh,ncl,nch);


    for (int i = nrl + 1; i < nrh - 1; i++) 
    {
        for (int j = ncl + 1; j < nch - 1; j++) 
        {
            ApplicationFiltres(image_ndg,Ix,Iy,i,j);
            norme = sqrt(pow(Ix[i][j],2) + pow(Iy[i][j],2));
            if(norme > SEUIL) In[i][j] = 255;
            else In[i][j] = 0;
        }
    }

    SavePGM_bmatrix(In,nrl,nrh,ncl,nch,"res.pgm");

    free_bmatrix(image_ndg,nrl,nrh,ncl,nch);
    free_bmatrix(Ix,nrl,nrh,ncl,nch);
    free_bmatrix(Iy,nrl,nrh,ncl,nch);    
    free_bmatrix(In,nrl,nrh,ncl,nch);

    return 0;
}


void ApplicationFiltres(byte **I, byte** Ix, byte** Iy, int i, int j)
{
    Ix[i][j] = floor((Sobel_vertical[0][0] * I[i - 1][j - 1] + Sobel_vertical[0][1] * I[i - 1][j] + Sobel_vertical[0][2] * I[i - 1][j + 1]
            + Sobel_vertical[1][0] * I[i][j - 1] + Sobel_vertical[1][1] * I[i][j] + Sobel_vertical[1][2] * I[i][j + 1]
            + Sobel_vertical[2][0] * I[i + 1][j - 1] + Sobel_vertical[2][1] * I[i + 1][j] + Sobel_vertical[2][2] * I[i + 1][j + 1]) / 9);

    Iy[i][j] = floor((Sobel_horizontal[0][0] * I[i - 1][j - 1] + Sobel_horizontal[0][1] * I[i - 1][j] + Sobel_horizontal[0][2] * I[i - 1][j + 1]
            + Sobel_horizontal[1][0] * I[i][j - 1] + Sobel_horizontal[1][1] * I[i][j] + Sobel_horizontal[1][2] * I[i][j + 1]
            + Sobel_horizontal[2][0] * I[i + 1][j - 1] + Sobel_horizontal[2][1] * I[i + 1][j] + Sobel_horizontal[2][2] * I[i + 1][j + 1]) / 9);
}