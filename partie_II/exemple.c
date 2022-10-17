#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

#define SEUIL 5
#define SEUIL_C 10
#define BLANC 255
#define NOIR 0
#define IMAGE_GREYSCALE 0
#define IMAGE_RGB 1
#define LAMBDA 0.5

const int reponse_impulsionnelle[3][3] = {{1,1,1}, {1,1,1}, {1,1,1}};
const int Sobel_horizontal[3][3] = {{-1,0,1} , {-2,0,2}, {-1,0,1}}; //Gradient horizontal Y
const int Sobel_vertical[3][3] = {{-1,-2,-1} , {0,0,0}, {1,2,1}}; //Gradient vertical X
const int filtre_gaussien[3][3] = {{0.3,0.5,0.3} , {0.5,1,0.5}, {0.3,0.5,0.3}};

void ApplicationFiltres(byte **I, byte** Ix, byte** Iy, int i, int j);
long detectionHarris(byte **Ix, byte** Iy, int i, int j);
void convert_rbg_to_greyscale(rgb8** image_rgb,byte** image_ndg,long nrl, long nrh, long ncl, long nch, double tauxRGB[3]);

int main(int argc, char** argv)
{
    byte **image_ndg; // Image en niveaux de gris
    byte **Ix;
	byte **Iy;
    byte **In;
    byte **Icoin;
    //byte **Harris;

    byte norme;
    long C;
    long nrh,nrl,nch,ncl;

    image_ndg=LoadPGM_bmatrix("route0.pgm",&nrl,&nrh,&ncl,&nch);

    Ix=bmatrix(nrl,nrh,ncl,nch);
    Iy=bmatrix(nrl,nrh,ncl,nch);
    In=bmatrix(nrl,nrh,ncl,nch);
    Icoin=bmatrix(nrl,nrh,ncl,nch);


    for (int i = nrl + 1; i < nrh - 1; i++) 
    {
        for (int j = ncl + 1; j < nch - 1; j++) 
        {
            ApplicationFiltres(image_ndg,Ix,Iy,i,j);
            C = detectionHarris(Ix,Iy, i, j);
            if( C > 1)
                printf("C : %ld \n", C);

            // norme = sqrt(pow(Ix[i][j],2) + pow(Iy[i][j],2));
            // if(norme > SEUIL) In[i][j] = 255;
            // else In[i][j] = 0;
           
        }
    }
    // <Ix^2><Iy^2> - <IxIy> - lambda (<Ix^2> + <Iy^2>)^2


    //SavePGM_bmatrix(In,nrl,nrh,ncl,nch,"res.pgm");

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

long detectionHarris(byte **Ix, byte** Iy, int i, int j)
{
    byte Ixtmp;
    byte Iytmp;
    byte Ixytmp;

    Ixtmp = floor((filtre_gaussien[0][0] * pow(Ix[i - 1][j - 1],2) + filtre_gaussien[0][1] * pow(Ix[i - 1][j],2) + filtre_gaussien[0][2] * pow(Ix[i - 1][j + 1],2)
            + filtre_gaussien[1][0] * pow(Ix[i][j - 1],2) + filtre_gaussien[1][1] * pow(Ix[i][j],2) + filtre_gaussien[1][2] * pow(Ix[i][j + 1],2)
            + filtre_gaussien[2][0] * pow(Ix[i + 1][j - 1],2) + filtre_gaussien[2][1] * pow(Ix[i + 1][j],2) + filtre_gaussien[2][2] * pow(Ix[i + 1][j + 1],2)) / 9);

    Iytmp = floor((filtre_gaussien[0][0] * pow(Iy[i - 1][j - 1],2) + filtre_gaussien[0][1] * pow(Iy[i - 1][j],2) + filtre_gaussien[0][2] * pow(Iy[i - 1][j + 1],2)
            + filtre_gaussien[1][0] * pow(Iy[i][j - 1],2) + filtre_gaussien[1][1] * pow(Iy[i][j],2) + filtre_gaussien[1][2] * pow(Iy[i][j + 1],2)
            + filtre_gaussien[2][0] * pow(Iy[i + 1][j - 1],2) + filtre_gaussien[2][1] * pow(Iy[i + 1][j],2) + filtre_gaussien[2][2] * pow(Iy[i + 1][j + 1],2)) / 9);

  
    Ixytmp = floor((filtre_gaussien[0][0] * (Ix[i - 1][j - 1]*Iy[i - 1][j - 1]) + filtre_gaussien[0][1] * (Ix[i - 1][j] * Iy[i - 1][j]) + filtre_gaussien[0][2] * (Ix[i - 1][j + 1] * Iy[i - 1][j + 1])
            + filtre_gaussien[1][0] * (Ix[i][j - 1] * Iy[i][j - 1]) + filtre_gaussien[1][1] * (Ix[i][j] * Iy[i][j]) + filtre_gaussien[1][2] * (Ix[i][j + 1] * Iy[i][j + 1])
            + filtre_gaussien[2][0] * (Ix[i + 1][j - 1] * Iy[i + 1][j - 1]) + filtre_gaussien[2][1] * (Ix[i + 1][j] * Iy[i + 1][j]) + filtre_gaussien[2][2] * (Ix[i + 1][j + 1] * Iy[i + 1][j + 1])) / 9);
    
    // if(Ixtmp > 0)
    //     printf("Ixtmp : %hhu\n", Ixtmp);
    // if(Iytmp > 0)
    //     printf("Iytmp : %hhu\n", Iytmp);
    // if(Ixytmp > 0)
    //     printf("Ixytmp : %hhu\n", Ixytmp);

    printf("Ixtmp * Iytmp : %hhu \n", Ixtmp);
    printf("Ixytmp : %hhu \n", Ixtmp);
    printf("Ixtmp * Iytmp : %hhu \n", Ixtmp);

    if(Iytmp > 0)
        printf("Iytmp : %hhu\n", Iytmp);
    if(Ixytmp > 0)
        printf("Ixytmp : %hhu\n", Ixytmp);

    if ((Ixtmp * Iytmp) - Ixytmp - LAMBDA * pow((Ixtmp + Iytmp),2) > 0 )
        printf("YES : %hhu", ((Ixtmp * Iytmp) - Ixytmp - LAMBDA * pow((Ixtmp + Iytmp),2)) );

    return (Ixtmp * Iytmp - Ixytmp - LAMBDA * pow((Ixtmp + Iytmp),2) );
}

    // <Ix^2><Iy^2> - <IxIy> - lambda (<Ix^2> + <Iy^2>)^2
