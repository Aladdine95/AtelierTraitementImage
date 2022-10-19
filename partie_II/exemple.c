#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"

/****************************/
/*      CONSTANTS           */
/****************************/

#define SEUIL_C 200
#define BLANC 255
#define NOIR 0
#define IMAGE_GREYSCALE 0
#define IMAGE_RGB 1
#define LAMBDA 0.1
#define USING_HARRIS 1


const int reponse_impulsionnelle[3][3] = {{1,1,1}, {1,1,1}, {1,1,1}};
const int lower_mask[3][3] = {{0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1}};
const int Sobel_horizontal[3][3] = {{-1,0,1} , {-2,0,2}, {-1,0,1}}; //Gradient horizontal Y
const int Sobel_vertical[3][3] = {{-1,-2,-1} , {0,0,0}, {1,2,1}}; //Gradient vertical X
const int filtre_gaussien[3][3] = {
    {0.2,0.1,0.2}, 
    {0.7,1,0.7}, 
    {0.5,0.1,0.5}
    };
const int gradient_mask_detector[3][3] = {{1, 1, 1}, {1, 0, 1}, {1, 1, 1}};

/****************************/
/*       PROTOTYPES         */
/****************************/

void ApplicationFiltres(byte **I, byte** Ix, byte** Iy, int i, int j);
long detectionHarris(byte **Ix, byte** Iy, int i, int j);
void convert_rbg_to_greyscale(rgb8** image_rgb,byte** image_ndg,long nrl, long nrh, long ncl, long nch, double tauxRGB[3]);
byte convolutionGradientMask(byte **I, int i, int j);
long gradient_detector(byte **Ix, byte **Iy, int i, int j);
void listePointsInteret(byte** image, byte** Ix, byte** Iy, byte **In, int* cptWhite, long nrl,long nrh,long ncl,long nch, int isHarris);
void vectorDisplacementEstimation(byte** In1, byte** In2, int** accumulationMat, int shift, long nrl, long nrh, long ncl, long nch);
void electedValueAccumulation(int** accumulationMat, int accuSize, int* vx, int* vy);

int main(int argc, char** argv)
{
    byte **image_ndg; // Image en niveaux de gris
    byte **image_ndg2;
    byte **Ix;
	byte **Iy;
    byte **In;
    byte **In2;

    long C;
    long nrh,nrl,nch,ncl;

    image_ndg=LoadPGM_bmatrix("rice.pgm",&nrl,&nrh,&ncl,&nch);
    image_ndg2=LoadPGM_bmatrix("rice2.pgm",&nrl,&nrh,&ncl,&nch);

    int cptWhiteIn;
    int cptWhiteIn2;

    cptWhiteIn = 0;
    cptWhiteIn2 = 0;

    Ix=bmatrix(nrl,nrh,ncl,nch);
    Iy=bmatrix(nrl,nrh,ncl,nch);
    In=bmatrix(nrl,nrh,ncl,nch);
    In2=bmatrix(nrl, nrh, ncl, nch);

    listePointsInteret(image_ndg, Ix, Iy, In, &cptWhiteIn, nrl, nrh, ncl, nch, USING_HARRIS);
    listePointsInteret(image_ndg2, Ix, Iy, In2, &cptWhiteIn2, nrl, nrh, ncl, nch, USING_HARRIS);

    int cptWhite;

    if(cptWhiteIn <= cptWhiteIn2){
        cptWhite = cptWhiteIn2;
    }
    else{
        cptWhite = cptWhiteIn;
    }

    int intervalAcc = 100; // Interval of accumulationMatrice is : [-50, 50]
    int shift = 50; // We are going to shift negative values from [-50, 50] to [0, intervalAcc], so shift is 50
    int** accumulationMat;
    int vx;
    int vy;

    vx = 0;
    vy = 0;

    accumulationMat = malloc(intervalAcc * sizeof(*accumulationMat));
    for(int i = 0; i < intervalAcc; i++){
        accumulationMat[i] = malloc(intervalAcc * sizeof(accumulationMat[0]));
    }

    for(int i = 0; i < intervalAcc; i++){
        for(int j = 0; j < intervalAcc; j++){
            accumulationMat[i][j] = 0;
        }
    }

    vectorDisplacementEstimation(In, In2, accumulationMat, shift, nrl, nrh, ncl, nch);

    electedValueAccumulation(accumulationMat, intervalAcc, &vx, &vy);
    // We unshift values
    vx = vx - shift;
    vy = vy - shift;
    printf("Vector components are: V(%d, %d)\n", vx, vy);

    SavePGM_bmatrix(In,nrl,nrh,ncl,nch,"res_rice.pgm");
    SavePGM_bmatrix(In2,nrl,nrh,ncl,nch,"res_rice2.pgm");

    free_bmatrix(image_ndg,nrl,nrh,ncl,nch);
    free_bmatrix(image_ndg2,nrl,nrh,ncl,nch);
    free_bmatrix(Ix,nrl,nrh,ncl,nch);
    free_bmatrix(Iy,nrl,nrh,ncl,nch);    
    free_bmatrix(In,nrl,nrh,ncl,nch);
    free_bmatrix(In2,nrl,nrh,ncl,nch);

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
            + filtre_gaussien[2][0] * pow(Ix[i + 1][j - 1],2) + filtre_gaussien[2][1] * pow(Ix[i + 1][j],2) + filtre_gaussien[2][2] * pow(Ix[i + 1][j + 1],2)) / (9 * 255));

    Iytmp = floor((filtre_gaussien[0][0] * pow(Iy[i - 1][j - 1],2) + filtre_gaussien[0][1] * pow(Iy[i - 1][j],2) + filtre_gaussien[0][2] * pow(Iy[i - 1][j + 1],2)
            + filtre_gaussien[1][0] * pow(Iy[i][j - 1],2) + filtre_gaussien[1][1] * pow(Iy[i][j],2) + filtre_gaussien[1][2] * pow(Iy[i][j + 1],2)
            + filtre_gaussien[2][0] * pow(Iy[i + 1][j - 1],2) + filtre_gaussien[2][1] * pow(Iy[i + 1][j],2) + filtre_gaussien[2][2] * pow(Iy[i + 1][j + 1],2)) / (9 * 255));

  
    Ixytmp = floor((filtre_gaussien[0][0] * (Ix[i - 1][j - 1]*Iy[i - 1][j - 1]) + filtre_gaussien[0][1] * (Ix[i - 1][j] * Iy[i - 1][j]) + filtre_gaussien[0][2] * (Ix[i - 1][j + 1] * Iy[i - 1][j + 1])
            + filtre_gaussien[1][0] * (Ix[i][j - 1] * Iy[i][j - 1]) + filtre_gaussien[1][1] * (Ix[i][j] * Iy[i][j]) + filtre_gaussien[1][2] * (Ix[i][j + 1] * Iy[i][j + 1])
            + filtre_gaussien[2][0] * (Ix[i + 1][j - 1] * Iy[i + 1][j - 1]) + filtre_gaussien[2][1] * (Ix[i + 1][j] * Iy[i + 1][j]) + filtre_gaussien[2][2] * (Ix[i + 1][j + 1] * Iy[i + 1][j + 1])) / (9 * 255));
    
    return (Ixtmp * Iytmp- Ixytmp) - (LAMBDA * pow((Ixtmp + Iytmp),2));
}

byte convolutionGradientMask(byte **I, int i, int j)
{
    return floor((gradient_mask_detector[0][0] * I[i - 1][j - 1] + gradient_mask_detector[0][1] * I[i - 1][j] + gradient_mask_detector[0][2] * I[i - 1][j + 1]
            + gradient_mask_detector[1][0] * I[i][j - 1] + gradient_mask_detector[1][1] * I[i][j] + gradient_mask_detector[1][2] * I[i][j + 1]
            + gradient_mask_detector[2][0] * I[i + 1][j - 1] + gradient_mask_detector[2][1] * I[i + 1][j] + gradient_mask_detector[2][2] * I[i + 1][j + 1]) / 9);
}

long gradient_detector(byte **Ix, byte **Iy, int i, int j){
    long IxConvolued;
    long IyConvolued;
    long IxIy;
    long IxIyConvolued;
    
    IxConvolued = convolutionGradientMask(Ix, i, j);
    IyConvolued = convolutionGradientMask(Iy, i, j);
    IxIy = floor(((Ix[i - 1][j - 1]*Iy[i - 1][j - 1]) + (Ix[i - 1][j] * Iy[i - 1][j]) + (Ix[i - 1][j + 1] * Iy[i - 1][j + 1])
            + (Ix[i][j - 1] * Iy[i][j - 1]) + (Ix[i][j] * Iy[i][j]) + (Ix[i][j + 1] * Iy[i][j + 1])
            + (Ix[i + 1][j - 1] * Iy[i + 1][j - 1]) + (Ix[i + 1][j] * Iy[i + 1][j]) + (Ix[i + 1][j + 1] * Iy[i + 1][j + 1])) / (9 * 255));
    IxIyConvolued = floor((gradient_mask_detector[0][0] * (Ix[i - 1][j - 1]*Iy[i - 1][j - 1]) + gradient_mask_detector[0][1] * (Ix[i - 1][j] * Iy[i - 1][j]) + gradient_mask_detector[0][2] * (Ix[i - 1][j + 1] * Iy[i - 1][j + 1])
            + gradient_mask_detector[1][0] * (Ix[i][j - 1] * Iy[i][j - 1]) + gradient_mask_detector[1][1] * (Ix[i][j] * Iy[i][j]) + gradient_mask_detector[1][2] * (Ix[i][j + 1] * Iy[i][j + 1])
            + gradient_mask_detector[2][0] * (Ix[i + 1][j - 1] * Iy[i + 1][j - 1]) + gradient_mask_detector[2][1] * (Ix[i + 1][j] * Iy[i + 1][j]) + gradient_mask_detector[2][2] * (Ix[i + 1][j + 1] * Iy[i + 1][j + 1])) / (9 * 255));

    long result = ((pow(Ix[i][j], 2) * pow(IyConvolued, 2)) + (pow(Iy[i][j], 2) * pow(IxConvolued, 2)) - ((2*IxIy) * IxIyConvolued)) / (pow(IxConvolued, 2) + pow(IyConvolued, 2));
    return result;
} 

void listePointsInteret(byte** image,byte** Ix,byte** Iy, byte **In, int* cptWhite, long nrl,long nrh,long ncl,long nch, int isHarris){

    int C;
    for (int i = nrl + 1; i < nrh - 1; i++) 
    {
        for (int j = ncl + 1; j < nch - 1; j++) 
        {
            ApplicationFiltres(image,Ix,Iy,i,j);
            if(isHarris){
                C = detectionHarris(Ix,Iy,i,j);
            }
            else{
                C = gradient_detector(Ix, Iy, i, j);
            }
            /* if(C != 0)
                printf("x : %d y : %d | C : %d | ABS : %d\n", i,j,C, abs(C)); */
            
            if(abs(C) > SEUIL_C){
                In[i][j] = BLANC;
                *cptWhite = *cptWhite + 1;
            }
            else{
                In[i][j] = NOIR;
            }   
        }
    }
}

void vectorDisplacementEstimation(byte** In1, byte** In2, int** accumulationMat, int shift, long nrl, long nrh, long ncl, long nch){

    for (int i = nrl + 1; i < nrh - 1; i++) 
    {
        for (int j = ncl + 1; j < nch - 1; j++) 
        {
            int condWhite = (In1[i][j] == BLANC);
            if(condWhite){
                if((In1[i][j] == In2[i][j])){
                    // Incrementing to the accumulation matrice for nil vector
                    accumulationMat[0][0] = accumulationMat[0][0] + 1;
                }
                else{
                    int found = 0;
                    int iterLoop = 1;
                    while(found != 1 && iterLoop <= 50){
                        for(int k = i - iterLoop; k < iterLoop; k++){
                            for(int l = j - iterLoop; l < iterLoop; l++){
                                if((k > nrl + 1 && k < nrh - 1) && (l > ncl + 1 && l < nch - 1)){
                                    if(In2[k][l] == BLANC){
                                        int vx = (k - i);
                                        int vy = (l - j);
                                        if((vx > -shift && vx < shift) && (vy > -shift && vy < shift)){
                                            //We shift the values
                                            vx = vx + shift;
                                            vy = vy + shift;
                                            accumulationMat[vx][vy] = accumulationMat[vx][vy] + 1;
                                            found = 1;
                                            goto out;
                                        }
                                    }
                                }
                            }
                        }
                        iterLoop++;
                    }
                }                    
            }
            out:
                //printf("went to\n");
        }
    }
}

void electedValueAccumulation(int** accumulationMat, int accuSize, int* vx, int* vy){
    int maxValue = 0;
    int maxX = 0;
    int maxY = 0;
    for(int i = 0; i < accuSize; i++){
        for(int j = 0; j < accuSize; j++){
            if(accumulationMat[i][j] != 0){
               printf("accumulationMat[%d][%d] : %ld\n",i, j, accumulationMat[i][j]);
            }
            if(maxValue < accumulationMat[i][j]){
                maxValue = accumulationMat[i][j];
                maxX = i;
                maxY = j;
            }
        }
    }
    vx = maxX;
    vy = maxY;
    printf("Max value is (%d, %d) : %d\n", vx, vy, maxValue);
}