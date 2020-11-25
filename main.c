//
//  main.c
//  PROJET_3ME008
//
//  Created by Gauthier LECLERCQ on 17/11/2020.
//  Copyright © 2020 Gauthier LECLERCQ. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define NX 40
#define nx 40.0



int main() {
    float h = 30/nx;      //h=L/nx=20/40
    float dt = 30/(nx*nx);
    float B[NX][NX];
    float b[NX];
    b[0]=dt/(h*h);
    for (int i=1; i<NX; i++) b[i]=0;
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            B[i][j]=0.0;
        }
    }
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            if (i==j) B[i][j] = 1-2*dt/(h*h);
            else if (j==i+1 || j==i-1) B[i][j] = dt/(h*h);
            else if (j==NX-2 && i==NX-1) B[i][j] = 2*dt/(h*h);
            else B[i][j] = 0;
        }
    }
    printf("mat B : \n\n");
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            printf(" %f ", B[i][j]);
        }
        printf("\n");
    }
    printf("vect b : \n\n");
    for (int i=0; i<NX; i++){
        printf(" %f \n", b[i]);
    }
        
    return 0;
}
