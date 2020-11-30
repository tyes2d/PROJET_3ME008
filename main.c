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


//Question 1
//permet de generer la matrice B
void Creation_Bb(float B[NX][NX], float b[NX], float h, float dt) {
    b[0]=dt/(h*h);
    for (int i=1; i<NX; i++) b[i]=0;
    for (int i=0; i<NX; i++){           //on initialise B, on lui met une valeur par défaut de 0.0
        for (int j=0; j<NX; j++){
            B[i][j]=0.0;
        }
    }
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            if (i==j) B[i][j] = 1-2*dt/(h*h);
            else if (j==i+1 || j==i-1) B[i][j] = dt/(h*h);
            else if (j==NX-2 && i==NX-1) B[i][j] = 2*dt/(h*h); //on cré notre matrice B avec les valeurs quil faut
            else B[i][j] = 0;
        }
    }
    printf("mat B : \n\n");     //affichage de B
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            printf(" %f ", B[i][j]);
        }
        printf("\n");
    }
    printf("vect b : \n\n");    //affichage de b
    for (int i=0; i<NX; i++){
        printf(" %f \n", b[i]);
    }
        
}
//Question 2/3
//Permet de calculer numériquement le profil de temperature dans la plaque au bout d'un temps donné Tstop
void solution_numerique(float B[NX][NX], float b[NX], float h, float dt, float Tinitiale, float Tfinale, int Tstop, float T[NX]){
    float res=0.0;
    for (int i=1; i<NX-1; i++) T[i]=Tinitiale;
    T[0]=Tfinale;
    T[NX-1]=Tfinale;
    for (int n=0; n*dt<Tstop; n++){ // n=Tstop/dt : nombre d'échantillon temporels sur lesquels il faut trouver T jusqu'au temps Tstop
        for (int i=1; i<NX-1;i++){  //On applique ici la formule T^n+1 = [B]*T^n + [b]
            res=0.0;
            for (int j=0; j<NX;j++) res += B[i][j]*T[j];
            T[i] = res+b[i];
        }
    }
}

int main(){
    float h = 30/nx;      //h=L/nx=20/40
    float dt = 30/(nx*nx);
    float dtprime=30/(5*nx*nx);
    float B[NX][NX], Bprime[NX][NX];
    float b[NX], bprime[NX];
    float Ti=750;
    float Tf=25;
    int Ts;
    float T[NX];
    for (int i=0; i<NX; i++) T[i]=0;
    
    printf("Au bout de combien de temps la plaque doit-elle etre retirée de l'eau (en seconde)? \n");
    scanf("%d", &Ts);
    Creation_Bb(B, b, h, dt);
    Creation_Bb(Bprime, bprime, h, dtprime);
    solution_numerique(B, b, h, dt, Ti, Tf, Ts, T);
    printf("vect T : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %f \n", T[i]);
    }
    solution_numerique(Bprime, bprime, h, dtprime, Ti, Tf, Ts, T);
    printf("vect Tprime : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %f \n", T[i]);
    }
}
