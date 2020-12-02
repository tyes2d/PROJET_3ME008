//
//  main.c
//  PROJET_3ME008
//
//  Created by Gautdxier LECLERCQ on 17/11/2020.
//  Copyrigdxt © 2020 Gautdxier LECLERCQ. All rigdxts reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define NX 40
#define nx 40.0



double norme_vect(double x[NX]) { //renvoie la norme d'un vecteur (tableau à 1 entrée)
    double var=0;
    for (int i=0;i<NX;i++) var += (x[i]*x[i]);
    return sqrt(var);
}

//Question 1
//permet de generer la matrice B
void Creation_Bb(double B[NX][NX], double b[NX], double dx, double dt, double a) {
    b[0]=a*dt/(dx*dx);
    for (int i=1; i<NX; i++) b[i]=0;
    for (int i=0; i<NX; i++){           //on initialise B, on lui met une valeur par défaut de 0.0
        for (int j=0; j<NX; j++){
            B[i][j]=0.0;
        }
    }
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            if (i==j) B[i][j] = 1-2*a*dt/(dx*dx);
            else if (j==i+1 || j==i-1) B[i][j] = a*dt/(dx*dx);
            else if (j==NX-2 && i==NX-1) B[i][j] = 2*a*dt/(dx*dx); //on cré notre matrice B avec les valeurs quil faut
            else B[i][j] = 0;
        }
    }
    /*printf("\nmat B : \n\n");     //affichage de B
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            printf(" %lf ", B[i][j]);
        }
        printf("\n");
    }
    printf("vect b : \n\n");    //affichage de b
    for (int i=0; i<NX; i++){
        printf(" %lf \n", b[i]);
    }*/
        
}
//Question 2/3
//Permet de calculer numériquement le profil de temperature dans la plaque au bout d'un temps donné Tstop
void solution_numerique(double B[NX][NX], double b[NX], double dx, double dt, double Tinitiale, double Tfinale, int Tstop, double T[NX]){
    double res=0.0;
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

//Question 4 : Calcul de la rapidité de convergence
void rapidite_conv(double B[NX][NX], double b[NX], double dx, double dt){
    //nous allons nous aider de la
}



int main(){
    double z=0;
    double a=0.00001;
    double dx = 0.03/nx;      //dx=L/nx=20/40
    double dt = 30/(nx*nx);
    double dtprime=30/(5*nx*nx);
    double dtQ5 ;
    double errT = 0.1; // on choisit  ici car une erreur de  degrés sur un centaine nous semble correct
    double normT=100;
    double B[NX][NX], Bprime[NX][NX], BQ5[NX][NX];
    double b[NX], bprime[NX], bQ5[NX];
    double Ti=750, Tf=25;
    int ts;
    double T[NX], Tprime[NX], TQ5[NX], deltaT[NX];
    for (int i=0; i<NX; i++) T[i]=0;
    
    printf("Au bout de combien de temps la plaque doit-elle etre retirée de l'eau (en seconde)? \n");
    scanf("%d", &ts);
    Creation_Bb(B, b, dx, dt, a);
    Creation_Bb(Bprime, bprime, dx, dtprime, a);
    solution_numerique(B, b, dx, dt, Ti, Tf, ts, T);
    printf("vect T : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %lf \n", T[i]);
    }
    solution_numerique(Bprime, bprime, dx, dtprime, Ti, Tf, ts, Tprime);
    printf("vect Tprime : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %lf \n", T[i]);
    }
    
    //Question 5
    
    dtQ5=dtprime;
    do {
        Creation_Bb(BQ5, bQ5, dx, dtQ5, a);
        solution_numerique(BQ5, bQ5, dx, dtQ5, Ti, Tf, ts, TQ5);
        for (int i=0; i<NX; i++) deltaT[i]=fabs(T[i]-TQ5[i]);
        normT = norme_vect(deltaT);
        dtQ5=dtQ5/5;
        z++;
    }
    while (normT>errT);
    printf("dt optimal : %lf\n", dtQ5);
}
