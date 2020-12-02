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



double norme_vect(double x[NX]) { //renvoie la norme d'un vecteur (tableau à 1 entrée)
    double var=0;
    int i;
    for (i=0;i<NX;i++) var += (x[i]*x[i]);
    return sqrt(var);
}

//Question 1
//permet de generer la matrice B
void Creation_Bb(double B[NX][NX], double b[NX], double h, double dt) {
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
void solution_numerique(double B[NX][NX], double b[NX], double h, double dt, double Tinitiale, double Tfinale, int Tstop, double T[NX]){
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
void rapidite_conv(double B[NX][NX], double b[NX], double h, double dt){
    //nous allons nous aider de la
}


/*double puissance_it( double A[MAX][MAX], double x0[MAX] , double x[MAX],double eps  , int n) {
    
    int i,j,k=0;
    double v[MAX], L1=1,L0=10000000, r, var ;
    
    
    for (i=0;i<n;i++) x[i] = x0[i];
    r = norme2_vect(x,n);
    for (i=0;i<n;i++) x[i] = x[i]/r;
    
    while ((L1 - L0) >= eps || (L1-L0) <= -eps ) {
        
        L0 = L1;

        for (i=0;i<n;i++) {
            var =0;
            for (j=0;j<n;j++) var+= A[i][j]*x[j];
            v[i] = var;
        }
                

        L1=0;
        for (i=0;i<n;i++) L1 +=    x[i]*v[i];
        

                
          for (i=0;i<n;i++) x[i] = v[i];
                
        r = norme2_vect(x,n);
        for (i=0;i<n;i++) x[i] = x[i]/r;
    
        k ++;
    }
            
    printf("%d itérations\n",k);

    return L1;
    
}*/


int main(){
    double a=0;
    double h = 30/nx;      //h=L/nx=20/40
    double dt = 30/(nx*nx);
    double dtprime=30/(5*nx*nx);
    double dtQ5 ;
    double errT = 40; // on choisit  ici car une erreur de  degrés sur un centaine nous semble correct
    double normT=100;
    double B[NX][NX], Bprime[NX][NX], BQ5[NX][NX];
    double b[NX], bprime[NX], bQ5[NX];
    double Ti=750, Tf=25;
    int ts;
    double T[NX], Tprime[NX], TQ5[NX], deltaT[NX];
    for (int i=0; i<NX; i++) T[i]=0;
    
    printf("Au bout de combien de temps la plaque doit-elle etre retirée de l'eau (en seconde)? \n");
    scanf("%d", &ts);
    Creation_Bb(B, b, h, dt);
    Creation_Bb(Bprime, bprime, h, dtprime);
    solution_numerique(B, b, h, dt, Ti, Tf, ts, T);
    printf("vect T : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %lf \n", T[i]);
    }
    solution_numerique(Bprime, bprime, h, dtprime, Ti, Tf, ts, Tprime);
    printf("vect Tprime : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %lf \n", T[i]);
    }
    
    //Question 5
    
    dtQ5=dtprime;
    while (normT>errT){
        Creation_Bb(BQ5, bQ5, h, dtQ5);
        solution_numerique(BQ5, bQ5, h, dtQ5, Ti, Tf, ts, TQ5);
        for (int i=0; i<NX; i++) deltaT[i]=fabs(T[i]-TQ5[i]);
        normT = norme_vect(deltaT);
        dtQ5=dtQ5/5;
        a++;
    }
    printf("dt optimal : %lf\n", dtQ5);
}
