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
    b[0]=25*a*dt/(dx*dx);               //ici T_0=25°C et T_L=25°C
    b[NX-1]=25*a*dt/(dx*dx);
    for (int i=1; i<NX-1; i++) b[i]=0;
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
void solution_numeriqueEXPLI(double B[NX][NX], double b[NX], double dx, double dt, double Tinitiale, double Tfinale, int Tstop, double T[NX]){
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

//-------------PARTIE 2------------------

//Question 6
void Creation_Cc(double C[NX][NX], double c[NX], double dx, double dt, double a) {
    c[0]=25*a*dt/(dx*dx);
    c[NX-1]=25*a*dt/(dx*dx);
       for (int i=1; i<NX-1; i++) c[i]=0;
       for (int i=0; i<NX; i++){           //on initialise B, on lui met une valeur par défaut de 0.0
           for (int j=0; j<NX; j++){
               C[i][j]=0.0;
           }
       }
       for (int i=0; i<NX; i++){
           for (int j=0; j<NX; j++){
               if (i==j) C[i][j] = 1+2*a*dt/(dx*dx);
               else if (j==i+1 || j==i-1) C[i][j] = -a*dt/(dx*dx);
               else C[i][j] = 0;
           }
       }
       /*printf("\nmat C : \n\n");     //affichage de C
       for (int i=0; i<NX; i++){
           for (int j=0; j<NX; j++){
               printf(" %lf ", C[i][j]);
           }
           printf("\n");
       }
       printf("vect c : \n\n");    //affichage de c
       for (int i=0; i<NX; i++){
           printf(" %lf \n", c[i]);
       }*/
}
//Question 7

void facto_LU( double C[NX][NX] , double L[NX][NX] , double U[NX][NX]) {
    int i,j,k;
    float l;
    for (j=0;j<NX;j++) {
        for (i=0;i<=j;i++) {
            l=0;
            for (k=0;k<i;k++) l += L[i][k]*U[k][j];
            U[i][j] = C[i][j] - l;
        }
        L[j][j] = 1;
        for (i=j+1;i<NX;i++) {
            l=0;
            for (k=0;k<=j-1;k++) l += L[i][k]*U[k][j];
            L[i][j] = (C[i][j] - l)/U[j][j];
        }
                
    }
    /*printf("\nmat L : \n\n");     //affichage de C
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            printf(" %lf ", L[i][j]);
        }
        printf("\n");
    }
    printf("\nmat U : \n\n");     //affichage de C
       for (int i=0; i<NX; i++){
           for (int j=0; j<NX; j++){
               printf(" %lf ", U[i][j]);
           }
           printf("\n");
       }*/
}



void resol_trig_inf( double A[NX][NX] , double x[NX], double b[NX]) {
    int i;
    float l;
    
    for (i=0;i<NX;i++) {
        l=0;
        for (int k=0;k<=i-1;k++) l += A[i][k]*x[k];
        x[i] = (b[i] -l)/A[i][i];
    }
}


void resol_trig_sup( double A[NX][NX] , double x[NX], double b[NX]) {
    int i;
    double l;
    
    for (i=NX-1;i>=0;i--) { //ici peut etre pb sur i=NX ou =NX-1
        l=0;
        for (int k=i+1;k<NX;k++) l += A[i][k]*x[k];
        x[i] = (b[i]-l)/A[i][i];
    }
}

void solution_numeriqueIMPLI(double C[NX][NX], double c[NX], double dx, double dt, double Tinitiale, double Tfinale, int Tstop, double T[NX]){
    double vect[NX];
    double y[NX];
    double L[NX][NX], U[NX][NX];
    T[0]=Tfinale;
    T[NX-1]=Tfinale;
    for (int i=1; i<NX-1; i++) T[i]=Tinitiale;
    facto_LU(C, L, U);
    for (int n=0; n*dt<Tstop; n++){ //On applique ici la formule [C]T^n+1 = T^n + [c]
        for (int i=0; i<NX; i++) vect[i]=T[i]+c[i];
        resol_trig_inf(L, y, vect);
        resol_trig_sup(U, T, y);
    }
}


int main(){
    double z=0;
    double a=0.00001;                   //Diffusivité du materiaux
    double dx = 0.03/nx;                //dx=L/nx=20/40
    double dt = 30/(nx*nx);             //Pas de temps
    double dtprime=30/(5*nx*nx);        //Pas de temps plus fin de la question 4
    double dtQ5 ;
    double errT = 1;                    // on choisit 1 ici, un choix qui nous semble cohérent avec l'ordre de grandeur des températures en jeu
    double normT=100;
    double B[NX][NX], Bprime[NX][NX], BQ5[NX][NX], C[NX][NX];
    double b[NX], bprime[NX], bQ5[NX], c[NX];
    double Ti=750, Tf=25;               //Condition initiale et finale
    int ts;
    double T[NX], Tprime[NX], TQ5[NX], TQ[NX], deltaT[NX];
    for (int i=0; i<NX; i++) T[i]=0;
    
    printf("Au bout de combien de temps la plaque doit-elle etre retirée de l'eau (en seconde, uniquement nombre entier)? \n"); //
    scanf("%d", &ts);                                                                                                           // PEUT ETRE LE MERTTRE DIRECT DANS LES FONCTIONS de  RESOLUTION ??
    Creation_Bb(B, b, dx, dt, a);
    Creation_Bb(Bprime, bprime, dx, dtprime, a);
    solution_numeriqueEXPLI(B, b, dx, dt, Ti, Tf, ts, T);
    printf("vect T : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %lf \n", T[i]);
    }
    solution_numeriqueEXPLI(Bprime, bprime, dx, dtprime, Ti, Tf, ts, Tprime);
    printf("vect Tprime : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %lf \n", Tprime[i]);
    }
    
    //Question 5
    
    dtQ5=dtprime;
    for (int i=0; i<NX; i++) TQ[i]=0;
    do {
        Creation_Bb(BQ5, bQ5, dx, dtQ5, a);
        solution_numeriqueEXPLI(BQ5, bQ5, dx, dtQ5, Ti, Tf, ts, TQ5);
        for (int i=0; i<NX; i++) deltaT[i]=fabs(TQ[i]-TQ5[i]);
        for (int i=0; i<NX; i++) TQ[i]=TQ5[i];
        normT = norme_vect(deltaT);
        dtQ5=dtQ5/5;
        z++;
    }
    while (normT>errT);
    printf("dt optimal : %lf\n", dtQ5);
    Creation_Cc(C, c, dx, dtQ5, a);
    
    //Question 7
    solution_numeriqueIMPLI(C, c, dx, dt, Ti, Tf, ts, T);
    printf("vect T : \n\n");    //affichage de T à Ts
    for (int i=0; i<NX; i++){
        printf(" %lf \n", T[i]);
    }
    
}
