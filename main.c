#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define NX 40           //Variable globale entiere pour la definition de tableau et pour les condition d'itération dans divers programmes
#define nx 40.0         //Variable globale flottante pour la definition des pas (temps et espace) qui sont de type 'double'



double norme_vect(double x[NX]) {               //Fonction qui renvoie la norme 2 d'un vecteur (tableau à 1 entrée)
    double var=0;
    for (int i=0;i<NX;i++) var += (x[i]*x[i]);
    return sqrt(var);
}

//Question 1

void Creation_Bb(double B[NX][NX], double b[NX], double dx, double dt, double a) {
    b[0]=25*a*dt/(dx*dx);                             //Ici T_0=25°C et T_L=25°C
    b[NX-1]=25*a*dt/(dx*dx);
    for (int i=1; i<NX-1; i++) b[i]=0;                //On initialise b, on lui met une valeur par défaut de 0.0
    for (int i=0; i<NX; i++){                         //On initialise B, on lui met une valeur par défaut de 0.0
        for (int j=0; j<NX; j++){
            B[i][j]=0.0;
        }
    }
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            if (i==j) B[i][j] = 1-2*a*dt/(dx*dx);
            else if (j==i+1 || j==i-1) B[i][j] = a*dt/(dx*dx);      //On remplit les diagonales au dessus et en dessous de la diagonale principale
            else if (j==NX-2 && i==NX-1) B[i][j] = 2*a*dt/(dx*dx);  //On remplit la diagonale principale de la matrice B
            else B[i][j] = 0;                                       //Le reste de la matrice est rempli de 0
        }
    }
    /*printf("\nmat B : \n\n");                                     //Affichage de B          DECOMMENTER CE BLOC POUR L'AFFICHAGE DE B et b
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            printf(" %lf ", B[i][j]);
        }
        printf("\n");
    }
    printf("vect b : \n\n");                                        //Affichage de b
    for (int i=0; i<NX; i++){
        printf(" %lf \n", b[i]);
    }*/
        
}

//Question 2/3

//Permet de calculer numériquement le profil de temperature dans la plaque au bout d'un temps donné Tstop
void solution_numeriqueEXPLI(double B[NX][NX], double b[NX], double dx, double dt, double Tinitiale, double Tfinale, int Tstop, double T[NX]){
    double res=0.0;
    int n;
    for (int i=0; i<NX; i++) T[i]=Tinitiale;
    for (n=0; n*dt<Tstop; n++){                               // n=Tstop/dt : nombre d'échantillon temporels sur lesquels il faut trouver T jusqu'au temps Tstop
        for (int i=0; i<NX;i++){                              //On applique ici la formule T^n+1 = [B]*T^n + [b]
            res=0.0;
            for (int j=0; j<NX;j++) res += B[i][j]*T[j];
            T[i] = res+b[i];
        }
    }
}
void Evo_temporelle_Question3(double B[NX][NX], double b[NX], double dx, double dt, double Tinitiale, double Tfinale, double T[NX], double EvoTemporelleT[6400]){
    //On utilise finalement quasiment la même fonction que 'solution_numeriqueEXPLI' 
    double res=0.0;
    for (int i=0; i<NX; i++) T[i]=Tinitiale;
    for (int n=0; n*dt<120; n++){                               //On fait varier n de 0 à 6399 pour remplir le tableau EvoTemporelleT[n]
        for (int i=0; i<NX;i++){
            res=0.0;
            for (int j=0; j<NX;j++) res += B[i][j]*T[j];
            T[i] = res+b[i];                                    //A chaque itération (à t=n*dt) on préleve T[19] pour le mettre dans EvoTemporelleT[n]
        }
        EvoTemporelleT[n]=T[19];
    }
}




//-------------PARTIE 2------------------




//Question 6

void Creation_Cc(double C[NX][NX], double c[NX], double dx, double dt, double a) {
    c[0]=25*a*dt/(dx*dx);                             //Ici T_0=25°C et T_L=25°C
    c[NX-1]=25*a*dt/(dx*dx); 
       for (int i=1; i<NX-1; i++) c[i]=0;             //On initialise c, on lui met une valeur par défaut de 0.0
       for (int i=0; i<NX; i++){                      //On initialise C, on lui met une valeur par défaut de 0.0
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
      /* printf("\nmat C : \n\n");                    //Affichage de C               DECOMMENTER CE BLOC POUR L'AFFICHAGE DE C et c
       for (int i=0; i<NX; i++){
           for (int j=0; j<NX; j++){
               printf(" %lf ", C[i][j]);
           }
           printf("\n");
       }
       printf("vect c : \n\n");                       //Affichage de c
       for (int i=0; i<NX; i++) printf(" %lf \n", c[i]);*/
}

//Question 7

void facto_LU( double C[NX][NX] , double L[NX][NX] , double U[NX][NX]) {
  //On a récupéré la fonction de factorisation LU des TP
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
    /*printf("\nmat L : \n\n");                      //affichage de L                DECOMMENTER CE BLOC POUR L'AFFICHAGE DE L ET U POUR LA RESOLUTION
    for (int i=0; i<NX; i++){
        for (int j=0; j<NX; j++){
            printf(" %lf ", L[i][j]);
        }
        printf("\n");
    }
    printf("\nmat U : \n\n");                        //affichage de U
       for (int i=0; i<NX; i++){
           for (int j=0; j<NX; j++){
               printf(" %lf ", U[i][j]);
           }
           printf("\n");
       }*/
}

void resol_trig_inf( double A[NX][NX] , double x[NX], double b[NX]) {
  //On a récupéré la fonction de la resolution triangulaire inférieure des TP
    int i;
    float l;
    
    for (i=0;i<NX;i++) {
        l=0;
        for (int k=0;k<=i-1;k++) l += A[i][k]*x[k];
        x[i] = (b[i] -l)/A[i][i];
    }
}

void resol_trig_sup( double A[NX][NX] , double x[NX], double b[NX]) {
  //On a récupéré la fonction de la resolution triangulaire supérieure des TP
    int i;
    double l;
    
    for (i=NX-1;i>=0;i--) {
        l=0;
        for (int k=i+1;k<NX;k++) l += A[i][k]*x[k];
        x[i] = (b[i]-l)/A[i][i];
    }
}

void solution_numeriqueIMPLI(double C[NX][NX], double c[NX], double dx, double dt, double Tinitiale, double Tfinale, int Tstop, double T[NX]){
    double vect[NX];
    double y[NX];                                             //Vecteur de résolution intermédiaire : LUx=b --> Ly=b --> Ux=y
    double L[NX][NX], U[NX][NX];
    for (int i=0; i<NX; i++) T[i]=Tinitiale;
    facto_LU(C, L, U);
    for (int n=0; n*dt<Tstop; n++){                           //On applique ici la formule [C]T^n+1 = T^n + [c]
        for (int i=0; i<NX; i++) vect[i]=0;
        for (int i=0; i<NX; i++) vect[i]=T[i]+c[i];
        resol_trig_inf(L, y, vect);
        resol_trig_sup(U, T, y);
    }
}



int main(){
    int z=0;                                //Variable quelconque
    double a=0.00001;                       //Diffusivité du materiaux
    double dx = 0.03/nx;                    //dx=L/nx=0.03/40
    double dt = 30/(nx*nx);                 //Pas de temps initial
    double evoT[6400];                      //Vecteur stockant la valeurs de T au milieu de la plaque à chaque pas de temps. Une case par seconde (6400*dt = 120)
    double dtprime=30/(5*nx*nx);            //Pas de temps pour la fin de la question 4
    double dtQ5, dtQ8 ;                     //Pas de temps pour la question 5 et la question 8
    double errT = 10;                       //On choisit 10 ici, une marge erreur qui nous semble correcte. Ca revient à une erreur moyenne sur x de 1.6°C
    double normT=100;                       //Variable stockant la norme du vecteur T (Température) à la question 5. On l'initialise suffisamment grand pour pouvoir lancer la boucle avec au moins une itération
    double normTQ8[6];                     //Pour stocker la différence entre les normes de la solution explicite et implicite à chaques itérations à la question 8
    double B[NX][NX], Bprime[NX][NX], BQ5[NX][NX];//Différentes matrices B pour les question 1, 4 et 5
    double C[NX][NX], CQ8[NX][NX];          //Différentes matrices C pour les question 6 et 8
    double b[NX], bprime[NX], bQ5[NX];      //Différentes vecteur b pour les question 1, 4 et 5
    double c[NX], cQ8[NX];                  //Différentes matrices C pour les question 6 et 8
    double Ti=750, Tf=25;                   //Condition de température initiale et finale
    int ts;                                 //Temps qui est saisi par l'utilisateur definissant quand on doit sortir la plaque de l'eau
    double T[NX], Tprime[NX], TQ5[NX], deltaT[NX], TQ8[NX];//Tous les vecteurs Températures qui vont nous etre utiles (respectivement) aux questions 1, 4, 5, 5 et 8
    for (int i=0; i<NX; i++) T[i]=0;        
    for(int i=0;i<6;i++) normTQ8[i]=0.0;    //On initialise T et normTQ8 pour eviter de se retrouver avec des valeurs ahurissantes plus tard
    
    //DEBUT DU MAIN
    
    printf("Au bout de combien de temps la plaque doit-elle etre retirée de l'eau (en seconde, uniquement nombre entier)?\n");
    scanf("%d", &ts);
    printf("\n------DEBUT DE LA PARTIE 1------\n\n");

    //Question 1
    
    Creation_Bb(B, b, dx, dt, a);                                                         //On crée B et b

    //Question 2

    solution_numeriqueEXPLI(B, b, dx, dt, Ti, Tf, ts, T);
    printf("\nVecteur T obtenu avec le pas de temps initial (méthode explicite):\n\n");   //affichage de T
    for (int i=0; i<NX; i++) printf(" %lf\n", T[i]); 

    //Question 3

    Creation_Bb(B, b, dx, dt, a);
    Evo_temporelle_Question3(B, b, dx, dt, Ti, Tf, T, evoT);
    /*printf("\nEvolution temporelle de T au milieu de la plaque :\n\n");                //DECOMMENTER CES DEUX LIGNES POUR L'EVOLUTION TEMPORELLE DE T AU MILIEU DE LA PLAQUE
    for (int i=0; i<6400; i=i+10) printf("%lf\n", evoT[i]);*/

    //Question 4

    Creation_Bb(Bprime, bprime, dx, dtprime, a);                                         //On crée B et b pour le pas de temps divisé par 5
    solution_numeriqueEXPLI(Bprime, bprime, dx, dtprime, Ti, Tf, ts, Tprime);
    printf("\nVecteur T avec le pas de temps initial divisé par 5 (méthode explicite):\n\n");    //affichage de T
    for (int i=0; i<NX; i++) printf(" %lf\n", Tprime[i]);
    //Pour le tracé de la différence entre T et Tprime, on a directement utilisé Excel en utilisant les valeurs que nous donnait le programme
    
    //Question 5
    
    dtQ5=dtprime;
    double TQ[NX];
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
    printf("\nVecteur T obtenu avec le pas de temps optimal (méthode explicite) :\n\n");
    for (int i=0; i<NX; i++) printf("%lf\n", TQ5[i]);
    printf("\nPas de temps optimal (méthode explicite) : %lf\n\n", dtQ5);

    printf("\n------DEBUT DE LA PARTIE 2------\n\n");

    //Question 6

    Creation_Cc(C, c, dx, dtQ5, a);
    
    //Question 7

    solution_numeriqueIMPLI(C, c, dx, dtQ5, Ti, Tf, ts, T);
    printf("Vecteur T question 7 (méthode implicite) :\n\n");    //affichage de T
    for (int i=0; i<NX; i++) printf(" %lf\n", T[i]);
    
    //Question 8

    dtQ8=dtQ5;
    z=0;
    do {
        Creation_Cc(CQ8, cQ8, dx, dtQ8, a);
        solution_numeriqueIMPLI(CQ8, cQ8, dx, dtQ8, Ti, Tf, ts, TQ8);
        for (int i=0; i<NX; i++) deltaT[i]=fabs(TQ5[i]-TQ8[i]);
        normTQ8[z] = norme_vect(deltaT);
        dtQ8=dtQ8*5;
        z++;
    }
    while (dtQ8<dt);
    printf("\nEvolution de le différence entre Solution Explicite et Implicite :\n\n");
    for(int i=0;i<6;i++) printf("%lf\n", normTQ8[i]);

return 0;
}
