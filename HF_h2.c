/* molecola di H2 con Hartree-Fock */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 20
#define ABS(a)     (((a) < 0) ? -(a) : (a))
#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))

 static double e2 = 2.0;
 static double pi = 3.14159265358979;
 double c[N], pos_e, alpha[N];
 int n_alpha;
 FILE *out;
 FILE *outf;
//dichiarazione subroutine
 double fi(double );
 double fi2(double );
 void int_sovrapp();
 void int_coulom();
 void int_cinet();
 void int_biell();
 void int_F();
 int diag  (int dimension, double f[N][N], double s[N][N], double *e, double v[N][N] );
 double h[N][N], f[N][N],g[N][N][N][N];
 double s[N][N], v[N][N], vcu[N][N], e[N];
 double  t[N][N];
 char filin[80];
 FILE *in, *out;
 double alpha_ij, alpha_sum,r_ij;
 double alpha_lm, r_lm;
 double Er, ER[600];
 double dr=0.1,r_a=0., r_b, r[N];
 double  aa, eold, enew; 
 int main(){
  out=fopen("energia","w");  
/* variables */
  double sqrt();   
  int ia,ib,ic,id, i, iter, maxiter=100;
  int j,k,l,m;
  double factor;
//inizializzo Er
  for(i=0;i<600;i++){
     ER[i]=0.;
    }
/*   Input data */
  fprintf(stdout, " Parameters of the Gaussians from file >> ");
  scanf("%80s", filin);
  in = fopen (filin,"r");
  fscanf(in,"%d", &n_alpha);
  if ( n_alpha > N)  {
    fprintf(stderr, "n_alpha (%d) > nalpx (%d)\n",n_alpha, N);
     exit (1);
  }
  /* Read parameters of the gaussians */
  for ( ia = 0 ; ia < n_alpha ; ia++ ) {
     fscanf(in,"%lf", &alpha[ia]);
     fprintf (stdout," Gaussian # %3d  coefficient >> %lf\n",ia, alpha[ia]);
  }
  fclose(in);
  /*starting solution very lousy guess*/
  for(i=1;i<N;i++){
     c[i]=0.0;
   }
  c[1]=1;
 /*1)loop sugli R */
  for(k=1;k<100;k++){
     r_b=dr*k;
 //riempio gli r
  for(i=0;i<n_alpha;i++){
     if(i<(n_alpha/2)){
        r[i]=0.;
      }
     else{
       r[i]=r_b;
      }
   }
//riempio le matrici
//integrale di sovrapposizione
  int_sovrapp();
//integrale di interazione coulombiana
  int_coulom();
//integrale cinetico
  int_cinet();
//integrale di interazione elettrone elettrone
  int_biell();
//hamiltoniano
  for(i=0;i<n_alpha;i++){
      for(j=0;j<n_alpha;j++){
         h[i][j]=t[i][j]+vcu[i][j];
          }
     }
//normalizzo gli autovettori
  factor=0.;
  for(i=0;i<n_alpha;i++){
     for(j=0;j<n_alpha;j++){
         factor+=s[i][j];
        } 
     }
  for(i=0;i<n_alpha;i++){
      c[i]=1./sqrt(factor);
     }
/*3) ciclo autoconsistente: calcolo F, diagonalizzo, trascrivo i coefficienti, calcolo energia*/
  enew=0.;
  for(iter=0;iter<maxiter; iter++){
     int_F();
     diag ( n_alpha,f,s,e,v);
     for(i=0;i<n_alpha;i++){
         c[i]=v[i][0];
         }
     eold=enew; 
     enew=0.0; 
     for(i=0;i<n_alpha;i++){
         for(j=0;j<n_alpha;j++){
             enew = enew + 2.0* (h[i][j])*c[i]*c[j];
             for(l=0;l<n_alpha;l++){
 	         for(m=0;m<n_alpha;m++){
 	             enew = enew + g[i][j][l][m]*c[i]*c[j]*c[l]*c[m];
                     }
                 }
             }
         }
  if ( ABS (enew-eold) < 5.0e-13){
     fprintf (stdout, "\n Convergence achieved at iter: %d, stopping\n",iter);
     goto L1; 
     }
    }//chiudo il ciclo sulle iterazioni
  fprintf (stdout, "\n Convergence not reached, stopping\n");
  exit(2);
  L1: 
/*4) aggiungo energia di repulsione dei nuclei e trascrivo sul file*/
  Er=1./r_b;
  ER[k]=Er+enew;
  fprintf(out, "%lf \t %lf \t %lf \n", r_b, enew+Er, enew);
  printf("# R=%lf \t E=%lf \n", r_b, enew);
  }//chiudo il ciclo sugli R
/*5)trovo L'R in corrispondenza di cui si ha il minimo di E e ricalcolo i coefficienti (questo per lo stato di equilibrio) altrimenti metto l'R che voglio e ricalcolo i coefficienti*/

//trovo il minimo di E e di conseguenza l'R
  double emin=0.;
  for(i=0;i<600;i++){
     if(ER[i]<emin){
        emin=ER[i];
        r_b=i*dr;
        }
     }
  printf(" Rmin: %g\t emin:%g\t\n", r_b,emin);
/*2)Riempio tutte le matrici: S,V,T,G e H per l'R desiderato */ 
 //r_b=40;
//integrale di sovrapposizione
  int_sovrapp();
//integrale di interazione coulombiana
  int_coulom();
//integrale cinetico
  int_cinet();
//integrale di interazione elettrone elettrone
  int_biell();
//hamiltoniano
  for(i=0;i<n_alpha;i++){
      for(j=0;j<n_alpha;j++){
          h[i][j]=t[i][j]+vcu[i][j];
          }
      }
//normalizzo gli autovettori
  factor=0.;
  for(i=0;i<n_alpha;i++){
     for(j=0;j<n_alpha;j++){
        factor+=s[i][j];
        }
     }
  for(i=0;i<n_alpha;i++){
      c[i]=1./sqrt(factor);
      }
/*3) ciclo autoconsistente: calcolo F, diagonalizzo, trascrivo i coefficienti, calcolo energia*/
  enew=0.;
  for(iter=0;iter<maxiter; iter++){
      int_F();
      diag ( n_alpha,f,s,e,v);
      for(i=0;i<n_alpha;i++){
          c[i]=v[i][0];
         }
      eold=enew; 
      enew=0.0; 
     for(i=0;i<n_alpha;i++){
         for(j=0;j<n_alpha;j++){
             enew = enew + 2.0* (h[i][j])*c[i]*c[j];
             for(l=0;l<n_alpha;l++){
 	         for(m=0;m<n_alpha;m++){
 	             enew = enew + g[i][j][l][m]*c[i]*c[j]*c[l]*c[m];
                     }
                 }
             }
          }
  
  if ( ABS (enew-eold) < 5.0e-13){
     fprintf (stdout, "\n Convergence achieved at iter: %d, stopping\n",iter);
     goto L2; 
      }
    }//chiudo il ciclo sulle iterazioni
  fprintf (stdout, "\n Convergence not reached, stopping\n");
  exit(2);
  printf("enew di r=%lg è pari a =%lg\n",r_b, enew);
  L2: 
/*6)trascrivo le funzioni d onda con i moduli quadri */
  printf("enew di r=%lg è pari a =%lg\n",r_b, enew);
  outf = fopen ("funzione", "w");
  double *y, *y2;
  double xmax=40, dx=0.1;
  int nx,iR;
  double fi(double x);
  double fi2(double x);
  nx=(int) (xmax/dx);
  y  = (double *) malloc (2*nx * sizeof(double));
  y2 = (double *) malloc (2*nx * sizeof(double));
  for (iR=0; iR<2*nx; iR++){
      y[iR]  = 0;
      y2[iR] = 0;
       }
  for (iR=0; iR<nx; iR++){
       y[nx+iR]  = fi(iR*dx);
       y[nx-iR]  = fi(-iR*dx);
       y2[nx+iR] = fi2(iR*dx);
       y2[nx-iR] = fi2(-iR*dx);
      }
  for (iR=1; iR<2*nx; iR++){
      fprintf(outf, "%lf \t %lf \t %lf\n", iR*dx - nx*dx, y[iR], y2[iR]);
      }
  FILE *outc;
  outc = fopen("coefficienti", "w");
  for (ia = 0; ia < n_alpha; ia++){
       fprintf(outc, "%d \t %f \t %f \n", ia, alpha[ia], c[ia]);
	}
  fclose(outc);
  fclose(outf);
  fclose(out);
  return 0;
  }

//subroutine funzioni d'onda
 double fi(double x){
        int a;
	double temp = 0;
	for (a=0; a<(n_alpha); a++){
	    temp += c[a]*exp(-alpha[a] * (x+(r_b/2)) * (x+(r_b/2)));
	    temp += c[a]*exp(-alpha[a] * (x-(r_b/2)) * (x-(r_b/2)));
	    }
	return temp;
 }

double fi2(double x){
        int a,b;
	double temp = 0;
	double Ra, Rb;
	for (a=0; a<n_alpha; a++){
		for (b=0; b<n_alpha; b++){
 	        temp += c[a]*c[b]*exp(-alpha[a] * (x+(r_b/2)) * (x+(r_b/2)))*exp(-alpha[b] * (x+(r_b/2)) * (x+(r_b/2)));
		temp += c[a]*c[b]*exp(-alpha[a] * (x-(r_b/2)) * (x-(r_b/2)))*exp(-alpha[b] * (x-(r_b/2)) * (x-(r_b/2)));
		}
	}
	return temp;
}

//subroutines integrali
void int_sovrapp(){
     int i,j;
     for(i=0;i<n_alpha;i++){
         for(j=0;j<n_alpha;j++){
             alpha_sum=alpha[i]+alpha[j];
             alpha_ij=(alpha[i]*alpha[j])/(alpha_sum);
             s[i][j]=pow((pi/alpha_sum),1.5)*exp(-alpha_ij*(r[i]-r[j])*(r[i]-r[j]));
             }
         }
     }

void int_coulom(){
     int i,j;
     double temp,arg_a,arg_b, f0_a, f0_b;
     for(i=0;i<n_alpha;i++){
         for(j=0;j<n_alpha;j++){
             alpha_sum=alpha[i]+alpha[j];
             alpha_ij=(alpha[i]*alpha[j])/(alpha_sum);
             r_ij=((alpha[i]*r[i])+(alpha[j]*r[j]))/(alpha_sum);
             temp=1.;
             arg_a=r_ij-r_a;
             arg_b=r_ij-r_b;
             if(arg_a==0){
                f0_a=-2*pi/alpha_sum*exp(-alpha_ij*(r[i]-r[j])*(r[i]-r[j]));
                }
             else{
                 f0_a=-s[i][j]/fabs(r_ij-r_a)*erf(sqrt(alpha_sum)*fabs(r_ij-r_a));
                 }
             if(arg_b==0){
                f0_b=-2*pi/alpha_sum*exp(-alpha_ij*(r[i]-r[j])*(r[i]-r[j]));
                }
             else{
                f0_b=-s[i][j]/fabs(r_ij-r_b)*erf(sqrt(alpha_sum)*fabs(r_ij-r_b));
               }
            temp=temp*(f0_a+f0_b);
            vcu[i][j]=temp;
          }
        } 
  }

  void int_cinet(){
       int i,j;
       for(i=0;i<n_alpha;i++){
           for(j=0;j<n_alpha;j++){
               alpha_sum=alpha[i]+alpha[j];
               alpha_ij=(alpha[i]*alpha[j])/(alpha_sum);
               t[i][j]=0.5*(alpha_ij*(6.-4.*alpha_ij*(r[i]-r[j])*(r[i]-r[j]))*s[i][j]);
                }
           }
   }

  void int_biell(){
       int i,j,l,m;
       double argo;
       for(i=0;i<n_alpha;i++){
           for(j=0;j<n_alpha;j++){
               for(l=0;l<n_alpha;l++){
                   for(m=0;m<n_alpha;m++){
                       alpha_ij=(alpha[i]+alpha[j]);
                       alpha_lm=(alpha[l]+alpha[m]);
                       r_ij=((alpha[i]*r[i])+(alpha[j]*r[j]))/(alpha_ij);
                       r_lm=((alpha[l]*r[l])+(alpha[m]*r[m]))/(alpha_lm);
                       argo=fabs(r_ij-r_lm);
                       if(argo==0){
                          g[i][j][l][m]=2*s[i][j]*s[l][m]*sqrt(alpha_ij*alpha_lm/(pi*(alpha_ij+alpha_lm)));
                          }
                       else{
                           g[i][j][l][m]=s[i][j]*s[l][m]/fabs(r_ij-r_lm)*erf(sqrt(alpha_ij*alpha_lm/(alpha_ij+alpha_lm))*fabs(r_ij-r_lm));
                          }
                      }
                   }
                }
           }
  }
//matrice F
 void int_F(){
      int i,j,l,m;
      for(i=0;i<n_alpha;i++){
          for(j=0;j<n_alpha;j++){
              f[i][j]=h[i][j];
              for(l=0;l<n_alpha;l++){
                  for(m=0;m<n_alpha;m++){
                      f[i][j]+=c[l]*g[i][j][l][m]*c[m];
                     }
                  }
               }
            }
      }
       
/* subroutine diag */
int diag(int n, double f[N][N], double s[N][N], double *e, double v[N][N])
{
  /*    Finds eigenvalues and eigenvectors of the generalized problem
	Hv=eSv, where H=hermitian matrix, S=overlap matrix */

  /* On input: n = dimension of the matrix to be diagonalized
               h = matrix to be diagonalized
               s = overlap matrix
     On output:
               e = eigenvalues
               v = eigenvectors
               s and h are unchanged */

  /* LOCAL variables */
  int lwork, i, j, k, nn, lda, info;
  /* lwork = dimension of workspace for lapack routine  */
  static double small = 1.0e-10;
  static char *V = "V";
  static char *U = "U";
  double *work, aux[N][N], uno[N][N];
 
  lwork=3*n;
  work = (double *) malloc( lwork * sizeof (double));

  /* Copy S into an auxiliary matrix (dsyev destroys the matrix) */
  /* The indices of aux are transposed for fortran compatibility */
  for ( i = 0; i < n; ++i ) {
    for ( j = 0; j < n; ++j ) {
      aux[j][i] = s[i][j];
    }
  }

  /*  Diagonalize S  */

  lda = N;

  dsyev_ ( V, U, &n, aux, &lda, e, work, &lwork, &info ) ;

  if ( info !=0 ) {
    fprintf (stderr, "S-matrix diagonalization failed\n");
    exit (1);
  }

  /*    Keep only linearly independent combinations
	(within a given threshold)  */

  nn = 0;
  for ( i = 0; i < n; ++i ) {
    if ( 

e[i] > small) {
      for ( j = 0; j < n; ++j ) {
	aux[nn][j] = aux[i][j] / sqrt(e[i]);
      }
      ++nn;
    }
  }

  if ( nn < n) fprintf (stdout, " # of linearly independent vectors = %d\n", nn);
  
  /*       Trasform H    */
  
  for ( i = 0; i < n; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      v[i][j] = 0.0;
      for  ( k = 0; k < n; ++k ) {
	v[i][j] = v[i][j] + f[i][k] * aux[j][k];
      }
    }
  }
   /* The indices of uno is transposed for fortran compatibility */
  for ( i = 0; i < nn; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      uno[j][i] = 0.0;
      for  ( k = 0; k < n; ++k ) {
	uno[j][i] = uno[j][i] + aux[i][k] * v[k][j];
      }
    }
  }
  
  /*  Diagonalize H  */
  
  dsyev_ ("V", "U", &nn, uno, &lda, e, work, &lwork, &info );

  if ( info !=0 ) {
    fprintf (stderr, "H-matrix diagonalization failed\n");
    exit (1);
  }

  for ( i = 0; i < n; ++i ) {
    for ( j = 0; j < nn; ++j ) {
      v[i][j] = 0.0;
      for  ( k = 0; k < n; ++k ) {
	v[i][j] = v[i][j] + aux[k][i] * uno[j][k];
      }
    }
  }
  free(work);
  return 0;
 }


