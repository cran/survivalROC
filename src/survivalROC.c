/* ---------------------------------------------------------------------- */
/*  survivalROC.c                                                        */
/* ---------------------------------------------------------------------- */
/* AUTHOR: P. HEAGERTY  */
/* COMMENTS BY: P. SAHA */
/* DATE: 06/03/08       */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define cfree free

/* ---------------------------------------------------------------------- */

void  survivalROC( S_time, S_status, S_ObsTimes, 
                   S_x, S_UniqueValues,
                   S_PredictTime, S_SurvT, S_span, 
                   S_TP, S_FP,
                   S_n, S_p, S_q )

  double *S_time, *S_status, *S_x, *S_UniqueValues, *S_PredictTime;
  double *S_ObsTimes, *S_SurvT, *S_span, *S_TP, *S_FP;
  int    *S_n, *S_p, *S_q;
{
  double *X, *UniqueValues, *SurvivalTime, *Status, PredictTime, span;
  double *ObsTimes, St, *SgivenX, *pX, *cdfX, STX, alt_St;
  double num, den, tj, ProdLimit, p1;
  int    n, p, q;
  int    i, j, k;
  int    index0, index1, IndexSpan, FirstIndex, LastIndex, TargetIndex;

  /* printf("[0]\n"); */
/* Step: INITIALIZATION -- setting the counters and initializing the variables */
  n = *S_n; /* total no of observations   */
  q = *S_q; /* no of unique failure times */
  p = *S_p; /* no of unique marker values */
  span = *S_span; /* span for the other method (from paper, Akritas) ?? */
  PredictTime = *S_PredictTime; /* this is the cut-off we want to use, e.g., 1yr, 2yrs, etc. */

  X = (double *)malloc((unsigned) (n)*sizeof(double));
  for( j=0; j<n; j++ ) X[j] = *(S_x+j);

  UniqueValues = (double *)malloc((unsigned) (q)*sizeof(double));
  for( j=0; j<q; j++ ) UniqueValues[j] = *(S_UniqueValues+j);

  SurvivalTime = (double *)malloc((unsigned) (n)*sizeof(double));
  for( j=0; j<n; j++ ) SurvivalTime[j] = *(S_time+j);

  Status = (double *)malloc((unsigned) (n)*sizeof(double));
  for( j=0; j<n; j++ ) Status[j] = *(S_status+j);

  ObsTimes = (double *)malloc((unsigned) (p)*sizeof(double));
  for( j=0; j<p; j++ ) ObsTimes[j] = *(S_ObsTimes+j);

  /* printf("[1]\n");*/

  /* *********************************** */
  /*   assume that data are passed with  */
  /*   X sorted                          */
  /*                                     */
  /*  BE SURE to sort in S!!!            */
  /* *********************************** */

  /* *********************************** */
  /*   assume that data are passed with  */
  /*   ObsTimes sorted                   */
  /*                                     */
  /*  BE SURE to sort in S!!!            */
  /* *********************************** */

  /* overall survival */

  ProdLimit = 1.0;

  for( j=0; j<p; j++ ){

    tj = ObsTimes[j];

    num = 0.0;
    den = 0.0;

    if( tj <= PredictTime ){
      for( k=0; k<n; k++ ){
	if( SurvivalTime[k] >= tj ){
	  den += 1.0;
	  if( SurvivalTime[k]==tj & Status[k]==1.0 ) num += 1.0;
	}
      }
    }
    if( den>0.0 ) ProdLimit = ProdLimit * ( 1.0 - num/den );
  }/* j */

  St = ProdLimit;

  /* printf("\n   St = %f\n\n", St );*/

  /* printf("[2]\n"); */

  /* conditional survival */

  SgivenX = (double *)malloc((unsigned) (q)*sizeof(double));
  for( j=0; j<q; j++ ) SgivenX[j] = 0.0;

  pX = (double *)malloc((unsigned) (q)*sizeof(double));
  for( j=0; j<q; j++ ) pX[j] = 0.0;
  cdfX = (double *)malloc((unsigned) (q)*sizeof(double));
  for( j=0; j<q; j++ ) cdfX[j] = 0.0;

  IndexSpan = (int)floor( ((double)(n)*span + 0.5 ) );

  /* printf("\n  IndexSpan = %i\n\n", IndexSpan ); */

  for( i=0; i<q; i++ ){

    for( j=0; j<n; j++ ){
      if( X[j]==UniqueValues[i] ) pX[i] += 1.0/( (double)(n) );
      if( X[j]<=UniqueValues[i] ) cdfX[i] += 1.0/( (double)(n) );
    }

    /* symmetric window */

    FirstIndex = n-1;    
    LastIndex = n-1;    
    for( j=0; j<(n-1); j++ ){
      if( X[j] < UniqueValues[i] & X[j+1] >= UniqueValues[i] ) FirstIndex=j;
      if( X[j] <= UniqueValues[i] & X[j+1] > UniqueValues[i] ) LastIndex=j;
    }
    TargetIndex = floor( ( (double)FirstIndex + 
                           (double)LastIndex )/2.0 + 0.5 );
 
    
    index0 = TargetIndex - IndexSpan;
    if( index0<0 ) index0=0;
    if( index0 > FirstIndex ) index0 = FirstIndex;

    index1 = TargetIndex + IndexSpan;
    if( index1>=n ) index1=(n-1);
    if( index1 < LastIndex ) index1 = LastIndex;

    ProdLimit = 1.0;

    for( j=0; j<p; j++ ){

      tj = ObsTimes[j];

      num = 0.0;
      den = 0.0;

      if( tj <= PredictTime ){
	for( k=0; k<=(index1-index0); k++ ){
	  if( SurvivalTime[index0+k] >= tj ){
	    den += 1.0;
	    if( SurvivalTime[index0+k]==tj & Status[index0+k] == 1.0 ) num += 1.0;
	  }
	}
      }
      if( den > 0.0 ) ProdLimit = ProdLimit * ( 1.0 - num/den );
	
    }/* j */

    SgivenX[i] = ProdLimit;

  }/* i */
    
  /* printf("[3]\n"); */ 

 alt_St = 0.0;
 p1 = 0.0;
 for( j=0; j<q; j++ ){
   if( j>0 ) pX[j] = cdfX[j] - cdfX[j-1];
   alt_St += SgivenX[j] * pX[j];
   p1 += pX[j];
 }

 /* printf("alt_St = %f, p1 = %f \n", alt_St, p1 ); */

 St = alt_St;
 
  for( i=0; i<(q-1); i++ ){
    STX = 0.0;
    p1 = 0.0;
    for( j=(i+1); j<q; j++ ){
      STX += SgivenX[j] * pX[j];
      p1 += pX[j];
    }
    *(S_TP+i) = (p1 - STX)/(1.0-St);
    *(S_FP+i) = STX/St;
  }

*S_SurvT = St;

/* printf("[4]\n");  */
  

cfree( X );
cfree( UniqueValues );
cfree( SurvivalTime );
cfree( Status );
cfree( ObsTimes );
cfree( SgivenX );
cfree( pX ); 
cfree( cdfX ); 


return;

  }
/* ---------------------------------------------------------------------- */





