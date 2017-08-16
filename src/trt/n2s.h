#ifndef N2S_SJM
#define N2S_SJM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


char digit_num[10]={'0','1','2','3','4','5','6','7','8','9'};
int NUM[5]={10,100,1000,10000,100000};

#ifdef __cplusplus
extern "C" {
#endif

char *n2s(int n);
void filename(char *file_tmp, char *CsgPre, int is, char *CsgSuffix);

#ifdef __cplusplus
}
#endif

void filename(char *file_tmp, char *CsgPre, int is, char *CsgSuffix)
{
   strcpy(file_tmp,CsgPre);
   if( is >= 0 ) {
      strcat(file_tmp,n2s(is));
   }
   if( strlen(CsgSuffix) > 0 ) {
      strcat(file_tmp,CsgSuffix);
   }
   return;
}

char *n2s(int n)
{
   int is_neg=0;
   if( n < 0 ) {
      is_neg=1;
      n=-1*n;
   }
   char *res;
   if( n>NUM[4] ) {
     n=(int)fmod(n,NUM[4]);
   }
   if( n < NUM[0] ) {
     res=(char *)malloc(2*sizeof(char)+is_neg);
     if( is_neg ) {
        res[0]='-';
     }
     res[is_neg]=digit_num[n];
     res[1+is_neg]='\0';
     return res;
   }
   int i,j,k;
   for( i=1;i<5;i++ ) {
      if( n<NUM[i] ) {
        res=(char *)malloc((i+2)*sizeof(char)+is_neg);
        if( is_neg ) {
           res[0]='-';
        }
        res[i+1+is_neg]='\0';
        k=(int)floor(n/NUM[i-1]);
        n=(int)fmod(n,NUM[i-1]);
        res[is_neg]=digit_num[k];
        for( j=1;j<i;j++ ) {
           k=(int)floor(n/NUM[i-j-1]);
           n=(int)fmod(n,NUM[i-j-1]);
           res[j+is_neg]=digit_num[k];
        }
	res[i+is_neg]=digit_num[n];
	break;
      }
    }
    return res;
}
        
   
#endif
