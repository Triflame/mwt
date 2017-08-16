#ifndef EALLOC_SJM
#define EALLOC_SJM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cwp_sjm.h"

#ifdef __cplusplus
extern "C" {
#endif

int *ealloc1int(size_t n);
int *realloc1int(int *v, size_t n1);
int **ealloc2int(size_t n1,size_t n2);
int ***ealloc3int(size_t n1,size_t n2,size_t n3);
void free1int(int *v);
void free2int(int **v);
void free3int(int ***v);
float *ealloc1float(size_t n);
float *realloc1float(float *v, size_t n1);
float **ealloc2float(size_t n1,size_t n2);
float ***ealloc3float(size_t n1,size_t n2,size_t n3);
void free1float(float *v);
void free2float(float **v);
void free3float(float ***v);
double *ealloc1double(size_t n);
double *realloc1double(double *v, size_t n1);
double **ealloc2double(size_t n1,size_t n2);
double ***ealloc3double(size_t n1,size_t n2,size_t n3);
void free1double(double *v);
void free2double(double **v);
void free3double(double ***v);
char *ealloc1char(size_t n);
char **ealloc2char(size_t n1,size_t n2);
void free1char(char *p);
void *ealloc1(size_t n, size_t size);
void *realloc1(void *v, size_t n, size_t size);
void free1(void *p);
void **ealloc2 (size_t n1, size_t n2, size_t size);
void free2 (void **p);
void ***ealloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void free3 (void ***p);
FILE *efopen(char *fname,char *type);

#ifdef __cplusplus
}
#endif


int *ealloc1int(size_t n)
{
   return (int *)ealloc1(n,sizeof(int));
}

int *realloc1int(int *v, size_t n1)
{
   return (int *)realloc1(v,n1,sizeof(int));
}

int **ealloc2int(size_t n1,size_t n2)
{
   return (int **)ealloc2(n1,n2,sizeof(int));
}

int ***ealloc3int(size_t n1,size_t n2,size_t n3)
{
   return (int ***)ealloc3(n1,n2,n3,sizeof(int));
}

void free1int(int *v)
{
   free1(v);
}
 
void free2int(int **v)
{
   free2((void**)v);
}

void free3int(int ***v)
{
   free3((void***)v);
}

float *ealloc1float(size_t n)
{
   return (float *)ealloc1(n,sizeof(float));
}

float *realloc1float(float *v, size_t n1)
{
   return (float *)realloc1(v,n1,sizeof(float));
}

float **ealloc2float(size_t n1,size_t n2)
{
   return (float **)ealloc2(n1,n2,sizeof(float));
}

float ***ealloc3float(size_t n1,size_t n2,size_t n3)
{
   return (float ***)ealloc3(n1,n2,n3,sizeof(float));
}

void free1float(float *v)
{
   free1(v);
}
 
void free2float(float **v)
{
   free2((void**)v);
}

void free3float(float ***v)
{
   free3((void***)v);
}

double *ealloc1double(size_t n)
{
   return (double *)ealloc1(n,sizeof(double));
}

double *realloc1double(double *v, size_t n1)
{
   return (double *)realloc1(v,n1,sizeof(double));
}

double **ealloc2double(size_t n1,size_t n2)
{
   return (double **)ealloc2(n1,n2,sizeof(double));
}

double ***ealloc3double(size_t n1,size_t n2,size_t n3)
{
   return (double ***)ealloc3(n1,n2,n3,sizeof(double));
}

void free1double(double *v)
{
   free1(v);
}
 
void free2double(double **v)
{
   free2((void**)v);
}

void free3double(double ***v)
{
   free3((void***)v);
}

char *ealloc1char(size_t n)
{
   return (char *)ealloc1(n,sizeof(char));
}

char **ealloc2char(size_t n1,size_t n2)
{
   return (char **)ealloc2(n1,n2,sizeof(char));
}

void free1char(char *v)
{
   free1(v);
}

void *ealloc1(size_t n, size_t size)
{
   void *res;
   res=malloc(n*size);
   if( !res ) {
     printf("memory request failed.\n");
     exit(1);
   }
   return res;
}

void *realloc1(void *v, size_t n, size_t size)
{
   void *res;
   res=realloc(v,n*size);
   if( !res ) {
     printf("memory request failed.\n");
     exit(1);
   }
   return res;
}

void free1(void *p)
{
    free(p);
}

void **ealloc2 (size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;

        if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) {
          printf("memory request failed.\n");
          exit(1);
        }
        if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
          free(p);
          printf("memory request failed.\n");
          exit(1);
        }
        for (i2=0; i2<n2; i2++)
          p[i2] = (char*)p[0]+size*n1*i2;
        return p;
}
void free2 (void **p)
{
    free(p[0]);
    free(p);
}

void ***ealloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
        size_t i3,i2;
        void ***p;

        if ((p=(void***)malloc(n3*sizeof(void**)))==NULL) {
          printf("memory request failed.\n");
          exit(1);
        }
        if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
          free(p);
          printf("memory request failed.\n");
          exit(1);
        }
        if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
          free(p[0]);
          free(p);
          printf("memory request failed.\n");
          exit(1);
        }

        for (i3=0; i3<n3; i3++) {
                p[i3] = p[0]+n2*i3;
                for (i2=0; i2<n2; i2++)
                        p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
        }
        return p;
}

void free3 (void ***p)
{
        free(p[0][0]);
        free(p[0]);
        free(p);
}

FILE *efopen(char *fname,char *type)
{
    FILE *fp;
    if( (fp=fopen(fname,type))==NULL ) {
      printf("open file %s failed.\n",fname);
      exit(1);
    }
    return fp;
}

#endif
