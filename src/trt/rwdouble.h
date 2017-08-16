#ifndef RWDOUBLE_SJM
#define RWDOUBLE_SJM

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void freaddouble(FILE *fp, size_t n, double *data);
void freaddouble_2(char *fname, char *type, size_t n, double *data);
void fwritedouble(FILE *fp, size_t n, double *data);
void fwritedouble_2(char *fname, char *type, size_t n, double *data);
void freaddouble1d(char *fname, size_t n, double *data);
void fwritedouble1d(char *fname, size_t n, double *data);
void freadfloat2d(char *fname, int nz, int nx, float **data);
void freadfloat2d_2(FILE *fp, int nz, int nx, float **data);
void freaddouble2d(char *fname, int nx, int nz, double **data);
void fwritedouble2d(char *fname, int nx, int nz, double **data);
void fwritefloat2d(char *fname, int nz, int nx, float **data);
void fwritefloat2d_2(FILE *fp, int nz, int nx, float **data);

#ifdef __cplusplus
}
#endif

void freadfloat2d_2(FILE *fp, int nz, int nx, float **data)
{
   int i,j;
   float *tmp=ealloc1float(nx*nz);
   fread(tmp,sizeof(float),nx*nz,fp);
   for( i=0;i<nz;i++ ) {
      for( j=0;j<nx;j++ ) {
         data[j][i]=tmp[i*nz+j];
      }
   }
   fclose(fp);
   free(tmp);
   return;
}

void freadfloat2d(char *fname, int nz, int nx, float **data)
{
   int i,j;
   FILE *fp=efopen(fname,"rb");
   float *tmp=ealloc1float(nx*nz);
   fread(tmp,sizeof(float),nx*nz,fp);
   for( i=0;i<nz;i++ ) {
      for( j=0;j<nx;j++ ) {
         data[j][i]=tmp[i*nz+j];
      }
   }
   fclose(fp);
   free(tmp);
   return;
}

void freaddouble2d(char *fname, int nx, int nz, double **data)
{
   int i,j;
   FILE *fp=efopen(fname,"rb");
   float *tmp=ealloc1float(nx*nz);
   fread(tmp,sizeof(float),nx*nz,fp);
   for( i=0;i<nx;i++ ) {
      for( j=0;j<nz;j++ ) {
         data[j][i]=(double)tmp[i*nz+j];
      }
   }
   fclose(fp);
   free(tmp);
   return;
}

void fwritefloat2d(char *fname, int nz, int nx, float **data)
{
   int i,j;
   FILE *fp=efopen(fname,"wb");
   float *tmp=ealloc1float(nx*nz);
   for( i=0;i<nz;i++ ) {
      for( j=0;j<nx;j++ ) {
         tmp[i*nz+j]=data[j][i];
      }
   }
   fwrite(tmp,sizeof(float),nx*nz,fp);
   fclose(fp);
   free(tmp);
   return;
}

void fwritefloat2d_2(FILE *fp, int nz, int nx, float **data)
{
   int i,j;
   float *tmp=ealloc1float(nx*nz);
   for( i=0;i<nz;i++ ) {
      for( j=0;j<nx;j++ ) {
         tmp[i*nz+j]=data[j][i];
      }
   }
   fwrite(tmp,sizeof(float),nx*nz,fp);
   fclose(fp);
   free(tmp);
   return;
}

void fwritedouble2d(char *fname, int nx, int nz, double **data)
{
   int i,j;
   FILE *fp=efopen(fname,"wb");
   float *tmp=ealloc1float(nx*nz);
   for( i=0;i<nx;i++ ) {
      for( j=0;j<nz;j++ ) {
         tmp[i*nz+j]=data[j][i];
      }
   }
   fwrite(tmp,sizeof(float),nx*nz,fp);
   fclose(fp);
   free(tmp);
   return;
}

void freaddouble1d(char *fname, size_t n, double *data)
{
   int i;
   FILE *fp=efopen(fname,"rb");
   float *tmp=ealloc1float(n);
   fread(tmp,sizeof(float),n,fp);
   for( i=0;i<n;i++ ) {
      data[i]=(double)tmp[i];
   }
   free(tmp);
   fclose(fp);
   return;
}

void freaddouble(FILE *fp, size_t n, double *data)
{
   int i;
   float *tmp=ealloc1float(n);
   fread(tmp,sizeof(float),n,fp);
   for( i=0;i<n;i++ ) {
      data[i]=(double)tmp[i];
   }
   free(tmp);
   return;
}

void fwritedouble1d(char *fname, size_t n, double *data)
{
   int i;
   FILE *fp=efopen(fname,"wb");
   float *tmp=ealloc1float(n);
   for( i=0;i<n;i++ ) {
      tmp[i]=data[i];
   }
   fwrite(tmp,sizeof(float),n,fp);
   free(tmp);
   fclose(fp);
   return;
}

void fwritedouble(FILE *fp, size_t n, double *data)
{
   int i;
   float *tmp=ealloc1float(n);
   for( i=0;i<n;i++ ) {
      tmp[i]=data[i];
   }
   fwrite(tmp,sizeof(float),n,fp);
   free(tmp);
   return;
}

void freaddouble_2(char *fname,  char *type, size_t n, double *data)
{
   int i;
   FILE *fp=efopen(fname,type);
   float *tmp=ealloc1float(n);
   fread(tmp,sizeof(float),n,fp);
   for( i=0;i<n;i++ ) {
      data[i]=(double)tmp[i];
   }
   free(tmp);fclose(fp);
   return;
}

void fwritedouble_2(char *fname, char *type, size_t n, double *data)
{
   int i;
   FILE *fp=efopen(fname,type);
   float *tmp=ealloc1float(n);
   for( i=0;i<n;i++ ) {
      tmp[i]=data[i];
   }
   fwrite(tmp,sizeof(float),n,fp);
   free(tmp);fclose(fp);
   return;
}

#endif
