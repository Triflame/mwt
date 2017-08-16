#ifndef EIK2D_SJM
#define EIK2D_SJM

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "cwp_sjm.h"

#ifdef __cplusplus
extern "C" {
#endif

void eik2d(int nx, int nz, double *ttb,
           double dx, double *slow,
           double xs, double zs);
void corners(double *ttb, double *slow,
             double hh, int nx, int nz);
void eband(int *ix_front2, int *iz_front2,
           double *ttb_front2, int *front_size2,
           int n1, int k1, int nx, int nz,
           double hh, double *slow, double *ttb);
void add(int *ix_front2, int *iz_front2,
         double *ttb_front2, int *front_size2,
         int n1, int k1, double trt);
void merge(int *ix_front, int *iz_front, double *ttb_front,
           int front_big, int *front_size,
           int *ix_front2, int *iz_front2,
           double *ttb_front2, int front_size2);
void update(int *ix_front2, int *iz_front2,
            double *ttb_front2, int front_size2,
            double *ttb, int nz);
void sort(int *ix_front, int *iz_front,
          double *ttb, int front_size);
void sort_old(int *ix_front, int *iz_front,
          int nz,
          double *ttb, int front_size);
void source(double xs, double zs,
            int ix0, int iz0, double dx,
            int nx, int nz,
            double *slow, double *ttb,
            int *ix_front, int *iz_front,
            double *ttb_front, int *front_size);
int nround(double a);
void fgets2(char *str, int num, FILE *stream) ;
double max1(double x, double y);
double max(double x, double y);
int max2(int x, int y);
double min1(double x, double y);
double min(double x, double y);
int min2(int x, int y);

#ifdef __cplusplus
}
#endif

void eik2d(int nx, int nz, double *ttb,
           double dx, double *slow,
	   double xs, double zs)
{
  int i,j;
  double vmax=0.0,ssmin=1.0;
  int k;
  if( slow[0] > 1.0 ){
    for( i=0;i<nx;i++ ){
      for( j=0;j<nz;j++ ){
	k=i*nz+j;
	vmax=max(vmax,slow[k]);
        slow[k]=1./slow[k];
        ttb[k]=-2.0;
      }
    }
    ssmin=1./vmax;
  } else {
    for( i=0;i<nx;i++ ){
      for( j=0;j<nz;j++ ){
        k=i*nz+j;
        ssmin=min(ssmin,slow[k]);
        ttb[k]=-2.0;
      }
    }
    vmax=1./ssmin;
  }
    
  double hh=dx*dx;
  int ix0=NINT(xs/dx);
  int iz0=NINT(zs/dx);
  ttb[ix0*nz+iz0]=sqrt(pow(ix0*dx-xs,2)
                    +pow(iz0*dx-zs,2))*slow[ix0*nz+iz0];
  double ddx=1.0;
  double dt=ddx*dx/vmax;
/*  printf("dt=%12.10f\n",dt);exit(5);*/
  int ncp=20*(nx+nz);
  int ix_front[ncp],iz_front[ncp],front_size; 
  double ttb_front[ncp];

  source(xs,zs,ix0,iz0,dx,nx,nz,slow,
         ttb,ix_front,iz_front,
         ttb_front,&front_size);
  sort(ix_front,iz_front,ttb_front,front_size);
  double tmin;

  int nt=front_size;
  int jj=0,kkk=0;
  while( nt >= 1 ) {
    int ix_front2[ncp],iz_front2[ncp];
    double ttb_front2[ncp];
    int front_size2=0;
    int front_big=front_size;
    tmin=ttb_front[0];
/*
      for( i=0;i<nt;i++ ) {
        printf("%d %d %d %12.9f\n",i+1,ix_front[i],iz_front[i],ttb_front[i]);
      }
*/
/*
      if( (tmin>0.7) && (tmin < 0.8) ) {
         int ind0=63*nz;
         for( j=19;j<=25;j++ ) {
           printf("%d:%f %f %f %f %f %f %f\n",j,
             ttb[ind0-3*nz+j],ttb[ind0-2*nz+j],
             ttb[ind0-nz+j],ttb[ind0+j],
             ttb[ind0+nz+j],ttb[ind0+2*nz+j],
             ttb[ind0+3*nz+j]);
         }
       }
*/
    for( i=0;i<nt;i++ ) {
      int n1=ix_front[i];
      int k1=iz_front[i];
      if( ttb[n1*nz+k1] > (tmin+dt) ) {
        front_big=i;
        break;
      }
      eband(ix_front2,iz_front2,ttb_front2,&front_size2,
            n1,k1,nx,nz,hh,slow,ttb);
/*
      printf("i=%d size2=%d\n",i,front_size2);
      for( int kk=0;kk<front_size2;kk++ ) {
        printf("%d:(%d,%d):%f\n",kk,ix_front2[kk],iz_front2[kk],ttb_front2[kk]);
      }
*/
    }
    sort(ix_front2,iz_front2,ttb_front2,front_size2);
/*
      printf("Front2:\n");
      for( int kk=0;kk<front_size2;kk++ ) {
        printf("%d:(%d,%d):%f\n",kk,ix_front2[kk],iz_front2[kk],ttb_front2[kk]);      }
      printf("Front1:\n");
      for( int kk=0;kk<front_size;kk++ ) {
        printf("%d:(%d,%d):%f\n",kk,ix_front[kk],iz_front[kk],ttb_front[kk]);      }
*/
/*
    update(ix_front2,iz_front2,ttb_front2,front_size2,ttb,nz);
*/
    merge(ix_front,iz_front,ttb_front,front_big,&front_size,
          ix_front2,iz_front2,ttb_front2,front_size2);
/*
     printf("After Merge:\n");
     for( int kk=0;kk<front_size;kk++ ) {
        printf("%d:(%d,%d):%f---%f:%f\n",kk,ix_front[kk],iz_front[kk],
           ttb_front[kk],tmin,tmin+dt);
     }
     if( jj > 20000 ) exit(4);
*/
    nt=front_size;
/*
    if( (nt==1) && (ix_front[0] == 643) ) {
         kkk=kkk+1;
         if(kkk>3)exit(4);
    }
*/
  } while( nt>= 1 );
  corners(ttb,slow,hh,nx,nz);
  return;      
}

void corners(double *ttb, double *slow,
             double hh, int nx, int nz)
{
  double s1,sq;
  s1=(slow[0]+slow[nz+1]+slow[nz]+slow[1])/4;
  sq=2*hh*s1*s1-pow(ttb[nz]-ttb[1],2);
  ttb[0]=ttb[nz+1]+sqrt(sq);
  s1=(slow[(nx-1)*nz]+slow[(nx-2)*nz+1]
         +slow[(nx-2)*nz]+slow[(nx-1)*nz+1])/4;
  sq=2*hh*s1*s1-pow(ttb[(nx-2)*nz]-ttb[(nx-1)*nz+1],2);
  ttb[(nx-1)*nz]=ttb[(nx-2)*nz+1]+sqrt(sq);
  s1=(slow[nx*nz-1]+slow[(nx-2)*nz+nz-2]
         +slow[(nx-1)*nz+nz-2]+slow[(nx-2)*nz+nz-1])/4;
  sq=2*hh*s1*s1-pow(ttb[(nx-1)*nz+nz-2]-ttb[(nx-2)*nz+nz-1],2);
  ttb[(nx-1)*nz+nz-1]=ttb[(nx-2)*nz+nz-2]+sqrt(sq);
  s1=(slow[nz-1]+slow[nz+nz-2]
         +slow[nz-2]+slow[nz+nz-1])/4;
  sq=2*hh*s1*s1-pow(ttb[nz-2]-ttb[nz+nz-1],2);
  ttb[nz-1]=ttb[nz+nz-2]+sqrt(sq);
  return;
}

void eband(int *ix_front2, int *iz_front2,
           double *ttb_front2, int *front_size2,
           int n1, int k1, int nx, int nz, 
           double hh, double *slow, double *ttb)
{
  int itrat=0;
  int ind=n1*nz+k1;
  if( ttb[ind+1] < 0.0 ) {  		/* (n1,k1+1) */
    itrat=itrat+8;
  }
  if( ttb[ind-1] < 0.0 ) { 		/* (n1,k1-1) */
    itrat=itrat+4;
  }
  if( ttb[ind+nz] < 0.0 ) { 		/* (n1+1,k1) */
    itrat=itrat+2;
  }
  if( ttb[ind-nz] < 0.0 ) {    		/* (n1-1,k1) */
    itrat=itrat+1;
  }
/*
  if( (n1==62) && ( (k1<=23)&&(k1>=21) ) ) {
    printf("(%d,%d):%d\n",n1,k1,itrat);
  }
*/
  double s1,sq,trt;
  int m1,m2,m3;
  switch( itrat ) {
    case 1:
      if( (ttb[ind-nz-1]>-1.0)
           ) {         /* (n1-1,k1-1) */
        s1=(slow[ind-1]+slow[ind-nz-1]
                   +slow[ind]+slow[ind-nz])/4.0;
        sq=hh*s1*s1*2.0-pow((ttb[ind-nz-1]-ttb[ind]),2);
        if( sq<0.0 ) {
          sq=2*s1*s1*hh;
        }
        trt=ttb[ind-1]+sqrt(sq);
        trt=max(trt,ttb[ind]);
      } else if( (ttb[ind-nz+1]>-1.0)
         ) {
        s1=(slow[ind+1]+slow[ind-nz]
                  +slow[ind-nz+1]+slow[ind])/4.0;
        sq=hh*s1*s1*2.0-pow((ttb[ind-nz+1]-ttb[ind]),2);
        if( sq<0.0 ) {
          sq=2*s1*s1*hh;
        }
        trt=ttb[ind+1]+sqrt(sq);
        trt=max(trt,ttb[ind]);
      } else {
        s1=(slow[ind-nz]+(slow[ind]
                  +slow[ind+1]+slow[ind-1])/3)/2;
        sq=hh*s1*s1-0.25*pow((ttb[ind+1]-ttb[ind-1]),2);
        if( sq<0.0 ) {
          sq=s1*s1*hh;
        }
        trt=ttb[ind]+sqrt(sq);
      }
      if( (n1-1) > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
      } 
      ttb[ind-nz]=trt;
      break;
    case 2:
      if( (ttb[ind+nz-1]>-1.0) ) {
        s1=(slow[ind-1]+slow[ind+nz]
                   +slow[ind+nz-1]+slow[ind])/4;
        sq=hh*s1*s1*2-pow(ttb[ind+nz-1]-ttb[ind],2);
        if( sq<0.0 ) {
          sq=2*s1*s1*hh;
        }
        trt=ttb[ind-1]+sqrt(sq);
        trt=max(trt,ttb[ind]);
      } else if( (ttb[ind+nz+1]>-1.0) ) {
        s1=(slow[ind+1]+slow[ind+nz]
                   +slow[ind]+slow[ind+nz+1])/4;
        sq=hh*s1*s1*2-pow(ttb[ind+nz+1]-ttb[ind],2);
        if( sq<0.0 ) {
          sq=2*s1*s1*hh;
        }
        trt=ttb[ind+1]+sqrt(sq);
        trt=max(trt,ttb[ind]);
      } else {
        s1=(slow[ind+nz]+(slow[ind]
                    +slow[ind+1]+slow[ind-1])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
        if( sq<0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
      }
      if( n1+1 < nx-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
      } 
      ttb[ind+nz]=trt;
      break;
    case 3:
      s1=(slow[ind-nz]+(slow[ind]
                  +slow[ind+1]+slow[ind-1])/3)/2;
      sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
      if( sq<0.0 ) {
        sq=hh*s1*s1;
      }
      trt=ttb[ind]+sqrt(sq); 
      if( n1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
      } 
      ttb[ind-nz]=trt;
      s1=(slow[ind+nz]+(slow[ind]
                      +slow[ind+1]+slow[ind-1])/3)/2;
      sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
      if( sq<0.0 ) {
        sq=hh*s1*s1;
      }
      trt=ttb[ind]+sqrt(sq);
      if( n1+1 < nx-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
      } 
      ttb[ind+nz]=trt;
      break;
    case 4:
      if( (ttb[ind-nz-1]>-1.0)
           ) {
        s1=(slow[ind-1]+slow[ind-nz]
                  +slow[ind]+slow[ind-nz-1])/4;
        sq=hh*s1*s1*2-pow(ttb[ind]-ttb[ind-nz-1],2);
        if( sq<0.0 ) {
          sq=hh*s1*s1*2;
        }
        trt=ttb[ind-nz]+sqrt(sq);
        trt=max(trt,ttb[ind]);
      } else if( (ttb[ind+nz-1]>-1.0)
          ) {
        s1=(slow[ind-1]+slow[ind+nz]
                   +slow[ind]+slow[ind+nz-1])/4;
        sq=hh*s1*s1*2-pow(ttb[ind]-ttb[ind+nz-1],2);
        if( sq<0.0 ) {
          sq=hh*s1*s1*2;
        }
        trt=ttb[ind+nz]+sqrt(sq);
        trt=max(trt,ttb[ind]);
      } else {
        s1=(slow[ind-1]+(slow[ind]
                   +slow[ind-nz]+slow[ind+nz])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind-nz]-ttb[ind+nz],2);
        if( sq<0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
      }
      if( k1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
      } 
      ttb[ind-1]=trt;
      break;
    case 5:
      m1=0;m2=0;
      if( ttb[ind+nz-1] > -1.0 ) {
        s1=(slow[ind+nz]+slow[ind-1]
               +slow[ind+nz-1]+slow[ind])/4;
        sq=2*hh*s1*s1-pow(ttb[ind+nz-1]-ttb[ind],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind+nz]+sqrt(sq);
        trt=max(trt,ttb[ind]);
        if( k1-1 > 0 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
        } 
        ttb[ind-1]=trt;
	m1=1;
      }
      if( ttb[ind-nz+1] > -1.0 ) {
        s1=(slow[ind-nz]+slow[ind+1]
                +slow[ind]+slow[ind-nz+1])/4;
        sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz+1],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind+1]+sqrt(sq);
        trt=max(trt,ttb[ind]);
        if( n1-1 > 0 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
        } 
        ttb[ind-nz]=trt;
        m2=1;
      }
      if( (m1==1) && (m2==1) ) {
        break;
      }
      if( (m1==0) && (m2==1) ) {
        s1=(slow[ind-1]+(slow[ind]
                +slow[ind+nz]+slow[ind-nz])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
        if(k1-1 > 0 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
        } 
        ttb[ind-1]=trt;
	break;
      }
      if( (m1==1) && (m2==0) ) {
        s1=(slow[ind-nz]+(slow[ind]
                 +slow[ind+1]+slow[ind-1])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
        if( n1-1 > 0 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
        } 
        ttb[ind-nz]=trt;
	break;
      }
      if( (m1==0) && (m2==0) ) {
        if( ttb[ind+2*nz] > -1.0 ) {
          s1=(slow[ind+nz-1]+(slow[ind+nz]
                  +slow[ind]+slow[ind+2*nz])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind]-ttb[ind+2*nz],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind+nz]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          ttb[ind+nz-1]=trt;
          if( ((n1+1) < nx-1) && ((k1-1) > 0) ) {
            add(ix_front2, iz_front2,
               ttb_front2, front_size2,
               n1+1, k1-1, trt);
          }
        }
        if( ttb[ind+2] > -1.0 ) {
          s1=(slow[ind-nz+1]+(slow[ind+1]
                 +slow[ind]+slow[ind+2])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind]-ttb[ind+2],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind+1]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          ttb[ind-nz+1]=trt;
          if( ((n1-1)>0) && ((k1+1)<(nz-1)) ) {
            add(ix_front2, iz_front2,
               ttb_front2, front_size2,
               n1-1, k1+1, trt);
          }
        }
        if( ttb[ind+nz-1] > -1.0 ) {
          s1=(slow[ind+nz]+slow[ind-1]
                 +slow[ind+nz-1]+slow[ind])/4;
          sq=2*hh*s1*s1-pow(ttb[ind+nz-1]-ttb[ind],2);
          if( sq < 0.0 ) {
            sq=2*hh*s1*s1;
          }
          trt=ttb[ind+nz]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          if( k1-1 > 0 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1, k1-1, trt);
          } 
          ttb[ind-1]=trt;
          m1=1;
        }
        if( ttb[ind-nz+1] > -1.0 ) {
          s1=(slow[ind-nz]+slow[ind+1]
                  +slow[ind]+slow[ind-nz+1])/4;
          sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz+1],2);
          if( sq < 0.0 ) {
            sq=2*hh*s1*s1;
          }
          trt=ttb[ind+1]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          if( n1-1 > 0 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1-1, k1, trt);
          } 
          ttb[ind-nz]=trt;
          m2=1;
        }
        if( (m1==1) && (m2==1) ) {
          break;
        }
        if( (m1==0) && (m2==1) ) {
          s1=(slow[ind-1]+(slow[ind]
                  +slow[ind+nz]+slow[ind-nz])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind]+sqrt(sq);
          if(k1-1 > 0 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1, k1-1, trt);
          } 
          ttb[ind-1]=trt;
          break;
        }
        if( (m1==1) && (m2==0) ) {
          s1=(slow[ind-nz]+(slow[ind]
                   +slow[ind+1]+slow[ind-1])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind]+sqrt(sq);
          if( n1-1 > 0 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1-1, k1, trt);
          }
          ttb[ind-nz]=trt;
          break;
        }
      }
      break;
    case 6:
      m1=0;
      m2=0;
      if( ttb[ind+nz+1] > -1.0 ) {
        s1=(slow[ind+1]+slow[ind+nz]
               +slow[ind+nz+1]+slow[ind])/4;
        sq=2*hh*s1*s1-pow(ttb[ind+nz+1]-ttb[ind],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind+1]+sqrt(sq);
        trt=max(trt,ttb[ind]);
        if( (n1+1) < nx-1 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
        } 
        ttb[ind+nz]=trt;
	m1=1;
      }
      if( ttb[ind-nz-1] > -1.0 ) {
        s1=(slow[ind-nz]+slow[ind-1]
                +slow[ind]+slow[ind-nz-1])/4;
        sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz-1],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind-nz]+sqrt(sq);
        trt=max(trt,ttb[ind]);
        if( k1-1 > 0 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
        }
        ttb[ind-1]=trt;
        m2=1;
      }
      if( (m1==1) && (m2==1) ) {
        break;
      }
      if( (m1==0) && (m2==1) ) {
        s1=(slow[ind+nz]+(slow[ind]
                +slow[ind+1]+slow[ind-1])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
        if( (n1+1) < nx-1 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
        }
        ttb[ind+nz]=trt;
	break;
      }
      if( (m1==1) && (m2==0) ) {
        s1=(slow[ind-1]+(slow[ind]
                 +slow[ind+nz]+slow[ind-nz])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
        if( k1-1 > 0 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
        } 
        ttb[ind-1]=trt;
	break;
      }
      if( (m1==0) && (m2==0) ) {
        if( ttb[ind+2] > -1.0 ) {
          s1=(slow[ind+nz+1]+(slow[ind+1]
                  +slow[ind]+slow[ind+2])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind]-ttb[ind+2],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind+1]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          ttb[ind+nz+1]=trt;
          if( ((n1+1) < nx-1) && ((k1+1) < nz-1) ) {
            add(ix_front2, iz_front2,
               ttb_front2, front_size2,
               n1+1, k1+1, trt);
          }
        }
        if( ttb[ind-2*nz] > -1.0 ) {
          s1=(slow[ind-nz-1]+(slow[ind-nz]
                 +slow[ind]+slow[ind-2*nz])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind]-ttb[ind-2*nz],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind-nz]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          ttb[ind-nz-1]=trt;
          if( ((n1-1)>0) && ((k1-1)>0) ) {
            add(ix_front2, iz_front2,
               ttb_front2, front_size2,
               n1-1, k1-1, trt);
          }
        }
        if( ttb[ind+nz+1] > -1.0 ) {
          s1=(slow[ind+1]+slow[ind+nz]
                 +slow[ind+nz+1]+slow[ind])/4;
          sq=2*hh*s1*s1-pow(ttb[ind+nz+1]-ttb[ind],2);
          if( sq < 0.0 ) {
            sq=2*hh*s1*s1;
          }
          trt=ttb[ind+1]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          if( (n1+1)<nx-1 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1+1, k1, trt);
          } 
          ttb[ind+nz]=trt;
          m1=1;
        }
        if( ttb[ind-nz-1] > -1.0 ) {
          s1=(slow[ind-nz]+slow[ind-1]
                  +slow[ind]+slow[ind-nz-1])/4;
          sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz-1],2);
          if( sq < 0.0 ) {
            sq=2*hh*s1*s1;
          }
          trt=ttb[ind-nz]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          if( k1-1 > 0 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1, k1-1, trt);
          } 
          ttb[ind-1]=trt;
          m2=1;
        }
        if( (m1==1) && (m2==1) ) {
          break;
        }
        if( (m1==0) && (m2==1) ) {
          s1=(slow[ind+nz]+(slow[ind]
                  +slow[ind+1]+slow[ind-1])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind]+sqrt(sq);
          if( (n1+1) < nx-1 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1+1, k1, trt);
          } 
          ttb[ind+nz]=trt;
          break;
        }
        if( (m1==1) && (m2==0) ) {
          s1=(slow[ind-1]+(slow[ind]
                   +slow[ind+nz]+slow[ind-nz])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind]+sqrt(sq);
          if( k1-1 > 0 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1, k1-1, trt);
          } 
          ttb[ind-1]=trt;
          break;
        }
      }
      break;
    case 7:
      s1=(slow[ind-nz]+slow[ind+1]
             +slow[ind]+slow[ind-nz+1])/4;
      sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz+1],2);
      if( sq < 0.0 ) {
        sq=2*hh*s1*s1;
      }
      trt=ttb[ind+1]+sqrt(sq);
      trt=max(trt,ttb[ind]);
      if( n1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
      } 
      ttb[ind-nz]=trt;
      s1=(slow[ind+1]+slow[ind+nz]
              +slow[ind]+slow[ind+nz+1])/4;
      sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind+nz+1],2);
      if( sq < 0.0 ) {
        sq=2*hh*s1*s1;
      }
      trt=ttb[ind+1]+sqrt(sq);
      trt=max(trt,ttb[ind]);
      if( n1+1 < nx-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
      } 
      ttb[ind+nz]=trt;
      s1=(slow[ind-1]+(slow[ind]
              +slow[ind+nz]+slow[ind-nz])/3)/2;
      sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
      if( sq < 0.0 ) {
        sq=hh*s1*s1;
      }
      trt=ttb[ind]+sqrt(sq);
      if( k1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
      } 
      ttb[ind-1]=trt;
      break;
    case 8:
      if( (ttb[ind-nz+1]>-1.0)
           ) {
        s1=(slow[ind+1]+slow[ind-nz]
             +slow[ind]+slow[ind-nz+1])/4;
        sq=hh*s1*s1*2-pow(ttb[ind]-ttb[ind-nz+1],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
	trt=ttb[ind-nz]+sqrt(sq);
        trt=max(trt,ttb[ind]);
      } else if( (ttb[ind+nz+1]>-1.0)
           ) {
        s1=(slow[ind+nz]+slow[ind+1]
              +slow[ind]+slow[ind+nz+1])/4;
        sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind+nz+1],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind+nz]+sqrt(sq);
        trt=max(trt,ttb[ind]);
      } else {
        s1=(slow[ind+1]+(slow[ind]
               +slow[ind-nz]+slow[ind+nz])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind-nz]-ttb[ind+nz],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
      }
      if( k1+1 < nz-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
      } 
      ttb[ind+1]=trt;
      break;
    case 9:
      m1=0;m2=0;
      if( ttb[ind+nz+1] > -1.0 ) {
        s1=(slow[ind+1]+slow[ind+nz]
               +slow[ind+nz+1]+slow[ind])/4;
        sq=2*hh*s1*s1-pow(ttb[ind+nz+1]-ttb[ind],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind+nz]+sqrt(sq);
        trt=max(trt,ttb[ind]);
        if( (k1+1)<nz-1  ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
        } 
        ttb[ind+1]=trt;
	m1=1;
      }
      if( ttb[ind-nz-1] > -1.0 ) {
        s1=(slow[ind-nz]+slow[ind-1]
                +slow[ind]+slow[ind-nz-1])/4;
        sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz-1],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind-1]+sqrt(sq);
        trt=max(trt,ttb[ind]);
        if( n1-1 > 0 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
        } 
        ttb[ind-nz]=trt;
        m2=1;
      }
      if( (m1==1) && (m2==1) ) {
        break;
      }
      if( (m1==0) && (m2==1) ) {
        s1=(slow[ind+1]+(slow[ind]
                +slow[ind+nz]+slow[ind-nz])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
        if( (k1+1) < nz-1 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
        } 
        ttb[ind+1]=trt;
	break;
      }
      if( (m1==1) && (m2==0) ) {
        s1=(slow[ind-nz]+(slow[ind]
                 +slow[ind+1]+slow[ind-1])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
        if( n1-1 > 0 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
        } 
        ttb[ind-nz]=trt;
	break;
      }
      if( (m1==0) && (m2==0) ) {
        if( ttb[ind+2*nz] > -1.0 ) {
          s1=(slow[ind+nz+1]+(slow[ind+nz]
                  +slow[ind]+slow[ind+2*nz])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind]-ttb[ind+2*nz],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind+nz]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          ttb[ind+nz+1]=trt;
          if( ((n1+1) < nx-1) && ((k1+1)<nz-1) ) {
            add(ix_front2, iz_front2,
               ttb_front2, front_size2,
               n1+1, k1+1, trt);
          }
        }
        if( ttb[ind-2] > -1.0 ) {
          s1=(slow[ind-nz-1]+(slow[ind-1]
                 +slow[ind]+slow[ind-2])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind]-ttb[ind-2],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind-1]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          ttb[ind-nz-1]=trt;
          if( ((n1-1)>0) && ((k1-1)>0) ) {
            add(ix_front2, iz_front2,
               ttb_front2, front_size2,
               n1-1, k1-1, trt);
          }
        }
        if( ttb[ind+nz+1] > -1.0 ) {
          s1=(slow[ind+1]+slow[ind+nz]
                 +slow[ind+nz+1]+slow[ind])/4;
          sq=2*hh*s1*s1-pow(ttb[ind+nz+1]-ttb[ind],2);
          if( sq < 0.0 ) {
            sq=2*hh*s1*s1;
          }
          trt=ttb[ind+nz]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          if( k1+1 < nz-1 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1, k1+1, trt);
          } 
          ttb[ind+1]=trt;
          m1=1;
        }
        if( ttb[ind-nz-1] > -1.0 ) {
          s1=(slow[ind-nz]+slow[ind-1]
                  +slow[ind]+slow[ind-nz-1])/4;
          sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz-1],2);
          if( sq < 0.0 ) {
            sq=2*hh*s1*s1;
          }
          trt=ttb[ind-1]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          if( n1-1 > 0 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1-1, k1, trt);
          } 
          ttb[ind-nz]=trt;
          m2=1;
        }
        if( (m1==1) && (m2==1) ) {
          break;
        }
        if( (m1==0) && (m2==1) ) {
          s1=(slow[ind+1]+(slow[ind]
                  +slow[ind+nz]+slow[ind-nz])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind]+sqrt(sq);
          if( k1+1 < nz-1 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1, k1+1, trt);
          } 
          ttb[ind+1]=trt;
          break;
        }
        if( (m1==1) && (m2==0) ) {
          s1=(slow[ind-nz]+(slow[ind]
                   +slow[ind-1]+slow[ind+1])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind]+sqrt(sq);
          if( n1-1 > 0 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1-1, k1, trt);
          } 
          ttb[ind-nz]=trt;
          break;
        }
      }
      break;
    case 10:
      m1=0;m2=0;
      if( ttb[ind+nz-1] > -1.0 ) {
        s1=(slow[ind+nz]+slow[ind-1]
               +slow[ind+nz-1]+slow[ind])/4;
        sq=2*hh*s1*s1-pow(ttb[ind+nz-1]-ttb[ind],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind-1]+sqrt(sq);
        trt=max(trt,ttb[ind]);
        if( n1+1 < nx-1 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
        } 
        ttb[ind+nz]=trt;
	m1=1;
      }
      if( ttb[ind-nz+1] > -1.0 ) {
        s1=(slow[ind-nz]+slow[ind+1]
                +slow[ind]+slow[ind-nz+1])/4;
        sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz+1],2);
        if( sq < 0.0 ) {
          sq=2*hh*s1*s1;
        }
        trt=ttb[ind-nz]+sqrt(sq);
        trt=max(trt,ttb[ind]);
        if( k1+1 < nz-1 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
        } 
        ttb[ind+1]=trt;
        m2=1;
      }
      if( (m1==1) && (m2==1) ) {
        break;
      }
      if( (m1==0) && (m2==1) ) {
        s1=(slow[ind+nz]+(slow[ind]
                +slow[ind+1]+slow[ind-1])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
        if( n1+1 < nx-1 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
        } 
        ttb[ind+nz]=trt;
	break;
      }
      if( (m1==1) && (m2==0) ) {
        s1=(slow[ind+1]+(slow[ind]
                 +slow[ind+nz]+slow[ind-nz])/3)/2;
        sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
        if( sq < 0.0 ) {
          sq=hh*s1*s1;
        }
        trt=ttb[ind]+sqrt(sq);
        if( k1+1 < nz-1 ) {
          add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
        } 
        ttb[ind+1]=trt;
	break;
      }
      if( (m1==0) && (m2==0) ) {
        if( ttb[ind-2] > -1.0 ) {
          s1=(slow[ind+nz-1]+(slow[ind-1]
                  +slow[ind]+slow[ind-2])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind]-ttb[ind-2],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind-1]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          ttb[ind+nz-1]=trt;
          if( ((n1+1) < nx-1) && ((k1-1) > 0) ) {
            add(ix_front2, iz_front2,
               ttb_front2, front_size2,
               n1+1, k1-1, trt);
          }
        }
        if( ttb[ind-2*nz] > -1.0 ) {
          s1=(slow[ind-nz+1]+(slow[ind-nz]
                 +slow[ind]+slow[ind-2*nz])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind]-ttb[ind-2*nz],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind-nz]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          ttb[ind-nz+1]=trt;
          if( ((n1-1)>0) && ((k1+1)<(nz-1)) ) {
            add(ix_front2, iz_front2,
               ttb_front2, front_size2,
               n1-1, k1+1, trt);
          }
        }
        if( ttb[ind+nz-1] > -1.0 ) {
          s1=(slow[ind+nz]+slow[ind-1]
                 +slow[ind+nz-1]+slow[ind])/4;
          sq=2*hh*s1*s1-pow(ttb[ind+nz-1]-ttb[ind],2);
          if( sq < 0.0 ) {
            sq=2*hh*s1*s1;
          }
          trt=ttb[ind-1]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          if( n1+1 < nx-1 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1+1, k1, trt);
          } 
          ttb[ind+nz]=trt;
          m1=1;
        }
        if( ttb[ind-nz+1] > -1.0 ) {
          s1=(slow[ind-nz]+slow[ind+1]
                  +slow[ind]+slow[ind-nz+1])/4;
          sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz+1],2);
          if( sq < 0.0 ) {
            sq=2*hh*s1*s1;
          }
          trt=ttb[ind-nz]+sqrt(sq);
          trt=max(trt,ttb[ind]);
          if( k1+1 < nz-1 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1, k1+1, trt);
          } 
          ttb[ind+1]=trt;
          m2=1;
        }
        if( (m1==1) && (m2==1) ) {
          break;
        }
        if( (m1==0) && (m2==1) ) {
          s1=(slow[ind+nz]+(slow[ind]
                  +slow[ind+1]+slow[ind-1])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind+1]-ttb[ind-1],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind]+sqrt(sq);
          if( n1+1 < nx-1 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1+1, k1, trt);
          } 
          ttb[ind+nz]=trt;
          break;
        }
        if( (m1==1) && (m2==0) ) {
          s1=(slow[ind+1]+(slow[ind]
                   +slow[ind+nz]+slow[ind-nz])/3)/2;
          sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
          if( sq < 0.0 ) {
            sq=hh*s1*s1;
          }
          trt=ttb[ind]+sqrt(sq);
          if( k1+1 < nz-1 ) {
            add(ix_front2, iz_front2,
             ttb_front2, front_size2,
             n1, k1+1, trt);
          } 
          ttb[ind+1]=trt;
          break;
        }
      }
      break;
    case 11:
      s1=(slow[ind+nz]+slow[ind-1]
              +slow[ind]+slow[ind+nz-1])/4;
      sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind+nz-1],2);
      if( sq < 0.0 ) {
        sq=2*hh*s1*s1;
      }
      trt=ttb[ind-1]+sqrt(sq);
      trt=max(trt,ttb[ind]);
      if( n1+1 < nx-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
      } 
      ttb[ind+nz]=trt;
      s1=(slow[ind-nz]+slow[ind-1]
              +slow[ind]+slow[ind-nz-1])/4;
      sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind-nz-1],2);
      if( sq < 0.0 ) {
        sq=2*hh*s1*s1;
      }
      trt=ttb[ind-1]+sqrt(sq);
      trt=max(trt,ttb[ind]);
      if( n1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
      } 
      ttb[ind-nz]=trt;
      s1=(slow[ind+1]+(slow[ind]
             +slow[ind+nz]+slow[ind-nz])/3)/2;
      sq=hh*s1*s1-0.25*pow(ttb[ind+nz]-ttb[ind-nz],2);
      if( sq < 0.0 ) {
        sq=hh*s1*s1;
      }
      trt=ttb[ind]+sqrt(sq);
      if( k1+1 < nz-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
      } 
      ttb[ind+1]=trt;
      break;
    case 12:
      s1=(slow[ind-1]+(slow[ind]
               +slow[ind+nz]+slow[ind-nz])/3)/2;
      sq=hh*s1*s1-0.25*pow(ttb[ind-nz]-ttb[ind+nz],2);
      if( sq < 0.0 ) {
        sq=hh*s1*s1;
      }
      trt=ttb[ind]+sqrt(sq);
      if( k1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
      } 
      ttb[ind-1]=trt;
      s1=(slow[ind+1]+(slow[ind]
               +slow[ind+nz]+slow[ind-nz])/3)/2;
      sq=hh*s1*s1-0.25*pow(ttb[ind-nz]-ttb[ind+nz],2);
      if( sq < 0.0 ) {
        sq=hh*s1*s1;
      }
      trt=ttb[ind]+sqrt(sq);
      if( k1+1 < nz-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
      } 
      ttb[ind+1]=trt;
      break;
    case 13:
      s1=(slow[ind+nz]+slow[ind-1]
              +slow[ind]+slow[ind+nz-1])/4;
      sq=2*hh*s1*s1-pow(ttb[ind]-ttb[ind+nz-1],2);
      if( sq < 0.0 ) {
        sq=2*hh*s1*s1;
      }
      trt=ttb[ind+nz]+sqrt(sq);
      trt=max(trt,ttb[ind]);
      if( k1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
      } 
      ttb[ind-1]=trt;
      s1=(slow[ind+1]+slow[ind+nz]
             +slow[ind+nz+1]+slow[ind])/4;
      sq=2*hh*s1*s1-pow(ttb[ind+nz+1]-ttb[ind],2);
      if( sq < 0.0 ) {
        sq=2*hh*s1*s1;
      }
      trt=ttb[ind+nz]+sqrt(sq);
      trt=max(trt,ttb[ind]);
      if( k1+1 < nz-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
      } 
      ttb[ind+1]=trt;
      s1=(slow[ind-nz]+(slow[ind]
              +slow[ind-1]+slow[ind+1])/3)/2;
      sq=hh*s1*s1-0.25*pow(ttb[ind-1]-ttb[ind+1],2);
      if( sq < 0.0 ) {
        sq=hh*s1*s1;
      }
      trt=ttb[ind]+sqrt(sq);
      if( n1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1-1, k1, trt);
      } 
      ttb[ind-nz]=trt;
      break;
    case 14:
      s1=(slow[ind-nz]+slow[ind+1]
                 +slow[ind-nz+1]+slow[ind])/4;
      sq=2*hh*s1*s1-pow(ttb[ind-nz+1]-ttb[ind],2);
      if( sq < 0.0 ) {
        sq=2*hh*s1*s1;
      }
      trt=ttb[ind-nz]+sqrt(sq);
      trt=max(trt,ttb[ind]);
      if( k1+1 < nz-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1+1, trt);
      } 
      ttb[ind+1]=trt;
      s1=(slow[ind-nz]+slow[ind-1]
                  +slow[ind]+slow[ind-nz-1])/4;
      sq=2*hh*s1*s1-pow(ttb[ind-nz-1]-ttb[ind],2);
      if( sq < 0.0 ) {
        sq=2*hh*s1*s1;
      }
      trt=ttb[ind-nz]+sqrt(sq);
      trt=max(trt,ttb[ind]);
      if( k1-1 > 0 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1, k1-1, trt);
      } 
      ttb[ind-1]=trt;
      s1=(slow[ind+nz]+(slow[ind]
             +slow[ind-1]+slow[ind+1])/3)/2;
      sq=hh*s1*s1-0.25*pow(ttb[ind-1]-ttb[ind+1],2);
      if( sq<0.0 ) {
        sq=hh*s1*s1;
      }
      trt=ttb[ind]+sqrt(sq);
      if( n1+1 < nx-1 ) {
        add(ix_front2, iz_front2,
           ttb_front2, front_size2,
           n1+1, k1, trt);
      } 
      ttb[ind+nz]=trt;
  }
  return;
}

void add(int *ix_front2, int *iz_front2,
         double *ttb_front2, int *front_size2,
         int n1, int k1, double trt)
{
/*
  for( int i=0;i<(*front_size2);i++ ) {
    if( (n1==ix_front2[i]) && (k1==iz_front2[i]) ) {
      if( trt < ttb_front2[i] ) {
        ttb_front2[i]=trt;
      }
      return;
    }
  }
*/
  ix_front2[(*front_size2)]=n1;
  iz_front2[(*front_size2)]=k1;
  ttb_front2[(*front_size2)]=trt;
  (*front_size2)=(*front_size2)+1;
  return;
}

void merge(int *ix_front, int *iz_front, double *ttb_front,
           int front_big, int *front_size,
           int *ix_front2, int *iz_front2,
           double *ttb_front2, int front_size2)
{
  int nt=front_size2+(*front_size)-front_big;
  int ix[nt],iz[nt];
  int i;
  double tmp[nt];
  int index1=front_big,index2=0,index0=0;

  while( (index1 < (*front_size))
         && (index2 < front_size2) ) {
/*
    if( (ttb_front[index1] < ttb_front2[index2])
        || (fabs(ttb_front[index1]-ttb_front2[index2]) < DBL_MIN) ) {
*/
    if( LE(ttb_front[index1],ttb_front2[index2]) ) {
      ix[index0]=ix_front[index1];
      iz[index0]=iz_front[index1];
      tmp[index0]=ttb_front[index1];
      index1++;
    } else {
      ix[index0]=ix_front2[index2];
      iz[index0]=iz_front2[index2];
      tmp[index0]=ttb_front2[index2];
      index2++;
    }
    index0++;
  }
  while( index1 < (*front_size) ) {
    ix[index0]=ix_front[index1];
    iz[index0]=iz_front[index1];
    tmp[index0]=ttb_front[index1];
    index0++;index1++;
  }
  while( index2 < front_size2 ) {
    ix[index0]=ix_front2[index2];
    iz[index0]=iz_front2[index2];
    tmp[index0]=ttb_front2[index2];
    index0++;index2++;
  }
  *front_size=nt;
  for(i=0;i<nt;i++ ) {
    ix_front[i]=ix[i];
    iz_front[i]=iz[i];
    ttb_front[i]=tmp[i];
  }
  return;
}

void update(int *ix_front2, int *iz_front2,
            double *ttb_front2, int front_size2,
            double *ttb, int nz)
{
  int i;
  for( i=0;i<front_size2;i++ ) {
    int n1=ix_front2[i];
    int k1=iz_front2[i];
    ttb[n1*nz+k1]=ttb_front2[i];
  }
  return;
}

void sort(int *ix_front, int *iz_front,
          double *ttb, int front_size)
{
  if( front_size<= 1 ) return;
  double aln2i=1./0.69314718;
  double tiny=0.000001;
  int lognb2=(int)floor(log(front_size)*aln2i+tiny);
  int m=front_size, nn, j;

  for(nn=0;nn<lognb2;nn++ ) {
    m=m/2;
    int k=front_size-m;
    for(j=0;j<k;j++ ) {
      int i=j;
/*      while( ((ttb[i+m]<ttb[i]) || (fabs(ttb[i+m]-ttb[i])<0.00000001))*/
      while( ((ttb[i+m]<ttb[i]))
           && (i>=0) ) {
        int n1=ix_front[i];
        int n2=iz_front[i];
        double tmp=ttb[i];
        ix_front[i]=ix_front[i+m];
        iz_front[i]=iz_front[i+m];
        ttb[i]=ttb[i+m];
        ix_front[i+m]=n1;
        iz_front[i+m]=n2;
        ttb[i+m]=tmp;
        i=i-m;
      }
    }
  }
  return;
}

void sort_old(int *ix_front, int *iz_front, 
          int nz,
          double *ttb, int front_size) 
{
  double aln2i=1./0.69314718;
  double tiny=0.00001;
  int lognb2=(int)floor(log(front_size)*aln2i+tiny);
  int m=front_size, nn, j;

  for(nn=0;nn<lognb2;nn++ ) {
    m=m/2;
    int k=front_size-m;
    for(j=0;j<k;j++ ) {
      int i=j;
      while( (ttb[ix_front[i+m]*nz+iz_front[i+m]]
                <ttb[ix_front[i]*nz+iz_front[i]]) && (i>=0) ) {
        int n1=ix_front[i];
        int n2=iz_front[i];
        ix_front[i]=ix_front[i+m];
        iz_front[i]=iz_front[i+m];
        ix_front[i+m]=n1;
        iz_front[i+m]=n2;
        i=i-m;
      }
    }
  }
  return;
}

void source(double xs, double zs, 
            int ix0, int iz0, double dx,
            int nx, int nz,
            double *slow, double *ttb,
            int *ix_front, int *iz_front,
            double *ttb_front, int *front_size)
{
  int lx1,lx2,lz1,lz2;
  lx1=ix0-2;lx2=ix0+2;
  lz1=iz0-2;lz2=iz0+2;
  lx1=max2(0,lx1);lx2=min2(nx-1,lx2);
  lz1=max2(0,lz1);lz2=min2(nz-1,lz2);

  int i,j;
  int kx,kz;
  for( i=lx1;i<=lx2;i++ ) {
    for( j=lz1;j<=lz2;j++ ) {
      double s1=0.0,ss1=0.0;
      int kx1,kx2,kz1,kz2,k2;
      kx1=min2(i,ix0);kx2=max2(i,ix0);
      kz1=min2(j,iz0);kz2=max2(j,iz0);
      for( kx=kx1;kx<=kx2;kx++ ) {
        for( kz=kz1;kz<=kz2;kz++ ) {
	  int k1;
          k1=kx*nz+kz;
          s1=s1+slow[k1];
          ss1=ss1+1.0;
        }
      }
      s1=s1/ss1;
      k2=i*nz+j;
      ttb[k2]=sqrt(pow(i*dx-xs,2)+pow(j*dx-zs,2))*s1;
    }
  }
  *front_size = 0;
  for( kz=iz0-2;kz<=iz0+2;kz++ ) {
    if( (kz<=(nz-2)) && (kz>=1) ) {
      if( (ix0-2)>=1 ) {
        ix_front[(*front_size)]=ix0-2;
        iz_front[(*front_size)]=kz;
	ttb_front[(*front_size)]=ttb[(ix0-2)*nz+kz];
        *front_size=(*front_size)+1;
      }
      if( (ix0+2)<= (nx-2) ) {
        ix_front[(*front_size)]=ix0+2;
        iz_front[(*front_size)]=kz;
	ttb_front[(*front_size)]=ttb[(ix0+2)*nz+kz];
        *front_size=(*front_size)+1;
      }
    }
  }
  for( kx=ix0-1;kx<=(ix0+1);kx++  ) {
    if( (kx<=(nx-2)) && (kx>=1) ) {
      if( (iz0-2)>=1 ) {
        ix_front[(*front_size)]=kx;
        iz_front[(*front_size)]=iz0-2;
	ttb_front[(*front_size)]=ttb[kx*nz+(iz0-2)];
        *front_size=(*front_size)+1;
      }
      if( (iz0+2)<=(nz-2) ) {
        ix_front[(*front_size)]=kx;
        iz_front[(*front_size)]=iz0+2;
	ttb_front[(*front_size)]=ttb[kx*nz+(iz0+2)];
        *front_size=(*front_size)+1;
      }
    }
  }
  return;
}
      
int nround(double a) 
{
  int res;
  if( a>0 ) {
    res=(int)floor(a+0.5);
  } else {
    res=(int)-floor(-a+0.5);
  }
  return res;
}

void fgets2(char *str, int num, FILE *stream)
{
    fgets(str,num,stream);
    int i=strlen(str);
    if( str[i-1] == '\n' ) {
      str[i-1]='\0';
    }
    return;
}

double max1(double x, double y) 
{
    return (x>=y)?x:y;
}

double max(double x, double y) 
{
    return (x>=y)?x:y;
}

int max2(int x, int y) 
{
    return (x>=y)?x:y;
}

double min1(double x, double y)
{
    return (x<=y)?x:y;
}

double min(double x, double y)
{
    return (x<=y)?x:y;
}

int min2(int x, int y)
{
    return (x<=y)?x:y;
}

#endif
