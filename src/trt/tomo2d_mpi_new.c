#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "eik2d.h"
#include "par_sjm.h"
#include "n2s.h"
#include "ealloc.h"
#include "rwdouble.h"

#define LINE_LENGTH 256
#define TOMO_LONG_DECK_NAME "./TOMO_MPI_NEW.DECK"

#ifdef __cplusplus
extern "C" {
#endif

void get_assigned(int *id1, int *id2, int njob1, int njob2, int numprocs, int myrank);
void get_slow_shot(double *slow, float xmin_shot, float xmin, float dx,
                   int nx_shot, int nz_shot, float *vw, int nz0);
void get_surf_shot(int *nsurf_shot, float xmin_shot, float xmin,
                   float dx, int nx_shot, int *nsurf);
void get_grad(double *gw, float *pw, double *ra, int *pre,
              float xmin_shot, float xmin, float dx,
              int nx_shot, int nz_shot, int nz0);
void output(int nx0, int nz0, int ns, int *nr, int ng, int *nsurf,
            float *vw, double *dump, char *VfinalFile, 
            float *pw, char *RayFile,
            float **residual, char *ResidualFile, 
            int ite, float slow_min);
int step(FILE *fplog, int ns, int *nr, float **trt,
         float *xs, float *zs, float **xr, float **zr, 
         int nx_shot, int nz_shot,
         int nx0, int nz0, float dx, float xmin, float xmax,
         float offset_min, float offset_max,
         double a, double pur, int nite,
         double *gw, float *vw, double *dump, double *slow, 
         double *ts, double *alp,
         float slow_min, float slow_max);
void ray2d_tomo(double xr, double zr,
		double xs, double zs,
		double h,  double *t1,
	        int nx,	   int nz,
                int *ind_ray, int *is_visited,
		double trt_diff, 
	        double trt,
	        double rrr,
		int *nsurf, double *slow,
		double *ra, int *pre);
void pad_1(double *slow,int nx,int nz);
double getvalue2d(double x,double z, double h,
                  int nz, double *data);
void sort_1(float *surfx, float *surfz, int k);
void sort_2(int *ind, int kn);
float getvalue1d(float x,
                 int ndata, float *xcor, float *fx);
double getamax1(double *slow, int nx, int nz, int *nsurf);
double getmax(float *slow, int nx, int nz, int *nsurf);
void get_middle(int *ix00, int *iz00,
                float *pre, int nx, int nz);
void expand(int ix0, int iz0, double *gw, int *nsurf, int nx, int nz);
void expand_2(int ix0, int iz0, float *gw, float *pw, int *nsurf, int nx, int nz);
void expand_3(int *iz00, float *gw, float *pw, int *nsurf, int nx, int nz);
void fillHole_1(int ix0, int iz0, int nx, int nz, int k,
                double *g, int *nsurf);
void fillHole_2(int ix0, int iz0, int nx, int nz, int k,
                float *g, float *pw, int *nsurf);
void smooth2d(double *g, double *dump, int nx, int nz, 
              int irsx, int irsz, int *nsurf,
              int ix1_node, int ix2_node);

#ifdef __cplusplus
}
#endif

int main(int argc, char **argv)
{
   int numprocs,myrank;
   MPI_Status status;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

/* input file name */
   char *CoordFile,*VinitFile;
/* output file name */
   char *VfinalFile,*RmsFile;
   char *ResidualFile,*RayFile,*LogFile;

   CoordFile=ealloc1char(LINE_LENGTH);
   VinitFile=ealloc1char(LINE_LENGTH);
   VfinalFile=ealloc1char(LINE_LENGTH);
   RmsFile=ealloc1char(LINE_LENGTH);
   ResidualFile=ealloc1char(LINE_LENGTH);
   RayFile=ealloc1char(LINE_LENGTH);
   LogFile=ealloc1char(LINE_LENGTH);

   float xmin,xmax;
   float zmin,zmax;
   float offset_min,offset_max;
   float vmin,vmax;
   float slow_min;
   float slow_max;
   int input_initial_vel;
   int numsch,nite;
   float dx;
   int ns,ng;
   int *iter_sch,*irsx_sch,*irsz_sch;
   int *nr;
   float *xs,*zs,**xr,**zr,**trt;
   int nx0,nz0;

   float *data_float=ealloc1float(100);
   int *data_int=ealloc1int(100);
   register int i,j,k;
   FILE *fp,*fplog;
   int IsWaterBottomConstraint;
   int npWater=0;
   float *xWaterBottom,*zWaterBottom;

   if( myrank == 0 ) {
     char *DeckFile;
     DeckFile = (char*) malloc(100*sizeof(char));
     DeckFile = argv[1];
     GetParInit(DeckFile);

//      char *deck_name=ealloc1char(LINE_LENGTH);
//      strcpy(deck_name,TOMO_LONG_DECK_NAME);
//      GetParInit(deck_name);
      EGetParStr("LogFile",LogFile,0,"Log");
      fplog=efopen(LogFile,"w");

      EGetParStr("CoordFile",CoordFile,1,"");
      fprintf(fplog,"CoordFile=%s\n",CoordFile);
      int input_initial_vel;
      input_initial_vel=GetParStr("VinitFile",VinitFile);
      if( input_initial_vel ) {
         fprintf(fplog,"VinitFile=%s\n",VinitFile);
      }

      EGetParStr("VfinalFile",VfinalFile,0,"Vfinal");
      fprintf(fplog,"VfinalFile=%s\n",VfinalFile);
      EGetParStr("RmsFile",RmsFile,0,"RMS");
      fprintf(fplog,"RmsFile=%s\n",RmsFile);
      EGetParStr("ResidualFile",ResidualFile,0,"Res");
      fprintf(fplog,"ResidualFile=%s\n",ResidualFile);
      EGetParStr("RayFile",RayFile,0,"Ray");
      fprintf(fplog,"RayFile=%s\n",RayFile);

      EGetParFloat("XMIN",1,&xmin,0,0.0);
      EGetParFloat("XMAX",1,&xmax,0,xmin-1.0);
      EGetParFloat("ZMIN",1,&zmin,0,0.0);
      EGetParFloat("ZMAX",1,&zmax,1,zmin-1.0);
      EGetParFloat("OFFSET_MIN",1,&offset_min,0,0.0);
      EGetParFloat("OFFSET_MAX",1,&offset_max,0,offset_min-1.0);
      data_float[0]=zmin;
      data_float[1]=zmax;

      EGetParFloat("VMIN",1,&vmin,1,1200.0);
      EGetParFloat("VMAX",1,&vmax,1,5000.0);
      slow_min=1.0/vmax;
      slow_max=1.0/vmin;
      data_float[2]=vmin;
      data_float[3]=vmax;
      data_float[4]=slow_min;
      data_float[5]=slow_max;

      char *line=ealloc1char(LINE_LENGTH);
      float x1,x2,x3,x4,x5;
      char *WaterBottomFile=ealloc1char(LINE_LENGTH);
      EGetParInt("IsWaterBottomConstraint",1,&IsWaterBottomConstraint,0,0);
      if(IsWaterBottomConstraint) {
         EGetParStr("WaterBottomFile",WaterBottomFile,1,"WaterBottomFile");
         fp=efopen(WaterBottomFile,"r");
         while( !feof(fp) ) {
            strcpy(line,"");
            if( fgets(line,LINE_LENGTH,fp) && (strlen(line) > 3)  ) {
               npWater++;
            }
         }
         xWaterBottom=ealloc1float(npWater);
         zWaterBottom=ealloc1float(npWater);
	 rewind(fp);
         i=0;
         while( !feof(fp) ) {
            strcpy(line,"");
            if( fgets(line,LINE_LENGTH,fp) && (strlen(line) > 3)  ) {
               sscanf(line,"%f%f",&x1,&x2);
               xWaterBottom[i]=x1;
               zWaterBottom[i]=x2;
               i++;
            }
         }
         fclose(fp);
         free(WaterBottomFile);
      }

      EGetParInt("NumbeOfSmooth",1,&numsch,1,3);
      EGetParInt("Nite",1,&nite,1,4);
      EGetParFloat("Dx",1,&dx,1,12.5);
 
      data_int[0]=input_initial_vel;
      data_int[1]=nite;
      data_float[6]=dx;

      irsx_sch=ealloc1int(numsch);
      irsz_sch=ealloc1int(numsch);
      iter_sch=ealloc1int(numsch);
      int *xtmp=ealloc1int(3);
      for( i=0;i<numsch;i++ ) {
         sprintf(line,"ITER, IRSX, IRSZ %d",i+1);
         EGetParInt(line,3,xtmp,1,5);
         iter_sch[i]=xtmp[0];irsx_sch[i]=xtmp[1];irsz_sch[i]=xtmp[2];
         if( i > 0 ) {
            iter_sch[i]=iter_sch[i]+iter_sch[i-1];
         }
         data_int[3*i+2]=iter_sch[i];
         data_int[3*i+3]=irsx_sch[i];
         data_int[3*i+4]=irsz_sch[i];
      } 
      free(xtmp);

      free(argtbl);

      fp=efopen(CoordFile,"r");
      float xmin0,xmax0;
      int i1,i2,isin,max_shots=2000;
      int *indshot=ealloc1int(max_shots);
      ns=0;
      fscanf(fp,"%d%d%f%f%f%f%f",&i1,&i2,&x1,&x2,&x3,&x4,&x5);
      if((xmax<xmin) || (BTE(xmin,x1,xmax) && BTE(xmin,x3,xmax))) {
         indshot[ns]=i1;ns++;ng=i2;xmin0=MIN(x1,x3);xmax0=MAX(x1,x3);
      }
      while( !feof(fp) ) {
         strcpy(line,"");
         if( fgets(line,LINE_LENGTH,fp) ) {
            sscanf(line,"%d%d%f%f%f%f%f",&i1,&i2,&x1,&x2,&x3,&x4,&x5);
            if((xmax<xmin) || (BTE(xmin,x1,xmax) && BTE(xmin,x3,xmax))) {
               isin=0;
	       for(i=0;i<ns;i++) {
                  if( i1 == indshot[i] ) {
                     isin=1;
                     break;
                  }
               }
               if(isin == 0) {
                  if(ns == max_shots) {
                     max_shots=2*max_shots;
                     indshot=realloc1int(indshot,max_shots);
	          }
                  indshot[ns]=i1;
                  ns++;
               }
               ng=MAX(ng,i2);
               xmin0=MIN(xmin0,x1);xmin0=MIN(xmin0,x3);
               xmax0=MAX(xmax0,x1);xmax0=MAX(xmax0,x3);
            }
         }
      }
      data_int[3*numsch+2]=ns;
      data_int[3*numsch+3]=ng;

      nr=ealloc1int(ns);
      xs=ealloc1float(ns);
      zs=ealloc1float(ns);
      xr=ealloc2float(ng,ns);
      zr=ealloc2float(ng,ns);
      trt=ealloc2float(ng,ns);
      rewind(fp);
      for( i=0;i<ns;i++ ) {
         nr[i] = 0;
         for( j=0;j<ng;j++ ) {
            xr[i][j]=xmin0-1000.0;
            trt[i][j]=-1.0;
         }
      }
      sscanf(line,"%d%d%f%f%f%f%f",&i1,&i2,&x1,&x2,&x3,&x4,&x5);
      if((xmax<xmin) || (BTE(xmin,x1,xmax) && BTE(xmin,x3,xmax))) {
         nr[0]=i2;xs[0]=x1;zs[0]=x2;
         xr[0][i2-1]=x3;zr[0][i2-1]=x4;trt[0][i2-1]=x5;
      }
      while( !feof(fp) ) {
         strcpy(line,"");
         if( fgets(line,LINE_LENGTH,fp) ) {
            sscanf(line,"%d%d%f%f%f%f%f",&i1,&i2,&x1,&x2,&x3,&x4,&x5);
            if((xmax<xmin) || (BTE(xmin,x1,xmax) && BTE(xmin,x3,xmax))) {
	       for(i=0;i<ns;i++) {
                  if( i1 == indshot[i] ) {
                     break;
                  }
               }
               nr[i]=MAX(nr[i],i2);
               xs[i]=x1;zs[i]=x2;
               xr[i][i2-1]=x3;zr[i][i2-1]=x4;trt[i][i2-1]=x5;
            }
         }
      }
      free(indshot);
      fclose(fp);
      free(CoordFile);free(line);

      if( xmax < xmin ) {
         xmax=xmax0;
         xmin=xmin0;
      }
      nx0=NINT((xmax-xmin)/dx)+1;
      nz0=NINT((zmax-zmin)/dx)+1;
      fprintf(fplog,"NX,NZ=%d,%d\n",nx0,nz0);
      data_int[3*numsch+4]=nx0;
      data_int[3*numsch+5]=nz0;
      data_float[7]=xmin;
      data_float[8]=xmax;
 
      float off1,off2;
      off1=-1.0;
      off2=20000.0;

      for( i=0;i<ns;i++ ) {
         for( j=0;j<nr[i];j++ ) {
            if(BTE(xmin,xs[i],xmax) && BTE(xmin,xr[i][j],xmax)) {
               off2=MAX(off2,ABS(xs[i]-xr[i][j]));
               off1=MIN(off1,ABS(xs[i]-xr[i][j]));
            }
         }
      }

      if( offset_max < offset_min ) {
	 offset_min=off1;
         offset_max=off2;
      }
      offset_max=MIN(offset_max,off2);
      offset_min=MAX(offset_min,off1);
      data_float[9]=offset_min;
      data_float[10]=offset_max;
      data_int[3*numsch+6]=IsWaterBottomConstraint;
      data_int[3*numsch+7]=npWater;

   }

   MPI_Bcast(data_float,11,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPI_Bcast(&numsch,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(data_int,3*numsch+8,MPI_INT,0,MPI_COMM_WORLD);
   
   zmin=data_float[0];
   zmax=data_float[1];
   vmin=data_float[2];
   vmax=data_float[3];
   slow_min=data_float[4];
   slow_max=data_float[5];
   dx=data_float[6];
   xmin=data_float[7];
   xmax=data_float[8];
   offset_min=data_float[9];
   offset_max=data_float[10];
   input_initial_vel=data_int[0];
   nite=data_int[1];
   
   irsx_sch=ealloc1int(numsch);
   irsz_sch=ealloc1int(numsch);
   iter_sch=ealloc1int(numsch);
   for( i=0;i<numsch;i++ ) {
      iter_sch[i]=data_int[3*i+2];
      irsx_sch[i]=data_int[3*i+3];
      irsz_sch[i]=data_int[3*i+4];
   }
   ns=data_int[3*numsch+2];
   ng=data_int[3*numsch+3];
   nx0=data_int[3*numsch+4];
   nz0=data_int[3*numsch+5];

   if (myrank == 0)
   {
     printf("numsch = %d\n", numsch);
     printf("nx = %d, nz = %d, ns = %d, ng = %d\n",
     nx0, nz0, ns, ng);
   }

   IsWaterBottomConstraint=data_int[3*numsch+6];
   npWater=data_int[3*numsch+7];
   if( IsWaterBottomConstraint ) {
      if( myrank != 0 ) {
         xWaterBottom=ealloc1float(npWater);
         zWaterBottom=ealloc1float(npWater);
      }
      MPI_Bcast(xWaterBottom,npWater,MPI_FLOAT,0,MPI_COMM_WORLD);
      MPI_Bcast(zWaterBottom,npWater,MPI_FLOAT,0,MPI_COMM_WORLD);
   }

   free(data_int);free(data_float);
 
   if( myrank == 0 ) {
      fprintf(fplog,"xmin,xmax=%f,%f\n",xmin,xmax);
      fprintf(fplog,"zmin,zmax=%f,%f\n",zmin,zmax);
      fprintf(fplog,"offset_min,offset_max=%f,%f\n",offset_min,offset_max);
      fprintf(fplog,"vmin,vmax=%f,%f\n",vmin,vmax);
      fprintf(fplog,"slow_min,slow_max=%f,%f\n",slow_min,slow_max);
      fprintf(fplog,"dx=%f\n",dx);
      fprintf(fplog,"input_initial_vel=%d\n",input_initial_vel);
      fprintf(fplog,"nite=%d\n",nite);
      fprintf(fplog,"numsch=%d\n",numsch);   
      for( i=0;i<numsch;i++ ) {
         fprintf(fplog,"ITER, IRSX, IRSZ %d=%d,%d,%d\n",i+1,iter_sch[i],irsx_sch[i],irsz_sch[i]);
      }
      fprintf(fplog,"ns,ng=%d,%d\n",ns,ng);
      fprintf(fplog,"nx0,nz0=%d,%d\n",nx0,nz0);
   }

   if( myrank == 0 ) {
      nr=realloc1int(nr,ns);
      xs=realloc1float(xs,ns);
      zs=realloc1float(zs,ns);
   } else {
      nr=ealloc1int(ns);
      xs=ealloc1float(ns);
      zs=ealloc1float(ns);
   }

   MPI_Bcast(nr,ns,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(xs,ns,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPI_Bcast(zs,ns,MPI_FLOAT,0,MPI_COMM_WORLD);
   
   float *data_dump1=ealloc1float(ns*ng);
   float *data_dump2=ealloc1float(ns*ng);
   float *data_dump3=ealloc1float(ns*ng);
   if( myrank == 0 ) {
      for( i=0;i<ns;i++ ) {
         for( j=0;j<nr[i];j++ ) {
            data_dump1[j*ns+i]=xr[i][j];
            data_dump2[j*ns+i]=zr[i][j];
            data_dump3[j*ns+i]=trt[i][j];
         }
      }
   }
   xr=ealloc2float(ng,ns);
   zr=ealloc2float(ng,ns);
   trt=ealloc2float(ng,ns);
   int ntrt=0;
   for( i=0;i<ns;i++ ) {
      ntrt=ntrt+nr[i];
   }
   MPI_Bcast(data_dump1,ntrt,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPI_Bcast(data_dump2,ntrt,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPI_Bcast(data_dump3,ntrt,MPI_FLOAT,0,MPI_COMM_WORLD);
   for( i=0;i<ns;i++ ) {
      for( j=0;j<nr[i];j++ ) {
         xr[i][j]=data_dump1[j*ns+i];
         zr[i][j]=data_dump2[j*ns+i];
         trt[i][j]=data_dump3[j*ns+i];
      }
   }
   free(data_dump1);free(data_dump2);free(data_dump3);

   int nx_shot=NINT(offset_max/dx)+11;
   int nz_shot=nz0+10;
   double *slow=ealloc1double(nx_shot*nz_shot);
   double *ra=ealloc1double(nx_shot*nz_shot);
   double *ts=ealloc1double(nx_shot*nz_shot);
   int *pre=ealloc1int(nx_shot*nz_shot);

   float *vw=ealloc1float(nx0*nz0);
   double *gw=ealloc1double(nx0*nz0);
   double *dump=ealloc1double(nx0*nz0);
   float *dumpf=ealloc1float(nx0*nz0);
   float *dumps=ealloc1float(ns);
   float *pw=ealloc1float(nx0*nz0);
   float *rms=ealloc1float(ns);
   float *RMS=ealloc1float(iter_sch[numsch-1]+1);
   float **residual=ealloc2float(ng,ns);
   
   if( input_initial_vel) {
      if( myrank == 0 ) {
         printf("Read initial model of size %dx%d from %s\n", nx0, nz0,VinitFile);
         fp=efopen(VinitFile,"rb");
         fread(vw,sizeof(float),nx0*nz0,fp);
         fclose(fp);
      }
      MPI_Bcast(vw,nx0*nz0,MPI_FLOAT,0,MPI_COMM_WORLD);
   } else {
      float tmp=(vmax-vmin)/(float)(nz0-1);
      float tmp2;
      for( i=0;i<nz0;i++ ) {
         tmp2=vmin+tmp*(float)i;
         for( j=0;j<nx0;j++ ) {
            vw[j*nz0+i]=tmp2;
         }
      }
   }
   free(VinitFile);

   if(myrank == 0)
   {
     printf("Output an initial velocity as file vel.dat\n");
     FILE *file = fopen("vel.dat","w");
     for (i=0;i<nx0*nz0;i++)
       fprintf(file, "%f\n", vw[i]);
     fclose(file);
   }

   for( i=0;i<nx0*nz0;i++ ) {
      vw[i]=1.0/vw[i];
      if (vw[i] <= 0.0) printf("Error\n");
   }

// Topography  
   float *surfx,*surfz;
   surfx=ealloc1float(ns*ng+ns);
   surfz=ealloc1float(ns*ng+ns);
   k=0;
   for( i=0;i<ns;i++ ) {
      surfx[k]=xs[i];
      surfz[k]=zs[i];
      k++;
      for( j=0;j<nr[i];j++ ) {
         surfx[k]=xr[i][j];
         surfz[k]=zr[i][j];
         k++;
      }
   }
   sort_1(surfx,surfz,k);
   int kkk=0;
   float zerodis=0.001*dx;
   for( i=1;i<k;i++ ) {
      if( ABS(surfx[kkk]-surfx[i]) > zerodis ) {
         kkk++;
         surfx[kkk]=surfx[i];
         surfz[kkk]=surfz[i];
      }
   }
   kkk++;
   int *nsurf=ealloc1int(nx0);
   for( i=0;i<nx0;i++ ) {
      float zz;
      zz=getvalue1d(xmin+(float)i*dx,kkk,surfx,surfz);
      nsurf[i]=NINT(zz/dx);
   }
   free(surfx);free(surfz);

   int *nWater=ealloc1int(nx0);
   if(IsWaterBottomConstraint) {
      for( i=0;i<nx0;i++ ) {
         float zz;
         zz=getvalue1d(xmin+(float)i*dx,npWater,xWaterBottom,zWaterBottom);
         nWater[i]=NINT(zz/dx);
      }
   }

   int *nsurf_shot=ealloc1int(nx_shot);

/* main loop */
   int ite,isch=0;
   float pad=5.0*dx;
   float ra0=3.0/sqrt((offset_max/dx+1.0)*0.5);	
   int irsx,irsz;
   double ccc;
   int *ind_ray=ealloc1int(10*(nx_shot+nz_shot));
   int *is_visited=ealloc1int(nx_shot*nz_shot);
   int *iz0_tmp=ealloc1int(nx0);
   int id1,id2;
   get_assigned(&id1,&id2,0,ns-1,numprocs,myrank);
   int ix1_node,ix2_node;
   get_assigned(&ix1_node,&ix2_node,0,nx0-1,numprocs,myrank);
   /*printf("myrank,id1,id2=%d,%d,%d\n",myrank,id1,id2);*/

   float water1=1.0/1490.0;
   float water2=1.0/1500.0;

   if (myrank == 0)
     printf("max iter = %d\n", iter_sch[numsch-1]);

   for( ite=0;ite<=iter_sch[numsch-1];ite++ )
   {
     if (myrank == 0) printf("iter = %d\n", ite+1);

      for( i=0;i<nx0*nz0;i++ ) {
         pw[i]=0;gw[i]=0;
      }
      for( i=0;i<ns;i++ ) {
         rms[i]=0.0;
      }
      if( IsWaterBottomConstraint ) {
         int idump;         
         for( i=0;i<nx0;i++ ) {
            for( j=0;j<nWater[i];j++ ) {
               idump=i*nz0+j;
               vw[idump]= MIN(vw[idump],water1);
               vw[idump]= MAX(vw[idump],water2);
            }
         }
      }
      for( i=id1;i<=id2;i++ ) {
         rms[i]=0.0;
         if( BTE(xmin,xs[i],xmax) ) {
            float xmin_shot=xs[i];
            float xmax_shot=xs[i];
            for( j=0;j<nr[i];j++ ) {
               if( BTE(xmin,xr[i][j],xmax) ) {
                  xmin_shot=MIN(xmin_shot,xr[i][j]);
                  xmax_shot=MAX(xmax_shot,xr[i][j]);
               }
            }
            int nx_shot0=NINT((xmax_shot-xmin_shot)/dx)+11;
            get_slow_shot(slow,xmin_shot,xmin,dx,nx_shot0,nz_shot,vw,nz0);
            get_surf_shot(nsurf_shot,xmin_shot,xmin,dx,nx_shot0,nsurf);
            for( j=0;j<nx_shot*nz_shot;j++ ) {
               ts[j]=99999.0;
               ra[j]=0.0;
               pre[j]=0;
            }

            eik2d(nx_shot0,nz_shot,ts,(double)dx,slow,
                   (double)(xs[i]-xmin_shot+pad),(double)(zs[i]+pad));

            for( j=0;j<nr[i];j++ ) {
               if( BTE(xmin,xr[i][j],xmax) && trt[i][j] > 0.0 ) {
                  float dis=sqrt((xr[i][j]-xs[i])*(xr[i][j]-xs[i])
                                +(zr[i][j]-zs[i])*(zr[i][j]-zs[i]));
                  if( (dis > zerodis) && BTE(offset_min,dis,offset_max) ) {
                     double tt;
                     tt=getvalue2d(xr[i][j]-xmin_shot+pad,zr[i][j]+pad,dx,nz_shot,ts);
                     residual[i][j]=tt-trt[i][j];
                     rms[i]=rms[i]+residual[i][j]*residual[i][j];
                     if( ite == iter_sch[numsch-1] ) {
                        ray2d_tomo(xr[i][j]-xmin_shot+pad,zr[i][j]+pad,
                                   xs[i]-xmin_shot+pad,zs[i]+pad,
                                   dx,ts,nx_shot0,nz_shot,
                                   ind_ray,is_visited,
                                   residual[i][j],tt,0.0,nsurf_shot,slow,
                                   ra,pre);
                     } else {
                        ray2d_tomo(xr[i][j]-xmin_shot+pad,zr[i][j]+pad,
                                   xs[i]-xmin_shot+pad,zs[i]+pad,
                                   dx,ts,nx_shot0,nz_shot,
                                   ind_ray,is_visited,
                                   residual[i][j],tt,ra0,nsurf_shot,slow,
                                   ra,pre);
                     }
                  } /* end if for offset check */
               } /* end if for receiver x-coord check */
            } /* end loop for receiver */
            get_grad(gw,pw,ra,pre,xmin_shot,xmin,dx,nx_shot0,nz_shot,nz0);
         } /* end if for shot x-coord check */
      } /* end loop for shot */

      for( i=0;i<ns;i++ ) {
         dumps[i]=0.0;
      }
      for( i=0;i<nx0*nz0;i++ ) {
         dump[i]=0.0;
         dumpf[i]=0.0;
      }

      MPI_Allreduce(rms,dumps,ns,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(gw,dump,nx0*nz0,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(pw,dumpf,nx0*nz0,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

      for( i=0;i<ns;i++ ) {
         rms[i]=dumps[i];
      }
      for( i=0;i<nx0*nz0;i++ ) {
         gw[i]=dump[i];
         pw[i]=dumpf[i];
      }

      ccc=0.0;
      for( i=0;i<ns;i++ ) {
         ccc=ccc+rms[i];
      }
      RMS[ite]=sqrt(ccc/(float)ntrt);
      if( myrank == 0 ) {
     /* printf("it=%d,rms=%f\n",ite,RMS[ite]);*/
         fprintf(fplog,"it=%d,rms=%f\n",ite,RMS[ite]);
         fflush(fplog);
      }

      if( ite == iter_sch[isch] ) {
         if( myrank == 0 ) {
            output(nx0,nz0,ns,nr,ng,nsurf,
                vw,dump,VfinalFile,pw,RayFile,
                residual,ResidualFile,
                isch,slow_min);
         }
      }
      
      if( ite < iter_sch[numsch-1] ) {
         for( i=0;i<nx0*nz0;i++ ) {
            gw[i]=gw[i]/(pw[i]+0.01);
         }
         for( i=0;i<nx0;i++ ) {
            int knz=nsurf[i]-1;
            for( j=nsurf[i];j<nz0;j++ ) {
               if( pw[i*nz0+j] > 0.5 ) {
                  knz=j;
               }
            }
            iz0_tmp[i]=knz;
         }
/*
         int ix00_tmp,iz00_tmp;
         ix00_tmp=nx0/2;iz00_tmp=(nsurf[ix00_tmp]+iz0_tmp[ix00_tmp])/2;
         expand(ix00_tmp,iz00_tmp,gw,nx0,nz0);
*/
/*
         for( i=0;i<nx0;i++ ) {
            for( j=MAX(1,nsurf[i]);j<nz0;j++ ) {
               if( pw[i*nz0+j] < 0.5 ) {
                  gw[i*nz0+j]=gw[i*nz0+j-1];
               }
            }
         }
*/
         if( ite < iter_sch[isch] ) {
            irsx=irsx_sch[isch];irsz=irsz_sch[isch];
         } else {
            isch++;
            irsx=irsx_sch[isch];irsz=irsz_sch[isch];
         }
         smooth2d(gw,dump,nx0,nz0,irsx,irsz,nsurf,ix1_node,ix2_node);
         for( i=0;i<nx0*nz0;i++ ) {
            gw[i]=0.0;
         }
         MPI_Allreduce(dump,gw,nx0*nz0,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
         double gmax=getamax1(gw,nx0,nz0,nsurf);
         double smax=getmax(vw,nx0,nz0,nsurf);
         double pur=smax/100.0/gmax;
         double alp;
         ccc=0.0;
	 int iiii;
         for( i=0;i<nite;i++ ) {
            iiii=(ite+i*(ns/nite))%ns;
            ccc=ccc+rms[iiii];
         }
   /* search a step */
         int nstop=0;
         double *aa=ealloc1double(2);
         double *al=ealloc1double(2);
         int ii;
         for( i=0;i<nx0*nz0;i++ ) {
            dump[i]=vw[i];
         }
         al[0]=5.0;al[1]=10.0;
         int doagain=1;
         while( doagain == 1 ) {
            for( ii=0;ii<2;ii++ ) {
               aa[ii]=0.0;
               for( i=0;i<nx0*nz0;i++ ) {
                  vw[i]=dump[i]-gw[i]*pur*al[ii];
                  vw[i]=MAX(vw[i],slow_min);
                  vw[i]=MIN(vw[i],slow_max);
               }
               if( IsWaterBottomConstraint ) {
                  int idump;         
                  for( i=0;i<nx0;i++ ) {
                     for( j=0;j<nWater[i];j++ ) {
                        idump=i*nz0+j;
                        vw[idump]= MIN(vw[idump],water1);
                        vw[idump]= MAX(vw[idump],water2);
                     }
                  }
               }
               for( iiii=0;iiii<nite;iiii++ ) {
	          i=(ite+iiii*(ns/nite))%ns;
                  if( (i>=id1) && (i<=id2) ) {
                     if( BTE(xmin,xs[i],xmax) ) {
                        float xmin_shot=xs[i],xmax_shot=xs[i];
                        for( j=0;j<nr[i];j++ ) {
                           if( BTE(xmin,xr[i][j],xmax) ) {
                              xmin_shot=MIN(xmin_shot,xr[i][j]);
                              xmax_shot=MAX(xmax_shot,xr[i][j]);
                           }
                        }
                        int nx_shot0=NINT((xmax_shot-xmin_shot)/dx)+11;
                        get_slow_shot(slow,xmin_shot,xmin,dx,nx_shot0,
                                      nz_shot,vw,nz0);
                        for( j=0;j<nx_shot*nz_shot;j++ ) {
                           ts[j]=99999.0;
                        }
                        eik2d(nx_shot0,nz_shot,ts,dx,slow,
                              xs[i]-xmin_shot+pad,zs[i]+pad);
                        for( j=0;j<nr[i];j++ ) {
                           if( BTE(xmin,xr[i][j],xmax) && trt[i][j] > 0.0 ) {
                              float dis=sqrt((xr[i][j]-xs[i])*(xr[i][j]-xs[i])
                                +(zr[i][j]-zs[i])*(zr[i][j]-zs[i]));
                              if( (dis > zerodis) && BTE(offset_min,dis,offset_max) ) {
                                 double tt;
                                 tt=getvalue2d(xr[i][j]-xmin_shot+pad,zr[i][j]+pad,dx,nz_shot,ts);
                                 aa[ii]=aa[ii]+(tt-trt[i][j])*(tt-trt[i][j]);
	                      }
                           } /* end check receiver x-coord */
                        } /* end loop receiver */
                     } /* end check shot x-coord */
                  } /* check shot range */
               } /* end loop shot */
               double rms_tmp=0.0;
               MPI_Allreduce(&aa[ii],&rms_tmp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
               aa[ii]=rms_tmp;
               if( myrank == 0 ) {
	          fprintf(fplog,"ii=%d,aa=%f,al=%f,a0=%f\n",ii+1,aa[ii],al[ii],ccc);
               }
            } /* end loop ii */
            double per=(ccc*al[0]*al[0]-aa[1]*al[0]*al[0]-ccc*al[1]*al[1]+
                     aa[0]*al[1]*al[1])
                    /(2.0*(ccc*al[0]-aa[1]*al[0]-ccc*al[1]+aa[0]*al[1]));
            per=MIN(per,al[1]);
            if( per < 0.05 ) {
               if( (aa[0]<ccc) && (aa[0]<aa[1]) ) {
                  per=al[0];
                  doagain=0;
               } else if( (aa[1]<ccc) && (aa[1]<aa[0]) ) {
                  per=al[1];
                  doagain=0;
               } else if( al[0]>0.3 ) {
                  al[0]=al[0]/2.0;
                  al[1]=al[1]/2.0;
               } else {
                  nstop=1;
                  per=0.0;
                  doagain=0;
               }
            } else {
               doagain=0;
            }
            if( doagain == 0 ) {
               if( per < 0.0005 ) {
                  nstop=1;
               } else {
                  alp=pur*per;
               }
               for( i=0;i<nx0*nz0;i++ ) {
                  vw[i]=dump[i];
               }
            }
         } /* end while */
         free(al);free(aa);
   
         nstop=0;
         if( nstop == 0 ) {
            for( i=0;i<nx0*nz0;i++ ) {
               vw[i]=vw[i]-alp*gw[i];
               vw[i]=MAX(slow_min,vw[i]);
               vw[i]=MIN(slow_max,vw[i]);
            }
            if( IsWaterBottomConstraint ) {
               int idump;         
               for( i=0;i<nx0;i++ ) {
                  for( j=0;j<nWater[i];j++ ) {
                     idump=i*nz0+j;
                     vw[idump]= MIN(vw[idump],water1);
                     vw[idump]= MAX(vw[idump],water2);
                  }
               }
            }
            for( i=0;i<nx0;i++ ) {
               for( j=iz0_tmp[i]+1;j<nz0;j++ ) {
                  vw[i*nz0+j]=vw[i*nz0+iz0_tmp[i]];
               }
            }
/*
            for( i=0;i<nx0;i++ ) {
               for( j=0;j<=iz0_tmp[i];j++ ) {
                  pw[i*nz0+j]=1.0;
               }
            }
            expand_3(iz0_tmp,vw,pw,nsurf,nx0,nz0);
            for( i=0;i<nx0*nz0;i++ ) {
               vw[i]=MIN(vw[i],slow_max);
               vw[i]=MAX(vw[i],slow_min);
	    }
*/
         } else {
            int itmp=iter_sch[isch]-ite-1;
            for( i=isch;i<numsch;i++ ) {
               iter_sch[i]=iter_sch[i]-itmp;
            }
         }
      }
   } /* end loop for iteration */

   if( myrank == 0 ) {
      fp=efopen(RmsFile,"w");
      for( i=0;i<=iter_sch[numsch-1];i++ ) {
         fprintf(fp,"%f\n",RMS[i]);
      }
      fclose(fp);
      fclose(fplog);
   }

   MPI_Finalize();
/*
   free(RmsFile);
   free(VfinalFile);free(RayFile);free(ResidualFile);free(LogFile);
   free(irsx_sch);free(irsz_sch);free(iter_sch);
   free(nr);free(xs);free(zs);free2float(xr);free2float(zr);free2float(trt);
   free(slow);free(ra);free(ts);free(pre);free(vw);free(gw);free(pw);
   free(dump);free(dumpf);free(dumps);free(rms);free(RMS); 
   free2float(residual);free(nsurf);free(nsurf_shot);
   free(ind_ray);free(is_visited);
   free(iz0_tmp);
*/
  return 0;
}

void get_slow_shot(double *slow, float xmin_shot, float xmin, float dx,
                   int nx_shot, int nz_shot, float *vw, int nz0)
{
  int ix0=NINT((xmin_shot-xmin)/dx);
  int i,j;
  for( i=5;i<nx_shot-5;i++ ) {
     for( j=5;j<nz_shot-5;j++ ) {
        slow[i*nz_shot+j]=vw[(ix0+i-5)*nz0+j-5];
     }
  }
  pad_1(slow,nx_shot,nz_shot);
  return;
}

void get_surf_shot(int *nsurf_shot, float xmin_shot, float xmin,
                   float dx, int nx_shot, int *nsurf)
{
  int ix0=NINT((xmin_shot-xmin)/dx);
  int i;
  for( i=5;i<nx_shot-5;i++ ) {
     nsurf_shot[i]=nsurf[ix0+i-5]+5;
  }
  for( i=0;i<5;i++ ) {
     nsurf_shot[i]=nsurf_shot[5];
  }
  for( i=nx_shot-5;i<nx_shot;i++ ) {
     nsurf_shot[i]=nsurf_shot[nx_shot-6];
  }
  return;
}

void get_grad(double *gw, float *pw, double *ra, int *pre,
              float xmin_shot, float xmin, float dx,
              int nx_shot, int nz_shot, int nz0)
{
  int ix0=NINT((xmin_shot-xmin)/dx);
  int i,j;
  for( i=5;i<nx_shot-5;i++ ) {
     for( j=5;j<nz_shot-5;j++ ) {
        gw[(ix0+i-5)*nz0+j-5]=gw[(ix0+i-5)*nz0+j-5]+ra[i*nz_shot+j];
        pw[(ix0+i-5)*nz0+j-5]=pw[(ix0+i-5)*nz0+j-5]+(float)pre[i*nz_shot+j];
     }
  }
  return;
}

void output(int nx0, int nz0, int ns, int *nr, int ng, int *nsurf,
            float *vw, double *dump, char *VfinalFile, 
            float *pw, char *RayFile,
            float **residual, char *ResidualFile, 
            int ite, float slow_min)
{
   FILE *fp;
   char *file_tmp=ealloc1char(LINE_LENGTH);
   strcpy(file_tmp,VfinalFile);
   strcat(file_tmp,n2s(ite+1)); 
   fp=efopen(file_tmp,"wb");
   int i,j;
   for( i=0;i<nx0*nz0;i++ ) {
      dump[i]=vw[i];
   }
/*
   for( i=0;i<nx0;i++ ) {
      for( j=nz0-1;j>=0;j-- ) {
         if( pw[i*nz0+j] < 0.5 ) {
            vw[i*nz0+j]=slow_min;
         } else {
            break;
         }
      }
   }
*/
   for( i=0;i<nx0;i++ ) {
      for( j=0;j<nsurf[i];j++ ) {
         vw[i*nz0+j]=0.0;
      }
      for( j=nsurf[i];j<nz0;j++ ) {
         vw[i*nz0+j]=1.0/vw[i*nz0+j];
      }
   }
   fwrite(vw,sizeof(float),nx0*nz0,fp);
   fclose(fp);
   for( i=0;i<nx0*nz0;i++ ) {
      vw[i]=dump[i];
   }

   strcpy(file_tmp,RayFile);
   strcat(file_tmp,n2s(ite+1));
   fp=efopen(file_tmp,"wb");
   fwrite(pw,sizeof(float),nx0*nz0,fp);
   fclose(fp);
   
/*
   strcpy(file_tmp,ResidualFile);
   strcat(file_tmp,n2s(ite+1)); 
   fp=efopen(file_tmp,"wb");
   for( i=0;i<ns;i++ ) {
      fwrite(&residual[i][0],sizeof(float),nr[i],fp);
   }
   fclose(fp);
*/

   free(file_tmp);
   return;
}

int step(FILE *fplog, int ns, int *nr, float **trt,
         float *xs, float *zs, float **xr, float **zr, 
         int nx_shot, int nz_shot,
         int nx0, int nz0, float dx, float xmin, float xmax,
         float offset_min, float offset_max,
         double a, double pur, int nite,
         double *gw, float *vw, double *dump, double *slow, 
         double *ts, double *alp,
         float slow_min, float slow_max)
{
  int nstop=0;
  double *aa=ealloc1double(2);
  double *al=ealloc1double(2);
  float pad=5.0*dx;
  float zerodis=0.001*dx;
  int i,j,ii;
  for( i=0;i<nx0*nz0;i++ ) {
     dump[i]=vw[i];
  }
  al[0]=5.0;al[1]=10.0;
  int doagain=1;
  while( doagain == 1 ) {
     aa[0]=0.0;aa[1]=0.0;
     for( ii=0;ii<2;ii++ ) {
        for( i=0;i<nx0*nz0;i++ ) {
           vw[i]=dump[i]-gw[i]*pur*al[ii];
           vw[i]=MAX(vw[i],slow_min);
           vw[i]=MIN(vw[i],slow_max);
        }
        for( i=0;i<ns;i+=ns/nite ) {
           if( BTE(xmin,xs[i],xmax) ) {
              float xmin_shot=xs[i],xmax_shot=xs[i];
              for( j=0;j<nr[i];j++ ) {
                 if( BTE(xmin,xr[i][j],xmax) ) {
                    xmin_shot=MIN(xmin_shot,xr[i][j]);
                    xmax_shot=MAX(xmax_shot,xr[i][j]);
                 }
              }
              int nx_shot0=NINT((xmax_shot-xmin_shot)/dx)+11;
              get_slow_shot(slow,xmin_shot,xmin,dx,nx_shot0,nz_shot,vw,nz0);
              for( j=0;j<nx_shot*nz_shot;j++ ) {
                 ts[j]=99999.0;
              }
              eik2d(nx_shot0,nz_shot,ts,dx,slow,xs[i]-xmin+pad,zs[i]+pad);
              for( j=0;j<nr[i];j++ ) {
                 if( BTE(xmin,xr[i][j],xmax) && trt[i][j] > 0.0 ) {
                    float dis=sqrt((xr[i][j]-xs[i])*(xr[i][j]-xs[i])
                                +(zr[i][j]-zs[i])*(zr[i][j]-zs[i]));
                    if( (dis > zerodis) && BTE(offset_min,dis,offset_max) ) {
                       double tt;
                       tt=getvalue2d(xr[i][j]-xmin+pad,zr[i][j]+pad,dx,nz_shot,ts);
                       aa[ii]=aa[ii]+(tt-trt[i][j])*(tt-trt[i][j]);
	            }
                 } /* end check receiver x-coord */
              } /* end loop receiver */
           } /* end check shot x-coord */
        } /* end loop shot */
	fprintf(fplog,"ii=%d,aa=%f,al=%f,a0=%f\n",ii+1,aa[ii],al[ii],a);
     } /* end loop ii */
     double per=(a*al[0]*al[0]-aa[1]*al[0]*al[0]-a*al[1]*al[1]+aa[0]*al[1]*al[1])
                  /(2.0*(a*al[0]-aa[1]*al[0]-a*al[1]+aa[0]*al[1]));
     per=MIN(per,al[1]);
     if( per < 0.05 ) {
        if( (aa[0]<a) && (aa[0]<aa[1]) ) {
           per=al[0];
           doagain=0;
        } else if( (aa[1]<a) && (aa[1]<aa[0]) ) {
           per=al[1];
           doagain=0;
        } else if( al[0]>0.3 ) {
           al[0]=al[0]/2.0;
           al[1]=al[1]/2.0;
        } else {
           nstop=1;
           for( i=0;i<nx0*nz0;i++ ) {
              slow[i]=dump[i];
           }
           free(al);free(aa);
           return nstop;
        }
     } else {
        doagain=0;
     }
     if( doagain == 0 ) {
        if( per < 0.0005 ) {
           nstop=1;
           for( i=0;i<nx0*nz0;i++ ) {
              slow[i]=dump[i];
           }
           free(al);free(aa);
           return nstop;
        }
        *alp=pur*per;
        for( i=0;i<nx0*nz0;i++ ) {
           slow[i]=dump[i];
        }
        free(al);free(aa);
        return 0;
     }
  } /* end while */
  free(al);free(aa);
  return 0;
}

void get_middle(int *ix00, int *iz00,
                float *pre, int nx, int nz)
{
   int isGet=0;
   int i,k,ix0,iz0;
   ix0=nx/2;iz0=ix0*nz;
   float tmp=-1.0;

   for( i=0;i<nz;i++ ) {
      if( pre[iz0+i] > tmp ) {
         ix0=i;
         tmp=pre[iz0+i];
      }
   }
   *ix00=nx/2;*iz00=ix0;
   return;
}

void expand(int ix0, int iz0, double *gw, int *nsurf, int nx, int nz)
{
   int i,j,k;
   int kmax=MAX(nx,nz);
   kmax=MAX(kmax,ix0);kmax=MAX(kmax,nx-1-ix0);
   kmax=MAX(kmax,iz0);kmax=MAX(kmax,nz-1-iz0);
   for(k=1;k<=kmax;k++) {
      if( iz0-k >= 0 ) {
         for(i=ix0;i<=MIN(nx-1,ix0+k);i++) {
            if( LE(ABS(gw[i*nz+iz0-k]),0.0) && ((iz0-k)>=nsurf[i]) ) {
               fillHole_1(i,iz0-k,nx,nz,5,gw,nsurf);
            }
         }
         for(i=ix0-1;i>=MAX(0,ix0-k);i--) {
            if( LE(ABS(gw[i*nz+iz0-k]),0.0) && ((iz0-k)>=nsurf[i]) ) {
               fillHole_1(i,iz0-k,nx,nz,5,gw,nsurf);
            }
         }
      }
      if( ix0+k < nx ) {
         for(i=iz0-1;i>=MAX(0,iz0-k);i--) {
            if( LE(ABS(gw[(ix0+k)*nz+i]),0.0) && (i>=nsurf[ix0+k]) ) {
               fillHole_1(ix0+k,i,nx,nz,5,gw,nsurf);
            }
         }
         for(i=iz0;i<=MIN(nz-1,iz0+k);i++) {
            if( LE(ABS(gw[(ix0+k)*nz+i]),0.0) && (i>=nsurf[ix0+k]) ) {
               fillHole_1(ix0+k,i,nx,nz,5,gw,nsurf);
            }
         }
      }
      if( iz0+k < nz ) {
         for(i=ix0;i<=MIN(nx-1,ix0+k);i++) {
            if( LE(ABS(gw[i*nz+iz0+k]),0.0) && ((iz0+k)>=nsurf[i]) ) {
               fillHole_1(i,iz0+k,nx,nz,5,gw,nsurf);
            }
         }
         for(i=ix0-1;i>=MAX(0,ix0-k);i--) {
            if( LE(ABS(gw[i*nz+iz0+k]),0.0) && ((iz0+k)>=nsurf[i]) ) {
               fillHole_1(i,iz0+k,nx,nz,5,gw,nsurf);
            }
         }
      }
      if( ix0-k >= 0 ) {
         for(i=iz0-1;i>=MAX(0,iz0-k);i--) {
            if( LE(ABS(gw[(ix0-k)*nz+i]),0.0) && (i>=nsurf[ix0-k]) ) {
               fillHole_1(ix0-k,i,nx,nz,5,gw,nsurf);
            }
         }
         for(i=iz0;i<=MIN(nz-1,iz0+k);i++) {
            if( LE(ABS(gw[(ix0-k)*nz+i]),0.0) && (i>=nsurf[ix0-k]) ) {
               fillHole_1(ix0-k,i,nx,nz,5,gw,nsurf);
            }
         }
      }
   }
   return;
}

void expand_2(int ix0, int iz0, float *gw, float *pw, int *nsurf, int nx, int nz)
{
   int i,j,k;
   int kmax=MAX(nx,nz);
   kmax=MAX(kmax,ix0);kmax=MAX(kmax,nx-1-ix0);
   kmax=MAX(kmax,iz0);kmax=MAX(kmax,nz-1-iz0);
   for(k=1;k<=kmax;k++) {
      if( iz0-k >= 0 ) {
         for(i=ix0;i<=MIN(nx-1,ix0+k);i++) {
            if( (pw[i*nz+iz0-k]<0.5) && ((iz0-k)>=nsurf[i]) ) {
               fillHole_2(i,iz0-k,nx,nz,5,gw,pw,nsurf);
            }
         }
         for(i=ix0-1;i>=MAX(0,ix0-k);i--) {
            if( (pw[i*nz+iz0-k]<0.5) && ((iz0-k)>=nsurf[i]) ) {
               fillHole_2(i,iz0-k,nx,nz,5,gw,pw,nsurf);
            }
         }
      }
      if( ix0+k < nx ) {
         for(i=iz0-1;i>=MAX(0,iz0-k);i--) {
            if( (pw[(ix0+k)*nz+i]<0.5) && (i>=nsurf[ix0+k]) ) {
               fillHole_2(ix0+k,i,nx,nz,5,gw,pw,nsurf);
            }
         }
         for(i=iz0;i<=MIN(nz-1,iz0+k);i++) {
            if( (pw[(ix0+k)*nz+i]<0.5) && (i>=nsurf[ix0+k]) ) {
               fillHole_2(ix0+k,i,nx,nz,5,gw,pw,nsurf);
            }
         }
      }
      if( iz0+k < nz ) {
         for(i=ix0;i<=MIN(nx-1,ix0+k);i++) {
            if( (pw[i*nz+iz0+k]<0.5) && ((iz0+k)>=nsurf[i]) ) {
               fillHole_2(i,iz0+k,nx,nz,5,gw,pw,nsurf);
            }
         }
         for(i=ix0-1;i>=MAX(0,ix0-k);i--) {
            if( (pw[i*nz+iz0+k]<0.5) && ((iz0+k)>=nsurf[i]) ) {
               fillHole_2(i,iz0+k,nx,nz,5,gw,pw,nsurf);
            }
         }
      }
      if( ix0-k >= 0 ) {
         for(i=iz0-1;i>=MAX(0,iz0-k);i--) {
            if( (pw[(ix0-k)*nz+i]<0.5) && (i>=nsurf[ix0-k]) ) {
               fillHole_2(ix0-k,i,nx,nz,5,gw,pw,nsurf);
            }
         }
         for(i=iz0;i<=MIN(nz-1,iz0+k);i++) {
            if( (pw[(ix0-k)*nz+i]<0.5) && (i>=nsurf[ix0-k]) ) {
               fillHole_2(ix0-k,i,nx,nz,5,gw,pw,nsurf);
            }
         }
      }
   }
   return;
}

void expand_3(int *iz00, float *gw, float *pw, int *nsurf, int nx, int nz)
{
   int i,j,k,ix0,iz0;
   int mppp=0;
   for( i=0;i<nx;i++ ) {
      if( iz00[i] < nz-1 ) {
         mppp=1;
         break;
      }
   }
   int ipp0;
   while( mppp > 0 ) {
      ipp0=nx/2;
      for( i=ipp0;i<nx;i++ ) {
      if( iz00[i] < nz-1 ) {
         ix0=i;
         iz0=iz00[i];
         if( (ix0>0) && (pw[(ix0-1)*nz+iz0]<0.5) && (iz0>=nsurf[ix0-1]) ) {
            fillHole_2(ix0-1,iz0,nx,nz,5,gw,pw,nsurf);
         }
         if( (iz0<nz-1) && (pw[ix0*nz+iz0+1]<0.5) && (iz0+1>=nsurf[ix0]) ) {
            fillHole_2(ix0,iz0+1,nx,nz,5,gw,pw,nsurf);
         }
         if( (ix0<nx-1) && (pw[(ix0+1)*nz+iz0]<0.5) && (iz0>=nsurf[ix0+1]) ) {
            fillHole_2(ix0+1,iz0,nx,nz,5,gw,pw,nsurf);
         }
      }
      }
      for( i=ipp0-1;i>=0;i-- ) {
      if( iz00[i] < nz-1 ) {
         ix0=i;
         iz0=iz00[i];
         if( (ix0>0) && (pw[(ix0-1)*nz+iz0]<0.5) && (iz0>=nsurf[ix0-1]) ) {
            fillHole_2(ix0-1,iz0,nx,nz,5,gw,pw,nsurf);
         }
         if( (iz0<nz-1) && (pw[ix0*nz+iz0+1]<0.5) && (iz0+1>=nsurf[ix0]) ) {
            fillHole_2(ix0,iz0+1,nx,nz,5,gw,pw,nsurf);
         }
         if( (ix0<nx-1) && (pw[(ix0+1)*nz+iz0]<0.5) && (iz0>=nsurf[ix0+1]) ) {
            fillHole_2(ix0+1,iz0,nx,nz,5,gw,pw,nsurf);
         }
      }
      }
      for( k=0;k<nx;k++ ) {
         int knz=iz00[k];
         for( j=knz+1;j<nz;j++ ) {
            if( pw[k*nz+j] > 0.5 ) {
               iz00[k]=j;
            }
         }
      }
      mppp=0;
      for( i=0;i<nx;i++ ) {
         if( iz00[i] < nz-1 ) {
            mppp=1;
            break;
         }
      }
   }
   return;
}

void fillHole_1(int ix0, int iz0, int nx, int nz, int k,
                double *g, int *nsurf)
{
    double tmp1,tmp2,tmp3;
    int i,j,n,m;
    int num=0;
    tmp1=0;tmp2=0;tmp3=0;
    for(i=-k;i<=k;i++) {
       if( ((i+ix0)>=0) && ((i+ix0) < nx) ) {
          for( j=-k;j<=k;j++) {
             if( ((j+iz0)>=0) && ((j+iz0) < nz) ) {
                n=i*i+j*j;
                m=(i+ix0)*nz+(j+iz0);
                if( (ABS(g[m])>0.0) && ((j+iz0)>=nsurf[i+ix0]) ) {
                   num++;
                   tmp1=1.0/sqrt(sqrt((double)(n)));
                   tmp2=tmp2+tmp1;
                   tmp3=tmp3+g[m]*tmp1;
                }
             }
          }
       }
    }
    g[ix0*nz+iz0]=tmp3/tmp2;
    return;
}

void fillHole_2(int ix0, int iz0, int nx, int nz, int k,
                float *g, float *pw, int *nsurf)
{
    double tmp1,tmp2,tmp3;
    int i,j,n,m;
    int num=0;
    tmp1=0;tmp2=0;tmp3=0;
    for(i=-k;i<=k;i++) {
       if( ((i+ix0)>=0) && ((i+ix0) < nx) ) {
          for( j=-k;j<=k;j++) {
             if( ((j+iz0)>=0) && ((j+iz0) < nz) ) {
                n=i*i+j*j;
                m=(i+ix0)*nz+(j+iz0);
                if( (pw[m]>0.5) && ((j+iz0)>=nsurf[i+ix0]) ) {
                   num++;
                   tmp1=1.0/sqrt(sqrt((double)(n)));
                   tmp2=tmp2+tmp1;
                   tmp3=tmp3+g[m]*tmp1;
                }
             }
          }
       }
    }
    g[ix0*nz+iz0]=tmp3/tmp2;
    pw[ix0*nz+iz0]=1.0;
    return;
}

void smooth2d(double *g, double *dump, int nx, int nz, 
              int irsx, int irsz, int *nsurf,
              int ix1_node, int ix2_node)
{
  double m1,aaa;
  int i,j,ind,l,m,im1,im2,l1,l2,il,jm,ind2;
  for( i=0;i<nx*nz;i++ ) {
     dump[i]=0.0;
  }
  for( i=ix1_node;i<=ix2_node;i++ ) {
     for( j=nsurf[i];j<nz;j++ ) {
        m1=0.0;
        ind=i*nz+j;
        l1=MAX(-i,-irsx);
        l2=MIN(nx-1-i,irsx);
        for( l=l1;l<=l2;l++ ) {
           il=i+l;
           im1=MAX(nsurf[il]-j,-irsz);
           im2=MIN(nz-1-j,irsz);
           for( m=im1;m<=im2;m++ ) {
              jm=j+m;
              ind2=il*nz+jm;
              if( ((m != 0) || (l != 0)) && NE(g[ind2],0.0) ) {
                 aaa=1.0/sqrt(sqrt(l*l*1.0+m*m*1.0));
                 dump[ind]=dump[ind]+g[ind2]*aaa;
                 m1=m1+aaa;
              }
           }
        }
        dump[ind]=(g[ind]+dump[ind])/(m1+1.0);
     }
  } 
  return;
}

double getamax1(double *slow, int nx, int nz, int *nsurf)
{
  double res=0.0;
  int i,j;
  for( i=0;i<nx;i++ ) {
     for( j=nsurf[i];j<nz;j++ ) {
        res=MAX(res,ABS(slow[i*nz+j]));
     }
  }
  return res;
}

double getmax(float *slow, int nx, int nz, int *nsurf)
{
  double res=0.0;
  int i,j;
  for( i=0;i<nx;i++ ) {
     for( j=nsurf[i];j<nz;j++ ) {
        res=MAX(res,slow[i*nz+j]);
     }
  }
  return res;
}

void sort_1(float *surfx, float *surfz, int front_size)
{
  if( front_size<= 1 ) return;
  double aln2i=1./0.69314718;
  double tiny=0.000001;
  int lognb2=(int)floor(log(front_size)*aln2i+tiny);
  int m=front_size;
  int nn, j;

  for(nn=0;nn<lognb2;nn++)
  {
    m=m/2;
    int k=front_size-m;
    for(j=0;j<k;j++)
    {
      int i=j;
      while( (surfx[i+m]<surfx[i])
           && (i>=0) ) {
        double tmp=surfx[i];
        double tmp2=surfz[i];
        surfx[i]=surfx[i+m];
        surfz[i]=surfz[i+m];
        surfx[i+m]=tmp;
        surfz[i+m]=tmp2;
        i=i-m;
      }
    }
  }
  return;
}

void sort_2(int *ind, int kn)
{
   if( kn <= 1 ) return;
   double aln2i=1./0.69314718;
   double tiny=0.000001;
   int lognb2=(int)floor(log(kn)*aln2i+tiny);
   int m=kn;
   register int i,j,nn;
   for( nn=0;nn<lognb2;nn++ ) {
      m=m/2;
      int k=kn-m;
      for( j=0;j<k;j++ ) {
         i=j;
         while( ((ind[i+m]<ind[i]))
           && (i>=0) ) {
           int n1=ind[i];
           ind[i]=ind[i+m];
           ind[i+m]=n1;
           i=i-m;
         }
      }
   }
   return;
}

float getvalue1d(float x,
                 int ndata, float *xcor, float *fx)
{
   float res;
   int i;

   if( LE(x,xcor[0]) ) {
     return fx[0];
   } else if( GE(x,xcor[ndata-1]) ) {
     return fx[ndata-1];
   } else {
     for(i=1;i<ndata;i++)
     {
        if( LE(x,xcor[i]) ) {
          double a=(xcor[i]-x)/(xcor[i]-xcor[i-1]);
          double b=1.0-a;
          return fx[i]*b+fx[i-1]*a;
        }
     }
   }
}
     
double getvalue2d(double x,double z, double h,
                  int nz, double *data)
{
    double ta=x/h;
    double ua=z/h;
    int ia=(int)floor(ta);
    int ja=(int)floor(ua);
    ta=ta-ia;
    ua=ua-ja;
    double tss;
    tss=(1.0-ta)*(1.0-ua)*data[ia*nz+ja]+
        ta*(1.0-ua)*data[(ia+1)*nz+ja]+
        ta*ua*data[(ia+1)*nz+ja+1]+
        (1.0-ta)*ua*data[ia*nz+ja+1];
    return tss;
}

void ray2d_tomo(double xr, double zr,
		double xs, double zs,
		double h,  double *t1,
	        int nx,	   int nz,
                int *ind_ray, int *is_visited,
		double trt_diff, 
	        double trt,
	        double ra1,
		int *nsurf, double *slow,
		double *ra, int *pre)
{
    int i,j;

    double zerodis=0.0001*h;
    double zerodis2=0.1;
    int nxr,nzr,nxs,nzs,i2i,nraymax;
    nxr=NINT(xr/h);nzr=NINT(zr/h);
    nxs=NINT(xs/h);nzs=NINT(zs/h);
    i2i=NINT(ra1*sqrt(ABS(nxr-nxs)*1.0));
    nraymax=50*(nx+nz);  
    
    double td1,td2,td3;
    td1=(xs<xr)?1.5:0.9;
    td2=(xs<xr)?0.9:1.5;
/*
    td1=1.5;td2=0.9;
*/
    td3=2.0;
    
    double x0=xr,z0=zr;
    double a=trt_diff/sqrt((x0-xs)*(x0-xs)
                   +(z0-zs)*(z0-zs));
    double t0=trt;

    int ii1=nxr,ii2=nzr;
    int ind=ii1*nz+ii2;
    int kn=0;
    ra[ind]=ra[ind]+a;
    pre[ind]=pre[ind]+1;

    double x01=0.0,z01=0.0;
    double x1,z1,x2,z2;
    double ta;
    int kkk=0,ia,ja;
    int nx1,nz1,nx2,nz2;

    do {
       x1=0.0,z1=0.0,x2=0.0,z2=0.0;
       if( kkk != 0 ) {
           if( (ABS(x0-x01)<zerodis) &&
               (ABS(z0-z01)<zerodis) ) {
               ta=100000.0;
               for( i=nx1;i<=nx2;i++ ) {
                    for( j=nz1;j<=nz2;j++ ) {
                         if( GE(ABS(i*h-x0),zerodis) ||
                             GE(ABS(j*h-z0),zerodis) ) {
                             ind=i*nz+j;
                             if( LE(t1[ind],ta) ) {
                                 ia=i;ja=j;
                                 ta=t1[ind];
                             }
                         }
                    }
               }
	       x0=ia*h;z0=ja*h;
               t0=ta;
           }
       } 
       x01=x0;z01=z0;
       kkk++;

       nx1=NINT((x0-td2*h)/h);
       nx2=NINT((x0+td1*h)/h);
       nz1=NINT((z0-td2*h)/h);
       nz2=NINT((z0+td1*h)/h);
       nx1=MAX(0,nx1);
       nx2=MIN(nx-1,nx2);
       nz1=MAX(0,nz1);
       nz2=MIN(nz-1,nz2);

       int ind1=nx1*nz+nz1;
       int ind2=nx2*nz+nz1;
       int ind3=nx1*nz+nz2;
       int ind4=nx2*nz+nz2;
       int flag_goto100=0;

       if( LE((t0-t1[ind1])*(t0-t1[ind2]), 0.0) ) {
           z1=nz1*h;
           double x11=nx1*h,x12=nx2*h,x13=(nx1+1)*h;
           double t11=t1[ind1];
           double t12=t1[ind2];
           double t13=t1[ind1+nz];
           if( NE((t11-t12)*(t11-t13)*(t12-t13),0.0) ) {
               x1=(t0-t12)*(t0-t13)*x11/((t11-t12)*(t11-t13))
                 +(t0-t11)*(t0-t12)*x13/((t13-t11)*(t13-t12))
                 +(t0-t11)*(t0-t13)*x12/((t12-t11)*(t12-t13));
           } else {
               x1=nx1*h+(nx2-nx1)*h/(t1[ind2]-t1[ind1])
                 *(t0-t1[ind1]);
           }
       }

       if( LE((t0-t1[ind1])*(t0-t1[ind3]),0.0) ) {
           double z11=nz1*h;
           double z12=(nz1+1)*h;
           double z13=(nz2+1)*h;
           double t11=t1[ind1];
           double t12=t1[ind1+1];
           double t13=t1[ind3];
	   if( EQ(z1,0.0) ) {
               x1=nx1*h;
	       if( NE((t11-t12)*(t11-t13)*(t12-t13),0.0) ) {
                   z1=(t0-t12)*(t0-t13)*z11/((t11-t12)*(t11-t13))
                     +(t0-t11)*(t0-t12)*z13/((t13-t11)*(t13-t12))
                     +(t0-t11)*(t0-t13)*z12/((t12-t11)*(t12-t13));
               } else {
                   z1=nz1*h+(nz2-nz1)*h*(t0-t1[ind1])/(t1[ind3]-t1[ind1]);
               }
           } else {
               x2=nx1*h;
               if( NE((t11-t12)*(t11-t13)*(t12-t13),0.0) ) {
                   z2=(t0-t12)*(t0-t13)*z11/((t11-t12)*(t11-t13))
                     +(t0-t11)*(t0-t12)*z13/((t13-t11)*(t13-t12))
                     +(t0-t11)*(t0-t13)*z12/((t12-t11)*(t12-t13));
               } else {
                   z2=nz1*h+(nz2-nz1)*h*(t0-t1[ind1])/(t1[ind3]-t1[ind1]);
               }
	       flag_goto100=1;
           }
       }

       if( LE((t0-t1[ind3])*(t0-t1[ind4]),0.0) && 
                   (flag_goto100 == 0) ) {
           double x11=nx1*h;
           double x12=nx2*h;
           double x13=(nx1+1)*h;
           double t11=t1[ind3];
           double t12=t1[ind4];
           double t13=t1[ind3+nz];
           if( EQ(z1,0.0) ) {
               z1=nz2*h;
               if( NE((t11-t12)*(t11-t13)*(t12-t13),0.0) ) {
                   x1=(t0-t12)*(t0-t13)*x11/((t11-t12)*(t11-t13))
                     +(t0-t11)*(t0-t12)*x13/((t13-t11)*(t13-t12))
                     +(t0-t11)*(t0-t13)*x12/((t12-t11)*(t12-t13));
               } else {
                   x1=nx1*h+(nx2-nx1)*h*(t0-t1[ind3])/(t1[ind4]-t1[ind3]);
               }
           } else {
               z2=nz2*h;
               if( NE((t11-t12)*(t11-t13)*(t12-t13),0.0) ) {
                   x2=(t0-t12)*(t0-t13)*x11/((t11-t12)*(t11-t13))
                     +(t0-t11)*(t0-t12)*x13/((t13-t11)*(t13-t12))
                     +(t0-t11)*(t0-t13)*x12/((t12-t11)*(t12-t13));
               } else {
                   x2=nx1*h+(nx2-nx1)*h*(t0-t1[ind3])/(t1[ind4]-t1[ind3]);
               }
               flag_goto100=1;
           }
       }

       if( LE((t0-t1[ind2])*(t0-t1[ind4]),0.0) &&
           (flag_goto100==0) ) {
           x2=nx2*h;
           double z11=nz1*h;
           double z12=(nz1+1)*h;
           double z13=(nz2+1)*h;
           double t11=t1[ind2];
           double t12=t1[ind2+1];
           double t13=t1[ind4];
           if( NE((t11-t12)*(t11-t13)*(t12-t13),0.0) ) {
               z2=(t0-t12)*(t0-t13)*z11/((t11-t12)*(t11-t13))
                 +(t0-t11)*(t0-t12)*z13/((t13-t11)*(t13-t12))
                 +(t0-t11)*(t0-t13)*z12/((t12-t11)*(t12-t13));
           } else {
               z2=nz1*h+(nz2-nz1)*h*(t0-t1[ind2])/(t1[ind4]-t1[ind2]);
           }
       }

       double x,z,t,zp;
       if( EQ((x1-x0)*(x2-x1)*(x2-x0),0.0) ) {
           if( NE(x2,x1) ) {
               zp=(z2-z1)/(x2-x1);
           } else if( NE(x1,x0) ) {
               zp=(z1-z0)/(x1-x0);
           } else if( NE(x2,x0) ) {
               zp=(z2-z0)/(x2-x0);
           }
       } else {
           double tmp1=x2-x0,tmp2=x1-x0;
	   tmp1=tmp1*tmp1;tmp2=tmp2*tmp2;
           zp=((z1*tmp1-z2*tmp2)-z0*(tmp1-tmp2))
              /((x1-x0)*(x2-x0)*(x2-x1));
       }
	   
       if( EQ(zp,0.0) ) {
           x=x0;
           t=t1[ind1]+(t1[ind2]-t1[ind1])*(x-nx1*h)/((nx2-nx1)*h);
           if( t < t0 ) {
               t0=t;
               x0=x;
               z0=nz1*h;
           } else {
               t=t1[ind3]+(t1[ind4]-t1[ind3])*(x-nx1*h)/((nx2-nx1)*h);
               if( t < t0 ) {
                   t0=t;
                   x0=x;
                   z0=nz2*h;
               }
           }
       } else {
           zp=1.0/zp;
           z=-zp*nx1*h+zp*x0+z0;
           t=t1[ind1]+(t1[ind3]-t1[ind1])*(z-nz1*h)/((nz2-nz1)*h);
	   if( LE((z-nz1*h)*(z-nz2*h),0.0) &&
               (t < t0) ) {
               t0=t;
               x0=nx1*h;
               z0=z;
           } else {
               z=-zp*nx2*h+zp*x0+z0;
               t=t1[ind2]+(t1[ind4]-t1[ind2])*(z-nz1*h)/((nz2-nz1)*h);
               if( LE((z-nz1*h)*(z-nz2*h),0.0) &&
                      (t < t0) ) {
                   t0=t;
                   x0=nx2*h;
                   z0=z;
               } else {
                   x=(zp*x0+z0-nz1*h)/zp;
                   t=t1[ind1]+(t1[ind2]-t1[ind1])*(x-nx1*h)
                           /((nx2-nx1)*h);
                   if( LE((x-nx1*h)*(x-nx2*h),0.0) &&
                       (t < t0) ) {
                       t0=t;
                       x0=x;
                       z0=nz1*h;
                   } else {
                       x=(zp*x0+z0-nz2*h)/zp;
                       t=t1[ind3]+(t1[ind4]-t1[ind3])*(x-nx1*h)
                               /((nx2-nx1)*h); 
	               if( LE((x-nx1*h)*(x-nx2*h),0.0) &&
                                   (t < t0) ) {
                           t0=t;
                           x0=x;
		           z0=nz2*h;
                       }
                   }
               }
           }
       }
       ii1=NINT(x0/h);
       ii2=NINT(z0/h);
       int ii3=MIN(nx-1,ii1+1);
       int jj1=MAX(0,ii2-i2i),jj2=MIN(nz-1,ii2+i2i+1);
       for( i=ii1;i<=ii3;i++ ) {
          for( j=MAX(jj1,nsurf[i]);j<=jj2;j++ ) {
             int ind_tmp=i*nz+j;
             ra[ind_tmp]=ra[ind_tmp]+a;
             pre[ind_tmp]=pre[ind_tmp]+1;
          }
       }
    } while( GE(sqrt((x0-xs)*(x0-xs)+(z0-zs)*(z0-zs)),td3*h) &&
                       (kkk <= nraymax) );

    ii1=nxs;
    ii2=nzs;
    int ind_tmp=ii1*nz+ii2;
    ra[ind_tmp]=ra[ind_tmp]+a;
    pre[ind_tmp]=pre[ind_tmp]+1;

    return;
}

void pad_1(double *slow,int nx,int nz)
{
   int i,j;
   for( i=5;i<nx-5;i++ ) {
      for( j=0;j<5;j++ ) {
         slow[i*nz+j]=slow[i*nz+5]*1.5;
      }
      for( j=nz-5;j<nz;j++ ) {
         slow[i*nz+j]=slow[i*nz+nz-6]*1.5;
      }
   }
   for( j=0;j<nz;j++ ) {
      for( i=0;i<5;i++ ) {
         slow[i*nz+j]=slow[5*nz+j]*1.5;
      }
      for( i=nx-5;i<nx;i++ ) {
         slow[i*nz+j]=slow[(nx-6)*nz+j]*1.5;
      }
   }
   return;
}
   
void get_assigned(int *id1, int *id2, int njob1, int njob2, int numprocs, int myrank)
{
   if( (numprocs < 1) || (myrank < 0) ) {
      printf("Something wrong! nump=%d, myrank=%d\n",numprocs,myrank);
      exit(1);
   }
   if( numprocs == 1 ) {
      *id1=njob1;
      *id2=njob2;
   } else {
      int nnn=njob2-njob1+1;
      int ii1=nnn/numprocs;
      int ii2=nnn-ii1*numprocs;
      if( myrank < (numprocs-ii2) ) {
         *id1=njob1+myrank*ii1;
         *id2=njob1+(myrank+1)*ii1-1;
        } else {
           *id1=njob1+ii2+myrank*ii1+myrank-numprocs;
           *id2=njob1+ii2+(myrank+1)*ii1+myrank-numprocs;
      }
   }
   return;
}
