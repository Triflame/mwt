#ifndef PAR_SJM
#define PAR_SJM

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef CWP_SJM
#include "cwp_sjm.h"
#endif

/* parameter table */
typedef struct {
        char *name;             /* external name of parameter   */
        char *asciival;         /* ascii value of parameter     */
} pointer_table;

/* global variables declared and used internally */
static pointer_table *argtbl;   /* parameter table              */
static int nargs;               /* number of args that parse    */
static int tabled = FALSE;      /* true when parameters tabled  */
#define LINE_LENGTH 256
static char FILENAME[LINE_LENGTH];

#ifdef __cplusplus
extern "C" {
#endif

int EGetParInt(char *name, int num, int *ptr, int IS_REQUIRED, int default_value);
int EGetParLongInt(char *name, int num, long int *ptr, int IS_REQUIRED, long int default_value);
int EGetParFloat(char *name, int num, float *ptr, int IS_REQUIRED, float default_value);
int EGetParDouble(char *name, int num, double *ptr, int IS_REQUIRED, double default_value);
int EGetParStr(char *name, char *ptr, int IS_REQUIRED, char *default_value);
int GetParInt (char *name, int num, int *ptr);
int GetParLongInt (char *name, int num, long int *ptr);
int GetParFloat (char *name, int num, float *ptr);
int GetParDouble (char *name, int num, double *ptr);
int GetParStr (char *name, char *ptr);
int getnpar (int n, char *name, char *type, void *ptr);
void GetParInit(char *fname);
static int getparindex (int n, char *name);
static void tabulate ();
char **getparse(char *str, int *num);
int *getints(char *str, int *num);
long int *getlongint(char *str, int *num);
float *getfloats(char *str, int *num);
double *getdoubles(char *str, int *num);

#ifdef __cplusplus
}
#endif

int EGetParInt(char *name, int num, int *ptr, int IS_REQUIRED, int default_value)
{
   if( !GetParInt(name,num,ptr) ) {
       if( IS_REQUIRED ) {
           printf("Please add item %s in the deck file.\n",name);
           exit(2);
       } else {
           if( num == 1 ) {
              *ptr=default_value;
           }
       }
       return 0;
   }
   printf("%s=%d\n",name,*ptr);
   return 1;
}

int EGetParLongInt(char *name, int num, long int *ptr, int IS_REQUIRED, long int default_value)
{
   if( !GetParLongInt(name,num,ptr) ) {
       if( IS_REQUIRED ) {
           printf("Please add item %s in the deck file.\n",name);
           exit(2);
       } else {
           if( num == 1 ) {
              *ptr=default_value;
           }
       }
       return 0;
   }
   printf("%s=%d\n",name,*ptr);
   return 1;
}

int EGetParFloat(char *name, int num, float *ptr, int IS_REQUIRED, float default_value)
{
   if( !GetParFloat(name,num,ptr) ) {
       if( IS_REQUIRED ) {
           printf("Please add item %s in the deck file.\n",name);
           exit(2);
       } else {
           if( num == 1 ) {
              *ptr=default_value;
           }
       }
       return 0;
   }
   printf("%s=%f\n",name,*ptr);
   return 1;
}

int EGetParDouble(char *name, int num, double *ptr, int IS_REQUIRED, double default_value)
{  
   if( !GetParDouble(name,num,ptr) ) {
       if( IS_REQUIRED ) {
           printf("Please add item %s in the deck file.\n",name);
           exit(2);
       } else {
           if( num == 1 ) {
              *ptr=default_value;
           }
       }
       return 0;
   }
   printf("%s=%f\n",name,*ptr);
   return 1;
}

int EGetParStr(char *name, char *ptr, int IS_REQUIRED, char *default_value)
{
   if( !GetParStr(name,ptr) ) {
       if( IS_REQUIRED ) {
           printf("Please add item %s in the deck file.\n",name);
           exit(2);
       } else {
           strcpy(ptr,default_value);
       }
       return 0;
   }
   printf("%s=%s\n",name,ptr);
   return 1;
}

int GetParInt (char *name, int num, int *ptr)
{
        int i;                  /* index of name in symbol table        */
	

        i = getparindex(0,name);/* Get parameter index */
        if (i < 0) return 0;    /* Not there */

        int n;
        int *tmp;
        tmp = getints(argtbl[i].asciival, &n);
        if( n < num ) {
          printf("there are only %d < %d numbers as required\n",n,num);
          exit(2);
        }
	if( num <= 1 ) {
          ptr[0]=tmp[0];
        } else {
          for( i=0;i<num;i++ ) {
             ptr[i]=tmp[i];
          }
	}
        free(tmp);
        return 1;
}

int GetParLongInt (char *name, int num, long int *ptr)
{
        int i;                  /* index of name in symbol table        */
	

        i = getparindex(0,name);/* Get parameter index */
        if (i < 0) return 0;    /* Not there */

        int n;
        long int *tmp;
        tmp = getlongint(argtbl[i].asciival, &n);
        if( n < num ) {
          printf("there are only %d < %d numbers as required\n",n,num);
          exit(2);
        }
	if( num <= 1 ) {
          *ptr=tmp[0];
        } else {
          *ptr=tmp[num-1];
	}
        free(tmp);
        return 1;
}

int GetParFloat (char *name, int num, float *ptr)
{
        int i;                  /* index of name in symbol table        */

        i = getparindex(0,name);/* Get parameter index */
        if (i < 0) return 0;    /* Not there */

        int n;
        float *tmp;
        tmp = getfloats(argtbl[i].asciival, &n);
        if( n < num ) {
          printf("there are only %d < %d numbers as required\n",n,num);
          exit(2);
        }
	if( num <= 1 ) {
          *ptr=tmp[0];
        } else {
          *ptr=tmp[num-1];
	}
        free(tmp);
        return 1;
}

int GetParDouble (char *name, int num, double *ptr)
{
        int i;                  /* index of name in symbol table        */

        i = getparindex(0,name);/* Get parameter index */
        if (i < 0) return 0;    /* Not there */

        int n;
	double *tmp;
        tmp = getdoubles(argtbl[i].asciival, &n);
        if( n < num ) {
          printf("there are only %d < %d numbers as required\n",n,num);
          exit(2);
        }
	if( num <= 1 ) {
          *ptr=tmp[0];
        } else {
          *ptr=tmp[num-1];
	}
        free(tmp);
        return 1;
}

int GetParStr (char *name, char *ptr)
{
        int i;                  /* index of name in symbol table        */

        i = getparindex(0,name);/* Get parameter index */
        if (i < 0) return 0;    /* Not there */

        strcpy(ptr,argtbl[i].asciival);
        return 1;

}

/*
 * Return the index of the n'th occurrence of a parameter name,
 * except if n==0, return the index of the last occurrence.
 * Return -1 if the specified occurrence does not exist.
 */
static int getparindex (int n, char *name)
{
        int i;
        if (n==0) {
                for (i=nargs-1; i>=0; --i) {
                        if (!strcmp(name,argtbl[i].name)) break;
                }
                return i;
        } else {
                for (i=0; i<nargs; ++i)
                        if (!strcmp(name,argtbl[i].name))
                                if (--n==0) break;
                if (i<nargs)
                        return i;
                else
                        return -1;
        }
}


void GetParInit(char *fname)
{
	if( tabled==TRUE ) {
          free(argtbl);
        }
        nargs=0;
        argtbl=(pointer_table*)malloc(100*sizeof(pointer_table));
	if( !argtbl ) {
          printf("Memory request failed\n");
          exit(1);
        }
        strcpy(FILENAME,fname);
        tabulate();
        tabled=TRUE;
        return;
}

/* Install symbol table */
static void tabulate ()
{

        int i;
        char *eqptr;
        char line[LINE_LENGTH];
        int debug=FALSE;
        
        FILE *fp;
        if( (fp=fopen(FILENAME,"r"))== NULL ) {
          printf("cannot open file %s\n",FILENAME);
          exit(0);
        }
        while( !feof(fp) ) {
          fgets(line,LINE_LENGTH,fp);
          if( eqptr=strchr( line, '#' )  ) {
            *eqptr='\0';
          }
          if( eqptr=strchr( line, '=' )  ) {
            *eqptr='\0';
            eqptr++;
            int n1=strlen(line),n2=strlen(eqptr);
            if( n1>0 ) {
              int i,i1,i2;
	      i1=n1;
              for( i=0;i<n1;i++ ) {
                if( isgraph(line[i]) ) {
                  i1=i;
                  break;
                }
              }
	      i2=-1;
              for( i=n1-1;i>=0;i-- ) {
                if( isgraph(line[i]) ) {
                  i2=i;
                  break;
                }
              }
	      char *tmp,*tmp2;
	      tmp=(char *)malloc((i2-i1+2)*sizeof(char));
              if( !tmp ) {
                printf("memory request failed.\n");
                exit(1);
              }
              for( i=i1;i<=i2;i++ ) {
                tmp[i-i1]=line[i];
              }
              tmp[i2+1-i1]='\0';
              argtbl[nargs].name=&tmp[0];
              if( n2 > 0 ) {
                i1=n2;
                for( i=0;i<n2;i++ ) {
                  eqptr[i] &= 0x7F;
                  if( isgraph(eqptr[i]) ) {
                    i1=i;
                    break;
                  }
                }
                i2=-1;
                for( i=n2-1;i>=0;i-- ) {
                  eqptr[i] &= 0x7F;
                  if( isgraph(eqptr[i]) ) {
                    i2=i;
                    break;
                  }
                }
                if( (i2-i1+2) > 0 ) {
	           tmp2=(char *)malloc((i2-i1+2)*sizeof(char));
                   if( !tmp2 ) {
                      printf("memory request failed.\n");
                      exit(1);
                   }
                   for( i=i1;i<=i2;i++ ) {
                      tmp2[i-i1]=eqptr[i];
                   }
                   tmp2[i2+1-i1]='\0';
                } else {
                   tmp2=(char *)malloc(sizeof(char));
                   tmp2[0]='\0';
                }
                argtbl[nargs].asciival=&tmp2[0];
              }
              nargs++;
            }
          }
        }
	fclose(fp);
        realloc(argtbl,nargs*sizeof(pointer_table));
        return;
}

int *getints(char *str, int *num)
{
   int *res;
   char **p1;
   int n,i,tmp;

   p1=getparse(str,&n);
   if( n==0 ) return NULL;
   res=(int *)malloc(sizeof(int));
   if( !res ) {
     printf("memory request failed\n");
     exit(1);
   }
   for( i=0;i<n;i++ ) {
     tmp=atoi(p1[i]);
     res[i]=tmp;
   }
   *num=n;
   free(p1);
   return res;
}

long int *getlongint(char *str, int *num)
{
   long int *res;
   char **p1;
   int n,i;

   p1=getparse(str,&n);
   if( n==0 ) return NULL;
   res=(long int *)malloc(sizeof(long int));
   if( !res ) {
     printf("memory request failed\n");
     exit(1);
   }
   for( i=0;i<n;i++ ) {
     res[i]=(long int)atof(p1[i]);
   }
   *num=n;
   free(p1);
   return res;
}

float *getfloats(char *str, int *num)
{
   float *res;
   char **p1;
   int n,i;

   p1=getparse(str,&n);
   if( n==0 ) return NULL;
   res=(float *)malloc(sizeof(float));
   if( !res ) {
     printf("memory request failed\n");
     exit(1);
   }
   for( i=0;i<n;i++ ) {
     res[i]=(float)atof(p1[i]);
   }
   *num=n;
   free(p1);
   return res;
}

double *getdoubles(char *str, int *num) 
{
   char **p1;
   int n,i;
   double tmp;
	
   p1=getparse(str,&n);
   if( n==0 ) return NULL;
   double *res;
   res=(double *)malloc(sizeof(double));
   if( !res ) {
     printf("memory request failed\n");
     exit(1);
   }
   for( i=0;i<n;i++ ) {
     tmp=atof(p1[i]);
     res[i]=tmp;
   }
   *num=n;
   free(p1);
   return &res[0];
}
 
char **getparse(char *str, int *num)
{
   char *str_tmp,**res;
   res=(char **)malloc(sizeof(char));
   if( !res ) {
     printf("memory request failed\n");
     exit(1);
   }
   int i=0,nstr=0,nlen;
   nlen=strlen(str);
   do {
     if( i==nlen ) break;
     while( ((str[i]==' ') || (str[i]==',')) && (i<nlen) ) {
       i++;
     }
     if( i==nlen ) break;
     str_tmp=(char *)malloc(sizeof(char));
     if( !str_tmp ) {
       printf("memory request failed\n");
       exit(1);
     }
     res[nstr]=str_tmp;
     int j=0;
     while( (str[i]!=' ') && (str[i]!=',') && (i<nlen) ) {
       str_tmp[j]=str[i];
       str_tmp[j] &= 0x7f;
       i++;j++;
     }
     str_tmp[j]='\0';
     nstr++;
     i++;
   } while( i<nlen );
   *num=nstr;
   return res;
}

#endif
