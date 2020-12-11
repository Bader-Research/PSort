/*                                                                    
 *
 * d_benchmarks.sc - Sorting Benchmarks for Double Precision 
 *               Floating Point Numbers (Split-C Code)
 *
 * 
 * "Copyright (c) 1996 The Regents of the University of Maryland.
 * All rights reserved.
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose, without fee, and without written agreement is
 * hereby granted, provided that the above copyright notice and the following
 * two paragraphs appear in all copies of this software.
 * 
 * IN NO EVENT SHALL THE UNIVERSITY OF MARYLAND BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
 * OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
 * MARYLAND HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * THE UNIVERSITY OF MARYLAND SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
 * ON AN "AS IS" BASIS, AND THE UNIVERSITY OF MARYLAND HAS NO OBLIGATION TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS."
 *
 * Authors:             David R. Helman  <helman@umiacs.umd.edu>
			David A. Bader   <dbader@umiacs.umd.edu>
 *                      Joseph F. Ja'Ja' <joseph@umiacs.umd.edu>
 *                      Institute for Advanced Computer Studies
 *                      Department of Electrical Engineering 
 *                      AV Williams Building
 *                      College Park, MD 20742
 *                      
 * Version:             1.0
 * Creation Date:       August 15, 1996
 * Filename:            d_benchmarks.sc
 * History:
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <split-c/split-c.h>
#include <split-c/control.h>


/* The largest allowable double precision floating point value: */

#define MAX_VAL 1.7976931348623E+308

/* The largest possible value returned by the random() function plus one: */

#define RANDOM_MAX 2147483648.0

#define PROCS_m1 (PROCS-1)



/*************************************************************************** 
*
*			  UNIFORM BENCHMARK [U]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*
***************************************************************************/


void U_d_Benchmark(int col_size, double *In)

{
  
int i;

double multiplier,offset;

offset = RANDOM_MAX/2.0;
multiplier = MAX_VAL/offset; 

for (i=0;i<col_size;i++) 
  In[i] = (((double) random()) - offset) * multiplier;

}



/*************************************************************************** 
*
*			  GAUSSIAN BENCHMARK [G]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*
***************************************************************************/


void G_d_Benchmark(int col_size, double *In)

{

int i;

double multiplier,offset;

offset = 2.0*RANDOM_MAX;
multiplier = MAX_VAL/offset;

for (i=0; i<col_size;i++)
  In[i] = (((double) random()) + ((double) random()) + ((double) random()) + 
	   ((double) random()) - offset) * multiplier;

}



/*************************************************************************** 
*
*			  g-GROUP BENCHMARK [g-G]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*       group    -> the size of the processor groups.
*       
***************************************************************************/


void gG_d_Benchmark(int col_size, double *In, int group)

{

int i,j,k,v;

double multiplier,offset,t1,t2;

offset = RANDOM_MAX/2.0;
multiplier = MAX_VAL/offset; 
t1 = RANDOM_MAX/((double) PROCS);
v = ((MYPROC - (MYPROC % group)) + (PROCS/2)) % PROCS;


k=0;
for (i=0;i<group;i++)
  {t2 = (((double) v) * t1) - offset;
   for (j=0;j<(col_size/group);j++)
     In[k++] = (t2 + (((double) random())/((double) PROCS))) * multiplier;
   v = (++v) % PROCS; 
  } 

}



/*************************************************************************** 
*
*			  BUCKET SORTED BENCHMARK [B]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*       
***************************************************************************/


void B_d_Benchmark(int col_size, double *In)

{

int i,j,k;

double multiplier,offset,t1,t2;

offset = RANDOM_MAX/2.0;
multiplier = MAX_VAL/offset; 
t1 = RANDOM_MAX/((double) PROCS);

k=0;

for (i=0;i<PROCS;i++)
  {t2 = (((double) i) * t1) - offset;
   for (j=0;j<(col_size/PROCS);j++)
     In[k++] = (t2 + (((double) random())/((double) PROCS))) * multiplier;
  } 

}



/*************************************************************************** 
*
*			  STAGGERED BENCHMARK [S]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*
***************************************************************************/


void S_d_Benchmark(int col_size, double *In)

{

int i,k,target;

double multiplier,offset,t1;

offset = RANDOM_MAX/2.0;
multiplier = MAX_VAL/offset; 

if (MYPROC < (PROCS/2))
  target = 2*MYPROC + 1;
else 
  target = (MYPROC - (PROCS/2))*2;

t1 = (((double) target) * (RANDOM_MAX/((double) PROCS))) - offset;

k=0;
for (i=0;i<col_size;i++)
  In[k++] = (t1 + (((double) random())/((double) PROCS))) * multiplier;

}



/*************************************************************************** 
*			  ZERO BENCHMARK [Z]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*       
***************************************************************************/


void Z_d_Benchmark(int col_size, double  *In)

{
int i;

for (i=0;i<col_size;i++)
  In[i] = 0.0;
}



/*************************************************************************** 
*
*			  WORST-LOAD REGULAR BENCHMARK [WR]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*       s        -> the number of samples collected from each processor.
*       
***************************************************************************/


void WR_d_Benchmark(int col_size, double *In, int s)

{

int i,j,k,t1;

double multiplier,offset,t2;

offset = RANDOM_MAX/2.0;
multiplier = MAX_VAL/offset; 
t1 = (int) (RANDOM_MAX/((double) PROCS));

k=0;

if (MYPROC == 0)
  {for (i=0;i<(PROCS-2);i+=2)
     {t2 = ((double) (i * t1)) - offset;
      for (j=1;j<(col_size/PROCS);j++)
	In[k++] = (t2 + 1.0 + ((double) (random() % (t1 - 1)))) * multiplier;
      In[k++] = (t2 += ((double) t1)) * multiplier;
      for (j=1;j<(col_size/s);j++)
	In[k++] = t2 * multiplier;
      for (j=(col_size/s);j<(col_size/PROCS);j++)
	In[k++] = (t2 + 1.0 + ((double) (random() % (t1 - 1)))) * multiplier;
      In[k++] = (t2 + ((double) t1)) * multiplier;
     } 
   t2 = ((double) ((PROCS - 2) * t1)) - offset;
   for (j=1;j<(col_size/PROCS);j++)
     In[k++] = (t2 + 1.0 + ((double) (random() % (t1 - 1)))) * multiplier;
   In[k++] = (t2 += ((double) t1)) * multiplier;
   for (j=1;j<(col_size/s);j++)
     In[k++] = t2 * multiplier;
   for (j=(col_size/s);j<=(col_size/PROCS);j++)
     In[k++] = (t2 + 1.0 + ((double) (random() % (t1 - 1)))) * multiplier;
  }
else
  {for (i=0;i<PROCS;i+=2)
     {t2 = ((double) (i * t1)) - offset;
      for (j=1;j<((col_size/PROCS) + (col_size/s));j++)
	In[k++] = (t2 + 1.0 + ((double) (random() % (t1 - 1)))) * multiplier;
      In[k++] = (t2 + ((double) (t1 + MYPROC))) * multiplier;
      for (j=0;j<((col_size/PROCS) - (col_size/s));j++)
	In[k++] = (t2 + ((double) (t1 + 1 + MYPROC)) + 
		   ((double) (random() % (t1 - 1 - MYPROC)))) * multiplier;
     } 
  }

}



/*************************************************************************** 
*
*			  DETERMINISTIC DUPLICATES BENCHMARK [DD]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*
***************************************************************************/


void DD_d_Benchmark(int col_size, double *In)

{

int i, tot, mask;

double *s_ptr, *f_ptr;

mask = (int) ((log((double) (col_size * PROCS)))/(log (2.00)));

if (MYPROC < PROCS_m1)
  {tot = 0;

   for (i=0;i<PROCSLOG;i++)
     if (MYPROC >= (tot += (PROCS/(2 << i)))) 
       mask--;

   for (i=0;i<col_size;i++)
     *(In + i) = (double) mask;
  }
else
  {tot = col_size/2;
   mask = mask - PROCSLOG;
   f_ptr = (s_ptr = In) + col_size;
   while (tot >= 1)
     {for (i=0;i<tot;i++) 
	*(s_ptr++) = (double) mask;
      tot = tot/2;
      mask--;
     }
   while (s_ptr < f_ptr)
     *(s_ptr++) = (double)mask;
  }

}



/*************************************************************************** 
*
*			  RANDOMIZED DUPLICATES BENCHMARK [RD]: 
*
*
* Called with:
*       col_size -> the number of elements per processor. 
*       *In      -> a local pointer to the array to be initialized.       
*       range    -> the range of random values.
*
***************************************************************************/


void RD_d_Benchmark(int col_size, double *In, int range)

{

int i,j,tot,target;
int *Temp = (int *) malloc(range*sizeof(int));

tot = 0;

for (i=0;i<range;i++)
  tot += ((*(Temp + i)) = 1 + ((random())%range)); 
  
for (i=0;i<range;i++)  
  *(Temp + i) = ((int) ((double) (*(Temp + i))*
			(((double) col_size)/((double) tot))));

j = 0;


for (i=0;i<range;i++)  
  {if ((j + (*(Temp + i))) <= col_size)
     tot = j + (*(Temp + i));
   else
     tot = col_size;
   target = ((random()) % range);
   while (j < tot)
     *(In + (j++)) = (double) target;
  }

while (j < col_size)
  *(In + (j++)) = (double) target;


free(Temp);

}
