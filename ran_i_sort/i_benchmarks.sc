

/*                                                                    
 *
 * i_benchmarks.sc - Sorting Benchmarks for Integers (Split-C Code)
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
 * Filename:            i_benchmarks.sc
 * History:
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <split-c/split-c.h>
#include <split-c/control.h>


/* The largest allowable integer value: */

#define MAX_VAL 2147483647 

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


void U_i_Benchmark(int col_size, int *In)

{

int i;

for (i=0;i<col_size;i++)
  In[i] = ((random()) & MAX_VAL); 
  
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


void G_i_Benchmark(int col_size, int *In)

{

int i;


for (i=0;i<col_size;i++)
  In[i] =  ((int) (random() & MAX_VAL)/4) +
	   ((int) (random() & MAX_VAL)/4) + 
	   ((int) (random() & MAX_VAL)/4) + 
	   ((int) (random() & MAX_VAL)/4);  

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


void gG_i_Benchmark(int col_size, int *In, int group)

{

int i,j,k,v,t1,t2;


v = ((MYPROC - (MYPROC % group)) + (PROCS/2)) % PROCS;
t1 = ((MAX_VAL >> PROCSLOG) + 1);

k=0;
for (i=0;i<group;i++)
  {t2 = v * t1;
   for (j=0;j<(col_size/group);j++)
     In[k++] = t2 + ((random() & MAX_VAL) >> PROCSLOG);
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


void B_i_Benchmark(int col_size, int *In)

{

int i,j,k,t1,t2;


t1 = ((MAX_VAL >> PROCSLOG) + 1);

k=0;
for (i=0;i<PROCS;i++)
  {t2 = i * t1;
   for (j=0;j<(col_size/PROCS);j++)
     In[k++] = t2 + ((random() & MAX_VAL) >> PROCSLOG);
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
       

void S_i_Benchmark(int col_size, int *In)

{

int i,k,t1,target;

if (MYPROC < (PROCS/2))
  target = 2*MYPROC + 1;
else 
  target = (MYPROC - (PROCS/2))*2;

t1 = target*((MAX_VAL >> PROCSLOG) + 1);

k=0;
for (i=0;i<col_size;i++)
  In[k++] = t1 + ((random() & MAX_VAL) >> PROCSLOG);
  
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


void Z_i_Benchmark(int col_size, int  *In)

{
int i;

for (i=0;i<col_size;i++)
  In[i] = 0;
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


void WR_i_Benchmark(int col_size, int *In, int s)

{

int i,j,k,t1,t2;

t1 = (MAX_VAL/PROCS) + 1;

k=0;

if (MYPROC == 0)
  {for (i=0;i<(PROCS-2);i+=2)
     {t2 = i * t1;
      for (j=1;j<(col_size/PROCS);j++)
	In[k++] = t2 + 1 + (random() % (t1 - 1));
      In[k++] = (t2 += t1);
      for (j=1;j<(col_size/s);j++)
	In[k++] = t2;
      for (j=(col_size/s);j<(col_size/PROCS);j++)
	In[k++] = t2 + 1 + (random() % (t1 - 1));
      In[k++] = t2 + t1;
     } 
   t2 = (PROCS - 2) * t1;
   for (j=1;j<(col_size/PROCS);j++)
     In[k++] = t2 + 1 + (random() % (t1 - 1));
   In[k++] = (t2 += t1);
   for (j=1;j<(col_size/s);j++)
     In[k++] = t2;
   for (j=(col_size/s);j<=(col_size/PROCS);j++)
     In[k++] = t2 + 1 + (random() % (t1 - 1));
  }
else
  {for (i=0;i<PROCS;i+=2)
     {t2 = i * t1;
      for (j=1;j<((col_size/PROCS) + (col_size/s));j++)
	In[k++] = t2 + 1 + (random() % (t1 - 1));
      In[k++] = t2 + t1 + MYPROC;
      for (j=0;j<((col_size/PROCS) - (col_size/s));j++)
	In[k++] = t2 + t1 + 1 + MYPROC + (random() % (t1 - 1 - MYPROC));
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
       

void DD_i_Benchmark(int col_size, int *In)

{

int i, tot, mask;

int *s_ptr, *f_ptr;

mask = (int) ((log((double) (col_size * PROCS)))/(log(2.00)));

if (MYPROC < PROCS_m1)
  {tot = 0;

   for (i=0;i<PROCSLOG;i++)
     if (MYPROC >= (tot += (PROCS/(2 << i)))) 
       mask--;

   for (i=0;i<col_size;i++)
     *(In + i) = mask;
  }
else
  {tot = col_size/2;
   mask = mask - PROCSLOG;
   f_ptr = (s_ptr = In) + col_size;
   while (tot >= 1)
     {for (i=0;i<tot;i++) 
	*(s_ptr++) = mask;
      tot = tot/2;
      mask--;
     }
   while (s_ptr < f_ptr)
     *(s_ptr++) = mask;
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
       

void RD_i_Benchmark(int col_size, int *In, int range)

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
     *(In + (j++)) = target;
  }

while (j < col_size)
  *(In + (j++)) = target;


free(Temp);

}











