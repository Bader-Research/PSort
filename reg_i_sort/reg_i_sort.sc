
/*                                                                    
 *
 * reg_i_sort.sc - Function to Sort Integers by a New Deterministic
 *                 Algorithm Based on the Regular Sampling Approach
 *                 (Split-C Code)
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
 * Filename:            reg_i_sort.sc
 * History:
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <split-c/split-c.h>
#include <split-c/control.h>

/* The largest allowable integer value: */

#define MAX_VAL 2147483647 

#define MAX_VAL_m1 (MAX_VAL-1)

#define PROCS_m1 (PROCS-1)

typedef struct split
  {int val;
   int count;
  } split_t;


void assert_spread_malloc(void *spread spreadPtr)
/* Check spreadPtr to make sure it is not NULL */
{
    if (spreadPtr == NULL) {
	fprintf(stderr,"ERROR: a SpreadPtr is NULL\n");
	fflush(stderr);
	exit(1);
    }
}

void assert_malloc(void *ptr) {
    if (ptr==NULL) {
	fprintf(stderr,"ERROR: PE%2d cannot malloc\n",MYPROC);
	fflush(stderr);
	exit(1);
    }
}


/**************************************************************************
*
*                     ALL_REGULAR_I_SORT
*
*
*      The function all_regular_i_sort is called with the following
* parameters.
*
*   (1) col_size -> the number of elements per processor in the input array, 
* 
*   (2) s -> the number of samples to be collected from each processor.
*
*   (3) Buffer1 & 
*   (4) Buffer2 ->  spread pointers to the base of two spread arrays at 
*                   processor P_0.  When all_regular_i_sort is called, 
*                   Buffer1 holds the input array to be sorted.  As sorting 
*                   progresses, the data is moved back and forth between 
*                   the two arrays until at completion either array may 
*                   hold the sorted output.  Thus, either array must be large enough 
*                   to hold the maximum number of elements possible
*                   at a processor at the completion of sorting, which is
*                   (col_size + col_size*PROCS/s - PROCS) elements.  But if 
*                   this worst case occurs, then the array receiveing the 
*                   transpose in Step (7) must actually be of size 
*                   PROCS::(col_size/PROCS + col_size/s + 5*PROCS).
*                               
*   (5) *out_array_ptr -> a local pointer to a global pointer to the
*                         array which holds the completely sorted output.
*
*   (6) *val_ptr-> a local pointer to the variable which describes how
*                  many items are present at that processor from the          
*                   completely sorted output.
*
*   (7) bits -> the number of bits used to represent the sorting radix  
*             (i.e. - the sorting radix will be 2^(bits)).
*
*   (8) value_range -> the largest possible integer value in the set to
*                      be sorted, where `value_range' can not exceed  
*                      MAX_VAL.
*
***************************************************************************/


void all_regular_i_sort(int col_size,int s, int *spread Buffer1, 
			int *spread Buffer2, int *spread *out_array_ptr,
			int *val_ptr, int bits, int value_range)


{


register int i,j,k,t,entries,values,v,times,stride,rem,dest,tot,tot1,tot2,
             adjustments,temp1,temp2,max,done,l,r,begin,digits,radix,
             radix_m1,shift1,shift2,log_value_range,tval,stride1,start_row,
             l_start_col,r_start_col,finish_row,l_finish_col,r_finish_col,
             pre_t,cur_c,values_m1,buck_size,ones,PROCS_f,over;

int *info_ptr,*A_ptr,*B_ptr,*C_ptr,*temp_ptr,**bin_ptr,*start_ptr,
    *finish_ptr, *s1_ptr, *s2_ptr, *D_ptr,*Count,*Next,*Offset,*t_ptr,
    *c_ptr,*s_ptr,*s_sample_ptr,*f_sample_ptr,**Bin,**Data,*Bin_Count,
    *Info,*Number,*Count1,*Count2,**Loc,*Sample1,*Sample2;

split_t *split_ptr;

int *spread start_sptr;
int *spread finish_sptr;
int *spread temp_sptr;
int *spread A;
int *spread B;
int *spread C;
int *spread D;

split_t *spread Splitters;

int *global dest_ptr;

Bin = (int **) malloc(PROCS*PROCS*sizeof(int *));

assert_malloc(Bin);

Data = (int **) malloc(PROCS*sizeof(int *));

assert_malloc(Data);

Bin_Count = (int *) malloc(PROCS*PROCS*sizeof(int));

assert_malloc(Bin_Count);

Info = (int *) malloc(4*PROCS*PROCS*sizeof(int));

assert_malloc(Info);

Number = (int *) malloc(PROCS*sizeof(int));

assert_malloc(Number);

Count1 = (int *) malloc(((1 << bits)+1)*sizeof(int));

assert_malloc(Count1);

Count2 = (int *) malloc(((1 << bits)+1)*sizeof(int));

assert_malloc(Count2);

Loc = (int **) malloc(PROCS*sizeof(int *));

assert_malloc(Loc);

Sample1 = (int *) malloc(((PROCS*s) + 2*PROCS)*sizeof(int));

assert_malloc(Sample1);

Sample2 = (int *) malloc(((PROCS*s) + 2*PROCS)*sizeof(int));

assert_malloc(Sample2);

Splitters = all_spread_malloc(PROCS,PROCS*sizeof(split_t));

assert_spread_malloc(Splitters);


/***************************************************************************
			       Step (1)
***************************************************************************/


over = col_size/PROCS + col_size/s + 5*PROCS;
buck_size = col_size/PROCS;
PROCS_f = 4*PROCS;

radix = (1 << bits);  
radix_m1 = radix - 1;
values = col_size;

log_value_range = (int) ceil(log((double) value_range)/log((double) 2));

digits = (int) ceil(((double) log_value_range)/((double) bits));

rem = log_value_range % bits;

if (rem == 0)
  rem = bits;

start_sptr = Buffer1;
finish_sptr = Buffer2;

start_ptr = ((int *) (start_sptr + MYPROC));
finish_ptr = ((int *) (finish_sptr + MYPROC));
Count = (Offset = Count1) + 1;
Next = Count2 + 1;


if (digits > 1)
  {t_ptr = (c_ptr = Count) + radix;
   while (c_ptr < t_ptr)
     *(c_ptr++) = 0;
   t_ptr = (c_ptr = Next) + radix;
   while (c_ptr < t_ptr)
     *(c_ptr++) = 0;
   t_ptr = (c_ptr = start_ptr) + values;
   while (c_ptr < t_ptr)
     (*(Count + ((*(c_ptr++)) & radix_m1)))++;
   *(Offset) = 0;
   t_ptr = (c_ptr = Offset) + radix;
   while ((++c_ptr) < t_ptr)
     *c_ptr +=  *(c_ptr - 1);
   t_ptr = (c_ptr = start_ptr) + values;
   while (c_ptr < t_ptr)
     {(*(Next + (((v = *(c_ptr++)) >> bits) & radix_m1)))++;
      *(finish_ptr + ((*(Offset + (v & radix_m1)))++)) = v;
     }
   temp_sptr = start_sptr;
   start_sptr = finish_sptr;
   finish_sptr = temp_sptr;
   start_ptr = ((int *) (start_sptr + MYPROC));
   finish_ptr = ((int *) (finish_sptr + MYPROC));

   temp_ptr = Count;
   Count = Next;
   Next = temp_ptr;
   Offset = Count - 1;
   for (j = 1;j < (digits - 1);j++)
     {shift1 = j * bits;
      shift2 = (j + 1)*bits;
      t_ptr = (c_ptr = Next) + radix;
      while (c_ptr < t_ptr)
	*(c_ptr++) = 0;
      *(Offset) = 0;
      t_ptr = (c_ptr = Offset) + radix;
      while ((++c_ptr) < t_ptr)
	*c_ptr +=  *(c_ptr - 1);
      t_ptr = (c_ptr = start_ptr) + values;
      while (c_ptr < t_ptr)
	{(*(Next + (((v = (*(c_ptr ++))) >> shift2) & radix_m1)))++;
	 *(finish_ptr + ((*(Offset + ((v >> shift1) & radix_m1)))++)) = v;
	
	}
      temp_sptr = start_sptr;
      start_sptr = finish_sptr;
      finish_sptr = temp_sptr;
      start_ptr = ((int *) (start_sptr + MYPROC));
      finish_ptr = ((int *) (finish_sptr + MYPROC));
 
      temp_ptr = Count;
      Count = Next;
      Next = temp_ptr;
      Offset = Count - 1;
     }
   t = buck_size + 2; 
   temp_ptr = finish_ptr - t + 1;
   for(i=0;i<PROCS;i++)
     *(Loc + i) = (temp_ptr += t);
   radix_m1 = (radix = (1 << rem)) - 1;
   shift1 = (digits - 1) * bits;
   *(Offset) = 0;
   t_ptr = (c_ptr = Offset) + radix;
   while ((++c_ptr) < t_ptr)
     *c_ptr +=  *(c_ptr - 1);
   t_ptr = (c_ptr = start_ptr) + values;
   while (c_ptr < t_ptr)
     {t = ((*(Offset + (((*(c_ptr)) >> shift1) & radix_m1)))++);
      *((*(Loc + (t & PROCS_m1))) + (t >> PROCSLOG)) = *(c_ptr++);
     }
  } 
else
  {t = buck_size + 2; 
   temp_ptr = finish_ptr - t + 1;
   for(i=0;i<PROCS;i++)
     *(Loc + i) = (temp_ptr += t);
   radix_m1 = (radix = (1 << rem)) - 1;
   t_ptr = (c_ptr = Count) + radix;
   while (c_ptr < t_ptr)
     *(c_ptr++) = 0;
   t_ptr = (c_ptr = start_ptr) + values;
   while (c_ptr < t_ptr)
     (*(Count + ((*(c_ptr++)) & radix_m1)))++;
   *(Offset) = 0;
   t_ptr = (c_ptr = Offset) + radix;
   while ((++c_ptr) < t_ptr)
     *c_ptr +=  *(c_ptr - 1);
   t_ptr = (c_ptr = start_ptr) + values;
   while (c_ptr < t_ptr)
     {t = (*(Offset + (((*(c_ptr)) >> shift1) & radix_m1)))++;
      *((*(Loc + (t & PROCS_m1))) + (t >> PROCSLOG)) = *(c_ptr++);
     }
  }

barrier();


/***************************************************************************
			       Step (2)
***************************************************************************/


#if (defined(SP2))

A = finish_sptr;
B = start_sptr;
max = buck_size + 2;
A_ptr = ((int *) (A + MYPROC)); 
B_ptr = ((int *) (B + MYPROC));

barrier();

mpc_index(A_ptr,B_ptr,(max*sizeof(int)),ALLGRP);

barrier();

#else

A = finish_sptr;
B = start_sptr;
max = buck_size + 2;
A_ptr = (*(Loc + MYPROC)); 
t = max*MYPROC + 1;
B_ptr = (((int *) (B + MYPROC)) + t);

for (i=0;i<buck_size;i++)
  *(B_ptr + i) = *(A_ptr + i);

for (i=1;i<PROCS;i++)
  {dest = (MYPROC + i) & PROCS_m1;
   dest_ptr = ((int *global) (B + dest)) + t;
   bulk_store(dest_ptr,(*(Loc + dest)),(buck_size*(sizeof(int))));
   all_store_sync();
  }

barrier();

#endif


/***************************************************************************
			       Step (3)
***************************************************************************/


start_ptr = ((int *) (B + MYPROC));  
finish_ptr = Sample1;

adjustments = 0;
stride = col_size/(PROCS*s);
stride1 = buck_size - stride;

if (MYPROC == PROCS_m1)
  {for (j=0;j<(PROCS >> 1);j++)
     {s2_ptr = (s1_ptr = (start_ptr + stride)) + max;
      (*(Number + j)) = (s << 1);
      v = stride1;
      while ((*(s1_ptr + v)) == MAX_VAL)
	{adjustments++;
	 (*(s1_ptr + (v--))) = MAX_VAL_m1; 
	} 
      v = stride1;
      while ((*(s2_ptr + v)) == MAX_VAL)
	{adjustments++;
	 (*(s2_ptr + (v--))) = MAX_VAL_m1; 
	} 
      temp_ptr = finish_ptr + (s << 1);
      l = r = s;
      while (finish_ptr < temp_ptr)
	{if ((l) && (r))
	   {if ((*s1_ptr) < (*s2_ptr))
	      {*(finish_ptr++) = *(s1_ptr);
	       l--;
	       s1_ptr += stride;
	      }
	    else 
	      {*(finish_ptr++) = *(s2_ptr);
	       r--;
	       s2_ptr += stride;  
	      } 
	   }   
	 else
	   {if (r)
	      {*(finish_ptr++) = (*s2_ptr);
	       s2_ptr += stride;
	      }
	    else
	      {*(finish_ptr++) = (*s1_ptr);
	       s1_ptr += stride;
	      }
	   }
	}
      start_ptr += (2*max);
      *(finish_ptr++) = MAX_VAL;  
     }

   if (PROCSLOG > 1)
     {s_sample_ptr = Sample1;
      f_sample_ptr = Sample2;
     }
   else
     f_sample_ptr = Sample1;
   for (i=1;i<PROCSLOG;i++)
     {start_ptr = s_sample_ptr;
      finish_ptr = f_sample_ptr;
      for (j=0;j<(PROCS >> i);j+=2)
	{tot1 = (*(Number+j));
	 tot2 = (*(Number+j+1));
	 s2_ptr = (s1_ptr = start_ptr) + tot1 + 1;
	 (*(Number+(j>>1))) = tot = tot1 + tot2;
	 temp_ptr = finish_ptr + tot;
	 while (finish_ptr < temp_ptr)
	   {if ((*s1_ptr) < (*s2_ptr))
	      *(finish_ptr++) = *(s1_ptr++);
	    else 
	      *(finish_ptr++) = *(s2_ptr++);
	   } 
	 start_ptr += (tot + 2);
	 *(finish_ptr++) = MAX_VAL;  
	}
      if (i < (PROCSLOG - 1))
	{temp_ptr = f_sample_ptr;
	 f_sample_ptr = s_sample_ptr;
	 s_sample_ptr = temp_ptr;
	}
     }

   finish_ptr = f_sample_ptr; 
   values = v = *Number;
   for (i=0;i<adjustments;i++)  
     *(finish_ptr + (--v)) = MAX_VAL;
  }



/***************************************************************************
			       Step (4)
***************************************************************************/


stride = col_size/(s*PROCS);

if (MYPROC == PROCS_m1)
  {temp_ptr = f_sample_ptr - 1;
   split_ptr = ((split_t *) (Splitters + PROCS_m1)) - 1;
   
   for (i=1;i<PROCS;i++)
     {tval = (++split_ptr)->val =  *(temp_ptr += s);
      if (tval > (*(temp_ptr - 1)))
	split_ptr->count = stride;
      else
	{if (tval == (*(temp_ptr - s + 1)))
	   split_ptr->count = buck_size;
	 else
	   {l = ((i - 1) * s) + 1; 
	    r = (i * s) - 2;
	    t = (l + r) >> 1;
	    done = 0;
	    while (! done)
	      {if (tval > (*(f_sample_ptr + t)))
		 {l = t + 1;
		  t = (l + r) >> 1;
		 }
	       else
		 {if (tval > (*(f_sample_ptr + t - 1)))
		    {split_ptr->count = ((i * s) - t) * stride;
		     done = 1;
		    }
		  else
		    {r = t - 1;
		     t = (l + r) >> 1;
		    }
		 }
	      }
	   }
	}
     }
  }

barrier();


/***************************************************************************
			       Step (5)
***************************************************************************/


v = MYPROC; 

ones = 0;

while ((v & 1))
  {ones++;
   v >>= 1;
  }

#if (defined(T3D))

ones--;

for (i=PROCSLOG-1;i>=0;i--)
  {if (ones >= i)
      bulk_store(((split_t *global) (Splitters + MYPROC - (1 << i))),
	     ((split_t *) (Splitters + MYPROC)),(PROCS_m1 * (sizeof(split_t))));
   all_store_sync();
  }

barrier();

#else

if (MYPROC < PROCS_m1)
  store_sync(PROCS_m1*(sizeof(split_t)));

reset_store_sync();

for (i=(ones-1);i>=0;i--)
  bulk_store(((split_t *global) (Splitters + MYPROC - (1 << i))),
	     ((split_t *) (Splitters + MYPROC)),(PROCS_m1 * (sizeof(split_t))));

barrier();

#endif


/***************************************************************************
			       Step (6)
***************************************************************************/


B_ptr = ((int *) (B + MYPROC)) - max + 1;
split_ptr = ((split_t *) (Splitters + MYPROC));
(split_ptr + PROCS_m1)->val = MAX_VAL;
bin_ptr = Bin;
info_ptr = Info;
values_m1 = (values = buck_size) - 1;

for (i=0;i<PROCS;i++)
  {begin =  -1;
   done = l = 0;
   t = (r = values_m1) >> 1; 
   *(bin_ptr++) = (B_ptr += max);
   *(info_ptr++) = 0;    
   *(info_ptr++) = 0;    
   for (j=0;j<PROCS;j++)
     {if ((j) && (tval == ((split_ptr + j)->val)))
	{if (j == PROCS_m1)
	   {*(info_ptr++) = 0;    
	    *(info_ptr++) = values_m1 - begin;  
	   } 
	 else  
	   {*(info_ptr++) = 0;    
	    rem = pre_t - begin;
	    cur_c = (split_ptr + j)->count;
	    if ((*B_ptr) > tval)
	      {*(info_ptr++) = 0; 
	       *(info_ptr++) = 0; 
	       *(info_ptr++) = 0; 
	       *(bin_ptr++) = B_ptr;
	      }
	    else
	      {if (rem == 0)
		 {*(info_ptr++) = 0;  
		  *(info_ptr++) = 0;  
		  *(info_ptr++) = l + 1;  
		  *(bin_ptr++) = B_ptr + l + 1;
		 }
	       else
		 {if (cur_c >= rem)
		    {*(info_ptr++) = rem; 
		     (split_ptr + j)->count -= rem;
		     t = pre_t;
		     if (t == values_m1)
		       {l = r = begin = values_m1;
			*(bin_ptr++) = B_ptr + values_m1;
			*(info_ptr++) = 0;  
			*(info_ptr++) = values;  
		       } 
		     else  
		       {l = begin = t;
			r = values_m1;
			t = (l + r) >> 1;
			*(bin_ptr++) = B_ptr + l + 1;
			*(info_ptr++) = 0;    
			*(info_ptr++) = l + 1; 
		       }
		    }
		  else  
		    {*(info_ptr++) = cur_c; 
		     (split_ptr + j)->count = 0;
		     if (tval != ((split_ptr + j + 1)->val))
		       t = pre_t;
		     else
		       t = begin + cur_c;
		     if (t == values_m1)
		       {if (tval != ((split_ptr + j + 1)->val))
			  {*(info_ptr++) = rem - cur_c; 
			   *(bin_ptr++) = B_ptr + begin + cur_c + 1; 
			  }
			else
			  {*(info_ptr++) = 0; 
			   *(bin_ptr++) = B_ptr + values_m1; 
			  }
			l = r = begin = values_m1;
			*(info_ptr++) = values; 
		       } 
		     else  
		       {if (tval != ((split_ptr + j + 1)->val))
			  {*(info_ptr++) = rem - cur_c; 
			   *(bin_ptr++) = B_ptr + begin + cur_c + 1; 
			  }
			else
			  {*(info_ptr++) = 0;   
			   *(bin_ptr++) = B_ptr + t + 1; 
			  }
			l = begin = t;
			r = values_m1;
			t = (l + r) >> 1;
			*(info_ptr++) = pre_t + 1; 
		       }
		    }
		 }
	      }
	   }
	}
      else  
	{if (j == PROCS_m1)
	   {
	    *(info_ptr++) = values_m1 - begin; 
	    *(info_ptr++) = 0;
	   }
	 else  
	   {tval = (split_ptr + j)->val; 
	    while (! done)
	      {if (tval <= (*(B_ptr + t)))
		 {if (! t)
		    done = 1;
		  else
		    {r = t - 1;
		     t = (l + r) >> 1;
		    }
		 }
	       else
		 {if ((t == (values_m1)) || (tval <= (*(B_ptr + t + 1)))) 
		    done = 1;
		  else
		    {l = t + 1;
		     t = (l + r) >> 1;
		    }
		 }
		}
	    if ((t == values_m1) || 
		((t) && (tval < (*(B_ptr + t + 1)))) || 
		((!t) && ((tval < (*B_ptr)) || ((tval > (*B_ptr)) && (tval < (*(B_ptr + 1)))))))  
	      {pre_t = t;
	       if (tval < (*B_ptr))
		 {*(info_ptr++) = 0;
		  *(info_ptr++) = 0;
		  *(info_ptr++) = 0;
		  *(info_ptr++) = 0;
		  *(bin_ptr++) = B_ptr;
		  begin = -1;
		  l = 0;
		  t = (r = values_m1) >> 1;
		  done = 0;
		 }
	       else
		 {if (t == (values_m1))
		    {*(info_ptr++) = values_m1 - begin;
		     *(info_ptr++) = 0; 
		     *(info_ptr++) = 0; 
		     *(info_ptr++) = values;  
		     *(bin_ptr++) = B_ptr + values_m1;
		     l = r = t = begin = values_m1;
		     done = 0;
		    }   
		  else        
		    {*(info_ptr++) = t - begin;
		     *(info_ptr++) = 0; 
		     *(info_ptr++) = 0; 
		     *(info_ptr++) = t + 1; 
		     *(bin_ptr++) = B_ptr + t + 1;
		     l = begin = t;
		     r = values_m1; 
		     t = (l + r) >> 1;
		     done = 0;
		    }
		 }
	      }
	    else 
	      {if (tval == (*B_ptr))
		 {*(info_ptr++) = 0; 
		  l = 0;
		 }
	       else
		 {*(info_ptr++) = t - begin;  
		  l = (begin = t) + 1;
		 }
	       if (((r = t + 4) < values) && 
					  ((*(B_ptr + r)) > tval));
	       else   
		 r = values_m1;
	       done = 0;
	       t = (l + r) >> 1;
	       while (! done)
		 {if ((*(B_ptr + t)) > tval)
		    {r = t - 1;
		     t = (l + r) >> 1;
		    } 
		  else  
		    {if ((t == values_m1) || ((*(B_ptr + t + 1)) > tval))
		       done = 1;
		     else
		       {l = t + 1;
			t = (l + r) >> 1;
		       }
		    }
		 }
	       pre_t = t;  
	       rem = t - begin;
	       cur_c = (split_ptr + j)->count;
	       done = 0;
	       if (cur_c >= rem)
		 {*(info_ptr++) = rem; 
		  (split_ptr + j)->count -= rem;
		  if (t == values_m1)
		    {l = r = begin = values_m1;
		     *(bin_ptr++) = B_ptr + values_m1;
		     *(info_ptr++) = 0;
		     *(info_ptr++) = values; 
		    } 
		  else  
		    {l = begin = t;
		     r = values_m1;
		     t = (l + r) >> 1;
		     *(bin_ptr++) = B_ptr + l + 1;
		     *(info_ptr++) = 0;   
		     *(info_ptr++) = l + 1; 
		    } 
		 }
	       else 
		 {*(info_ptr++) = cur_c;
		  (split_ptr + j)->count = 0;
		  if (tval == ((split_ptr + j + 1)->val))
		    t = begin + cur_c;
		  done = 0;
		  if (t == values_m1)
		    {if (tval != ((split_ptr + j + 1)->val))
		       {*(info_ptr++) = rem - cur_c; 
			*(bin_ptr++) = B_ptr + begin + cur_c + 1; 
		       }
		     else
		       {*(info_ptr++) = 0;   
			*(bin_ptr++) = B_ptr + values_m1; 
		       }
		     l = r = begin = values_m1;
		     *(info_ptr++) = values; 
		    } 
		  else  
		    {if (tval != ((split_ptr + j + 1)->val))
		       {*(info_ptr++) = rem - cur_c;   
			*(bin_ptr++) = B_ptr + begin + cur_c + 1; 
		       }
		     else
		       {*(info_ptr++) = 0;  
			*(bin_ptr++) = B_ptr + t + 1; 
		       }
		     l = begin = t;
		     r = values_m1;
		     t = (l + r) >> 1;
		     *(info_ptr++) = pre_t + 1; 
		    }
		 }
	      }
	   }

	}
     }          
  }             

barrier();


/***************************************************************************
			       Step (7)
***************************************************************************/


#if (defined(SP2))

for (i=0;i<PROCS;i++)
  *(Bin_Count + i) = 0;
info_ptr = Info - 4;

for (i=0;i<PROCS;i++)
  for (j=0;j<PROCS;j++)
    {info_ptr += 4;
     *(Bin_Count + j) += (*info_ptr + *(info_ptr + 2) + *(info_ptr + 3));
    }

max = *(Bin_Count);

for (i=1;i<PROCS;i++)
  if ((*(Bin_Count + i)) > max)
     max = (*(Bin_Count + i));
  
barrier();

max = all_reduce_to_all_max(max);

C = finish_sptr;
max += (2 + PROCS_f);
 
finish_ptr = ((int *) (C + MYPROC)) - max;

for (i=0;i<PROCS;i++)
  {bin_ptr = (Bin + i - PROCS);
   info_ptr = (Info + (i << 2) - PROCS_f);
   C_ptr = (finish_ptr += max);
   for (j=0;j<PROCS;j++)
     {B_ptr = *(bin_ptr += PROCS);  
      info_ptr += PROCS_f;
      *(C_ptr + 1) = *(info_ptr + 1); 
      entries = ((*(C_ptr) = *(info_ptr)) + 
		(*(C_ptr + 2) = *(info_ptr + 2)) + 
		(*(C_ptr + 3) = *(info_ptr + 3)));  
      C_ptr += 4; 
      for (k=0;k<entries;k++)
	*(C_ptr + k) = *(B_ptr + k);
      C_ptr += entries;
     }
  }    
  
D = start_sptr;
C_ptr = ((int *) (C + MYPROC));
D_ptr = ((int *) (D + MYPROC));

barrier();

mpc_index(C_ptr,D_ptr,(max*sizeof(int)),ALLGRP);


barrier();

D = finish_sptr;
C = start_sptr;

#else

max = over;

B = start_sptr;
C = finish_sptr;
t = max*MYPROC;

bin_ptr = (Bin + MYPROC - PROCS);
info_ptr = (Info + (MYPROC << 2) - PROCS_f);
values = 0;
C_ptr = ((int *) (C + MYPROC)) + t;

for (i=0;i<PROCS;i++) 
  {B_ptr = *(bin_ptr += PROCS);
   info_ptr += PROCS_f;
   *(C_ptr + 1) = *(info_ptr + 1); 
   entries = ((*(C_ptr) = *(info_ptr)) + (*(C_ptr + 2) = *(info_ptr + 2)) + 
					 (*(C_ptr + 3) = *(info_ptr + 3)));  
   C_ptr += 4; 
   for (j=0;j<entries;j++)
     *(C_ptr + j) = *(B_ptr + j);
   C_ptr += entries;
   values+=entries; 
  }

for (i=1;i<PROCS;i++)
  {dest = (MYPROC + i) & PROCS_m1;
   dest_ptr = ((int *global) (C + dest)) + t;
   bin_ptr = (Bin + dest - PROCS);
   info_ptr = (Info + (dest << 2) - PROCS_f);
   for (j=0;j<PROCS;j++)
     {bulk_store(dest_ptr,(info_ptr += PROCS_f),
				     (4*(sizeof(int))));
      entries = (*(info_ptr)) + (*(info_ptr + 2)) + (*(info_ptr + 3));  
      bulk_store((dest_ptr + 4),(*(bin_ptr += PROCS)),
				     (entries*(sizeof(int))));
      dest_ptr += (entries + 4);
      values+=entries;
     }
   all_store_sync();
  }

C = finish_sptr;
D = start_sptr;

barrier();

#endif


/***************************************************************************
			       Step (8)
***************************************************************************/


C_ptr = ((int *) (C + MYPROC)) - max;
finish_ptr = ((int *) (D + MYPROC));
split_ptr = ((split_t *) (Splitters + MYPROC));  

for (i=0;i<PROCS;i++)
  *(Data + i) = (*(Bin + i) = (C_ptr += max)) + 4;

for (i=0;i<PROCS;i++)
  {entries = 0;
     for (j=0;j<PROCS;j++)
       if (times = *((*(Bin + j)))) 
	 {temp_ptr = finish_ptr + times;
	  while (finish_ptr < temp_ptr)
	    *(finish_ptr ++) = *((*(Data + j)) ++);
	  entries += times;
	 }
      start_row = *((*Bin) + 1);
      entries += (temp2 = *((*Bin) + 2));
      l_start_col = l_finish_col = 0;
      if (temp2 > 1)
	{finish_row = start_row + temp2 - 1;
	 r_start_col = r_finish_col = 1;
	}
      else
	{if (temp2)
	   {finish_row = start_row;
	    r_finish_col = 1;
	    r_start_col = 0;
	   }
	 else
	   {finish_row = start_row;
	    r_finish_col = 0;
	    r_start_col = 0;
	   }
	}   
      for (j=1;j<PROCS;j++)
	{temp1 = *((*(Bin + j)) + 1);
	 entries += (temp2 = *((*(Bin + j)) + 2));
	 if (temp2)
	   {if (((temp1 + temp2 - 1) == finish_row) && (r_finish_col))
	      r_finish_col++;
	    if ((temp1 == start_row)  && (r_start_col))
	      r_start_col++;
	    if (temp1 < start_row)
	      {start_row = temp1;
	       l_start_col = j;
	       r_start_col = j + 1;
	      }
	   }
	} 
      for (j=l_start_col;j<r_start_col;j++)
	*(finish_ptr++) = *((*(Data + j))++); 
      if ((times = (finish_row - start_row - 1)) > 0)
	{for (k=0;k<PROCS;k++)
	   {temp_ptr = finish_ptr + k;       
	    s_ptr = *(Data + k);
	    for (t=0;t<times;t++)
	      *(temp_ptr + (t << PROCSLOG)) = *(s_ptr + t);
	    (*(Data + k)) += times;
	   }
	 finish_ptr += (PROCS*times);
	}
      for (j=l_finish_col;j<r_finish_col;j++)
	*(finish_ptr++) = *((*(Data + j))++); 
      
   for (j=0;j<PROCS;j++)
     if (times = *((*(Bin + j)) + 3)) 
       {temp_ptr = finish_ptr + times;
	while (finish_ptr < temp_ptr)
	  *(finish_ptr ++) = *((*(Data + j)) ++);
	entries += times;
       }
   finish_ptr++;
   *(Number + i) = entries;   
   for (j=0;j<PROCS;j++)     
     *(Data + j) = (*(Bin + j) = *(Data + j)) + 4;
  }

barrier();


/***************************************************************************
			       Step (9)
***************************************************************************/


start_sptr = D;
finish_sptr =  C;

adjustments = 0;

start_ptr = ((int *) (start_sptr + MYPROC));
finish_ptr = ((int *) (finish_sptr + MYPROC));
for (j=0;j<PROCS;j+=2)
  {tot1 = *(Number + j);
   s2_ptr = (s1_ptr = start_ptr) + tot1 + 1;
   tot2 = *(Number + j + 1);
   (*(Number + (j >> 1))) = tot = tot1 + tot2;
   (*(s1_ptr + tot1)) = (*(s2_ptr + tot2)) = MAX_VAL;
   v = tot1 - 1;
   if (tot1)
     while ((*(s1_ptr + v)) == MAX_VAL)
       {adjustments++;
	(*(s1_ptr + (v--))) = MAX_VAL_m1; 
       } 
   v = tot2 - 1;
   if (tot2)
     while ((*(s2_ptr + v)) == MAX_VAL)
       {adjustments++;
	(*(s2_ptr + (v--))) = MAX_VAL_m1; 
       } 
   temp_ptr = finish_ptr + tot;
   while (finish_ptr < temp_ptr)
     {if ((*s1_ptr) < (*s2_ptr))
	*(finish_ptr++) = *(s1_ptr++);
      else 
	*(finish_ptr++) = *(s2_ptr++);
     } 
   start_ptr += (tot+2);
   *(finish_ptr++) = MAX_VAL;  
  }

if (PROCSLOG > 1)   
  {temp_sptr = finish_sptr;  
   finish_sptr = start_sptr;
   start_sptr = temp_sptr;
  }

for (i=1;i<PROCSLOG;i++)
  {start_ptr = ((int *) (start_sptr + MYPROC));
   finish_ptr = ((int *) (finish_sptr + MYPROC));
   for (j=0;j<(PROCS >> i);j+=2)
     {tot1 = (*(Number+j));
      tot2 = (*(Number+j+1));
      s2_ptr = (s1_ptr = start_ptr) + tot1 + 1;
      (*(Number+(j>>1))) = tot = tot1 + tot2;
      temp_ptr = finish_ptr + tot;
      while (finish_ptr < temp_ptr)
	{if ((*s1_ptr) < (*s2_ptr))
	   *(finish_ptr++) = *(s1_ptr++);
	 else 
	   *(finish_ptr++) = *(s2_ptr++);
	} 
      start_ptr += (tot + 2);
      *(finish_ptr++) = MAX_VAL;  
     }
   if (i < (PROCSLOG - 1))
     {temp_sptr = finish_sptr;
      finish_sptr = start_sptr;
      start_sptr = temp_sptr;
     }
  }

finish_ptr = ((int *) (finish_sptr + MYPROC)); 
values = v = *Number;
for (i=0;i<adjustments;i++)  
  *(finish_ptr + (v--)) = MAX_VAL;

*val_ptr = values;
*out_array_ptr = finish_sptr;

barrier();
 
free(Bin);
free(Bin_Count);
free(Info);
free(Data);
free(Loc);
free(Number);
free(Count1);
free(Count2);
free(Sample1);
free(Sample2);

all_spread_free(Splitters);

}
