

/*                                                                    
 *
 * ran_i_sort.sc - Function to Sort Integers by a Novel Variation of 
 *                 Sample Sort (Split-C Code)
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
 * Filename:            ran_i_sort.sc
 * History:
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <split-c/split-c.h>
#include <split-c/control.h>

/* The largest allowable integer value: */

#define MAX_VAL 2147483647 

#define MAX_VAL_m1 MAX_VAL-1

/* See the introduction below for an explanation of BUNDLE: */ 

#define BUNDLE 4

#define PROCS_m1 (PROCS-1)


typedef struct record
  {int counts;
   struct record *next; 
   int *bin;
  } record_t;

typedef struct split
  {int val;
   double frac;
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

/*****************************************************************************
 *                           RRANDOM()
 *
 *               (The 4.2 BSD version of random( ))
 *
 * Copyright (c) 1983, 1993
 *	The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *	This product includes software developed by the University of
 *	California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * An improved random number generation package.  In addition to the standard
 * rand()/srand() like interface, this package also has a special state info
 * interface.  The iinitstate() routine is called with a seed, an array of
 * bytes, and a count of how many bytes are being passed in; this array is
 * then initialized to contain information for random number generation with
 * that much state information.  Good sizes for the amount of state
 * information are 32, 64, 128, and 256 bytes.  The state can be switched by
 * calling the ssetstate() routine with the same array as was initiallized
 * with iinitstate().  By default, the package runs with 128 bytes of state
 * information and generates far better random numbers than a linear
 * congruential generator.  If the amount of state information is less than
 * 32 bytes, a simple linear congruential R.N.G. is used.
 *
 * Internally, the state information is treated as an array of longs; the
 * zeroeth element of the array is the type of R.N.G. being used (small
 * integer); the remainder of the array is the state information for the
 * R.N.G.  Thus, 32 bytes of state information will give 7 longs worth of
 * state information, which will allow a degree seven polynomial.  (Note:
 * the zeroeth word of state information also has some other information
 * stored in it -- see ssetstate() for details).
 * 
 * The random number generation technique is a linear feedback shift register
 * approach, employing trinomials (since there are fewer terms to sum up that
 * way).  In this approach, the least significant bit of all the numbers in
 * the state table will act as a linear feedback shift register, and will
 * have period 2^deg - 1 (where deg is the degree of the polynomial being
 * used, assuming that the polynomial is irreducible and primitive).  The
 * higher order bits will have longer periods, since their values are also
 * influenced by pseudo-random carries out of the lower bits.  The total
 * period of the generator is approximately deg*(2**deg - 1); thus doubling
 * the amount of state information has a vast influence on the period of the
 * generator.  Note: the deg*(2**deg - 1) is an approximation only good for
 * large deg, when the period of the shift register is the dominant factor.
 * With deg equal to seven, the period is actually much longer than the
 * 7*(2**7 - 1) predicted by this formula.
 */

/*
 * For each of the currently supported random number generators, we have a
 * break value on the amount of state information (you need at least this
 * many bytes of state info to support this random number generator), a degree
 * for the polynomial (actually a trinomial) that the R.N.G. is based on, and
 * the separation between the two lower order coefficients of the trinomial.
 */
#define TYPE_0          0               /* linear congruential */
#define BREAK_0         8
#define DEG_0           0
#define SEP_0           0

#define TYPE_1          1               /* x**7 + x**3 + 1 */
#define BREAK_1         32
#define DEG_1           7
#define SEP_1           3

#define TYPE_2          2               /* x**15 + x + 1 */
#define BREAK_2         64
#define DEG_2           15
#define SEP_2           1

#define TYPE_3          3               /* x**31 + x**3 + 1 */
#define BREAK_3         128
#define DEG_3           31
#define SEP_3           3

#define TYPE_4          4               /* x**63 + x + 1 */
#define BREAK_4         256
#define DEG_4           63
#define SEP_4           1

/*
 * Array versions of the above information to make code run faster --
 * relies on fact that TYPE_i == i.
 */
#define MAX_TYPES       5               /* max number of types above */

static int degrees[MAX_TYPES] = { DEG_0, DEG_1, DEG_2, DEG_3, DEG_4 };
static int seps [MAX_TYPES] =   { SEP_0, SEP_1, SEP_2, SEP_3, SEP_4 };

/*
 * Initially, everything is set up as if from:
 *
 *      iinitstate(1, &randtbl, 128);
 *
 * Note that this initialization takes advantage of the fact that srandom()
 * advances the front and rear pointers 10*rand_deg times, and hence the
 * rear pointer which starts at 0 will also end up at zero; thus the zeroeth
 * element of the state information, which contains info about the current
 * position of the rear pointer is just
 *
 *      MAX_TYPES * (rptr - state) + TYPE_3 == TYPE_3.
 */

static long randtbl[DEG_3 + 1] = {
	TYPE_3,
	0x9a319039, 0x32d9c024, 0x9b663182, 0x5da1f342, 0xde3b81e0, 0xdf0a6fb5,
	0xf103bc02, 0x48f340fb, 0x7449e56b, 0xbeb1dbb0, 0xab5c5918, 0x946554fd,
	0x8c2e680f, 0xeb3d799f, 0xb11ee0b7, 0x2d436b86, 0xda672e2a, 0x1588ca88,
	0xe369735d, 0x904f35f7, 0xd7158fd6, 0x6fa6f051, 0x616e6b96, 0xac94efdc,
	0x36413f93, 0xc622c298, 0xf5a42ab8, 0x8a88d77b, 0xf5ad9d0e, 0x8999220b,
	0x27fb47b9,
};

/*
 * fptr and rptr are two pointers into the state info, a front and a rear
 * pointer.  These two pointers are always rand_sep places aparts, as they
 * cycle cyclically through the state information.  (Yes, this does mean we
 * could get away with just one pointer, but the code for random() is more
 * efficient this way).  The pointers are left positioned as they would be
 * from the call
 *
 *      iinitstate(1, randtbl, 128);
 *
 * (The position of the rear pointer, rptr, is really 0 (as explained above
 * in the initialization of randtbl) because the state table pointer is set
 * to point to randtbl[1] (as explained below).
 */
static long *fptr = &randtbl[SEP_3 + 1];
static long *rptr = &randtbl[1];

/*
 * The following things are the pointer to the state information table, the
 * type of the current generator, the degree of the current polynomial being
 * used, and the separation between the two pointers.  Note that for efficiency
 * of random(), we remember the first location of the state information, not
 * the zeroeth.  Hence it is valid to access state[-1], which is used to
 * store the type of the R.N.G.  Also, we remember the last location, since
 * this is more efficient than indexing every time to find the address of
 * the last element to see if the front and rear pointers have wrapped.
 */
static long *state = &randtbl[1];
static int rand_type = TYPE_3;
static int rand_deg = DEG_3;
static int rand_sep = SEP_3;
static long *end_ptr = &randtbl[DEG_3 + 1];


/*
 * random:
 *
 * If we are using the trivial TYPE_0 R.N.G., just do the old linear
 * congruential bit.  Otherwise, we do our fancy trinomial stuff, which is
 * the same in all the other cases due to all the global variables that have
 * been set up.  The basic operation is to add the number at the rear pointer
 * into the one at the front pointer.  Then both pointers are advanced to
 * the next location cyclically in the table.  The value returned is the sum
 * generated, reduced to 31 bits by throwing away the "least random" low bit.
 *
 * Note: the code takes advantage of the fact that both the front and
 * rear pointers can't wrap on the same call by not testing the rear
 * pointer if the front one has wrapped.
 *
 * Returns a 31-bit random number.
 */
inline long
rrandom()
{
	long i;

		*fptr += *rptr;
		i = (*fptr >> 1) & 0x7fffffff;  /* chucking least random bit */
		if (++fptr >= end_ptr) {
			fptr = state;
			++rptr;
		} else if (++rptr >= end_ptr)
			rptr = state;

	return(i);
}


/*
 * srandom:
 *
 * Initialize the random number generator based on the given seed.  If the
 * type is the trivial no-state-information type, just remember the seed.
 * Otherwise, initializes state[] based on the given "seed" via a linear
 * congruential generator.  Then, the pointers are set to known locations
 * that are exactly rand_sep places apart.  Lastly, it cycles the state
 * information a given number of times to get rid of any initial dependencies
 * introduced by the L.C.R.N.G.  Note that the initialization of randtbl[]
 * for default usage relies on values produced by this routine.
 */

inline void 
srrandom(x)
	unsigned int x;
{
	register int i, j;

	if (rand_type == TYPE_0)
		state[0] = x;
	else {
		j = 1;
		state[0] = x;
		for (i = 1; i < rand_deg; i++)
			state[i] = 1103515245 * state[i - 1] + 12345;
		fptr = &state[rand_sep];
		rptr = &state[0];
		for (i = 0; i < 10 * rand_deg; i++)
			(void)rrandom();
	}
}

/*
 * iinitstate:
 *
 * Initialize the state information in the given array of n bytes for future
 * random number generation.  Based on the number of bytes we are given, and
 * the break values for the different R.N.G.'s, we choose the best (largest)
 * one we can and set things up for it.  srandom() is then called to
 * initialize the state information.
 * 
 * Note that on return from srandom(), we set state[-1] to be the type
 * multiplexed with the current value of the rear pointer; this is so
 * successive calls to iinitstate() won't lose this information and will be
 * able to restart with ssetstate().
 * 
 * Note: the first thing we do is save the current state, if any, just like
 * ssetstate() so that it doesn't matter when iinitstate is called.
 *
 * Returns a pointer to the old state.
 */
char *
iinitstate(seed, arg_state, n)
	unsigned int seed;              /* seed for R.N.G. */
	char *arg_state;                /* pointer to state array */
	int n;                          /* # bytes of state info */
{
	register char *ostate = (char *)(&state[-1]);

	if (rand_type == TYPE_0)
		state[-1] = rand_type;
	else
		state[-1] = MAX_TYPES * (rptr - state) + rand_type;
	if (n < BREAK_0) {
		(void)fprintf(stderr,
		    "rrandom: not enough state (%d bytes); ignored.\n", n);
		return(0);
	}
	if (n < BREAK_1) {
		rand_type = TYPE_0;
		rand_deg = DEG_0;
		rand_sep = SEP_0;
	} else if (n < BREAK_2) {
		rand_type = TYPE_1;
		rand_deg = DEG_1;
		rand_sep = SEP_1;
	} else if (n < BREAK_3) {
		rand_type = TYPE_2;
		rand_deg = DEG_2;
		rand_sep = SEP_2;
	} else if (n < BREAK_4) {
		rand_type = TYPE_3;
		rand_deg = DEG_3;
		rand_sep = SEP_3;
	} else {
		rand_type = TYPE_4;
		rand_deg = DEG_4;
		rand_sep = SEP_4;
	}
	state = &(((long *)arg_state)[1]);      /* first location */
	end_ptr = &state[rand_deg];     /* must set end_ptr before srandom */
	srrandom(seed);
	if (rand_type == TYPE_0)
		state[-1] = rand_type;
	else
		state[-1] = MAX_TYPES*(rptr - state) + rand_type;
	return(ostate);
}

/*
 * ssetstate:
 *
 * Restore the state from the given state array.
 *
 * Note: it is important that we also remember the locations of the pointers
 * in the current state information, and restore the locations of the pointers
 * from the old state information.  This is done by multiplexing the pointer
 * location into the zeroeth word of the state information.
 *
 * Note that due to the order in which things are done, it is OK to call
 * ssetstate() with the same state as the current state.
 *
 * Returns a pointer to the old state information.
 */
char *
ssetstate(arg_state)
	char *arg_state;
{
	register long *new_state = (long *)arg_state;
	register int type = new_state[0] % MAX_TYPES;
	register int rear = new_state[0] / MAX_TYPES;
	char *ostate = (char *)(&state[-1]);

	if (rand_type == TYPE_0)
		state[-1] = rand_type;
	else
		state[-1] = MAX_TYPES * (rptr - state) + rand_type;
	switch(type) {
	case TYPE_0:
	case TYPE_1:
	case TYPE_2:
	case TYPE_3:
	case TYPE_4:
		rand_type = type;
		rand_deg = degrees[type];
		rand_sep = seps[type];
		break;
	default:
		(void)fprintf(stderr,
		    "random: state info corrupted; not changed.\n");
	}
	state = &new_state[1];
	if (rand_type != TYPE_0) {
		rptr = &state[rear];
		fptr = &state[(rear + rand_sep) % rand_deg];
	}
	end_ptr = &state[rand_deg];             /* set end_ptr too */
	return(ostate);
}


/**************************************************************************
*
*                     ALL_RANDOM_I_SORT
*
*
*      The function all_random_i_sort is called with the following
* parameters.
*
*   (1) col_size -> the number of elements per processor in the input array, 
* 
*   (2) bits -> the number of bits used to represent the sorting radix  
*             (i.e. - the sorting radix will be 2^(bits)).
*
*   (3) over -> (over*col_size/PROCS) forms the estimate of the maximum 
*               number of elements exchanged by any two processors in 
*               either Step(2) or Step (7).  
*
*   (4) Buffer1 & 
*   (5) Buffer2 ->  spread pointers to the base of two spread arrays at 
*                   processor P_0.  When all_random_i_sort is called, 
*                   Buffer1 holds the input array to be sorted.  As sorting 
*                   progresses, the data is moved back and forth between 
*                   the two arrays until at completion either array may 
*                   hold the sorted output.  Thus, either array must be large enough 
*                   to hold the maximum number of elements possible
*                   at a processor at any point during sorting, which we
*                   estimate to be (over*col_size) elements.
*                               
*   (6) *out_array_ptr -> a local pointer to a global pointer to the
*                         array which holds the completely sorted output.
*
*   (7) *val_ptr-> a local pointer to the variable which describes how
*                  many items are present at that processor from the          
*                  completely sorted output.
*
*   (8) value_range -> the largest possible integer value in the set to
*                      be sorted, where `value_range' can not exceed  
*                      MAX_VAL.
*
*   Finally, the preprocessor definition BUNDLE needs to be explained.
*   The idea is that in Step (1), as we place the input elements into 
*   buckets, it really isn't necessary to randomly place each and every
*   element.   Rather, we might place the elements into bucket in bundles
*   of size BUNDLE.  The larger we make BUNDLE, the less calls we need to 
*   random number generator rrandom().  On the other hand, if BUNDLE is too
*   large, the performance will begin to suffer because of poor load
*   balance.  We have found that setting BUNDLE to 4 seems to work nicely,
*   but the user can certainly reset this value as he chooses.   
*
***************************************************************************/



void all_random_i_sort(int col_size, int bits, int over, int *spread Buffer1, 
		     int *spread Buffer2, int *spread *out_array_ptr,
		     int *val_ptr, int value_range)

{

int buck_size = (over*col_size/PROCS);

register int i,j,k,t,entries,values,v,times,stride,stride_p1,stride_m1,result,
             rem,dest,tot,tot1,tot2,adjustments,max,done,l,r,begin,digits,
             radix,radix_m1,shift,shift1,shift2,log_value_range,pre_begin,
             pre_r,size,values_m1,tval,ones;

int *in_ptr,*A_ptr,*B_ptr,*C_ptr,*temp_ptr, *E_ptr, *bin_ptr,*start_ptr,
    *finish_ptr, *s1_ptr, *s2_ptr, *D_ptr,*Count,*Next,*Offset,*t_ptr,
    *c_ptr,*s_ptr, *record_ptr,**Bin, *Bin_Count, *Number, *Count1, 
    *Count2, **Loc;


int *spread start_sptr;
int *spread finish_sptr;
int *spread temp_sptr;
int *spread A;
int *spread B;
int *spread C;
int *spread D;
int *spread E;

int *global dest_ptr;

split_t *split_ptr;

split_t *spread Splitters;


Bin = (int **) malloc(PROCS*sizeof(int *));
assert_malloc(Bin);

Bin_Count = (int *) malloc(PROCS*sizeof(int));
assert_malloc(Bin_Count);

Number = (int *) malloc(PROCS*sizeof(int));
assert_malloc(Number);

Count1 = (int *) malloc(((1 << bits)+1)*sizeof(int));
assert_malloc(Count1);

Count2 = (int *) malloc(((1 << bits)+1)*sizeof(int));
assert_malloc(Count2);

Loc = (int **) malloc(PROCS*sizeof(int *));
assert_malloc(Loc);

Splitters = all_spread_malloc(PROCS,PROCS*sizeof(split_t));
assert_spread_malloc(Splitters);


/***************************************************************************
			       Step (1)
***************************************************************************/


A = Buffer2;
in_ptr = ((int *) (Buffer1 + MYPROC));
A_ptr = ((int *) (A+MYPROC)) - buck_size - BUNDLE + 1;
srrandom(21+1001*MYPROC);

for (i=0;i<PROCS;i++)
  *(Loc + i) = (A_ptr += buck_size);

#if BUNDLE == 1

times = (int) 31/PROCSLOG;
shift = 31 - times*PROCSLOG;
result = col_size/times;
rem = col_size - times*result;
for (i=0;i<result;i++)
  {*(++(*(Loc + ((v = (rrandom() >> shift)) & PROCS_m1)))) = *(in_ptr++); 
   for (j=1;j<times;j++)
     *(++(*(Loc + ((v >>= PROCSLOG) & PROCS_m1)))) = *(in_ptr++); 
  }

v = rrandom();
for (i=0;i<rem;i++)
  *(++(*(Loc + ((v >>= PROCSLOG) & PROCS_m1)))) = *(in_ptr++); 

#endif

#if BUNDLE > 1

times = (int) (31)/PROCSLOG;
shift = 31 - times*PROCSLOG;
result = col_size/(times*BUNDLE);
rem = (col_size - BUNDLE*times*result)/BUNDLE;
for (i=0;i<result;i++)
  {record_ptr = (*(Loc + ((v = (rrandom() >> shift)) & PROCS_m1)) += BUNDLE); 
   for (k=0;k<BUNDLE;k++)
     *(record_ptr + k) = *(in_ptr++);
   for (j=1;j<times;j++)
     {record_ptr = (*(Loc + ((v >>= PROCSLOG) & PROCS_m1)) += BUNDLE);
      for (k=0;k<BUNDLE;k++)
        *(record_ptr + k) = *(in_ptr++);
     }
  }

v = rrandom();
for (i=0;i<rem;i++)
  {record_ptr = (*(Loc + ((v >>= PROCSLOG) & PROCS_m1)) += BUNDLE); 
   for (k=0;k<BUNDLE;k++)
     *(record_ptr + k) = *(in_ptr++);
  }
  
#endif

A_ptr = ((int *) (A + MYPROC)); 
max = *(Bin_Count) = *(Loc) + BUNDLE - 1 - A_ptr; 
for (i=1;i<PROCS;i++)
  if (((*(Bin_Count + i)) = (*(Loc+i) + BUNDLE - 1 - (A_ptr+=buck_size)))
	 > max)
       max = (*(Bin_Count + i)); 

max = all_reduce_to_all_max(max);

if ((max-2) > buck_size)
  {fprintf(stderr,"ERROR: Too many elements to route in Step(7)\n");
   fflush(stderr);
   exit(1);
  }

barrier();


/***************************************************************************
			       Step (2)
***************************************************************************/


#if (defined(SP2))

B = Buffer1;

max += 2;
start_ptr = ((int *) (A + MYPROC)) - buck_size;
finish_ptr = ((int *) (B + MYPROC)) - max; 

for (i=0;i<PROCS;i++)
  {entries = *(finish_ptr += max) = *(Bin_Count + i);
   start_ptr += buck_size;
   for (j=1;j<=entries;j++)
     *(finish_ptr+j) = *(start_ptr+j);
  }

C = Buffer2;
B_ptr = ((int *) (B + MYPROC));
C_ptr = ((int *) (C + MYPROC));


barrier();

mpc_index(B_ptr,C_ptr,(max*sizeof(int)),ALLGRP);

barrier();

start_sptr = Buffer2;
finish_sptr = Buffer1;


#else

B = Buffer1;
A_ptr = ((int *) (A + MYPROC)); 
t = buck_size*MYPROC;
B_ptr = (((int *) (B + MYPROC)) + t);
*B_ptr = entries = *(Bin_Count);

for (i=1;i<=entries;i++)
  *(B_ptr + i) = *(A_ptr + i);

for (i=1;i<PROCS;i++)
  {dest = (MYPROC + i) & PROCS_m1;
   dest_ptr = ((int *global) (B + dest)) + t;
   *(A_ptr += buck_size) = (*(Bin_Count + i)); 
   bulk_store(dest_ptr,A_ptr,((*(Bin_Count + i) + 1)*(sizeof(int))));
   barrier();
  }

all_store_sync();

start_sptr = Buffer1;
finish_sptr = Buffer2;
max = buck_size;

barrier();

#endif


/***************************************************************************
			       Step (3)
***************************************************************************/


radix = (1 << bits);  
radix_m1 = radix - 1;
values = 0;

log_value_range = (int) ceil(log((double) value_range)/log((double) 2));

digits = (int) ceil(((double) log_value_range)/((double) bits));

rem = log_value_range % bits;

if (rem == 0)
  rem = bits;

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
   s_ptr = start_ptr - max;
   for (j=0;j<PROCS;j++)
     {values += (entries = *(c_ptr = (s_ptr += max)));
      t_ptr = c_ptr + entries;
      while ((++c_ptr) <= t_ptr)
	(*(Count + ((*c_ptr) & radix_m1)))++;
     }
   *(Offset) = 0;
   t_ptr = (c_ptr = Offset) + radix;
   while ((++c_ptr) < t_ptr)
     *c_ptr +=  *(c_ptr - 1);
   s_ptr = start_ptr - max;
   for (j=0;j<PROCS;j++)
     {entries = *(c_ptr = (s_ptr += max));
      t_ptr = c_ptr + entries;
      while ((++c_ptr) <= t_ptr)
	{(*(Next + (((v = *c_ptr) >> bits) & radix_m1)))++;
	 *(finish_ptr + ((*(Offset + (v & radix_m1)))++)) = v;
	}
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
   
   radix_m1 = (radix = (1 << rem)) - 1;
   shift1 = (digits - 1) * bits;
   *(Offset) = 0;
   t_ptr = (c_ptr = Offset) + radix;
   while ((++c_ptr) < t_ptr)
     *c_ptr +=  *(c_ptr - 1);
   t_ptr = (c_ptr = start_ptr) + values;
   while (c_ptr < t_ptr)
     {v = *(c_ptr++);
      *(finish_ptr + ((*(Offset + ((v >> shift1) & radix_m1)))++)) = v;
     }
  } 
else
  {radix_m1 = (radix = (1 << rem)) - 1;
   t_ptr = (c_ptr = Count) + radix;
   while (c_ptr < t_ptr)
     *(c_ptr++) = 0;
   s_ptr = start_ptr - max;
   for (j=0;j<PROCS;j++)
     {values += (entries = *(c_ptr = (s_ptr += max)));
      t_ptr = c_ptr + entries;
      while ((++c_ptr) <= t_ptr)
	(*(Count + ((*c_ptr) & radix_m1)))++;
     }
   *(Offset) = 0;
   t_ptr = (c_ptr = Offset) + radix;
   while ((++c_ptr) < t_ptr)
     *c_ptr +=  *(c_ptr - 1);
   s_ptr = start_ptr - max;
   for (j=0;j<PROCS;j++)
     {entries = *(c_ptr = (s_ptr += max));
      t_ptr = c_ptr + entries;
      while ((++c_ptr) <= t_ptr)
	*(finish_ptr + ((*(Offset + ((v = *c_ptr) & radix_m1)))++)) = v;
     }
  }
  
C = finish_sptr;

barrier();


/***************************************************************************
			       Step (4)
***************************************************************************/


if (MYPROC == 0)
  {C_ptr = ((int *) C);
   stride = (int) (values/PROCS);
   rem = values - PROCS*stride;
   stride_p1 = stride + 1;
   stride_m1 = stride - 1;
   split_ptr = ((split_t *) Splitters);
   values_m1 = values - 1;
   tot = -1;
   for (i=0;i<rem;i++)
     {tval = (split_ptr + i)->val = *(C_ptr + (tot += stride_p1));
      if ((!(i)) || ((split_ptr + i - 1)->val != tval))
	{if ((tval > (*(C_ptr + tot - 1))) && (tval < (*(C_ptr + tot + 1))))
	   (split_ptr + i)->frac = 1.0;
	 else
	   {l = tot - stride;
	    r = tot;
	    t = (l + r) >> 1;
	    done = 0;
	    while (! done)
	      {if (tval > (*(C_ptr + t)))
		 {l = t + 1;
		  t = (l + r) >> 1;
		 }
	       else
		 {if ((! t) || (tval > (*(C_ptr + t - 1))))
		    {pre_begin = t - 1;
		     done = 1;
		    }
		  else
		    {r = t - 1;
		     t = (l + r) >> 1;
		    }
		 }
	      }
	    l = tot;
	    r = values_m1;
	    if ((tot + PROCS) < values)
	      t = tot + PROCS;
	    else
	      t = (l + r) >> 1;
	    done = 0;
	    while (! done)
	      {if ((*(C_ptr + t)) > tval)
		 {r = t - 1;
		  t = (l + r) >> 1;
		 } 
	       else  
		 {if ((t == values_m1) || ((*(C_ptr + t + 1)) > tval))
		    done = 1;
		  else
		    {l = t + 1;
		     t = (l + r) >> 1;
		    }
		 }
	      }
	    pre_r = t;  
	    size = pre_r - pre_begin;
	    (split_ptr + i)->frac = (((double) (tot - pre_begin))/ 
				     ((double) size)); 
	   }  
	}
      else      
	(split_ptr + i)->frac = (((double) stride_p1)/((double) size)); 
     }
   tot++;
   for (i=rem;i<PROCS_m1;i++)
     {tval = (split_ptr + i)->val = *(C_ptr + (tot += stride));
      if ((!(i)) || ((split_ptr + i - 1)->val != tval))
	{if ((tval > (*(C_ptr + tot - 1))) && (tval < (*((C_ptr + tot + 1)))))
	   (split_ptr + i)->frac = 1.0;
	 else
	   {l = tot - stride_m1;
	    r = tot;
	    t = (l + r) >> 1;
	    done = 0;
	    while (! done)
	      {if (tval > (*(C_ptr + t)))
		 {l = t + 1;
		  t = (l + r) >> 1;
		 }
	       else
		 {if ((! t) || (tval > (*(C_ptr + t - 1))))
		    {pre_begin = t - 1;
		     done = 1;
		    }
		  else
		    {r = t - 1;
		     t = (l + r) >> 1;
		    }
		 }
	      }
	    l = tot;
	    r = values_m1;
	    if ((tot + PROCS) < values)
	      t = tot + PROCS;
	    else
	      t = (l + r) >> 1;
	    done = 0;
	    while (! done)
	      {if ((*(C_ptr + t)) > tval)
		 {r = t - 1;
		  t = (l + r) >> 1;
		 } 
	       else  
		 {if ((t == values_m1) || ((*(C_ptr + t + 1)) > tval))
		    done = 1;
		  else
		    {l = t + 1;
		     t = (l + r) >> 1;
		    }
		 }
	      }
	    pre_r = t;  
	    size = pre_r - pre_begin;
	    (split_ptr + i)->frac = (((double) (tot - pre_begin))/ 
				     ((double) size)); 
	   }  
	}
      else      
	(split_ptr + i)->frac = (((double) stride)/((double) size)); 
     }
  }


/***************************************************************************
			       Step (5)
***************************************************************************/


barrier();

v = PROCS - 1 - MYPROC; 

ones = 0;

while ((v & 1))
  {ones++;
   v >>= 1;
  }

#if (defined(T3D))

ones--;

for (i=PROCSLOG-1;i>=0;i--)
  {if (ones >= i)
      bulk_put(((split_t *global) (Splitters + MYPROC + (1 << i))),
	     ((split_t *) (Splitters + MYPROC)),(PROCS_m1 * (sizeof(split_t))));
   sync();
   barrier();
  }

#else


if (MYPROC > 0)
  store_sync(PROCS_m1*(sizeof(split_t)));

reset_store_sync();

for (i=(ones-1);i>=0;i--)
  bulk_store(((split_t *global) (Splitters + MYPROC + (1 << i))),
	     ((split_t *) (Splitters + MYPROC)),(PROCS_m1 * (sizeof(split_t))));

barrier();
	    
#endif


/***************************************************************************
			       Step (6)
***************************************************************************/


C_ptr = ((int *) (C + MYPROC));
split_ptr = ((split_t *) (Splitters + MYPROC));
(split_ptr + PROCS_m1)->val = MAX_VAL;
values_m1 = values - 1;

begin = -1;
l = 0;
r = values_m1; 
if (PROCS > 2)
  t = shift = 2 * ((int) values/PROCS);
else
  t = shift = values_m1;
done = 0;
*Bin = C_ptr;

if (MYPROC)
  {for (i=0;i<PROCS;i++)
     {if (i == PROCS_m1)
	*(Bin_Count + i) = values_m1 - begin;
      else
	{tval = (split_ptr + i)->val; 
	 if (((i)) && (tval == ((split_ptr + i - 1)->val)))
	   {if ((pre_begin += ((int) ((split_ptr + i)->frac * 
			  ((double) size)))) <= pre_r);
	    else
	      pre_begin = pre_r;
	    *(Bin_Count + i) = pre_begin - begin;
	    done = 0;
	    if (pre_begin == values_m1)
	      {l = r = begin = values_m1;
	       *(Bin + i + 1) = C_ptr + values_m1;
	      } 
	    else  
	      {l = begin = pre_begin;
	       r = values_m1;
	       t = (l + r) >> 1;
	       *(Bin + i + 1) = C_ptr + pre_begin + 1;
	      } 
	   }
	 else  
	   {while (! done)
	      {if (tval <= (*(C_ptr + t)))
		 {if (! t)
		    done = 1;
		  else
		    {r = t - 1;
		     t = (l + r) >> 1;
		    }
		 }
	       else
		 {if ((t == (values_m1)) || (tval <= (*(C_ptr + t + 1)))) 
		    done = 1;
		  else
		    {l = t + 1;
		     t = (l + r) >> 1;
		    }
		 }
	      }
	    if ((t == values_m1) || ((t) && (tval < (*(C_ptr + t + 1)))) || 
		     ((!t) && ((tval < (*C_ptr)) || 
		     ((tval > (*C_ptr)) && (tval < (*(C_ptr + 1)))))))  
	      {if (tval < (*C_ptr))
		 {*(Bin_Count + i) = 0;  
		  *(Bin + i + 1) = C_ptr;
		  begin = -1;
		  l = 0;
		  r = values_m1;
		  t = shift;
		  done = 0;
		 }
	       else
		 {if (t == (values_m1))
		    {*(Bin_Count + i) = values_m1 - begin;  
		     *(Bin + i + 1) = C_ptr + values_m1;
		     l = r = t = begin = values_m1;
		     done = 0;
		    }   
		  else        
		    {*(Bin_Count + i) = t - begin;  
		     *(Bin + i + 1) = C_ptr + t + 1;
		     l = begin = t;
		     r = values_m1; 
		     if ((t += shift) < values_m1);
		     else
		       t = (l + r) >> 1;
		     done = 0;
		    }
		 }
	       size = 0;
	       pre_r = pre_begin = begin;
	      }
	    else  
	      {if (tval == (*C_ptr))
		 {pre_begin = -1;
		  l = 0;
		 }
	       else
		 {pre_begin = t;
		  l = t;
		 }
	       r = values_m1; 
	       if ((t += PROCS) < values);
	       else
		 t = (l + r) >> 1;
	       done = 0;
	       while (! done)
		 {if ((*(C_ptr + t)) > tval)
		    {r = t - 1;
		     t = (l + r) >> 1;
		    } 
		  else  
		    {if ((t == values_m1) || ((*(C_ptr + t + 1)) > tval))
		       done = 1;
		     else
		       {l = t + 1;
			t = (l + r) >> 1;
		       }
		    }
		 }
	       pre_r = t;  
	       size = pre_r - pre_begin;
	       if ((pre_begin += ((int) ((split_ptr + i)->frac * 
			  ((double) size)))) <= pre_r);
	       else
		 pre_begin = pre_r;
	       *(Bin_Count + i) = pre_begin - begin;
	       done = 0;
	       if (pre_begin == values_m1)
		 {l = r = begin = values_m1;
		  *(Bin + i + 1) = C_ptr + values_m1;
		 } 
	       else  
		 {l = begin = pre_begin;
		  r = values_m1;
		  t = (l + r) >> 1;
		  *(Bin + i + 1) = C_ptr + pre_begin + 1;
		 } 
	      }
	   }
	}
     }             
  }
else
  {stride = (int) (values/PROCS);
   rem = values - PROCS*stride;
   stride_p1 = stride + 1;
   temp_ptr = C_ptr - stride_p1;
   for (i=0;i<rem;i++)
     {*(Bin + i) =  (temp_ptr += stride_p1);
      *(Bin_Count + i) = stride_p1;
     }
   temp_ptr ++;
   for (i=rem;i<PROCS_m1;i++)
     {*(Bin + i) =  (temp_ptr += stride);
      *(Bin_Count + i) = stride;
     }
   if (rem == PROCS_m1)
     {*(Bin + i) =  (temp_ptr + stride_p1);
      *(Bin_Count + i) = stride;
     }
   else
     {*(Bin + i) =  (temp_ptr + stride);
      *(Bin_Count + i) = stride;
     }
  }   


/***************************************************************************
			       Step (7)
***************************************************************************/


#if (defined(SP2))

max = *(Bin_Count);

for (i=1;i<PROCS;i++)
  if ((*(Bin_Count + i)) > max)
     max = (*(Bin_Count + i));
  
barrier();

max = all_reduce_to_all_max(max);

if ((max-2) > buck_size)
  {fprintf(stderr,"ERROR: Too many elements to route in Step(7)\n");
   fflush(stderr);
   exit(1);
  }


D = start_sptr;
max += 2;
start_ptr = ((int *) (C + MYPROC));
finish_ptr = bin_ptr = ((int *) (D + MYPROC)); 

for (i=0;i<PROCS;i++)
  {entries = *(finish_ptr++) = *(Bin_Count + i);
   start_ptr = *(Bin + i);
   for (j=0;j<entries;j++)
     *(finish_ptr+j) = *(start_ptr+j);
   finish_ptr = (bin_ptr += max);
  }

E = finish_sptr;
E_ptr = ((int *) (E + MYPROC));
D_ptr = ((int *) (D + MYPROC));

barrier();

mpc_index(D_ptr,E_ptr,(max*sizeof(int)),ALLGRP);

barrier();

start_sptr = E;
finish_sptr =  D;

#else

max = buck_size - 2;

for (i=0;i<PROCS;i++)
  if ((*(Bin_Count)) > max)
    {fprintf(stderr,"ERROR: Too many elements to route in Step(7)\n");
     fflush(stderr);
     exit(1);
    }

max = buck_size;

D = start_sptr;
t = max*MYPROC;
entries = (*(Bin_Count + MYPROC)); 
C_ptr = (*(Bin + MYPROC));
D_ptr = ((int *) (D + MYPROC)) + t;

*(D_ptr ++) = entries; 

for (i=0;i<entries;i++)
  *(D_ptr + i) = *(C_ptr + i);

for (i=1;i<PROCS;i++)
  {dest = (MYPROC + i) & PROCS_m1;
   dest_ptr = ((int *global) (D + dest)) + t;
   *dest_ptr :- (*(Bin_Count + dest));
   bulk_store((dest_ptr + 1),(*(Bin + dest)),
				     ((*(Bin_Count + dest))*(sizeof(int))));
   barrier();
  }
 
all_store_sync();

start_sptr = D;
finish_sptr =  C;

barrier();

#endif


/***************************************************************************
			       Step (8)
***************************************************************************/


adjustments = 0;

start_ptr = ((int *) (start_sptr + MYPROC));
finish_ptr = ((int *) (finish_sptr + MYPROC));
for (j=0;j<(PROCS >> 1);j++)
  {tot1 = (*(start_ptr));
   s2_ptr = (s1_ptr = (start_ptr + 1)) + max;
   tot2 =  (*(s2_ptr - 1));
   (*(Number + j)) = tot = tot1 + tot2;
   (*(s1_ptr + tot1)) = (*(s2_ptr + tot2)) = MAX_VAL;
   v = tot1 - 1;
   while ((*(s1_ptr + v)) == MAX_VAL)
     {adjustments++;
      (*(s1_ptr + (v--))) = MAX_VAL_m1; 
     } 
   v = tot2 - 1;
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
   start_ptr += (2*max);
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

free(Bin);
free(Bin_Count);
free(Loc);
free(Number);
free(Count1);
free(Count2);

all_spread_free(Splitters);

}

