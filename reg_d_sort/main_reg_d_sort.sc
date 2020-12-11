
/*                                                                    
 *
 * main_reg_d_sort.sc - Main Function for Calling Function to Sort Double 
 *                      Precision Floating Point Numbers by a New 
 *                      Deterministic Algorithm Based on the Regular
 *                      Sampling Approach - Includes Timing and Error 
 *                      Checking Routines (Split-C Code)
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
 * Filename:            main_reg_d_sort.sc
 * History:
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <split-c/split-c.h>
#include <split-c/control.h>
#include "reg_d_sort.h"
#include "d_benchmarks.h"

#define PROCS_m1 (PROCS - 1)

/*************************************************************************** 
*
*                         CHECK FUNCTION: 
*
*
*
* The function check verifies that the data has been sorted correctly.  If 
* an error is detected, the value of tot_time (the time required for 
* sorting) is set to -1. 
*
* Called with:
*       Out           -> a spread pointer to the spread array containing
*                        the sorted data.
*       col_size      -> the number of elements per processor before sorting. 
*       values        -> the number of elements at a  particular processor
*                        after sorting.
*       *tot_time_ptr -> a pointer to the location of the variable tot_time.
*       s             -> the number of samples collected from each processor.
*        
***************************************************************************/   



void Check(double *spread Out, int col_size, int values, double *tot_time_ptr, int s)

{

int sum, result1, result2, result3, result4, i;

double temp;

double *out_ptr;


/* First, check that the output values are correctly sorted at each
processor: */

sum = 0;
result1  = 0;
out_ptr = ((double *) (Out + MYPROC));
for (i=1;i<values;i++)
  if (out_ptr[i] >= out_ptr[i-1])
    sum++;  
if (sum == (values - 1))    
  result1++;

barrier();

result1 = all_reduce_to_one_add(result1); 


/* Second, check that the largest value at each processor is less than or
equal to the first value at the next processor */

result2 = 0;

if (MYPROC == PROCS_m1)
  result2++;
else 
  {temp = *(Out + MYPROC + 1);
   if (out_ptr[values-1] <= temp)
	result2++;
      
  }        
  
barrier();

result2 = all_reduce_to_one_add(result2); 


/* Third, check that the total number of elements at all the processors at the completion
is identical to the total number of elements at all the processors at the inception. */


result3 = all_reduce_to_one_add(values);


/* Lastly, check that the load balance is witin the gauranteed range */

result4 = 0;

if (values <= (col_size + PROCS*(col_size/s - 1)))
  result4++;

result4 = all_reduce_to_one_add(result4);

barrier();

if ((result1 != PROCS) || (result2 != PROCS) || 
    (result3 != (col_size*PROCS)) || (result4 != PROCS))
  *tot_time_ptr = -1;

}


/*************************************************************************** 
*
*			 SPLITC_MAIN FUNCTION: 
*
*
*
*      The function splitc_main provides a convenient framework from which 
* to call the sorting routine all_regular_d_sort.  Using the command 
* line parameters provided by the user, it calls the function 
* all_regular_d_sort four seperate times and then verifies the correctness 
* of each sort.  If the input is sorted correctly all four times, then the 
* fastest of the four sorting times is reported; otherwise, the value -1 
* is returned indicating that something has not been done correctly.
*  The reason why we automatically include four trials is to minimize 
* the effect on timing of such disruptions as lazy loading and context 
* switching.
*
***************************************************************************/



splitc_main(int argc, char *argv[])

{ 

int col_size,benchmark,t,values,group,range,s,w;

int *val_ptr;

double min_time, tot_time, *tot_time_ptr;

double *In;

double *spread Buffer1; 
double *spread Buffer2;
double *spread Out;
double *spread *out_array_ptr;


/*************************************************************************** 
*
*                 CALLING THE SPLITC_ MAIN FUNCTION: 
*
*
*
*      In order to call all_regular_d_sort, there are restrictions 
* on the possible values of its parameters which are enforced by
* splitc_main.  These restrictions do not reflect any intrinsic 
* limitations of the algorithm.  Rather, they made the initial construction 
* of the code simpler, and we have not as yet had the opportunity to 
* modify the code to run correctly for values outside these limits.
*
*     Note that the number of processors `PROCS' must be a power of two.  
*
* Splitc_main can be called with the following command line parameters: 
*
* (1) col_size -> the number of elements per processor in the input array,
*                 where there are three restrictions on its value: 
*
*                  <1> -> The value of `col_size' must be the same for
*                         every processor.
*                  <2> -> The value of `col_size' must be at least 
*                         max(PROCS*PROCS,1024).
*                  <3> -> PROCS must divide `col_size' evenly.                    
*
* (2) benchmark -> the benchmark used to initiate the input.  For a 
*                  description of these benchmarks, see our accompanying
*                  papers.  To choose a particular benchmark, enter 
*                  the corresponding number below.  The choices are:
*              
*                    <1> -> Uniform [U].
*                    <2> -> Gaussian [G].
*                    <3> -> g-Group [g-G].
*                    <4> -> Bucket Sorted [B].
*                    <5> -> Staggered [S].
*                    <6> -> Zero [Z].
*                    <7> -> Worst-Load Regular [WR].
*                    <8> -> Deterministic Duplicates [DD].
*                    <9> -> Randomized Duplicates [RD].
*
*
* (3) s -> the number of samples to be collected from each processor, 
*          where there are two restrictions on its value:
*
*           <1> The value of s must be at least PROCS.
*           <2> s must evenly divide `col_size/PROCS'.
*
* (4) group -> the size of the processor groups used for the g-Group [g-G] 
*              Benchmark, where `group' must divide PROCS evenly.
*
* (5) range -> the range of random values used for the Randomized 
*              Duplicates [RD] Benchmark.
*
   
***************************************************************************/

 
if (PROCSLOG < 0)
  {fprintf(stderr,"ERROR: PROCS must be a power of two!\n");
   fflush(stderr);
   exit(1);
  }

if (argc > 1)
  col_size = atoi(argv[1]);
else
  col_size = (PROCS*PROCS) > 1024 ? (PROCS*PROCS) : 1024;

if (all_reduce_to_all_min(col_size) != all_reduce_to_all_max(col_size))
  {fprintf(stderr,"ERROR: col_size must be the same for every processor!\n");
   fflush(stderr);
   exit(1);
  }

if ((col_size < (PROCS*PROCS)) || (col_size < 1024))
  {fprintf(stderr,"ERROR: col_size must be at least max(PROCS*PROCS,1024)!\n");
   fflush(stderr);
   exit(1);
  }

if ((col_size % PROCS) > 0)
  {fprintf(stderr,"ERROR: col_size must be evenly divided by PROCS!\n");
   fflush(stderr);
   exit(1);
  }

if (argc > 2)
  benchmark = atoi(argv[2]);
else
  benchmark = 1; 

if ((benchmark > 9) || (benchmark < 1))
  {fprintf(stderr,"ERROR: benchmark choice must be an integer between 1 and 9!\n");
   fflush(stderr);
   exit(1);
  }

if (argc > 3)
  s = atoi(argv[3]);
else
  s = PROCS;

if (s < PROCS)
  {fprintf(stderr,"ERROR: s must be at least PROCS!\n");
   fflush(stderr);
   exit(1);
  }

if (((col_size/PROCS) < s) || (((col_size/PROCS) % s) > 0))
  {fprintf(stderr,"ERROR: col_size/PROCS must be evenly divided by s!\n");
   fflush(stderr);
   exit(1);
  }

if (argc > 4)
  group = atoi(argv[4]);
else
  group = 2;

if (s < PROCS)
  {fprintf(stderr,"ERROR: s must be at least PROCS!\n");
   fflush(stderr);
   exit(1);
  }

if ((PROCS % group) > 0)
  {fprintf(stderr,"ERROR: group must evenly divide PROCS!\n");
   fflush(stderr);
   exit(1);
  }

if (argc > 5)
  range = atoi(argv[5]);
else
  range = 32;


min_time = 0;

srandom(21+1001*MYPROC); 

for (t=0;t<4;t++)
  {if (((col_size/PROCS + col_size/s + 5*PROCS)*PROCS) > (col_size + col_size/2 + 2))
     {Buffer1 = all_spread_malloc(PROCS,
		    ((col_size/PROCS + col_size/s + 5*PROCS)*PROCS*sizeof(double)));
      assert_spread_malloc(Buffer1);

      Buffer2 = all_spread_malloc(PROCS,
		   ((col_size/PROCS + col_size/s + 5*PROCS )*PROCS*sizeof(double)));
      assert_spread_malloc(Buffer2);
     }
   else
     {Buffer1 = all_spread_malloc(PROCS,
			     ((col_size + col_size/2 + 2)*sizeof(double)));
      assert_spread_malloc(Buffer1);

      Buffer2 = all_spread_malloc(PROCS,
			     ((col_size + col_size/2 + 2)*sizeof(double)));
      assert_spread_malloc(Buffer2);
     }

   out_array_ptr = &Out;
   val_ptr = &values;
      
   In =  ((double *) (Buffer1+MYPROC));  
   switch(benchmark)
     {
      case 1:
	U_d_Benchmark(col_size, In);
	break;

      case 2:
        G_d_Benchmark(col_size, In);
	break;

      case 3:
        gG_d_Benchmark(col_size, In, group);	
	break;
   
      case 4:
        B_d_Benchmark(col_size, In);
	break;

      case 5:
	S_d_Benchmark(col_size, In);
	break;

      case 6:
	Z_d_Benchmark(col_size, In);
	break;

      case 7:
	WR_d_Benchmark(col_size, In, s);
	break;

      case 8:
	DD_d_Benchmark(col_size, In);
	break;
   
      case 9:
	RD_d_Benchmark(col_size, In, range);
	break;

     }

      barrier();
      tot_time = get_seconds();
      
      all_regular_d_sort(col_size,s,Buffer1,Buffer2,out_array_ptr,
		      val_ptr);
      
      barrier();
      tot_time = get_seconds() - tot_time;

      tot_time_ptr = &tot_time;
      Check(Out,col_size,values,tot_time_ptr,s);

      if (t)
	{if (tot_time < min_time)
	   min_time = tot_time;
	}
      else 
	min_time = tot_time;
      all_spread_free(Buffer1);
      all_spread_free(Buffer2);

  }
  
on_one printf("For col_size = %d, benchmark = %d, number of samples = %d,\n
group = %d, and range = %d\n
the time required for sorting is %f seconds.\n",
	       col_size,benchmark,s,group,range,min_time);
on_one printf("---------------------------------------------------------------------------\n");
     


}


