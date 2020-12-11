/*                                                                    
 *
 * i_benchmarks.h - Include File for Sorting Benchmarks for Integers
 *                  (Split-C Code)
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
 * Filename:            i_benchmarks.h
 * History:
 */



#ifndef I_BENCHMARKS_H
#define I_BENCHMARKS_H


void U_i_Benchmark(int col_size, int *In);
void G_i_Benchmark(int col_size, int *In);
void gG_i_Benchmark(int col_size, int *In, int group);
void B_i_Benchmark(int col_size, int *In);
void S_i_Benchmark(int col_size, int *In);
void Z_i_Benchmark(int col_size, int  *In);
void WR_i_Benchmark(int col_size, int *In, int s);
void DD_i_Benchmark(int col_size, int *In);
void RD_i_Benchmark(int col_size, int *In, int range);

#endif
