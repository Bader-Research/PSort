# PSort: Parallel Sorting using Sampling with Randomized and Deterministic Approaches

Previous schemes for sorting on general-purpose parallel machines have
had to choose between poor load balancing and irregular communication
or multiple rounds of all-to-all personalized communication. This code
provides a sample sort which uses only two rounds of regular
all-to-all personalized communication in a scheme that yields very
good load balancing with virtually no overhead. Moreover, unlike
previous variations, our algorithm efficiently handles the presence of
duplicate values without the overhead of tagging each element with a
unique identifier. This algorithm was implemented in Split-C and run
on a variety of platforms, including the Thinking Machines CM-5, the
IBM SP-2, and the Cray Research T3D. This implementation of our new
algorithm seems to outperform all similar algorithms known to the
authors on these platforms, and its performance is invariant over the
set of input distributions unlike previous efficient algorithms. Our
results also compare favorably with those reported for the simpler
ranking problem posed by the NAS Integer Sorting (IS) Benchmark.

The code also provides a new deterministic parallel sorting algorithm
based on the regular sampling approach. The performance compares
closely to that of our random sample sort algorithm, which seems to
outperform all similar algorithms known to the authors on these
platforms. Together, their performance is nearly invariant over the
set of input distributions, unlike previous efficient
algorithms. However, unlike our randomized sorting algorithm, the
performance and memory requirements of our regular sorting algorithm
can be deterministically guaranteed.

References:

D. R. Helman, D. A. Bader, and J. Já Já . "A Randomized Parallel Sorting Algorithm With an Experimental Study," Journal of Parallel and Distributed Computing, 52(1):1-23, 1998.

D.R. Helman, J. Já Já , D.A. Bader. "A New Deterministic Parallel Sorting Algorithm With an Experimental Evaluation," ACM Journal of Experimental Algorithmics, 3(4):1-24, 1998.
