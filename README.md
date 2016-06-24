# BLPhi
Sequence alignment is becoming increasingly important in our current day and age, and with the rise of coprocessors, it is important to adapt sequence alignment algorithms to the new architecture. Parallelization using SIMD technology has previously been achieved that implement alignment algorithms efficiently such as SWIPE, described by Rognes in 2011.

The Intel Xeon Phi provides a solid architecture which can be used and exploited to maximize the speed in sequence alignment. It is therefore important, to develop algorithms that are able to use efficiently the coprocessor to maximize throughput.


A different approach and implementation is described and benchmarked. The new program, called BLPhi, aims to reduce execution time by using a filtering method. BLPhi takes advantage of the architecture of the Intel Xeon Phi and provides a parallel solution to sequence alignment using SIMD technology.


While exact alignment methods such as the Smith-Waterman are too sensitive for database searching, BLPhi uses a filtering technique that can be adapted to any length given. This gives the user the ability to adapt the search for his needs and may, perhaps, make the search more restrictive.  


An analysis with current state of the art alignment methods is also presented. We perform speed and accuracy analysis to see the potential that BLPhi has among its competitors.
