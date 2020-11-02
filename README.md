# comp90025-pmc-notes
COMP90025 - Parallel and Multicore Computing - 2020S2 - Exam review/summary sheet

## terminology
- knowledge [location in xuliny's lecture slide]
- {xxx} link current knowledge with xxx (you can just search xxx in this document)
- def: definition

## 01 introduction
### definition
- Work: w(n) = t(n) × p(n) [problem with input size n]
  - t(n) = parallel steps
  - p(n) = #processors
- Size: 
  - total number of operations that the parallel algorithm undertakes
  - T(n) ≤ size(n) ≤ w(n) 
  - prove at [24]
- Brent's Principle: proof and def [26]
- Optimality: def [27]
  - A PRAM algorithm is optimal if ...
- Efficient: def [27]
- Speedup = S(p) = T(n) / t(n) **(complexity)**
    - def [28]: reduced the complexity of the problem by a factor of ...
    - T(n) = the **best sequential complexity** for a given problem
    - t(n) = the **parallel complexity** of a given parallel algorithm
- Speedup = ts / tp **(in wall-clock time)**
  - def [36]
- Amdahl's Law: S(p) = 1/f | large p
  - f = fraction of this time that cannot be parallelized
  - predict the maximum achievable speedup for a given program
- Gustafson's Law: S(p) = p + (1 − p) s.
  - s = sequential part
- efficiency = E = S(p) / p
  - max p for optimal processor allocation? [29]
- feasibility: def [32]
  - 如果我们algorithm demand的资源随n增长的太快的话, 氪金也解决不了问题   
  - feasible [33]
  - feasible highly parallel [33]
  - inherently sequential [33]
  - Nick’s Class: def [34]
### PRAM
- PRAM [20]
- 4 PRAM sub-categories [21]
  - C: concurrent, E: exclusive, R: read, W: write
  - EREW
  - CREW
  - CRCW
  - ERCW (not commonly used)
- 4 CW conflicts model -> CRCW is used [22]
  - ||description
    |---|---|
    |COMMON|processors都agree了才write
    |ARBITRARY|随便选一个processor去执行write
    |PRIORITY|有指定priority(processor with lowest id)的processor执行write
    |COMBINING|被write的是"a combination of what processors are going to write?"<br /> 例如: max(values written), min(), sum() <br />十分的Powerful, 比如sort in O(1) time, 但machine也不好build, 所以人们不常用
- Power of PRAM model variations
  - The models are listed in increasing order of "power":
    1. any algorithm that runs on a EREW PRAM will run on a CREW PRAM,
    2. any algorithm that runs on a CREW PRAM will run on a COMMON CRCW PRAM, and so on. 
    3. [Small PRAMs can simulate larger PRAMs](http://pages.cs.wisc.edu/~tvrdik/2/html/Section2.html#Simulation1) 
    - However it is not true in the other direction.
      - impossible example [44]
  - Any algorithm for a CRCW PRAM in the PRIORITY model 
    - can be "simulated" by 
      - an EREW PRAM with the same number of processors  
        and 
      - with the parallel time increased by only a factor of O(log p).
    - { EREW SimulatePriority }
  - Any algorithm for a PRIORITY PRAM 
    - can be simulated by
      - a COMMON PRAM with no loss in parallel time  
       and
      - provided sufficiently many processors are available.
    - { COMMON SimulatePriority }
  - CR PRAMs 
    - can be simulated by
      - EREW PRAM by
        - without using any additional processor/time complexity
    - { Optimal EREW Broadcast }
### EREW algorithm
- Suboptimal EREW Lambda
  - |||
    |---|---|
    |input size|n
    |t(n)|O(log n) steps
    |p(n)|p = n/2
    |T(n)|O(n)
    |Lambda|:: Ord a => (a, a) -> a
  - ```
    j = n/2

    while j ≥ 1
        for each processor i from 0 to j-1:
            input[i] = Lambda(input[2i], input[2i + 1])
        j /= 2
    
    return input[0]
    ```
- Optimal EREW Summation [42] 
    - |||
      |---|---|
      |input size|n
      |t(n)|O(n/p + log p) = O(log n) steps
      |p(n)|p = n/log n
      |T(n)|O(n)
    - optimal EREW pattern (1)
- Optimal EREW Maximum
  - |||
    |---|---|
    |input size|n^2
    |t(n)|O(n/p + log n) <br /> = O(log n / n + log n^2 - log log n) steps by substitute in p(n) <br /> = O(log n / n + 2 log n) steps <br /> = O(log n) steps
    |p(n)|p = n^2/log n
    |T(n)|O(n^2)
   - optimal EREW pattern (1)
   - ```
     # initialize
     B[p]
     for each processor i from 0 to p-1:
        B[i] = MIN
     
     for each processor i from 0 to p: # in O(log n) steps with n / log n processors
        # input size / p(n) = log n
        for j from i*(input size/p) to (i+1)*(input size/p) - 1:
            B[i] = max(B[i], input[j])
     
     # Suboptimal EREW Maximum
     #   in O(log p) steps with p/2 processors
     return Suboptimal EREW Lambda=Max (input = B, input size = p)
     ```
- Suboptimal EREW Replication [53]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(log n)
    |p(n)|p = n
    |T(n)|O(n)
- Optimal EREW Replication [54]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(n/p + log n)
    |p(n)|p = n/log n
    |T(n)|O(n)
  - input an arr [a1, a2, ..., an] return [a1, a1, ..., a1]
- Optimal EREW Broadcast [55]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(log p) = O(log n + log log n) = O(log n) steps
    |p(n)|p = n/log n
    |T(n)|O(n)
  - processor 0 has val, have a [val, val, ..., val].size=n in memory
    - by using Suboptimal EREW Replication
- EREW SimulatePriority [59]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(log n) steps
    |p(n)|p = n
    |T(n)|O(n)
  - How to execute:
    ```
    Input: W[n]: W[i] is the address that PRIORITY processor i wants to write to
    # initialize A[n] with 
    for i in range (0, p):
        A[i] = (i, w[i], False) # (index, address, whether we are allowed to write)
        
    sort(A, lambda = \x, y, z -> cmp y, then cmp x then, increase=True) # Sort by address W[i] then id[i]

    for i in range (0, p):
        if i == 0:  # Lowest id for the smallest address always wins
            (a, b, _) <- A[0]
            A[0] = (a, b, True) 
        else:
            (a, b, _) <- A[i]
            c <- A[i] 跟 A[i-1] 比较 snd; if 相等 True otherwise False
            A[i] <- (a, b, c)

    for i in range (0, p):
        if A[i].c:
            # processor i performs write
            PRIORITY processor i can write to W[i]
        else:
            Discard write operation
    ```
#### optimal EREW pattern
- (1)
  - ```
      t(n): O(n/p + log n)
      Input: array
      initialize localArr[p]

      for each processor_i: 
          sequentially cal subarray to localArr[i]
      
      # reduce
      parallelize reducing localArr to finalResult by Suboptimal EREW algorithm (LAMBDA)
      ```
### COMMON algorithm
- Optimal COMMON Logical_OR [44]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(n/p) steps
    |p(n)|n
    |T(n)|O(n^2)
  - similar to optimal EREW pattern (1)
    - except no O(log p) because we use localVar instead of localArr (no race condition beacuse we have COMMON), so no parallelize reduce just return
- Suboptimal COMMON Maximum [47]
  - |||
    |---|---|
    |input size|n^2
    |t(n)|O(log log n) steps
    |p(n)|n^2
    |T(n)|O(n)
- COMMON SimulatePriority [57]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(1) steps
    |p(n)|n^2
    |T(n)|O(n) on Priority
  - How to execute:
    ```
    Input: W[n]: W[i] is the address that PRIORITY processor i wants to write to
    M[n] = [True, ..., True]
    
    check whether W has duplicate, if any, M[higher p_i] = False
      # only M[lowest i] = True if CW Conflict
    
    for i in range(0, p):
        if M[p_i]:
            # processor i performs write
            PRIORITY processor p_i can write to W[p_i]
        else:
            Discard Write operation
    ``` 
     
### PRIORITY algorithm
- Optimal PRIORITY ElementUnique
  - |||
    |---|---|
    |input size|n
    |t(n)|O(1) steps
    |p(n)|n
    |T(n)|O(n)
  - optimal Priority pattern (1)
#### optimal Priority pattern
- (1)
  - ```
      t(n): O(1)
      Input: array
      initialize localArr[p], localVar

      for each processor_i:
        localArr[input[i]] <- i  # As PRIORITY, only lowest i will be written

      for each processor_i: 
        if localArr[input[i]] != i:  # implies other processor has duplicate value
          change localArr
      
      return localArr
      ```


## 02 architecture
- PE: processing element
- criticism of the PRAM [2] = why we need architecture
- Flynn's taxnomy [4]
  - SISD [5]
  - SIMD [6, 7]
    - good for data parallelism
    - GPUs contain many small SIMD modules
  - MIMD [8]
  - SIMD v.s. MIMD
    - SIMD
      - adv
        1. simple to code
        2. PRAM algorithms are close to SIMD
    - MIMD
      - adv
        1. more flexible | given a number of PE
      - disadv
        1. more expensive
        2. have multiple control unit here, which takes more space on chip and cost more time for instruction stream to arrive **(force the clock speed to be lower)**
  - MISD: There are no machines widely accepted to be MISD
- Schwartz's parallel machine classes [10]
  - paracomputers = shared memory multiprocessor
  - ultracomputer = distributed memory multiprocessor
  - Shared versus distributed memory [12]
    - shared memory system
      - = symmetric multiprocessor (SMP) | all processors are identical
      - adv: easier to code
    - distributed memory system
      - = message passing machine
      - adv: scale up with less cost
    - distributed shared memory (DSM)
      - = has distributed memory but a single address space
- Uniformity of shared memory access [13]
  - |a shared memory system can be either|||adv|disadv
    |---|---|---|---|---
    |UMA|uniform memory access|[14]|Caches improve performance| but need protocols to keep cache coherency
    |NUMA|non-uniform memory access|[15]|Access to local memory is fast| Access to remote memory is slower
    |COMA|cache-only memory architecture|[]|
- Simultaneous memory access [16]
  - Real memory access is different form theory
  - Memory Access Patterns:
    1. Memory is divided into one or more banks (roughly, sticks of RAM).
    2. Once a bank is accessed, it has a recharge time before it can be accessed again.
  - Addressing memory banks [17]   
    <img width="80%" src="./docs/1.jpg"/>
  - So, parallel memory access, rather than parallel processing
- Coprocessors [18]
  - GPU and Xeon Phi
  - |adv|disadv|
    |---|---|
    |low cost parallelism|architecturally they present more challenges
- Implicit versus Explicit [19]
  - parallelism implicit/explicit
  - decomposition implicit/explicit
  - mapping implicit/explicit
  - communication implicit/explicit
  - ||adv||
    |---|---|---|
    |implicit|easier to code
    |explicit|high performance|OpenMP
- SPMD and MPMD [20]
  |||def|Framework
  |---|---|---|---|
  |SPMD|Single Program Multiple Data|a single program that executes on a different portion of the data|OpenMP #pragma for<br />OpenMPI<br />OpenCL
  |MPMD|Multiple Program Multiple Data|1. different programs are written, for different threads<br /> 2. one thread may be a "master" thread, with different code from "worker" threads|OpenMP #pragma sections/single
- |parallelism model||suitable when||
  |---|---|---|---|
  |Thread parallelism model|[21]|running on<br /> 1. a single UMA machine or <br />2. a distributed shared memory architecture that implicitly distributes threads|OpenMP
  |Process parallelism model|[22]|running on multiple machines, i.e., when there is no shared memory (NUMA).
  |Hybrid parallelism model|[23]|running a single process on each machine and having each process use multiple threads

## 03 OpenMP
- What is OpenMP? [3]
  - Thread Parallelism model
  - Explicit Parallelism
  - Fork-Join Model [4]
- Compiler Directives [7]
  - ``` #pragma omp parallel default ( shared ) private ( beta ) ```
- Run-time Library Routines [8, 9]
  - ``` int omp_get_num_threads ( void ); ```
  - [大全](https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/optimization-and-programming-guide/openmp-support/openmp-library-support/openmp-run-time-library-routines.html)
- Environment Variables [16, 17]
  - ``` export OMP_NUM_THREADS=10 ```
  - ```OMP_DYNAMIC```
  - ```OMP_NESTED```
  - ```OMP_THREAD_LIMIT```
  - ```OMP_STACKSIZE```
  - [大全](https://docs.oracle.com/cd/E19205-01/819-5270/aewcb/index.html)
- C/C++ general code structure [10, 11, 12]
- OpenMP Directives [18, 19, 20, 21, 22, 23, 24]
  - ``` 
    # pragma omp parallel [ clause [ clause ...]] 
    structured_block 
    ```
  - Work-Sharing Constructs 
    - NOWAIT
      - If specified, then threads do not synchronize at the end of the parallel loop.
    - ``` #pramga omp parallel for ``` [23]
      - SCHEDULE (Dynamic vs Static [24, 28]) 
        - Describes how iterations of the loop are divided among the threads in the team.
        - |scheduler|adv|disadv
          |---|---|---|
          |Dynamic|better when the iterations may take very different amounts of time|has scheduler overhead <br />After each iteration, the threads must stop and receive a new value of the loop variable to use for its next iteration.
          |Guided|allows a tradeoff: <br /> 1. less overhead than dynamic <br />2. better splitting of workload among threads than static.
      - ORDERED clause [30]
        - Specifies that the iterations of the loop must be executed as they would be in a serial program.
      - COLLAPSE clause [31]
        - Specifies how many loops in a nested loop should be collapsed into one large iteration space and divided according to the schedule clause.
    - ``` 
      #pramga omp sections [ nowait ] 
          # pragma omp section
              structured_block
          # pragma omp section
              structured_block
          # pragma omp section
              structured_block
      ``` 
      [32, 33]
    - ```
      # pragma omp single [ nowait ]
          structured_block
      ```
      [34, 35]
- Synchronization and sharing variables
  1. critical [37]: ensures executions of the critical block to one thread at a time.
  2. barrier [38]: ensures that all threads in the team have reach the barrier before any
thread is allowed to continue beyond the barrier
  3. atomic [39]: ensures that updates to a variable are free from conflict.
     - atomic v.s. critical [39]
  4. Example - parallel pi program
- Memory Consistency and flush
  1. flush [40, 41, 42, 43]
- Data environment
  1. stack and heap [44]
     - In a multi-threaded situation each thread will have its own completely independent stack,
       - so the threads will copy the stack variables but they will share the heap variables.
  2. shared and private [45, 46, 47]
  3. firstprivate(x) [48]
     - initializes each private copy with the corresponding value from the master thread
  4. lastprivate(x) [49]
     - passes the value of a private from the last thread to a global variable.
  5. reduction [50, 51]
     - ``` # pragma omp parallel for reduction(+:x) default (shared) ```
       - add will initialize x with default 0
       - perform add in each thread
       - and add all threads' results and global variable to global variable
- Dijkstra algorithm [52, 53, 54, 55]
## 04a communication
- Time to send a message: t_m
  - t_m = t_f + t_b * l
    - t_f = a fixed time
    - t_b = a time per byte
    - l = number of bytes in the message
  - efficient message passing = t_b * l >> t_f
- Parallel run time: t_p
  - t_p = t_comp + t_comm
    - t_comp: computation time
    - t_comm: communication time
  - speedup: S
    - S = t_s / t_p
      - t_s: sequential time
- Delay hiding [6]
  - = when t_p not valid?
  - Communication and computation can happen at the same time
  - Part of the art of parallel programming is making sure that all processors/processes have enough data to keep them busy
- Granularity [7]
  - granularity = t_comp / t_comm
  - tradeoff between parallelism and communication
- Communication primitives [8]
  - what are the basic tasks happen when send data
- send() and receive() [9, 10]
  - send/receive(destination identifier, message type, data type and contents)
  - blocking or or local-blocking non-blocking [11, 12, 13, 14]
    - non-blocking(): send will keep sending in the background
    - non-blocking(): receive will typically return whatever data has been received so far, and indicate to the caller how much has been read [15]
  - Blocking and non-blocking connote (意味着) synchronous and asynchronous communication. [16]
- Communication patterns [17-26]
  - Why we need communication pattern?
    - Distributed memory architectures require message passing and this in turn leads to communication patterns.
    - Shared memory parallel programs don’t explicitly use communication patterns
    - Distributed shared memory architecture needs
## 04b OpenMPI

## 05 prefix sum
- sequential prefix sum (or other dyadic operation)
  - |||
    |---|---|
    |input size|n
    |T(n)|O(n)
  - ```
    int X[n];
    for(i=1; i<n; i++) {
        X[i] = X[i-1] + X[i];
    }
    ```
- Suboptimal CREW Upper/Lower parallel prefix algorithm [6]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(log n)
    |p(n)|p = n
    |T(n)|O(n)
  - can use { Optimal EREW Broadcast } to make CREW to EREW for operation on line 12
  - Execution on X = [5, 2, 3, 8, 2, 3, 4, 5], n = 8
    - -> means the return of the function call
    - line 7, UpperLower([5, 2, 3, 8], 4) -> S[0:4] = [5, 7, 10, 18]
      - line 7, UpperLower([5, 2], 2) -> S[0:2] = [5, 7]
        - line 7, UpperLower([5], 1) -> S[0:1] = 5
        - line 9, UpperLower([2], 1) -> S[1:2] = 2
        - // now S = [5, 7]
        - line 10, i from 1 to 1:
          - S[1] += S[0] -> S[1] = 5+2 = 7
        - return S = [5, 7]
      - line 9, UpperLower([3, 8], 2) -> S[2:4] = [3, 11]
        - line 7, UpperLower([3], 1) -> S[0:1] = 3
        - line 9, UpperLower([8], 1) -> S[1:2] = 8
        - // now S = [3, 8]
        - line 10, i from 1 to 1:
          - S[1] += S[0] -> S[1] = 3+8 = 11
        - return S = [3, 11]
      - // now S = [5, 7, 3, 11]
      - line 10, i from 2 to 3:
        - S[2] += S[1] -> S[2] = 7+3 = 10
        - S[3] += S[1] -> S[3] = 11+7 = 18
      - return S = [5, 7, 10, 18]
    - line 9, UpperLower([2, 3, 4, 5], 4) -> S[4:8] = [2, 5, 9, 14]
      - line 7, UpperLower([2, 3], 2) -> S[0:2] = [2, 5]
        - line 7, UpperLower([2], 1) -> S[0:1] = 2
        - line 9, UpperLower([3], 1) -> S[1:2] = 3
        - // now S = [2, 3]
        - line 10, i from 1 to 1:
          - S[1] += S[0] -> S[1] = 2+3 = 5
        - return S = [2, 5]
      - line 9, UpperLower([4, 5], 2) -> S[2:4] = [4, 9]
        - line 7, UpperLower([4], 1) -> S[0:1] = 4
        - line 9, UpperLower([5], 1) -> S[1:2] = 5
        - // now S = [4, 5]
        - line 10, i from 1 to 1:
          - S[1] += S[0] -> S[1] = 4+5 = 9
        - return S = [4, 9]
      - // now S = [2, 5, 4, 9]
      - line 10, i from 2 to 3:
        - S[2] += S[1] -> S[2] = 5+4 = 9
        - S[3] += S[1] -> S[3] = 5+9 = 14
      - return S = [2, 5, 9, 14]
    - // now S = [5, 7, 10, 18, 2, 5, 4, 9]
      - line 10, i from 4 to 7:
        - S[4] += S[3] -> S[4] = 2+18 = 20
        - S[5] += S[3] -> S[5] = 5+18 = 23
        - S[6] += S[3] -> S[6] = 9+18 = 27
        - S[7] += S[3] -> S[7] = 14+18 = 32
    - return S = [5, 7, 10, 18, 20, 23, 27, 32]
- Suboptimal EREW Odd/Even parallel prefix algorithm [12, 13]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(log n)
    |p(n)|p = n
    |T(n)|O(n)
  - Execution on X = [5, 2, 3, 8, 2, 3, 4, 5], n = 8
    - -> means the return of the function call
    - line 2, False
    - line 8, index [0, 2, 4, 6] gives
      - S = [5, _, 3, _, 2, _, 4, _]
    - line 11, index [1, 3, 5, 7] gives
      - P's index floor(i/2) = [0, 1, 2, 3]
      - P = [7, 11, 5, 9]
    - line 13, OddEven(X=P, 4) -> A = [7, 9, 12, 20, 22]
      - line 2, True
        - S = [7, ...]
        - for i from 1 to 4
          - S[1] = X[1] + S[0] = 7+7 = 14
          - S[2] = X[2] + S[1] = 11+14 = 28
          - S[3] = X[3] + S[2] = 5+28 = 33
          - S[4] = X[4] + S[3] = 9+20 = 22
        - return S = [7, 9, 12, 20, 22]
    - line 14, for i from 2 to 7
      - line 15, if i mod 2 == 0:
        - S[2] += A[0] = 3+7 = 10
        - S[4] += A[1] = 2+9 = 11
        - S[6] += A[2] = 4+12 = 16
    - // now S = [5, _, 10, _, 11, _, 16, _]
    - return S = [5, _, 10, _, 11, _, 16, _]
- Ladner and Fischer's parallel prefix algorithm [18]
- Optimal EREW Prefix sum
  - |||
    |---|---|
    |input size|n
    |t(n)|O(log n)
    |p(n)|p = n / log n
    |T(n)|O(n)
  - ```
    input[n]
    # do sequential prefix sum on sub-array size=log(n) in parallel
    # O(n/p) = O(log n) steps
    for processor_i in range(0, p) do in parallel:
        for j in range(ceil( i * (n/p) ) + 1, ceil( (i+1) * (n/p) ) :  # +1 as first element in subarray doesn't need to be prefixed
            input[j] += input[j-1]

    # O(p) steps
    for k in range(1, p):
        last <- input[ceil(k * (n/p)) - 1]  # take sub-array (k-1)'s last prefix sum

        # O(1) step
        for processor_i in range(0, p) do in parallel:
            input[ceil(k*p) + i] += last  # update previous subarray's prefix to cur subarray
    ```
## 06

## 07

## 08

## 09

## 10
