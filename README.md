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
    - However it is not true in the other direction.
      - impossible example [44]
  - Any algorithm for a CRCW PRAM in the PRIORITY model 
    - can be "simulated" by 
      - an EREW PRAM with the same number of processors  
        and 
      - with the parallel time increased by only a factor of O(log p).
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
  - ||adv|
    |---|---
    |implicit|easier to code
    |explicit|high performance
- SPMD and MPMD [20]
  |||def|Framework
  |---|---|---|---|
  |SPMD|Single Program Multiple Data|a single program that executes on a different portion of the data|OpenMP<br />OpenMPI<br />OpenCL
  |MPMD|Multiple Program Multiple Data|1. different programs are written, for different threads<br /> 2. one thread may be a "master" thread, with different code from "worker" threads
- |parallelism model|suitable when|
  |---|---|
  |Thread model [21]|running on<br /> 1. a single UMA machine or <br />2. a distributed shared memory architecture that implicitly distributes threads
  |Process model [22]|running on multiple machines, i.e., when there is no shared memory (NUMA).
  |Hybrid model [23]|running a single process on each machine and having each process use multiple threads

## 03

## 04

## 05

## 06

## 07

## 08

## 09

## 10
