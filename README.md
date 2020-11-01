# comp90025-pmc-notes
COMP90025 - Parallel and Multicore Computing - 2020S2 - Exam review/summary sheet

## terminology
- knowledge [location in xuliny's lecture slide]
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
### algorithm
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
- Optimal COMMON Logical_OR [44]
  - |||
    |---|---|
    |input size|n
    |t(n)|O(n/p) steps
    |p(n)|n
    |T(n)|O(n^2)
  - optimal EREW pattern (1)
- Suboptimal COMMON Maximum [47]
  - |||
    |---|---|
    |input size|n^2
    |t(n)|O(log log n) steps
    |p(n)|n^2
    |T(n)|O(n)
- Optimal EREW Maximum
  - |||
    |---|---|
    |input size|n^2
    |t(n)|O(n/p + log n) <br /> = O(log n / n + log n^2 - log log n) steps by substitute in p(n) <br /> = O(log n / n + 2 log n) steps <br /> = O(log n) steps
    |p(n)|p = n^2/log n
    |T(n)|O(n^2)
   - optimal EREW pattern
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
### optimal EREW pattern
- (1)
  - ```
      Input: array
      initialize localArr[p]

      for each processor_i: 
          sequentially cal subarray to localArr[i]
      
      parallelize reducing localArr to finalResult by Suboptimal EREW algorithm (LAMBDA)
      ```

## 02 

## 03

## 04

## 05

## 06

## 07

## 08

## 09

## 10
