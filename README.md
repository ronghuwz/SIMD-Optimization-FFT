# 1. Explanation of the code (functionality of the program, program flow, and how you have used SIMD to accelerate)

## a. Functionality:
The program is designed to perform a 2D FFT on a square matrix of complex numbers (represented as std::complex<double>). FFT is a widely used algorithm in digital signal processing for converting signals between time and frequency domains.

## b. Program Flow:
### Define the FFT function: 
The program first defined the FFT function.
### Initialization: 
A 2D square array (data) of size N x N is created and initialized with complex numbers.
### 1D FFT on Rows: 
The program first applies the FFT row-wise. It iterates through each row of the matrix and performs FFT on it.
### 1D FFT on Columns: 
Then, it extracts each column from the 2D array, performs FFT on these columns, and places the transformed column back into the original matrix. This two-step process effectively accomplishes the 2D FFT.

## c. SIMD Acceleration Explained:
### AVX2 Intrinsics: 
The program uses AVX2 intrinsics to accelerate parts of the FFT computation by: 
#### #include <immintrin.h> 
#### #pragma GCC target("avx2") 
AVX2 instructions allow performing operations on multiple data points simultaneously, leveraging the CPU's vector processing capabilities.
### Manual Vectorization: 
The manual vectorization involves explicitly loading data into AVX2 registers, performing vectorized operations, and then storing the results back. This method bypasses the need for the compiler to auto-vectorize the code, offering more direct control over the optimizations.
### Complex Number Handling: 
AVX2 does not natively support complex numbers, so the program separately handles real and imaginary parts as double-precision floating-point arrays. \
This separation allows the use of AVX2's floating-point arithmetic instructions.
For each group of four elements, it: \
Loads real and imaginary parts of "even" and "odd" elements separately into AVX2 registers. \
Calculates the "twiddle factors" (complex exponential terms) using AVX2 instructions for cosine and sine operations. \
Performs complex multiplication of the "odd" elements with the twiddle factors, followed by adding or subtracting these from the "even" elements, effectively combining them into the FFT result. 
### Performance Gains: 
By processing four double-precision elements in parallel, the SIMD-optimized sections of the code can achieve significant speedups, especially on large datasets (in this case N= 10000) where such parallelizable operations dominate the computation time.

# 2. Estimated speed up with explanation
AVX2 instructions operate on 256-bit wide registers. For double-precision floating-point operations, each AVX2 register can hold four double-precision values. Therefore, in an ideal scenario where the computation is purely CPU-bound, using AVX2 could theoretically provide up to a 4x speedup for the operations that have been vectorized, compared to scalar operations. However, due to factors like memory bandwidth limitations and overhead, the effective speedup for the entire FFT computation (including non-vectorized parts) might be lower, potentially in the range of 1.5x to 3x, depending on the specific hardware and runtime conditions.

# 3. Compilation steps and flags
Compile the C++ file named fft_SIMD_program_final.cpp into an executable named fft_SIMD_final.exe using gcc, then run the executable. 
### gcc Flag used: 
#### g++ -O2 -mavx2 fft_SIMD_program_final.cpp -o fft_SIMD_final
