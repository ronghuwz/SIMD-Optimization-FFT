// Fast Fourier Transform (FFT) using 2D Cooley-Tukey algorithm GCC

// Original:
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <complex>
#include <ctime>
#include <math.h>
#include <immintrin.h>   // AVX intrinsics
#pragma GCC target("avx2")

const int N = 10000; // Size of the 2D square array

// Define a complex number
using Complex = std::complex<double>;

// Cooley-Tukey FFT algorithm
void simd_fft(std::vector<Complex>& data) {

    const int n = data.size();
    if (n <= 1)
    {
       return; 
    }
    

    std::vector<Complex> even(n/2);
    std::vector<Complex> odd(n/2);
    
    for (int i = 0; i < n / 2; ++i) {
        even[i] = data[i * 2];
        odd[i] = data[i * 2 + 1];
    }
    
    simd_fft(even);
    simd_fft(odd);

    
    std::vector<double> even_real_part (n / 2); // Allocate memory for real parts
    std::vector<double> even_imag_part (n / 2); // Allocate memory for imaginary parts
    std::vector<double> odd_real_part (n / 2); // Allocate memory for real parts
    std::vector<double> odd_imag_part (n / 2); // Allocate memory for imaginary parts

    // Decompose the data into real + imaginary
    for (int i = 0; i < n / 2; ++i) {
        even_real_part[i] = even[i].real();
        even_imag_part[i] = even[i].imag();
        odd_real_part[i] = odd[i].real();
        odd_imag_part[i] = odd[i].imag();
    }
    

    
    __m256d even_real, even_imag, odd_real, odd_imag, twiddle_real, twiddle_imag, odd_temp_real, odd_temp_imag,
      plus_temp_real, plus_temp_imag, sub_temp_real, sub_temp_imag;
    for (int k = 0;  k < n / 2; k += 4) {
        if (k + 4 >= n / 2)
        {
            //break;
             for (int ii = k;  ii < n / 2; ++ ii) {
                 Complex t = std::polar(1.0, -2.0 * M_PI * ii / n) * odd[ii];
                 data[ii] = even[ii] + t;
                 data[ii + n / 2] = even[ii] - t;  
             }
             continue;
        }
        // Load real and imaginary parts of even and odd elements
        // std::cout << "Here 2.1, k: " << k << std::endl;
        even_real = _mm256_loadu_pd(&even_real_part[k]); // Load real parts
        even_imag = _mm256_loadu_pd(&even_imag_part[k]); // Load real parts
        // std::cout << "Here 2.1.0, k: " << k << std::endl;
        odd_real = _mm256_loadu_pd(&odd_real_part[k]); // Load imaginary parts
        odd_imag = _mm256_loadu_pd(&odd_imag_part[k]); // Load imaginary parts
        // std::cout << "Here 2.2, k: " << k << std::endl;

        // Normally, you would precompute or calculate these
        twiddle_real = _mm256_set1_pd(cos(-2.0 * M_PI * k / n));
        twiddle_imag = _mm256_set1_pd(sin(-2.0 * M_PI * k / n));

        // Complex multiplication: (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
        odd_temp_real = _mm256_sub_pd(_mm256_mul_pd(odd_real, twiddle_real), _mm256_mul_pd(odd_imag, twiddle_imag));
        odd_temp_imag = _mm256_add_pd(_mm256_mul_pd(odd_real, twiddle_imag), _mm256_mul_pd(odd_imag, twiddle_real));
        // std::cout << "Here 2.4, k: " << k << std::endl;
        
        plus_temp_real = _mm256_add_pd(even_real, odd_temp_real);
        plus_temp_imag = _mm256_add_pd(even_imag, odd_temp_imag);
        sub_temp_real = _mm256_sub_pd(even_real, odd_temp_real);
        sub_temp_imag = _mm256_sub_pd(even_imag, odd_temp_imag);
        

        for (int m = 0; m < 4; m ++){
            data[k + m] = std::complex<double>(plus_temp_real[m], plus_temp_imag[m]);
            data[k + m + n / 2] = std::complex<double>(sub_temp_real[m], sub_temp_imag[m]);   
        }
    }
  
}

int main() {
    // Create a large 2D array (N x N)
    std::vector<std::vector<Complex>> data(N, std::vector<Complex>(N));

    // Initialize data with some values (e.g., a pattern)

    // Start the timer
    std::time_t start_time = std::time(nullptr);

    // Perform 2D FFT on the data
    for (int i = 0; i < N; i++) {
        // std::cout <<"i: " << i << ", data size: " << data[i].size() << std::endl;
        simd_fft(data[i]);
    }
    for (int j = 0; j < N; j++) {
        std::vector<Complex> column(N);
        for (int i = 0; i < N; i++) {
            column[i] = data[i][j];
        }
        simd_fft(column);
        for (int i = 0; i < N; i++) {
            data[i][j] = column[i];
        }
    }

    // End the timer
    std::time_t end_time = std::time(nullptr);

    std::cout << "2D FFT completed in " << difftime(end_time, start_time) << " seconds." << std::endl;

    return 0;
}