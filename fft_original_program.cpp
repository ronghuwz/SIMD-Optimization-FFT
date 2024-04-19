// Fast Fourier Transform (FFT) using 2D Cooley-Tukey algorithm

// Original:

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <ctime>

const int N = 10000; // Size of the 2D square array

// Define a complex number
using Complex = std::complex<double>;

// Cooley-Tukey FFT algorithm
void fft(std::vector<Complex>& data) {
    const int n = data.size();
    if (n <= 1) return;

    std::vector<Complex> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; i ++) {
        even[i] = data[i * 2];
        odd[i] = data[i * 2 + 1];
    }

    fft(even);
    fft(odd);

    for (int k = 0; k < n / 2; k++) {
        Complex t = std::polar(1.0, -2.0 * M_PI * k / n) * odd[k];
        data[k] = even[k] + t;
        data[k + n / 2] = even[k] - t;
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
        fft(data[i]);
    }
    for (int j = 0; j < N; j++) {
        std::vector<Complex> column(N);
        for (int i = 0; i < N; i++) {
            column[i] = data[i][j];
        }
        fft(column);
        for (int i = 0; i < N; i++) {
            data[i][j] = column[i];
        }
    }

    // End the timer
    std::time_t end_time = std::time(nullptr);

    std::cout << "2D FFT completed in " << difftime(end_time, start_time) << " seconds." << std::endl;

    return 0;
}