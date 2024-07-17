#include <omp.h>
#include <iostream>

// [[Rcpp::export]]
void test_openmp() {
    #pragma omp parallel
    {
        std::cout << "Hello from thread " << omp_get_thread_num() << std::endl;
    }
}
