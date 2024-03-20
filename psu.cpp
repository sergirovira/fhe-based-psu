#include <random>
#include "openfhe.h"
#include <sstream>
#include <iomanip>
#include <string>
#include <cassert>
#include <openssl/sha.h>
#include <chrono>
#include "scheme/ckksrns/ckksrns-fhe.h"
#include "typedef.h"

// header files needed for serialization
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <unordered_map>

#include "cryptocontext.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
#include "scheme/ckksrns/ckksrns-parametergeneration.h"

using namespace lbcrypto;

int bin_size;
int num_multiplications_getlist;
int num_multiplications_getrot;
int num_multiplications_getnrot;
int num_multiplications_server;

int spliting = 16;

class Benchmark {
public:
  Benchmark() { start_point = std::chrono::high_resolution_clock::now(); }
  ~Benchmark() { Stop(); }
  void Stop() {
    std::chrono::time_point<std::chrono::high_resolution_clock> end_point =
        std::chrono::high_resolution_clock::now();

    auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(start_point)
                     .time_since_epoch()
                     .count();
    auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(end_point)
                   .time_since_epoch()
                   .count();
    std::cout << " " << (end - start) << " milliseconds" << std::endl;
  }

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start_point;
};

// Function to calculate the binomial coefficient (n choose k)
long long binomialCoefficient(int n, int k) {
    if (k == 0 || k == n) {
        return 1;
    }
    long long result = 1;
    for (int i = 1; i <= k; ++i) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

// Function to calculate the coefficients of the polynomial
void calculateCoefficients(int n) {
    for (int i = 0; i <= n; ++i) {
        double coefficient = 1.0 / std::pow(4, i) * binomialCoefficient(2 * i, i);
        std::cout << "Coefficient for x^" << i << ": " << coefficient << std::endl;
    }
}

void removeZeros(std::vector<double>& vec) {
    vec.erase(std::remove_if(vec.begin(), vec.end(), [](double val) {
        return val == 0.0;
    }), vec.end());
}

void removeZeros(std::vector<int>& vec) {
    vec.erase(std::remove_if(vec.begin(), vec.end(), [](int val) {
        return val == 0.0;
    }), vec.end());
}

bool areVectorsEqual(const std::vector<int>& vec1, const std::vector<double>& vec2) {
    // Check if the sizes are the same
    if (vec1.size() != vec2.size()) {
        std::cout << "Psu NOT correct (duplicates)" << std::endl;
        return 0;
    }

    // Check each element for equality
    for (size_t i = 0; i < vec1.size(); ++i) {
        if (log2(fabs(vec1[i] - vec2[i])) >= -3) {
            std::cout << "Psu NOT correct (wrong removed)" << vec1[i] << ", " << vec2[i] << ", " << log2(fabs(vec1[i] - vec2[i])) << std::endl;
            return 0;
        }
    }
    // Vectors are equal
    std::cout << "Psu correct" << std::endl;
    return 1;
}

bool areVectorsEqual(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    // Check if the sizes are the same
    if (vec1.size() != vec2.size()) {
        std::cout << "Psu NOT correct (duplicates)" << std::endl;
        return 0;
    }

    // Check each element for equality
    for (size_t i = 0; i < vec1.size(); ++i) {
        if (fabs(vec1[i] - vec2[i]) >= 1e-5) {
            std::cout << "Psu NOT correct (wrong removed)" << vec1[i] << ", " << vec2[i] << ", " << log2(fabs(vec1[i] - vec2[i])) << std::endl;
            return 0;
        }
    }
    // Vectors are equal
    std::cout << "Psu correct" << std::endl;
    return 1;
}

template <typename T> void print_vec1d(const vec1d_t<T> &v) {
  std::cout << "[";
  for (const auto &i : v) {
    //std::cout << std::to_string(i) << ", ";
    std::cout << i << ", ";
  }
  std::cout << "]\n";
}

template <typename T> void print_vec2d(const vec2d_t<T> &v) {
  std::cout << "[\n";
  for (const auto &i : v) {
    std::cout << "\t";
    print_vec1d(i);
    std::cout << ",";
  }
  std::cout << "]\n";
}


template <typename T> void print_vec3d(const vec3d_t<T> &v) {
  std::cout << "[\n";
  for (const auto &i : v) {
    std::cout << "\t";
    print_vec2d(i);
    std::cout << ",";
  }
  std::cout << "]\n";
}

double CalculateApproximationError(const std::vector<std::complex<double>>& result,
                                   const std::vector<std::complex<double>>& expectedResult) {
    if (result.size() != expectedResult.size())
        OPENFHE_THROW(config_error, "Cannot compare vectors with different numbers of elements");

    // using the Euclidean norm
    double avrg = 0;
    for (size_t i = 0; i < result.size(); ++i) {
        avrg += std::pow(std::abs(result[i].real() - expectedResult[i].real()), 2);
    }

    avrg = std::sqrt(avrg) / result.size();  // get the average
    return std::abs(std::log2(avrg));
}

void printBeautifulMatrix(const std::vector<std::vector<complex_t>>& matrix) {
    // Set the width and precision for formatting
    const int width = 8;
    const int precision = 2;

    // Iterate over the matrix and print each element with proper formatting
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            std::cout << std::setw(width) << std::fixed << std::setprecision(precision) << real(element);
        }
        std::cout << std::endl;
    }
}

std::vector<double> extract_vector(plaintext_t ptxt) {
  size_t slots = ptxt->GetCKKSPackedValue().size();

  std::vector<double> ptxt_vec;

  for (size_t i = 0; i < slots; i++) {
    ptxt_vec.push_back(real(ptxt->GetCKKSPackedValue()[i]));
  }

  return ptxt_vec;

}

std::vector<double> multiplyMatrixVector(const std::vector<std::vector<complex_t>>& matrix, const std::vector<double>& vector) {
    std::vector<double> result;

    if (matrix[0].size() != vector.size()) {
        std::cerr << "Error: Matrix and vector dimensions are not compatible for multiplication." << std::endl;
        return result;  
    }

    for (size_t i = 0; i < matrix.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < vector.size(); ++j) {
            sum += real(matrix[i][j]) * vector[j];
        }
        result.push_back(sum);
    }

    return result;
}

plaintext_t multiplyMatrixPlaintext(const std::vector<std::vector<complex_t>>& matrix, plaintext_t ptxt, cryptocontext_t cc, int depth, int N) {

    auto ptxt_vec = extract_vector(ptxt);

    auto res = multiplyMatrixVector(matrix,ptxt_vec);

    auto res_ptxt = cc->MakeCKKSPackedPlaintext(res, depth, 0, nullptr, N/2);

    return res_ptxt;
}

std::vector<double> multiplyVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    std::vector<double> result;
/*
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vectors must be of the same size." << std::endl;
        return result;  
    }
*/
    for (size_t i = 0; i < std::min(vector1.size(), vector2.size()); i++) {
        result.push_back(vector1[i] * vector2[i]);
    }

    return result;
}

std::vector<double> subVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    std::vector<double> result;

    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vectors must be of the same size." << std::endl;
        return result;  
    }

    for (size_t i = 0; i < vector1.size(); i++) {
        result.push_back(vector1[i] - vector2[i]);
    }

    return result;
}

std::vector<double> addVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    std::vector<double> result;

    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vectors must be of the same size." << std::endl;
        return result;  
    }

    for (size_t i = 0; i < vector1.size(); i++) {
        result.push_back(vector1[i] + vector2[i]);
    }

    return result;
}

std::vector<double> decodeVector(const std::vector<double>& vector1) {
    std::vector<double> result;

    for (size_t i = 0; i < vector1.size(); i++) {
        if(fabs(vector1[i]) > 1e-8) {
            result.push_back(1/vector1[i]);
        } else {
            result.push_back(0);
        }
    }

    return result;
}

std::vector<double> compareVectors(const std::vector<double>& vector1, const std::vector<double>& vector2, size_t bound) {
    std::vector<double> result;

    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vectors must be of the same size." << std::endl;
        return result;  
    }

    for (size_t i = 0; i < bound; i++) {
        if (std::fabs(vector1[i] - vector2[i]) < 1e-5) {
            result.push_back(1);
        } else {
            result.push_back(0);
        }
    }

    for (size_t i = bound; i < vector1.size(); i++) {
        result.push_back(0);
    }

    return result;
}

std::vector<double> signVector(const std::vector<double>& v) {
    std::vector<double> result(v.size());

    for (size_t i = 0; i < v.size(); i++) {
        if (fabs(v[i]) < 1e-8) {
            result[i] = 0;
        }
        else if (v[i] > 0) {
            result[i] = 1;
        }
        else {
            result[i] = -1;
        }
    }
    return result;
}

plaintext_t multiplyPlaintexts(plaintext_t ptxt1, plaintext_t ptxt2, cryptocontext_t cc, int depth, int N) {
    auto ptxt1_vec = extract_vector(ptxt1);
    auto ptxt2_vec = extract_vector(ptxt2);
    //std::cout << ptxt1_vec.size() << ' ' << ptxt2_vec.size() << std::endl;
    auto res = multiplyVectors(ptxt1_vec, ptxt2_vec);

    auto res_ptxt = cc->MakeCKKSPackedPlaintext(res, depth, 0, nullptr, N/2);

    return res_ptxt;
}

plaintext_t decodePlaintext(plaintext_t ptxt1, cryptocontext_t cc, int depth, int N) {
    auto ptxt1_vec = extract_vector(ptxt1);

    auto res = decodeVector(ptxt1_vec);

    auto res_ptxt = cc->MakeCKKSPackedPlaintext(res, depth, 0, nullptr, N/2);

    return res_ptxt;
}

plaintext_t subPlaintexts(plaintext_t ptxt1, plaintext_t ptxt2, cryptocontext_t cc, int depth, int N) {
    auto ptxt1_vec = extract_vector(ptxt1);
    auto ptxt2_vec = extract_vector(ptxt2);

    auto res = subVectors(ptxt1_vec, ptxt2_vec);

    auto res_ptxt = cc->MakeCKKSPackedPlaintext(res, depth, 0, nullptr, N/2);

    return res_ptxt;
}

plaintext_t addPlaintexts(plaintext_t ptxt1, plaintext_t ptxt2, cryptocontext_t cc, int depth, int N) {
    auto ptxt1_vec = extract_vector(ptxt1);
    auto ptxt2_vec = extract_vector(ptxt2);

    auto res = addVectors(ptxt1_vec, ptxt2_vec);

    auto res_ptxt = cc->MakeCKKSPackedPlaintext(res, depth, 0, nullptr, N/2);

    return res_ptxt;
}

plaintext_t comparePlaintexts(plaintext_t ptxt1, plaintext_t ptxt2, cryptocontext_t cc, int depth, int N, int bound) {
    auto ptxt1_vec = extract_vector(ptxt1);
    auto ptxt2_vec = extract_vector(ptxt2);

    auto res = compareVectors(ptxt1_vec, ptxt2_vec, bound);

    auto res_ptxt = cc->MakeCKKSPackedPlaintext(res, depth, 0, nullptr, N/2);

    return res_ptxt;
}

plaintext_t copyPlaintext(plaintext_t ptxt1, cryptocontext_t cc, int depth, int N) {
    auto ptxt1_vec = extract_vector(ptxt1);
    auto res_ptxt = cc->MakeCKKSPackedPlaintext(ptxt1_vec, depth, 0, nullptr, N/2);
    return res_ptxt;
}


template <typename T> void print_ciphertext(ciphertext_t c, cryptocontext_t& cc, keypair_t keyPair, int length){
    Plaintext output;
    cc->Decrypt(keyPair.secretKey, c, &output);
    output->SetLength(length);

    int num_values = length;
    vec1d_t<double> r;

    for (int i = 0; i < num_values; i++){
        auto a = real(output->GetCKKSPackedValue()[i]);
        if(fabs(a) <= 1e-5) {
            r.push_back(0);
        } else {
            r.push_back(a);
        }
    }
    print_vec1d(r);
}

template <typename T> void print_ciphertext2(ciphertext_t c, cryptocontext_t& cc, keypair_t keyPair, int length){
    Plaintext output;
    cc->Decrypt(keyPair.secretKey, c, &output);
    output->SetLength(length);

    int num_values = length;
    vec1d_t<double> r;

    for (int i = 0; i < num_values; i++){
        auto a = real(output->GetCKKSPackedValue()[i]);
        if(fabs(a) <= 1e-5) {
            r.push_back(0);
        } else {
            r.push_back(1/a);
        }
    }
    print_vec1d(r);
}

template <typename T> void print_ciphertext_raw(ciphertext_t c, cryptocontext_t& cc, keypair_t keyPair, int length){
    Plaintext output;
    cc->Decrypt(keyPair.secretKey, c, &output);
    output->SetLength(length);
    std::cout << output << std::endl;
}

double add_elements_ciphertext(ciphertext_t c, cryptocontext_t& cc, keypair_t keyPair, int length){
    Plaintext output;
    cc->Decrypt(keyPair.secretKey, c, &output);
    output->SetLength(length);

    int num_values = length;
    double r = 0;

    for (int i = 0; i < num_values; i++){
        auto a = real(output->GetCKKSPackedValue()[i]);
        if(fabs(a) <= 1e-5) {
            r += 0;
        } else {
            r += a;
        }
        //std::cout << "i: " << i << "r: " << r << "a: " << a << std::endl;
    }
    return r;
}

double add_elements_ciphertext_decoded(ciphertext_t c, cryptocontext_t& cc, keypair_t keyPair, int length){
    Plaintext output;
    cc->Decrypt(keyPair.secretKey, c, &output);
    output->SetLength(length);

    int num_values = length;
    double r = 0;

    for (int i = 0; i < num_values; i++){
        auto a = real(output->GetCKKSPackedValue()[i]);
        if(fabs(a) <= 1e-5) {
            r += 0;
        } else {
            r += pow(2,1/a);
        }
        //std::cout << "i: " << i << "r: " << r << "a: " << a << std::endl;
    }
    return r;
}

template <typename T> vec1d_t<double> print_ciphertext_decode(ciphertext_t c, cryptocontext_t& cc, keypair_t keyPair, int length){
    Plaintext output;
    cc->Decrypt(keyPair.secretKey, c, &output);
    output->SetLength(length);

    int num_values = length;
    vec1d_t<double> r;

    for (int i = 0; i < num_values; i++){
        auto a = real(output->GetCKKSPackedValue()[i]);
        if(a < 1.0/pow(2,spliting)) { //2^splitting is the maximum value that we encode, if something is smaller it means it is a 0.
            r.push_back(0);
        } else {
            r.push_back(1/a);
        }
    }
    //print_vec1d(r);
    return r;
}

int compute_diff_pos(ciphertext_t ctxt1, ciphertext_t ctxt2, int len,
                     cryptocontext_t cc, keypair_t keys) {
  plaintext_t ptxt1;
  plaintext_t ptxt2;
  double max = -1;

  cc->Decrypt(keys.secretKey, ctxt1, &ptxt1);
  ptxt1->SetLength(len);
  cc->Decrypt(keys.secretKey, ctxt2, &ptxt2);
  ptxt2->SetLength(len);

  size_t slots = ptxt1->GetCKKSPackedValue().size();

  for (size_t i = 0; i < slots; i++) {
    auto diff = fabs(real(ptxt1->GetCKKSPackedValue()[i]) -
                     real(ptxt2->GetCKKSPackedValue()[i]));
    if (diff > max) {
      max = diff;
    }
  }

  return max;
}

std::vector<int> computePSU(const std::vector<std::vector<int>>& vecOfVecs) {
    std::unordered_set<int> unionSet;

    // Iterate through each vector in the vector of vectors
    for (const auto& vec : vecOfVecs) {
        // Insert all elements of the vector into the set
        unionSet.insert(vec.begin(), vec.end());
    }

    // Convert the set to a vector
    std::vector<int> result(unionSet.begin(), unionSet.end());

    return result;
}

std::vector<double> computePSU(const std::vector<std::vector<double>>& vecOfVecs) {
    std::unordered_set<double> unionSet;

    // Iterate through each vector in the vector of vectors
    for (const auto& vec : vecOfVecs) {
        // Insert all elements of the vector into the set
        unionSet.insert(vec.begin(), vec.end());
    }

    // Convert the set to a vector
    std::vector<double> result(unionSet.begin(), unionSet.end());

    return result;
}

int compute_diff_pos_ptxt(ciphertext_t ctxt, plaintext_t ptxt_in, int len,
                     cryptocontext_t cc, keypair_t keys) {
    plaintext_t ptxt;
    //double max = -1;

    cc->Decrypt(keys.secretKey, ctxt, &ptxt);
    //ptxt->SetLength(len);

    ptxt->SetLength(std::min(ptxt->GetLength(), ptxt_in->GetLength()));

    auto result = ptxt->GetCKKSPackedValue();
    auto expectedResult = ptxt_in->GetCKKSPackedValue();
/*
    if (result.size() != expectedResult.size()) {
        OPENFHE_THROW(config_error, "Cannot compare vectors with different numbers of elements");
    }
    */
    // using the Euclidean norm
    double avrg = 0;
    for (size_t i = 0; i < result.size(); ++i) {
        avrg += std::pow(std::abs(result[i].real() - expectedResult[i].real()), 2);
        if(std::pow(std::abs(result[i].real() - expectedResult[i].real()), 2) > 0.001) {
            std::cout << "Mitja: " << i << " " << std::pow(std::abs(result[i].real() - expectedResult[i].real()), 2) << " result: " << result[i].real()  << " expected: " << expectedResult[i].real() << std::endl;
        }
        if(i < 10){
            //std::cout << "ALPHA: " << log2(fabs(result[i].real() - expectedResult[i].real())) << " " << result[i].real() << " " <<  expectedResult[i].real() << " " << result[i].real() - expectedResult[i].real() << std::endl;
        }
    }

    avrg = std::sqrt(avrg) / result.size();  // get the average

    std::cout << "Precision: " << std::abs(std::log2(avrg)) << std::endl;

    return std::abs(std::log2(avrg));

int compute_diff(vec1d_t<ciphertext_t> &ctxt1, vec1d_t<ciphertext_t> &ctxt2,
                 int len, cryptocontext_t cc, keypair_t keys) {
  for (int i = 0; i < len; i++) {
    int diff = compute_diff_pos(ctxt1[i], ctxt2[i], len, cc, keys);
    if (diff > 0) {
      std::cout << "NO PASS" << std::endl;
      //exit(0);
    }
  }
  std::cout << "PASS" << std::endl;
  return true;
}

int compute_diff_ptxt(vec1d_t<ciphertext_t> &ctxt, vec1d_t<plaintext_t> &ptxt_in,
                 int len, cryptocontext_t cc, keypair_t keys) {
  for (int i = 0; i < len; i++) {
    int diff = compute_diff_pos_ptxt(ctxt[i], ptxt_in[i], cc->GetRingDimension()/2, cc, keys);
    if (diff < 10) {
      std::cout << "NO PASS" << std::endl;
      //exit(0);
    }
  }
  std::cout << "PASS" << std::endl;
  return true;
}

int compute_diff_ptxt_single(ciphertext_t &ctxt, plaintext_t &ptxt_in,
                 int len, cryptocontext_t cc, keypair_t keys) {
  for (int i = 0; i < len; i++) {
    int diff = compute_diff_pos_ptxt(ctxt, ptxt_in, cc->GetRingDimension()/2, cc, keys);
    if (diff < 10) {
      std::cout << "NO PASS" << std::endl;
      //exit(0);
    }
  }
  std::cout << "PASS" << std::endl;
  return true;
}


// -------------------------------------------------------
/// @brief extracts the psu-list of the ciphertext a by multiplying it by a mask
/// @param a ciphertext containing the list
/// @param M size of the resulting list
/// @param j bin to which this list corresponds
/// @param cryptoContext crypto context of the used scheme
/// @return returns a ciphertext containing the psu-list in the adequate slots
Ciphertext<DCRTPoly> GetListPSU(ciphertext_t& a, const size_t M, size_t j, cryptocontext_t& cc, keypair_t keyPair) {
    std::vector<double> v((j + 1) * M, 0);
    for (size_t k = 0; k < M; ++k) {
        v[j * M + k] = 1;
    }
    size_t N = cc->GetRingDimension();
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(v, a->GetNoiseScaleDeg(), a->GetLevel(), nullptr, N/2);

    std::vector<double> w((j + 1) * M, 1.0);
    Plaintext ptxt_ones = cc->MakeCKKSPackedPlaintext(w, a->GetNoiseScaleDeg(), a->GetLevel(), nullptr, N/2);

    a = cc->EvalSub(ptxt_ones, a); //ends in nan if we perform sign in the EvalLinearTransform function
    
    Ciphertext<DCRTPoly> list = cc->EvalMult(a, ptxt);
    num_multiplications_getlist += 1;
    //std::cout << "GetListPsu" << num_multiplications_getlist << std::endl;

    return list;
}

plaintext_t GetListPSU_plaintext(plaintext_t a, const size_t M, size_t j, cryptocontext_t cc, int depth, int N, int bound) {
    auto a_vec = extract_vector(a);
    auto s = a_vec.size();

    std::vector<double> w((j + 1) * M, 1.0);

    while(w.size() < s){
        w.push_back(0);
    }

    Plaintext ptxt_ones = cc->MakeCKKSPackedPlaintext(w, depth, 0, nullptr, N/2);

    a = subPlaintexts(ptxt_ones, a, cc, depth, N);

    auto a_aux = extract_vector(a);

    for(int i = bound; i < N/2; i++){
        a_aux[i] = 0;
    }

    a = cc->MakeCKKSPackedPlaintext(a_aux, depth, 0, nullptr, N/2);
    std::vector<double> v(N/2, 0);
    for (size_t k = 0; k < M; ++k) {
        v[j * M + k] = 1;
    }

    plaintext_t ptxt = cc->MakeCKKSPackedPlaintext(v, depth, 0, nullptr, N/2);

    plaintext_t list = multiplyPlaintexts(a, ptxt, cc, depth, N);

    return list;
}

std::vector<double> evaluatePolynomial(const std::vector<double>& coefficients, const std::vector<double>& x_values) {
    std::vector<double> results;

    for (double x : x_values) {
        int degree = coefficients.size() - 1;  // Degree of the polynomial

        double result = coefficients[degree];
        for (int i = degree - 1; i >= 0; --i) {
            result = result * x + coefficients[i];
        }

        results.push_back(result);
    }

    return results;
}

Ciphertext<DCRTPoly> HomomorphicSign(Ciphertext<DCRTPoly>& x, size_t dg, size_t df, const cryptocontext_t& cc, keypair_t keyPair) {
    //std::vector<double> coef_g = {0.0, 2126.0/1024, 0.0, -1359.0/1024}; //n=1
    //std::vector<double> coef_g = {0.0, 3334.0/1024, 0.0, -6108.0/1024, 0.0, 3796.0/1024}; //n=2
    //std::vector<double> coef_g = {0.0, 4589.0/1024, 0.0, -16577.0/1024, 0.0, 25614.0/1024, 0.0, -12860.0/1024}; //n=3
    std::vector<double> coef_g = {0.0, 5850.0/1024, 0.0, -34974.0/1024, 0.0, 97015.0/1024, 0.0, -113492.0/1024, 0.0, 46623.0/1024}; //n=4
    
    // Plaintext ptxt;
    // auto N = cc->GetRingDimension();
    //#pragma omp parallel for
    for (size_t i = 0; i < dg; ++i) {
        
        // cc->Decrypt(keyPair.secretKey, x, &ptxt);
        // ptxt->SetLength(N/2);
        // auto output_vec = extract_vector(ptxt);
        // auto polyeval = evaluatePolynomial(coef_g, output_vec);
        // auto polyeval_ptxt = cc->MakeCKKSPackedPlaintext(polyeval, 63, 0, nullptr, N/2); 

        Ciphertext<DCRTPoly> d = cc->EvalPoly(x, coef_g);
        //std::cout << "Depth x before " << d->GetLevel() << ", i: " << i << std::endl;
        if(d->GetLevel() >= 21) {
            d = cc->EvalBootstrap(d,2,17);
        }
        x = d;
        //std::cout << "Depth x after" << x->GetLevel() << ", i: " << i << std::endl;
        // compute_diff_ptxt_single(x, polyeval_ptxt, 1, cc, keyPair); 

    //     // print_ciphertext<double>(x, cc, keyPair, N/2);
    //     // std::cout << polyeval_ptxt << std::endl;
        

    }
    //x = cc->EvalBootstrap(x,2,17);
    //std::vector<double> coef_f = {0.0, 3.0/2, 0.0, -1.0/2}; //n=1
    //std::vector<double> coef_f = {0.0, 15.0/8, 0.0, -10.0/8, 0.0, 3.0/8}; //n=2
    //std::vector<double> coef_f = {0.0, 35.0/16, 0.0, -35.0/16, 0.0, 21.0/16, 0.0, -5.0/16}; //n=3
    std::vector<double> coef_f = {0.0, 315.0/128, 0.0, -420.0/128, 0.0, 378.0/128, 0.0, -180.0/128, 0.0, 35.0/128}; //n=4
    //std::vector<double> coef_f = {0.0, 693.0/256, 0.0, -1155.0/256, 0.0, 1386.0/256, 0.0, -990.0/256,  0.0, 385.0/256, 0.0, -63.0/256}; //n=5
    // std::vector<double> coef_f = {0.0, 3003.0/1024, 0.0, -6006.0/1024, 0.0, 9009.0/1024, 0.0, -8580.0/1024, 0.0, 5005.0/1024, 0.0, -1638.0/1024, 0.0, 231.0/1024}; //n=6
    //#pragma omp parallel for
    for (size_t i = 0; i < df; ++i) {
        
        // cc->Decrypt(keyPair.secretKey, x, &ptxt);
        // ptxt->SetLength(N/2);
        // auto output_vec = extract_vector(ptxt);
        // auto polyeval = evaluatePolynomial(coef_f, output_vec);
        // auto polyeval_ptxt = cc->MakeCKKSPackedPlaintext(polyeval, 63, 0, nullptr, N/2); 
        
        Ciphertext<DCRTPoly> d = cc->EvalPoly(x, coef_f);
        if(d->GetLevel() >= 21) {
            d = cc->EvalBootstrap(d,2,17);
        }
        x = d; 
        //std::cout << "Depth x " << x->GetLevel() << ", i: " << i << std::endl;
        // compute_diff_ptxt_single(x, polyeval_ptxt, 1, cc, keyPair);  
        // print_ciphertext<double>(x, cc, keyPair, N/2);
        // std::cout << polyeval_ptxt << std::endl;  
        
    }
    return x;
}

plaintext_t signPlaintext(plaintext_t& ptxt, cryptocontext_t& cc, int depth, int N) {
    auto ptxt_vec = extract_vector(ptxt);

    auto res = signVector(ptxt_vec);

    auto res_ptxt = cc->MakeCKKSPackedPlaintext(res, depth, 0, nullptr, N/2);

    return res_ptxt;
}

/// @brief evaluates a Homomorphic Comparison
/// @param a first ciphertext to compare
/// @param b second ciphertext to compare
/// @param dg order of the composition of the g dfunction in the sigmoid approximation
/// @param df order of the composition of the f dfunction in the sigmoid approximation
/// @param cryptoContext crypto context of the used scheme
/// @return returns a ciphertext containing encryption of 1 in the k-th slot if a[k] == b[k], encryption of 0 otherwise
Ciphertext<DCRTPoly> HomomorphicComparison(Ciphertext<DCRTPoly>& x, Ciphertext<DCRTPoly>& y, const size_t& dg, const size_t& df, cryptocontext_t cc, keypair_t keyPair, int N, size_t i, size_t j, size_t n, size_t b) {
    Ciphertext<DCRTPoly> c = cc->EvalSub(x, y);

    // Plaintext out;
    // cc->Decrypt(keyPair.secretKey, c, &out);
    // out->SetLength(i*std::ceil((int)n/b)*std::ceil((int)n/b));
    // std::cout << "c = " << out << std::endl;

    //std::cout << "C depth" << c->GetLevel() << std::endl;

    Ciphertext<DCRTPoly> d = HomomorphicSign(c, dg, df, cc, keyPair);

    //std::cout << "D depth" << d->GetLevel() << std::endl;

    c = cc->EvalSquare(d);

    //std::cout << "C^2 depth" << c->GetLevel() << std::endl;

    // std::cout << "after square" << std::endl;
    // cc->Decrypt(keyPair.secretKey, c, &out);
    // out->SetLength(i*std::ceil((int)n/b)*std::ceil((int)n/b));
    // std::cout << "c = " << out << std::endl;
    std::vector<double> v(i*bin_size*bin_size, 1.0);
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(v, c->GetNoiseScaleDeg(), c->GetLevel(), nullptr, N/2);
    Ciphertext<DCRTPoly> e = cc->EvalSub(ptxt, c);

    //std::cout << "e depth" << e->GetLevel() << std::endl;

    return e;
}


size_t sha256_mod_b(int input, size_t range) {
    // Extract the first 5 bits of 'input'
    //int first5Bits = input & 0b11111;
    int first5Bits = input & ((1 << spliting) - 1);

    // Extract the remaining bits
    //int remainingBits = input >> spliting;

    // Extract first log2(range) bits
    int resultBits = (input >> spliting) & ((1 << static_cast<int>(log2(range))) - 1);

    // Hash the first 5 bits
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_CTX sha256Context;
    std::string first5BitsString = std::to_string(first5Bits);
    SHA256_Init(&sha256Context);
    SHA256_Update(&sha256Context, first5BitsString.c_str(), first5BitsString.length());
    SHA256_Final(hash, &sha256Context);

    // Convert the first sizeof(size_t) bytes of the hash to a size_t value
    size_t hashValue;
    std::memcpy(&hashValue, hash, sizeof(size_t));

    hashValue = hashValue % range;

    // XOR with the remaining bits
    size_t xorResult = hashValue ^ static_cast<size_t>(resultBits);

    return xorResult % range;
    // Apply modulo to limit the value within the specified range
    //size_t result = xorResult % range;

    //return result;
}

std::vector<std::vector<std::vector<double>>> ClientSetup_Clear_Bins(vec1d_t<int> list, int i, size_t B, int n, int N, vec2d_t<double>& ages_bins, vec2d_t<double>& bmi_bins, vec2d_t<double>& ones_bins) {
    
    std::vector<std::vector<std::vector<double>>> w;
    std::cout << std::endl;
    std::cout << "Party " << i << std::endl;
    std::vector<std::vector<double>> v(B);
    std::vector<std::vector<double>> v2(B);

    std::uniform_int_distribution<> dist(2, 90); // for age
    std::uniform_int_distribution<> dist2(13, 38); // for BMI
    std::random_device rd;
    std::mt19937 gen(rd());

    //std::cout << "size of list " << list.size() << std::endl;

    double aa = 2;
    double b = pow(2,spliting);
    double c = 0.0100001;
    double d = 1;

    for (int& a : list) {
        size_t h = sha256_mod_b(a, B);
        v[h].push_back(((c - d)/(aa - b))*a + (aa*d - b*c)/(aa-b));
        v2[h].push_back(1.0/a);

        ages_bins[h].push_back(1.0/log2(dist(gen)));
        bmi_bins[h].push_back(dist2(gen));
        ones_bins[h].push_back(1);
    }

    w.push_back(v);
    w.push_back(v2);

    return w;
}


Ciphertext<DCRTPoly> ClientSetup(std::vector<double> a, int depth, int N, keypair_t keys, cryptocontext_t cc) {
    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(a, depth, 0, nullptr, N/2);
    auto ctxt = cc->Encrypt(keys.publicKey, ptxt);
    return ctxt;
}

Ciphertext<DCRTPoly> GetCopy(Ciphertext<DCRTPoly>& a, size_t j, const size_t n, const size_t b, usint& depth, const size_t N, const CryptoContext<DCRTPoly>& cryptoContext) {
    std::vector<double> u(j * ceil((double)n / b), 0);
    std::vector<double> v(ceil((double)n / b), 1.0);
    u.reserve(u.size() + v.size());
    u.insert(u.end(), v.begin(), v.end());
    Plaintext mask = cryptoContext->MakeCKKSPackedPlaintext(u, a->GetNoiseScaleDeg(), a->GetLevel(), nullptr, N/2);
    Ciphertext<DCRTPoly> ctxt = cryptoContext->EvalMult(a, mask);
    return ctxt;
}

void vecToFullPacking(cryptocontext_t cc, vec1d_t<int> indexlist, ciphertext_t c, keypair_t keyPair, int N) {
    //Plaintext output;

    for(int i = 0; i < (int)indexlist.size(); i++){
        //std::cout << "rotation" << indexlist[i] << std::endl;
        auto a = cc->EvalRotate(c,indexlist[i]);
        cc->EvalAddInPlace(c,a);
        //cc->Decrypt(keyPair.secretKey, c, &output);
        //output->SetLength(N/2);
        //std::cout << output << std::endl;
    }
}

ciphertext_t GetRot(cryptocontext_t cc, keypair_t keys, vec1d_t<int> packingindex, vec1d_t<int> rotationindex, vec1d_t<int> partiesindex, ciphertext_t c, int N, int n, int depth, int index) {
    
    //Precomputation for the hoisted automorphisms
    auto cPrecomp = cc->EvalFastRotationPrecompute(c);

    Ciphertext<DCRTPoly> rot;

    //Generate rot using fast rotations
    for(int i = 0; i < (int)rotationindex.size(); i++){
        std::vector<double> mask(N/2, 0.0);
        //std::cout << "rotation" << rotationindex[i] << std::endl;
        auto a = cc->EvalFastRotation(c, rotationindex[i], 2*N,cPrecomp);
        for(int j = 0; j < n; j++){
            mask[n*rotationindex[i] + j] = 1.0;
        }
        Plaintext ptxt_mask = cc->MakeCKKSPackedPlaintext(mask, a->GetNoiseScaleDeg(), a->GetLevel(), nullptr, N/2);
        if(i == 0){
            rot = cc->EvalMult(a, ptxt_mask);
            num_multiplications_getrot += 1;
            //std::cout << "GetRot: " << num_multiplications_getrot << std::endl;
        } else {
            auto b = cc->EvalMult(a, ptxt_mask);
            num_multiplications_getrot += 1;
            //std::cout << "GetRot: " << num_multiplications_getrot << std::endl;
            cc->EvalAddInPlace(rot, b);
            //Plaintext output;
            //cc->Decrypt(keys.secretKey, b, &output);
            //output->SetLength(n*n);
            //std::cout << output << std::endl;
        }
    }
    
    rot = cc->EvalRotate(rot, partiesindex[index]);

    return rot;
}

ciphertext_t GetNonRot(cryptocontext_t cc, keypair_t keys, vec1d_t<int> packingindex, vec1d_t<int> rotationindex, vec1d_t<int> partiesindex, ciphertext_t c, int N, int n, int depth, int index) {
    
    //Precomputation for the hoisted automorphisms
    auto cPrecomp = cc->EvalFastRotationPrecompute(c);

    Ciphertext<DCRTPoly> rot;

    //Generate rot using fast rotations
    for(int i = 0; i < (int)rotationindex.size(); i++){
        std::vector<double> mask(N/2, 0.0);
        auto a = cc->EvalFastRotation(c, rotationindex[i], 2*N, cPrecomp);
        for(int j = 0; j < n; j++){
            mask[n*rotationindex[i] + j] = 1.0;
        }
        Plaintext ptxt_mask = cc->MakeCKKSPackedPlaintext(mask, a->GetNoiseScaleDeg(), a->GetLevel(), nullptr, N/2);
        if(i == 0){
            rot = cc->EvalMult(a, ptxt_mask);
            num_multiplications_getnrot += 1;
            std::cout << "GetNonRot: " << num_multiplications_getnrot << std::endl;
        } else {
            auto b = cc->EvalMult(a, ptxt_mask);
            num_multiplications_getnrot += 1;
            std::cout << "GetNonRot: " << num_multiplications_getnrot << std::endl;
            cc->EvalAddInPlace(rot, b);
        }
    }
    
    rot = cc->EvalRotate(rot, partiesindex[index]);

    return rot;
}

/// @brief computes b
/// @param n size of the sets
/// @param P number of sets (parties)
/// @param N ring dimension
/// @return returns the minimum b such that bin_size^2 * (P-1) <= N/2
size_t find_B(const size_t n, const size_t P, const size_t N) {
    int B = 1;
    //std::cout << "find b:" << pow(pow(2, ceil(log2(ceil((double)n/B)))), 2) * (P - 1) << std::endl;
    while (pow(pow(2, ceil(log2(ceil((double)n/B)))), 2) * (P - 1) >= N/2) {
        ++B;
    }
    return B;
}

std::vector<int> createArray(int n) {
    std::vector<int> arr;
    for (int i = 1; i <= n; ++i) {
        arr.push_back(i);
    }
    return arr;
}

const vec2d_t<int> genRotationArrays(int num_Blocks, int cutting_point) {
  auto closest_power2 = 1 << (int)floor(log2(cutting_point * num_Blocks)) + 1;

  vec2d_t<int> rotation_arrays;

  for (int i = 1; i < num_Blocks; i++) {
    vec1d_t<int> array(closest_power2);

    std::iota(array.begin(), array.end(), 1);

    for (int j = 0; j < cutting_point; j++) {
      std::swap(array[j], array[j + i * cutting_point]);
    }

    rotation_arrays.push_back(array);
  }

  return rotation_arrays;
}



const vec1d_t<int> getAllRotations(const std::vector<int>& arr, int total) {
    vec1d_t<int> rotation;
    int n = arr.size();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            rotation.push_back(arr[(i + j) % n] + i*n);
        }
    }

    for (int i = n*n; i < total; ++i) {
    rotation.push_back(i);
    }

    return rotation;
}

std::vector<double> getAllRotations_plaintext(const std::vector<double>& v, int i, int n, int B, int N) {
    std::vector<double> rotations;
    int vector_size = v.size();

    //TODO: FIX THIS
    if((int)vector_size > bin_size) {
        //std::cout << "Error in size (n is too small compared to the range of input values)" << std::endl;
        std::cout << "real size " << vector_size << " estimated size= " << bin_size << std::endl;
        exit(1);
    }

    //std::cout << "Vector size: " << vector_size << std::endl;

    int zeros = (i)*pow(bin_size,2);

    //std::cout << "Number of zeros: " << (i)*pow(bin_size,2) << std::endl;

    for (int i = 0; i < zeros; i++) {
        rotations.push_back(0);
    }

    for (int i = 0; i < vector_size; i++) {
        for (int j = 0; j < vector_size; j++) {
            rotations.push_back(v[(i + j) % vector_size]);
        }
    }

    int current_length = rotations.size();

    //std::cout << "current length: " << current_length << std::endl;

    for(int i = current_length; i < N/2; i++){
        rotations.push_back(0);
    }

    while ((int)rotations.size() > N/2) { // TODO: REVISE, should not be needed
        rotations.pop_back();
    }

    //std::cout << "GetRotPtxt size = " << rotations.size() << std::endl;

    return rotations;
}

std::vector<double> getAllRotations_plaintext2(const std::vector<double>& v, int i, int n, int B, int N) {
    std::vector<double> rotations;
    int vector_size = v.size();

    std::cout << "Vector size: " << vector_size << std::endl;

    int zeros = (i)*pow(bin_size,2);

    std::cout << "Number of zeros: " << (i)*pow(bin_size,2) << std::endl;

    for (int i = 0; i < zeros; i++) {
        rotations.push_back(0);
    }

    for (int i = 0; i < vector_size; i++) {
        for (int j = 0; j < vector_size; j++) {
            rotations.push_back(v[(i + j) % vector_size]);
        }
    }

    int current_length = rotations.size();

    for(int i = current_length; i < N/2; i++){
        rotations.push_back(0);
    }

    return rotations;
}

const vec2d_t<complex_t> genPermutation(vec1d_t<int> array, size_t n) {
  vec2d_t<complex_t> matrix(n, vec1d_t<complex_t>(n));
  for (unsigned int t = 0; t < n; t++) {
    matrix[t][array[t] - 1] = 1;
  }

  return matrix;
}

int compute_error(cryptocontext_t cc, keypair_t keys, plaintext_t ptxt,
                  ciphertext_t ctxt, int len) {
  plaintext_t out;
  double max = -1;

  cc->Decrypt(keys.secretKey, ctxt, &out);
  out->SetLength(len);

  for (int i = 0; i < len; i++) {
    auto diff = fabs(real(ptxt->GetCKKSPackedValue()[i]) -
                    real(out->GetCKKSPackedValue()[i]));
    if (diff > max) {
      max = diff;
    }
  }

  return (int)log2(ceil(max));
}

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

bool checkBins(const std::vector<std::vector<double>>& bins) {
    for (const auto& vec : bins) {
        if(vec.size() > (size_t)bin_size) {
            return false; // At least one element is not smaller than bin_size
        }
    }
    return true; // All elements are smaller than bin_size
}

void print_moduli_chain(const DCRTPoly& poly){
    int num_primes = poly.GetNumOfElements();
    double total_bit_len = 0.0;
    for (int i = 0; i < num_primes; i++) {
        auto qi = poly.GetParams()->GetParams()[i]->GetModulus();
        std::cout << "q_" << i << ": " 
                    << qi
                    << ",  log q_" << i <<": " << log(qi.ConvertToDouble()) / log(2)
                    << std::endl;
        total_bit_len += log(qi.ConvertToDouble()) / log(2);
    }   
    std::cout << "Total bit length: " << total_bit_len << std::endl;
}

#include <omp.h>
int main() {


  fheckkrns_t ckksrns;

  num_multiplications_getlist = 0;
num_multiplications_getrot = 0;
num_multiplications_getnrot = 0;
num_multiplications_server = 0;
    
  vec1d_t<uint32_t> dim1 = {0, 0};

  // Level budget = {e,d} controls the computational complexity of the
  // homomorphic encoding and decoding steps in CKKS bootstrapping, both are
  // homomorphic evaluations of linear transforms. A higher budget allows faster
  // computation of these steps at the expense of using a deeper circuit
  // (consuming more levels). On the other hand, lower budget would be slower
  // but it uses a shallow circuit. Recommended values, found experimentally,
  // are e or d = {1,2,3, or 4} (they do not need to be equal)
  vec1d_t<uint32_t> levelBudget = {4, 4};

#if NATIVEINT == 128
  lbcrypto::ScalingTechnique rescaleTech = FLEXIBLEAUTOEXT;
  usint dcrtBits = 78;
  usint firstMod = 89; /*firstMod*/
#else
  lbcrypto::ScalingTechnique rescaleTech = lbcrypto::FLEXIBLEAUTOEXT;
  usint dcrtBits = 59;
  usint firstMod = 60; /*firstMod*/
#endif

  // lbcrypto::FLEXIBLEAUTOEXT

  auto secretKeyDist =
      lbcrypto::UNIFORM_TERNARY; // Check section 3 to see why uniform ternary
                                 // makes sense:
                                 // https://eprint.iacr.org/2020/1118.pdf
  //auto n = 1 << 15;

  //usint depth = 57;
    int dg = 8;
    int df = 2;

    const size_t P = 2; //#parties
    const size_t n = 1 << 12; // 1 << 3
    auto size_of_permutation = n;

    std::random_device rd;
    std::mt19937 gen(rd());
    //const size_t A = 1 << 16; //range of the sets // 1 << 16
    //std::uniform_int_distribution<> distrib(1, A);

    const size_t A = 1 << int(ceil(P*log2(n))); //range of the sets, we assume all elements are different
    //const size_t A = 1 << 20; //range of the sets, we assume all elements are different
    spliting = log2(A);
    std::cout << "A: " << log2(A) << std::endl;
    std::cout << "split: " << spliting << std::endl;
    std::uniform_int_distribution<> distrib(2, A); //We assume the indices are mapped to the range [2, A] 

    uint32_t levelsAvailableAfterBootstrap = 7;
    uint32_t margin = 0;
    usint depth = levelsAvailableAfterBootstrap + (FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist)) + margin;
    std::cout << "BootstrapDepth: " << (FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist)) << std::endl;
    std::cout << "Depth: " << depth << std::endl;

    CryptoContext<DCRTPoly> cc;

    
  lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> parameters;

  parameters.SetMultiplicativeDepth(depth);
  parameters.SetScalingModSize(dcrtBits);
  parameters.SetFirstModSize(firstMod);
  parameters.SetScalingTechnique(rescaleTech);
  parameters.SetSecretKeyDist(secretKeyDist);
  parameters.SetNumLargeDigits(3);
  parameters.SetKeySwitchTechnique(lbcrypto::HYBRID);

  parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);

  cc = GenCryptoContext(parameters);

  // Turn on features
  cc->Enable(lbcrypto::PKE);
  cc->Enable(lbcrypto::KEYSWITCH);
  cc->Enable(lbcrypto::LEVELEDSHE);
  cc->Enable(lbcrypto::ADVANCEDSHE);
  cc->Enable(lbcrypto::FHE);

  

 auto N = cc->GetRingDimension();
 std::cout << "N = " << N << std::endl;
 //std::cout << *cc->GetCryptoParameters()  <<std::endl;
 NativeInteger q = cc->GetElementParams()->GetParams()[0]->GetModulus().ConvertToInt();
 std::cout << "q = " << q << std::endl;

    size_t B = find_B(n, P, N);
    //minimum B
    std::cout << "B min: " << B << std::endl;
    bin_size = pow(2, ceil(log2(ceil((double)n/B))));
    B = 52;
    std::cout << "B used: " << B << std::endl;

    uint32_t numSlots = N/2;
    uint32_t M     = cc->GetCyclotomicOrder();
    uint32_t slots = (numSlots == 0) ? M / 4 : numSlots;
    ckksrns.m_bootPrecomMap[slots] = std::make_shared<CKKSBootstrapPrecom>();

    keypair_t keyPair;

    std::cout << "Key Generation";
    {
    Benchmark bench;
    keyPair = cc->KeyGen();
    }

    ////////////////////
    // CLIENT WORK
    ////////////////////

    //generate one random set for each party
    std::vector<std::vector<int>> parties(P);
    std::vector<std::set<int>> set_parties(P);

    for (size_t i = 0; i < P; ++i) {
        for (size_t j = 0; j < (size_t)size_of_permutation; ++j) {
            set_parties[i].insert(distrib(gen));
        }
    }
    
    for (size_t i = 0; i < P; ++i) {
        for (auto it = set_parties[i].begin(); it != set_parties[i].end(); ++it) {
            parties[i].push_back((int)*it);
        }
    }

    // for (size_t i = 0; i < P; ++i) {
    //     std::cout << "party " << i << " " << parties[i] << std::endl;
    // }

    Plaintext output;

    std::vector<ciphertext_t> parties_ciphertexts;
    std::vector<ciphertext_t> parties_ciphertexts2;
    std::vector<plaintext_t> parties_plaintexts;
    std::vector<plaintext_t> parties_plaintexts2;
    std::vector<std::vector<double>> parties_plaintexts_vec;
    std::vector<std::vector<std::vector<double>>> parties_cleartexts; 
    //every party encrypts its dataset into one ciphertext

    //for the PCSU
    vec2d_t<double> ages(P);
    vec2d_t<double> bmi(P);
    vec2d_t<double> ones_parties(P);
    vec1d_t<plaintext_t> ages_plaintext(P);
    vec1d_t<ciphertext_t> ages_ciphertext(P);

    vec1d_t<plaintext_t> bmi_plaintext(P);
    vec1d_t<ciphertext_t> bmi_ciphertext(P);

    vec1d_t<plaintext_t> ones_parties_plaintext(P);
    vec1d_t<ciphertext_t> ones_parties_ciphertext(P);

std::cout << "Client setup (all parties)";
{
  Benchmark bench;

        for(int i = 0; i < (int)P; i++){
            vec2d_t<double> ages_bins(B);
            vec2d_t<double> bmi_bins(B);
            vec2d_t<double> ones_bins(B);

            auto split = ClientSetup_Clear_Bins(parties[i],i,B,size_of_permutation,N, ages_bins, bmi_bins, ones_bins); //REMOVE THIS WHEN BENCHMARKING

            auto bins = split[0];
            auto bins2 = split[1];

            size_t totalSize = B * bin_size;

            std::vector<double> a;
            a.reserve(totalSize);
            std::vector<double> b;
            b.reserve(totalSize);
            std::vector<double> c;
            c.reserve(totalSize);
            std::vector<double> d;
            d.reserve(totalSize);
            for (int j = 0; j < (int)B; ++j) {
                //std::cout << "size of bin: " << j << " " << bins[j].size() << std::endl;
                //std::cout << bins[j] << std::endl;
                while ((int)bins[j].size() < bin_size) {
                    bins[j].push_back(0);
                    bmi_bins[j].push_back(0);
                    ages_bins[j].push_back(0);
                    ones_bins[j].push_back(0);
                }
                //std::cout << "size of bins" << bins[j].size() << std::endl;
                a.insert(a.end(), bins[j].begin(), bins[j].end());
                b.insert(b.end(), ages_bins[j].begin(), ages_bins[j].end());
                c.insert(c.end(), bmi_bins[j].begin(), bmi_bins[j].end());
                d.insert(d.end(), ones_bins[j].begin(), ones_bins[j].end());
            }
            ages[i] = b;
            bmi[i] = c;
            ones_parties[i] = d;
            std::vector<double> aa;
            aa.reserve(totalSize);
            for (int j = 0; j < (int)B; ++j) {
                while ((int)bins2[j].size() < bin_size) {
                    bins2[j].push_back(0);
                }
                //std::cout << "size of bins" << bins2[j].size() << std::endl;
                aa.insert(aa.end(), bins2[j].begin(), bins2[j].end());
            }

            parties_cleartexts.push_back(bins);
            auto ctp = ClientSetup(a, depth, N, keyPair, cc);
            parties_ciphertexts.push_back(ctp);

            auto ctp2 = ClientSetup(aa, depth, N, keyPair, cc);
            parties_ciphertexts2.push_back(ctp2);

            //std::cout << "A: " << a << std::endl;

            parties_plaintexts_vec.push_back(decodeVector(a));

            Plaintext ptxt = cc->MakeCKKSPackedPlaintext(a, depth, 0, nullptr, N/2);
            parties_plaintexts.push_back(ptxt);

            Plaintext ptxt2 = cc->MakeCKKSPackedPlaintext(aa, depth, 0, nullptr, N/2);
            parties_plaintexts2.push_back(ptxt2);

            ages_plaintext[i] = cc->MakeCKKSPackedPlaintext(ages[i], depth, 0, nullptr, N/2);
            ages_ciphertext[i] = cc->Encrypt(keyPair.publicKey, ages_plaintext[i]);

            bmi_plaintext[i] = cc->MakeCKKSPackedPlaintext(bmi[i], depth, 0, nullptr, N/2);
            bmi_ciphertext[i] = cc->Encrypt(keyPair.publicKey, bmi_plaintext[i]);        

            ones_parties_plaintext[i] = cc->MakeCKKSPackedPlaintext(ones_parties[i], depth, 0, nullptr, N/2);
            ones_parties_ciphertext[i] = cc->Encrypt(keyPair.publicKey, ones_parties_plaintext[i]);

        }
}
std::cout << "\t - Total number ciphertexts to server: " << parties_ciphertexts.size() + 3 << std::endl;
 
    std::cout << "Size of databases: " << n << std::endl;
    std::cout << "Bits per index: " << log2(A) << std::endl;
    std::cout << "Number of bins (B): " << B << std::endl;
    std::cout << "Elements per bin: " << bin_size << std::endl;

    //GET APPROXIMATE APPROXIMATE TIMMINGS BASED ON MICROBENCHMARKS
    auto num_rot_keys = P + log2(N) - log2(bin_size) - 1 + bin_size + 2*B;
    auto mult_by_mask = 100;
    auto rot = 495;
    //auto precomp = 120;

    std::cout << "Number of rotation keys: " << num_rot_keys << std::endl;
    std::cout << "Time nrot: " << P*(mult_by_mask + (B-1)*(rot+mult_by_mask)) << std::endl;

    std::vector<int32_t> partiesindex = {};

    for(int i = 1; i <= (int)P; i++){
        partiesindex.push_back((1-i)*bin_size*bin_size);
    }

    auto numRotations = log2(N) - log2(bin_size) - 1; //We assume that bin_size is a power of 2
    std::vector<int32_t> packingindex = {};

    for(int i = 0; i < numRotations; i++){
        packingindex.push_back(bin_size*(1 << i));
    }

    std::vector<int32_t> rotationindex = {};

    for(int i = 0; i < bin_size; i++){
        rotationindex.push_back(i);
            //addIfNotPresent(packingindex,i);
    }

    std::vector<int32_t> V = {};
    for(int i = 1; i <= (int)B; i++){
        V.push_back(i*bin_size);
        V.push_back(-i*bin_size);
    }
    

size_t num_rot = partiesindex.size() + packingindex.size() + rotationindex.size() + V.size();
std::cout << "Number of rotation keys: " << num_rot << std::endl;

std::cout << "Bootstrapping and rotation key generation";
{
  Benchmark bench;
  cc->EvalMultKeyGen(keyPair.secretKey);

  //int evalfunction = 1; //if 0 we don't evaluate a function over PSU

  //if(P > 2 || evalfunction == 1) { //we do not need bootstrapping if P = 2 and we do not evaluate any function after PSU
    cc->EvalBootstrapSetup(levelBudget, dim1, N/2);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, N/2);
  //}

    cc->EvalRotateKeyGen(keyPair.secretKey, partiesindex);

    cc->EvalRotateKeyGen(keyPair.secretKey, packingindex);

    cc->EvalRotateKeyGen(keyPair.secretKey, rotationindex);

    cc->EvalRotateKeyGen(keyPair.secretKey, V);

}

    ////////////////////
    // SERVER WORK IN THE CLEAR
    ////////////////////
    std::vector<std::vector<plaintext_t>> rotated_plaintext;

    for(int i = 0; i < (int)P; i++){
        
        std::vector<plaintext_t> bins;
        bins.resize(B*sizeof(plaintext_t));

        #pragma omp parallel for
        for(int j = 0; j < (int)B; j++) {
            bins[j] = cc->MakeCKKSPackedPlaintext(getAllRotations_plaintext(parties_cleartexts[i][j], i, n, B, N), depth, 0, nullptr, N/2);
        }

        rotated_plaintext.push_back(bins);

    }

    std::vector<std::vector<plaintext_t>> nrotated_plaintext;

    for(int i = 0; i < (int)P; i++){
        
        std::vector<plaintext_t> bins;

        for(int j = 0; j < (int)B; j++) {

            std::vector<double> vec;
            uint s = vec.size();
            while(s < N/2) {
                vec.insert(vec.end(), parties_cleartexts[i][j].begin(), parties_cleartexts[i][j].end());
                //std::cout << "size of vec" << vec.size() << " " << parties_cleartexts[i][j].size() << std::endl;
                s = vec.size();
                while (s%((int)(bin_size)) != 0) { 
                    vec.push_back(0);
                    s++;
                }

            } 

            auto nrotated_bin_plaintext = cc->MakeCKKSPackedPlaintext(vec, depth, 0, nullptr, N/2);
            bins.push_back(nrotated_bin_plaintext);
        }

        nrotated_plaintext.push_back(bins);
    }

    ////////////////////
    // SERVER WORK
    ////////////////////

    std::vector<double> mask(N/2,0.0);
    for(int i = 0; i < bin_size; i++){
        mask[i] = 1.0;
    }
    Plaintext ptxt_mask = cc->MakeCKKSPackedPlaintext(mask, depth, 0, nullptr, N/2);

    vec1d_t<plaintext_t> nrot_masks;

    for(int i = 1; i < (int)P; i++){
        std::vector<double> mask(N/2,0.0);
        //std::cout << "bound: " << i*bin_size*bin_size << std::endl;
        for(int j = 0; j < i*bin_size*bin_size; j++){
            mask[j] = 1.0;
        }
        Plaintext ptxt_mask_nrot = cc->MakeCKKSPackedPlaintext(mask, depth, 0, nullptr, N/2);
        nrot_masks.push_back(ptxt_mask_nrot);
    }
    
    //extract non rotated copies of each ciphertext
    vec2d_t<ciphertext_t> nrotated;

    int i = 0;
std::cout << "Extract nrot"; //TODO: Remove nrotated from 1st party
{
  Benchmark bench;
    for(ciphertext_t& ctp : parties_ciphertexts) {
        //auto ctpPrecomp = cc->EvalFastRotationPrecompute(ctp);
        vec1d_t<ciphertext_t> encrypted_bins;

        auto bin0 = cc->EvalMult(ctp,ptxt_mask);

        encrypted_bins.push_back(bin0);

        for(int j = 1; j < (int)B; j++){
            //std::cout << "Server: " << j << std::endl;
            //auto a = cc->EvalFastRotation(ctp,j*bin_size,2*N,ctpPrecomp);
            auto a = cc->EvalRotate(ctp,j*bin_size);
            auto bin = cc->EvalMult(a,ptxt_mask);

            encrypted_bins.push_back(bin);
        }
        encrypted_bins.resize(bin_size*bin_size);
        nrotated.push_back(encrypted_bins);
        i++;
    }
}



std::cout << "\t Size of nrots: " << nrotated.size() << " x " << nrotated[0].size() << std::endl;
    
    //extract rotated copies of each ciphertext
    vec2d_t<ciphertext_t> rotated;
std::cout << "Extract rotated: " << std::endl;
{
  Benchmark bench;
    for(int i = 0; i < (int)P; i++){
        vec1d_t<ciphertext_t> rot;
        //rot.resize(B*sizeof(ciphertext_t));
        rot.resize(bin_size*bin_size);
        #pragma omp parallel for
        for(int j = 0; j < (int)B; j++) {
            vecToFullPacking(cc, packingindex, nrotated[i][j], keyPair, N);
            //GetRot(cc, keyPair, packingindex, rotationindex, partiesindex, nrotated[i][j], N, bin_size, depth, i);
            rot[j] = GetRot(cc, keyPair, packingindex, rotationindex, partiesindex, nrotated[i][j], N, bin_size, depth, i);
        }
        rotated.push_back(rot);
    }
}
    std::cout << "\t Size of rots: " << rotated.size() << " x " << rotated[0].size() << std::endl;

    vec1d_t<ciphertext_t> lists(P);
    vec1d_t<plaintext_t> lists_plaintext(P);
    vec2d_t<ciphertext_t> pre_lists(P);
    vec2d_t<plaintext_t> pre_lists_plaintext(P);
    vec1d_t<ciphertext_t> psu(P);
    vec1d_t<plaintext_t> psu_plaintext(P);
    //psu[0] = parties_ciphertexts[0];
    psu_plaintext[0] = parties_plaintexts[0];


{
  Benchmark bench;

    for (int i = 1; i < (int)P; i++) {
        
        #pragma omp parallel for
        for(int j = 0; j < (int)B; j++) {

            //std::cout << "nrotated level before: " << nrotated[i][j]->GetLevel() << std::endl;
            
            nrotated[i][j] = cc->EvalMult(nrotated[i][j], nrot_masks[i - 1]);

            //std::cout << "nrotated level after: " << nrotated[i][j]->GetLevel() << std::endl;

            num_multiplications_server += 1;
            //std::cout << "Server: " << num_multiplications_server << std::endl;

            nrotated_plaintext[i][j] = multiplyPlaintexts(nrotated_plaintext[i][j], nrot_masks[i-1], cc, depth, N);
            auto bound =  i*bin_size*bin_size;

            auto rotated_sum = rotated[0][j]->Clone();
            auto rotated_sum_plaintext = rotated_plaintext[0][j];

           // std::cout << "rotated_sum j: " << "z : " << 0 << std::endl;
            //compute_diff_ptxt_single(rotated[0][j], rotated_plaintext[0][j], 1, cc, keyPair); //REMOVE WHEN BENCHMARKING
            //print_ciphertext<double>(rotated[0][j], cc, keyPair, N/2);
            //std::cout << "ROTATED PLAINTEXT TOP" << rotated_plaintext[0][j] << std::endl;


            for(int z = 1; z <= i-1; z++){
                cc->EvalAddInPlace(rotated_sum,rotated[z][j]);
                rotated_sum_plaintext = addPlaintexts(rotated_sum_plaintext, rotated_plaintext[z][j], cc, depth, N);
                std::cout << "rotated_sum: " << "z : " << z << std::endl;
                compute_diff_ptxt_single(rotated_sum, rotated_sum_plaintext, 1, cc, keyPair); //REMOVE WHEN BENCHMARKING
            }

            auto x_plaintext = comparePlaintexts(rotated_sum_plaintext, nrotated_plaintext[i][j], cc, depth, N, bound);

            if(i > 1) {
                 auto rot_after = cc->EvalBootstrap(rotated_sum, 2, 13);
                 rotated_sum = rot_after->Clone();
            }

            ciphertext_t comp;
            std::cout << "- Homomorphic Comparison " << " P: " << i+1 << ", B: " << j+1 << " ";
            {
                Benchmark bench;
                comp = HomomorphicComparison(rotated_sum, nrotated[i][j], dg, df, cc, keyPair, N, i, j, n, B);
                //print_ciphertext<double>(rotated_sum, cc, keyPair, n);
                //print_ciphertext<double>(nrotated[i][j], cc, keyPair, n);
            }
            
            
            if(1) {
                std::cout << "x_" << i+1 << "_" << j+1 << std::endl;
                compute_diff_ptxt_single(comp, x_plaintext, 1, cc, keyPair); //REMOVE WHEN BENCHMARKING
            }
            
        
            int t = bin_size;
            auto comp_rot = cc->EvalRotate(comp, t);

            //std::cout << "comp_rot depth: " << comp_rot->GetLevel() << std::endl;

            auto res = cc->EvalAdd(comp_rot, comp);

            for(size_t mu = 2; mu < (size_t)t; mu++) {
                //rotPrecomp = cc->EvalFastRotationPrecompute(comp_rot);
                comp_rot = cc->EvalRotate(comp_rot, t);
                cc->EvalAddInPlace(res, comp_rot);
            }

            res = cc->EvalRotate(res, -j*bin_size);

            auto x_vec = x_plaintext->GetCKKSPackedValue();
            std::vector<double> z(N/2, 0.0);
            for (size_t k = 0; k < N/2; k++) {
                for (size_t l = 0; l < (N/2)/t; l++) {
                    z[k] = z[k] + x_vec[(k+l*t)%(N/2)].real();
                }
            }
            Plaintext y_plaintext = cc->MakeCKKSPackedPlaintext(z, depth, 0, nullptr, N/2);
            
            auto r = GetListPSU(res, bin_size, j, cc, keyPair); 

            auto r_plaintext = GetListPSU_plaintext(y_plaintext, bin_size, j, cc, depth, N, bound);

            pre_lists[i].push_back(r);
            pre_lists_plaintext[i].push_back(r_plaintext);

        }

        lists[i] = pre_lists[i][0];
        lists_plaintext[i] = pre_lists_plaintext[i][0];

        for (int k = 1; k < (int)B; k++) {
            lists_plaintext[i] = addPlaintexts(lists_plaintext[i],pre_lists_plaintext[i][k], cc, depth, N);
            cc->EvalAddInPlace(lists[i], pre_lists[i][k]);
        }
        
        psu[i] = cc->EvalMult(parties_ciphertexts2[i], lists[i]); //here we remove the ID's that are duplicated among the parties
        psu_plaintext[i] = multiplyPlaintexts(parties_plaintexts[i], lists_plaintext[i], cc, depth, N);
    }
    std::cout << "PSU computation ";
}

    psu[0] = parties_ciphertexts2[0];

    vec2d_t<double> decrypted_psu;

    plaintext_t output2;
    
    for(int i = 0; i < (int)P; i++){
        decrypted_psu.push_back(print_ciphertext_decode<double>(psu[i], cc, keyPair, bin_size*bin_size));
    }

    std::cout << "Necessary depth: " << psu[P-1]->GetLevel() << " (current: " << depth << ")" << std::endl;

    //std::cout << "True PSU: " << std::endl;
    auto psu_easy = computePSU(parties);
    std::sort(psu_easy.begin(), psu_easy.end());
    removeZeros(psu_easy);

    std::vector<double> mergedPSU;

    for (const auto& innerVector : decrypted_psu) {
        mergedPSU.insert(mergedPSU.end(), innerVector.begin(), innerVector.end());
    }

    std::sort(mergedPSU.begin(), mergedPSU.end());
    removeZeros(mergedPSU);

    std::cout << "Length psu easy: " << psu_easy.size() << std::endl;
    std::cout << "Length psu encrypted: " << mergedPSU.size() << std::endl;

    //print_vec1d(psu_easy);
    //print_vec1d(mergedPSU);

    areVectorsEqual(psu_easy,mergedPSU);

    ////////////////////
    // COMPUTE OVER THE UNION
    ////////////////////

    vec2d_t<double> ages_aux(P);

    double sum_bmi = 0;
    double sum_card = 0;

    for (size_t i = 0; i < P; i++) {
        for (size_t j = 0; j < B*bin_size; j++) {
            if (fabs(ages[i][j]) < 1e-8) {
                ages_aux[i].push_back(0);
            }
            else {
                ages_aux[i].push_back(pow(2, 1.0/ages[i][j]));
            }
         }
    }

    for (size_t i = 0; i < P; i++) {
        
        vec1d_t<double> ag = ages_aux[i];
        vec1d_t<double> bm = bmi[i];

        vec1d_t<double> signs(ag.size(), 0);

        for(size_t j = 0; j < ag.size(); j++) {
            if(ag[j] - 59.5 > 0){
                signs[j] = 1;
            }
            if(ag[j] - 59.5 < 0){
                signs[j] = -1;
            }
        }

        for(size_t j = 0; j < ag.size(); j++) {
            signs[j] += 1;
            signs[j] /= 2;
            
            if(i > 0) {
                signs[j] *= extract_vector(lists_plaintext[i])[j];
            }

            signs[j] *= extract_vector(ones_parties_plaintext[i])[j];

            bm[j] *= signs[j];
            sum_bmi += bm[j];
            sum_card += signs[j];
        }
        //print_vec1d(signs);

    }


    //now we filter the people aged over age 60
    vec1d_t<double> lower_bound_vec(B*bin_size, 1.0/log2(59.5));
    plaintext_t lower_bound_plaintext = cc->MakeCKKSPackedPlaintext(lower_bound_vec, depth, 0, nullptr, N/2);

    vec1d_t<double> ones(B*bin_size, 1);
    plaintext_t ones_plaintext = cc->MakeCKKSPackedPlaintext(ones, depth, 0, nullptr, N/2);

    vec1d_t<double> one_half(B*bin_size, 0.5);
    plaintext_t one_half_plaintext = cc->MakeCKKSPackedPlaintext(one_half, depth, 0, nullptr, N/2);

    vec1d_t<double> ini(B*bin_size, 0);
    plaintext_t bmi_sum_plaintext = cc->MakeCKKSPackedPlaintext(ini, depth, 0, nullptr, N/2);
    ciphertext_t bmi_sum_ciphertext = cc->Encrypt(keyPair.publicKey, bmi_sum_plaintext);

    vec1d_t<ciphertext_t> filtered_age;
    vec1d_t<plaintext_t> filtered_age_p;

std::cout << "Compute g: " << std::endl;
{
  Benchmark bench;

    for (size_t i = 0; i < P; i++) {
        plaintext_t older_60_ptxt = subPlaintexts(lower_bound_plaintext, ages_plaintext[i], cc, depth, N);
        ciphertext_t older_60_ctxt = cc->EvalSub(lower_bound_plaintext, ages_ciphertext[i]);

        older_60_ptxt = signPlaintext(older_60_ptxt, cc, depth, N);
        older_60_ctxt = HomomorphicSign(older_60_ctxt, dg, df, cc, keyPair);
        std::cout << "Precision for older_60_ptxt and ctxt";
        compute_diff_ptxt_single(older_60_ctxt, older_60_ptxt, 1, cc, keyPair);

        older_60_ptxt = addPlaintexts(older_60_ptxt, ones_plaintext, cc, depth, N);
        cc->EvalAddInPlace(older_60_ctxt, ones_plaintext);

        auto filtered_age_plaintext = multiplyPlaintexts(older_60_ptxt, one_half_plaintext, cc, depth, N);
        auto filtered_age_ciphertext = cc->EvalMult(older_60_ctxt, one_half_plaintext);

        filtered_age.push_back(filtered_age_ciphertext);
        filtered_age_p.push_back(filtered_age_plaintext);


        compute_diff_ptxt_single(filtered_age_ciphertext, filtered_age_plaintext, 1, cc, keyPair);

        plaintext_t filtered_bmi_plaintext = multiplyPlaintexts(bmi_plaintext[i], filtered_age_plaintext, cc, depth, N);
        ciphertext_t filtered_bmi_ciphertext = cc->EvalMult(bmi_ciphertext[i], filtered_age_ciphertext);

        bmi_sum_plaintext = addPlaintexts(bmi_sum_plaintext, filtered_bmi_plaintext, cc, depth, N);
        cc->EvalAddInPlace(bmi_sum_ciphertext, filtered_bmi_ciphertext);

    }
    
    plaintext_t cardinality_sum_plaintext = ones_parties_plaintext[0];
    ciphertext_t cardinality_sum_ciphertext = cc->EvalMult(ones_parties_ciphertext[0], filtered_age[0]);
    
    std::cout << "filtered_ages: " << std::endl;
    for (size_t i = 1; i < P; i++) {
        auto a = cc->EvalMult(ones_parties_ciphertext[i], filtered_age[i]);        
        cc->EvalAddInPlace(cardinality_sum_ciphertext, a);

        auto a_p = multiplyPlaintexts(ones_parties_plaintext[i], filtered_age_p[i], cc, depth, N);
        cardinality_sum_plaintext = addPlaintexts(cardinality_sum_plaintext, a_p, cc, depth, N);
    }


    double numerator = add_elements_ciphertext(bmi_sum_ciphertext, cc, keyPair, B*bin_size);
    double denominator = add_elements_ciphertext(cardinality_sum_ciphertext, cc, keyPair, B*bin_size);

    std::cout << "Computation encrypted: " << numerator/denominator << std::endl;
    std::cout << "Computation clear: " << sum_bmi/sum_card << std::endl;

    std::cout << "ctxt numerator then denominator = " << numerator << ' ' << denominator << std::endl;
    std::cout << "ptxt numerator then denominator = " << sum_bmi << ' ' << sum_card << std::endl;
}
}

    