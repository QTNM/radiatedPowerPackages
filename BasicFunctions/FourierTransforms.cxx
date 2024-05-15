/*
  FourierTransforms.cxx
*/

#include "BasicFunctions/FourierTransforms.h"

namespace rad {
FFTWComplex *doFFT(int length, double *theInput) {
  const int numFreqs = (length / 2) + 1;
  fftw_complex *out;
  out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * numFreqs);
  fftw_plan p = fftw_plan_dft_r2c_1d(length, theInput, out, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);

  rad::FFTWComplex *result = new rad::FFTWComplex[numFreqs];
  for (int i = 0; i < numFreqs; i++) {
    result[i] = rad::FFTWComplex(out[i][0], out[i][1]);
  }
  fftw_free(out);
  return result;
}

double *doInverseFFT(int length, const rad::FFTWComplex *theInput) {
  const int numFreqs = (length / 2) + 1;
  fftw_complex *in;
  in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * numFreqs);
  for (int i = 0; i < numFreqs; i++) {
    in[i][0] = theInput[i].re;
    in[i][1] = theInput[i].im;
  }
  double *result = new double[length];
  fftw_plan p = fftw_plan_dft_c2r_1d(length, in, result, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_free(in);
  return result;
}
}  // namespace rad