//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_H
#define QMCPLUSPLUS_DELAYED_UPDATE_H

#include "config.h"
#include <Numerics/OhmmsPETE/OhmmsVector.h>
#include <Numerics/OhmmsPETE/OhmmsMatrix.h>
#include "Numerics/OhmmsBlas.h"
#include "QMCWaveFunctions/DiracMatrix.h"
#include "Numerics/BlasThreadingEnv.h"
#include "RAJA/RAJA.hpp"
namespace qmcplusplus
{
/** implements delayed update on CPU using BLAS
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class DelayedUpdate
{
  /// define real type
  using real_type = typename scalar_traits<T>::real_type;
  /// orbital values of delayed electrons
  qmcplusplus::Matrix<T> U;
  RAJA::View<T,RAJA::Layout<2>> Uview;
  /// rows of Ainv corresponding to delayed electrons
  qmcplusplus::Matrix<T> V;
  /// Matrix inverse of B, at maximum KxK
  qmcplusplus::Matrix<T> Binv;
  /// scratch space, used during inverse update
  qmcplusplus::Matrix<T> tempMat;
  /// temporal scratch space used by SM-1
  Vector<T> temp;
  /// new column of B
  Vector<T> p;
  /// list of delayed electrons
  std::vector<int> delay_list;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count;
  /// matrix inversion engine
  DiracMatrix<T_FP, T> detEng;

public:
  /// default constructor
  DelayedUpdate() : delay_count(0) {}

  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, int delay)
  {
    std::array<long int, 2> rowMajor{{0,1}};
    V.resize(delay, norb);
    U.resize(delay, norb);

    p.resize(delay);
    temp.resize(norb);
    tempMat.resize(norb, delay);
    Binv.resize(delay, delay);
    delay_list.resize(delay);
  }

  /** compute the inverse of the transpose of matrix A
   * @param logdetT orbital value matrix
   * @param Ainv inverse matrix
   */
  inline void invert_transpose(const qmcplusplus::Matrix<T>& logdetT, qmcplusplus::Matrix<T>& Ainv, real_type& LogValue, real_type& PhaseValue)
  {
    detEng.invert_transpose(logdetT, Ainv, LogValue, PhaseValue);
    // safe mechanism
    delay_count = 0;
  }

  /** initialize internal objects when Ainv is refreshed
   * @param Ainv inverse matrix
   */
  inline void initializeInv(const qmcplusplus::Matrix<T>& Ainv)
  {
    // safe mechanism
    delay_count = 0;
  }

  /** compute the row of up-to-date Ainv
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   */
  template<typename VVT>
  inline void getInvRow(const qmcplusplus::Matrix<T>& Ainv, int rowchanged, VVT& invRow)
  {
    if (delay_count == 0)
    {
      // Ainv is fresh, directly access Ainv
      std::copy_n(Ainv[rowchanged], invRow.size(), invRow.data());
      return;
    }
    const T cone(1);
    const T czero(0);
    const int norb     = Ainv.rows();
    const int lda_Binv = Binv.cols();
    // save Ainv[rowchanged] to invRow
    std::copy_n(Ainv[rowchanged], norb, invRow.data());
    // multiply V (NxK) Binv(KxK) U(KxN) invRow right to the left
    BLAS::gemv('T', norb, delay_count, cone, U.data(), norb, invRow.data(), 1, czero, p.data(), 1);
    BLAS::gemv('N', delay_count, delay_count, cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
    BLAS::gemv('N', norb, delay_count, -cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
  }

  /** accept a move with the update delayed
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   * @param psiV new orbital values
   *
   * Before delay_count reaches the maximum delay, only Binv is updated with a recursive algorithm
   */
  template<typename VVT>
  inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
  {
    const T cminusone(-1);
    const T czero(0);
    const int norb     = Ainv.rows();
    const int lda_Binv = Binv.cols();
    std::copy_n(Ainv[rowchanged], norb, V[delay_count]);
    std::copy_n(psiV.data(), norb, U[delay_count]);
    delay_list[delay_count] = rowchanged;
    // the new Binv is [[X Y] [Z x]]
    BLAS::gemv('T', norb, delay_count + 1, cminusone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
    // x
    T y = -p[delay_count];
    for (int i = 0; i < delay_count; i++)
      y += Binv[delay_count][i] * p[i];
    Binv[delay_count][delay_count] = y = T(1) / y;
    // Y
    BLAS::gemv('T', delay_count, delay_count, y, Binv.data(), lda_Binv, p.data(), 1, czero, Binv.data() + delay_count,
               lda_Binv);
    // X
    BLAS::ger(delay_count, delay_count, cminusone, Binv[delay_count], 1, Binv.data() + delay_count, lda_Binv,
              Binv.data(), lda_Binv);
    // Z
    for (int i = 0; i < delay_count; i++)
      Binv[delay_count][i] *= -y;
    delay_count++;
    // update Ainv when maximal delay is reached
    if (delay_count == lda_Binv)
      updateInvMat(Ainv);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  inline void updateInvMat(Matrix<T>& Ainv)
  {
    if (delay_count == 0)
      return;
    // update the inverse matrix
    const T cone(1);
    const T czero(0);
    const int norb = Ainv.rows();
    if (delay_count == 1)
    {
      // this is a special case invoking the Fahy's variant of Sherman-Morrison update.
      // Only use the first norb elements of tempMat as a temporal array
      BLAS::gemv('T', norb, norb, cone, Ainv.data(), norb, U[0], 1, czero, temp.data(), 1);
      temp[delay_list[0]] -= cone;
      BLAS::ger(norb, norb, -Binv[0][0], V[0], 1, temp.data(), 1, Ainv.data(), norb);
    }
    else
    {

#define USING_RAJA
#ifdef USING_RAJA
      const int lda_Binv = Binv.cols();
      using View2 = RAJA::View<T, RAJA::Layout<2>>;
      View2 Ainv_(Ainv.data(), Ainv.rows(), Ainv.cols());
      View2 U_(U.data(), U.rows(), U.cols());
      View2 V_(V.data(), V.rows(), V.cols());
      View2 Binv_(Binv.data(), Binv.rows(), Binv.cols());
      View2 tempMat_(tempMat.data(), tempMat.rows(), tempMat.cols());
      using namespace RAJA;
      using POL=KernelPolicy<
        statement::Tile<0, tile_fixed<64>, RAJA::seq_exec,
        statement::Tile<1, tile_fixed<64>, RAJA::seq_exec,
        statement::Tile<2, tile_fixed<64>, RAJA::seq_exec,
        statement::For<0,omp_parallel_for_exec,
          statement::For<1, seq_exec,
            statement::Lambda<0,Segs<0,1>>,
              statement::For<2, seq_exec,
                statement::Lambda<1,Segs<0,1,2>>
      >>>>>>>;
      
      auto dc_seg = RangeSegment(0,delay_count);
      auto norb_seg = RangeSegment(0,norb);

      auto knl1 = make_kernel<POL>(make_tuple(RangeSegment(0,delay_count), RangeSegment(0,norb), RangeSegment(0,norb)),
        [&](auto i, auto j) {tempMat_(j,i) = 0;},
        [&](auto i, auto j, auto k) {
          tempMat_(j,i) += U_(i,k) * Ainv_(j,k);
        });

      auto knl2 = make_forall<omp_parallel_for_exec>(dc_seg, [&](auto i) {
        tempMat_(delay_list[i], i) -= 1;
      });

      auto knl3 = make_kernel<POL>(make_tuple(norb_seg, dc_seg, dc_seg),
        [&](auto i, auto j) {U_(j,i) = 0;},
        [&](auto i, auto j, auto k) {
          U_(j,i) += V_(k,i) * Binv_(j,k);
        });

      auto knl4 = make_kernel<POL>(make_tuple(norb_seg, norb_seg, dc_seg),
        [&](auto i, auto j) {},
        [&](auto i, auto j, auto k) {
          Ainv_(j,i) -= U_(k,i) * tempMat_(j,k);
        });

#ifdef USING_FORMAT_DECISIONS1
      auto dec = format_decisions(tie(U_, Ainv_, tempMat_, V_), knl1, knl2, knl3, knl4);
      dec.set_format_for(Ainv_, {0,1}, knl1);
      dec.set_format_for(Ainv_, {1,0}, knl4);
      dec.set_format_for(U_, {0,1}, knl1);
      dec.set_format_for(U_, {1,0}, knl3);
      dec.set_format_for(U_, {1,0}, knl4);
      dec.set_format_for(tempMat_, {1,0}, knl1);
      dec.set_format_for(tempMat_, {0,1}, knl4);
      dec.set_format_for(V_, {1,0}, knl3);
      auto comp = dec.finalize();
      comp();
#elif defined(USING_FORMAT_DECISIONS2)
      auto dec = format_decisions(tie(U_, V_), knl1, knl2, knl3, knl4);
      dec.set_format_for(U_, {0,1}, knl1);
      dec.set_format_for(U_, {1,0}, knl4);
      dec.set_format_for(V_, {1,0}, knl3);
      auto comp = dec.finalize();
      comp();
#elif defined(USING_FORMAT_DECISIONS0)
      auto dec = format_decisions(tie(U_, V_, Ainv_, tempMat_), knl1, knl2, knl3, knl4);
      auto comp = dec.finalize();
      comp();
#elif defined(USING_FORMAT_DECISIONS3)
      auto dec = format_decisions(tie(U_), knl1, knl2, knl3, knl4);
      dec.set_format_for(U_, {0,1}, knl1);
      dec.set_format_for(U_, {1,0}, knl4);
      auto comp = dec.finalize();
      comp();
#elif defined(USING_FORMAT_DECISIONS4)
      auto dec = format_decisions(tie(V_), knl1, knl2, knl3, knl4);
      dec.set_format_for(V_, {1,0}, knl3);
      auto comp = dec.finalize();
      comp();
#elif defined(USING_FORMAT_DECISIONS5)
      auto dec = format_decisions(tie(V_), knl1, knl2, knl3, knl4);
      dec.set_format_for(V_, {1,0}, knl1);
      dec.set_format_for(V_, {1,0}, knl2);
      dec.set_format_for(V_, {1,0}, knl3);
      dec.set_format_for(V_, {1,0}, knl4);
      dec.set_output_format(V_, {1,0});
      auto comp = dec.finalize();
      comp();

#elif defined(USING_FORMAT_DECISIONS6)
      auto dec = format_decisions(tie(Ainv_), knl1, knl2, knl3, knl4);
      dec.set_format_for(Ainv_, {1,0}, knl1);
      dec.set_format_for(Ainv_, {1,0}, knl2);
      dec.set_format_for(Ainv_, {1,0}, knl3);
      dec.set_format_for(Ainv_, {1,0}, knl4);
      dec.set_output_format(Ainv_, {1,0});
      auto comp = dec.finalize();
      comp();

#else
std::cout << "default raja\n";
      knl1();
      knl2();
      knl3();
      knl4();
#endif

#else


std::cout << "nonraja execution.\n";

      const int lda_Binv     = Binv.cols();
      int num_threads_nested = getNextLevelNumThreads();
      // always use serial when norb is small or only one second level thread
      bool use_serial(norb <= 256 || num_threads_nested == 1);
      if (use_serial || BlasThreadingEnv::NestedThreadingSupported())
      {
        // threading depends on BLAS
        BlasThreadingEnv knob(use_serial ? 1 : num_threads_nested);
        BLAS::gemm('T', 'N', delay_count, norb, norb, cone, U.data(), norb, Ainv.data(), norb, czero, tempMat.data(),
                   lda_Binv);
        for (int i = 0; i < delay_count; i++)
          tempMat(delay_list[i], i) -= cone;
        BLAS::gemm('N', 'N', norb, delay_count, delay_count, cone, V.data(), norb, Binv.data(), lda_Binv, czero,
                   U.data(), norb);
        BLAS::gemm('N', 'N', norb, norb, delay_count, -cone, U.data(), norb, tempMat.data(), lda_Binv, cone,
                   Ainv.data(), norb);
      }
      else
      {
        // manually threaded version of the above GEMM calls
#pragma omp parallel
        {
          const int block_size = getAlignedSize<T>((norb + num_threads_nested - 1) / num_threads_nested);
          int num_block        = (norb + block_size - 1) / block_size;
#pragma omp for
          for (int ix = 0; ix < num_block; ix++)
          {
            int x_offset = ix * block_size;
            BLAS::gemm('T', 'N', delay_count, std::min(norb - x_offset, block_size), norb, cone, U.data(), norb,
                       Ainv[x_offset], norb, czero, tempMat[x_offset], lda_Binv);
          }
#pragma omp master
          for (int i = 0; i < delay_count; i++)
            tempMat(delay_list[i], i) -= cone;
#pragma omp for
          for (int iy = 0; iy < num_block; iy++)
          {
            int y_offset = iy * block_size;
            BLAS::gemm('N', 'N', std::min(norb - y_offset, block_size), delay_count, delay_count, cone,
                       V.data() + y_offset, norb, Binv.data(), lda_Binv, czero, U.data() + y_offset, norb);
          }
#pragma omp for collapse(2) nowait
          for (int iy = 0; iy < num_block; iy++)
            for (int ix = 0; ix < num_block; ix++)
            {
              int x_offset = ix * block_size;
              int y_offset = iy * block_size;
              BLAS::gemm('N', 'N', std::min(norb - y_offset, block_size), std::min(norb - x_offset, block_size),
                         delay_count, -cone, U.data() + y_offset, norb, tempMat[x_offset], lda_Binv, cone,
                         Ainv[x_offset] + y_offset, norb);
            }
        }
      }

#endif

    }
    delay_count = 0;
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DELAYED_UPDATE_H
