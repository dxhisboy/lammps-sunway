/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ATOM_VEC_KOKKOS_H
#define LMP_ATOM_VEC_KOKKOS_H

#include "atom_vec.h"
#include "kokkos_type.h"
#include <type_traits>

namespace LAMMPS_NS {

union d_ubuf {
  double d;
  int64_t i;
  KOKKOS_INLINE_FUNCTION
  d_ubuf(double arg) : d(arg) {}
  KOKKOS_INLINE_FUNCTION
  d_ubuf(int64_t arg) : i(arg) {}
  KOKKOS_INLINE_FUNCTION
  d_ubuf(int arg) : i(arg) {}
};

class AtomVecKokkos : public AtomVec {
 public:
  AtomVecKokkos(class LAMMPS *);
  virtual ~AtomVecKokkos() {}

  virtual void sync(ExecutionSpace space, unsigned int mask) = 0;
  virtual void modified(ExecutionSpace space, unsigned int mask) = 0;
  virtual void sync_overlapping_device(ExecutionSpace space, unsigned int mask) {};

  virtual int
    pack_comm_self(const int &n, const DAT::tdual_int_2d &list,
                   const int & iswap, const int nfirst,
                   const int &pbc_flag, const int pbc[]) = 0;
  //{return 0;}
  virtual int
    pack_comm_kokkos(const int &n, const DAT::tdual_int_2d &list,
                     const int & iswap, const DAT::tdual_xfloat_2d &buf,
                     const int &pbc_flag, const int pbc[]) = 0;
  //{return 0;}
  virtual void
    unpack_comm_kokkos(const int &n, const int &nfirst,
                       const DAT::tdual_xfloat_2d &buf) = 0;
  virtual int
    pack_border_kokkos(int n, DAT::tdual_int_2d k_sendlist,
                       DAT::tdual_xfloat_2d buf,int iswap,
                       int pbc_flag, int *pbc, ExecutionSpace space) = 0;
  //{return 0;};
  virtual void
    unpack_border_kokkos(const int &n, const int &nfirst,
                         const DAT::tdual_xfloat_2d &buf,
                         ExecutionSpace space) = 0;

  virtual int
    pack_exchange_kokkos(const int &nsend, DAT::tdual_xfloat_2d &buf,
                         DAT::tdual_int_1d k_sendlist,
                         DAT::tdual_int_1d k_copylist,
                         ExecutionSpace space, int dim, X_FLOAT lo, X_FLOAT hi) = 0;
  //{return 0;};
  virtual int
    unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf, int nrecv,
                           int nlocal, int dim, X_FLOAT lo, X_FLOAT hi,
                           ExecutionSpace space) = 0;
  //{return 0;};

 protected:

  class CommKokkos *commKK;
  size_t buffer_size;
  void* buffer;

  #ifdef KOKKOS_HAVE_CUDA
  template<class ViewType>
  Kokkos::View<typename ViewType::data_type,
               typename ViewType::array_layout,
               Kokkos::CudaHostPinnedSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged> >
  create_async_copy(const ViewType& src) {
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same<typename ViewType::execution_space,LMPDeviceType>::value,
                   Kokkos::CudaHostPinnedSpace,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer = Kokkos::kokkos_malloc<Kokkos::CudaHostPinnedSpace>(src.capacity());
       buffer_size = src.capacity();
    } else if (buffer_size < src.capacity()) {
       buffer = Kokkos::kokkos_realloc<Kokkos::CudaHostPinnedSpace>(buffer,src.capacity());
       buffer_size = src.capacity();
    }
    return mirror_type( buffer ,
                             src.dimension_0() ,
                             src.dimension_1() ,
                             src.dimension_2() ,
                             src.dimension_3() ,
                             src.dimension_4() ,
                             src.dimension_5() ,
                             src.dimension_6() ,
                             src.dimension_7() );
  }

  template<class ViewType>
  void perform_async_copy(const ViewType& src, unsigned int space) {
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same<typename ViewType::execution_space,LMPDeviceType>::value,
                   Kokkos::CudaHostPinnedSpace,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer = Kokkos::kokkos_malloc<Kokkos::CudaHostPinnedSpace>(src.capacity()*sizeof(typename ViewType::value_type));
       buffer_size = src.capacity();
    } else if (buffer_size < src.capacity()) {
       buffer = Kokkos::kokkos_realloc<Kokkos::CudaHostPinnedSpace>(buffer,src.capacity()*sizeof(typename ViewType::value_type));
       buffer_size = src.capacity();
    }
    mirror_type tmp_view( (typename ViewType::value_type*)buffer ,
                             src.dimension_0() ,
                             src.dimension_1() ,
                             src.dimension_2() ,
                             src.dimension_3() ,
                             src.dimension_4() ,
                             src.dimension_5() ,
                             src.dimension_6() ,
                             src.dimension_7() );
    if(space == Device) {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.h_view),
      Kokkos::deep_copy(LMPHostType(),src.d_view,tmp_view);
      src.modified_device() = src.modified_host();
    } else {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.d_view),
      Kokkos::deep_copy(LMPHostType(),src.h_view,tmp_view);
      src.modified_device() = src.modified_host();
    }
  }
  #else
  template<class ViewType>
  void perform_async_copy(ViewType& src, unsigned int space) {
    if(space == Device)
      src.template sync<LMPDeviceType>();
    else
      src.template sync<LMPHostType>();
  }
  #endif
};

}

#endif

/* ERROR/WARNING messages:

*/
