// Copyright (c) 2020 Mario Brcic, Robert Vaser

#ifndef SPOA_ARCHITECTURES_HPP_
#define SPOA_ARCHITECTURES_HPP_

namespace spoa {

enum class Architecture {
  kAVX2,
  kSSE4_1,
  kSSE2,
  kAutomatic
};

}  // namespace spoa

#endif  // SPOA_ARCHITECTURES_HPP_
