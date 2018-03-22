#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>

namespace aa {

enum class AminoAcid : uint8_t {
  ALANINE = 0,     // A
  CYSTEINE,        // C
  ASPARTIC_ACID,   // D
  GLUTAMIC_ACID,   // E
  PHENYLALANINE,   // F
  GLYCINE,         // G
  HISTIDINE,       // H
  ISOLEUCINE,      // I
  LYSINE,          // K
  LEUCINE,         // L
  METHIONINE,      // M
  ASPARAGINE,      // N
  PROLINE,         // P
  GLUTAMINE,       // Q
  ARGININE,        // R
  SERINE,          // S
  THREONINE,       // T
  VALINE,          // V
  TRYPTOPHAN,      // W
  TYROSINE,        // Y
  STOP
};

// From https://stackoverflow.com/questions/8357240/how-to-automatically-convert-strongly-typed-enum-into-int
template <typename E>
constexpr typename std::underlying_type<E>::type to_underlying(E e) noexcept {
  return static_cast<typename std::underlying_type<E>::type>(e);
}

constexpr char one_letter_codes[] = "ACDEFGHIKLMNPQRSTVWY*";

constexpr char to_one_letter(AminoAcid aa) {
  return one_letter_codes[to_underlying(aa)];
}

constexpr uint8_t dignucl(char ch) {
  switch (ch) {
    default:
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
  }
}

constexpr size_t codon_to_idx(char codon0, char codon1, char codon2) {
  return
    dignucl(codon0) << (2*2) |
    dignucl(codon1) << (2*1) |
    dignucl(codon2) << (2*0);
}

constexpr size_t codon_to_idx(const char codon[3]) {
  return codon_to_idx(codon[0], codon[1], codon[2]);
}

constexpr uint8_t aa_table[64] = {
  0x08, 0x0B, 0x08, 0x0B, 0x10, 0x10, 0x10, 0x10, 0x0E, 0x0F, 0x0E, 0x0F, 0x07, 0x07, 0x0A, 0x07,
  0x0D, 0x06, 0x0D, 0x06, 0x0C, 0x0C, 0x0C, 0x0C, 0x0E, 0x0E, 0x0E, 0x0E, 0x09, 0x09, 0x09, 0x09,
  0x03, 0x02, 0x03, 0x02, 0x00, 0x00, 0x00, 0x00, 0x05, 0x05, 0x05, 0x05, 0x11, 0x11, 0x11, 0x11,
  0x14, 0x13, 0x14, 0x13, 0x0F, 0x0F, 0x0F, 0x0F, 0x14, 0x01, 0x12, 0x01, 0x09, 0x04, 0x09, 0x04
};

constexpr AminoAcid to_aa(const char codon[3]) {
  size_t idx = codon_to_idx(codon);
  return AminoAcid(aa_table[idx]);
}

constexpr AminoAcid to_aa(char codon0, char codon1, char codon2) {
  size_t idx = codon_to_idx(codon0, codon1, codon2);
  return AminoAcid(aa_table[idx]);
}

inline std::string translate(const char *nts) {
  size_t len_nts = strlen(nts);
  size_t len_aas = len_nts / 3;  // floor

  std::string aas(len_aas, ' ');
  for (size_t i = 0; i < len_aas; ++i) {
    aas[i] = to_one_letter(to_aa(nts + 3*i));
  }

  return aas;
}

inline std::string translate(const std::string &nts) {
  return translate(nts.c_str());
}

}  // namespace aa

using aa::to_aa;
using aa::translate;
using aa::to_one_letter;
