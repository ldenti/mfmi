#ifndef PANPP_FMD_H
#define PANPP_FMD_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stack>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "bits/deltavector.h"
#include "bits/nibblevector.h"
#include "bits/rlevector.h"
#include "bits/succinctvector.h"
#include "rlcsa.h"

#include "utils.h"

#ifdef USE_NIBBLE_VECTORS
// Use Nibble Vectors to encode our range endpoint bitmaps
typedef CSA::NibbleVector RangeVector;
typedef CSA::NibbleEncoder RangeEncoder;
#else
// Use RLEVectors to encode our range endpoint bitmaps
typedef CSA::RLEVector RangeVector;
typedef CSA::RLEEncoder RangeEncoder;
#endif

typedef std::pair<CSA::usint, CSA::usint> pair_type;

static const CSA::usint NUM_BASES = 6;

/**
 * Represents the state (or result) of an FMD-index search, which is two ranges
 * (one for the forward sequence, and one for the reverse complement) of equal
 * length. The ranges are stored as two start indices and a length. They can be
 * in either SA space (not counting the text start symbols at the beginning of
 * the BWT) or in BWT space.
 *
 * Range semantics are inclusive, so a length = 0 range holds 1 thing and its
 * reverse complement.
 */
struct FMDPosition {
  CSA::usint forward_start;
  CSA::usint reverse_start;
  // Offset 0 = only the entry at start/end. -1 = empty.
  CSA::sint end_offset;
  FMDPosition();
  FMDPosition(CSA::usint forward_start, CSA::usint reverse_start,
              CSA::usint end_offset);
  /**
   * Flip the FMDPosition around so the reverse complement interval is the
   * forward interval and visa versa.
   */
  FMDPosition flip() const;

  /**
   * Are two FMDPositions equal?
   */
  bool operator==(const FMDPosition &other) const;

  /**
   * Is an FMDPosition empty?
   */
  inline bool isEmpty() const { return end_offset < 0; }

  /**
   * Return the actual number of matches represented by an FMDPosition.
   */
  inline CSA::usint getLength() const { return end_offset + 1; }
};

/**
 * Provide pretty-printing for FMDPositions. See
 * <http://www.parashift.com/c++-faq/output-operator.html>
 */
std::ostream &operator<<(std::ostream &o, FMDPosition const &position);

const FMDPosition EMPTY_FMD_POSITION = FMDPosition(0, 0, -1);

/**
 * Defines an RLCSA index derivative that represents an FMD-index: an index of
 * DNA sequences (over the alphabet {A, C, G, T, N}) where all texts are present
 * with their reverse complements.
 *
 * In such an index, an ongoing search can be extended or retracted at either
 * end in O(1) time.
 *
 * See the paper "Exploring single-sample SNP and INDEL calling with whole-
 * genome de novo assembly" (2012), by Heng Li, which defines the FMD-index.
 */
class FMD : public CSA::RLCSA {

public:
  // We can only be constructed on a previously generated RLCSA index that
  // just happens to meet our requirements.
  explicit FMD(const std::string &base_name, bool print = false);

  /**
   * Extend a search by a character, either backward or forward. Ranges are in
   * BWT coordinates.
   */
  FMDPosition extend(FMDPosition range, CSA::usint c, bool backward) const;

  /**
   * Count occurrences of a pattern using the FMD search algorithm, iterating
   * through the pattern either forward or backward.
   */
  FMDPosition fmdCount(const std::string &pattern, bool backward = true) const;

  /**
   * Get an FMDPosition for the part of the BWT for things starting with the
   * given character.
   */
  FMDPosition getCharPosition(CSA::usint c) const;

protected:
  // How many times did we extend a search while mapping?
  static CSA::usint extends;

  // How many times did we restart a search while mapping?
  static CSA::usint restarts;

private:
  /**
   * Get an FMDPosition covering the whole SA.
   */
  FMDPosition getSAPosition() const;

  /**
   * Convert an FMDPosition in BWT coordinates to one in SA coordinates, in
   * place.
   */
  void convertToSAPosition(FMDPosition &bwt_position) const;

  // These are not allowed.
  FMD();
  FMD(const FMD &);
  FMD &operator=(const FMD &);
};

#endif
