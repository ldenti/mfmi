#include "fmd_simple.hpp"

FMDPosition::FMDPosition(CSA::usint forward_start, CSA::usint reverse_start,
                         CSA::usint end_offset)
    : forward_start(forward_start), reverse_start(reverse_start),
      end_offset(end_offset) {}

FMDPosition::FMDPosition()
    : forward_start(0), reverse_start(0), end_offset(-1) {}

FMDPosition FMDPosition::flip() const {
  // Swap the two intervals of the bi-interval
  return FMDPosition(reverse_start, forward_start, end_offset);
}

bool FMDPosition::operator==(const FMDPosition &other) const {
  // Compare all the fields.
  return forward_start == other.forward_start &&
         reverse_start == other.reverse_start && end_offset == other.end_offset;
}

std::ostream &operator<<(std::ostream &o, FMDPosition const &position) {
  // Report both the ranges that we represent.
  return o << position.forward_start << "-"
           << (position.forward_start + position.end_offset) << "|"
           << position.reverse_start << "-"
           << (position.reverse_start + position.end_offset);
}

FMD::FMD(const std::string &base_name, bool print) : RLCSA(base_name, print) {}

FMDPosition FMD::extend(FMDPosition range, CSA::usint c, bool backward) const {

  // More or less directly implemented off of algorithms 2 and 3 in "Exploring
  // single-sample SNP and INDEL calling with whole-genome de novo assembly"
  // (Li, 2012). However, our character indices are one less, since we don't
  // allow search patterns to include the end-of-text symbol. We also use
  // alphabetical ordering instead of the paper's N-last ordering in the FM-
  // index, and consequently need to assign reverse ranges in alphabetical order
  // by reverse complement.

  if (backward) {

    // Only allow characters in the index
    if (c >= CSA::CHARS || this->array[c] == 0) {
      return EMPTY_FMD_POSITION;
    }
    // Only allow DNA bases
    if (!isBase(c)) {
      return EMPTY_FMD_POSITION;
    }

    // We have an array of FMDPositions, one per base, that we will fill in by a
    // tiny dynamic programming.
    FMDPosition answers[NUM_BASES];

    for (CSA::usint base = 0; base < NUM_BASES; base++) {
      // Go through the bases in arbitrary order.

      // Count up the number of characters < this base, including sequence stop
      // characters.
      CSA::usint start =
          this->alphabet->cumulative(base) + this->number_of_sequences - 1;

      // Get a pointer to the bit vector for this letter, which might be NULL if
      // this base never appeared.
      CSA::PsiVector *vector = this->array[base];

      if (vector == NULL) {

        // Fill in forward_start and length with the knowledge that this
        // character doesn't exist. forward_start should never get used, but
        // end_offset will get used and probably needs to be -1 for empty.
        answers[base].end_offset = -1;

      } else {
        // Get an iterator for the bit vector for this character, for
        // calculating ranks/occurrences.
        CSA::PsiVector::Iterator iter(*vector);

        // Fill in the forward-strand start positions and range end_offsets for
        // each base's answer. TODO: do we want at_least set or not? What does
        // it do?

        // First cache the forward_start rank we re-use
        CSA::usint forward_start_rank = iter.rank(range.forward_start, true);

        answers[base].forward_start = start + forward_start_rank;
        answers[base].end_offset =
            iter.rank(range.forward_start + range.end_offset, false) -
            forward_start_rank;
      }
    }

    // Since we don't keep an FMDPosition for the non-base end-of-text
    // character, we need to track its length separately in order for the DP
    // algorithm given in the paper to be implementable. We calculate
    // occurrences of the text end character (i.e. how much of the current range
    // is devoted to things where an endOfText comes next) implicitly: it's
    // whatever part of the length of the range is unaccounted-for by the other
    // characters. We need to use the length accessor because ranges with one
    // thing have the .end_offset set to 0.
    CSA::usint endOfTextLength = range.getLength();

    for (CSA::usint base = 0; base < NUM_BASES; base++) {
      // Go through the bases in order and account for their lengths.
      endOfTextLength -= answers[base].getLength();
    }

    // The endOfText character is the very first character we need to account
    // for when subdividing the reverse range and picking which subdivision to
    // take.

    // Next, allocate the range for the base that comes first in alphabetical
    // order by reverse complement.
    answers[0].reverse_start = range.reverse_start + endOfTextLength;

    for (CSA::usint base = 1; base < NUM_BASES; base++) {
      // For each subsequent base in alphabetical order by reverse complement
      // (as stored in BASES), allocate it the next part of the reverse range.

      answers[base].reverse_start = answers[fm6_comp(base)].reverse_start +
                                    answers[fm6_comp(base)].getLength();
    }

    // Now all the per-base answers are filled in.

    for (CSA::usint base = 0; base < NUM_BASES; ++base) {
      // For each base in arbitrary order
      if (base == c) {
        // This is the base we're actually supposed to be extending with. Return
        // its answer.
        return answers[base];
      }
    }

    // If we get here, they gave us something not in BASES somehow.
    throw "Unrecognized base";
  } else {
    // Flip the interval, do backwards search with the reverse complement of the
    // base, and then flip back.
    return this->extend(range.flip(), fm6_comp(c), true).flip();
  }
}

FMDPosition FMD::fmdCount(const std::string &pattern, bool backward) const {
  if (pattern.length() == 0) {
    return this->getSAPosition();
  }

  // Keep an FMDPosition to store our intermediate result in.
  FMDPosition index_position;

  if (backward) {
    // Start at the end of the pattern and work towards the front

    std::string::const_reverse_iterator iter = pattern.rbegin();
    index_position = this->getCharPosition((CSA::uchar)*iter);
    if (index_position.isEmpty()) {
      return index_position;
    }

    for (++iter; iter != pattern.rend(); ++iter) {
      // Backwards extend with subsequent characters.
      index_position = this->extend(index_position, *iter, true);
      if (index_position.isEmpty()) {
        return EMPTY_FMD_POSITION;
      }
    }
  } else {
    // Start at the front of the pattern and work towards the end.

    std::string::const_iterator iter = pattern.begin();
    index_position = this->getCharPosition((CSA::uchar)*iter);
    if (index_position.isEmpty()) {
      return index_position;
    }

    for (++iter; iter != pattern.end(); ++iter) {
      // Forwards extend with subsequent characters.
      index_position = this->extend(index_position, *iter, false);
      if (index_position.isEmpty()) {
        return EMPTY_FMD_POSITION;
      }
    }
  }

  this->convertToSAPosition(index_position);

  return index_position;
}

CSA::usint FMD::extends = 0;
CSA::usint FMD::restarts = 0;

FMDPosition FMD::getCharPosition(CSA::usint c) const {
  if (c >= CSA::CHARS || this->array[c] == 0) {
    return EMPTY_FMD_POSITION;
  }

  if (!isBase(c)) {
    return EMPTY_FMD_POSITION;
  }

  pair_type forward_range = this->alphabet->getRange(c);
  this->convertToBWTRange(forward_range);

  pair_type reverse_range = this->alphabet->getRange(fm6_comp(c));
  this->convertToBWTRange(reverse_range);

  // Make the FMDPosition rolling together both ranges.
  // TODO: Make sure both ranges are the same length, as they should be.
  FMDPosition position(forward_range.first, reverse_range.first,
                       forward_range.second - forward_range.first);

  return position;
}

FMDPosition FMD::getSAPosition() const {
  return FMDPosition(0, 0, this->data_size - 1);
}

void FMD::convertToSAPosition(FMDPosition &bwt_position) const {
  bwt_position.forward_start -= this->number_of_sequences;
  bwt_position.reverse_start -= this->number_of_sequences;
}
