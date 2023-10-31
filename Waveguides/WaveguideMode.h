/*
  WaveguideMode.h

  Class representing modes in waveguides
*/

#include <ostream>
#include <string>

namespace rad {
enum ModeType { kTE, kTM, kTEM };

class WaveguideMode {
 private:
  // Meanings of i and j will vary depending on waveguide type
  unsigned int i;
  unsigned int j;
  ModeType type;

 public:
  /// @brief Default constructor
  WaveguideMode() : i(0), j(0), type(kTEM) {}

  /// @brief Parametrised constructor
  /// @param ind1
  /// @param ind2
  /// @param mType Mode type, kTE, kTM or kTEM
  WaveguideMode(unsigned int ind1, unsigned int ind2, ModeType mType)
      : i(ind1), j(ind2), type(mType) {}

  /// @brief Destructor
  ~WaveguideMode() {}

  // Getters for class members
  /// @brief Getter for mode type
  /// @return Mode type enum
  ModeType GetModeType() { return type; }

  /// @brief Get first mode index
  /// @return Mode index 1
  unsigned int GetModeIndex1() { return i; }

  /// @brief Get second mode index
  /// @return Mode index 2
  unsigned int GetModeIndex2() { return j; }

  /// @brief Overload output function
  /// @param os
  /// @param wm
  friend std::ostream& operator<<(std::ostream& os, WaveguideMode const& wm) {
    std::string msg;
    if (wm.type == kTE) {
      msg = "TE";
    } else if (wm.type == kTM) {
      msg = "TM";
    } else {
      msg = "TEM";
    }
    msg = msg + std::to_string(wm.i) + std::to_string(wm.j);
    return os << msg;
  }
};
}  // namespace rad