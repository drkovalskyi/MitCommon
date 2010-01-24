//--------------------------------------------------------------------------------------------------
// $Id: Angle.h,v 1.3 2009/07/20 03:12:22 loizides Exp $
//
// Angle classes, with automatic range checking.  SignedAngle ranges from -PI to PI, UnsignedAngle
// from 0 to 2PI.  The two types of angle convert implicitly into double.
//
// The Angle classes reimplement (inline) most arithmetic operators, instead of relying on the
// implicit conversion to double and converting the result back when it is assigned to an Angle
// variable. There are two reasons for this:
//
// 1) Temporary values are kept in range (example: std::cout << alpha+beta)
// 2) Optimization (the difference of two Angles is guaranteed to lie less
//    than 2PI outside the allowed range, for example)
//
// The only exception is multiplication: the "2" in "2 * alpha" could be interpreted as a length. If
// you want the result of a multiplication to be bounded, assign it to a temporary Angle, or use *=.
//
// Due to operator ambiguity, some mixed-type arithmetics might fail:
// for example, (int + Angle) or (Angle / long)
// The most common operations are supported:
// double +- angle; angle +-/ double; angle / int
//
// Comparison between Angles is done via conversion to double. Make sure you do not compare a
// SignedAngle to an UnsignedAngle.
//
// Since the bounds and return types for the two classes are different, operators cannot be factored
// in a single superclass.
//
// In C++, implicit conversion operators are shallow - they are not nested.  Therefore, if you want
// to convert a SignedAngle into an UnsignedAngle, you need an intermediate cast into double. For
// example:
//
//     SignedAngle s; UnsignedAngle u = double(s);
//
// Adding a direct conversion between the two types of angle would result in several ambiguities.
//
// Author: C.Paus (stolen from CDF implementation of Paolo Gatti, University of Padova / INFN,
//                 therefore not all our coding conventions fulfilled) 
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_MATHTOOLS_ANGLE_H
#define MITCOMMON_MATHTOOLS_ANGLE_H

#include <Riostream.h>
#include <Rtypes.h>
#include <math.h>
#include <limits>

namespace mithep
{
//--------------------------------------------------------------------------------------------------
//
// Signed angles, ranging from -PI to +PI
//
//--------------------------------------------------------------------------------------------------
  class SignedAngle
  {
    public:
      // Constructor
      SignedAngle(double angle = 0.0);

      // No destructor; use standard copy constructor and assignment operator

      // Implicit conversion to double
      operator             double() const;

      // Self-modifying operations
      SignedAngle&         operator += (const SignedAngle& other);
      SignedAngle&         operator -= (const SignedAngle& other);
      SignedAngle&         operator *= (double other);
      SignedAngle&         operator /= (double other);

      // Other operations
      friend SignedAngle   operator + (double me, const SignedAngle& other);
      SignedAngle          operator + (const SignedAngle& other) const;
      SignedAngle          operator + (double other) const;
      friend SignedAngle   operator - (double me, const SignedAngle& other);
      SignedAngle          operator - (const SignedAngle& other) const;
      SignedAngle          operator - (double other) const;
      SignedAngle          operator / (double other) const;
      SignedAngle          operator / (int other) const;

      // I/O
      //friend std::ostream& operator << (std::ostream& out, const SignedAngle& me)
      //  { out << me.fValue;  return out; }

      //friend std::istream& operator >> (std::istream& in, SignedAngle& me)
      //  { in >> me.fValue; me.FixRangeSlow(); return in; }

    private:
      // Enforce the correct range for the angle's value. The first version is faster, but assumes 
      // the current value is within 2PI of the correct range.  The second version is slower but 
      // always works
      void                 FixRangeFast();
      void                 FixRangeSlow();

      double               fValue;

    ClassDef(SignedAngle, 0) // Signed angle class
  };


//--------------------------------------------------------------------------------------------------
//
// Unsigned angles, ranging from 0 to 2PI
//
//--------------------------------------------------------------------------------------------------
  class UnsignedAngle
  {
    public:
      // Constructor
      UnsignedAngle(double angle = 0.0);

      // No destructor; use standard copy constructor and assignment operator

      // Implicit conversion to double
      operator             double() const;

      // Self-modifying operations
      UnsignedAngle&       operator += (const UnsignedAngle& other);
      UnsignedAngle&       operator -= (const UnsignedAngle& other);
      UnsignedAngle&       operator *= (double other);
      UnsignedAngle&       operator /= (double other);

      // Other operations
      friend UnsignedAngle operator + (double me, const UnsignedAngle& other);
      UnsignedAngle        operator + (const UnsignedAngle& other) const;
      UnsignedAngle        operator + (double other) const;
      friend UnsignedAngle operator - (double me, const UnsignedAngle& other);
      UnsignedAngle        operator - (const UnsignedAngle& other) const;
      UnsignedAngle        operator - (double other) const;
      UnsignedAngle        operator / (double other) const;
      UnsignedAngle        operator / (int other) const;

      // I/O
//       friend std::ostream& operator << (std::ostream& out, const UnsignedAngle& me)
//         { out << me.fValue; return out; }
//       friend std::istream& operator >> (std::istream& in, UnsignedAngle& me)
//         { in >> me.fValue; me.FixRangeSlow(); return in; }

    private:
      // Enforce the correct range for the angle's value. The first version is faster, 
      // but assumes the current value is within 2PI of the correct range. 
      // The second version is slower but always works.

      void                 FixRangeFast();
      void                 FixRangeSlow();

      double               fValue;

    ClassDef(UnsignedAngle, 0) // Unsigned angle class
  };

// By default, "Angles" are unsigned (CDF convention: 0-2PI range).
  typedef UnsignedAngle Angle;

//--------------------------------------------------------------------------------------------------
  inline
  void mithep::SignedAngle::FixRangeFast()
  {
    if      (fValue >= M_PI)
    fValue -= 2 * M_PI;
    else if (fValue < -M_PI)
      fValue += 2 * M_PI;
  }

//--------------------------------------------------------------------------------------------------
  inline
  void mithep::SignedAngle::FixRangeSlow()
  {
    if (fValue < -M_PI) {
      if (fValue < double(std::numeric_limits<int>::min()) * 2 * M_PI)
        std::cout << "@SUB=SignedAngle::FixRangeSlow"
                  << "Angle out of range:  " << fValue << std::endl;
      else {
        int shift = int((M_PI - fValue) / (2 * M_PI));
        fValue += shift * 2 * M_PI;
      }
    }
    // The previous step might have brought -3PI to +PI, so no 'else' here.
    if (fValue >= M_PI) {
      if (fValue > double(std::numeric_limits<int>::max()) * 2 * M_PI) {
        std::cout << "@SUB=SignedAngle::FixRangeSlow"
                  << "Angle out of range:  " << fValue << std::endl;
      }
      else {
        int shift = int((fValue + M_PI) / (2 * M_PI));
        fValue -= shift * 2 * M_PI;
      }
    }
  }


//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle::SignedAngle(double angle)
    : fValue(angle)
  {
    // Constructor

    FixRangeSlow();
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle::operator double() const
  {
    // Implicit conversion to double
    return fValue;
  }

// Self-modifying operations
//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle& mithep::SignedAngle::operator += (const mithep::SignedAngle& other)
  {
    fValue += other.fValue;
    FixRangeFast(); // other.fValue is guaranteed to be in range
    return *this;
  }

//--------------------------------------------------------------------------------------------------
inline
  mithep::SignedAngle& mithep::SignedAngle::operator -= (const mithep::SignedAngle& other)
  {
    fValue -= other.fValue;
    FixRangeFast(); // other.fValue is guaranteed to be in range
    return *this;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle& mithep::SignedAngle::operator *= (double other)
  {
    fValue *= other;
    FixRangeSlow();
    return *this;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle& mithep::SignedAngle::operator /= (double other)
  {
    fValue /= other;
    FixRangeSlow();
    return *this;
  }

// Other operations
//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle operator + (double me, const mithep::SignedAngle& other)
  {
    mithep::SignedAngle result(me);
    result += other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle mithep::SignedAngle::operator + (const mithep::SignedAngle& other) const
  {
    mithep::SignedAngle result(*this);
    result += other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle mithep::SignedAngle::operator + (double other) const
  {
    return mithep::SignedAngle(fValue + other);
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle operator - (double me, const mithep::SignedAngle& other)
  {
    mithep::SignedAngle result(me);
    result -= other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle mithep::SignedAngle::operator - (const mithep::SignedAngle& other) const
  {
    mithep::SignedAngle result(*this);
    result -= other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle mithep::SignedAngle::operator - (double other) const
  {
    return mithep::SignedAngle(fValue - other);
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle mithep::SignedAngle::operator / (double other) const
  {
    mithep::SignedAngle result(*this);
    result /= other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::SignedAngle mithep::SignedAngle::operator / (int other) const
  {
    mithep::SignedAngle result(*this);
    result /= other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  void mithep::UnsignedAngle::FixRangeFast()
  {
    if      (fValue >= 2 * M_PI)
      fValue -= 2 * M_PI;
    else if (fValue < 0)
      fValue += 2 * M_PI;
  }

//--------------------------------------------------------------------------------------------------
  inline
  void mithep::UnsignedAngle::FixRangeSlow()
  {
    if (fValue < 0) {
      if (fValue < double(std::numeric_limits<int>::min()) * 2 * M_PI)
        std::cout << "@SUB=UnsignedAngle::FixRangeSlow"
                  << "Angle out of range:  " << fValue << std::endl;
      else  {
        int shift = 1 - int(fValue / (2 * M_PI));
        fValue += shift * 2 * M_PI;
      }
    }
  // The previous step might have brought -2PI to +2PI, so no 'else' here.
    if (fValue >= 2 * M_PI) {
      if (fValue > double(std::numeric_limits<int>::max()) * 2 * M_PI)
        std::cout << "@SUB=UnsignedAngle::FixRangeSlow"
                  << "Angle out of range:  " << fValue << std::endl;
      else {
        int shift = int(fValue / (2 * M_PI));
        fValue -= shift * 2 * M_PI;
      }
    }
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle::UnsignedAngle(double angle)
    : fValue(angle)
  {
    // Constructor
    FixRangeSlow();
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle::operator double() const
  {
    // Implicit conversion to double  
    
    return fValue;
  }

// Self-modifying operations
//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle& mithep::UnsignedAngle::operator += (const mithep::UnsignedAngle& other)
  {
    fValue += other.fValue;
    FixRangeFast(); // other.fValue is guaranteed to be in range
    return *this;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle& mithep::UnsignedAngle::operator -= (const mithep::UnsignedAngle& other)
  {
    fValue -= other.fValue;
    FixRangeFast(); // other.fValue is guaranteed to be in range
    return *this;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle& mithep::UnsignedAngle::operator *= (double other)
  {
    fValue *= other;
    FixRangeSlow();
    return *this;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle& mithep::UnsignedAngle::operator /= (double other)
  {
    fValue /= other;
    FixRangeSlow();
    return *this;
  }

// Other operations
//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle operator + (double me, const mithep::UnsignedAngle& other)
  {
    mithep::UnsignedAngle result(me);
    result += other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle mithep::UnsignedAngle::operator + (const mithep::UnsignedAngle& other) const
  {
    mithep::UnsignedAngle result(*this);
    result += other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle mithep::UnsignedAngle::operator + (double other) const
  {
    return mithep::UnsignedAngle(fValue + other);
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle operator - (double me, const mithep::UnsignedAngle& other)
  {
    mithep::UnsignedAngle result(me);
    result -= other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle mithep::UnsignedAngle::operator - (const mithep::UnsignedAngle& other) const
  {
    mithep::UnsignedAngle result(*this);
    result -= other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle mithep::UnsignedAngle::operator - (double other) const
  {
    return mithep::UnsignedAngle(fValue - other);
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle mithep::UnsignedAngle::operator / (double other) const
  {
    mithep::UnsignedAngle result(*this);
    result /= other;
    return result;
  }

//--------------------------------------------------------------------------------------------------
  inline
  mithep::UnsignedAngle mithep::UnsignedAngle::operator / (int other) const
  {
    mithep::UnsignedAngle result(*this);
    result /= other;
    return result;
  }
} //close namespace here because of friend declarations
#endif
