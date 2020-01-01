/*
(**************************************************************************)
(*                                                                        *)
(*                                Schifra                                 *)
(*                Reed-Solomon Error Correcting Code Library              *)
(*                                                                        *)
(* Release Version 0.0.1                                                  *)
(* http://www.schifra.com                                                 *)
(* Copyright (c) 2000-2020 Arash Partow, All Rights Reserved.             *)
(*                                                                        *)
(* The Schifra Reed-Solomon error correcting code library and all its     *)
(* components are supplied under the terms of the General Schifra License *)
(* agreement. The contents of the Schifra Reed-Solomon error correcting   *)
(* code library and all its components may not be copied or disclosed     *)
(* except in accordance with the terms of that agreement.                 *)
(*                                                                        *)
(* URL: http://www.schifra.com/license.html                               *)
(*                                                                        *)
(**************************************************************************)
*/


#ifndef INCLUDE_SCHIFRA_GALOIS_FIELD_ELEMENT_HPP
#define INCLUDE_SCHIFRA_GALOIS_FIELD_ELEMENT_HPP


#include <iostream>
#include <vector>

#include "schifra_galois_field.hpp"


namespace schifra
{

   namespace galois
   {

      class field_element
      {
      public:

         field_element(const field& gfield)
         : field_(gfield),
           poly_value_(-1)
         {}

         field_element(const field& gfield,const field_symbol& v)
         : field_(const_cast<field&>(gfield)),
           poly_value_(v)
         {}

         field_element(const field_element& gfe)
         : field_(const_cast<field&>(gfe.field_)),
           poly_value_(gfe.poly_value_)
         {}

        ~field_element()
         {}

         inline field_element& operator = (const field_element& gfe)
         {
            if ((this != &gfe) && (&field_ == &gfe.field_))
            {
               poly_value_ = gfe.poly_value_;
            }

            return *this;
         }

         inline field_element& operator = (const field_symbol& v)
         {
            poly_value_ = v & field_.size();
            return *this;
         }

         inline field_element& operator += (const field_element& gfe)
         {
            poly_value_ ^= gfe.poly_value_;
            return *this;
         }

         inline field_element& operator += (const field_symbol& v)
         {
            poly_value_ ^= v;
            return *this;
         }

         inline field_element& operator -= (const field_element& gfe)
         {
            *this += gfe;
            return *this;
         }

         inline field_element& operator -= (const field_symbol& v)
         {
            *this += v;
            return *this;
         }

         inline field_element& operator *= (const field_element& gfe)
         {
            poly_value_ = field_.mul(poly_value_, gfe.poly_value_);
            return *this;
         }

         inline field_element& operator *= (const field_symbol& v)
         {
            poly_value_ = field_.mul(poly_value_, v);
            return *this;
         }

         inline field_element& operator /= (const field_element& gfe)
         {
            poly_value_ = field_.div(poly_value_, gfe.poly_value_);
            return *this;
         }

         inline field_element& operator /= (const field_symbol& v)
         {
            poly_value_ = field_.div(poly_value_, v);
            return *this;
         }

         inline field_element& operator ^= (const int& n)
         {
            poly_value_ = field_.exp(poly_value_,n);
            return *this;
         }

         inline bool operator == (const field_element& gfe) const
         {
            return ((field_  == gfe.field_) && (poly_value_ == gfe.poly_value_));
         }

         inline bool operator == (const field_symbol& v) const
         {
            return (poly_value_ == v);
         }

         inline bool operator != (const field_element& gfe) const
         {
            return ((field_ != gfe.field_) || (poly_value_ != gfe.poly_value_));
         }

         inline bool operator != (const field_symbol& v) const
         {
            return (poly_value_ != v);
         }

         inline bool operator < (const field_element& gfe)
         {
            return (poly_value_ < gfe.poly_value_);
         }

         inline bool operator < (const field_symbol& v)
         {
            return (poly_value_ < v);
         }

         inline bool operator > (const field_element& gfe)
         {
            return (poly_value_ > gfe.poly_value_);
         }

         inline bool operator > (const field_symbol& v)
         {
            return (poly_value_ > v);
         }

         inline field_symbol index() const
         {
            return field_.index(poly_value_);
         }

         inline field_symbol poly() const
         {
            return poly_value_;
         }

         inline field_symbol& poly()
         {
            return poly_value_;
         }

         inline const field& galois_field() const
         {
            return field_;
         }

         inline field_symbol inverse() const
         {
            return field_.inverse(poly_value_);
         }

         inline void normalize()
         {
            poly_value_ &= field_.size();
         }

         friend std::ostream& operator << (std::ostream& os, const field_element& gfe);

      private:

         const field& field_;
         field_symbol poly_value_;

      };

      inline field_element operator + (const field_element& a, const field_element& b);
      inline field_element operator - (const field_element& a, const field_element& b);
      inline field_element operator * (const field_element& a, const field_element& b);
      inline field_element operator * (const field_element& a, const field_symbol&  b);
      inline field_element operator * (const field_symbol&  a, const field_element& b);
      inline field_element operator / (const field_element& a, const field_element& b);
      inline field_element operator ^ (const field_element& a, const int&           b);

      inline std::ostream& operator << (std::ostream& os, const field_element& gfe)
      {
         os << gfe.poly_value_;
         return os;
      }

      inline field_element operator + (const field_element& a, const field_element& b)
      {
         field_element result = a;
         result += b;
         return result;
      }

      inline field_element operator - (const field_element& a, const field_element& b)
      {
         field_element result = a;
         result -= b;
         return result;
      }

      inline field_element operator * (const field_element& a, const field_element& b)
      {
         field_element result = a;
         result *= b;
         return result;
      }

      inline field_element operator * (const field_element& a, const field_symbol& b)
      {
         field_element result = a;
         result *= b;
         return result;
      }

      inline field_element operator * (const field_symbol& a, const field_element& b)
      {
         field_element result = b;
         result *= a;
         return result;
      }

      inline field_element operator / (const field_element& a, const field_element& b)
      {
         field_element result = a;
         result /= b;
         return result;
      }

      inline field_element operator ^ (const field_element& a, const int& b)
      {
         field_element result = a;
         result ^= b;
         return result;
      }

   } // namespace galois

} // namespace schifra

#endif
