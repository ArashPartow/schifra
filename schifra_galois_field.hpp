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


#ifndef INCLUDE_SCHIFRA_GALOIS_FIELD_HPP
#define INCLUDE_SCHIFRA_GALOIS_FIELD_HPP


#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>
#include <string>


namespace schifra
{

   namespace galois
   {

      typedef int          field_symbol;
      const   field_symbol GFERROR = -1;

      class field
      {
      public:

         field(const int  pwr, const std::size_t primpoly_deg, const unsigned int* primitive_poly);
        ~field();

         bool operator==(const field& gf) const;
         bool operator!=(const field& gf) const;

         inline field_symbol index(const field_symbol value) const
         {
            return index_of_[value];
         }

         inline field_symbol alpha(const field_symbol value) const
         {
            return alpha_to_[value];
         }

         inline unsigned int size() const
         {
            return field_size_;
         }

         inline unsigned int pwr() const
         {
            return power_;
         }

         inline unsigned int mask() const
         {
            return field_size_;
         }

         inline field_symbol add(const field_symbol& a, const field_symbol& b) const
         {
            return (a ^ b);
         }

         inline field_symbol sub(const field_symbol& a, const field_symbol& b) const
         {
            return (a ^ b);
         }

         inline field_symbol normalize(field_symbol x) const
         {
            while (x < 0)
            {
               x += static_cast<field_symbol>(field_size_);
            }

            while (x >= static_cast<field_symbol>(field_size_))
            {
               x -= static_cast<field_symbol>(field_size_);
               x  = (x >> power_) + (x & field_size_);
            }

            return x;
         }

         inline field_symbol mul(const field_symbol& a, const field_symbol& b) const
         {
            #if !defined(NO_GFLUT)
               return mul_table_[a][b];
            #else
               if ((a == 0) || (b == 0))
                  return 0;
               else
                  return alpha_to_[normalize(index_of_[a] + index_of_[b])];
            #endif
         }

         inline field_symbol div(const field_symbol& a, const field_symbol& b) const
         {
            #if !defined(NO_GFLUT)
               return div_table_[a][b];
            #else
               if ((a == 0) || (b == 0))
                  return 0;
               else
                  return alpha_to_[normalize(index_of_[a] - index_of_[b] + field_size_)];
            #endif
         }

         inline field_symbol exp(const field_symbol& a, int n) const
         {
            #if !defined(NO_GFLUT)
               if (n >= 0)
                  return exp_table_[a][n & field_size_];
               else
               {
                  while (n < 0) n += field_size_;

                  return (n ? exp_table_[a][n] : 1);
               }
            #else
               if (a != 0)
               {
                  if (n < 0)
                  {
                     while (n < 0) n += field_size_;
                     return (n ? alpha_to_[normalize(index_of_[a] * n)] : 1);
                  }
                  else if (n)
                     return alpha_to_[normalize(index_of_[a] * static_cast<field_symbol>(n))];
                  else
                     return 1;
               }
               else
                  return 0;
            #endif
         }

         #ifdef LINEAR_EXP_LUT
         inline field_symbol* const linear_exp(const field_symbol& a) const
         {
            #if !defined(NO_GFLUT)
               static const field_symbol upper_bound = 2 * field_size_;
               if ((a >= 0) && (a <= upper_bound))
                  return linear_exp_table_[a];
               else
                  return reinterpret_cast<field_symbol*>(0);
            #else
               return reinterpret_cast<field_symbol*>(0);
            #endif
         }
         #endif

         inline field_symbol inverse(const field_symbol& val) const
         {
            #if !defined(NO_GFLUT)
               return mul_inverse_[val];
            #else
               return alpha_to_[normalize(field_size_ - index_of_[val])];
            #endif
         }

         inline unsigned int prim_poly_term(const unsigned int index) const
         {
            return prim_poly_[index];
         }

         friend std::ostream& operator << (std::ostream& os, const field& gf);

      private:

         field();
         field(const field& gfield);
         field& operator=(const field& gfield);

         void         generate_field(const unsigned int* prim_poly_);
         field_symbol gen_mul       (const field_symbol& a, const field_symbol& b) const;
         field_symbol gen_div       (const field_symbol& a, const field_symbol& b) const;
         field_symbol gen_exp       (const field_symbol& a, const std::size_t&  n) const;
         field_symbol gen_inverse   (const field_symbol& val) const;

         std::size_t create_array(char buffer_[],
                                  const std::size_t& length,
                                  const std::size_t offset,
                                  field_symbol** array);

         std::size_t create_2d_array(char buffer_[],
                                     std::size_t row_cnt, std::size_t col_cnt,
                                     const std::size_t offset,
                                     field_symbol*** array);
         unsigned int   power_;
         std::size_t    prim_poly_deg_;
         unsigned int   field_size_;
         unsigned int   prim_poly_hash_;
         unsigned int*  prim_poly_;
         field_symbol*  alpha_to_;    // aka exponential or anti-log
         field_symbol*  index_of_;    // aka log
         field_symbol*  mul_inverse_; // multiplicative inverse
         field_symbol** mul_table_;
         field_symbol** div_table_;
         field_symbol** exp_table_;
         field_symbol** linear_exp_table_;
         char*          buffer_;
      };

      field::field(const int  pwr, const std::size_t primpoly_deg, const unsigned int* primitive_poly)
      : power_(pwr),
        prim_poly_deg_(primpoly_deg),
        field_size_((1 << power_) - 1)
      {
         alpha_to_    = new field_symbol [field_size_ + 1];
         index_of_    = new field_symbol [field_size_ + 1];

         #if !defined(NO_GFLUT)

         #ifdef LINEAR_EXP_LUT
         static const std::size_t buffer_size = ((6 * (field_size_ + 1) * (field_size_ + 1)) + ((field_size_ + 1) * 2)) * sizeof(field_symbol);
         #else
         static const std::size_t buffer_size = ((4 * (field_size_ + 1) * (field_size_ + 1)) + ((field_size_ + 1) * 2)) * sizeof(field_symbol);
         #endif

         buffer_ = new char[buffer_size];
         std::size_t offset = 0;
         offset = create_2d_array(buffer_,(field_size_ + 1),(field_size_ + 1),offset,&mul_table_);
         offset = create_2d_array(buffer_,(field_size_ + 1),(field_size_ + 1),offset,&div_table_);
         offset = create_2d_array(buffer_,(field_size_ + 1),(field_size_ + 1),offset,&exp_table_);

         #ifdef LINEAR_EXP_LUT
         offset = create_2d_array(buffer_,(field_size_ + 1),(field_size_ + 1) * 2,offset,&linear_exp_table_);
         #else
         linear_exp_table_ = 0;
         #endif

         offset = create_array(buffer_,(field_size_ + 1) * 2,offset,&mul_inverse_);

         #else

           buffer_      = 0;
           mul_table_   = 0;
           div_table_   = 0;
           exp_table_   = 0;
           mul_inverse_ = 0;
           linear_exp_table_ = 0;

         #endif

         prim_poly_ = new unsigned int [prim_poly_deg_ + 1];

         for (unsigned int i = 0; i < (prim_poly_deg_ + 1); ++i)
         {
            prim_poly_[i] = primitive_poly[i];
         }

         prim_poly_hash_ = 0xAAAAAAAA;

         for (std::size_t i = 0; i < (prim_poly_deg_ + 1); ++i)
         {
            prim_poly_hash_ += ((i & 1) == 0) ? (  (prim_poly_hash_ <<  7) ^  primitive_poly[i] * (prim_poly_hash_ >> 3)) :
                                                (~((prim_poly_hash_ << 11) + (primitive_poly[i] ^ (prim_poly_hash_ >> 5))));
         }

         generate_field(primitive_poly);
      }

      field::~field()
      {
         if (0 !=  alpha_to_) { delete [] alpha_to_;  alpha_to_  = 0; }
         if (0 !=  index_of_) { delete [] index_of_;  index_of_  = 0; }
         if (0 != prim_poly_) { delete [] prim_poly_; prim_poly_ = 0; }

         #if !defined(NO_GFLUT)

         if (0 != mul_table_) { delete [] mul_table_; mul_table_ = 0; }
         if (0 != div_table_) { delete [] div_table_; div_table_ = 0; }
         if (0 != exp_table_) { delete [] exp_table_; exp_table_ = 0; }

         #ifdef LINEAR_EXP_LUT
         if (0 != linear_exp_table_) { delete [] linear_exp_table_; linear_exp_table_ = 0; }
         #endif

         if (0 != buffer_) { delete [] buffer_; buffer_ = 0; }

         #endif
      }

      inline bool field::operator==(const field& gf) const
      {
         return (
                  (this->power_ == gf.power_) &&
                  (this->prim_poly_hash_ == gf.prim_poly_hash_)
                );
      }

      inline bool field::operator!=(const field& gf) const
      {
         return !field::operator ==(gf);
      }

      inline void field::generate_field(const unsigned int* prim_poly)
      {
         /*
            Note: It is assumed that the degree of the primitive
                  polynomial will be equivelent to the m value as
                  in GF(2^m)
         */

         field_symbol mask = 1;

         alpha_to_[power_] = 0;

         for (field_symbol i = 0; i < static_cast<field_symbol>(power_); ++i)
         {
            alpha_to_[i]            = mask;
            index_of_[alpha_to_[i]] = i;

            if (prim_poly[i] != 0)
            {
               alpha_to_[power_] ^= mask;
            }

            mask <<= 1;
         }

         index_of_[alpha_to_[power_]] = power_;

         mask >>= 1;

         for (field_symbol i = power_ + 1; i < static_cast<field_symbol>(field_size_); ++i)
         {
            if (alpha_to_[i - 1] >= mask)
              alpha_to_[i] = alpha_to_[power_] ^ ((alpha_to_[i - 1] ^ mask) << 1);
            else
              alpha_to_[i] = alpha_to_[i - 1] << 1;

            index_of_[alpha_to_[i]] = i;
         }

         index_of_[0] = GFERROR;
         alpha_to_[field_size_] = 1;

         #if !defined(NO_GFLUT)

           for (field_symbol i = 0; i < static_cast<field_symbol>(field_size_ + 1); ++i)
           {
              for (field_symbol j = 0; j < static_cast<field_symbol>(field_size_ + 1); ++j)
              {
                 mul_table_[i][j] = gen_mul(i,j);
                 div_table_[i][j] = gen_div(i,j);
                 exp_table_[i][j] = gen_exp(i,j);
              }
           }

           #ifdef LINEAR_EXP_LUT
           for (field_symbol i = 0; i < static_cast<field_symbol>(field_size_ + 1); ++i)
           {
              for (int j = 0; j < static_cast<field_symbol>(2 * field_size_); ++j)
              {
                 linear_exp_table_[i][j] = gen_exp(i,j);
              }
           }
           #endif

           for (field_symbol i = 0; i < static_cast<field_symbol>(field_size_ + 1); ++i)
           {
              mul_inverse_[i] = gen_inverse(i);
              mul_inverse_[i + (field_size_ + 1)] = mul_inverse_[i];
           }

         #endif
      }

      inline field_symbol field::gen_mul(const field_symbol& a, const field_symbol& b) const
      {
         if ((a == 0) || (b == 0))
           return 0;
         else
           return alpha_to_[normalize(index_of_[a] + index_of_[b])];
      }

      inline field_symbol field::gen_div(const field_symbol& a, const field_symbol& b) const
      {
         if ((a == 0) || (b == 0))
           return 0;
         else
           return alpha_to_[normalize(index_of_[a] - index_of_[b] + field_size_)];
      }

      inline field_symbol field::gen_exp(const field_symbol& a, const std::size_t& n) const
      {
         if (a != 0)
           return ((n == 0) ? 1 : alpha_to_[normalize(index_of_[a] * static_cast<field_symbol>(n))]);
         else
           return 0;
      }

      inline field_symbol field::gen_inverse(const field_symbol& val) const
      {
         return alpha_to_[normalize(field_size_ - index_of_[val])];
      }

      std::size_t field::create_array(char buffer[],
                                      const std::size_t& length,
                                      const std::size_t offset,
                                      field_symbol** array)
      {
         const std::size_t row_size = length * sizeof(field_symbol);
         (*array) = new(buffer + offset)field_symbol[length];
         return row_size + offset;
      }

      std::size_t field::create_2d_array(char buffer[],
                                         std::size_t row_cnt, std::size_t col_cnt,
                                         const std::size_t offset,
                                         field_symbol*** array)
      {
         const std::size_t row_size = col_cnt * sizeof(field_symbol);
         char* buffer__offset = buffer + offset;
         (*array) = new field_symbol* [row_cnt];
         for (std::size_t i = 0; i < row_cnt; ++i)
         {
            (*array)[i] = new(buffer__offset + (i * row_size))field_symbol[col_cnt];
         }
         return (row_cnt * row_size) + offset;
      }

      inline std::ostream& operator << (std::ostream& os, const field& gf)
      {
         for (std::size_t i = 0; i < (gf.field_size_ + 1); ++i)
         {
            os << i << "\t" << gf.alpha_to_[i] << "\t" << gf.index_of_[i] << std::endl;
         }

         return os;
      }

      /* 1x^0 + 1x^1 + 0x^2 + 1x^3 */
      const unsigned int primitive_polynomial00[]    = {1, 1, 0, 1};
      const unsigned int primitive_polynomial_size00 = 4;

      /* 1x^0 + 1x^1 + 0x^2 + 0x^3 + 1x^4*/
      const unsigned int primitive_polynomial01[]    = {1, 1, 0, 0, 1};
      const unsigned int primitive_polynomial_size01 = 5;

      /* 1x^0 + 0x^1 + 1x^2 + 0x^3 + 0x^4 + 1x^5 */
      const unsigned int primitive_polynomial02[]    = {1, 0, 1, 0, 0, 1};
      const unsigned int primitive_polynomial_size02 = 6;

      /* 1x^0 + 1x^1 + 0x^2 + 0x^3 + 0x^4 + 0x^5 + 1x^6 */
      const unsigned int primitive_polynomial03[]    = {1, 1, 0, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size03 = 7;

      /* 1x^0 + 0x^1 + 0x^2 + 1x^3 + 0x^4 + 0x^5 + 0x^6 + 1x^7 */
      const unsigned int primitive_polynomial04[]    = {1, 0, 0, 1, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size04 = 8;

      /* 1x^0 + 0x^1 + 1x^2 + 1x^3 + 1x^4 + 0x^5 + 0x^6 + 0x^7 + 1x^8 */
      const unsigned int primitive_polynomial05[]    = {1, 0, 1, 1, 1, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size05 = 9;

      /* 1x^0 + 1x^1 + 1x^2 + 0x^3 + 0x^4 + 0x^5 + 0x^6 + 1x^7 + 1x^8 */
      const unsigned int primitive_polynomial06[]    = {1, 1, 1, 0, 0, 0, 0, 1, 1};
      const unsigned int primitive_polynomial_size06 = 9;

      /* 1x^0 + 0x^1 + 0x^2 + 0x^3 + 1x^4 + 0x^5 + 0x^6 + 0x^7 + 0x^8 + 1x^9 */
      const unsigned int primitive_polynomial07[]    = {1, 0, 0, 0, 1, 0, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size07 = 10;

      /* 1x^0 + 0x^1 + 0x^2 + 1x^3 + 0x^4 + 0x^5 + 0x^6 + 0x^7 + 0x^8 + 0x^9 + 1x^10 */
      const unsigned int primitive_polynomial08[]    = {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size08 = 11;

      /* 1x^0 + 0x^1 + 1x^2 + 0x^3 + 0x^4 + 0x^5 + 0x^6 + 0x^7 + 0x^8 + 0x^9 + 0x^10 + 1x^11 */
      const unsigned int primitive_polynomial09[]    = {1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size09 = 12;

      /* 1x^0 + 1x^1 + 0x^2 + 0x^3 + 1x^4 + 0x^5 + 1x^6 + 0x^7 + 0x^8 + 0x^9 + 0x^10 + 0x^11 + 1x^12 */
      const unsigned int primitive_polynomial10[]    = {1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size10 = 13;

      /* 1x^0 + 1x^1 + 0x^2 + 1x^3 + 1x^4 + 0x^5 + 0x^6 + 0x^7 + 0x^8 + 0x^9 + 0x^10 + 0x^11 + 0x^12 + 1x^13 */
      const unsigned int primitive_polynomial11[]    = {1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size11 = 14;

      /* 1x^0 + 1x^1 + 0x^2 + 0x^3 + 0x^4 + 0x^5 + 1x^6 + 0x^7 + 0x^8 + 0x^9 + 1x^10 + 0x^11 + 0x^12 + 0x^13 + 1x^14 */
      const unsigned int primitive_polynomial12[]    = {1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size12 = 15;

      /* 1x^0 + 1x^1 + 0x^2 + 0x^3 + 0x^4 + 0x^5 + 0x^6 + 0x^7 + 0x^8 + 0x^9 + 0x^10 + 0x^11 + 0x^12 + 0x^13 + 0x^14 + 1x^15 */
      const unsigned int primitive_polynomial13[]    = {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size13 = 16;

      /* 1x^0 + 1x^1 + 0x^2 + 1x^3 + 0x^4 + 0x^5 + 0x^6 + 0x^7 + 0x^8 + 0x^9 + 0x^10 + 0x^11 + 1x^12 + 0x^13 + 0x^14 + 0x^15 + 1x^16 */
      const unsigned int primitive_polynomial14[]    = {1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1};
      const unsigned int primitive_polynomial_size14 = 17;

   } // namespace galois

} // namespace schifra

#endif
