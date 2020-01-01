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


#ifndef INCLUDE_SCHIFRA_GALOIS_FIELD_POLYNOMIAL_HPP
#define INCLUDE_SCHIFRA_GALOIS_FIELD_POLYNOMIAL_HPP


#include <cassert>
#include <iostream>
#include <vector>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_element.hpp"


namespace schifra
{

   namespace galois
   {

      class field_polynomial
      {
      public:

         field_polynomial(const field& gfield);
         field_polynomial(const field& gfield, const unsigned int& degree);
         field_polynomial(const field& gfield, const unsigned int& degree, const field_element element[]);
         field_polynomial(const field_polynomial& polynomial);
         field_polynomial(const field_element& gfe);
        ~field_polynomial() {}

         bool valid() const;
         int deg() const;
         const field& galois_field() const;
         void set_degree(const unsigned int& x);
         void simplify();

         field_polynomial& operator  =  (const field_polynomial& polynomial);
         field_polynomial& operator  =  (const field_element&       element);
         field_polynomial& operator +=  (const field_polynomial&    element);
         field_polynomial& operator +=  (const field_element&       element);
         field_polynomial& operator -=  (const field_polynomial&    element);
         field_polynomial& operator -=  (const field_element&       element);
         field_polynomial& operator *=  (const field_polynomial& polynomial);
         field_polynomial& operator *=  (const field_element&       element);
         field_polynomial& operator /=  (const field_polynomial&    divisor);
         field_polynomial& operator /=  (const field_element&       element);
         field_polynomial& operator %=  (const field_polynomial&    divisor);
         field_polynomial& operator %=  (const unsigned int&          power);
         field_polynomial& operator ^=  (const unsigned int&              n);
         field_polynomial& operator <<= (const unsigned int&              n);
         field_polynomial& operator >>= (const unsigned int&              n);

         field_element&    operator[] (const std::size_t&            term);
         field_element     operator() (const field_element&         value);
         field_element     operator() (field_symbol                 value);

         const field_element& operator[](const std::size_t&    term) const;
         const field_element  operator()(const field_element& value) const;
         const field_element  operator()(field_symbol         value) const;

         bool operator==(const field_polynomial& polynomial) const;
         bool operator!=(const field_polynomial& polynomial) const;

         bool monic() const;

         field_polynomial derivative() const;

         friend std::ostream& operator << (std::ostream& os, const field_polynomial& polynomial);

      private:

         typedef std::vector<field_element>::iterator       poly_iter;
         typedef std::vector<field_element>::const_iterator const_poly_iter;

         void simplify(field_polynomial& polynomial) const;

         field& field_;
         std::vector<field_element> poly_;
      };

      field_polynomial operator + (const field_polynomial& a, const field_polynomial& b);
      field_polynomial operator + (const field_polynomial& a, const field_element&    b);
      field_polynomial operator + (const field_element&    a, const field_polynomial& b);
      field_polynomial operator + (const field_polynomial& a, const field_symbol&     b);
      field_polynomial operator + (const field_symbol&     a, const field_polynomial& b);
      field_polynomial operator - (const field_polynomial& a, const field_polynomial& b);
      field_polynomial operator - (const field_polynomial& a, const field_element&    b);
      field_polynomial operator - (const field_element&    a, const field_polynomial& b);
      field_polynomial operator - (const field_polynomial& a, const field_symbol&     b);
      field_polynomial operator - (const field_symbol&     a, const field_polynomial& b);
      field_polynomial operator * (const field_polynomial& a, const field_polynomial& b);
      field_polynomial operator * (const field_element&    a, const field_polynomial& b);
      field_polynomial operator * (const field_polynomial& a, const field_element&    b);
      field_polynomial operator / (const field_polynomial& a, const field_polynomial& b);
      field_polynomial operator / (const field_polynomial& a, const field_element&    b);
      field_polynomial operator % (const field_polynomial& a, const field_polynomial& b);
      field_polynomial operator % (const field_polynomial& a, const unsigned int& power);
      field_polynomial operator ^ (const field_polynomial& a, const int&              n);
      field_polynomial operator <<(const field_polynomial& a, const unsigned int&     n);
      field_polynomial operator >>(const field_polynomial& a, const unsigned int&     n);
      field_polynomial         gcd(const field_polynomial& a, const field_polynomial& b);

      inline field_polynomial::field_polynomial(const field& gfield)
      : field_(const_cast<field&>(gfield))
      {
         poly_.clear();
         poly_.reserve(256);
      }

      inline field_polynomial::field_polynomial(const field& gfield, const unsigned int& degree)
      : field_(const_cast<field&>(gfield))
      {
         poly_.reserve(256);
         poly_.resize(degree + 1,field_element(field_,0));
      }

      inline field_polynomial::field_polynomial(const field& gfield, const unsigned int& degree, const field_element element[])
      : field_(const_cast<field&>(gfield))
      {
         poly_.reserve(256);

         if (element != NULL)
         {
            /*
               It is assumed that element is an array of field elements
               with size/element count of degree + 1.
            */
            for (unsigned int i = 0; i <= degree; ++i)
            {
               poly_.push_back(element[i]);
            }
         }
         else
            poly_.resize(degree + 1, field_element(field_, 0));
      }

      inline field_polynomial::field_polynomial(const field_polynomial& polynomial)
      : field_(const_cast<field&>(polynomial.field_)),
        poly_ (polynomial.poly_)
      {}

      inline field_polynomial::field_polynomial(const field_element& element)
      : field_(const_cast<field&>(element.galois_field()))
      {
         poly_.resize(1,element);
      }

      inline bool field_polynomial::valid() const
      {
         return (poly_.size() > 0);
      }

      inline int field_polynomial::deg() const
      {
         return static_cast<int>(poly_.size()) - 1;
      }

      inline const field& field_polynomial::galois_field() const
      {
         return field_;
      }

      inline void field_polynomial::set_degree(const unsigned int& x)
      {
         poly_.resize(x - 1,field_element(field_,0));
      }

      inline field_polynomial& field_polynomial::operator = (const field_polynomial& polynomial)
      {
         if ((this != &polynomial) && (&field_ == &(polynomial.field_)))
         {
            poly_ = polynomial.poly_;
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator = (const field_element& element)
      {
         if (&field_ == &(element.galois_field()))
         {
            poly_.resize(1,element);
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator += (const field_polynomial& polynomial)
      {
         if (&field_ == &(polynomial.field_))
         {
            if (poly_.size() < polynomial.poly_.size())
            {
               const_poly_iter it0 = polynomial.poly_.begin();

               for (poly_iter it1 = poly_.begin(); it1 != poly_.end(); ++it0, ++it1)
               {
                  (*it1) += (*it0);
               }

               while (it0 != polynomial.poly_.end())
               {
                  poly_.push_back(*it0);
                  ++it0;
               }
            }
            else
            {
               poly_iter it0 = poly_.begin();

               for (const_poly_iter it1 = polynomial.poly_.begin(); it1 != polynomial.poly_.end(); ++it0, ++it1)
               {
                  (*it0) += (*it1);
               }
            }

            simplify(*this);
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator += (const field_element& element)
      {
         poly_[0] += element;
         return *this;
      }

      inline field_polynomial& field_polynomial::operator -= (const field_polynomial& element)
      {
         return (*this += element);
      }

      inline field_polynomial& field_polynomial::operator -= (const field_element& element)
      {
         poly_[0] -= element;
         return *this;
      }

      inline field_polynomial& field_polynomial::operator *= (const field_polynomial& polynomial)
      {
         if (&field_ == &(polynomial.field_))
         {
            field_polynomial product(field_,deg() + polynomial.deg() + 1);

            poly_iter result_it = product.poly_.begin();

            for (poly_iter it0 = poly_.begin(); it0 != poly_.end(); ++it0)
            {
               poly_iter current_result_it = result_it;

               for (const_poly_iter it1 = polynomial.poly_.begin(); it1 != polynomial.poly_.end(); ++it1)
               {
                  (*current_result_it) += (*it0) * (*it1);
                  ++current_result_it;
               }

               ++result_it;
            }

            simplify(product);
            poly_ = product.poly_;
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator *= (const field_element& element)
      {
         if (field_ == element.galois_field())
         {
            for (poly_iter it = poly_.begin(); it != poly_.end(); ++it)
            {
               (*it) *= element;
            }
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator /= (const field_polynomial& divisor)
      {
         if (
             (&field_       == &divisor.field_) &&
             (deg()         >=   divisor.deg()) &&
             (divisor.deg() >=               0)
            )
         {
            field_polynomial quotient (field_, deg() - divisor.deg() + 1);
            field_polynomial remainder(field_, divisor.deg() - 1);

            for (int i = static_cast<int>(deg()); i >= 0; i--)
            {
               if (i <= static_cast<int>(quotient.deg()))
               {
                  quotient[i] = remainder[remainder.deg()] / divisor[divisor.deg()];

                  for (int j = static_cast<int>(remainder.deg()); j > 0; --j)
                  {
                     remainder[j] = remainder[j - 1] + (quotient[i] * divisor[j]);
                  }

                  remainder[0] = poly_[i] + (quotient[i] * divisor[0]);
               }
               else
               {
                  for (int j = static_cast<int>(remainder.deg()); j > 0; --j)
                  {
                     remainder[j] = remainder[j - 1];
                  }

                  remainder[0] = poly_[i];
               }
            }

            simplify(quotient);
            poly_ = quotient.poly_;
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator /= (const field_element& element)
      {
         if (field_ == element.galois_field())
         {
            for (poly_iter it = poly_.begin(); it != poly_.end(); ++it)
            {
               (*it) /= element;
            }
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator %= (const field_polynomial& divisor)
      {
         if (
              (field_        == divisor.field_) &&
              (deg()         >= divisor.deg() ) &&
              (divisor.deg() >=             0 )
            )
         {
            field_polynomial quotient (field_, deg() - divisor.deg() + 1);
            field_polynomial remainder(field_, divisor.deg() - 1);

            for (int i = static_cast<int>(deg()); i >= 0; i--)
            {
               if (i <= static_cast<int>(quotient.deg()))
               {
                  quotient[i] = remainder[remainder.deg()] / divisor[divisor.deg()];

                  for (int j = static_cast<int>(remainder.deg()); j > 0; --j)
                  {
                     remainder[j] = remainder[j - 1] + (quotient[i] * divisor[j]);
                  }

                  remainder[0] = poly_[i] + (quotient[i] * divisor[0]);
               }
               else
               {
                  for (int j = static_cast<int>(remainder.deg()); j > 0; --j)
                  {
                     remainder[j] = remainder[j - 1];
                  }

                  remainder[0] = poly_[i];
               }
            }

            poly_ = remainder.poly_;
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator %= (const unsigned int& power)
      {
         if (poly_.size() >= power)
         {
            poly_.resize(power,field_element(field_,0));
            simplify(*this);
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator ^= (const unsigned int& n)
      {
         field_polynomial result = *this;

         for (std::size_t i = 0; i < n; ++i)
         {
            result *= *this;
         }

         *this = result;

         return *this;
      }

      inline field_polynomial& field_polynomial::operator <<= (const unsigned int& n)
      {
         if (poly_.size() > 0)
         {
            size_t initial_size = poly_.size();

            poly_.resize(poly_.size() + n, field_element(field_,0));

            for (size_t i = initial_size - 1; static_cast<int>(i) >= 0; --i)
            {
               poly_[i + n] = poly_[i];
            }

            for (unsigned int i = 0; i < n; ++i)
            {
               poly_[i] = 0;
            }
         }

         return *this;
      }

      inline field_polynomial& field_polynomial::operator >>= (const unsigned int& n)
      {
         if (n <= poly_.size())
         {
            for (unsigned int i = 0; i <= deg() - n; ++i)
            {
               poly_[i] = poly_[i + n];
            }

            poly_.resize(poly_.size() - n,field_element(field_,0));
         }
         else if (static_cast<int>(n) >= (deg() + 1))
         {
            poly_.resize(0,field_element(field_,0));
         }

         return *this;
      }

      inline const field_element& field_polynomial::operator [] (const std::size_t& term) const
      {
         assert(term < poly_.size());
         return poly_[term];
      }

      inline field_element& field_polynomial::operator [] (const std::size_t& term)
      {
         assert(term < poly_.size());
         return poly_[term];
      }

      inline field_element field_polynomial::operator () (const field_element& value)
      {
         field_element result(field_,0);

         if (!poly_.empty())
         {
            int i = 0;
            field_symbol total_sum = 0 ;
            field_symbol value_poly_form = value.poly();

            for (poly_iter it = poly_.begin(); it != poly_.end(); ++it, ++i)
            {
               total_sum ^= field_.mul(field_.exp(value_poly_form,i), (*it).poly());
            }

            result = total_sum;
         }

         return result;
      }

      inline const field_element field_polynomial::operator () (const field_element& value) const
      {
         if (!poly_.empty())
         {
            int i = 0;
            field_symbol total_sum = 0 ;
            field_symbol value_poly_form = value.poly();

            for (const_poly_iter it = poly_.begin(); it != poly_.end(); ++it, ++i)
            {
               total_sum ^= field_.mul(field_.exp(value_poly_form,i), (*it).poly());
            }

            return field_element(field_,total_sum);
         }

         return field_element(field_,0);
      }

      inline field_element field_polynomial::operator () (field_symbol value)
      {
         if (!poly_.empty())
         {
            int i = 0;
            field_symbol total_sum = 0 ;

            for (const_poly_iter it = poly_.begin(); it != poly_.end(); ++it, ++i)
            {
               total_sum ^= field_.mul(field_.exp(value,i), (*it).poly());
            }

            return field_element(field_,total_sum);
         }

         return field_element(field_,0);
      }

      inline const field_element field_polynomial::operator () (field_symbol value) const
      {
         if (!poly_.empty())
         {
            int i = 0;
            field_symbol total_sum = 0 ;

            for (const_poly_iter it = poly_.begin(); it != poly_.end(); ++it, ++i)
            {
               total_sum ^= field_.mul(field_.exp(value, i), (*it).poly());
            }

            return field_element(field_,total_sum);
         }

         return field_element(field_,0);
      }

      inline bool field_polynomial::operator == (const field_polynomial& polynomial) const
      {
         if (field_ == polynomial.field_)
         {
            if (poly_.size() != polynomial.poly_.size())
              return false;
            else
            {
               const_poly_iter it0 = polynomial.poly_.begin();

               for (const_poly_iter it1 = poly_.begin(); it1 != poly_.end(); ++it0, ++it1)
               {
                  if ((*it0) != (*it1))
                    return false;
               }

               return true;
            }
         }
         else
           return false;
      }

      inline bool field_polynomial::operator != (const field_polynomial& polynomial) const
      {
         return !(*this == polynomial);
      }

      inline field_polynomial field_polynomial::derivative() const
      {
         if ((*this).poly_.size() > 1)
         {
            field_polynomial deriv(field_,deg());

            const std::size_t upper_bound = poly_.size() - 1;

            for (std::size_t i = 0; i < upper_bound; i += 2)
            {
               deriv.poly_[i] = poly_[i + 1];
            }

            simplify(deriv);
            return deriv;
         }

         return field_polynomial(field_,0);
      }

      inline bool field_polynomial::monic() const
      {
         return (poly_[poly_.size() - 1] == static_cast<galois::field_symbol>(1));
      }

      inline void field_polynomial::simplify()
      {
         simplify(*this);
      }

      inline void field_polynomial::simplify(field_polynomial& polynomial) const
      {
         std::size_t poly_size = polynomial.poly_.size();

         if ((poly_size > 0) && (polynomial.poly_.back() == 0))
         {
            poly_iter it    = polynomial.poly_.end  ();
            poly_iter begin = polynomial.poly_.begin();

            std::size_t count = 0;

            while ((begin != it) && (*(--it) == 0))
            {
               ++count;
            }

            if (0 != count)
            {
               polynomial.poly_.resize(poly_size - count, field_element(field_,0));
            }
         }
      }

      inline field_polynomial operator + (const field_polynomial& a, const field_polynomial& b)
      {
         field_polynomial result = a;
         result += b;
         return result;
      }

      inline field_polynomial operator + (const field_polynomial& a, const field_element& b)
      {
         field_polynomial result = a;
         result += b;
         return result;
      }

      inline field_polynomial operator + (const field_element& a, const field_polynomial& b)
      {
         field_polynomial result = b;
         result += a;
         return result;
      }

      inline field_polynomial operator + (const field_polynomial& a, const field_symbol& b)
      {
         return a + field_element(a.galois_field(),b);
      }

      inline field_polynomial operator + (const field_symbol& a, const field_polynomial& b)
      {
         return b + field_element(b.galois_field(),a);
      }

      inline field_polynomial operator - (const field_polynomial& a, const field_polynomial& b)
      {
         field_polynomial result = a;
         result -= b;
         return result;
      }

      inline field_polynomial operator - (const field_polynomial& a, const field_element& b)
      {
         field_polynomial result = a;
         result -= b;
         return result;
      }

      inline field_polynomial operator - (const field_element& a, const field_polynomial& b)
      {
         field_polynomial result = b;
         result -= a;
         return result;
      }

      inline field_polynomial operator - (const field_polynomial& a, const field_symbol& b)
      {
         return a - field_element(a.galois_field(),b);
      }

      inline field_polynomial operator - (const field_symbol& a, const field_polynomial& b)
      {
         return b - field_element(b.galois_field(),a);
      }

      inline field_polynomial operator * (const field_polynomial& a, const field_polynomial& b)
      {
         field_polynomial result = a;
         result *= b;
         return result;
      }

      inline field_polynomial operator * (const field_element& a, const field_polynomial& b)
      {
         field_polynomial result = b;
         result *= a;
         return result;
      }

      inline field_polynomial operator * (const field_polynomial& a, const field_element& b)
      {
         field_polynomial result = a;
         result *= b;
         return result;
      }

      inline field_polynomial operator / (const field_polynomial& a, const field_polynomial& b)
      {
         field_polynomial result = a;
         result /= b;
         return result;
      }

      inline field_polynomial operator / (const field_polynomial& a, const field_element& b)
      {
         field_polynomial result = a;
         result /= b;
         return result;
      }

      inline field_polynomial operator % (const field_polynomial& a, const field_polynomial& b)
      {
         field_polynomial result = a;
         result %= b;
         return result;
      }

      inline field_polynomial operator % (const field_polynomial& a, const unsigned int& n)
      {
         field_polynomial result = a;
         result %= n;
         return result;
      }

      inline field_polynomial operator ^ (const field_polynomial& a, const int& n)
      {
         field_polynomial result = a;
         result ^= n;
         return result;
      }

      inline field_polynomial operator << (const field_polynomial& a, const unsigned int& n)
      {
         field_polynomial result = a;
         result <<= n;
         return result;
      }

      inline field_polynomial operator >> (const field_polynomial& a, const unsigned int& n)
      {
         field_polynomial result = a;
         result >>= n;
         return result;
      }

      inline field_polynomial gcd(const field_polynomial& a, const field_polynomial& b)
      {
         if (&a.galois_field() == &b.galois_field())
         {
            if ((!a.valid()) && (!b.valid()))
            {
               field_polynomial error_polynomial(a.galois_field());
               return error_polynomial;
            }

            if (!a.valid()) return b;
            if (!b.valid()) return a;

            field_polynomial x = a % b;
            field_polynomial y = b;
            field_polynomial z = x;

            while ((z = (y % x)).valid())
            {
               y = x;
               x = z;
            }
            return x;
         }
         else
         {
            field_polynomial error_polynomial(a.galois_field());
            return error_polynomial;
         }
      }

      inline field_polynomial generate_X(const field& gfield)
      {
         const field_element xgfe[2] = {
                                         galois::field_element(gfield, 0),
                                         galois::field_element(gfield, 1)
                                       };

         field_polynomial X_(gfield,1,xgfe);

         return X_;
      }

      inline std::ostream& operator << (std::ostream& os, const field_polynomial& polynomial)
      {
         if (polynomial.deg() >= 0)
         {
            /*
            for (unsigned int i = 0; i < polynomial.poly_.size(); ++i)
            {
               os << polynomial.poly[i].index()
                  << ((i != (polynomial.deg())) ? " " : "");
            }

            std::cout << " poly form: ";
            */

            for (unsigned int i = 0; i < polynomial.poly_.size(); ++i)
            {
               os << polynomial.poly_[i].poly()
                  << " "
                  << "x^"
                  << i
                  << ((static_cast<int>(i) != (polynomial.deg())) ? " + " : "");
            }
         }

         return os;
      }

   } // namespace galois

} // namespace schifra

#endif
