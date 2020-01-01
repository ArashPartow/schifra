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


#ifndef INCLUDE_SCHIFRA_GALOIS_UTILITIES_HPP
#define INCLUDE_SCHIFRA_GALOIS_UTILITIES_HPP


#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <string>
#include <sstream>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"


namespace schifra
{

   namespace galois
   {

      inline std::string convert_to_string(const unsigned int& value, const unsigned int& width)
      {
         std::stringstream stream;
         stream << std::setw(width) << std::setfill('0') << value;
         return stream.str();
      }

      inline std::string convert_to_string(const int& value, const unsigned int& width)
      {
         std::stringstream stream;
         stream << std::setw(width) << std::setfill('0') << value;
         return stream.str();
      }

      inline std::string convert_to_bin(const unsigned int& value, const unsigned int& field_descriptor)
      {
         std::string output = std::string(field_descriptor, ' ');

         for (unsigned int i = 0; i < field_descriptor; ++i)
         {
            output[i] = ((((value >> (field_descriptor - 1 - i)) & 1) == 1) ? '1' : '0');
         }

         return output;
      }

      inline void alpha_table(std::ostream& os, const field& gf)
      {
         std::vector<std::string> str_list;

         for (unsigned int i = 0; i < gf.size() + 1; ++i)
         {
            str_list.push_back("alpha^" + convert_to_string(gf.index(i),2) + "\t" +
                                          convert_to_bin   (i,gf.pwr())    + "\t" +
                                          convert_to_string(gf.alpha(i),2));
         }

         std::sort(str_list.begin(),str_list.end());
         std::copy(str_list.begin(),str_list.end(),std::ostream_iterator<std::string>(os,"\n"));
      }

      inline void polynomial_alpha_form(std::ostream& os, const field_polynomial& polynomial)
      {
         for (int i = 0; i < (polynomial.deg() + 1); ++i)
         {
            field_symbol alpha_power = polynomial.galois_field().index(polynomial[i].poly());

            if (alpha_power != 0)
               os << static_cast<unsigned char>(224) << "^" << convert_to_string(alpha_power,2);
            else
               os << 1;

            os << " * "
               << "x^"
               << i
               << ((i != (polynomial.deg())) ? " + " : "");
         }
      }

      inline void polynomial_alpha_form(std::ostream& os, const std::string& prepend, const field_polynomial& polynomial)
      {
         os << prepend;
         polynomial_alpha_form(os,polynomial);
         os << std::endl;
      }

   } // namespace reed_solomon

} // namespace schifra

#endif
