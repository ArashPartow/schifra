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


#ifndef INCLUDE_SCHIFRA_SEQUENTIAL_ROOT_GENERATOR_POLYNOMIAL_CREATOR_HPP
#define INCLUDE_SCHIFRA_SEQUENTIAL_ROOT_GENERATOR_POLYNOMIAL_CREATOR_HPP


#include <cstddef>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_element.hpp"
#include "schifra_galois_field_polynomial.hpp"


namespace schifra
{

   inline bool make_sequential_root_generator_polynomial(const galois::field& field,
                                                         const std::size_t initial_index,
                                                         const std::size_t num_elements,
                                                         galois::field_polynomial& generator_polynomial)
   {
      if (
           (initial_index >= field.size()) ||
           ((initial_index + num_elements) >  field.size())
         )
      {
         return false;
      }

      galois::field_element alpha(field, 2);
      galois::field_polynomial X = galois::generate_X(field);
      generator_polynomial = galois::field_element(field, 1);

      for (std::size_t i = initial_index; i < (initial_index + num_elements); ++i)
      {
         generator_polynomial *= (X + (alpha ^ static_cast<galois::field_symbol>(i)));
      }

      return true;
   }

} // namespace schifra

#endif
