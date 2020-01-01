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


#include <cstddef>
#include <string>

#include "schifra_galois_field.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_file_decoder.hpp"


int main()
{
   const std::size_t field_descriptor    =   8;
   const std::size_t gen_poly_index      = 120;
   const std::size_t code_length         = 255;
   const std::size_t fec_length          =   6;
   const std::string input_file_name     = "input.schifra";
   const std::string output_file_name    = "output.decoded";

   typedef schifra::reed_solomon::decoder<code_length,fec_length> decoder_t;
   typedef schifra::reed_solomon::file_decoder<code_length,fec_length> file_decoder_t;

   const schifra::galois::field field(field_descriptor,
                                      schifra::galois::primitive_polynomial_size06,
                                      schifra::galois::primitive_polynomial06);

   const decoder_t rs_decoder(field,gen_poly_index);

   file_decoder_t(rs_decoder, input_file_name, output_file_name);

   return 0;
}
