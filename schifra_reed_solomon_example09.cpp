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


/*
   Description: This example will demonstrate how to instantiate a 16-bit per
                symbol Reed-Solomon encoder and decoder, add the full amount
                of possible errors, correct the errors, and output the various
                pieces of relevant information. Furthermore this example will
                demonstrate the use of the Reed-Solomon codec without the use
                of LUTs for the finite field computations.
*/


#include <cstddef>
#include <iostream>
#include <string>

#define NO_GFLUT
#include "schifra_galois_field.hpp"
#undef NO_GFLUT
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_error_processes.hpp"


int main()
{
   /* Finite Field Parameters */
   const std::size_t field_descriptor                =   16;
   const std::size_t generator_polynomial_index      =    0;
   const std::size_t generator_polynomial_root_count = 1000;

   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 65535;
   const std::size_t fec_length  =  1000;
   const std::size_t data_length = code_length - fec_length;

   /* Instantiate Finite Field and Generator Polynomials */
   const schifra::galois::field field(field_descriptor,
                                      schifra::galois::primitive_polynomial_size14,
                                      schifra::galois::primitive_polynomial14);

   schifra::galois::field_polynomial generator_polynomial(field);

   if (
        !schifra::make_sequential_root_generator_polynomial(field,
                                                            generator_polynomial_index,
                                                            generator_polynomial_root_count,
                                                            generator_polynomial)
      )
   {
      std::cout << "Error - Failed to create sequential root generator!" << std::endl;
      return 1;
   }

   /* Instantiate Encoder and Decoder (Codec) */
   typedef schifra::reed_solomon::encoder<code_length,fec_length,data_length> encoder_t;
   typedef schifra::reed_solomon::decoder<code_length,fec_length,data_length> decoder_t;

   const encoder_t encoder(field, generator_polynomial);
   const decoder_t decoder(field, generator_polynomial_index);

   std::vector<unsigned int> data(code_length, 0x0000 & field.mask());

   /* Instantiate RS Block For Codec */
   schifra::reed_solomon::block<code_length,fec_length> block;
   schifra::reed_solomon::block<code_length,fec_length> original_block;

   schifra::reed_solomon::copy(&data[0], data.size(), block);

   original_block = block;

   /* Transform message into Reed-Solomon encoded codeword */
   if (!encoder.encode(block))
   {
      std::cout << "Error - Critical encoding failure! "
                << "Msg: " << block.error_as_string()  << std::endl;
      return 1;
   }

   /* Add errors at every 3rd location starting at position zero */
   schifra::corrupt_message_all_errors_wth_mask(block, 0, field.mask(), 3);

   if (!decoder.decode(block))
   {
      std::cout << "Error - Critical decoding failure! "
                << "Msg: " << block.error_as_string()  << std::endl;
      return 1;
   }
   else if (!schifra::are_blocks_equivelent(block, original_block, true, true))
   {
      std::cout << "Error - Error correction failed!" << std::endl;
      return 1;
   }

   std::cout << "Encoder Parameters [" << encoder_t::trait::code_length << ","
                                       << encoder_t::trait::data_length << ","
                                       << encoder_t::trait::fec_length  << "]" << std::endl;

   std::cout << "Decoder Parameters [" << decoder_t::trait::code_length << ","
                                       << decoder_t::trait::data_length << ","
                                       << decoder_t::trait::fec_length  << "]" << std::endl;

   return 0;
}
