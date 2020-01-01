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
   Description: This example will demonstrate how to instantiate a 12-bit per
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
#include "schifra_reed_solomon_bitio.hpp"


int main()
{
   /* Finite Field Parameters */
   const std::size_t field_descriptor                =  12;
   const std::size_t generator_polynomial_index      = 120;
   const std::size_t generator_polynomial_root_count = 100;

   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 4095; //(2^12 - 1)
   const std::size_t fec_length  =  100;
   const std::size_t data_length = code_length - fec_length;

   /* 12-bit Symbol Mask Parameter */
   const int mask = 0x00000FFF;

   /*
     Note: One must make sure to be using primitive polynomials
           of correct degree for generating elements in the specified
           field ie: a primitive polynomial of degree 12 for GF(2^12).
   */

   /* Instantiate Finite Field and Generator Polynomials */
   const schifra::galois::field field(field_descriptor,
                                      schifra::galois::primitive_polynomial_size10,
                                      schifra::galois::primitive_polynomial10);

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
   typedef schifra::reed_solomon::encoder<code_length,fec_length> encoder_t;
   typedef schifra::reed_solomon::decoder<code_length,fec_length> decoder_t;

   const encoder_t encoder(field, generator_polynomial);
   const decoder_t decoder(field, generator_polynomial_index);

   /*
     Note: The data length represents the number of code symbols that will be used.
           The effective data length is then the number of code symbols multipled
           by the number of bits per symbol, in this case it is 12-bits per code
           symbol.
   */

   /* Instantiate RS Block For Codec */
   schifra::reed_solomon::block<code_length,fec_length> block;
   schifra::reed_solomon::block<code_length,fec_length> original_block;

   /* Populate RS Blocks with 12-bit data symbols */
   for (std::size_t i = 0; i < data_length; ++i)
   {
      block[i] = static_cast<int>(i & mask);
      original_block[i] = static_cast<int>(i & mask);
   }

   /* Transform message into Reed-Solomon encoded codeword */
   if (!encoder.encode(block))
   {
      std::cout << "Error - Critical encoding failure! "
                << "Msg: " << block.error_as_string()  << std::endl;
      return 1;
   }

   /* Add errors at every 3rd location starting at position zero */
   schifra::corrupt_message_all_errors_wth_mask(block, 0, mask, 3);

   if (!decoder.decode(block))
   {
      std::cout << "Error - Critical decoding failure! "
                << "Msg: " << block.error_as_string()  << std::endl;
      return 1;
   }
   else if (!schifra::are_blocks_equivelent(block, original_block, data_length))
   {
      std::cout << "Error - Error correction failed!" << std::endl;
      return 1;
   }

   std::cout << "Encoder Parameters [" << schifra::reed_solomon::encoder<code_length,fec_length>::trait::code_length << ","
                                       << schifra::reed_solomon::encoder<code_length,fec_length>::trait::data_length << ","
                                       << schifra::reed_solomon::encoder<code_length,fec_length>::trait::fec_length  << "]" << std::endl;

   std::cout << "Decoder Parameters [" << schifra::reed_solomon::decoder<code_length,fec_length>::trait::code_length << ","
                                       << schifra::reed_solomon::decoder<code_length,fec_length>::trait::data_length << ","
                                       << schifra::reed_solomon::decoder<code_length,fec_length>::trait::fec_length  << "]" << std::endl;

   return 0;
}
