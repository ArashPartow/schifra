/*
(**************************************************************************)
(*                                                                        *)
(*                                Schifra                                 *)
(*                Reed-Solomon Error Correcting Code Library              *)
(*                                                                        *)
(* Release Version 0.0.1                                                  *)
(* http://www.schifra.com                                                 *)
(* Copyright (c) 2000-2015 Arash Partow, All Rights Reserved.             *)
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
   Description: This example will demonstrate how to instantiate a Reed-Solomon
                encoder and decoder, add the full amount of possible erasures,
                correct the erasures, and output the various pieces of relevant
                information.
*/


#include <cstddef>
#include <iostream>
#include <string>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_error_processes.hpp"


int main()
{
   /* Finite Field Parameters */
   const std::size_t field_descriptor                 =   8;
   const std::size_t generator_polynommial_index      = 120;
   const std::size_t generator_polynommial_root_count =  32;

   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 255;
   const std::size_t fec_length  =  32;
   const std::size_t data_length = code_length - fec_length;

   /* Instantiate Finite Field and Generator Polynomials */
   schifra::galois::field field(field_descriptor,
                                schifra::galois::primitive_polynomial_size06,
                                schifra::galois::primitive_polynomial06);

   schifra::galois::field_polynomial generator_polynomial(field);

   schifra::sequential_root_generator_polynomial_creator(field,
                                                         generator_polynommial_index,
                                                         generator_polynommial_root_count,
                                                         generator_polynomial);

   /* Instantiate Encoder and Decoder (Codec) */
   schifra::reed_solomon::encoder<code_length,fec_length> encoder(field,generator_polynomial);
   schifra::reed_solomon::decoder<code_length,fec_length> decoder(field,generator_polynommial_index);

   std::string message = "A professional is a person who knows more and more about less and less until they know everything about nothing";
               message = message + std::string(data_length - message.length(),static_cast<unsigned char>(0x00));

   std::cout << "Original Message:   [" << message << "]" << std::endl;

   /* Instantiate RS Block For Codec */
   schifra::reed_solomon::block<code_length,fec_length> block;

   /* Transform message into Reed-Solomon encoded codeword */
   if (!encoder.encode(message,block))
   {
      std::cout << "Error - Critical encoding failure!" << std::endl;
      return 1;
   }

   schifra::reed_solomon::erasure_locations_t erasure_location_list;

   /* Add erasures at every 2nd location starting at position zero */
   schifra::corrupt_message_all_erasures00(block,erasure_location_list,0,2);

   std::cout << "Corrupted Codeword: [" << block << "]" << std::endl;

   if (!decoder.decode(block,erasure_location_list))
   {
      std::cout << "Error - Critical decoding failure!" << std::endl;
      return 1;
   }
   else if (!schifra::is_block_equivelent(block,message))
   {
      std::cout << "Error - Error correction failed!" << std::endl;
      return 1;
   }

   block.data_to_string(message);

   std::cout << "Corrected Message:  [" << message << "]" << std::endl;

   std::cout << "Encoder Parameters [" << schifra::reed_solomon::encoder<code_length,fec_length>::trait::code_length << ","
                                       << schifra::reed_solomon::encoder<code_length,fec_length>::trait::data_length << ","
                                       << schifra::reed_solomon::encoder<code_length,fec_length>::trait::fec_length  << "]" << std::endl;

   std::cout << "Decoder Parameters [" << schifra::reed_solomon::decoder<code_length,fec_length>::trait::code_length << ","
                                       << schifra::reed_solomon::decoder<code_length,fec_length>::trait::data_length << ","
                                       << schifra::reed_solomon::decoder<code_length,fec_length>::trait::fec_length  << "]" << std::endl;

   return 0;
}
