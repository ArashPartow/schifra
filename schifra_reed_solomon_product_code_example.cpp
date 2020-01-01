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
   Description: This example demonstrates how to implement a simple Reed-Solomon
                based product-code encoder and decoder.
*/


#include <cstddef>
#include <cstdio>
#include <iostream>
#include <vector>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_product_code.hpp"
#include "schifra_utilities.hpp"


int main()
{
   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 255;
   const std::size_t fec_length  =  32;
   const std::size_t data_length = code_length - fec_length;

   /* Finite Field Parameters */
   const std::size_t field_descriptor                =   8;
   const std::size_t generator_polynomial_index      = 120;
   const std::size_t generator_polynomial_root_count = fec_length;

   /* Input/ Output Data Lengths */
   const std::size_t input_data_length  = data_length * data_length;
   const std::size_t output_data_length = code_length * code_length;

   /* Instantiate Finite Field and Generator Polynomials */
   const schifra::galois::field field(field_descriptor,
                                      schifra::galois::primitive_polynomial_size06,
                                      schifra::galois::primitive_polynomial06);

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

   schifra::reed_solomon::square_product_code_encoder<code_length,fec_length> spc_encoder(encoder);
   schifra::reed_solomon::square_product_code_decoder<code_length,fec_length> spc_decoder(decoder);

   unsigned char* input_data_01 = new unsigned char[ input_data_length];
   unsigned char* input_data_02 = new unsigned char[ input_data_length];
   unsigned char* output_data   = new unsigned char[output_data_length];

   const std::size_t max_iterations = 1000;

   for (std::size_t round = 0; round < max_iterations; ++round)
   {
      for (std::size_t i = 0; i < input_data_length; ++i)
      {
         input_data_01[i] = static_cast<unsigned char>(i);
         input_data_02[i] = 0;
      }

      spc_encoder.clear();
      spc_decoder.clear();

      schifra::utils::timer block_timer;
      block_timer.start();

      spc_encoder.encode_and_interleave(input_data_01);
      spc_encoder.output(output_data);

      const std::size_t max_error_count = (code_length * (fec_length >> 1));
      std::size_t error_count   = 0;
      std::size_t error_modulus = (static_cast<std::size_t>(2.0 * code_length / (fec_length))) + 1;

      for (std::size_t i = 0; i < output_data_length; ++i)
      {
         if ((0 == (i % error_modulus)) && (error_count < max_error_count))
         {
            output_data[i] ^= 0xFF;
            error_count++;
         }
      }

      spc_decoder.deinterleave_and_decode(output_data);
      spc_decoder.output(input_data_02);

      block_timer.stop();
      double block_time = block_timer.time();

      for (std::size_t i = 0; i < input_data_length; ++i)
      {
         if (input_data_01[i] != input_data_02[i])
         {
            std::cout << "Error: input01 and input02 at " << i << " differ." << std::endl;
            return 1;
         }
      }

      printf("Round %lu Data: %5.3fMB\tTime: %5.3fsec\tRate: %5.3fMbps\r",
             static_cast<unsigned long>(round),
             ((round + 1.0) * input_data_length) / 1048576.0,
             block_time,
             (8.0 * input_data_length) / (1048576.0 * block_time));
   }

   delete[] input_data_01;
   delete[] input_data_02;
   delete[] output_data;

   return 0;
}
