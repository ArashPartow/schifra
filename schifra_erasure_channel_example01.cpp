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
   Description: This example will demonstrate how to implement a simple Reed-Solomon based
                channel code for Erasure Channels. This is done by creating a square stack
                of Reed-Solomon blocks, encoding the blocks, then interleaving the square
                stack by transposing it. Whole rows are then erased to simulate the packet
                loss that might occur over an erasure channel such as UDP or some form of
                wireless data communications. Finally a mapping is calculated in order to
                ascertain the erasure locations for the deinterleaved stack rows based on
                the "missing" rows in the interleaved stack. Once the incomplete interleaved
                stack has been dinterleaved, each of the rows or in this case the Reed-Solomon
                blocks are decoded. The decoding process is coupled with the designated
                erasures from the previous steps.

                Possible relevant use-cases for the functionality of this demonstration are:
                1. Overcoming UDP packet loss
                2. Data reconstruction from distributed data storage systems
*/


#include <cstddef>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_erasure_channel.hpp"
#include "schifra_utilities.hpp"


int main()
{
   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 255;
   const std::size_t fec_length  =  20;
   const std::size_t data_length = code_length - fec_length;
   const std::size_t stack_size  = code_length;
   const std::size_t data_size   = stack_size * data_length;

   /* Finite Field Parameters */
   const std::size_t field_descriptor                =   8;
   const std::size_t generator_polynomial_index      = 120;
   const std::size_t generator_polynomial_root_count = fec_length;

   /* Instantiate Finite Field and Generator Polynomials */
   schifra::galois::field field(field_descriptor,
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

   encoder_t rs_encoder(field,generator_polynomial);
   decoder_t rs_decoder(field,generator_polynomial_index);

   schifra::reed_solomon::erasure_code_decoder<code_length,fec_length> rs_erasure_decoder(field,generator_polynomial_index);

   schifra::reed_solomon::block<code_length,fec_length> block_stack[stack_size];

   schifra::utils::timer timer;
   timer.start();

   unsigned char send_data[data_size];
   unsigned char recv_data[data_size];

   const std::size_t max_iterations = 100000;

   for (std::size_t iteration = 0; iteration < max_iterations; ++iteration)
   {
      /* Populate block stack with data */
      for (std::size_t i = 0; i < data_size; ++i)
      {
         send_data[i] = static_cast<unsigned char>((i * 3 + 7 * iteration) & 0xFF);
      }

      schifra::utils::timer block_timer;
      block_timer.start();

      schifra::reed_solomon::copy<unsigned char,code_length,fec_length,stack_size>(send_data,data_size,block_stack);

      schifra::reed_solomon::erasure_channel_stack_encode<code_length,fec_length>(rs_encoder,block_stack);

      /* Add Erasures - Simulate network packet loss (e.g: UDP) */
      schifra::reed_solomon::erasure_locations_t missing_row_index;
      missing_row_index.clear();

      for (std::size_t i = 0; i < fec_length; ++i)
      {
         std::size_t missing_index = (iteration + (i * 3)) % stack_size;
         block_stack[missing_index].clear();
         missing_row_index.push_back(missing_index);
      }

      schifra::reed_solomon::erasure_channel_stack_decode<code_length,fec_length>(rs_decoder,missing_row_index,block_stack);

      schifra::reed_solomon::copy<unsigned char,code_length,fec_length,stack_size>(block_stack,recv_data);

      block_timer.stop();

      double block_time = block_timer.time();

      for (std::size_t i = 0; i < data_size; ++i)
      {
         if (recv_data[i] != send_data[i])
         {
            std::cout << std::endl << std::endl;
            std::cout << "Error: Final block stack comparison failed! stack: " << i << std::endl;

            return 1;
         }
      }

      printf("Round %5lu Blocks Decoded:%7lu Data:%7.3fMB\tTime:%7.3fsec\tRate:%7.3fMbps\tBRate:%7.3fBlk/sec\r",
             static_cast<unsigned long>(iteration),
             static_cast<unsigned long>((iteration + 1) * stack_size),
             (1.0 * iteration * stack_size * data_length) / 1048576.0,
             block_time,
             (8.0 * stack_size * data_length) / (1048576.0 * block_time),
             (1.0 * stack_size) / block_time);
   }

   timer.stop();
   double time = timer.time();

   std::cout << "Blocks Decoded: " << max_iterations * stack_size                                   << "\t"
                "Data: "           << (1.0 * max_iterations * stack_size * data_length) / 1048576.0 << "MB\t"
                "Time: "           << time                                                          << "sec\t"
                "Rate: "           << (8.0 * max_iterations * stack_size * data_length) / (1048576.0 * time)  << "Mbps" << std::endl;

   return 0;
}
