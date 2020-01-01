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
   Description: This example will demonstrate how one can instantiate a
                generalised codec for a specific finite field and generator
                polynomial. Also how encoding and decoding can be done without
                having to specify seperate instances of the Reed-Solomon encoders
                and decoders.
*/


#include <cstddef>
#include <iostream>
#include <vector>

#include "schifra_galois_field.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_general_codec.hpp"
#include "schifra_error_processes.hpp"
#include "schifra_utilities.hpp"


template <std::size_t code_length, std::size_t fec_length>
bool do_encode_decode(const schifra::reed_solomon::general_codec<code_length>& codec,
                            schifra::reed_solomon::block<code_length,fec_length>& block)
{
   for (int i = 0; i < static_cast<int>(code_length - fec_length); ++i) block[i] = i;

   if (!codec.encode(block))
   {
      std::cout << "do_encode_decode() - Error: Failed to encode block!" << std::endl;
      return false;
   }

   schifra::corrupt_message_all_errors(block,fec_length >> 1,0,3);

   if (!codec.decode(block))
   {
      std::cout << "do_encode_decode() - Error: Failed to decode block!" << std::endl;
      return false;
   }

   for (int i = 0; i < static_cast<int>(code_length - fec_length); ++i)
   {
      if (block[i] != i)
      {
         std::cout << "do_encode_decode() - Error: Decoding failure!" << std::endl;
      }
   }

   return true;
}

int main()
{
   /* Finite Field Parameters */
   const std::size_t field_descriptor =   8;
   const std::size_t gen_poly_index   = 120;

   /* Instantiate Finite Field and Generator Polynomials */
   const schifra::galois::field field(field_descriptor,
                                      schifra::galois::primitive_polynomial_size06,
                                      schifra::galois::primitive_polynomial06);

   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 255;

   /* Reed Solomon Code Parameters */
   schifra::reed_solomon::general_codec<code_length> codec(field, gen_poly_index);

   schifra::reed_solomon::block<code_length,  2>  block_fec02;
   schifra::reed_solomon::block<code_length,  4>  block_fec04;
   schifra::reed_solomon::block<code_length,  6>  block_fec06;
   schifra::reed_solomon::block<code_length,  8>  block_fec08;
   schifra::reed_solomon::block<code_length, 10>  block_fec10;
   schifra::reed_solomon::block<code_length, 12>  block_fec12;
   schifra::reed_solomon::block<code_length, 14>  block_fec14;
   schifra::reed_solomon::block<code_length, 16>  block_fec16;
   schifra::reed_solomon::block<code_length, 18>  block_fec18;
   schifra::reed_solomon::block<code_length, 20>  block_fec20;
   schifra::reed_solomon::block<code_length, 22>  block_fec22;
   schifra::reed_solomon::block<code_length, 24>  block_fec24;
   schifra::reed_solomon::block<code_length, 26>  block_fec26;
   schifra::reed_solomon::block<code_length, 28>  block_fec28;
   schifra::reed_solomon::block<code_length, 30>  block_fec30;
   schifra::reed_solomon::block<code_length, 32>  block_fec32;
   schifra::reed_solomon::block<code_length, 64>  block_fec64;
   schifra::reed_solomon::block<code_length, 80>  block_fec80;
   schifra::reed_solomon::block<code_length, 96>  block_fec96;
   schifra::reed_solomon::block<code_length,128> block_fec128;

   do_encode_decode(codec, block_fec02);
   do_encode_decode(codec, block_fec04);
   do_encode_decode(codec, block_fec06);
   do_encode_decode(codec, block_fec08);
   do_encode_decode(codec, block_fec10);
   do_encode_decode(codec, block_fec12);
   do_encode_decode(codec, block_fec14);
   do_encode_decode(codec, block_fec16);
   do_encode_decode(codec, block_fec18);
   do_encode_decode(codec, block_fec20);
   do_encode_decode(codec, block_fec22);
   do_encode_decode(codec, block_fec24);
   do_encode_decode(codec, block_fec26);
   do_encode_decode(codec, block_fec28);
   do_encode_decode(codec, block_fec30);
   do_encode_decode(codec, block_fec32);
   do_encode_decode(codec, block_fec64);
   do_encode_decode(codec, block_fec80);
   do_encode_decode(codec, block_fec96);
   do_encode_decode(codec,block_fec128);

   return 0;
}
