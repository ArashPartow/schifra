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
#include <iostream>
#include <string>

#include "schifra_galois_field.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_file_encoder.hpp"
#include "schifra_reed_solomon_file_decoder.hpp"
#include "schifra_reed_solomon_file_interleaver.hpp"
#include "schifra_fileio.hpp"
#include "schifra_error_processes.hpp"


/*
   Note:
   This example demonstrates how integrating RS ECC and interleaving on a large
   enough block of data allows for the correction of burst errors of length many
   times longer than what can normally be corrected. In this case 255 times longer.
   The example will do the following:
   1. Create a file
   2. Reed-Solomon encode the file with 2 symbol error detecting 1 symbol error correcting capability
   3. Interleave the file
   4. Corrupt the file with a burst error of code length
   5. Deinterleave the file
   6. Reed-Solomon decode the file attempting to correct all encountered errors
   7. Compare the original file to the final output file
*/

void create_file(const std::string& file_name, const std::size_t file_size)
{
   std::string buffer = std::string(file_size,0x00);

   for (std::size_t i = 0; i < buffer.size(); ++i)
   {
      buffer[i] = static_cast<unsigned char>(i & 0xFF);
   }

   schifra::fileio::write_file(file_name,buffer);
}

int main()
{
   const std::size_t code_length         = 255;
   const std::size_t fec_length          =   2;
   const std::size_t data_length         = code_length - fec_length;
   const std::size_t field_descriptor    =   8;
   const std::size_t gen_poly_index      = 120;
   const std::size_t gen_poly_root_count = fec_length;
   const std::size_t stack_size          = 255;

   const std::string input_file_name                = "input.dat";
   const std::string rsencoded_output_file_name     = "output.rsenc";
   const std::string interleaved_output_file_name   = "output.intr";
   const std::string deinterleaved_output_file_name = "output.deintr";
   const std::string rsdecoded_file_name            = "output.rsdec";

   const schifra::galois::field field(field_descriptor,
                                      schifra::galois::primitive_polynomial_size06,
                                      schifra::galois::primitive_polynomial06);

   schifra::galois::field_polynomial generator_polynomial(field);

   if (
        !schifra::make_sequential_root_generator_polynomial(field,
                                                            gen_poly_index,
                                                            gen_poly_root_count,
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
   const decoder_t decoder(field, gen_poly_index);

   create_file(input_file_name,data_length * stack_size - 3);

   schifra::reed_solomon::file_encoder<code_length,fec_length>
                          (
                            encoder,
                            input_file_name,
                            rsencoded_output_file_name
                          );

   schifra::reed_solomon::file_interleaver<code_length,stack_size>
                          (
                            rsencoded_output_file_name,
                            interleaved_output_file_name
                          );

   schifra::corrupt_file_with_burst_errors
            (
              interleaved_output_file_name,
              10,
              code_length * (fec_length >> 1)
            );

   schifra::reed_solomon::file_deinterleaver<code_length,stack_size>
                          (
                            interleaved_output_file_name,
                            deinterleaved_output_file_name
                          );

   schifra::reed_solomon::file_decoder<code_length,fec_length>
                          (
                            decoder,
                            deinterleaved_output_file_name,
                            rsdecoded_file_name
                          );

   if (!schifra::fileio::files_identical(input_file_name, rsdecoded_file_name))
   {
      std::cout << "ERROR - Input file and decoded deinterleaved files are not equivelent!" << std::endl;
      return 1;
   }

   return 0;
}
