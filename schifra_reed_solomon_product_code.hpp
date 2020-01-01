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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_PRODUCT_CODE_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_PRODUCT_CODE_HPP


#include <cstddef>
#include <iostream>
#include <fstream>

#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_interleaving.hpp"
#include "schifra_reed_solomon_bitio.hpp"
#include "schifra_ecc_traits.hpp"


namespace schifra
{

   namespace reed_solomon
   {
      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length = code_length - fec_length>
      class square_product_code_encoder
      {
      public:

         typedef encoder<code_length,fec_length> encoder_type;
         typedef block<code_length,fec_length> block_type;
         typedef traits::reed_solomon_triat<code_length,fec_length,data_length> trait;
         typedef unsigned char data_type;
         typedef data_type* data_ptr_type;

         enum { data_size  = data_length * data_length };
         enum { total_size = code_length * code_length };

         square_product_code_encoder(const encoder_type& enc)
         : encoder_(enc)
         {}

         bool encode(data_ptr_type data)
         {
            data_ptr_type curr_data_ptr = data;

            for (std::size_t row = 0; row < data_length; ++row, curr_data_ptr += data_length)
            {
               copy(curr_data_ptr, data_length, block_stack_[row]);

               if (!encoder_.encode(block_stack_[row]))
               {
                  return false;
               }
            }

            block_type vertical_block;

            for (std::size_t col = 0; col < code_length; ++col)
            {
               for (std::size_t row = 0; row < data_length; ++row)
               {
                  vertical_block[row] = block_stack_[row][col];
               }

               if (!encoder_.encode(vertical_block))
               {
                  return false;
               }

               for (std::size_t fec_index = 0; fec_index < fec_length; ++fec_index)
               {
                  block_stack_[data_length + fec_index].fec(fec_index) = vertical_block.fec(fec_index);
               }
            }

            return true;
         }

         bool encode_and_interleave(data_ptr_type data)
         {
            if (!encode(data))
            {
               return false;
            }

            interleave<code_length,fec_length>(block_stack_);

            return true;
         }

         void output(data_ptr_type output_data)
         {
            for (std::size_t row = 0; row < code_length; ++row, output_data += code_length)
            {
               bitio::convert_symbol_to_data<traits::symbol<code_length>::size>(block_stack_[row].data,output_data,code_length);
            }
         }

         void clear()
         {
            for (std::size_t i = 0; i < code_length; ++i)
            {
               block_stack_[i].clear();
            }
         }

      private:

         square_product_code_encoder(const square_product_code_encoder& spce);
         square_product_code_encoder& operator=(const square_product_code_encoder& spce);

         block_type block_stack_[code_length];
         const encoder_type& encoder_;
      };

      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length = code_length - fec_length>
      class square_product_code_decoder
      {
      public:

         typedef decoder<code_length,fec_length> decoder_type;
         typedef block<code_length,fec_length> block_type;
         typedef traits::reed_solomon_triat<code_length,fec_length,data_length> trait;
         typedef unsigned char data_type;
         typedef data_type* data_ptr_type;

         enum { data_size  = data_length * data_length };
         enum { total_size = code_length * code_length };

         square_product_code_decoder(const decoder_type& decoder)
         : decoder_(decoder)
         {}

         void decode(data_ptr_type data)
         {
            copy_proxy(data);
            decode_proxy();
         }

         void deinterleave_and_decode(data_ptr_type data)
         {
            copy_proxy(data);
            interleave<code_length,fec_length>(block_stack_);
            decode_proxy();
         }

         void output(data_ptr_type output_data)
         {
            for (std::size_t row = 0; row < data_length; ++row, output_data += data_length)
            {
               bitio::convert_symbol_to_data<traits::symbol<code_length>::size>(block_stack_[row].data,output_data,data_length);
            }
         }

         void clear()
         {
            for (std::size_t i = 0; i < code_length; ++i)
            {
               block_stack_[i].clear();
            }
         }

      private:

         square_product_code_decoder(const square_product_code_decoder& spcd);
         square_product_code_decoder& operator=(const square_product_code_decoder& spcd);

         void copy_proxy(data_ptr_type data)
         {
            for (std::size_t row = 0; row < code_length; ++row, data += code_length)
            {
               bitio::convert_data_to_symbol<traits::symbol<code_length>::size>(data,code_length,block_stack_[row].data);
            }
         }

         void decode_proxy()
         {
            bool first_iteration_failure = false;

            for (std::size_t row = 0; row < data_length; ++row)
            {
               if (!decoder_.decode(block_stack_[row]))
               {
                  first_iteration_failure = true;
               }
            }

            if (!first_iteration_failure)
            {
               /*
                 Either no errors detected or all errors have
                 been detected and corrected.
               */
               return;
            }

            block_type vertical_block;

            for (std::size_t col = 0; col < code_length; ++col)
            {
               for (std::size_t row = 0; row < data_length; ++row)
               {
                  vertical_block[row] = block_stack_[row][col];
               }

               decoder_.decode(vertical_block);
            }
         }

         block_type block_stack_[code_length];
         const decoder_type& decoder_;
      };

   } // namespace reed_solomon

} // namespace schifra

#endif
