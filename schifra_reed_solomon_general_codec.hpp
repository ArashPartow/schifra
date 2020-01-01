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


#ifndef INCLUDE_SCHIFRA_REED_GENERAL_CODEC_HPP
#define INCLUDE_SCHIFRA_REED_GENERAL_CODEC_HPP


#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_ecc_traits.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t code_length, std::size_t fec_length>
      void* create_encoder(const galois::field& field,
                           const std::size_t& gen_poly_index)
      {
         const std::size_t data_length = code_length - fec_length;
         traits::validate_reed_solomon_code_parameters<code_length,fec_length,data_length>();
         galois::field_polynomial gen_polynomial(field);

         if (
              !make_sequential_root_generator_polynomial(field,
                                                         gen_poly_index,
                                                         fec_length,
                                                         gen_polynomial)
            )
         {
            return reinterpret_cast<void*>(0);
         }

         return new encoder<code_length,fec_length>(field,gen_polynomial);
      }

      template <std::size_t code_length, std::size_t fec_length>
      void* create_decoder(const galois::field& field,
                           const std::size_t& gen_poly_index)
      {
         const std::size_t data_length = code_length - fec_length;
         traits::validate_reed_solomon_code_parameters<code_length,fec_length,data_length>();
         return new decoder<code_length,fec_length>(field,static_cast<unsigned int>(gen_poly_index));
      }

      template <std::size_t code_length, std::size_t max_fec_length = 128>
      class general_codec
      {
      public:

         general_codec(const galois::field& field,
                       const std::size_t& gen_poly_index)
         {
            for (std::size_t i = 0; i < max_fec_length; ++i)
            {
               encoder_[i] = 0;
               decoder_[i] = 0;
            }

            encoder_[  2] = create_encoder<code_length,  2>(field, gen_poly_index);
            encoder_[  4] = create_encoder<code_length,  4>(field, gen_poly_index);
            encoder_[  6] = create_encoder<code_length,  6>(field, gen_poly_index);
            encoder_[  8] = create_encoder<code_length,  8>(field, gen_poly_index);
            encoder_[ 10] = create_encoder<code_length, 10>(field, gen_poly_index);
            encoder_[ 12] = create_encoder<code_length, 12>(field, gen_poly_index);
            encoder_[ 14] = create_encoder<code_length, 14>(field, gen_poly_index);
            encoder_[ 16] = create_encoder<code_length, 16>(field, gen_poly_index);
            encoder_[ 18] = create_encoder<code_length, 18>(field, gen_poly_index);
            encoder_[ 20] = create_encoder<code_length, 20>(field, gen_poly_index);
            encoder_[ 22] = create_encoder<code_length, 22>(field, gen_poly_index);
            encoder_[ 24] = create_encoder<code_length, 24>(field, gen_poly_index);
            encoder_[ 26] = create_encoder<code_length, 26>(field, gen_poly_index);
            encoder_[ 28] = create_encoder<code_length, 28>(field, gen_poly_index);
            encoder_[ 30] = create_encoder<code_length, 30>(field, gen_poly_index);
            encoder_[ 32] = create_encoder<code_length, 32>(field, gen_poly_index);
            encoder_[ 64] = create_encoder<code_length, 64>(field, gen_poly_index);
            encoder_[ 80] = create_encoder<code_length, 80>(field, gen_poly_index);
            encoder_[ 96] = create_encoder<code_length, 96>(field, gen_poly_index);
            encoder_[128] = create_encoder<code_length,128>(field, gen_poly_index);

            decoder_[  2] = create_decoder<code_length,  2>(field, gen_poly_index);
            decoder_[  4] = create_decoder<code_length,  4>(field, gen_poly_index);
            decoder_[  6] = create_decoder<code_length,  6>(field, gen_poly_index);
            decoder_[  8] = create_decoder<code_length,  8>(field, gen_poly_index);
            decoder_[ 10] = create_decoder<code_length, 10>(field, gen_poly_index);
            decoder_[ 12] = create_decoder<code_length, 12>(field, gen_poly_index);
            decoder_[ 14] = create_decoder<code_length, 14>(field, gen_poly_index);
            decoder_[ 16] = create_decoder<code_length, 16>(field, gen_poly_index);
            decoder_[ 18] = create_decoder<code_length, 18>(field, gen_poly_index);
            decoder_[ 20] = create_decoder<code_length, 20>(field, gen_poly_index);
            decoder_[ 22] = create_decoder<code_length, 22>(field, gen_poly_index);
            decoder_[ 24] = create_decoder<code_length, 24>(field, gen_poly_index);
            decoder_[ 26] = create_decoder<code_length, 26>(field, gen_poly_index);
            decoder_[ 28] = create_decoder<code_length, 28>(field, gen_poly_index);
            decoder_[ 30] = create_decoder<code_length, 30>(field, gen_poly_index);
            decoder_[ 32] = create_decoder<code_length, 32>(field, gen_poly_index);
            decoder_[ 64] = create_decoder<code_length, 64>(field, gen_poly_index);
            decoder_[ 80] = create_decoder<code_length, 80>(field, gen_poly_index);
            decoder_[ 96] = create_decoder<code_length, 96>(field, gen_poly_index);
            decoder_[128] = create_decoder<code_length,128>(field, gen_poly_index);
         }

        ~general_codec()
         {
            delete static_cast<reed_solomon::encoder<code_length,  2>*>(encoder_[  2]);
            delete static_cast<reed_solomon::encoder<code_length,  4>*>(encoder_[  4]);
            delete static_cast<reed_solomon::encoder<code_length,  6>*>(encoder_[  6]);
            delete static_cast<reed_solomon::encoder<code_length,  8>*>(encoder_[  8]);
            delete static_cast<reed_solomon::encoder<code_length, 10>*>(encoder_[ 10]);
            delete static_cast<reed_solomon::encoder<code_length, 12>*>(encoder_[ 12]);
            delete static_cast<reed_solomon::encoder<code_length, 14>*>(encoder_[ 14]);
            delete static_cast<reed_solomon::encoder<code_length, 16>*>(encoder_[ 16]);
            delete static_cast<reed_solomon::encoder<code_length, 18>*>(encoder_[ 18]);
            delete static_cast<reed_solomon::encoder<code_length, 20>*>(encoder_[ 20]);
            delete static_cast<reed_solomon::encoder<code_length, 22>*>(encoder_[ 22]);
            delete static_cast<reed_solomon::encoder<code_length, 24>*>(encoder_[ 24]);
            delete static_cast<reed_solomon::encoder<code_length, 26>*>(encoder_[ 26]);
            delete static_cast<reed_solomon::encoder<code_length, 28>*>(encoder_[ 28]);
            delete static_cast<reed_solomon::encoder<code_length, 30>*>(encoder_[ 30]);
            delete static_cast<reed_solomon::encoder<code_length, 32>*>(encoder_[ 32]);
            delete static_cast<reed_solomon::encoder<code_length, 64>*>(encoder_[ 64]);
            delete static_cast<reed_solomon::encoder<code_length, 80>*>(encoder_[ 80]);
            delete static_cast<reed_solomon::encoder<code_length, 96>*>(encoder_[ 96]);
            delete static_cast<reed_solomon::encoder<code_length,128>*>(encoder_[128]);

            delete static_cast<reed_solomon::decoder<code_length,  2>*>(decoder_[  2]);
            delete static_cast<reed_solomon::decoder<code_length,  4>*>(decoder_[  4]);
            delete static_cast<reed_solomon::decoder<code_length,  6>*>(decoder_[  6]);
            delete static_cast<reed_solomon::decoder<code_length,  8>*>(decoder_[  8]);
            delete static_cast<reed_solomon::decoder<code_length, 10>*>(decoder_[ 10]);
            delete static_cast<reed_solomon::decoder<code_length, 12>*>(decoder_[ 12]);
            delete static_cast<reed_solomon::decoder<code_length, 14>*>(decoder_[ 14]);
            delete static_cast<reed_solomon::decoder<code_length, 16>*>(decoder_[ 16]);
            delete static_cast<reed_solomon::decoder<code_length, 18>*>(decoder_[ 18]);
            delete static_cast<reed_solomon::decoder<code_length, 20>*>(decoder_[ 20]);
            delete static_cast<reed_solomon::decoder<code_length, 22>*>(decoder_[ 22]);
            delete static_cast<reed_solomon::decoder<code_length, 24>*>(decoder_[ 24]);
            delete static_cast<reed_solomon::decoder<code_length, 26>*>(decoder_[ 26]);
            delete static_cast<reed_solomon::decoder<code_length, 28>*>(decoder_[ 28]);
            delete static_cast<reed_solomon::decoder<code_length, 30>*>(decoder_[ 30]);
            delete static_cast<reed_solomon::decoder<code_length, 32>*>(decoder_[ 32]);
            delete static_cast<reed_solomon::decoder<code_length, 64>*>(decoder_[ 64]);
            delete static_cast<reed_solomon::decoder<code_length, 80>*>(decoder_[ 80]);
            delete static_cast<reed_solomon::decoder<code_length, 96>*>(decoder_[ 96]);
            delete static_cast<reed_solomon::decoder<code_length,128>*>(decoder_[128]);
         }

         template <typename Block>
         bool encode(Block& block) const
         {
            /*
               cl : code length
               fl : fec length
            */
            typedef reed_solomon::encoder<Block::trait::code_length,Block::trait::fec_length> encoder_type;
            traits::__static_assert__<(Block::trait::fec_length <= max_fec_length)>();
            if (encoder_[Block::trait::fec_length] == 0)
               return false;
            else
               return static_cast<encoder_type*>(encoder_[Block::trait::fec_length])->encode(block);
         }

         template <typename Block>
         bool decode(Block& block) const
         {
            typedef reed_solomon::decoder<Block::trait::code_length,Block::trait::fec_length> decoder_type;
            traits::__static_assert__<(Block::trait::fec_length <= max_fec_length)>();
            if (decoder_[Block::trait::fec_length] == 0)
               return false;
            else
               return static_cast<decoder_type*>(decoder_[Block::trait::fec_length])->decode(block);
         }

      private:

         void* encoder_[max_fec_length + 1];
         void* decoder_[max_fec_length + 1];
      };

   } // namespace reed_solomon

} // namespace schifra

#endif
