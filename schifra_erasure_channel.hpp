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


#ifndef INCLUDE_SCHIFRA_ERASURE_CHANNEL_HPP
#define INCLUDE_SCHIFRA_ERASURE_CHANNEL_HPP


#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_interleaving.hpp"
#include "schifra_utilities.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t block_length, std::size_t fec_length>
      inline void interleaved_stack_erasure_mapper(const std::vector<std::size_t>& missing_row_index,
                                                         std::vector<erasure_locations_t>& erasure_row_list)
      {
         erasure_row_list.resize(block_length);

         for (std::size_t i = 0; i < block_length; ++i)
         {
            erasure_row_list[i].reserve(fec_length);
         }

         for (std::size_t i = 0; i < missing_row_index.size(); ++i)
         {
            for (std::size_t j = 0; j < block_length; ++j)
            {
               erasure_row_list[j].push_back(missing_row_index[i]);
            }
         }
      }

      template <std::size_t code_length, std::size_t fec_length>
      inline bool erasure_channel_stack_encode(const encoder<code_length,fec_length>& encoder,
                                                     block<code_length,fec_length> (&output)[code_length])
      {
         for (std::size_t i = 0; i < code_length; ++i)
         {
            if (!encoder.encode(output[i]))
            {
               std::cout << "erasure_channel_stack_encode() - Error: Failed to encode block[" << i <<"]" << std::endl;

               return false;
            }
         }

         interleave<code_length,fec_length>(output);

         return true;
      }

      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length = code_length - fec_length>
      class erasure_code_decoder : public decoder<code_length,fec_length,data_length>
      {
      public:

         typedef decoder<code_length,fec_length,data_length> decoder_type;
         typedef typename decoder_type::block_type block_type;
         typedef std::vector<galois::field_polynomial> polynomial_list_type;

         erasure_code_decoder(const galois::field& gfield,
                              const unsigned int& gen_initial_index)
         : decoder<code_length,fec_length,data_length>(gfield, gen_initial_index)
         {
            for (std::size_t i = 0; i < code_length; ++i)
            {
               received_.push_back(galois::field_polynomial(decoder_type::field_, code_length - 1));
               syndrome_.push_back(galois::field_polynomial(decoder_type::field_));
            }
         };

         bool decode(block_type rsblock[code_length], const erasure_locations_t& erasure_list) const
         {
            if (
                 (!decoder_type::decoder_valid_) ||
                 (erasure_list.size() != fec_length)
               )
            {
               return false;
            }

            for (std::size_t i = 0; i < code_length; ++i)
            {
               decoder_type::load_message    (received_[i], rsblock  [i]);
               decoder_type::compute_syndrome(received_[i], syndrome_[i]);
            }

            erasure_locations_t erasure_locations;
            decoder_type::prepare_erasure_list(erasure_locations,erasure_list);

            galois::field_polynomial gamma(galois::field_element(decoder_type::field_, 1));

            decoder_type::compute_gamma(gamma,erasure_locations);

            std::vector<int> gamma_roots;

            find_roots_in_data(gamma,gamma_roots);

            polynomial_list_type omega;

            for (std::size_t i = 0; i < code_length; ++i)
            {
               omega.push_back((gamma * syndrome_[i]) % fec_length);
            }

            galois::field_polynomial gamma_derivative = gamma.derivative();

            for (std::size_t i = 0; i < gamma_roots.size(); ++i)
            {
               int error_location                  = static_cast<int>(gamma_roots[i]);
               galois::field_symbol  alpha_inverse = decoder_type::field_.alpha(error_location);
               galois::field_element denominator   = gamma_derivative(alpha_inverse);

               if (denominator == 0)
               {
                  return false;
               }

               for (std::size_t j = 0; j < code_length; ++j)
               {
                  galois::field_element numerator = (omega[j](alpha_inverse) * decoder_type::root_exponent_table_[error_location]);
                  /*
                    A minor optimization can be made in the event the
                    numerator is equal to zero by not executing the
                    following line.
                  */
                  rsblock[j][error_location - 1] ^= decoder_type::field_.div(numerator.poly(),denominator.poly());
               }
            }

            return true;
         }

      private:

         void find_roots_in_data(const galois::field_polynomial& poly, std::vector<int>& root_list) const
         {
            /*
               Chien Search, as described in parent, but only
               for locations within the data range of the message.
            */
            root_list.reserve(fec_length << 1);
            root_list.resize(0);

            std::size_t polynomial_degree = poly.deg();
            std::size_t root_list_size = 0;

            for (int i = 1; i <= static_cast<int>(data_length); ++i)
            {
               if (0 == poly(decoder_type::field_.alpha(i)).poly())
               {
                  root_list.push_back(i);
                  root_list_size++;

                  if (root_list_size == polynomial_degree)
                  {
                     break;
                  }
               }
            }
         }

         mutable polynomial_list_type received_;
         mutable polynomial_list_type syndrome_;

      };

      template <std::size_t code_length, std::size_t fec_length>
      inline bool erasure_channel_stack_decode(const decoder<code_length,fec_length>& general_decoder,
                                               const erasure_locations_t& missing_row_index,
                                                     block<code_length,fec_length> (&output)[code_length])
      {
         if (missing_row_index.empty())
         {
            return true;
         }

         interleave<code_length,fec_length>(output);

         for (std::size_t i = 0; i < code_length; ++i)
         {
            if (!general_decoder.decode(output[i],missing_row_index))
            {
               std::cout << "[2] erasure_channel_stack_decode() - Error: Failed to decode block[" << i <<"]" << std::endl;

               return false;
            }
         }

         return true;
      }

      template <std::size_t code_length, std::size_t fec_length>
      inline bool erasure_channel_stack_decode(const erasure_code_decoder<code_length,fec_length>& erasure_decoder,
                                               const erasure_locations_t& missing_row_index,
                                                     block<code_length,fec_length> (&output)[code_length])
      {
         /*
           Note: 1. Missing row indicies must be unique.
                 2. Missing row indicies must exist within
                    the stack's size.
                 3. There will be NO errors in the rows (aka output)
                 4. The information members of the blocks will
                    not be utilized.
                 There are NO exceptions to these rules!
         */
         if (missing_row_index.empty())
         {
            return true;
         }
         else if (missing_row_index.size() == fec_length)
         {
            interleave<code_length,fec_length>(output);

            return erasure_decoder.decode(output,missing_row_index);
         }
         else
            return erasure_channel_stack_decode<code_length,fec_length>(
                      static_cast<const decoder<code_length,fec_length>&>(erasure_decoder),
                      missing_row_index,
                      output);
      }

   } // namespace reed_solomon

} // namepsace schifra


#endif
