/*
(**************************************************************************)
(*                                                                        *)
(*                                Schifra                                 *)
(*                Reed-Solomon Error Correcting Code Library              *)
(*                                                                        *)
(* Release Version 0.0.1                                                  *)
(* http://www.schifra.com                                                 *)
(* Copyright (c) 2000-2017 Arash Partow, All Rights Reserved.             *)
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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_DECODER_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_DECODER_HPP


#include "schifra_galois_field.hpp"
#include "schifra_galois_field_element.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_ecc_traits.hpp"


namespace schifra
{
   namespace reed_solomon
   {

      template <std::size_t code_length, std::size_t fec_length, std::size_t data_length = code_length - fec_length>
      class decoder
      {
      public:

         typedef traits::reed_solomon_triat<code_length,fec_length,data_length> trait;
         typedef block<code_length,fec_length> block_type;

         decoder(const galois::field& gfield, const unsigned int& gen_initial_index)
         : decoder_valid_(false),
           field_(gfield),
           X_(galois::generate_X(field_)),
           gen_initial_index_(0)
         {
            if (field_.size() == code_length)
            {
               //Note: code_length and field size can be used interchangeably
               gen_initial_index_ = gen_initial_index;
               create_lookup_tables();
               decoder_valid_ = true;
            }
         };

         bool decode(block_type& rsblock) const
         {
            std::vector<std::size_t> erasure_list;
            return decode(rsblock,erasure_list);
         }

         bool decode(block_type& rsblock, const erasure_locations_t& erasure_list) const
         {
            if ((!decoder_valid_) || (erasure_list.size() > fec_length))
            {
               rsblock.errors_detected  = 0;
               rsblock.errors_corrected = 0;
               rsblock.zero_numerators  = 0;
               rsblock.unrecoverable    = true;

               return false;
            }

            galois::field_polynomial received(field_,code_length - 1);
            load_message(received,rsblock);

            galois::field_polynomial syndrome(field_);

            if (compute_syndrome(received,syndrome) == 0)
            {
               rsblock.errors_detected  = 0;
               rsblock.errors_corrected = 0;
               rsblock.zero_numerators  = 0;
               rsblock.unrecoverable    = false;

               return true;
            }

            erasure_locations_t erasure_locations;
            prepare_erasure_list(erasure_locations,erasure_list);

            galois::field_polynomial lambda(field_);

            compute_gamma(lambda,erasure_locations);

            if (erasure_list.size() < fec_length)
            {
               modified_berlekamp_massey_algorithm(lambda,syndrome,erasure_list.size());
            }

            std::vector<int> error_locations;

            find_roots(lambda,error_locations);

            if (0 == error_locations.size())
            {
               /*
                 Syndrome is non-zero yet no error locations have
                 been obtained, conclusion:
                 It is possible that there are MORE errrors in the
                 message than can be detected and corrected for this
                 particular code.
               */

               rsblock.errors_detected  = 0;
               rsblock.errors_corrected = 0;
               rsblock.zero_numerators  = 0;
               rsblock.unrecoverable    = true;

               return false;
            }
            else if (((2 * error_locations.size()) - erasure_list.size()) > fec_length)
            {
               /*
                  Too many errors\erasures! 2E + S <= fec_length
                   L =  E + S
                   E =  L - S
                  2E = 2L - 2S
                  2E + S = 2L - 2S + S
                         = 2L - S
                 Where:
                  L : Error Locations
                  E : Errors
                  S : Erasures

               */

               rsblock.errors_detected  = error_locations.size();
               rsblock.errors_corrected = 0;
               rsblock.zero_numerators  = 0;
               rsblock.unrecoverable    = true;

               return false;
            }
            else
               rsblock.errors_detected  = error_locations.size();

            return forney_algorithm(error_locations,lambda,syndrome,rsblock);
         }

      private:

         decoder();
         decoder(const decoder& dec);
         decoder& operator=(const decoder& dec);

      protected:

         void load_message(galois::field_polynomial& received, const block_type& rsblock) const
         {
            /*
              Load message data into received polynomial in reverse order.
            */

            for (std::size_t i = 0; i < code_length; ++i)
            {
               received[code_length - 1 - i] = rsblock[i];
            }
         }

         void create_lookup_tables()
         {
            root_exponent_table_.reserve(field_.size() + 1);

            for (int i = 0; i < static_cast<int>(field_.size() + 1); ++i)
            {
               root_exponent_table_.push_back(field_.exp(field_.alpha(code_length - i),(1 - gen_initial_index_)));
            }

            syndrome_exponent_table_.reserve(fec_length);

            for (int i = 0; i < static_cast<int>(fec_length); ++i)
            {
               syndrome_exponent_table_.push_back(field_.alpha(gen_initial_index_ + i));
            }

            gamma_table_.reserve(field_.size() + 1);

            for (int i = 0; i < static_cast<int>(field_.size() + 1); ++i)
            {
               gamma_table_.push_back((1 + (X_ * galois::field_element(field_,field_.alpha(i)))));
            }
         }

         void prepare_erasure_list(erasure_locations_t& erasure_locations, const erasure_locations_t& erasure_list) const
         {
            /*
              Note: 1. Erasure positions must be unique.
                    2. Erasure positions must exist within the code block.
                    There are NO exceptions to these rules!
            */

            erasure_locations.resize(erasure_list.size());

            for (std::size_t i = 0; i < erasure_list.size(); ++i)
            {
               erasure_locations[i] = (code_length - 1 - erasure_list[i]);
            }
         }

         int compute_syndrome(const galois::field_polynomial& received,
                                    galois::field_polynomial& syndrome) const
         {
            int error_flag = 0;
            syndrome = galois::field_polynomial(field_,fec_length - 1);

            for (std::size_t i = 0; i < fec_length; ++i)
            {
               syndrome[i]  = received(syndrome_exponent_table_[i]);
               error_flag  |= syndrome[i].poly();
            }

            return error_flag;
         }

         void compute_gamma(galois::field_polynomial& gamma, const erasure_locations_t& erasure_locations) const
         {
            gamma = galois::field_element(field_,1);

            for (std::size_t i = 0; i < erasure_locations.size(); ++i)
            {
               gamma *= gamma_table_[erasure_locations[i]];
            }
         }

         void find_roots(const galois::field_polynomial& poly, std::vector<int>& root_list) const
         {
            /*
               Chien Search: Find the roots of the error locator polynomial
               via an exhaustive search over all non-zero elements in the
               given finite field.
            */

            root_list.reserve(fec_length << 1);
            root_list.resize(0);

            std::size_t polynomial_degree = poly.deg();
            std::size_t root_list_size = 0;

            for (int i = 1; i <= static_cast<int>(code_length); ++i)
            {
               if (0 == poly(field_.alpha(i)).poly())
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

         void compute_discrepancy(galois::field_element&          discrepancy,
                                  const galois::field_polynomial& lambda,
                                  const galois::field_polynomial& syndrome,
                                  const std::size_t&              l,
                                  const std::size_t&              round) const
         {
            /*
               Compute the lambda discrepancy at the current round of BMA
            */

            int upper_bound = std::min(static_cast<int>(l),lambda.deg());
            discrepancy = 0;

            for (int i = 0; i <= upper_bound; ++i)
            {
               discrepancy += lambda[i] * syndrome[round - i];
            }
         }

         void modified_berlekamp_massey_algorithm(galois::field_polynomial&       lambda,
                                                  const galois::field_polynomial& syndrome,
                                                  const std::size_t               erasure_count) const
         {
            /*
               Modified Berlekamp-Massey Algorithm
               Identify the shortest length linear feed-back shift register (LFSR)
               that will generate the sequence equivalent to the syndrome.
            */

            int i = -1;
            std::size_t l = erasure_count;

            galois::field_element discrepancy(field_,0);
            galois::field_polynomial previous_lambda = lambda << 1;

            for (std::size_t round = erasure_count; round < fec_length; ++round)
            {
               compute_discrepancy(discrepancy,lambda,syndrome,l,round);

               if (discrepancy != 0)
               {
                  galois::field_polynomial tau = lambda - discrepancy * previous_lambda;

                  if (static_cast<int>(l) < (static_cast<int>(round) - i))
                  {
                     std::size_t tmp = round - i;
                     i = static_cast<int>(round - l);
                     l = tmp;
                     previous_lambda = lambda / discrepancy;
                  }

                  lambda = tau;
               }

               previous_lambda <<= 1;
            }
         }

         bool forney_algorithm(const std::vector<int>&         error_locations,
                               const galois::field_polynomial& lambda,
                               const galois::field_polynomial& syndrome,
                               block_type&                     rsblock) const
         {
            /*
               The Forney algorithm for computing the error magnitudes
            */
            galois::field_polynomial omega = (lambda * syndrome) % fec_length;
            galois::field_polynomial lambda_derivative = lambda.derivative();

            rsblock.errors_corrected = 0;
            rsblock.zero_numerators  = 0;

            galois::field_element numerator(field_);
            galois::field_element denominator(field_);

            for (std::size_t i = 0; i < error_locations.size(); ++i)
            {
               int error_location                  = error_locations[i];
               galois::field_symbol  alpha_inverse = field_.alpha(error_location);
               numerator                           = (omega(alpha_inverse) * root_exponent_table_[error_location]);
               denominator                         = lambda_derivative(alpha_inverse);

               if (numerator != 0)
               {
                  if (denominator != 0)
                  {
                     rsblock[error_location - 1] ^= field_.div(numerator.poly(),denominator.poly());
                     rsblock.errors_corrected++;
                  }
                  else
                  {
                     rsblock.unrecoverable = true;
                     return false;
                  }
               }
               else
                  ++rsblock.zero_numerators;
            }

            if (lambda.deg() == static_cast<int>(rsblock.errors_detected))
               return true;
            else
            {
               rsblock.unrecoverable = true;
               return false;
            }
         }

      protected:

         bool                                  decoder_valid_;
         const galois::field&                  field_;
         std::vector<galois::field_symbol>     root_exponent_table_;
         std::vector<galois::field_symbol>     syndrome_exponent_table_;
         std::vector<galois::field_polynomial> gamma_table_;
         galois::field_polynomial              X_;
         unsigned int                          gen_initial_index_;
      };

      template <std::size_t code_length,
                std::size_t fec_length,
                std::size_t data_length    = code_length - fec_length,
                std::size_t natural_length = 255,  // Needs to be in-sync with field size
                std::size_t padding_length = natural_length - data_length - fec_length >
      class shortened_decoder
      {
      public:

         typedef traits::reed_solomon_triat<code_length,fec_length,data_length> trait;
         typedef block<code_length,fec_length> block_type;

         shortened_decoder(const galois::field& field, const unsigned int gen_initial_index)
         : decoder_(field,gen_initial_index)
         {
            for (std::size_t i = 0; i < natural_length; ++i)
            {
               block_[i] = 0;
            }
         }

         inline bool decode(block_type& rsblock, const erasure_locations_t& erasure_list)
         {
            for (std::size_t i = 0; i < code_length; ++i)
            {
               block_.data[padding_length + i] = rsblock.data[i];
            }

            erasure_locations_t shifted_position_erasure_list(erasure_list.size(),0);

            for (std::size_t i = 0; i < erasure_list.size(); ++i)
            {
               shifted_position_erasure_list[i] = erasure_list[i] + padding_length;
            }

            if (decoder_.decode(block_,shifted_position_erasure_list))
            {
               for (std::size_t i = 0; i < code_length; ++i)
               {
                  rsblock.data[i] = block_.data[padding_length + i];
               }

               rsblock.errors_detected  = block_.errors_detected;
               rsblock.errors_corrected = block_.errors_corrected;
               rsblock.zero_numerators  = block_.zero_numerators;
               rsblock.unrecoverable    = block_.unrecoverable;

               return true;
            }
            else
            {
               rsblock.errors_detected  = block_.errors_detected;
               rsblock.errors_corrected = block_.errors_corrected;
               rsblock.zero_numerators  = block_.zero_numerators;
               rsblock.unrecoverable    = block_.unrecoverable;

               return false;
            }
         }

         inline bool decode(block_type& rsblock)
         {
            for (std::size_t i = 0; i < code_length; ++i)
            {
               block_.data[padding_length + i] = rsblock.data[i];
            }

            if (decoder_.decode(block_))
            {
               for (std::size_t i = 0; i < code_length; ++i)
               {
                  rsblock.data[i] = block_.data[padding_length + i];
               }

               rsblock.errors_detected  = block_.errors_detected;
               rsblock.errors_corrected = block_.errors_corrected;
               rsblock.zero_numerators  = block_.zero_numerators;
               rsblock.unrecoverable    = block_.unrecoverable;

               return true;
            }
            else
            {
               rsblock.errors_detected  = block_.errors_detected;
               rsblock.errors_corrected = block_.errors_corrected;
               rsblock.zero_numerators  = block_.zero_numerators;
               rsblock.unrecoverable    = block_.unrecoverable;

               return false;
            }
         }

      private:

         typedef decoder<natural_length,fec_length> natural_decoder_type;
         natural_decoder_type decoder_;
         typename natural_decoder_type::block_type block_;
      };

   } // namespace reed_solomon

} // namespace schifra

#endif
