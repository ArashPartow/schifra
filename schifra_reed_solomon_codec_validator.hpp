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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_CODEC_VALIDATOR_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_CODEC_VALIDATOR_HPP


#include <cstddef>
#include <iostream>
#include <string>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_ecc_traits.hpp"
#include "schifra_error_processes.hpp"
#include "schifra_utilities.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t code_length,
                std::size_t fec_length,
                typename encoder_type = encoder<code_length,fec_length>,
                typename decoder_type = decoder<code_length,fec_length>,
                std::size_t data_length = code_length - fec_length>
      class codec_validator
      {
      public:

         typedef block<code_length,fec_length> block_type;

         codec_validator(const galois::field& gf,
                         const unsigned int gpii,
                         const std::string& msg)
         : field_(gf),
           generator_polynomial_(galois::field_polynomial(field_)),
           rs_encoder_(reinterpret_cast<encoder_type*>(0)),
           rs_decoder_(reinterpret_cast<decoder_type*>(0)),
           message(msg),
           genpoly_initial_index_(gpii),
           blocks_processed_(0),
           block_failures_(0)
         {
            traits::equivalent_encoder_decoder<encoder_type,decoder_type>();

            if (
                 !make_sequential_root_generator_polynomial(field_,
                                                            genpoly_initial_index_,
                                                            fec_length,
                                                            generator_polynomial_)
               )
            {
               return;
            }

            rs_encoder_ = new encoder_type(field_,generator_polynomial_);
            rs_decoder_ = new decoder_type(field_,genpoly_initial_index_);

            if (!rs_encoder_->encode(message,rs_block_original))
            {
               std::cout << "codec_validator() - ERROR: Encoding process failed!" << std::endl;
               return;
            }
         }

         bool execute()
         {
            schifra::utils::timer timer;
            timer.start();

            bool result = stage1() &&
                          stage2() &&
                          stage3() &&
                          stage4() &&
                          stage5() &&
                          stage6() &&
                          stage7() &&
                          stage8() &&
                          stage9() &&
                         stage10() &&
                         stage11() &&
                         stage12() ;

            timer.stop();

            double time = timer.time();

            print_codec_properties();
            std::cout << "Blocks decoded: "       << blocks_processed_ <<
                         "\tDecoding Failures: "  << block_failures_   <<
                         "\tRate: "               << ((blocks_processed_ * data_length) * 8.0) / (1048576.0 * time) << "Mbps" << std::endl;
            /*
              Note: The throughput rate is not only the throughput of reed solomon
                    encoding and decoding, but also that of the steps needed to add
                    simulated transmission errors to the reed solomon block such as
                    the calculation of the positions and additions of errors and
                    erasures to the reed solomon block, which normally in a true
                    data transmission medium would not be taken into consideration.
            */
            return result;
         }

        ~codec_validator()
         {
            delete rs_encoder_;
            delete rs_decoder_;
         }

         void print_codec_properties()
         {
            std::cout << "Codec: RS(" << code_length << "," << data_length << "," << fec_length <<") ";
         }

      private:

         bool stage1()
         {
            /* Burst Error Only Combinations */

            const std::size_t initial_failure_count = block_failures_;

            for (std::size_t error_count = 1; error_count <= (fec_length >> 1); ++error_count)
            {
               for (std::size_t start_position = 0; start_position < code_length; ++start_position)
               {
                  block_type rs_block = rs_block_original;

                  corrupt_message_all_errors
                  (
                    rs_block,
                    error_count,
                    start_position,
                    1
                  );

                  if (!rs_decoder_->decode(rs_block))
                  {
                     print_codec_properties();
                     std::cout << "stage1() - Decoding Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (!is_block_equivelent(rs_block,message))
                  {
                     print_codec_properties();
                     std::cout << "stage1() - Error Correcting Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != rs_block.errors_corrected)
                  {
                     print_codec_properties();
                     std::cout << "stage1() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != error_count)
                  {
                     print_codec_properties();
                     std::cout << "stage1() - Error In The Number Of Detected Errors! Errors Detected: " << rs_block.errors_detected << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_corrected != error_count)
                  {
                     print_codec_properties();
                     std::cout << "stage1() - Error In The Number Of Corrected Errors! Errors Corrected: " << rs_block.errors_corrected << std::endl;
                     ++block_failures_;
                  }

                  ++blocks_processed_;
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage2()
         {
            /* Burst Erasure Only Combinations */

            const std::size_t initial_failure_count = block_failures_;

            erasure_locations_t erasure_list;

            for (std::size_t erasure_count = 1; erasure_count <= fec_length; ++erasure_count)
            {
               for (std::size_t start_position = 0; start_position < code_length; ++start_position)
               {
                  block_type rs_block = rs_block_original;

                  corrupt_message_all_erasures
                  (
                    rs_block,
                    erasure_list,
                    erasure_count,
                    start_position,
                    1
                  );

                  if (!rs_decoder_->decode(rs_block,erasure_list))
                  {
                     print_codec_properties();
                     std::cout << "stage2() - Decoding Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (!is_block_equivelent(rs_block,message))
                  {
                     std::cout << "stage2() - Error Correcting Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != rs_block.errors_corrected)
                  {
                     print_codec_properties();
                     std::cout << "stage2() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != erasure_count)
                  {
                     print_codec_properties();
                     std::cout << "stage2() - Error In The Number Of Detected Errors! Errors Detected: " << rs_block.errors_detected << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_corrected != erasure_count)
                  {
                     print_codec_properties();
                     std::cout << "stage2() - Error In The Number Of Corrected Errors! Errors Corrected: " << rs_block.errors_corrected << std::endl;
                     ++block_failures_;
                  }

                  ++blocks_processed_;
                  erasure_list.clear();
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage3()
         {
            /* Consecutive Burst Erasure and Error Combinations */

            const std::size_t initial_failure_count = block_failures_;

            erasure_locations_t erasure_list;

            for (std::size_t erasure_count = 1; erasure_count <= fec_length; ++erasure_count)
            {
               for (std::size_t start_position = 0; start_position < code_length; ++start_position)
               {
                  block_type rs_block = rs_block_original;

                  corrupt_message_errors_erasures
                  (
                    rs_block,
                    error_mode::erasures_errors,
                    start_position,erasure_count,
                    erasure_list
                  );

                  if (!rs_decoder_->decode(rs_block,erasure_list))
                  {
                     print_codec_properties();
                     std::cout << "stage3() - Decoding Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (!is_block_equivelent(rs_block,message))
                  {
                     print_codec_properties();
                     std::cout << "stage3() - Error Correcting Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != rs_block.errors_corrected)
                  {
                     print_codec_properties();
                     std::cout << "stage3() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                     ++block_failures_;
                  }

                  ++blocks_processed_;
                  erasure_list.clear();
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage4()
         {
            /* Consecutive Burst Error and Erasure Combinations */

            const std::size_t initial_failure_count = block_failures_;

            erasure_locations_t erasure_list;

            for (std::size_t erasure_count = 1; erasure_count <= fec_length; ++erasure_count)
            {
               for (std::size_t start_position = 0; start_position < code_length; ++start_position)
               {
                  block_type rs_block = rs_block_original;

                  corrupt_message_errors_erasures
                  (
                    rs_block,
                    error_mode::errors_erasures,
                    start_position,
                    erasure_count,
                    erasure_list
                  );

                  if (!rs_decoder_->decode(rs_block,erasure_list))
                  {
                     print_codec_properties();
                     std::cout << "stage4() - Decoding Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (!is_block_equivelent(rs_block,message))
                  {
                     print_codec_properties();
                     std::cout << "stage4() - Error Correcting Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != rs_block.errors_corrected)
                  {
                     print_codec_properties();
                     std::cout << "stage4() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                     ++block_failures_;
                  }

                  ++blocks_processed_;
                  erasure_list.clear();
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage5()
         {
            /* Distanced Burst Erasure and Error Combinations */

            const std::size_t initial_failure_count = block_failures_;

            erasure_locations_t erasure_list;

            for (std::size_t between_distance = 1; between_distance <= 10; ++between_distance)
            {
               for (std::size_t erasure_count = 1; erasure_count <= fec_length; ++erasure_count)
               {
                  for (std::size_t start_position = 0; start_position < code_length; ++start_position)
                  {
                     block_type rs_block = rs_block_original;

                     corrupt_message_errors_erasures
                     (
                       rs_block,
                       error_mode::erasures_errors,
                       start_position,
                       erasure_count,
                       erasure_list,
                       between_distance
                     );

                     if (!rs_decoder_->decode(rs_block,erasure_list))
                     {
                        print_codec_properties();
                        std::cout << "stage5() - Decoding Failure! start position: " << start_position << std::endl;
                        ++block_failures_;
                     }
                     else if (!is_block_equivelent(rs_block,message))
                     {
                        print_codec_properties();
                        std::cout << "stage5() - Error Correcting Failure! start position: " << start_position << std::endl;
                        ++block_failures_;
                     }
                     else if (rs_block.errors_detected != rs_block.errors_corrected)
                     {
                        print_codec_properties();
                        std::cout << "stage5() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                        ++block_failures_;
                     }

                     ++blocks_processed_;
                     erasure_list.clear();
                  }
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage6()
         {
            /* Distanced Burst Error and Erasure Combinations */

            const std::size_t initial_failure_count = block_failures_;

            erasure_locations_t erasure_list;

            for (std::size_t between_distance = 1; between_distance <= 10; ++between_distance)
            {
               for (std::size_t erasure_count = 1; erasure_count <= fec_length; ++erasure_count)
               {
                  for (std::size_t start_position = 0; start_position < code_length; ++start_position)
                  {
                     block_type rs_block = rs_block_original;

                     corrupt_message_errors_erasures
                     (
                       rs_block,
                       error_mode::errors_erasures,
                       start_position,
                       erasure_count,
                       erasure_list,between_distance
                     );

                     if (!rs_decoder_->decode(rs_block,erasure_list))
                     {
                        print_codec_properties();
                        std::cout << "stage6() - Decoding Failure! start position: " << start_position << std::endl;
                        ++block_failures_;
                     }
                     else if (!is_block_equivelent(rs_block,message))
                     {
                        print_codec_properties();
                        std::cout << "stage6() - Error Correcting Failure! start position: " << start_position << std::endl;
                        ++block_failures_;
                     }
                     else if (rs_block.errors_detected != rs_block.errors_corrected)
                     {
                        print_codec_properties();
                        std::cout << "stage6() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                        ++block_failures_;
                     }

                     ++blocks_processed_;
                     erasure_list.clear();
                  }
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage7()
         {
            /*  Intermittent Error Combinations */

            const std::size_t initial_failure_count = block_failures_;

            for (std::size_t error_count = 1; error_count < (fec_length >> 1); ++error_count)
            {
               for (std::size_t start_position = 0; start_position < code_length; ++start_position)
               {
                  for (std::size_t scale = 1; scale < 5; ++scale)
                  {
                     block_type rs_block = rs_block_original;

                     corrupt_message_all_errors
                     (
                       rs_block,
                       error_count,
                       start_position,
                       scale
                     );

                     if (!rs_decoder_->decode(rs_block))
                     {
                        print_codec_properties();
                        std::cout << "stage7() - Decoding Failure! start position: " << start_position << std::endl;
                        ++block_failures_;
                     }
                     else if (!is_block_equivelent(rs_block,message))
                     {
                        print_codec_properties();
                        std::cout << "stage7() - Error Correcting Failure! start position: " << start_position << std::endl;
                        ++block_failures_;
                     }
                     else if (rs_block.errors_detected != rs_block.errors_corrected)
                     {
                        print_codec_properties();
                        std::cout << "stage7() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                        ++block_failures_;
                     }
                     else if (rs_block.errors_detected != error_count)
                     {
                        print_codec_properties();
                        std::cout << "stage7() - Error In The Number Of Detected Errors! Errors Detected: " << rs_block.errors_detected << std::endl;
                        ++block_failures_;
                     }
                     else if (rs_block.errors_corrected != error_count)
                     {
                        print_codec_properties();
                        std::cout << "stage7() - Error In The Number Of Corrected Errors! Errors Corrected: " << rs_block.errors_corrected << std::endl;
                        ++block_failures_;
                     }

                     ++blocks_processed_;
                  }
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage8()
         {
            /* Intermittent Erasure Combinations */

            const std::size_t initial_failure_count = block_failures_;

            erasure_locations_t erasure_list;

            for (std::size_t erasure_count = 1; erasure_count <= fec_length; ++erasure_count)
            {
               for (std::size_t start_position = 0; start_position < code_length; ++start_position)
               {
                  for (std::size_t scale = 4; scale < 5; ++scale)
                  {
                     block_type rs_block = rs_block_original;

                     corrupt_message_all_erasures
                     (
                       rs_block,
                       erasure_list,
                       erasure_count,
                       start_position,
                       scale
                     );

                     if (!rs_decoder_->decode(rs_block,erasure_list))
                     {
                        print_codec_properties();
                        std::cout << "stage8() - Decoding Failure! start position: " << start_position << "\t scale: " << scale << std::endl;
                        ++block_failures_;
                     }
                     else if (!is_block_equivelent(rs_block,message))
                     {
                        print_codec_properties();
                        std::cout << "stage8() - Error Correcting Failure! start position: " << start_position << "\t scale: " << scale <<std::endl;
                        ++block_failures_;
                     }
                     else if (rs_block.errors_detected != (rs_block.errors_corrected + rs_block.zero_numerators))
                     {
                        print_codec_properties();
                        std::cout << "stage8() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                        ++block_failures_;
                     }
                     else if (rs_block.errors_detected > erasure_count)
                     {
                        print_codec_properties();
                        std::cout << "stage8() - Error In The Number Of Detected Errors! Errors Detected: " << rs_block.errors_detected << std::endl;
                        ++block_failures_;
                     }
                     else if (rs_block.errors_corrected > erasure_count)
                     {
                        print_codec_properties();
                        std::cout << "stage8() - Error In The Number Of Corrected Errors! Errors Corrected: " << rs_block.errors_corrected << std::endl;
                        ++block_failures_;
                     }
                     ++blocks_processed_;
                     erasure_list.clear();
                  }
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage9()
         {
            /* Burst Interleaved Error and Erasure Combinations */

            const std::size_t initial_failure_count = block_failures_;

            erasure_locations_t erasure_list;

            for (std::size_t erasure_count = 1; erasure_count <= fec_length; ++erasure_count)
            {
               for (std::size_t start_position = 0; start_position < code_length; ++start_position)
               {
                  block_type rs_block = rs_block_original;

                  corrupt_message_interleaved_errors_erasures
                  (
                    rs_block,
                    start_position,
                    erasure_count,
                    erasure_list
                  );

                  if (!rs_decoder_->decode(rs_block,erasure_list))
                  {
                     print_codec_properties();
                     std::cout << "stage9() - Decoding Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (!is_block_equivelent(rs_block,message))
                  {
                     print_codec_properties();
                     std::cout << "stage9() - Error Correcting Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != rs_block.errors_corrected)
                  {
                     print_codec_properties();
                     std::cout << "stage9() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                     ++block_failures_;
                  }
                  ++blocks_processed_;
                  erasure_list.clear();
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage10()
         {
            /* Segmented Burst Errors */

            const std::size_t initial_failure_count = block_failures_;

            for (std::size_t start_position = 0; start_position < code_length; ++start_position)
            {
               for (std::size_t distance_between_blocks = 0; distance_between_blocks < 5; ++distance_between_blocks)
               {
                  block_type rs_block = rs_block_original;

                  corrupt_message_all_errors_segmented
                  (
                    rs_block,
                    start_position,
                    distance_between_blocks
                  );

                  if (!rs_decoder_->decode(rs_block))
                  {
                     print_codec_properties();
                     std::cout << "stage10() - Decoding Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (!is_block_equivelent(rs_block,message))
                  {
                     print_codec_properties();
                     std::cout << "stage10() - Error Correcting Failure! start position: " << start_position << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != rs_block.errors_corrected)
                  {
                     print_codec_properties();
                     std::cout << "stage10() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                     ++block_failures_;
                  }

                  ++blocks_processed_;
               }
            }

            return (block_failures_ == initial_failure_count);
         }

         bool stage11()
         {
            /* No Errors */

            const std::size_t initial_failure_count = block_failures_;

            block_type rs_block = rs_block_original;

            if (!rs_decoder_->decode(rs_block))
            {
               print_codec_properties();
               std::cout << "stage11() - Decoding Failure!" << std::endl;
               ++block_failures_;
            }
            else if (!is_block_equivelent(rs_block,message))
            {
               print_codec_properties();
               std::cout << "stage11() - Error Correcting Failure!" << std::endl;
               ++block_failures_;
            }
            else if (rs_block.errors_detected != 0)
            {
               print_codec_properties();
               std::cout << "stage11() - Error Correcting Failure!" << std::endl;
               ++block_failures_;
            }
            else if (rs_block.errors_corrected != 0)
            {
               print_codec_properties();
               std::cout << "stage11() - Error Correcting Failure!" << std::endl;
               ++block_failures_;
            }
            else if (rs_block.unrecoverable)
            {
               print_codec_properties();
               std::cout << "stage11() - Error Correcting Failure!" << std::endl;
               ++block_failures_;
            }

            ++blocks_processed_;

            return (block_failures_ == initial_failure_count);
         }

         bool stage12()
         {
            /* Random Errors Only */

            const std::size_t initial_failure_count = block_failures_;

            std::vector<std::size_t> random_error_index;
            generate_error_index((fec_length >> 1),random_error_index,0xA5A5A5A5);

            for (std::size_t error_count = 1; error_count <= (fec_length >> 1); ++error_count)
            {
               for (std::size_t error_index = 0; error_index < error_index_size; ++error_index)
               {
                  block_type rs_block = rs_block_original;

                  corrupt_message_all_errors_at_index
                  (
                    rs_block,
                    error_count,
                    error_index,
                    random_error_index
                  );

                  if (!rs_decoder_->decode(rs_block))
                  {
                     print_codec_properties();
                     std::cout << "stage12() - Decoding Failure! error index: " << error_index << std::endl;
                     ++block_failures_;
                  }
                  else if (!is_block_equivelent(rs_block,message))
                  {
                     print_codec_properties();
                     std::cout << "stage12() - Error Correcting Failure! error index: " << error_index << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != rs_block.errors_corrected)
                  {
                     print_codec_properties();
                     std::cout << "stage12() - Discrepancy between the number of errors detected and corrected. [" << rs_block.errors_detected << "," << rs_block.errors_corrected << "]" << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_detected != error_count)
                  {
                     print_codec_properties();
                     std::cout << "stage12() - Error In The Number Of Detected Errors! Errors Detected: " << rs_block.errors_detected << std::endl;
                     ++block_failures_;
                  }
                  else if (rs_block.errors_corrected != error_count)
                  {
                     print_codec_properties();
                     std::cout << "stage12() - Error In The Number Of Corrected Errors! Errors Corrected: " << rs_block.errors_corrected << std::endl;
                     ++block_failures_;
                  }

                  ++blocks_processed_;
               }
            }

            return (block_failures_ == initial_failure_count);
         }

      protected:

         codec_validator() {}

      private:

         codec_validator(const codec_validator&);
         const codec_validator& operator=(const codec_validator&);

         const galois::field& field_;
         galois::field_polynomial generator_polynomial_;
         encoder_type* rs_encoder_;
         decoder_type* rs_decoder_;
         block_type rs_block_original;
         const std::string&  message;
         const unsigned int genpoly_initial_index_;
         unsigned int blocks_processed_;
         unsigned int block_failures_;
      };

      template <std::size_t data_length>
      void create_messages(std::vector<std::string>& message_list, const bool full_test_set = false)
      {
         /* Various message bit patterns */

         message_list.clear();

         if (full_test_set)
         {
            for (std::size_t i = 0; i < 256; ++i)
            {
               message_list.push_back(std::string(data_length, static_cast<unsigned char>(i)));
            }
         }
         else
         {
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x00)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xAA)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xA5)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xAC)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xCA)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x5A)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xCC)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xF0)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x0F)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xFF)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x92)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x6D)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x77)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0x7A)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xA7)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xE5)));
            message_list.push_back(std::string(data_length,static_cast<unsigned char>(0xEB)));
         }

         std::string tmp_str = std::string(data_length,static_cast<unsigned char>(0x00));

         for (std::size_t i = 0; i < data_length; ++i)
         {
            tmp_str[i] = static_cast<unsigned char>(i);
         }

         message_list.push_back(tmp_str);

         for (int i = data_length - 1; i >= 0; --i)
         {
            tmp_str[i] = static_cast<unsigned char>(i);
         }

         message_list.push_back(tmp_str);

         for (std::size_t i = 0; i < data_length; ++i)
         {
            tmp_str[i] = (((i & 0x01) == 1) ? static_cast<unsigned char>(i) : 0x00);
         }

         message_list.push_back(tmp_str);

         for (std::size_t i = 0; i < data_length; ++i)
         {
            tmp_str[i] = (((i & 0x01) == 0) ? static_cast<unsigned char>(i) : 0x00);
         }

         message_list.push_back(tmp_str);

         for (int i = data_length - 1; i >= 0; --i)
         {
            tmp_str[i] = (((i & 0x01) == 1) ? static_cast<unsigned char>(i) : 0x00);
         }

         message_list.push_back(tmp_str);

         for (int i = data_length - 1; i >= 0; --i)
         {
            tmp_str[i] = (((i & 0x01) == 0) ? static_cast<unsigned char>(i) : 0x00);
         }

         message_list.push_back(tmp_str);

         tmp_str = std::string(data_length,static_cast<unsigned char>(0x00));

         for (std::size_t i = 0; i < (data_length >> 1); ++i)
         {
               tmp_str[i] = static_cast<unsigned char>(0xFF);
         }

         message_list.push_back(tmp_str);

         tmp_str = std::string(data_length,static_cast<unsigned char>(0xFF));

         for (std::size_t i = 0; i < (data_length >> 1); ++i)
         {
            tmp_str[i] = static_cast<unsigned char>(0x00);
         }

         message_list.push_back(tmp_str);
      }

      template <std::size_t field_descriptor, std::size_t gen_poly_index, std::size_t code_length, std::size_t fec_length>
      inline bool codec_validation_test(const std::size_t prim_poly_size,const unsigned int prim_poly[])
      {
         const unsigned int data_length = code_length - fec_length;

         galois::field field(field_descriptor,prim_poly_size,prim_poly);
         std::vector<std::string> message_list;
         create_messages<data_length>(message_list);

         for (std::size_t i = 0; i < message_list.size(); ++i)
         {
            codec_validator<code_length,fec_length>
               validator(field, gen_poly_index, message_list[i]);

            if (!validator.execute())
            {
               return false;
            }
         }

         return true;
      }

      template <std::size_t field_descriptor,
                std::size_t gen_poly_index,
                std::size_t code_length,
                std::size_t fec_length>
      inline bool shortened_codec_validation_test(const std::size_t prim_poly_size,const unsigned int prim_poly[])
      {
         typedef shortened_encoder<code_length,fec_length> encoder_type;
         typedef shortened_decoder<code_length,fec_length> decoder_type;

         const unsigned int data_length = code_length - fec_length;

         galois::field field(field_descriptor,prim_poly_size,prim_poly);
         std::vector<std::string> message_list;
         create_messages<data_length>(message_list);

         for (std::size_t i = 0; i < message_list.size(); ++i)
         {
            codec_validator<code_length,fec_length,encoder_type,decoder_type>
               validator(field,gen_poly_index,message_list[i]);

            if (!validator.execute())
            {
               return false;
            }

         }

         return true;
      }

      inline bool codec_validation_test00()
      {
         return codec_validation_test<8,120,255,  2>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255,  4>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255,  6>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 10>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 12>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 14>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 16>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 18>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 20>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 22>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 24>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 32>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 64>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 80>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255, 96>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) &&
                codec_validation_test<8,120,255,128>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) ;
      }

      inline bool codec_validation_test01()
      {
         return shortened_codec_validation_test<8,120,126,14>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) && /* Intelsat 1 RS Code */
                shortened_codec_validation_test<8,120,194,16>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) && /* Intelsat 2 RS Code */
                shortened_codec_validation_test<8,120,219,18>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) && /* Intelsat 3 RS Code */
                shortened_codec_validation_test<8,120,225,20>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) && /* Intelsat 4 RS Code */
                shortened_codec_validation_test<8,  1,204,16>(galois::primitive_polynomial_size05,galois::primitive_polynomial05) && /* DBV/MPEG-2 TSP RS Code */
                shortened_codec_validation_test<8,  1,104,27>(galois::primitive_polynomial_size05,galois::primitive_polynomial05) && /* Magnetic Storage Outer RS Code */
                shortened_codec_validation_test<8,  1,204,12>(galois::primitive_polynomial_size05,galois::primitive_polynomial05) && /* Magnetic Storage Inner RS Code */
                shortened_codec_validation_test<8,120, 72,10>(galois::primitive_polynomial_size06,galois::primitive_polynomial06) ;  /* VDL Mode 3 RS Code */
      }

   } // namespace reed_solomon

} // namespace schifra

#endif
