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
   Description: This example will demonstrate the use of the Reed-Solomon
                erasure channel coding capabilities in a threaded context.
                One must note that the number of threads should not exceed
                the architecture's ability to efficiently and productively
                run the threads. A simple limiting strategy would be not to
                have more threads than the number of available cores on the
                processor.
*/


#include <cstddef>
#include <iostream>
#include <string>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_erasure_channel.hpp"
#include "schifra_ecc_traits.hpp"
#include "schifra_utilities.hpp"


const std::size_t round_count = 100;

template <typename Encoder, typename Decoder>
class erasure_process
{
public:

   erasure_process(const unsigned int& process_id,
                   schifra::galois::field& field,
                   schifra::galois::field_polynomial& generator_polynomial,
                   const std::size_t& generator_polynomial_index)
   : process_id_(process_id),
     total_time_(0.0),
     encoder_(field,generator_polynomial),
     decoder_(field,generator_polynomial_index)
   {}

   erasure_process& operator=(const erasure_process& ep)
   {
      process_id_ = ep.process_id_;
      total_time_ = ep.total_time_;
      return *this;
   }

   double time() { return total_time_; }

   void execute()
   {
      schifra::traits::equivalent_encoder_decoder<Encoder,Decoder>();

      const std::size_t code_length = Encoder::trait::code_length;
      const std::size_t fec_length  = Encoder::trait::fec_length;
      const std::size_t data_length = Encoder::trait::data_length;
      const std::size_t stack_size  = Encoder::trait::code_length;
      const std::size_t data_size   = stack_size * data_length;

      typedef schifra::reed_solomon::block<code_length,fec_length> block_type;

      block_type block_stack[stack_size];
      unsigned char send_data[data_size];
      unsigned char recv_data[data_size];

      schifra::utils::timer timer;
      total_time_ = 0.0;

      for (std::size_t k = 0; k < round_count; ++k)
      {
         /* Populate block stack with data */
         for (std::size_t i = 0; i < data_size; ++i)
         {
            send_data[i] = static_cast<unsigned char>((i * 3 + 7 * k) & 0xFF);
         }

         schifra::reed_solomon::copy<unsigned char,code_length,fec_length,stack_size>(send_data,data_size,block_stack);

         timer.start();

         schifra::reed_solomon::erasure_channel_stack_encode<code_length,fec_length>(encoder_,block_stack);

         /* Add Erasures - Simulate network packet loss (e.g: UDP) */
         schifra::reed_solomon::erasure_locations_t missing_row_index;
         missing_row_index.clear();

         for (std::size_t i = 0; i < fec_length; ++i)
         {
            std::size_t missing_index = (k + (i * 4)) % stack_size;
            block_stack[missing_index].clear();
            missing_row_index.push_back(missing_index);
         }

         schifra::reed_solomon::erasure_channel_stack_decode<code_length,fec_length>(decoder_,missing_row_index,block_stack);

         schifra::reed_solomon::copy<unsigned char,code_length,fec_length,stack_size>(block_stack,recv_data);

         timer.stop();
         total_time_ += timer.time();

         for (std::size_t i = 0; i < data_size; ++i)
         {
            if (recv_data[i] != send_data[i])
            {
               std::cout << "[" << process_id_ << "] Error: Final block stack comparison failed! stack: " << i << std::endl;
               return;
            }
         }
      }
   }

private:

   unsigned int process_id_;
   double total_time_;
   Encoder encoder_;
   Decoder decoder_;
};

int main()
{
   /* Finite Field Parameters */
   const std::size_t field_descriptor                =   8;
   const std::size_t generator_polynomial_index      = 120;
   const std::size_t generator_polynomial_root_count = 128;

   /* Reed Solomon Code Parameters */
   const std::size_t code_length = 255;
   const std::size_t fec_length  = 128;
   const std::size_t data_length = code_length - fec_length;

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

   typedef schifra::reed_solomon::encoder<code_length,fec_length>              encoder_type;
   typedef schifra::reed_solomon::erasure_code_decoder<code_length,fec_length> decoder_type;

   typedef erasure_process<encoder_type,decoder_type> erasure_process_type;
   typedef boost::shared_ptr<erasure_process_type>    erasure_process_ptr_type;

   const unsigned int max_thread_count = 4; // number of functional cores.

   std::vector<erasure_process_ptr_type> erasure_process_list;

   boost::thread_group threads;

   for (unsigned int i = 0; i < max_thread_count; ++i)
   {
      erasure_process_list.push_back(erasure_process_ptr_type(new
                                     erasure_process_type
                                     (
                                       i,
                                       field,
                                       generator_polynomial,
                                       generator_polynomial_index
                                     )));

      threads.create_thread(boost::bind(&erasure_process_type::execute,erasure_process_list[i]));
   }

   threads.join_all();

   double time = -1.0;

   /* Determine the process with the longest running time. */
   for (std::size_t i = 0; i < erasure_process_list.size(); ++i)
   {
      time = ((time < erasure_process_list[i]->time()) ? erasure_process_list[i]->time() : time);
   }

   double mbps = (max_thread_count * round_count * 8.0 * code_length * data_length) / (1048576.0 * time);

   std::cout << "Blocks decoded: " << max_thread_count * round_count * code_length << "\tTime: " << time <<"sec\tRate: " << mbps << "Mbps" << std::endl;

   return 0;
}

