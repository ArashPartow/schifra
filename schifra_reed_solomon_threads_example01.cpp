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
                encoder and decoder capabilities in a threaded context. One
                must note that the number of threads should not exceed the
                architecture's ability to efficiently and productively run the
                threads. A simple limiting strategy would be not to have more
                threads than the number of available cores on the processor.
*/


#include <cstddef>
#include <iostream>
#include <string>
#include <limits>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/thread.hpp>

#include "schifra_galois_field.hpp"
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_error_processes.hpp"
#include "schifra_ecc_traits.hpp"
#include "schifra_utilities.hpp"


const std::size_t round_count = 1000;

template <typename Encoder, typename Decoder>
class process
{
public:

   process(const unsigned int& process_id,
           const Encoder& encoder,
           const Decoder& decoder,
           const std::vector<std::string>& message_list)
   : process_id_(process_id),
     total_time_(0.0),
     encoder_(encoder),
     decoder_(decoder),
     message_list_(message_list)
   {}

   process& operator=(const process& proc)
   {
      process_id_ = proc.process_id_;
      total_time_ = proc.total_time_;
      return *this;
   }

   double time() { return total_time_; }

   inline void execute()
   {
      schifra::traits::equivalent_encoder_decoder<Encoder,Decoder>();
      typedef schifra::reed_solomon::block<Encoder::trait::code_length,Encoder::trait::fec_length> block_type;

      std::vector<block_type> block_list(message_list_.size());

      for (std::size_t i = 0; i < message_list_.size(); ++i)
      {
         if (!encoder_.encode(message_list_[i],block_list[i]))
         {
            std::cout << "[" << process_id_ << "] (0)Error - Critical encoding failure!" << std::endl;
            return;
         }
         schifra::corrupt_message_all_errors00(block_list[i],0,3);
      }

      schifra::utils::timer timer;
      timer.start();

      for (std::size_t k = 0; k < round_count; ++k)
      {
         for (std::size_t i = 0; i < message_list_.size(); ++i)
         {
            if (!decoder_.decode(block_list[i]))
            {
               std::cout << "[" << process_id_ << "] (1)Error - Critical decoding failure!" << std::endl;
               return;
            }
            else if (!schifra::is_block_equivelent(block_list[i],message_list_[i]))
            {
               std::cout << "[" << process_id_ << "] (2)Error - Error correction failed!" << std::endl;
               return;
            }
         }
      }

      timer.stop();
      total_time_ = timer.time();
   }

private:

   unsigned int process_id_;
   double total_time_;
   const Encoder& encoder_;
   const Decoder& decoder_;
   const std::vector<std::string>& message_list_;
};

void generate_messages(const std::size_t data_length, std::vector<std::string>& message_list)
{
   for (unsigned int c = 0; c < 256; ++c)
   {
      message_list.push_back(std::string(data_length,static_cast<unsigned char>(c)));
   }
}

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

   typedef schifra::reed_solomon::encoder<code_length,fec_length> encoder_type;
   typedef schifra::reed_solomon::decoder<code_length,fec_length> decoder_type;
   typedef process<encoder_type,decoder_type>                     process_type;
   typedef boost::shared_ptr<process_type>                        process_ptr_type;

   /* Instantiate Encoder and Decoder (Codec) */
   encoder_type encoder(field,generator_polynomial);
   decoder_type decoder(field,generator_polynomial_index);

   std::vector<std::string> message_list;

   generate_messages(data_length,message_list);

   const unsigned int max_thread_count = 4; // number of functional cores.
   std::vector<process_ptr_type> process_list;

   boost::thread_group threads;

   for (unsigned int i = 0; i < max_thread_count; ++i)
   {
      process_list.push_back(process_ptr_type(new process_type(i,encoder,decoder,message_list)));
      threads.create_thread(boost::bind(&process_type::execute,process_list[i]));
   }

   threads.join_all();

   double time = -1.0;

   /* Determine the process with the longest running time. */
   for (std::size_t i = 0; i < process_list.size(); ++i)
   {
      time = ((time < process_list[i]->time()) ? process_list[i]->time() : time);
   }

   double mbps = (max_thread_count * round_count * message_list.size() * data_length * 8.0) / (1048576.0 * time);

   std::cout << "Blocks decoded: " << max_thread_count * round_count * message_list.size() << "\tTime: " << time <<"sec\tRate: " << mbps << "Mbps" << std::endl;

   return 0;
}

