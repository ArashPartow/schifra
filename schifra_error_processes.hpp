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


#ifndef INCLUDE_SCHIFRA_ERROR_PROCESSES_HPP
#define INCLUDE_SCHIFRA_ERROR_PROCESSES_HPP


#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <deque>
#include <vector>

#include "schifra_reed_solomon_block.hpp"
#include "schifra_fileio.hpp"


namespace schifra
{

   template <std::size_t code_length, std::size_t fec_length>
   inline void add_erasure_error(const std::size_t& position, reed_solomon::block<code_length,fec_length>& block)
   {
      block[position] = (~block[position]) & 0xFF; // Or one can simply equate to zero
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void add_error(const std::size_t& position, reed_solomon::block<code_length,fec_length>& block)
   {
      block[position] = (~block[position]) & 0xFF;
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void add_error_4bit_symbol(const std::size_t& position, reed_solomon::block<code_length,fec_length>& block)
   {
      block[position] = (~block[position]) & 0x0F;
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_all_errors00(reed_solomon::block<code_length,fec_length>& rsblock,
                                           const std::size_t& start_position,
                                           const std::size_t& scale = 1)
   {
      for (std::size_t i = 0; i < (fec_length >> 1); ++i)
      {
         add_error((start_position + scale * i) % code_length,rsblock);
      }
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_all_errors_wth_mask(reed_solomon::block<code_length,fec_length>& rsblock,
                                                   const std::size_t& start_position,
                                                   const int& mask,
                                                   const std::size_t& scale = 1)
   {
      for (std::size_t i = 0; i < (fec_length >> 1); ++i)
      {
         std::size_t position = (start_position + scale * i) % code_length;
         rsblock[position] = (~rsblock[position]) & mask;

      }
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_all_errors(schifra::reed_solomon::block<code_length,fec_length>& rsblock,
                                          const std::size_t error_count,
                                          const std::size_t& start_position,
                                          const std::size_t& scale = 1)
   {
      for (std::size_t i = 0; i < error_count; ++i)
      {
         add_error((start_position + scale * i) % code_length,rsblock);
      }
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_all_erasures00(reed_solomon::block<code_length,fec_length>& rsblock,
                                              reed_solomon::erasure_locations_t& erasure_list,
                                              const std::size_t& start_position,
                                              const std::size_t& scale = 1)
   {
      std::size_t erasures[code_length];

      for (std::size_t i = 0; i < code_length; ++i) erasures[i] = 0;

      for (std::size_t i = 0; i < fec_length; ++i)
      {
         std::size_t error_position = (start_position + scale * i) % code_length;
         add_erasure_error(error_position,rsblock);
         erasures[error_position] = 1;
      }

      for (std::size_t i = 0; i < code_length; ++i)
      {
         if (erasures[i] == 1) erasure_list.push_back(i);
      }
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_all_erasures(reed_solomon::block<code_length,fec_length>& rsblock,
                                            reed_solomon::erasure_locations_t& erasure_list,
                                            const std::size_t erasure_count,
                                            const std::size_t& start_position,
                                            const std::size_t& scale = 1)
   {
      std::size_t erasures[code_length];

      for (std::size_t i = 0; i < code_length; ++i) erasures[i] = 0;

      for (std::size_t i = 0; i < erasure_count; ++i)
      {
         /* Note: Must make sure duplicate erasures are not added */
         std::size_t error_position = (start_position + scale * i) % code_length;
         add_erasure_error(error_position,rsblock);
         erasures[error_position] = 1;
      }

      for (std::size_t i = 0; i < code_length; ++i)
      {
         if (erasures[i] == 1) erasure_list.push_back(i);
      }
   }

   namespace error_mode
   {
      enum type
      {
         errors_erasures, // Errors first then erasures
         erasures_errors  // Erasures first then errors
      };
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_errors_erasures(reed_solomon::block<code_length,fec_length>& rsblock,
                                               const error_mode::type& mode,
                                               const std::size_t& start_position,
                                               const std::size_t& erasure_count,
                                               reed_solomon::erasure_locations_t& erasure_list,
                                               const std::size_t between_space = 0)
   {
      std::size_t error_count = (fec_length - erasure_count) >> 1;

      if ((2 * error_count) + erasure_count > fec_length)
      {
         std::cout << "corrupt_message_errors_erasures() - ERROR Too many erasures and errors!" << std::endl;
         std::cout << "Error Count:   " << error_count << std::endl;
         std::cout << "Erasure Count: " << error_count << std::endl;

         return;
      }

      std::size_t erasures[code_length];

      for (std::size_t i = 0; i < code_length; ++i) erasures[i] = 0;

      std::size_t error_position = 0;

      switch (mode)
      {
         case error_mode::erasures_errors : {
                                               for (std::size_t i = 0; i < erasure_count; ++i)
                                               {
                                                  error_position = (start_position + i) % code_length;
                                                  add_erasure_error(error_position,rsblock);
                                                  erasures[error_position] = 1;
                                               }

                                               for (std::size_t i = 0; i < error_count; ++i)
                                               {
                                                  error_position = (start_position + erasure_count + between_space + i) % code_length;
                                                  add_error(error_position,rsblock);
                                               }
                                            }
                                            break;

         case error_mode::errors_erasures : {
                                               for (std::size_t i = 0; i < error_count; ++i)
                                               {
                                                  error_position = (start_position + i) % code_length;
                                                  add_error(error_position,rsblock);
                                               }

                                               for (std::size_t i = 0; i < erasure_count; ++i)
                                               {
                                                  error_position = (start_position + error_count + between_space + i) % code_length;
                                                  add_erasure_error(error_position,rsblock);
                                                  erasures[error_position] = 1;
                                               }
                                            }
                                            break;
      }

      for (std::size_t i = 0; i < code_length; ++i)
      {
         if (erasures[i] == 1) erasure_list.push_back(i);
      }

   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_interleaved_errors_erasures(reed_solomon::block<code_length,fec_length>& rsblock,
                                                           const std::size_t& start_position,
                                                           const std::size_t& erasure_count,
                                                           reed_solomon::erasure_locations_t& erasure_list)
   {
      std::size_t error_count = (fec_length - erasure_count) >> 1;

      if ((2 * error_count) + erasure_count > fec_length)
      {
         std::cout << "corrupt_message_interleaved_errors_erasures() - [1] ERROR Too many erasures and errors!" << std::endl;
         std::cout << "Error Count:   " << error_count << std::endl;
         std::cout << "Erasure Count: " << error_count << std::endl;

         return;
      }

      std::size_t erasures[code_length];

      for (std::size_t i = 0; i < code_length; ++i) erasures[i] = 0;

      std::size_t e = 0;
      std::size_t s = 0;
      std::size_t i = 0;

      while ((e < error_count) || (s < erasure_count) || (i < (error_count + erasure_count)))
      {
        std::size_t error_position = (start_position + i) % code_length;

        if (((i & 0x01) == 0) && (s < erasure_count))
        {
           add_erasure_error(error_position,rsblock);
           erasures[error_position] = 1;
           s++;
        }
        else if (((i & 0x01) == 1) && (e < error_count))
        {
           e++;
           add_error(error_position,rsblock);
        }
        ++i;
      }

      for (std::size_t j = 0; j < code_length; ++j)
      {
         if (erasures[j] == 1) erasure_list.push_back(j);
      }

      if ((2 * e) + erasure_list.size() > fec_length)
      {
         std::cout << "corrupt_message_interleaved_errors_erasures() - [2] ERROR Too many erasures and errors!" << std::endl;
         std::cout << "Error Count:   " << error_count << std::endl;
         std::cout << "Erasure Count: " << error_count << std::endl;

         return;
      }
   }

   namespace details
   {
      template <std::size_t code_length, std::size_t fec_length, bool t>
      struct corrupt_message_all_errors_segmented_impl
      {
         static void process(reed_solomon::block<code_length,fec_length>& rsblock,
                             const std::size_t& start_position,
                             const std::size_t& distance_between_blocks = 1)
         {
            std::size_t block_1_error_count = (fec_length >> 2);
            std::size_t block_2_error_count = (fec_length >> 1) - block_1_error_count;

            for (std::size_t i = 0; i < block_1_error_count; ++i)
            {
               add_error((start_position + i) % code_length,rsblock);
            }

            std::size_t new_start_position = (start_position + (block_1_error_count)) + distance_between_blocks;

            for (std::size_t i = 0; i < block_2_error_count; ++i)
            {
               add_error((new_start_position + i) % code_length,rsblock);
            }
         }
      };

      template <std::size_t code_length, std::size_t fec_length>
      struct corrupt_message_all_errors_segmented_impl<code_length,fec_length,false>
      {
         static void process(reed_solomon::block<code_length,fec_length>&,
                             const std::size_t&, const std::size_t&)
         {}
      };
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_all_errors_segmented(reed_solomon::block<code_length,fec_length>& rsblock,
                                                    const std::size_t& start_position,
                                                    const std::size_t& distance_between_blocks = 1)
   {
      details::corrupt_message_all_errors_segmented_impl<code_length,fec_length,(fec_length > 2)>::
                  process(rsblock,start_position,distance_between_blocks);
   }

   inline bool check_for_duplicate_erasures(const std::vector<int>& erasure_list)
   {
      for (std::size_t i = 0; i < erasure_list.size(); ++i)
      {
         for (std::size_t j = i + 1; j < erasure_list.size(); ++j)
         {
            if (erasure_list[i] == erasure_list[j])
            {
               return false;
            }
         }
      }

      return true;
   }

   inline void dump_erasure_list(const schifra::reed_solomon::erasure_locations_t& erasure_list)
   {
      for (std::size_t i = 0; i < erasure_list.size(); ++i)
      {
         std::cout << "[" << i << "," << erasure_list[i] << "] ";
      }

      std::cout << std::endl;
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline bool is_block_equivelent(const reed_solomon::block<code_length,fec_length>& rsblock,
                                   const std::string& data,
                                   const bool display    = false,
                                   const bool all_errors = false)
   {
      std::string::const_iterator it = data.begin();

      bool error_found = false;

      for (std::size_t i = 0; i < code_length - fec_length; ++i, ++it)
      {
         if (static_cast<char>(rsblock.data[i] & 0xFF) != (*it))
         {
            error_found = true;

            if (display)
            {
               printf("is_block_equivelent() - Error at loc : %02d\td1: %02X\td2: %02X\n",
                      static_cast<unsigned int>(i),
                      rsblock.data[i],
                      static_cast<unsigned char>(*it));
            }

            if (!all_errors)
               return false;
         }
      }

      return !error_found;
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline bool are_blocks_equivelent(const reed_solomon::block<code_length,fec_length>& block1,
                                     const reed_solomon::block<code_length,fec_length>& block2,
                                     const std::size_t span = code_length,
                                     const bool display     = false,
                                     const bool all_errors  = false)
   {
      bool error_found = false;

      for (std::size_t i = 0; i < span; ++i)
      {
         if (block1[i] != block2[i])
         {
            error_found = true;

            if (display)
            {
               printf("are_blocks_equivelent() - Error at loc : %02d\td1: %04X\td2: %04X\n",
                      static_cast<unsigned int>(i),
                      block1[i],
                      block2[i]);
            }

            if (!all_errors)
               return false;
         }
      }

      return !error_found;
   }

   template <std::size_t code_length, std::size_t fec_length, std::size_t stack_size>
   inline bool block_stacks_equivelent(const reed_solomon::block<code_length,fec_length> block_stack1[stack_size],
                                       const reed_solomon::block<code_length,fec_length> block_stack2[stack_size])
   {
      for (std::size_t i = 0; i < stack_size; ++i)
      {
         if (!are_blocks_equivelent(block_stack1[i],block_stack2[i]))
         {
            return false;
         }
      }

      return true;
   }

   template <std::size_t block_length, std::size_t stack_size>
   inline bool block_stacks_equivelent(const reed_solomon::data_block<std::size_t,block_length> block_stack1[stack_size],
                                       const reed_solomon::data_block<std::size_t,block_length> block_stack2[stack_size])
   {
      for (std::size_t i = 0; i < stack_size; ++i)
      {
         for (std::size_t j = 0; j < block_length; ++j)
         {
            if (block_stack1[i][j] != block_stack2[i][j])
            {
               return false;
            }
         }
      }

      return true;
   }

   inline void corrupt_file_with_burst_errors(const std::string& file_name,
                                              const long& start_position,
                                              const long& burst_length)
   {
      if (!schifra::fileio::file_exists(file_name))
      {
         std::cout << "corrupt_file() - Error: " << file_name << " does not exist!" << std::endl;
         return;
      }

      if (static_cast<std::size_t>(start_position + burst_length) >= schifra::fileio::file_size(file_name))
      {
         std::cout << "corrupt_file() - Error: Burst error out of bounds." << std::endl;
         return;
      }

      std::vector<char> data(burst_length);

      std::ifstream ifile(file_name.c_str(), std::ios::in | std::ios::binary);

      if (!ifile)
      {
         return;
      }

      ifile.seekg(start_position,std::ios_base::beg);
      ifile.read(&data[0],burst_length);
      ifile.close();

      for (long i = 0; i < burst_length; ++i)
      {
         data[i] = ~data[i];
      }

      std::ofstream ofile(file_name.c_str(), std::ios::in | std::ios::out | std::ios::binary);

      if (!ofile)
      {
         return;
      }

      ofile.seekp(start_position,std::ios_base::beg);
      ofile.write(&data[0],burst_length);
      ofile.close();
   }

   static const std::size_t global_random_error_index[] =
                         {
                            13,  170,  148,   66,  228,  208,  182,   92,
                             4,  137,   97,   99,  237,  151,   15,    0,
                           119,  243,   41,  222,   33,  211,  188,    5,
                            44,   30,  210,  111,   54,   79,   61,  223,
                           239,  149,   73,  115,  201,  234,  194,   62,
                           147,   70,   19,   49,   72,   52,  164,   29,
                           102,  225,  203,  153,   18,  205,   40,  217,
                           165,  177,  166,  134,  236,   68,  231,  154,
                           116,  136,   47,  240,   46,   89,  120,  183,
                           242,   28,  161,  226,  241,  230,   10,  131,
                           207,  132,   83,  171,  202,  195,  227,  206,
                           112,   88,   90,  146,  117,  180,   26,   78,
                           118,  254,  107,  110,  220,    7,  192,  187,
                            31,  175,  127,  209,   32,   12,   84,  128,
                           190,  156,   95,  105,  104,  246,   91,  215,
                           219,  142,   36,  186,  247,  233,  167,  133,
                           160,   16,  140,  169,   23,   96,  155,  235,
                           179,   76,  253,  103,  238,   67,   35,  121,
                           100,   27,  213,   58,   77,  248,  174,   39,
                           214,   56,   42,  200,  106,   21,  129,  114,
                           252,  113,  168,   53,   25,  216,   64,  232,
                            81,   75,    2,  224,  250,   60,  135,  204,
                            48,  196,   94,   63,  244,  191,   93,  126,
                           138,  159,    9,   85,  249,   34,  185,  163,
                            17,   65,  184,   82,  109,  172,  108,   69,
                           150,    3,   20,  221,  162,  212,  152,   59,
                           198,   74,  229,   55,   87,  178,  141,  199,
                            57,  130,   80,  173,  101,  122,  144,   51,
                           139,   11,    8,  125,  158,  124,  123,   37,
                            14,   24,   22,   43,  197,   50,   98,    6,
                           176,  251,   86,  218,  193,   71,  145,    1,
                            45,   38,  189,  143,  245,  157,  181
                         };

   static const std::size_t error_index_size = sizeof(global_random_error_index) / sizeof(std::size_t);

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_all_errors_at_index(schifra::reed_solomon::block<code_length,fec_length>& rsblock,
                                                   const std::size_t error_count,
                                                   const std::size_t& error_index_start_position,
                                                   const bool display_positions = false)
   {
      schifra::reed_solomon::block<code_length,fec_length> tmp_rsblock = rsblock;

      for (std::size_t i = 0; i < error_count; ++i)
      {
         std::size_t error_position = (global_random_error_index[(error_index_start_position + i) % error_index_size]) % code_length;

         add_error(error_position,rsblock);

         if (display_positions)
         {
            std::cout << "Error index: " << error_position << std::endl;
         }
      }
   }

   template <std::size_t code_length, std::size_t fec_length>
   inline void corrupt_message_all_errors_at_index(schifra::reed_solomon::block<code_length,fec_length>& rsblock,
                                                   const std::size_t error_count,
                                                   const std::size_t& error_index_start_position,
                                                   const std::vector<std::size_t>& random_error_index,
                                                   const bool display_positions = false)
   {
      for (std::size_t i = 0; i < error_count; ++i)
      {
         std::size_t error_position = (random_error_index[(error_index_start_position + i) % random_error_index.size()]) % code_length;

         add_error(error_position,rsblock);

         if (display_positions)
         {
            std::cout << "Error index: " << error_position << std::endl;
         }
      }
   }

   inline void generate_error_index(const std::size_t index_size,
                                    std::vector<std::size_t>& random_error_index,
                                    std::size_t seed)
   {
      if (0 == seed)
      {
         seed = 0xA5A5A5A5;
      }

      ::srand(static_cast<unsigned int>(seed));

      std::deque<std::size_t> index_list;

      for (std::size_t i = 0; i < index_size; ++i)
      {
         index_list.push_back(i);
      }

      random_error_index.reserve(index_size);
      random_error_index.resize(0);

      while (!index_list.empty())
      {
         // possibly the worst way of doing this.
         std::size_t index = ::rand() % index_list.size();

         random_error_index.push_back(index_list[index]);
         index_list.erase(index_list.begin() + index);
      }
   }

} // namespace schifra

#endif
