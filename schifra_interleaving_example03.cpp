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

#include "schifra_reed_solomon_interleaving.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_utilities.hpp"


template <std::size_t block_length, std::size_t stack_size>
bool block_stacks_equivelent(const schifra::reed_solomon::data_block<std::size_t,block_length> block_stack1[stack_size],
                             const schifra::reed_solomon::data_block<std::size_t,block_length> block_stack2[stack_size])
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

int main()
{
   const std::size_t block_length   =   10;
   const std::size_t stack_size     =  255;
   const std::size_t max_iterations = 1000;

   schifra::reed_solomon::data_block<std::size_t,block_length> block_stack   [stack_size];
   schifra::reed_solomon::data_block<std::size_t,block_length> or_block_stack[stack_size];

   for (std::size_t i = 0; i < stack_size; ++i)
   {
      for (std::size_t j = 0; j < block_length; ++j)
      {
         block_stack[i][j] = i;
         or_block_stack[i][j] = i;
      }
   }

   for (std::size_t i = block_length - 4; i < block_length; ++i)
   {
      or_block_stack[stack_size - 1][i] = 0xA5;
      block_stack   [stack_size - 1][i] = 0xA5;
   }

   schifra::utils::timer timer;
   timer.start();

   for (std::size_t i = 0; i < max_iterations; ++i)
   {
      schifra::reed_solomon::interleave  <std::size_t,block_length>(block_stack,stack_size,block_length - 4);
      schifra::reed_solomon::deinterleave<std::size_t,block_length>(block_stack,stack_size,block_length - 4);
   }

   timer.stop();

   if (!block_stacks_equivelent<block_length,stack_size>(block_stack,or_block_stack))
   {
      std::cout <<"ERROR!" << std::endl;
      return 1;
   }

   double mbps = ((max_iterations * stack_size * block_length) * 8.0) / (1048576.0 * timer.time());

   std::cout << "Interleave Rate: " << mbps << "Mbps" << std::endl;

   return 0;
}
