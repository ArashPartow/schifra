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


int main()
{
   const std::size_t block_length   = 255;
   const std::size_t stack_size     = 255;
   const std::size_t max_iterations = 10000;

   schifra::reed_solomon::data_block<std::size_t,block_length> block_stack[stack_size];

   for (std::size_t i = 0; i < stack_size; ++i)
   {
      for (std::size_t j = 0; j < block_length; ++j)
      {
         block_stack[i][j] = i;
      }
   }

   schifra::utils::timer timer;
   timer.start();

   for (std::size_t i = 0; i < max_iterations; ++i)
   {
      schifra::reed_solomon::interleave  <std::size_t,block_length,stack_size>(block_stack);
      schifra::reed_solomon::deinterleave<std::size_t,block_length,stack_size>(block_stack);
   }

   timer.stop();

   double mbps = ((max_iterations * stack_size * block_length) * 8.0) / (1048576.0 * timer.time());

   std::cout << "Interleave Rate: " << mbps << "Mbps" << std::endl;

   return 0;
}
