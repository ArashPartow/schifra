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


#ifndef INCLUDE_SCHIFRA_REED_SOLOMON_INTERLEAVING_HPP
#define INCLUDE_SCHIFRA_REED_SOLOMON_INTERLEAVING_HPP


#include <cstddef>
#include <iostream>
#include <string>

#include "schifra_reed_solomon_block.hpp"


namespace schifra
{

   namespace reed_solomon
   {

      template <std::size_t code_length, std::size_t fec_length>
      inline void interleave(block<code_length,fec_length> (&block_stack)[code_length])
      {
         for (std::size_t i = 0; i < code_length; ++i)
         {
            for (std::size_t j = i + 1; j < code_length; ++j)
            {
               typename block<code_length,fec_length>::symbol_type tmp = block_stack[i][j];
               block_stack[i][j] = block_stack[j][i];
               block_stack[j][i] = tmp;
           }
         }
      }

      template <std::size_t code_length, std::size_t fec_length, std::size_t row_count>
      inline void interleave(block<code_length,fec_length> (&block_stack)[row_count])
      {
         block<code_length,fec_length> auxiliary_stack[row_count];

         std::size_t aux_row   = 0;
         std::size_t aux_index = 0;

         for (std::size_t index = 0; index < code_length; ++index)
         {
            for (std::size_t row = 0; row < row_count; ++row)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_index == code_length)
               {
                  aux_index = 0;
                  aux_row++;
               }
            }
         }

         copy<code_length,fec_length,row_count>(auxiliary_stack,block_stack);
      }

      template <std::size_t code_length, std::size_t fec_length, std::size_t row_count>
      inline void interleave(block<code_length,fec_length,row_count> (&block_stack)[row_count],
                             const std::size_t partial_code_length)
      {
         if (partial_code_length == code_length)
         {
            interleave<code_length,fec_length,row_count>(block_stack);
         }
         else
         {
            block<code_length,fec_length,row_count> auxiliary_stack[row_count];

            std::size_t aux_row   = 0;
            std::size_t aux_index = 0;

            for (std::size_t index = 0; index < partial_code_length; ++index)
            {
               for (std::size_t row = 0; row < row_count; ++row)
               {
                  auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

                  if (++aux_index == code_length)
                  {
                     aux_index = 0;
                     aux_row++;
                  }
               }
            }

            for (std::size_t index = partial_code_length; index < code_length; ++index)
            {
               for (std::size_t row = 0; row < row_count - 1; ++row)
               {
                  auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

                  if (++aux_index == code_length)
                  {
                     aux_index = 0;
                     aux_row++;
                  }
               }
            }

            for (std::size_t row = 0; row < row_count - 1; ++row)
            {
               for (std::size_t index = 0; index < code_length - fec_length; ++index)
               {
                  block_stack[row].data[index] = auxiliary_stack[row].data[index];
               }
               for (std::size_t index = 0; index < fec_length; ++index)
               {
                  block_stack[row].fec[index] = auxiliary_stack[row].fec[index];
               }
            }

            for (std::size_t index = 0; index < partial_code_length; ++index)
            {
               block_stack[row_count - 1][index] = auxiliary_stack[row_count - 1][index];
            }
         }
      }

      template <typename T, std::size_t block_length>
      inline void interleave(data_block<T,block_length> (&block_stack)[block_length])
      {
         for (std::size_t i = 0; i < block_length; ++i)
         {
            for (std::size_t j = i + 1; j < block_length; ++j)
            {
               T tmp = block_stack[i][j];
               block_stack[i][j] = block_stack[j][i];
               block_stack[j][i] = tmp;
           }
         }
      }

      template <typename T, std::size_t block_length, std::size_t row_count>
      inline void interleave(data_block<T,block_length> (&block_stack)[row_count])
      {
         data_block<T,block_length> auxiliary_stack[row_count];

         std::size_t aux_row   = 0;
         std::size_t aux_index = 0;

         for (std::size_t index = 0; index < block_length; ++index)
         {
            for (std::size_t row = 0; row < row_count; ++row)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_index == block_length)
               {
                  aux_index = 0;
                  aux_row++;
               }
            }
         }

         copy<T,block_length,row_count>(auxiliary_stack,block_stack);
      }

      template <typename T, std::size_t block_length, std::size_t row_count>
      inline void interleave(data_block<T,block_length> (&block_stack)[row_count],
                             const std::size_t partial_block_length)
      {
         if (partial_block_length == block_length)
         {
            interleave<T,block_length,row_count>(block_stack);
         }
         else
         {
            data_block<T,block_length> auxiliary_stack[row_count];

            std::size_t aux_row   = 0;
            std::size_t aux_index = 0;

            for (std::size_t index = 0; index < partial_block_length; ++index)
            {
               for (std::size_t row = 0; row < row_count; ++row)
               {
                  auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

                  if (++aux_index == block_length)
                  {
                     aux_index = 0;
                     aux_row++;
                  }
               }
            }

            for (std::size_t index = partial_block_length; index < block_length; ++index)
            {
               for (std::size_t row = 0; row < row_count - 1; ++row)
               {
                  auxiliary_stack[aux_row][aux_index] = block_stack[row][index];
                  if (++aux_index == block_length)
                  {
                     aux_index = 0;
                     aux_row++;
                  }
               }
            }

            for (std::size_t row = 0; row < row_count - 1; ++row)
            {
               for (std::size_t index = 0; index < block_length; ++index)
               {
                  block_stack[row][index] = auxiliary_stack[row][index];
               }
            }

            for (std::size_t index = 0; index < partial_block_length; ++index)
            {
               block_stack[row_count - 1][index] = auxiliary_stack[row_count - 1][index];
            }
         }
      }

      template <typename T, std::size_t block_length>
      inline void interleave(data_block<T,block_length> block_stack[],
                             const std::size_t row_count)
      {
         data_block<T,block_length>* auxiliary_stack = new data_block<T,block_length>[row_count];

         std::size_t aux_row   = 0;
         std::size_t aux_index = 0;

         for (std::size_t index = 0; index < block_length; ++index)
         {
            for (std::size_t row = 0; row < row_count; ++row)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_index == block_length)
               {
                  aux_index = 0;
                  aux_row++;
               }
            }
         }

         for (std::size_t row = 0; row < row_count; ++row)
         {
            for (std::size_t index = 0; index < block_length; ++index)
            {
               block_stack[row][index] = auxiliary_stack[row][index];
            }
         }

         delete[] auxiliary_stack;
      }

      template <typename T, std::size_t block_length>
      inline void interleave(data_block<T,block_length> block_stack[],
                             const std::size_t row_count,
                             const std::size_t partial_block_length)
      {
         data_block<T,block_length>* auxiliary_stack = new data_block<T,block_length>[row_count];

         std::size_t aux_row   = 0;
         std::size_t aux_index = 0;

         for (std::size_t index = 0; index < partial_block_length; ++index)
         {
            for (std::size_t row = 0; row < row_count; ++row)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_index == block_length)
               {
                  aux_index = 0;
                  aux_row++;
               }
            }
         }

         for (std::size_t index = partial_block_length; index < block_length; ++index)
         {
            for (std::size_t row = 0; row < row_count - 1; ++row)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_index == block_length)
               {
                  aux_index = 0;
                  aux_row++;
               }
            }
         }

         for (std::size_t row = 0; row < row_count - 1; ++row)
         {
            for (std::size_t index = 0; index < block_length; ++index)
            {
               block_stack[row][index] = auxiliary_stack[row][index];
            }
         }

         for (std::size_t index = 0; index < partial_block_length; ++index)
         {
            block_stack[row_count - 1][index] = auxiliary_stack[row_count - 1][index];
         }

         delete[] auxiliary_stack;
      }

      template <std::size_t code_length, std::size_t fec_length, std::size_t row_count>
      inline void deinterleave(block<code_length,fec_length> (&block_stack)[row_count])
      {
         block<code_length,fec_length> auxiliary_stack[row_count];

         std::size_t aux_row   = 0;
         std::size_t aux_index = 0;

         for (std::size_t row = 0; row < row_count; ++row)
         {
            for (std::size_t index = 0; index < code_length; ++index)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_row == row_count)
               {
                  aux_row = 0;
                  aux_index++;
               }
            }
         }

         copy<code_length,fec_length,row_count>(auxiliary_stack,block_stack);
      }

      template <std::size_t code_length, std::size_t fec_length, std::size_t row_count>
      inline void deinterleave(block<code_length,fec_length> (&block_stack)[row_count],
                               const std::size_t partial_code_length)
      {
         if (partial_code_length == code_length)
         {
            deinterleave<code_length,fec_length,row_count>(block_stack);
         }
         else
         {
            block<code_length,fec_length> auxiliary_stack[row_count];

            std::size_t aux_row1   = 0;
            std::size_t aux_index1 = 0;

            std::size_t aux_row2   = 0;
            std::size_t aux_index2 = 0;

            for (std::size_t i = 0; i < partial_code_length * row_count; ++i)
            {
               auxiliary_stack[aux_row1][aux_index1] = block_stack[aux_row2][aux_index2];

               if (++aux_row1 == row_count)
               {
                  aux_row1 = 0;
                  aux_index1++;
               }

               if (++aux_index2 == code_length)
               {
                  aux_index2 = 0;
                  aux_row2++;
               }
            }

            for (std::size_t i = 0; aux_index1 < code_length; ++i)
            {
               auxiliary_stack[aux_row1][aux_index1] = block_stack[aux_row2][aux_index2];

               if (++aux_row1 == (row_count - 1))
               {
                  aux_row1 = 0;
                  aux_index1++;
               }

               if (++aux_index2 == code_length)
               {
                  aux_index2 = 0;
                  aux_row2++;
               }
            }

            for (std::size_t row = 0; row < row_count - 1; ++row)
            {
               for (std::size_t index = 0; index < code_length; ++index)
               {
                  block_stack[row][index] = auxiliary_stack[row][index];
               }
            }

            for (std::size_t index = 0; index < partial_code_length; ++index)
            {
               block_stack[row_count - 1][index] = auxiliary_stack[row_count - 1][index];
            }
         }
      }

      template <typename T, std::size_t block_length>
      inline void deinterleave(data_block<T,block_length> (&block_stack)[block_length])
      {
         data_block<T,block_length> auxiliary_stack[block_length];

         for (std::size_t row = 0; row < block_length; ++row)
         {
            for (std::size_t index = 0; index < block_length; ++index)
            {
               auxiliary_stack[index][row] = block_stack[row][index];
            }
         }

         copy<T,block_length,block_length>(auxiliary_stack,block_stack);
      }

      template <typename T, std::size_t block_length, std::size_t row_count>
      inline void deinterleave(data_block<T,block_length> (&block_stack)[row_count])
      {
         data_block<T,block_length> auxiliary_stack[row_count];

         std::size_t aux_row   = 0;
         std::size_t aux_index = 0;

         for (std::size_t row = 0; row < row_count; ++row)
         {
            for (std::size_t index = 0; index < block_length; ++index)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_row == row_count)
               {
                  aux_row = 0;
                  aux_index++;
               }
            }
         }

         copy<T,block_length,row_count>(auxiliary_stack,block_stack);
      }

      template <typename T, std::size_t block_length>
      inline void deinterleave(data_block<T,block_length> block_stack[],
                               const std::size_t row_count)
      {
         data_block<T,block_length>* auxiliary_stack = new data_block<T,block_length>[row_count];

         std::size_t aux_row   = 0;
         std::size_t aux_index = 0;

         for (std::size_t row = 0; row < row_count; ++row)
         {
            for (std::size_t index = 0; index < block_length; ++index)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_row == row_count)
               {
                  aux_row = 0;
                  aux_index++;
               }
            }
         }

         for (std::size_t row = 0; row < row_count; ++row)
         {
            for (std::size_t index = 0; index < block_length; ++index)
            {
               block_stack[row][index] = auxiliary_stack[row][index];
            }
         }

         delete[] auxiliary_stack;
      }

      template <typename T, std::size_t block_length>
      inline void deinterleave(data_block<T,block_length> block_stack[],
                               const std::size_t row_count,
                               const std::size_t partial_block_length)
      {
         if (row_count == 1) return;

         data_block<T,block_length>* auxiliary_stack = new data_block<T,block_length>[row_count];

         std::size_t aux_row1   = 0;
         std::size_t aux_index1 = 0;

         std::size_t aux_row2   = 0;
         std::size_t aux_index2 = 0;

         for (std::size_t i = 0; i < partial_block_length * row_count; ++i)
         {
            auxiliary_stack[aux_row1][aux_index1] = block_stack[aux_row2][aux_index2];

            if (++aux_row1 == row_count)
            {
               aux_row1 = 0;
               aux_index1++;
            }

            if (++aux_index2 == block_length)
            {
               aux_index2 = 0;
               aux_row2++;
            }
         }

         for (std::size_t i = 0; aux_index1 < block_length; ++i)
         {
            auxiliary_stack[aux_row1][aux_index1] = block_stack[aux_row2][aux_index2];

            if (++aux_row1 == (row_count - 1))
            {
               aux_row1 = 0;
               aux_index1++;
            }

            if (++aux_index2 == block_length)
            {
               aux_index2 = 0;
               aux_row2++;
            }
         }

         for (std::size_t row = 0; row < row_count - 1; ++row)
         {
            for (std::size_t index = 0; index < block_length; ++index)
            {
               block_stack[row][index] = auxiliary_stack[row][index];
            }
         }

         for (std::size_t index = 0; index < partial_block_length; ++index)
         {
            block_stack[row_count - 1][index] = auxiliary_stack[row_count - 1][index];
         }

         delete[] auxiliary_stack;
      }

      template <typename T, std::size_t block_length, std::size_t skip_columns>
      inline void interleave_columnskip(data_block<T,block_length>* block_stack)
      {
         for (std::size_t i = 0; i < block_length; ++i)
         {
            for (std::size_t j = i + 1; j < block_length; ++j)
            {
               std::size_t x1 = i + skip_columns;
               std::size_t x2 = j + skip_columns;

               T tmp = block_stack[i][x2];
               block_stack[i][x2] = block_stack[j][x1];
               block_stack[j][x1] = tmp;
            }
         }
      }

      template <typename T, std::size_t block_length, std::size_t skip_columns>
      inline void interleave_columnskip(data_block<T,block_length>* block_stack, const std::size_t& row_count)
      {
         data_block<T,block_length>* auxiliary_stack = new data_block<T,block_length>[row_count];

         std::size_t aux_row   = 0;
         std::size_t aux_index = skip_columns;

         for (std::size_t index = skip_columns; index < block_length; ++index)
         {
            for (std::size_t row = 0; row < row_count; ++row)
            {
               auxiliary_stack[aux_row][aux_index] = block_stack[row][index];

               if (++aux_index == block_length)
               {
                  aux_index = skip_columns;
                  aux_row++;
               }
            }
         }

         for (std::size_t row = 0; row < row_count; ++row)
         {
            for (std::size_t index = skip_columns; index < block_length; ++index)
            {
               block_stack[row][index] = auxiliary_stack[row][index];
            }
         }

         delete[] auxiliary_stack;
      }

      template <typename T, std::size_t data_length>
      inline void interleave(T* block_stack[data_length])
      {
         for (std::size_t i = 0; i < data_length; ++i)
         {
            for (std::size_t j = i + 1; j < data_length; ++j)
            {
               T tmp = block_stack[i][j];
               block_stack[i][j] = block_stack[j][i];
               block_stack[j][i] = tmp;
            }
         }
      }

      template <typename T, std::size_t data_length, std::size_t skip_columns>
      inline void interleave_columnskip(T* block_stack[data_length])
      {
         for (std::size_t i = skip_columns; i < data_length; ++i)
         {
            for (std::size_t j = i + 1; j < data_length; ++j)
            {
               T tmp = block_stack[i][j];
               block_stack[i][j] = block_stack[j][i];
               block_stack[j][i] = tmp;
            }
         }
      }

   } // namespace reed_solomon

} // namespace schifra

#endif
