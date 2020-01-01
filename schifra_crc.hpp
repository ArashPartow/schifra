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


#ifndef INCLUDE_SCHIFRA_CRC_HPP
#define INCLUDE_SCHIFRA_CRC_HPP


#include <iostream>
#include <string>


namespace schifra
{

   class crc32
   {
   public:

      typedef std::size_t crc32_t;

      crc32(const crc32_t& _key, const crc32_t& _state = 0x00)
      : key(_key),
        state(_state),
        initial_state(_state)
      {
         initialize_crc32_table();
      }

      void reset()
      {
         state = initial_state;
      }

      void update_1byte(const unsigned char data)
      {
         state = (state >> 8) ^ table[data];
      }

      void update(const unsigned char data[], const std::size_t& count)
      {
         for (std::size_t i = 0; i < count; ++i)
         {
            update_1byte(data[i]);
         }
      }

      void update(char data[], const std::size_t& count)
      {
         for (std::size_t i = 0; i < count; ++i)
         {
            update_1byte(static_cast<unsigned char>(data[i]));
         }
      }

      void update(const std::string& data)
      {
         for (std::size_t i = 0; i < data.size(); ++i)
         {
            update_1byte(static_cast<unsigned char>(data[i]));
         }
      }

      void update(const std::size_t& data)
      {
         update_1byte(static_cast<unsigned char>((data      ) & 0xFF));
         update_1byte(static_cast<unsigned char>((data >>  8) & 0xFF));
         update_1byte(static_cast<unsigned char>((data >> 16) & 0xFF));
         update_1byte(static_cast<unsigned char>((data >> 24) & 0xFF));
      }

      crc32_t crc()
      {
         return state;
      }

   private:

      crc32& operator=(const crc32&);

      void initialize_crc32_table()
      {
         for (std::size_t i = 0; i < 0xFF; ++i)
         {
            crc32_t reg = i;

            for (int j = 0; j < 0x08; ++j)
            {
               reg = ((reg & 1) ? (reg >> 1) ^ key : reg >> 1);
            }

            table[i] = reg;
         }
      }

   protected:

      crc32_t key;
      crc32_t state;
      const crc32_t initial_state;
      crc32_t table[256];
   };

   class schifra_crc : public crc32
   {
   public:

      schifra_crc(const crc32_t _key)
      : crc32(_key,0xAAAAAAAA)
      {}

      void update(const unsigned char& data)
      {
         state = ((state >> 8) ^ table[data]) ^ ((state << 8) ^ table[~data]);
      }

      void update(const unsigned char data[], const std::size_t& count)
      {
         for (std::size_t i = 0; i < count; ++i)
         {
            update_1byte(data[i]);
         }
      }

      void update(const char data[], const std::size_t& count)
      {
         for (std::size_t i = 0; i < count; ++i)
         {
            update_1byte(static_cast<unsigned char>(data[i]));
         }
      }

      void update(const std::string& data)
      {
         for (std::size_t i = 0; i < data.size(); ++i)
         {
            update_1byte(static_cast<unsigned char>(data[i]));
         }
      }

      void update(const std::size_t& data)
      {
         update_1byte(static_cast<unsigned char>((data      ) & 0xFF));
         update_1byte(static_cast<unsigned char>((data >>  8) & 0xFF));
         update_1byte(static_cast<unsigned char>((data >> 16) & 0xFF));
         update_1byte(static_cast<unsigned char>((data >> 24) & 0xFF));
      }

   };

} // namespace schifra


#endif
