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


#ifndef INCLUDE_SCHIFRA_FILEIO_HPP
#define INCLUDE_SCHIFRA_FILEIO_HPP


#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>

#include "schifra_crc.hpp"


namespace schifra
{

   namespace fileio
   {

      inline void read_into_vector(const std::string& file_name, std::vector<std::string>& buffer)
      {
         std::ifstream file(file_name.c_str());
         if (!file) return;
         std::string line;
         while (std::getline(file,line))
         {
            buffer.push_back(line);
         }
         file.close();
      }

      inline void write_from_vector(const std::string& file_name, const std::vector<std::string>& buffer)
      {
         std::ofstream file(file_name.c_str());
         if (!file) return;
         std::ostream_iterator <std::string> os(file,"\n");
         std::copy(buffer.begin(),buffer.end(), os);
         file.close();
      }

      inline bool file_exists(const std::string& file_name)
      {
         std::ifstream file(file_name.c_str(), std::ios::binary);
         return ((!file) ? false : true);
      }

      inline std::size_t file_size(const std::string& file_name)
      {
         std::ifstream file(file_name.c_str(),std::ios::binary);
         if (!file) return 0;
         file.seekg (0, std::ios::end);
         return static_cast<std::size_t>(file.tellg());
      }

      inline void load_file(const std::string& file_name, std::string& buffer)
      {
         std::ifstream file(file_name.c_str(), std::ios::binary);
         if (!file) return;
         buffer.assign(std::istreambuf_iterator<char>(file),std::istreambuf_iterator<char>());
         file.close();
      }

      inline void load_file(const std::string& file_name, char** buffer, std::size_t& buffer_size)
      {
         std::ifstream in_stream(file_name.c_str(),std::ios::binary);
         if (!in_stream) return;
         buffer_size = file_size(file_name);
         *buffer = new char[buffer_size];
         in_stream.read(*buffer,static_cast<std::streamsize>(buffer_size));
         in_stream.close();
      }

      inline void write_file(const std::string& file_name, const std::string& buffer)
      {
         std::ofstream file(file_name.c_str(),std::ios::binary);
         file << buffer;
         file.close();
      }

      inline void write_file(const std::string& file_name, char* buffer, const std::size_t& buffer_size)
      {
         std::ofstream out_stream(file_name.c_str(),std::ios::binary);
         if (!out_stream) return;
         out_stream.write(buffer,static_cast<std::streamsize>(buffer_size));
         out_stream.close();
      }

      inline bool copy_file(const std::string& src_file_name, const std::string& dest_file_name)
      {
         std::ifstream src_file(src_file_name.c_str(),std::ios::binary);
         std::ofstream dest_file(dest_file_name.c_str(),std::ios::binary);
         if (!src_file) return false;
         if (!dest_file) return false;

         const std::size_t block_size = 1024;
         char buffer[block_size];

         std::size_t remaining_bytes = file_size(src_file_name);

         while (remaining_bytes >= block_size)
         {
            src_file.read(&buffer[0],static_cast<std::streamsize>(block_size));
            dest_file.write(&buffer[0],static_cast<std::streamsize>(block_size));
            remaining_bytes -= block_size;
         }

         if (remaining_bytes > 0)
         {
            src_file.read(&buffer[0],static_cast<std::streamsize>(remaining_bytes));
            dest_file.write(&buffer[0],static_cast<std::streamsize>(remaining_bytes));
            remaining_bytes = 0;
         }

         src_file.close();
         dest_file.close();

         return true;
      }

      inline bool files_identical(const std::string& file_name1, const std::string& file_name2)
      {
         std::ifstream file1(file_name1.c_str(),std::ios::binary);
         std::ifstream file2(file_name2.c_str(),std::ios::binary);
         if (!file1) return false;
         if (!file2) return false;
         if (file_size(file_name1) != file_size(file_name2)) return false;

         const std::size_t block_size = 1024;
         char buffer1[block_size];
         char buffer2[block_size];

         std::size_t remaining_bytes = file_size(file_name1);

         while (remaining_bytes >= block_size)
         {
            file1.read(&buffer1[0],static_cast<std::streamsize>(block_size));
            file2.read(&buffer2[0],static_cast<std::streamsize>(block_size));

            for (std::size_t i = 0; i < block_size; ++i)
            {
               if (buffer1[i] != buffer2[i])
               {
                  return false;
               }
            }

            remaining_bytes -= block_size;
         }

         if (remaining_bytes > 0)
         {
            file1.read(&buffer1[0],static_cast<std::streamsize>(remaining_bytes));
            file2.read(&buffer2[0],static_cast<std::streamsize>(remaining_bytes));

            for (std::size_t i = 0; i < remaining_bytes; ++i)
            {
               if (buffer1[i] != buffer2[i])
               {
                  return false;
               }
            }

            remaining_bytes = 0;
         }

         file1.close();
         file2.close();

         return true;
      }

      inline std::size_t file_crc(crc32& crc_module, const std::string& file_name)
      {
         std::ifstream file(file_name.c_str(),std::ios::binary);
         if (!file) return 0;

         const std::size_t block_size = 1024;
         char buffer[block_size];

         std::size_t remaining_bytes = file_size(file_name);

         crc_module.reset();

         while (remaining_bytes >= block_size)
         {
            file.read(&buffer[0],static_cast<std::streamsize>(block_size));
            crc_module.update(buffer,block_size);
            remaining_bytes -= block_size;
         }

         if (remaining_bytes > 0)
         {
            file.read(&buffer[0],static_cast<std::streamsize>(remaining_bytes));
            crc_module.update(buffer,remaining_bytes);
            remaining_bytes = 0;
         }

         return crc_module.crc();
      }

   } // namespace fileio

} // namespace schifra

#endif
