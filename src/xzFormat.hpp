/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of YAX.

    YAX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must 
    be included in all copies or substantial portions of YAX.
*/


#include <iostream>
#include <thread>
#include <algorithm>
#include <unistd.h>
#include <fcntl.h>
#include <lzma.h>
#include <string.h>


// Class to record recent int-value function calls
//    for FIFO-type caching / temporary memoization.
class call_tracker
{

	private:
		std::vector<int> most_recent;
		int n;
	public:
		call_tracker(const int& size) : n(size) {};
		call_tracker() {};
		void set_size(const int& n_recent){
			n = n_recent;
		}

		int push_new(const int& x){
			most_recent.push_back(x);
			if( most_recent.size() < n ){
				// Don't remove anything from cache. 
				return -1;
			}
			int out = most_recent.front();
			most_recent.erase(most_recent.begin());
			return out;
		}
		
		int check_new( const int& x ){
			if( most_recent.size() < n ){
				most_recent.push_back(x);
				return -1;
			}
			int x_i = -1;
			for(int i = n; i >= 0; i-- ){
				if( most_recent[i] == x ){
					x_i = i;
					break;
				}
			}
			if( x_i > 0 ){
				std::swap(most_recent[most_recent.size()-1], most_recent[x_i]);
				return -1;
			}else{
				most_recent.push_back(x);
				int out = most_recent.front();
				most_recent.erase(most_recent.begin());
				return out;
			}
		}	
};



// The code for the xzReader class is partly  
// derived from the following GitHub repository:
//     https://github.com/CTSRD-CHERI/cheritrace.git
// The license issued by the original authors 
// (Alfredo Mazzinghi and David T. Chisnall) for this code is here: 
//     http://www.beri-open-systems.org/legal/license-1-0.txt

class xzReader 
{
	protected:
		struct offsets
		{
			off_t c_start; size_t c_size;
			off_t u_start; size_t u_size;
		};
		size_t read_compressed(void *buffer, off_t start, size_t n_bytes){
			if ( start < 0 || start > compressed_file_size ){
				return 0;
			}
			if (start + n_bytes > compressed_file_size){
				n_bytes = compressed_file_size - start;
			}
			size_t n_completed = 0;
			while (n_bytes > 0)
			{
				ssize_t n_bytes_i = pread(compressed_file, buffer, n_bytes, start);
				if (n_bytes_i < 0) break;
				n_completed += n_bytes_i;
				start += n_bytes_i;
				n_bytes -= n_bytes_i;
				buffer = (void*)((char*)buffer + n_bytes_i);
			}
			return n_completed;
		}
		
		call_tracker cache_czar;
		std::vector<std::unique_ptr<uint8_t>> cache;
		
		int compressed_file;
		
		std::vector<offsets> block_offsets;
		size_t compressed_file_size = 0;

		lzma_stream_flags stream_flags;

	public:
		
		void open(const std::string& filepath)
		{

			compressed_file = ::open(filepath.c_str(), O_RDONLY);
			compressed_file_size = lseek(compressed_file, 0, SEEK_END);

			uint8_t buffer[12];
			read_compressed((void*)buffer, compressed_file_size-12, 12);
			
			if ( lzma_stream_footer_decode(&stream_flags, buffer) != LZMA_OK ){
				std::cerr << "ERROR: Cannot read xz footer.\n";
				close(compressed_file);
				abort();
			}
			
			std::unique_ptr<uint8_t> index_buffer(new uint8_t[stream_flags.backward_size]);
			read_compressed((void*)index_buffer.get(),
								  compressed_file_size - stream_flags.backward_size - 12,
								  stream_flags.backward_size);
			lzma_index *idx;
			uint64_t mem = UINT64_MAX;
			size_t pos = 0;
			
			if ( lzma_index_buffer_decode(&idx, &mem, nullptr,index_buffer.get(), &pos, stream_flags.backward_size) != LZMA_OK ){
				std::cerr << "ERROR: Cannot read xz index.\n";
				close(compressed_file);
				abort();
			}
			lzma_index_iter iter;
			lzma_index_iter_init(&iter, idx);
			while (!lzma_index_iter_next(&iter, LZMA_INDEX_ITER_ANY))
			{
				struct offsets block;
				block.c_start = iter.block.compressed_file_offset;
				block.c_size = iter.block.total_size;
				block.u_start = iter.block.uncompressed_file_offset;
				block.u_size = iter.block.uncompressed_size;
				block_offsets.push_back(block);
			}
			lzma_index_end(idx, nullptr);
			offsets &b = block_offsets.back();
			
			cache_czar.set_size(100);
			cache.resize(block_offsets.size());
		}
		xzReader(const std::string& filepath){
			open(filepath);
		}
		xzReader() {};
		
		size_t read(void *buffer, off_t start, size_t length) 
		{
			int block_idx = get_block_from_offset(start);
			if (block_idx < 0)
			{
				return 0;
			}
			size_t copied = 0;
			while (length > 0)
			{
				if (block_idx >= (int)block_offsets.size())
				{
					break;
				}
				if( cache_block(block_idx) < 0 ){
					std::cerr << "ERROR: Cannot read block " << block_idx << "\n";
				}
				const auto& data = cache[block_idx];
				auto &b = block_offsets[block_idx++];
				size_t copy_start = start - b.u_start;
				size_t copy_length = b.u_size - copy_start;
				copy_length = std::min(copy_length, length);
				memcpy(buffer, data.get()+copy_start, copy_length);
				copied += copy_length;
				start += copy_length;
				length -= copy_length;
				buffer = (void*)((char*)buffer + copy_length);
			}
			return copied;
		}
		
		int get_block_from_offset(off_t off)
		{
			int i = 0;
			for(const offsets& bi : block_offsets){
				if( (bi.u_start <= off) && (bi.u_start + bi.u_size > off) ){
					return i;
				}
				i++;
			}
			return -1;
		}

		int cache_block(const int& i)
		{
			if( cache[i] != nullptr ){
				return 1;
			}else{
				int prune_cache = cache_czar.push_new(i);
				if( prune_cache > 0 ){
					cache[prune_cache] = nullptr;
				}
			}
			
			offsets& b = block_offsets[i];
			
			std::unique_ptr<uint8_t> buffer(new uint8_t[b.c_size]);
			
			read_compressed((void*)buffer.get(), b.c_start, b.c_size);
			lzma_block block;
			lzma_filter filters[LZMA_FILTERS_MAX + 1];
			
			filters[0].id = LZMA_VLI_UNKNOWN;
			block.filters = filters;
			block.version = 1;
			block.check = stream_flags.check;
			block.header_size = lzma_block_header_size_decode(*buffer);
			
			if ( lzma_block_header_decode(&block, nullptr, buffer.get()) != LZMA_OK )
			{
				return -1;
			}
			
			cache[i] = std::unique_ptr<uint8_t>(new uint8_t[b.u_size]);
			size_t in_pos = block.header_size;
			size_t out_pos = 0;
			
			if ( lzma_block_buffer_decode(&block, nullptr, buffer.get(),&in_pos, b.c_size, cache[i].get(), &out_pos, b.u_size) != LZMA_OK)
			{
				return -1;
			}
			
			return 1; 
		}
};


// Class for writing xz files from output stream. 
// We want to pass small chunks of data at a time 
// (for example, results calculated for a given gene or 
//  genomic region) to be written to file via writeBytes, 
// and take care of blocking, indexing, and other
// details internally. 

// Uncompressed data is kept in buffer. Once buffer reaches
// capacity based on the specified bytes_per_block, we 
// compress the data in buffer and write a new xz block
// to the file.  Then, we clear buffer and continue. 

// At the end, we want to write any remaining data as the final 
// xz block, and then write the xz block index to file. 

/*
class xzWriter
{
	protected:
		// Buffer storing uncompressed data.
		void* buffer;
		// How many bytes have been copied into the buffer.
		size_t buffer_n;
		
		// How many (uncompressed) bytes will be stored per block.
		size_t bytes_per_block;
		
		// Output file descriptor. 
		int xz_out_file;
		
		int writeStreamHeader();
		int writeBlock();
		int writeIndex();
		int writeStreamFooter();
	public:
		xzWriter() {};
		xzWriter( const std::string& fp, const int& nb) { open(fp, nb); };
		int open(const std::string& filepath, const int& nb)
		{
			bytes_per_block = nb;
			buffer = malloc(bytes_per_block);
			xz_out_file = ::open(filepath.c_str(), O_WRONLY);
			writeStreamHeader();
		}
		int writeBytes(void* data, size_t input_size){
			// // Copy data to buffer.
			// size_t bytes_written = 0;
			// while( bytes_written < input_size ){
			//   if( input_size + buffer_n >= bytes_per_block ){
			//     // A. Copy first (bytes_per_block - buffer_n)
			//     //    bytes to buffer. Update the counts.
			//         bytes_written += bytes_per_block - buffer_n;
			//         buffer_n = bytes_per_block;
			//          
			//     // B. Compress buffer and write block to file.
			//         writeBlock();
			//     // C. Clear buffer. Update buffer_n. 
			//         buffer_n = 0;
			//   }else{
            //     // Copy data to buffer. Update the counts.
			//	       bytes_written += input_size;
			//         buffer_n += input_size;
			//   }
			// } 
			return 0;
		}
		int close(){
			// Write any remaining data to final block. 
			writeBlock();
			// Free the buffer. 
			free(buffer);
			// Write the XZ index and footer. 
			writeIndex();
			writeStreamFooter();
			// Close the file. 
			if( ::close(xz_out_file) < 0 ){
				std::cerr << "ERROR: Failed to close xz file.\n";
			}
			return 0;
		}
};
*/

