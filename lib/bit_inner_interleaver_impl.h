/* -*- c++ -*- */
/* 
 * Copyright 2013 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifndef INCLUDED_DVBT_BIT_INNER_INTERLEAVER_IMPL_H
#define INCLUDED_DVBT_BIT_INNER_INTERLEAVER_IMPL_H

#include <dvbt/bit_inner_interleaver.h>
#include <dvbt/dvbt_config.h>

namespace gr {
  namespace dvbt {

    class bit_inner_interleaver_impl : public bit_inner_interleaver
    {
    private:
      const dvbt_config config;

      // constellation
      unsigned char d_v;
      // Bit interleaver block size
      static const unsigned char d_bsize;

      //Table to keep interleaved indices
      unsigned char * d_perm;

      // permutation function
      unsigned char H(unsigned char e, unsigned char w);

    public:
      bit_inner_interleaver_impl(int ninput, int noutput, \
        dvbt_constellation_t constellation, dvbt_hierarchy_t hierarchy);
      ~bit_inner_interleaver_impl();

      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      // Where all the action really happens
      int general_work(int noutput_items,
		       gr_vector_int &ninput_items,
		       gr_vector_const_void_star &input_items,
		       gr_vector_void_star &output_items);
    };

  } // namespace dvbt
} // namespace gr

#endif /* INCLUDED_DVBT_BIT_INNER_INTERLEAVER_IMPL_H */

