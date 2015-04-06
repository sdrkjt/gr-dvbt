#!/usr/bin/env python
##################################################
# Gnuradio Python Flow Graph
# Title: Top Block
# Generated: Mon Apr  6 10:07:39 2015
##################################################

from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import fft
from gnuradio import filter
from gnuradio import gr
from gnuradio import uhd
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from gnuradio.wxgui import fftsink2
from gnuradio.wxgui import scopesink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import dvbt
import time
import wx

class top_block(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="Top Block")

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 32000

        ##################################################
        # Blocks
        ##################################################
        self.wxgui_scopesink2_0 = scopesink2.scope_sink_c(
        	self.GetWin(),
        	title="Scope Plot",
        	sample_rate=samp_rate,
        	v_scale=0,
        	v_offset=0,
        	t_scale=0,
        	ac_couple=False,
        	xy_mode=True,
        	num_inputs=1,
        	trig_mode=wxgui.TRIG_MODE_AUTO,
        	y_axis_label="Counts",
        )
        self.Add(self.wxgui_scopesink2_0.win)
        self.wxgui_fftsink2_0 = fftsink2.fft_sink_c(
        	self.GetWin(),
        	baseband_freq=570e6,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=10e6,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title="FFT Plot",
        	peak_hold=False,
        )
        self.Add(self.wxgui_fftsink2_0.win)
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	",".join(("addr=192.168.10.2", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_source_0.set_samp_rate(10e6)
        self.uhd_usrp_source_0.set_center_freq(1280e6, 0)
        self.uhd_usrp_source_0.set_gain(19, 0)
        self.rational_resampler_xxx_0 = filter.rational_resampler_ccc(
                interpolation=64,
                decimation=70,
                taps=None,
                fractional_bw=None,
        )
        self.fft_vxx_0 = fft.fft_vcc(2048, True, (window.rectangular(2048)), True, 1)
        self.dvbt_viterbi_decoder_0 = dvbt.viterbi_decoder(dvbt.QAM16, dvbt.NH, dvbt.C1_2, 512, 0, -1)
        self.dvbt_symbol_inner_interleaver_0 = dvbt.symbol_inner_interleaver(1512, dvbt.T2k, 0)
        self.dvbt_reed_solomon_dec_0 = dvbt.reed_solomon_dec(2, 8, 0x11d, 255, 239, 8, 51, 8)
        self.dvbt_ofdm_sym_acquisition_0 = dvbt.ofdm_sym_acquisition(1, 2048, 1705, 64, 30)
        self.dvbt_energy_descramble_0 = dvbt.energy_descramble(8)
        self.dvbt_dvbt_demap_0 = dvbt.dvbt_demap(1512, dvbt.QAM16, dvbt.NH, dvbt.T2k, 1)
        self.dvbt_demod_reference_signals_0 = dvbt.demod_reference_signals(gr.sizeof_gr_complex, 2048, 1512, dvbt.QAM16, dvbt.NH, dvbt.C1_2, dvbt.C1_2, dvbt.G1_32, dvbt.T2k, 0, 0)
        self.dvbt_convolutional_deinterleaver_0 = dvbt.convolutional_deinterleaver(136, 12, 17)
        self.dvbt_bit_inner_deinterleaver_0 = dvbt.bit_inner_deinterleaver(1512, dvbt.QAM16, dvbt.NH, dvbt.T2k)
        self.blocks_vector_to_stream_1 = blocks.vector_to_stream(gr.sizeof_char*1, 1512)
        self.blocks_vector_to_stream_0 = blocks.vector_to_stream(gr.sizeof_gr_complex*1, 1512)
        self.blocks_multiply_const_vxx_0 = blocks.multiply_const_vcc((0.0022097087, ))
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_char*1, "test_out.ts", False)
        self.blocks_file_sink_0.set_unbuffered(False)
        self.blocks_add_xx_0 = blocks.add_vcc(1)
        self.analog_noise_source_x_0 = analog.noise_source_c(analog.GR_GAUSSIAN, 0.000, 0)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_noise_source_x_0, 0), (self.blocks_add_xx_0, 1))    
        self.connect((self.blocks_add_xx_0, 0), (self.rational_resampler_xxx_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.dvbt_ofdm_sym_acquisition_0, 0))    
        self.connect((self.blocks_multiply_const_vxx_0, 0), (self.wxgui_fftsink2_0, 0))    
        self.connect((self.blocks_vector_to_stream_0, 0), (self.wxgui_scopesink2_0, 0))    
        self.connect((self.blocks_vector_to_stream_1, 0), (self.dvbt_viterbi_decoder_0, 0))    
        self.connect((self.dvbt_bit_inner_deinterleaver_0, 0), (self.blocks_vector_to_stream_1, 0))    
        self.connect((self.dvbt_convolutional_deinterleaver_0, 0), (self.dvbt_reed_solomon_dec_0, 0))    
        self.connect((self.dvbt_demod_reference_signals_0, 0), (self.blocks_vector_to_stream_0, 0))    
        self.connect((self.dvbt_demod_reference_signals_0, 0), (self.dvbt_dvbt_demap_0, 0))    
        self.connect((self.dvbt_dvbt_demap_0, 0), (self.dvbt_symbol_inner_interleaver_0, 0))    
        self.connect((self.dvbt_energy_descramble_0, 0), (self.blocks_file_sink_0, 0))    
        self.connect((self.dvbt_ofdm_sym_acquisition_0, 0), (self.fft_vxx_0, 0))    
        self.connect((self.dvbt_reed_solomon_dec_0, 0), (self.dvbt_energy_descramble_0, 0))    
        self.connect((self.dvbt_symbol_inner_interleaver_0, 0), (self.dvbt_bit_inner_deinterleaver_0, 0))    
        self.connect((self.dvbt_viterbi_decoder_0, 0), (self.dvbt_convolutional_deinterleaver_0, 0))    
        self.connect((self.fft_vxx_0, 0), (self.dvbt_demod_reference_signals_0, 0))    
        self.connect((self.rational_resampler_xxx_0, 0), (self.blocks_multiply_const_vxx_0, 0))    
        self.connect((self.uhd_usrp_source_0, 0), (self.blocks_add_xx_0, 0))    


    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.wxgui_scopesink2_0.set_sample_rate(self.samp_rate)

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"
    parser = OptionParser(option_class=eng_option, usage="%prog: [options]")
    (options, args) = parser.parse_args()
    tb = top_block()
    tb.Start(True)
    tb.Wait()
