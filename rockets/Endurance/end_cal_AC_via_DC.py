"""
Calibrate the AC data to exact mV/m based on the DC data in the overlap region around 10 Hz.

From Paulo: The DC data have a very accurate calibration from counts to mV/m. The AC data are not accurate.
Here we use the fact that the two channels overlap in frequency around 10 Hz to better calibrate the AC data.

"""
import sys
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
#import plot_spectrogram as ps




