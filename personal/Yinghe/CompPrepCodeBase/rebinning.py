#from numarray import *
import numpy as np

def womply(a, i, x):
#
#             I+1
#  EVALUATES  SUM A(J) * X**(I+1-J)
#             J=1
#
#  I IS THUS THE ORDER OF THE POLY, WHICH HAS I+1 COEFFICIENTS,
#  AND A(I+1) IS THE CONSTANT TERM.
#
#   n_params = 3
   
   if (i < 0):   
      womply = 0.0e0
      print '*** VALUE OF I OUT OF BOUNDS IN WOMPLY ***'
      return womply
   
   womply = a[0]
   if (i == 0):   
      return womply
   for j in np.arange(2, (i + 1)+(1)):
      womply = x * womply + a[j - 1]
   return womply
   
"""
##########################################
##########################################
NOTE:
##########################################
##########################################
This file created by putting "foo.pro" into the same folder as the file "idl2python".
while in the folder containing both "foo.pro" and "idl2python", simply type into the 
command line (not in python): "idl2python foo.pro"
doing "idl2python -d foo.pro" will print out all the code in the terminal, leaving out
the "-d" creates a file "foo.py"

HOPEFULLY THIS WORKS...


Original IDL code for following code block
      N=1

      RL=WOMPLY(ARCTWO,IN21,RX)
;print, rl
jump1:     L=WOMPLY(ARCONE,IN11,X)
;       print, 'l ', l
      DX=(RL-L)*RSLOPE
      X=X+DX
;      print, 'x ', x
      IF(ABS(DX) LT ACC) then RETURN
      IF(N GT 100) then RETURN
      N=N+1
      RSLOPE=1/WOMPLY(DARCONE,IN12,X)
      GOTO, jump1
      END

idl2python can't convert the "GOTO" command so I had to conver this snippet by myself.
I think what they were trying to do was a loop until one of the if statement conditions
were met, so i removed the GOTO line so that it would convert the whole code, and added
a while loop that wraps where the GOTO statement would go
"""

def womrebin_tform(rx, x, arctwo, darcone, in21, acc, rslope, in11, in12, arcone):
#
# COPYRIGHT 1982 M. ASHLEY. See comments with SUBROUTINE REBIN.
#
# Version 1.2     23 Apr 1982    - More arguments put into COMMON.
# Version 1.1     29 Mar 1982    - First working version
#

#   n_params = 10
   def _ret():  return (rx, x, arctwo, darcone, in21, acc, rslope, in11, in12, arcone)
   n = 1          
   rl = womply(arctwo, in21, rx)
#   print rl
   # jump1:
   while(x>0):
       l = womply(arcone, in11, x)
       print 'l ', l
       dx = (rl - l) * rslope
       x = x + dx
       print 'x ', x
       if (abs(dx) < acc):   
           return _ret()
       if (n > 100):   
         return _ret()
         exit
       else:
           n = n + 1
           rslope = 1 / womply(darcone, in12, x)


def womashrebin(wave, flux, nwave, nflux):

#  COMMENTS BELOW FROM ASHREBIN OF SARAH AND LOLITA
#  I have translated this into IDL
#  TM  2000-05-12
#***************************************************************
#    Please do not make changes to this program without documenting them below.
#
#
# Version 1.4      8 Jul 1982 - Minor bug fixed (NSGN had the wrong value if
#                               you were rebinning from decreasing to increasing
#                               wavelength scales (or vice versa) and if, at the
#                               same time the slopes were less than a certain
#                               small number).
# Version 1.3     23 Apr 1982 - Documentation added; some code changes.
# Version 1.2     22 Apr 1982 - MACRO PLY replaced by FORTRAN version.
# Version 1.1     29 Mar 1982 - First "working version".
#
#
# COPYRIGHT 1982 M. ASHLEY.  Use of this program at this stage of
# development is entirely the user's responsibility. Unauthorized
# use or copying of this program is strictly prohibited.
#
#
# HISTORY OF DEVELOPMENT
# ----------------------
#
#      This subroutine was written by Michael Ashley as a fast
# replacement for Peter Young's subroutine of the same name, which had
# been the bottleneck when running the data reduction package LOLITA.
#
#
# QUICK SUMMARY OF OPERATION
# --------------------------
#
#      REBIN will rebin data from the first NPIX entries of array WDATA
# into the array RBIN. If IQUAD equals zero the rebinning will be done
# using linear interpolation, otherwise quadratic interpolation will be
# used. The program can be used in two modes: if IMODE equals zero then
# the output bins will be shifted by SSKEW relative to the input bins,
# and each output bin will consis of the sum of NRBIN input bins; otherwise
# the rebinnning will be done from the wavelength scale ARC1 to the scale
# ARC2, with NRBIN now meaning the number of output bins. In this later
# mode the absolute magnitudes of NARC1 and NARC2 are the number of
# coefficients in ARC1 and ARC2 respectively, and their signs indicate
# whether logarithmic scales are to be used.
#
#
# DESCRIPTION OF ARGUMENTS
# ------------------------
#
#   IMODE   - INTEGER -  If equal to zero, selects mode 0, in which input
#                        bins can be grouped and shifted to form the output
#                        bins.
#
#                        If not equal to zero, selects mode 1, in which the
#                        input bins are tranferred to the output bins in
#                        accordance with the wavelength scales ARC1 and ARC2.
#
#   IQUAD   - INTEGER -  If equal to zero, selects linear interpolation when
#                        rebinning.
#
#                        If not equal to zero, selects quadratic interpolation.
#
#   WDATA   - REAL*4  -  Array containing input data.
#
#   NPIX    - INTEGER -  The number of bins (elements) in WDATA.
#
#   RBIN    - REAL*4  -  Output Array, into which the rebinned data will be
#                        placed.
#
#   NRBIN   - INTEGER -  MODE 0: The number of input bins to be added together
#                                to form one output bin.
#
#                        MODE 1: The number of bins in RBIN.
#
#   SSKEW   - REAL*4  -  The number of bins the input array is to shifted (used
#                        in mode 0 only).
#
#   NARC1   - INTEGER -  Used in MODE 1 only. The absolute value of NARC1 is
#                        the number of coefficients in ARC1. If NARC1 is
#                        negative, then a logarithmic scale is assumed.
#
#   ARC1    - REAL*8  -  Used in MODE 1 only. The array of polynomial coeff-
#                        icients defining the wavelength scale of the input
#                        data.
#
#   NARC2   - INTEGER -  As for NARC1 except refers to the output wavelength
#                        scale ARC2.
#
#   ARC2    - REAL*8  -  As for ARC2 except it is the output wavelength scale.
#
#
#
# POTENTIAL PROBLEMS
# ------------------
#
#      REBIN assumes a certain degree of responsibility on the part of the
# program that calls it. Only minimal checking is done to ensure sensible
# arguments. Needless to say, if REBIN crashes, a large amount of time can
# be wasted. It is up to the user to ensure that inappropriate arguments
# can never be sent to REBIN.
#
#      REBIN may crash in the following ways -
#
# [1] Floating point overflow in REBIN_POLY, due to invalid arc coefficients.
#
# [2] Floating point overflow, or argument out of range in WOMREBIN_TFORM when
#     taking the logarithm and/or exponent of a result returned from
#     REBIN_POLY. This would be caused by having invalid arc coefficients.
#
# [3] In the Newton Raphson solution for output bin edges in terms of the
#     input bins it is assumed that the wavelength scales are fairly smooth,
#     and single valued over the interval of interest. If bizarre coefficients
#     are used, all sorts of difficulties may occur.
#
# [4] Divide by zero in REBIN and/or WOMREBIN_TFORM if the derivative of the
#     ARC1 polynomial equals zero at a tested point.
#
# [5] The user must ensure tha the arrays WDATA and RBIN are at least as
#     large as NPIX and NRBIN respectively (except in the case where IMODE
#     equals zero, when RBIN must be as large as necessary to fit the
#     rebinned data).
#
# [6] ARC1 and ARC2 must be dimensioned at least as large as the absolute
#     values of NARC1 and NARC2 respectively.
#
#
#
# LINEAR/QUADRATIC REBINNING
# --------------------------
#
#      IQUAD is used as a flag to set either linear or quadratic rebinning.
# This refers to the way in which the input data is interpolated to arrive
# at the number of counts in a fraction of a bin. Suppose that one of the
# output bins completely covers one of the input bins, and partially covers
# the bins on either side. Bin <n> is treated as though it collects all the
# counts from <n-0.5> to <n+0.5>. It is not treated as sampling the data
# a position <n> only. In linear rebinning it is assumed that the count
# rate was uniform across each bin. Thus, for example, if an output bin covers
# one quarter of an input bin, then it receives one quarter of the counts
# in that bin.
#
#      In quadratic rebinning it is assumed that, locally, the original data
# followed a parabola, and that the number of counts in bin <n> results
# from integrating this parabola from <n-0.5> to <n+0.5>. When a fraction
# of an input bin is required, the equation of the parabola which fits this
# bin and the neighbouring bins is found, and then integrated appropriately.
#
#      In both linear and quadratic interpolation, after the number of counts
# in each output bin has been found, these are divided by the width of
# each output bin in terms of input bins. This ensures that the mean level
# of the data doess not change. Hence, after rebinning it is not necessarily
# valid to say that the noise in a bin is the square root of the number of
# counts in that bin.
#
#
#
#print, size(nwave)
#   n_params = 4
   def _ret():  return (wave, flux, nwave, nflux)
   
   wdata = np.array(flux, copy=0,dtype = float)
#   npix = (size(wave))[1]
   npix = len(wave)
#   print npix
#   nrbin = (size(nwave))[1]
   nrbin = len(nwave)
#   rbin = dblarr(nrbin)
   rbin = []
#   arc1 = dblarr(2)
   arc1= []
#   arcone = dblarr(2)
   arcone = [None]*2
#   arctwo = dblarr(2)
   arctwo = [None]*2
#   darcone = dblarr(2)
   darcone = [None]*2
   #print, wave[npix-1]
   #print, double(wave[0])
   #print, double(wave[npix-1])
   #print, double(wave[npix-2])
   arc1.append(np.array(wave[0], copy=0).astype(float) - (np.array(wave[npix - 1], copy=0).astype(float) - np.array(wave[0], copy=0).astype(float)) / (npix - 1))
   arc1.append(((np.array(wave[npix - 1], copy=0).astype(float) - np.array(wave[0], copy=0).astype(float)) * npix) / (npix - 1))
#   arc2 = dblarr(2)
   arc2 = []
   arc2.append(np.array(nwave[0], copy=0).astype(float) - (np.array(nwave[nrbin - 1], copy=0).astype(float) - np.array(nwave[0], copy=0).astype(float)) / (nrbin - 1))
   arc2.append(((np.array(nwave[nrbin - 1], copy=0).astype(float) - np.array(nwave[0], copy=0).astype(float)) * nrbin) / (nrbin - 1))
   
   # Declaration of parameters.
   maxcoeff = 11
   accbin = 0.0001e0
   narc1 = 2
   narc2 = 2
   # Declaration of variables.
   #     REAL*8 RX2,X1,X2,RSLOPE,DX,A,B,C,D,
   #    X DD,DDD,Y
   
   #print, 'arc1 ', arc1
   #print, 'arc2 ', arc2
   #print, 'nrbin ', nrbin
   inarc1 = abs(narc1)
   inarc2 = abs(narc2)
   in11 = inarc1 - 1
   in12 = inarc1 - 2
   in21 = inarc2 - 1
   lstart = 0
   acc = accbin / npix * 1.0e0
   print 'npix', npix
   rnpix = 1.0e0 / npix
   #print, rnpix
   rnrbin = 1.0e0 / nrbin
   
   rx2 = 0.5e0 * rnrbin
#   print 'inarc1', inarc1
   for i in np.arange(0, (inarc1 - 1)+(1)):
      arcone[inarc1 - 1 - i] = arc1[i]
      darcone[inarc1 - 1 - i] = (i) * arc1[i]
   for i in np.arange(0, (inarc2 - 1)+(1)):
      arctwo[inarc2 - 1 - i] = arc2[i]
#   print 'darcone ', darcone, in12
#   print 'arctwo',arctwo   
   rslope = 1 / womply(darcone, in12, 0.2e0)
#   print 'rslope ', rslope
   x1 = 0.2
   rx2, x1, arctwo, darcone, in21, acc, rslope, in11, in12, arcone = womrebin_tform(rx2, x1, arctwo, darcone, in21, acc, rslope, in11, in12, arcone)
   print 'x1', x1
   x1 = x1 * npix
#   print x1
   dx = 0
   nsgn = 1
   if ((womply(arctwo, in21, 1.0e0) - arc2[1]) * rslope < 0.0e0):   
      nsgn = -1
   nstop = nrbin
   
   j1 = round(x1) - 1
   print 'j1' , j1
   for k in np.arange(0, (nstop - 1)+(1)):
   
      rx2 = rx2 + rnrbin
      x2 = (x1 + dx) * rnpix
      rx2, x2, arctwo, darcone, in21, acc, rslope, in11, in12, arcone = womrebin_tform(rx2, x2, arctwo, darcone, in21, acc, rslope, in11, in12, arcone)
      x2 = x2 * npix
      #print, 'x2 ', x2
#      print 'npix' , npix
      dx = x2 - x1
      #       print, 'dx ', dx
      j2 = round(x2) - 1
      d = 0
      
      if (lstart == 0):   
         lstart = 1
         m1 = max(np.concatenate([np.array(np.concatenate([j1 - 1, npix - 1]), copy=0).min(), 1]))
         m2 = max(np.concatenate([np.array(np.concatenate([j1, npix - 1]), copy=0).min(), 1]))
         m3 = max(np.concatenate([np.array(np.concatenate([j1 + 1, npix - 1]), copy=0).min(), 1]))
         #print, 'aj1,j2 ', j1, j2
         
         #                print, 'm', m1, m2, m3
         #              print, 'w ', WDATA(M1), WDATA(M3), WDATA(M2)
         a = (wdata(m1) + wdata(m3)) * 0.5e0
         b = (a - wdata(m1)) * 0.5e0
         c = (13.0e0 / 12.0e0) * wdata(m2) - a / 12.0e0
         a = (a - wdata(m2)) / 3.0e0
         y = x1 - j1 - 1
         #print, 'y,x1,j1,', y, x1, j1
         dd = nsgn * ((((a * y) + b) * y + c) * y - b * 0.25e0) + a * 0.125e0 + c * 0.5e0
         #print, a, b, c, y, dd
      #print, 'nsgn', nsgn
      #print, 'j1,j2 ', j1, j2
      m1 = max(np.concatenate([np.array(np.concatenate([j2 - 1, npix - 1]), copy=0).min(), 1]))
      m2 = max(np.concatenate([np.array(np.concatenate([j2, npix - 1]), copy=0).min(), 1]))
      m3 = max(np.concatenate([np.array(np.concatenate([j2 + 1, npix - 1]), copy=0).min(), 1]))
      #print, 'm1,m2.m3', m1, m2, m3
      a = (wdata(m1) + wdata(m3)) * 0.5e0
      b = (a - wdata(m1)) * 0.5e0
      c = 1.083333333333333e0 * wdata(m2) - a * 0.08333333333333333e0
      a = (a - wdata(m2)) * 0.3333333333333333e0
      #print, 'x2,j2,', x2, j2
      
      y = x2 - j2 - 1
      #print, 'y', y
      d = d - dd
      dd = nsgn * ((((a * y) + b) * y + c) * y - b * 0.25e0)
      ddd = a * 0.125e0 + c * 0.5e0
      d = d + dd - ddd
      dd = dd + ddd
      
      for kk in np.arange(j1, (j2)+(nsgn), nsgn):
         d = d + wdata(max(np.concatenate([np.array(np.concatenate([kk, npix - 1]), copy=0).min(), 1])))
      #print, 'd', d
      rbin[k] = d / abs(dx)
      # print, rbin(k)
      x1 = x2
      j1 = j2
   nflux = np.array(rbin, copy=0).astype(float)
   # stop
   return _ret()
   

