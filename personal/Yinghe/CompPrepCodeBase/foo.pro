FUNCTION WOMPLY, a, i, x
;
;             I+1
;  EVALUATES  SUM A(J) * X**(I+1-J)
;             J=1
;
;  I IS THUS THE ORDER OF THE POLY, WHICH HAS I+1 COEFFICIENTS,
;  AND A(I+1) IS THE CONSTANT TERM.
;

 
      IF(I LT 0) THEN begin
      WOMPLY=0.D0
      print, '*** VALUE OF I OUT OF BOUNDS IN WOMPLY ***'
      RETURN, womply
      ENDIF
 
      WOMPLY=A(0)
        if (i eq 0) then return, womply
      for j = 2, i+1 do begin
       WOMPLY=X*WOMPLY+A(J-1)
        endfor
      RETURN, womply
  
      END
pro   WOMREBIN_TFORM, rx, x, arctwo, darcone, in21, acc, rslope, in11, in12, arcone
;
; COPYRIGHT 1982 M. ASHLEY. See comments with SUBROUTINE REBIN.
;
; Version 1.2     23 Apr 1982    - More arguments put into COMMON.
; Version 1.1     29 Mar 1982    - First working version
;

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

pro womashrebin, wave, flux, nwave, nflux

;  COMMENTS BELOW FROM ASHREBIN OF SARAH AND LOLITA
;  I have translated this into IDL
;  TM  2000-05-12
;***************************************************************
;    Please do not make changes to this program without documenting them below.
;
;
; Version 1.4      8 Jul 1982 - Minor bug fixed (NSGN had the wrong value if
;                               you were rebinning from decreasing to increasing
;                               wavelength scales (or vice versa) and if, at the
;                               same time the slopes were less than a certain
;                               small number).
; Version 1.3     23 Apr 1982 - Documentation added; some code changes.
; Version 1.2     22 Apr 1982 - MACRO PLY replaced by FORTRAN version.
; Version 1.1     29 Mar 1982 - First "working version".
;
;
; COPYRIGHT 1982 M. ASHLEY.  Use of this program at this stage of
; development is entirely the user's responsibility. Unauthorized
; use or copying of this program is strictly prohibited.
;
;
; HISTORY OF DEVELOPMENT
; ----------------------
;
;      This subroutine was written by Michael Ashley as a fast
; replacement for Peter Young's subroutine of the same name, which had
; been the bottleneck when running the data reduction package LOLITA.
;
;
; QUICK SUMMARY OF OPERATION
; --------------------------
;
;      REBIN will rebin data from the first NPIX entries of array WDATA
; into the array RBIN. If IQUAD equals zero the rebinning will be done
; using linear interpolation, otherwise quadratic interpolation will be
; used. The program can be used in two modes: if IMODE equals zero then
; the output bins will be shifted by SSKEW relative to the input bins,
; and each output bin will consis of the sum of NRBIN input bins; otherwise
; the rebinnning will be done from the wavelength scale ARC1 to the scale
; ARC2, with NRBIN now meaning the number of output bins. In this later
; mode the absolute magnitudes of NARC1 and NARC2 are the number of
; coefficients in ARC1 and ARC2 respectively, and their signs indicate
; whether logarithmic scales are to be used.
;
;
; DESCRIPTION OF ARGUMENTS
; ------------------------
;
;   IMODE   - INTEGER -  If equal to zero, selects mode 0, in which input
;                        bins can be grouped and shifted to form the output
;                        bins.
;
;                        If not equal to zero, selects mode 1, in which the
;                        input bins are tranferred to the output bins in
;                        accordance with the wavelength scales ARC1 and ARC2.
;
;   IQUAD   - INTEGER -  If equal to zero, selects linear interpolation when
;                        rebinning.
;
;                        If not equal to zero, selects quadratic interpolation.
;
;   WDATA   - REAL*4  -  Array containing input data.
;
;   NPIX    - INTEGER -  The number of bins (elements) in WDATA.
;
;   RBIN    - REAL*4  -  Output Array, into which the rebinned data will be
;                        placed.
;
;   NRBIN   - INTEGER -  MODE 0: The number of input bins to be added together
;                                to form one output bin.
;
;                        MODE 1: The number of bins in RBIN.
;
;   SSKEW   - REAL*4  -  The number of bins the input array is to shifted (used
;                        in mode 0 only).
;
;   NARC1   - INTEGER -  Used in MODE 1 only. The absolute value of NARC1 is
;                        the number of coefficients in ARC1. If NARC1 is
;                        negative, then a logarithmic scale is assumed.
;
;   ARC1    - REAL*8  -  Used in MODE 1 only. The array of polynomial coeff-
;                        icients defining the wavelength scale of the input
;                        data.
;
;   NARC2   - INTEGER -  As for NARC1 except refers to the output wavelength
;                        scale ARC2.
;
;   ARC2    - REAL*8  -  As for ARC2 except it is the output wavelength scale.
;
;
;
; POTENTIAL PROBLEMS
; ------------------
;
;      REBIN assumes a certain degree of responsibility on the part of the
; program that calls it. Only minimal checking is done to ensure sensible
; arguments. Needless to say, if REBIN crashes, a large amount of time can
; be wasted. It is up to the user to ensure that inappropriate arguments
; can never be sent to REBIN.
;
;      REBIN may crash in the following ways -
;
; [1] Floating point overflow in REBIN_POLY, due to invalid arc coefficients.
;
; [2] Floating point overflow, or argument out of range in WOMREBIN_TFORM when
;     taking the logarithm and/or exponent of a result returned from
;     REBIN_POLY. This would be caused by having invalid arc coefficients.
;
; [3] In the Newton Raphson solution for output bin edges in terms of the
;     input bins it is assumed that the wavelength scales are fairly smooth,
;     and single valued over the interval of interest. If bizarre coefficients
;     are used, all sorts of difficulties may occur.
;
; [4] Divide by zero in REBIN and/or WOMREBIN_TFORM if the derivative of the
;     ARC1 polynomial equals zero at a tested point.
;
; [5] The user must ensure tha the arrays WDATA and RBIN are at least as
;     large as NPIX and NRBIN respectively (except in the case where IMODE
;     equals zero, when RBIN must be as large as necessary to fit the
;     rebinned data).
;
; [6] ARC1 and ARC2 must be dimensioned at least as large as the absolute
;     values of NARC1 and NARC2 respectively.
;
;
;
; LINEAR/QUADRATIC REBINNING
; --------------------------
;
;      IQUAD is used as a flag to set either linear or quadratic rebinning.
; This refers to the way in which the input data is interpolated to arrive
; at the number of counts in a fraction of a bin. Suppose that one of the
; output bins completely covers one of the input bins, and partially covers
; the bins on either side. Bin <n> is treated as though it collects all the
; counts from <n-0.5> to <n+0.5>. It is not treated as sampling the data
; a position <n> only. In linear rebinning it is assumed that the count
; rate was uniform across each bin. Thus, for example, if an output bin covers
; one quarter of an input bin, then it receives one quarter of the counts
; in that bin.
;
;      In quadratic rebinning it is assumed that, locally, the original data
; followed a parabola, and that the number of counts in bin <n> results
; from integrating this parabola from <n-0.5> to <n+0.5>. When a fraction
; of an input bin is required, the equation of the parabola which fits this
; bin and the neighbouring bins is found, and then integrated appropriately.
;
;      In both linear and quadratic interpolation, after the number of counts
; in each output bin has been found, these are divided by the width of
; each output bin in terms of input bins. This ensures that the mean level
; of the data doess not change. Hence, after rebinning it is not necessarily
; valid to say that the noise in a bin is the square root of the number of
; counts in that bin.
;
;
;
;print, size(nwave)
wdata = double(flux)
npix = (size(wave))[1]
;print, npix
nrbin = (size(nwave))[1]
rbin = dblarr(nrbin)
arc1 = dblarr(2)
arcone = dblarr(2)
arctwo = dblarr(2)
darcone = dblarr(2)
;print, wave[npix-1]
;print, double(wave[0])
;print, double(wave[npix-1])
;print, double(wave[npix-2])
arc1(0) = double(wave[0])-(double(wave[npix-1])-double(wave[0]))/(npix-1)
arc1(1) = ((double(wave[npix-1])-double(wave[0]))*npix)/(npix-1)
arc2 = dblarr(2)
arc2(0) = double(nwave[0])-(double(nwave[nrbin-1])-double(nwave[0]))/(nrbin-1)
arc2(1) = ((double(nwave[nrbin-1])-double(nwave[0]))*nrbin)/(nrbin-1)
; Declaration of parameters.
MAXCOEFF=11
ACCBIN=0.0001D0
narc1 = 2
narc2 = 2
; Declaration of variables.
;     REAL*8 RX2,X1,X2,RSLOPE,DX,A,B,C,D,
;    X DD,DDD,Y

;print, 'arc1 ', arc1
;print, 'arc2 ', arc2
;print, 'nrbin ', nrbin
      INARC1=ABS(NARC1)
      INARC2=ABS(NARC2)
      IN11=INARC1-1
      IN12=INARC1-2
      IN21=INARC2-1
lstart = 0
      ACC=ACCBIN/NPIX*1.0D0
;print, npix
      RNPIX=1.0D0/NPIX
;print, npix
      RNRBIN=1.0D0/NRBIN
 
            RX2=0.5D0*RNRBIN
;print, 'inarc1', inarc1
        for i = 0, inarc1-1 do begin
               ARCONE(INARC1-1-I)=ARC1(I)
               DARCONE(INARC1-1-I)=(I)*ARC1(I)
            ENDfor
         for i = 0, inarc2-1 do begin  
               ARCTWO(INARC2-1-I)=ARC2(I)
         endfor
;print, 'darcone ', darcone, in12
            RSLOPE=1/WOMPLY(DARCONE,IN12,0.2D0)
        ;   print, 'rslope ', rslope
            X1=0.2
            WOMREBIN_TFORM, RX2,X1, arctwo, darcone, in21, acc, rslope, $
         in11, in12, arcone
;        print, x1
            X1=X1*NPIX
;print, x1
            DX=0
            NSGN=1
            IF((WOMPLY(ARCTWO,IN21,1.0D0)-ARC2(1))*RSLOPE LT 0.0D0) then  NSGN=-1
            NSTOP=NRBIN

      J1=round(X1)-1L
;print, j1
      for k = 0L, nstop-1 do begin
        
               RX2=RX2+RNRBIN
               X2=(X1+DX)*RNPIX
               WOMREBIN_TFORM, RX2,X2, arctwo, darcone, in21, acc, $
      rslope, in11, in12, arcone
               X2=X2*NPIX
;print, 'x2 ', x2
;print, npix
               DX=X2-X1
;       print, 'dx ', dx
         J2=round(X2)-1L
         D=0
         
               IF(LSTART eq 0) THEN begin
                     LSTART=1
                     M1=MAX([MIN([J1-1,NPIX-1]),1])
                     M2=MAX([MIN([J1,NPIX-1]),1])
                     M3=MAX([MIN([J1+1,NPIX-1]),1])
;print, 'aj1,j2 ', j1, j2

;                print, 'm', m1, m2, m3
;              print, 'w ', WDATA(M1), WDATA(M3), WDATA(M2)
                     A=(WDATA(M1)+WDATA(M3))*0.5D0
                     B=(A-WDATA(M1))*0.5D0
                     C=(13.0D0/12.0D0)*WDATA(M2)-A/12.0D0
                     A=(A-WDATA(M2))/3.0D0
                     Y=X1-J1-1
;print, 'y,x1,j1,', y, x1, j1
                     DD=NSGN*((((A*Y)+B)*Y+C)*Y-B*0.25D0)+ A*0.125D0+C*0.5D0
;print, a, b, c, y, dd    
ENDIF
;print, 'nsgn', nsgn
;print, 'j1,j2 ', j1, j2
               M1=MAX([MIN([J2-1,NPIX-1]),1])
               M2=MAX([MIN([J2,NPIX-1]),1])
               M3=MAX([MIN([J2+1,NPIX-1]),1])
;print, 'm1,m2.m3', m1, m2, m3
               A=(WDATA(M1)+WDATA(M3))*0.5D0
               B=(A-WDATA(M1))*0.5D0
               C=1.083333333333333D0*WDATA(M2)-A*0.08333333333333333D0
               A=(A-WDATA(M2))*0.3333333333333333D0
;print, 'x2,j2,', x2, j2

               Y=X2-J2-1
;print, 'y', y
               D=D-DD
               DD=NSGN*((((A*Y)+B)*Y+C)*Y-B*0.25D0)
               DDD=A*0.125D0+C*0.5D0
               D=D+DD-DDD
               DD=DD+DDD
         
         for kk = j1, j2, nsgn do begin
            D=D+WDATA(MAX([MIN([KK,NPIX-1]),1]))
         ENDFOR
;print, 'd', d
         RBIN(K)=D/ABS(DX)
        ; print, rbin(k)
         X1=X2
         J1=J2
      ENDfor
      nflux = float(rbin)
; stop
      RETURN

      END