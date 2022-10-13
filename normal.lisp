#|
Time-stamp: <2020-12-10 10:28:11 jpp>
File:       normal.lisp
Author:     JP Prin
License:    MIT
Purpose:    normal distribution

This file contains an implementation of the normal inverse cdf function, which
can be used in a lot of other applications, particularly ones that need to
calculate implied volatility.

The method I use below was published by Peter Acklam to approximate the quartile
function (aka the inverse cumulative distribution function) of the normal
distribution.  He made the algorithm freely available, but unfortunately the 
link to his website no longer works.

Peter Acklam estimates that "the absolute value of the relative error is less 
than 1.15 x 10 ^ -9 in the entire region.  Relative error is defined as

  x_approx - x_exact / x_exact

He says that the error could be greater if x < -38, but since these values 
almost never occur in practice he disregards that as being a problem since it 
would require an input of 2.885428351 x 10^-316.

He further comments that using IEEE double precision arithmetic we can not even
represent numbers that small in full precision, so the error range is fine.
|#

(in-package :godel.statistics)

;;; Coefficients in rational approximations

(defparameter a-params
  '(-3.969683028665376e+01   ; a1
     2.209460984245205e+02   ; a2
    -2.759285104469687e+02   ; a3
     1.383577518672690e+02   ; a4
    -3.066479806614716e+01   ; a5
    2.506628277459239e+00))  ; a6


(defparameter b-params
  '(-5.447609879822406e+01    ; b1
     1.615858368580409e+02    ; b2
    -1.556989798598866e+02    ; b3
     6.680131188771972e+01    ; b4
    -1.328068155288572e+01))  ; b5

(defparameter c-params
  '(-7.784894002430293e-03    ; c1
    -3.223964580411365e-01    ; c2
    -2.400758277161838e+00    ; c3
    -2.549732539343734e+00    ; c4
     4.374664141464968e+00    ; c5
     2.938163982698783e+00))  ; c6

(defparameter d-params
  '(7.784695709041462e-03    ; d1
    3.224671290700398e-01    ; d2
    2.445134137142996e+00    ; d3
    3.754408661907416e+00))  ; d4

;;; Break points

(defparameter p-low   0.02425)
(defparameter p-high  (1- p-low))

;;; the approximation

(defun normal-expand (val lst)
  (reduce #'(lambda (x y) (+ (* x val) y)) lst))

(defun inverse-normal-cdf (p)

  (assert (and (plusp p) (< p 1))

  (cond
    ((< p p-low)
     (let ((q (sqrt (* -2 (log p)))))
       (/ (normal-expand q c-parems) (1+ (normal-expand q d-parems)))))

    ((<= p p-high)
     (let* ((q (- p 0.5))
	    (r (* q q)))
       (/ (* q (normal-expand r a-parems))
	  (1+ (* r (normal-expand r b-parems))))))

    (t (let ((q (sqrt (* -2 (log (- 1 p))))))
	 (/ (normal-expand q c-parems) (1+ (* q (normal-expand q d-parems)))))))))

