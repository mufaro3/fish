;; We'll need NumCL to simplify the vector calculations.
;; This is equivalent to `import numpy as np`
(ql:quickload :numcl)
(ql:quickload :alexandria)
(defpackage :fish
  (:use :cl)
  (:local-nicknames (:nc :numcl)))
(in-package :fish)

;; First, we build the fish model based on the specificaitons.
;; Defining the constants for the simulation

;; For this simluation, we assume that we are modelling Pink
;; Salmon in the Fraser River of British Columbia.

;; Average length of Pink Salmon is about 55 cm.
(defparameter *fish-length* 0.55) ; meters

;; The Fraser River has an annual average volumetric flow rate
;; at its mouth of about 3550 m^3/s.
(defparameter *volumetric-flow-rate* 3550) ; meter^3 / second

;; The number of fish to model in this school
(defparameter *num-fish* 20)

;; Fish Model

(defstruct fish
  (position nil :type nc:array)
  (orientation nil :type nc:array))
(defparameter *cs-area* (* 4 pi (expt *fish-length* 2)))
(defparameter *self-propelled-speed* (/ *volumetric-flow-rate* *cs-area*))

;; Simulation Calculation Functions

(defun half-length-vector (fish)
  "The distance vector from the center-of-mass to the front or back"
  (nc:* 1/2 *fish-length* (fish-orientation fish)))

(defun source-pos (fish)
  "The front of the fish, x_f,i"
  (nc:+ (fish-position fish)
	(half-length-vector fish)))

(defun sink-pos (fish)
  "The back of the fish, x_b,i"
  (nc:- (fish-position)
	(half-length-vector fish)))

(defun norm (v)
  "Distance vector calculation"
  (nc:sqrt (nc:dot v v)))

(defun compute-decrease-vector (a b)
  "Same thing as computing the vector representation of 1/r^2"
  ((let ((distance-vec (nc:- a b)))
     (nc:/ distance-vec (expt (norm distance-vec) 3)))))

(defun source-decrease-vector (fish-a fish-b)
  "Computes the decrease vector 1/r^2 between the sources and sinks"
  (nc:-
   (compute-decrease-vector (source-pos fish-a) (source-pos fish-b))
   (compute-decrease-vector (source-pos fish-a) (sink-pos   fish-b))))

(defun sink-decrease-vector (fish-a fish-b)
  "Computes the decrease vector 1/r^2 for computing the sink velocity"
  (nc:-
   (compute-decrease-vector (sink-pos fish-a) (source-pos fish-b))
   (compute-decrease-vector (sink-pos fish-a) (sink-pos   fish-b))))

(defun compute-velocity (fish other-fishes contibution-fn)
  "Computes the total contribution"
  (let ((compute-contibution (alexandria:curry #'contribution-fn fish))
	(total-distance-contributions
	  (nc:sum (nc:asarray (mapcar #'compute-contribution other-fishes)))))
    (nc:* *self-propelled-speed*
	  (nc:+ (fish-orientation fish)
		(nc:* *fish-length* total-distance-contributions)))))

(defun source-vel (fish other-fishes)
  "Computes the source velocity for a given fish"
  (compute-velocity fish other-fishes source-decrease-vector))

(defun sink-vel (fish other-fishes)
  "Computes the sink velocity for a given fish"
  (compute-velocity fish other-fishes sink-decrease-vector))

(defun translational-vel (fish other-fishes)
  "Computes the translational velocity for a given fish"
  (nc:/ (nc:+ (source-vel fish other-fishes) (sink-vel fish other-fishes)) 2))

(defun rotational-vel (fish other-fishes)
  "Computes the rotational velocity for a given fish"
  (let ((relative-vel (nc:- (source-vel fish other-fishes)
			    (sink-vel fish other-fishes)))
	(lagrange-mult
	  (nc:* -1/2 (nc:dot relative-vel (fish-orientation fish)))))
    (nc:/
     (nc:+ relative-vel
	   (nc:* 2 lagrange-mult (fish-orientation fish)))
     *fish-length*)))

;; Simulation Function

(defun simulate (initial-state stop-time time-step) nil)
