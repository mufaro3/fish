;; We'll need NumCL to simplify the vector calculations.
;; This is equivalent to `import numpy as np`
(ql:quickload :numcl)
(ql:quickload :alexandria)
(ql:quickload :cl-opengl) ;; for plotting, equivalent to matplotlib
(ql:quickload :cl-glut)
(defpackage :fish
  (:use :cl)
  (:local-nicknames (:nc :numcl))
  (:export :simulate :orientation-vector))
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
(defparameter *initial-generation-bounds* (nc:asarray '(10 10 10)))

;; Fish Model

(defstruct fish
  (position nil :type nc:array)
  (orientation nil :type nc:array))
(defparameter *cs-area* (* 4 pi (expt *fish-length* 2)))
(defparameter *self-propelled-speed* (/ *volumetric-flow-rate* *cs-area*))

;; Linear Algebra Functions

(defun dot (a b)
  "Dot Product"
  (declare (type nc:array a b))
  (nc:sum (nc:* a b)))

(defun norm (v)
  "Distance vector calculation"
  (declare (type nc:array v))
  (nc:sqrt (dot v v)))

;; Simulation Calculation Functions

(defun half-length-vector (fish)
  "The distance vector from the center-of-mass to the front or back"
  (declare (type fish fish))
  (nc:* 1/2 *fish-length* (fish-orientation fish)))

(defun source-pos (fish)
  "The front of the fish, x_f,i"
  (declare (type fish fish))
  (nc:+ (fish-position fish)
	(half-length-vector fish)))

(defun sink-pos (fish)
  "The back of the fish, x_b,i"
  (declare (type fish fish))
  (nc:- (fish-position fish)
	(half-length-vector fish)))

(defun compute-decrease-vector (a b)
  "Same thing as computing the vector representation of 1/r^2"
  (declare (type nc:array a b))
  (let ((distance-vec (nc:- a b)))
     (nc:/ distance-vec (expt (norm distance-vec) 3))))

(defun source-decrease-vector (fish-a fish-b)
  "Computes the decrease vector 1/r^2 between the sources and sinks"
  (declare (type fish fish-a fish-b))
  (nc:-
   (compute-decrease-vector (source-pos fish-a) (source-pos fish-b))
   (compute-decrease-vector (source-pos fish-a) (sink-pos   fish-b))))

(defun sink-decrease-vector (fish-a fish-b)
  "Computes the decrease vector 1/r^2 for computing the sink velocity"
  (declare (type fish fish-a fish-b))
  (nc:-
   (compute-decrease-vector (sink-pos fish-a) (source-pos fish-b))
   (compute-decrease-vector (sink-pos fish-a) (sink-pos   fish-b))))

(defun compute-velocity (fish other-fishes contribution-fn)
  "Computes the total contribution"
  (declare (type fish     fish)
	   (type vector   other-fishes)
	   (type function contribution-fn))
  (let* ((compute-contribution
	   (alexandria:curry contribution-fn fish))
	 (total-distance-contributions
	   (nc:sum (nc:asarray
		    (map 'vector compute-contribution other-fishes)))))
    (nc:* *self-propelled-speed*
	  (nc:+ (fish-orientation fish)
		(nc:* *fish-length* total-distance-contributions)))))

(defun source-vel (fish other-fishes)
  "Computes the source velocity for a given fish"
  (declare (type fish fish)
	   (type vector other-fishes))
  (compute-velocity fish other-fishes #'source-decrease-vector))

(defun sink-vel (fish other-fishes)
  "Computes the sink velocity for a given fish"
  (declare (type fish fish)
	   (type vector other-fishes))
  (compute-velocity fish other-fishes #'sink-decrease-vector))

(defun translational-vel (fish other-fishes)
  "Computes the translational velocity for a given fish"
  (declare (type fish fish)
	   (type vector other-fishes))
  (nc:/ (nc:+ (source-vel fish other-fishes)
	      (sink-vel fish other-fishes))
	2))

(defun rotational-vel (fish other-fishes)
  "Computes the rotational velocity for a given fish"
  (declare (type fish fish)
	   (type vector other-fishes))
  (let* ((relative-vel
	   (nc:- (source-vel fish other-fishes)
		 (sink-vel fish other-fishes)))
	(lagrange-mult
	  (nc:* -1/2 (dot relative-vel (fish-orientation fish)))))
    (nc:/
     (nc:+ relative-vel
	   (nc:* 2 lagrange-mult (fish-orientation fish)))
     *fish-length*)))

;; Simulation Functions

;; Theta [0, 2PI) and Phi [0,PI)
(defparameter *angle-bounds* (vector (* 2 pi) pi))
(defun generate-random-fish ()
  "Computes a random fish with a random heading direction"
  (make-fish :position (map 'vector #'random *initial-generation-bounds*)
	     :orientation (map 'vector #'random *angle-bounds*)))

(defun generate-random-distribution (n)
  "Generates n random autofish within a maximum starting bounds"
  (declare (type fixnum n))
  (make-array n :initial-contents
	      (loop :for i :below n
		    :collect (generate-random-fish))))

(defparameter *orientation-vector-length* 2)
(defun spherical->cartesian (theta phi)
  "Converts the spherical orientation to a cartesian unit vector"
  (nc:asarray (vector
	       (* (sin phi) (cos theta))
	       (* (sin phi) (sin theta))
	       (cos phi))))

(defun orientation-vector (fish)
  "Generates the orientation vector for a given fish"
  (declare (type fish fish))
  (nc:* *orientation-vector-length*
	(apply #'spherical->cartesian (fish-orientation fish))))

(defun display-state (state)
  "Displays the state in a three-dimensional plot"
  (declare (type vector state))
  (gl:clear :color-buffer-git :depth-buffer-bit)

  (loop :for fish :across state :do
    ;; Display the fish as a point
    (gl:begin :points)
    (apply #'gl:vertex (fish-position fish))
    (gl:end)

    ;; Display the orientation as a small line
    (gl:begin :lines)
    (apply #'gl:vertex (fish-position fish))
    (apply #'gl:vertex (nc:+ (fish-position fish)
			     (orientation-vector fish)))
    (gl:end))
  
  (gl:flush)
  (swap-buffers)
  (sleep 1/60)) ; 60 FPS

(defstruct octree-node xyz xy-z x-yz x-y-z -x-y-z -x-yz -xy-z -xyz)
(defun generate-barnes-hut-octree (all-fish)
  (labels (())
    (recurse-octree-gen )))

(defun barnes-hut-simplify (other-fish-tree fish)
  "Applies the Barnes-Hut Approximation to cluster farther fish"
  ())

(defun calculate-next-state (state time-step use-barnes-hut?)
  "Computes the next state based on the previous state and the time-step"
  (declare (type vector state)
	   (type number time-step))
  (let ((barnes-hut-tree
	  (if use-barnes-hut?
	      (generate-barnes-hut-octree (coerce state 'list)))))
    (coerce
     (loop :for i :below (length state)
	   :for fish := (aref state i)
	   :for other-fish := (concatenate 'vector
					   (subseq state 0 i)
					   (subseq state (1+ i)))
	   :collect
	   (progn
	     (if use-barnes-hut?
		 (setf other-fish
		       (barnes-hut-simplify other-fish fish)))
	     (let* ((trans-vel (translational-vel fish other-fishes))
		    (rot-vel (rotational-vel fish other-fishes))
		    (pos-delta (nc:* trans-vel time-step))
		    (orientation-delta (nc:* rot-vel time-step)))
	       (make-fish
		(nc:+ (fish-position fish) pos-delta)
		(nc:+ (fish-orientation fish) orientation-delta)))))
     'vector)))

(defun simulate (stop-time time-step use-barnes-hut?)
  "The main simulation function"
  (declare (type number stop-time time-step)
	   (type bool   use-barnes-hut?))
  (let ((time 0)
	(state (generate-random-distribution *num-fish*)))
    (loop :while (< time stop-time)
	  :do (display-state state)
	      (setf state (calculate-next-state state time-step
						use-barnes-hut?))
	      (incf time time-step))))
