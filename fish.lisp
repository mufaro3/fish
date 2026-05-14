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
(defparameter *initial-generation-box-length* 10)
(defparameter *barnes-hut-threshold* 1/4) ; theta

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
  (declare (type fish fish)
	   (type (vector fish) other-fishes)
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
	   (type (vector fish) other-fishes))
  (compute-velocity fish other-fishes #'source-decrease-vector))

(defun sink-vel (fish other-fishes)
  "Computes the sink velocity for a given fish"
  (declare (type fish fish)
	   (type (vector fish) other-fishes))
  (compute-velocity fish other-fishes #'sink-decrease-vector))

(defun translational-vel (fish other-fishes)
  "Computes the translational velocity for a given fish"
  (declare (type fish fish)
	   (type (vector fish) other-fishes))
  (nc:/ (nc:+ (source-vel fish other-fishes)
	      (sink-vel fish other-fishes))
	2))

(defun rotational-vel (fish other-fishes)
  "Computes the rotational velocity for a given fish"
  (declare (type fish fish)
	   (type (vector fish) other-fishes))
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
(defparameter *initial-generation-bounds*
  (nc:asarray (make-array 3 :initial-element *initial-generation-box-length*)))

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
  (declare (type number theta phi))
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
  (declare (type (vector fish) state))
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

;; State Update Functions

(defstruct octree-node
  (center (nc:/ *initial-generation-bounds* 2) :type nc:array)
  (half-size (/ *initial-generation-box-length* 2) :type number)
  (children nil :type (vector octree-node 8))
  (object nil :type fish))

(defun octree-get-child-index (node next)
  "Obtains the `child index` of the fish, or the octant to place it in"
  (declare (type octree-node node)
	   (type fish next))
  (let* ((center (octree-node-center node))
	 (pos (fish-position next))
	 (x-bit (if > (aref pos 0) (aref center 0)) 1 0)
	 (y-bit (if > (aref pos 1) (aref center 1)) 2 0)
	 (z-bit (if > (aref pos 2) (aref center 2)) 4 0))
    (+ x-bit y-bit z-bit)))

(defun octree-insert (node next)
  "Inserts a fish into the octree"
  (declare (type octree-node node)
	   (type fish next))
  (let ((children (octree-node-children node))
	(index (child-index node next)))
    (if children (octree-insert (aref children index) next)
	(setf (octree-node-object node)
	      (make-octree-node :center (fish-position next)
				:half-size (/ (octree-node-half-size node) 2)
				:object next)))))

(defun generate-barnes-hut-octree (state)
  "Generates an octree representation of the state"
  (declare (optimize (safety 3) (speed 0))
	   (type (list fish) state))
  (labels ((recurse-octree-gen (tree remaining)
	     (declare (type list tree remaining))
	     (if (null remaining) tree
		 (recurse-octree-gen
		  (octree-insert tree (first remaining))
		  (rest remaining)))))
    (recurse-octree-gen (make-octree-node) state)))

(defun barnes-hut-simplify (node fish)
  "Produces a list of simplified fish interactions based on the Octree."
  (declare (type octree-node node)
	   (type fish fish))
  ;; Nodes are too far if the *fish-length* / distance < *barnes-hut-threshold* 
  (flet ((node-too-far? (node fish)
	   (declare (type octree-node node)
		    (type fish fish)
		    (optimize (speed 3) (safety 0)))
	   (let ((distance (norm (nc:- (octree-node-center node) (fish-pos fish)))))
	     (< (/ *fish-length* distance)
		*barnes-hut-threshold*))))
    (cond (((null node) nil)
	   ;; The node is too far away, so cluster the tree
	   (node-too-far? node fish)
	   (list (octree->cluster node)))
	  
	  ;; Leaf, so return the exact fish
	  ((null (octree-node-children node))
	   (octree-node-object node))

	  ;; Otherwise, return this fish and recurse across all children
	  (t (append (octree-node-object node)
		     (loop :for child :across (octree-node-children node) :append
			   (barnes-hut-simplify child fish))))))))

(defun calculate-next-state-barnes-hut (state time-step)
  "Calculates the next state using Barnes-Hut approximation"
  (declare (type vector state) (type number time-step))
  (let ((octree (generate-barnes-hut-octree (coerce 'list state))))
    (loop :for fish :across state :collect
	  ;; Perform a Barnes-Hut simplification to the octree
	  (let* ((other-fish-approximated (coerce (barnes-hut-simplify octree fish) 'vector))

		 ;; Compute translational and rotational velocity
		 (trans-vel (translational-vel fish other-fish-approximated))
		 (rot-vel (rotational-vel fish other-fish-approximated))

		 ;; Compute the change in translational position and orientation
		 (pos-delta (nc:* trans-vel time-step))
		 (orientation-delta (nc:* rot-vel time-step)))

	    ;; Produce next fish
	    (make-fish
	     (nc:+ (fish-position fish) pos-delta)
	     (nc:+ (fish-orientation fish) orientation-delta))))))

(defun calculate-next-state-iterative (state time-step)
  "Calculates the next state using standard iteration"
  (declare (type vector state) (type number time-step))

  ;; Loops through all of the autofish
  (loop :for i :below (length state)
	:for fish := (aref state i)
	:for other-fish := (concatenate 'vector
					(subseq state 0 i)
					(subseq state (1+ i)))
	:collect
	;; Compute the Translational Velocity
	(let* ((trans-vel (translational-vel fish other-fish))
	       
	       ;; Compute the rotational velocity
	       (rot-vel (rotational-vel fish other-fish))

	       ;; Compute the change in translational position
	       (pos-delta (nc:* trans-vel time-step))

	       ;; Compute the change in orientation
	       (orientation-delta (nc:* rot-vel time-step)))

	  ;; Add the changes to the fish and update 
	  (make-fish
	   (nc:+ (fish-position fish) pos-delta)
	   (nc:+ (fish-orientation fish) orientation-delta)))))

(defun calculate-next-state (state time-step use-barnes-hut?)
  "Computes the next state based on the previous state and the time-step"
  (declare (type vector state)
	   (type number time-step))
  (if use-barnes-hut?
      (calculate-next-state-barnes-hut state time-step)
      (calculate-next-state-iterative  state time-step)))

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
