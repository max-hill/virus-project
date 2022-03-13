;; Auxillary file to generate-partition-file.sh

(defun count-digits (n)
  "return the number of digits of n"
  (1+ (floor (/ (log n) (log 10)))))

(defparameter *A* (parse-integer (second *posix-argv*)))
(defparameter *B* (parse-integer (third *posix-argv*)))
;; (defparameter *position-padding* (count-digits *A*))
;; (defparameter *gene-number-padding* (count-digits (ceiling (/ *A* *B*))))
(defparameter *position-padding* 0)
(defparameter *gene-number-padding* 0)

(defun pad-zeros (n k)
  "left-pad n with k zeros"
  (format nil (format nil "~a~a,'0d"  "~" k) n))

(defun generate-partition-file (alignment-length block-length &optional (block-starting-position 1) (gene-number 1))
  "INPUT: two positive integers: (1) total alignment length and (2) desired block length. OUTPUT: the contents of an msa
partition file. USAGE: (make-partition-file 14125 500). The optional variables are internal (dummy) variables and are
not to be specified by the user. 2022-03-07 max-hill"
    (if (> (+ block-starting-position block-length) alignment-length)
        (format t "DNA, gene_~a=~a-~a~%"
		(pad-zeros gene-number *gene-number-padding*)
		(pad-zeros block-starting-position *position-padding*)
		alignment-length)
      (progn
	(format t "DNA, gene_~a=~a-~a~%"
		(pad-zeros gene-number *gene-number-padding*)
		(pad-zeros block-starting-position *position-padding*)
		(pad-zeros (+ block-starting-position block-length -1) *position-padding*))
	(generate-partition-file alignment-length block-length
				 (+ block-starting-position block-length)
				 (1+ gene-number)))))

