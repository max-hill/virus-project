(defun generate-partition-file (alignment-length block-length &optional (block-starting-position 1) (gene-number 1))
  "INPUT: two positive integers: (1) total alignment length and (2) desired
block length. OUTPUT: the contents of an msa partition file. USAGE: (make-partition-file 14125
500). The optional variables are internal (dummy) variables and are not to be
specified by the user. 2022-03-07 max-hill"
  (if (> (+ block-starting-position block-length) alignment-length)
      (format t "~%gene-~A ~A-~A" gene-number block-starting-position alignment-length)
      (progn (format t "~%gene-~A ~A-~A" gene-number block-starting-position (+ block-starting-position block-length -1))
	     (generate-partition-file alignment-length block-length (+ block-starting-position block-length) (1+ gene-number)))))
      
