;;
;; 
;; Matrix utility operations.
;;
;;
;; Copyright 2007-2017 Ivan Raikov.
;;
;;
;; This program is free software: you can redistribute it and/or
;; modify it under the terms of the GNU General Public License as
;; published by the Free Software Foundation, either version 3 of the
;; License, or (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful, but
;; WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;; General Public License for more details.
;;
;; A full copy of the GPL license can be found at
;; <http://www.gnu.org/licenses/>.
;;
;;

(module matrix-utils

        (make-numvector
         numvector-new
         list->numvector
         matrix-utils:error
         matrix-map
         matrix-fold
         matrix-fold-partial
         fill-matrix!
         matrix-ones
         matrix-zeros
         matrix-eye
         matrix-diag
         linspace
         logspace
         repmat
         meshgrid
         )

   (import scheme chicken foreign data-structures srfi-4 srfi-1)
   (require-extension interfaces srfi-4 srfi-42 srfi-4-comprehensions blas)

   
(interface numvector
  (define numvector-itemsize)
  (define (numvector-new n))
  (define (numvector-set! vect i val))
  (define (numvector-get vect i))
  (define (numvector-length vect))
  (define (list->numvector lst))
  )

#>

/* repeat a block of memory rep times */
void cmemrep(void *dest, size_t sz, size_t chunk, int rep)
{
  if(rep == 1) return;
  memcpy(dest + chunk, dest, sz*chunk); 
  if(rep & 1)
  {
    dest += chunk;
    memcpy(dest + chunk, dest, sz*chunk);
  }
  /* now repeat using a block twice as big */
  cmemrep(dest, sz, chunk<<1, rep>>1);
}

void crepmat(void *dest, const void *src, size_t sz, int ndim, 
	     int *destdimsize, int *dimsize, const int *dims, int *rep) 
{
  int i, chunk;
  int d = ndim-1;
  printf("crepmat: d= %d\n", d);
  /* copy the first repetition into dest */
  if(d == 0) 
  {
    chunk = dimsize[0];
    printf("crepmat: chunk = %d dest = %p src = %p\n", chunk, dest, src);
    memcpy(dest,src,sz*chunk);
  }
  else 
  {
    /* recursively repeat each slice of src */
       printf("crepmat: dims[%d]= %d\n", d, dims[d]);
    for(i=0; i<dims[d]; i++)
    {
       printf("crepmat: i = %d\n", i);
      crepmat(dest + i*destdimsize[d-1]*sz, src + i*dimsize[d-1]*sz, sz,
 	      ndim-1, destdimsize, dimsize, dims, rep);
    }
    chunk = destdimsize[d-1]*dims[d];
  }
  /* copy the result rep-1 times */
  cmemrep(dest,sz,chunk,rep[d]);
}

<#

(define vrepmat (foreign-safe-lambda void "crepmat"  nonnull-blob nonnull-blob unsigned-int int s32vector s32vector s32vector s32vector))

(define vmemrep (foreign-safe-lambda void "cmemrep"  nonnull-blob unsigned-int int int))


(define (matrix-utils:error x . rest)
  (let ((port (open-output-string)))
    (let loop ((objs (cons x rest)))
      (if (null? objs)
	  (begin
	    (newline port)
	    (error 'matrix-utils (get-output-string port)))
	  (begin (display (car objs) port)
		 (display " " port)
		 (loop (cdr objs)))))))



;;
;; FILL-MATRIX!:: VOBJ * ORDER * M * N * A * F * F0 * [IX * IY * EX * EY] -> A
;;
;; FILL-MATRIX! fills matrix A of size M x N with the values returned
;; by applying procedure F to each pair of indices in the matrix. VOBJ
;; is an implementation of the vector interface.
;;
;; Procedure F is of the form LAMBDA I J AX -> VAL * AX1, where I and
;; J are matrix indices, and AX is accumulator value (like in
;; fold). The initial value of AX is given by F0. Procedure F is
;; expected to return two values: the value to be placed in matrix A
;; at position [I,J], and the new accumulator value (or #f). 
;;
;; Optional arguments IX IY EX EY may specify a sub-matrix in matrix A
;; to be filled. These arguments are checked to make sure they specify
;; a valid sub-matrix.
;;
;; Procedure F must ensure that it returns values that are within the
;; range of the underlying vector type used.
;;
(define (fill-matrix! vobj order M N A f f0 #!key (ix 0) (iy 0) (ex M) (ey N))
  
  ;; optional arguments to specify a sub-matrix to be filled
  (if (not (and (fx>= ix 0) (fx>= iy 0) (fx<= ex M) (fx<= ey N)
                (fx<= ix ex) (fx<= iy ey)))
      (matrix-utils:error 'fill-matrix! "invalid sub-matrix dimensions: " (list ix iy ex ey)))

  (cond ((= order RowMajor)
         (fold-ec f0 (:parallel (:range b (fx* N ix) (fx* N M) N) (:range x ix ex)) (cons x b)
                  (lambda (x+b ax)
                    (fold-ec ax (:range y iy ey) y
                             (lambda (y ax)
                               (let ((i (fx- (car x+b) ix)) (j (fx- y iy)))
                                 (let-values (((val ax1) (f i j ax)))
                                   ((numvector-set! vobj) A (fx+ (cdr x+b) y) val)
                                   ax1)))))))
	   
        ((= order ColMajor)
         (fold-ec f0 (:parallel (:range b (fx* N ix) (fx* N M) M) (:range y iy ey)) (cons y b)
                  (lambda (y+b ax)
                    (fold-ec ax (:range x ix ex) x
                             (lambda (x ax)
                               (let ((i (fx- x ix)) (j (fx- (car y+b) iy)))
                                 (let-values (((val ax1) (f j i ax)))
                                   ((numvector-set! vobj) A (fx+ (cdr y+b) x) val)
                                   ax1)))))))
        
        (else (matrix-utils:error 'fill-matrix! "unknown order " order)))
  A)

;;
;; ONES:: VOBJ * M * N [* ORDER] -> A
;;
;; Procedure ONES returns a matrix A of size M x N, in which all
;; elements are 1.  Optional argument ORDER specifies the matrix
;; layout: ColMajor or RowMajor, default is RowMajor.
;; VOBJ is an implementation of the vector interface.
;;
(define (matrix-ones vobj m n #!key (order RowMajor))
  (let ((A ((numvector-new vobj) (* m n))))
    (fill-matrix! vobj order m n A (lambda (i j ax) (values 1.0 #f)) #f)))

;;
;; ZEROS:: VOBJ * M * N [* ORDER] -> A
;;
;; Procedure ZEROS returns a matrix A of size M x N, in which all
;; elements are 0.  Optional argument ORDER specifies the matrix
;; layout: ColMajor or RowMajor, default is RowMajor.
;; VOBJ is an implementation of the vector interface.
;;
(define (matrix-zeros vobj m n #!key (order RowMajor))
  (let ((A ((numvector-new vobj) (* m n))))
    (fill-matrix! vobj order m n A (lambda (i j ax) (values 0.0 #f)) #f)))


;;
;; EYE:: VOBJ * M * N [* ORDER] -> I
;;
;; Procedure EYE returns an identity matrix of size M x N.
;; Optional argument ORDER specifies the matrix layout: ColMajor
;; or RowMajor, default is RowMajor.
;;
;; MAKE-VECTOR is one of the homogeneous vector creation procedures
;; from SRFI-4, and FILL-MATRIX! is a procedure created by
;; MAKE-FILL-MATRIX, above. FILL-MATRIX! must operate on the same type
;; of vector as MAKE-VECTOR.
;;
(define (matrix-eye vobj m n #!key (order RowMajor))
  (let ((A ((numvector-new vobj) (* m n))))
    (fill-matrix! vobj order m n A (lambda (i j ax) (values (if (fx= i j) 1.0 0.0) #f)) #f)))

;;
;; DIAG:: M * N * V [* K * ORDER] -> D
;;
;; Procedure DIAG returns a diagonal matrix D of size M x N, with
;; vector V on diagonal K.  Argument K is optional.  If it is
;; positive, the vector is placed on the K-th super-diagonal of matrix
;; D.  If it is negative, it is placed on the -K-th sub-diagonal of
;; matrix D.  The default value of K is 0, and the vector is placed on
;; the main diagonal of matrix D. Optional argument ORDER specifies
;; the matrix layout: ColMajor or RowMajor, default is
;; RowMajor. Vector V is always assumed to be a row vector.
;;
(define (matrix-diag vobj m n v #!key (k 0) (order RowMajor))
  (let ((A  ((numvector-new vobj) (* m n)))
        (k  (if (eq? order RowMajor) k (- k))))
    (fill-matrix! vobj order m n A
                  (lambda (i j vi) 
                    (if (fx= k (fx- j i))
                        (values ((numvector-get vobj) v vi) (fx+ 1 vi))
                        (values 0.0 vi)))
                  0)))

;;
;; LINSPACE:: N * BASE * LIMIT -> V
;;
;; Procedure LINSPACE returns a row vector with N linearly spaced elements
;; between BASE and LIMIT.  The number of elements, N, must be greater
;; than 1.  The BASE and LIMIT are always included in the range.  If
;; BASE is greater than LIMIT, the elements are stored in decreasing
;; order.
;;
(define (linspace vobj n base limit)
    (if (not (> n 1)) 
	(matrix-utils:error 'linspace "vector size N must be greater than 1"))
    (let ((step  (/ (- limit base) (fx- n 1)))
	  (a     ((numvector-new vobj) n)))
      (fill-matrix! vobj RowMajor 1 n a 
		    (lambda (i j ax) (values (+ base (* 1.0 j step)) ax))  #f)))
      
;;
;; LOGSPACE:: N * BASE * LIMIT -> V
;;
;; Procedure LOGSPACE returns a row vector with elements that are
;; logarithmically spaced from 10^BASE to 10^LIMIT. The number of
;; elements, N, must be greater than 1.  The BASE and LIMIT are always
;; included in the range.  If BASE is greater than LIMIT, the elements
;; are stored in decreasing order.
;;
(define (logspace vobj n base limit)
    (if (not (> n 1)) 
	(matrix-utils:error 'logspace "vector size N must be greater than 1"))
    (let ((step  (/ (- limit base) (fx- n 1)))
	  (a     ((numvector-new vobj) n))
	  (b     (linspace vobj n base limit)))
      (fill-matrix! vobj RowMajor 1 n a 
		    (lambda (i j ax) 
		      (let ((v  (expt 10 ((numvector-get vobj) b j))))
			(values v  ax))) #f)))



;;
;; MATRIX-FOLD-PARTIAL:: M * N * A * F * P * X0 * [IX * IY * EX * EY] -> XN
;;
;; Procedure MATRIX-FOLD-PARTIAL applies the fold operation on a
;; matrix A of size M x N with the values returned by applying
;; procedure F to each pair of indices and the corresponding value at
;; that position in the matrix. MATRIX-FOLD-PARTIAL first applies the
;; predicate P to the indices, and if P returns true, then F is
;; applied.
;;
;; Procedure F is of the form LAMBDA V AX -> AX1, where V is a
;; matrix element at position (I,J) and AX is accumulator value. The
;; initial value of AX is given by X0. Procedure F is expected to
;; return the new accumulator value.
;;
;; Procedure P is of the form LAMBDA I J -> boolean, where I and J are
;; matrix indices.
;;
;; Optional arguments IX IY EX EY may specify a sub-matrix in matrix A
;; to be folded. These arguments are checked to make sure they specify
;; a valid sub-matrix.
;;

(define (matrix-fold-partial vobj order M N A f p x0 #!key (ix 0) (iy 0) (ex M) (ey N))

    ;; optional arguments to specify a sub-matrix 
  (if (not (and (fx>= ix 0) (fx>= iy 0) (fx<= ex M) (fx<= ey N)
                (fx<= ix ex) (fx<= iy ey)))
      (matrix-utils:error 'matrix-fold-partial "invalid sub-matrix dimensions: " (list ix iy ex ey)))
  
  (cond ((= order RowMajor)
         (fold-ec x0 (:parallel (:range b (fx* N ix) (fx* N M) N) (:range x ix ex)) (cons x b)
                  (lambda (x+b ax)
                    (fold-ec ax (:range y iy ey) y
                             (lambda (y ax)
                               (let ((i (fx- (car x+b) ix)) (j (fx- y iy)))
                                 (if (p i j)
                                     (f ((numvector-get vobj) A (fx+ (cdr x+b) y)) ax) ax)))))))
        
        ((= order ColMajor)
         (fold-ec x0 (:parallel (:range b (fx* N ix) (fx* N M) M) (:range y iy ey)) (cons y b)
                  (lambda (y+b ax)
                    (fold-ec ax (:range x ix ex) x
                             (lambda (x ax)
                               (let ((i (fx- x ix)) (j (fx- (car y+b) iy)))
                                 (if (p j i)
                                     (f ((numvector-get vobj) A (fx+ (cdr y+b) x)) ax) ax)))))))
        
        (else (matrix-utils:error 'matrix-fold-partial "unknown order " order))
        
        ))
     
  

;;
;; MATRIX-FOLD:: M * N * A * F * X0 * [IX * IY * EX * EY] -> XN
;;
;; Procedure MATRIX-FOLD applies the fold operation on a matrix A of
;; size M x N with the values returned by applying procedure F to each
;; pair of indices and the corresponding value at that position in the
;; matrix.
;;
;; Procedure F is of the form LAMBDA I J V AX -> AX1, where V is a
;; matrix element at position (I,J) and AX is accumulator value. The
;; initial value of AX is given by X0. Procedure F is expected to
;; return the new accumulator value.
;;
;; Optional arguments IX IY EX EY may specify a sub-matrix in matrix A
;; to be filled. These arguments are checked to make sure they specify
;; a valid sub-matrix.
;;

(define (matrix-fold vobj order M N A f x0 #!key (ix 0) (iy 0) (ex M) (ey N))

  (if (not (and (fx>= ix 0) (fx>= iy 0) (fx<= ex M) (fx<= ey N)
                (fx<= ix ex) (fx<= iy ey)))
      (matrix-utils:error 'matrix-fold "invalid sub-matrix dimensions: " (list ix iy ex ey)))

  (cond ((= order RowMajor)
         (fold-ec x0 (:parallel (:range b (fx* N ix) (fx* N M) N) (:range x ix ex)) (cons x b)
                  (lambda (x+b ax)
                    (fold-ec ax (:range y iy ey) y
                             (lambda (y ax)
                               (let ((i (fx- (car x+b) ix)) (j (fx- y iy)))
                                 (f i j ((numvector-get vobj) A (fx+ (cdr x+b) y)) ax) ax))))))
        
        ((= order ColMajor)
         (fold-ec x0 (:parallel (:range b (fx* N ix) (fx* N M) M) (:range y iy ey)) (cons y b)
                  (lambda (y+b ax)
                    (fold-ec ax (:range x ix ex) x
                             (lambda (x ax)
                               (let ((i (fx- x ix)) (j (fx- (car y+b) iy)))
                                 (f j i ((numvector-get vobj) A (fx+ (cdr y+b) x)) ax) ax))))))
        
        (else (matrix-utils:error 'matrix-fold "unknown order " order))
        ))
     

  
(define (matrix-map vobj order M N A f  #!key (ix 0) (iy 0) (ex M) (ey N))

     (if (not (and (fx>= ix 0) (fx>= iy 0) (fx<= ex M) (fx<= ey N)
		   (fx<= ix ex) (fx<= iy ey)))
	 (matrix-utils:error 'matrix-map "invalid sub-matrix dimensions: " (list ix iy ex ey)))
     
     (cond ((= order RowMajor)
	    (fill-matrix! vobj order M N A (lambda (i j ax) 
                                             (let ((v  ((numvector-get vobj) A (fx+ j (fx* i M))))) (values (f i j v) #f))) #f))
	   ((= order ColMajor)
	    (fill-matrix! vobj order M N A (lambda (i j ax) 
                                             (let ((v  ((numvector-get vobj) A (fx+ i (fx* j N))))) (values (f j i v) #f))) #f))
           
	   (else (matrix-utils:error 'matrix-map "unknown order " order))
           
           ))
			    




(define (repmat vobj src dims reps)
  
  (let* ((sz (numvector-itemsize vobj))
         (ndim (length dims))
	 (dimsize (make-s32vector ndim))
	 (vdims (list->s32vector dims)))

    (print "src = " src)
    (print "dims = " dims)
    (print "reps = " reps)
    (assert (every positive? reps))
    
    ;; compute dimension sizes
    (s32vector-set! dimsize 0 (car dims))
    (let recur ((i 1) (dims (cdr dims)))
      (if (< i ndim)
	  (begin 
	    (s32vector-set! dimsize i (* (car dims) (s32vector-ref dimsize (- i 1))))
	    (recur (+ 1 i) (cdr dims)))
	  ))

    ;; determine repetition vector 
    (let* ((nrep (length reps))
	   (ndimdest (if (> nrep ndim) nrep ndim))
	   (rep (make-s32vector ndimdest)))
      (let recur ((i 0) (reps reps))
	(if (< i nrep)
	    (let ((repv (car reps)))
	      (s32vector-set! rep i repv)
	      (recur (+ 1 i) (cdr reps)))))

      ;; compute output size
      (let ((destdims (make-s32vector ndimdest 0)))
	(let recur ((i 0))
	  (if (< i ndim)
	      (begin (s32vector-set! destdims i 
				    (* (s32vector-ref vdims i)
				       (s32vector-ref rep i )))
		     (recur (+ 1 i)))))

	(let ((extrarep
	       (let recur ((i ndim) (extrarep 1))
		 (if (< i ndimdest)
		     (let ((ri (s32vector-ref rep i)))
		       (s32vector-set! destdims i ri)
		       (recur (+ 1 i) (* extrarep ri)))
		     extrarep
		     ))))

	  (let ((destdimsize (make-s32vector ndimdest)))
	    (s32vector-set! destdimsize 0 (s32vector-ref destdims 0))
	    (let recur ((i 1))
	      (if (< i ndimdest)
		  (begin
		    (s32vector-set! destdimsize i 
				    (* (s32vector-ref destdimsize (- i 1))
				       (s32vector-ref destdims i)))
		    (recur (+ 1 i)))))

	    ;; return array
	    (let* ((destlen (s32vector-ref destdimsize (- ndimdest 1)))
		   (dest ((numvector-new vobj) destlen)))
              (print "destlen = " destlen)
	      (vrepmat dest src sz ndim destdimsize dimsize vdims rep)
	      (if (> ndimdest ndim)
		  (let ((n (s32vector-ref destdimsize (- ndim 1))))
		    (vmemrep dest sz n extrarep)))
	      dest
	      ))
	  ))
      ))
  )


	  

(define (meshgrid vobj x y  . rest)
  (let-optionals rest ((z #f))
  
  (let ((sz   (numvector-itemsize vobj))
        (lenx ((numvector-length vobj) x))
        (leny ((numvector-length vobj) y))
        (lenz (and z ((numvector-length vobj) z))))
    
    (let ((dimx (list 1 lenx))
          (dimy (list 1 leny))
          (dimz (and z (list 1 lenz))))
          
      (let ((xx (if z (repmat vobj
                              (repmat vobj x dimx (list leny 1))
                              (list leny lenx)
                              (list 1 1 lenz))
                    (repmat vobj x dimx (list leny 1))))
            (yy (if z (repmat vobj
                              (repmat vobj y dimy (list 1 lenx))
                              (list lenx leny)
                              (list 1 1 lenz))
                (repmat vobj y dimy (list 1 lenx))))
            (zz (if z (repmat vobj z dimz (list (* lenx leny) 1)) #f)))
        
        (filter identity (list xx yy zz)))
      
      ))
  ))

      

      

)
