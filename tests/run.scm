(use typeclass srfi-1 srfi-4 blas matrix-utils test)

(define order RowMajor)

(define f64-Vector
  (let (
        (vect-itemsize 8)
        (vect-new      make-f64vector)
        (vect-update   f64vector-set!)
        (vect-get      f64vector-ref)
        (vect-size     f64vector-length)
        (from-list     list->f64vector)
        )
    (make-<Vector> vect-itemsize
                   vect-new from-list
                   vect-get vect-update 
                   vect-size )
    ))

(define f32-Vector
  (let (
        (vect-itemsize 4)
        (vect-new      make-f32vector)
        (vect-update   f32vector-set!)
        (vect-get      f32vector-ref)
        (vect-size     f32vector-length)
        (from-list     list->f32vector)
        )
    (make-<Vector> vect-itemsize
                   vect-new from-list
                   vect-get vect-update 
                   vect-size )
    ))


(define a ((matrix-ones  f64-Vector) 10 8))
(define b ((matrix-zeros f64-Vector) 10 8))
(define c ((matrix-eye   f64-Vector) 10 8))
(define d ((matrix-diag  f64-Vector) 10 10 a))

(print a)
(print b)
(print c)
(print d)

(define aa ((linspace f64-Vector) 10 0 8))
(define bb ((logspace f64-Vector) 10 0 8))

(print aa)
(print bb)

(print (with-instance
        ((<Vector> f64-Vector))
        (get aa 0)))

(define x (with-instance
           ((<Vector> f32-Vector))
           (from-list (list 1 2 3))))

(test-group "repmat"
            (test ((repmat f32-Vector) x '(3 1) '())
                  ((repmat f32-Vector) x 3 '())))

