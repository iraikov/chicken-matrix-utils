(use interfaces srfi-1 srfi-4 blas matrix-utils test)

(define order RowMajor)

(define f64numvector
  (implementation
   numvector
   (define numvector-itemsize 8)
   (define numvector-new    make-f64vector)
   (define numvector-set!   f64vector-set!)
   (define numvector-get    f64vector-ref)
   (define numvector-length f64vector-length)
   (define list->numvector list->f64vector)
   ))

(define f32numvector
  (implementation
   numvector
   (define numvector-itemsize 4)
   (define numvector-new    make-f32vector)
   (define numvector-set!   f32vector-set!)
   (define numvector-get    f32vector-ref)
   (define numvector-length f32vector-length)
   (define list->numvector list->f32vector)
   ))


(define a (matrix-ones  f64numvector 10 8))
(define b (matrix-zeros f64numvector 10 8))
(define c (matrix-eye   f64numvector 10 8))
(define d (matrix-diag  f64numvector 10 10 a))

(print a)
(print b)
(print c)
(print d)

(define aa (linspace f32numvector 10 0 8))
(define bb (logspace f32numvector 10 0 8))

(print aa)
(print bb)

(define x ((list->numvector f32numvector) '(1 2 3)))

(test-group "repmat"
            (test (repmat f32numvector x '(3 1) '())
                  (repmat f32numvector x 3 '())))

