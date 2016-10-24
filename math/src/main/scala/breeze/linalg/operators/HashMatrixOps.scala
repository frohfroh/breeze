package breeze.linalg
package operators

import breeze.math.Ring
import breeze.math.Semiring
import breeze.storage.Zero
import breeze.collection.mutable.OpenAddressHashArray
import breeze.generic.UFunc
import breeze.generic.UFunc.InPlaceImpl2
import breeze.math.Field
import breeze.macros.expand

import scalaxy.debug._
//import breeze.generic.UFunc

import scala.reflect.ClassTag



//OpAdd, OpSub, OpMulScalar, OpMulMatrix, OpDiv, OpSet, OpMod, OpPow)

//operam elemento a elemento numa matriz e são caracterizadas pelos seguintes símbolos

//{_ + _},  {_ - _}, {_ * _}, {_ * _}, {_ / _}, {(a,b) => b}, {_ % _}, {_ pow _}

trait HashMatrixOps_Ring extends HashMatrixOpsLowPrio /*with SerializableLogging*/ {
  //this: CSCMatrixOps =>
//implements - HashMatrix 
implicit def hash_OpNeg[T:Ring]: OpNeg.Impl[HashMatrix[T], HashMatrix[T]] = {
    new OpNeg.Impl[HashMatrix[T], HashMatrix[T]] {
      val ring = implicitly[Ring[T]]
      def apply(a: HashMatrix[T]): HashMatrix[T] = {
        val acp = a.copy
		val harray : OpenAddressHashArray[T] = acp.harray
		val data : Array[T] = harray.data
		var i = 0
		while(i<data.length){
			data(i) = ring.negate(data(i))
			i=i+1
		}
        acp
      }
    }
  }


  //used for multiplying a HashMatrix with a vector
  implicit def canMulM_V_Semiring[T:Semiring:Zero:ClassTag]: BinaryRegistry[HashMatrix[T], Vector[T],OpMulMatrix.type, Vector[T]] =
    new BinaryRegistry[HashMatrix[T], Vector[T], OpMulMatrix.type, Vector[T]] {
    implicit val ring = implicitly[Semiring[T]]

    override def bindingMissing(a: HashMatrix[T], b: Vector[T]) = {
      require(a.cols == b.length, "Dimension Mismatch!")
      val res = HashVector.zeros[T](a.rows)
	  for(((f,t),v) <- a.activeIterator){
			b(f) = ring.+(b(f) , ring.*(b(t),v)) 
	  }
      res
    }
  }


  //probably not good for multiplying DenseMatrix 
  //implements HashMatrix * Matrix 
  implicit def canMulM_M_Semiring[T: Semiring : Zero : ClassTag]: OpMulMatrix.Impl2[HashMatrix[T], Matrix[T], HashMatrix[T]] =
    new OpMulMatrix.Impl2[HashMatrix[T], Matrix[T], HashMatrix[T]] {
      def apply(a: HashMatrix[T], b: Matrix[T]) = {
        val ring = implicitly[Semiring[T]]
        require(a.cols == b.rows, "HashMatrix Multiplication Dimension Mismatch")
        val res = HashMatrix.zeros[T](a.rows, b.cols)

		val byCol = a.activeIterator.toList.groupBy(_._1._2)

		for(((f_b,t_b),v_b) <- b.activeIterator){
			for(((f_a,_),v_a) <- byCol.getOrElse(f_b , Nil)){
				res(f_a,t_b) = ring.+(res(f_a,t_b) , ring.*(v_b,v_a))
			}
		}
        
		res
      }
    }

//implements HashMatrix + scalar
//If you're adding something, you accept loosing sparcity
//In the rare case where you're adding zero, if you want to keep sparcity, you must test it for yourself
  implicit def canAddM_S_Semiring[T: Semiring : ClassTag :Field:Zero]: OpAdd.Impl2[HashMatrix[T], T, DenseMatrix[T]] =
    new OpAdd.Impl2[HashMatrix[T], T, DenseMatrix[T]] {
        val s = implicitly[Semiring[T]]
        val zero = s.zero
      override def apply(v: HashMatrix[T], v2: T): DenseMatrix[T] = {
		val dm = v.toDenseMatrix
		dm + v2
      }
    }

//Same thing as above for HashMatrix - scalar
  implicit def canSubM_S_Ring[T: Semiring : ClassTag :Field:Zero]: OpSub.Impl2[HashMatrix[T], T, DenseMatrix[T]] =
    new OpSub.Impl2[HashMatrix[T], T, DenseMatrix[T]] {
        val s = implicitly[Semiring[T]]
        val zero = s.zero
      override def apply(v: HashMatrix[T], v2: T): DenseMatrix[T] = {
		val dm = v.toDenseMatrix
		dm - v2
      }
    }





  //used implicitly in Hash_T_UpdateOp
  implicit def canSetM_S_Semiring[T: Semiring : ClassTag]: OpSet.Impl2[HashMatrix[T], T, HashMatrix[T]] =
    new OpSet.Impl2[HashMatrix[T], T, HashMatrix[T]] {
      val r = implicitly[Semiring[T]]
      val zero = r.zero
      def apply(v: HashMatrix[T], v2: T): HashMatrix[T] = {
        val zm = HashMatrix.zeros[T](v.rows,v.cols)
        if (v2 == zero)
          return zm
        //we shouldn't be here for a sparse matrix, so we don't care about performance
        for(i <- 0 until v.rows){
          for(j <- 0 until v.cols){
            zm(i,j)   = v2
          }
        }
        zm
      }
    }

//implements HashMatrix * scalar
  @expand
  implicit def canMulM_S_Ring[@expand.args(OpMulMatrix,OpMulScalar) Op <: OpType, T:Ring:ClassTag]: Op.Impl2[HashMatrix[T],T,HashMatrix[T]] = {
    val ring = implicitly[Ring[T]]
    new Op.Impl2[HashMatrix[T], T, HashMatrix[T]] {
      def apply(hm: HashMatrix[T], v2: T): HashMatrix[T] = {
		val res = HashMatrix.zeros[T](hm.rows,hm.cols)
		for(((f,t),v) <- hm.activeIterator){
			res(f,t) = ring.*(v, v2) 
	  }
		res
      }
    }
  }
//implements HashMatrix *:* Matrix
  @expand
  implicit def HashMatrixCanMulScalarM_M_Semiring[@expand.args(Int, Long, Float, Double) A: Semiring : ClassTag : Zero , MM <: Matrix[A]]: OpMulScalar.Impl2[HashMatrix[A], MM, HashMatrix[A]] =
    new OpMulScalar.Impl2[HashMatrix[A], MM, HashMatrix[A]] {
      val ring = implicitly[Semiring[A]]
      final def apply(a: HashMatrix[A], b: MM): HashMatrix[A] = {
        val rows = a.rows
        val cols = a.cols
        require(rows == b.rows, "Matrices must have same number of rows!")
        require(cols == b.cols, "Matrices must have same number of cols!")
		val res = HashMatrix.zeros[A](a.rows,a.cols)
		for(((f,t),v) <- a.activeIterator){
			res(f,t) = ring.*(v, b(f,t)) 
	  	}
        res
      }
    }

//implements HashMatrix + HashMatrix
 implicit def HashMatrixCanAdd_M_M_Semiring[A:Semiring:Zero:ClassTag]: OpAdd.Impl2[HashMatrix[A], HashMatrix[A], HashMatrix[A]] =
    new OpAdd.Impl2[HashMatrix[A], HashMatrix[A], HashMatrix[A]] {
    val ring = implicitly[Semiring[A]]
    def apply(a: HashMatrix[A], b: HashMatrix[A]): HashMatrix[A] = {
      require(a.rows == b.rows, "Matrix dimensions must match")
      require(a.cols == b.cols, "Matrix dimensions must match")
	val res = a.copy
	for(((f,t),v) <- b.activeIterator){
		res(f,t) = ring.+(v, a(f,t)) 
	}
	res
    }
  }

  //implements HashMatrix - HashMatrix
  implicit def HashMatrixCanSub_M_M_Ring[A:Ring:Zero:ClassTag]: OpSub.Impl2[HashMatrix[A], HashMatrix[A], HashMatrix[A]] =
  new OpSub.Impl2[HashMatrix[A], HashMatrix[A], HashMatrix[A]] {
    val ring = implicitly[Ring[A]]
    def apply(a: HashMatrix[A], b: HashMatrix[A]): HashMatrix[A] = {
      require(a.rows == b.rows, "Matrix dimensions must match")
      require(a.cols == b.cols, "Matrix dimensions must match")
      val res = a.copy
      for(((f,t),v) <- b.activeIterator){
        res(f,t) = ring.-(v, a(f,t))
      }
      res
    }
  }

  //used for HashMatrix / scalar , HashMatrix % scalar  and HashMatrix ^:^ scalar
  @expand
  implicit def hash_T_Op[@expand.args(OpDiv, OpMod, OpPow) Op <: OpType, T:Field:ClassTag]
  (implicit @expand.sequence[Op]({f./(_,_)}, {f.%(_,_)},{f.pow(_,_)}) op: Op.Impl2[T,T,T]):
  Op.Impl2[HashMatrix[T], T, HashMatrix[T]] = {
    val f = implicitly[Field[T]]
    new Op.Impl2[HashMatrix[T], T, HashMatrix[T]] {
      def apply(a: HashMatrix[T], b: T): HashMatrix[T] = {
        if (b == f.zero) {
          // degenerate case, creates effectively dense matrix
          val default = op(f.zero, b)
          val res =  HashMatrix.zeros[T](a.rows,a.cols)
          for(i <- 0 until a.rows ; j <- 0 until a.cols){
            res(i,j) = default
          }
          //don't care about performance, one shouldn't use HashMatrix with this ops and zero anyway
          for(((fr,t),v) <- a.activeIterator){
            res(fr,t) = op(v, b)
          }
          res
        } else {
          val res =  HashMatrix.zeros[T](a.rows,a.cols)
          for(((fr,t),v) <- a.activeIterator){
            res(fr,t) = op(v, b)
          }
          res
        }
      }
    }
  }


  implicit def HashMatrixCanSetM_M_Semiring[T:Semiring:ClassTag]: OpSet.Impl2[HashMatrix[T], HashMatrix[T], HashMatrix[T]] = {
    val f = implicitly[Semiring[T]]
    new OpSet.Impl2[HashMatrix[T], HashMatrix[T], HashMatrix[T]] {
      def apply(a: HashMatrix[T], b: HashMatrix[T]): HashMatrix[T] = {
        val rows  = a.rows
        val cols  = a.cols
        require(rows == b.rows, "Matrices must have same number of rows!")
        require(cols == b.cols, "Matrices must have same number of cols!")
        b.copy
      }
    }
  }

  //should be used for HashMatrix := Matrix
  //which is used in turn for HashMatrix *= scalar & other similar ones
  //unfourtunately a new matrix is created and all the non-zero values are copied unnecessarily
  protected def updateFromPure_Hash_T[T:Zero, Op<:OpType, Other , MT <: Matrix[T]](implicit op: UFunc.UImpl2[Op, HashMatrix[T], Other, MT])
  : UFunc.InPlaceImpl2[Op, HashMatrix[T], Other] = {
    val zero = implicitly[Zero[T]]
    new UFunc.InPlaceImpl2[Op, HashMatrix[T], Other] {
      def apply(a: HashMatrix[T], b: Other) {
        val result: MT = op(a, b)
        //inefficient, but as written above, without a legitimate use case...
        //We could use var for the arguments of a HashMatrix so to avoid copying from results to a...
        //which is more or less what is done in CSCMatrix
        for(((fr,t),v) <- a.activeIterator){
          a(fr,t) = zero.zero
        }
        for(((fr,t),v) <- result.activeIterator){
          a(fr,t) = v
        }
      }
    }
  }

 //used for HashMatrix *= scalar & other similar ones
  //for sparsity-losing operations, the not in place versions gives a DenseMatrix (eg Hashmatrix + Scalar =  DenseMatrix)
  //the inplace operations must "keep" the HashMatrix, but use the not in place operations, copying the elements
  //The zeros that may appear are not recovered and keep on wasting memory and time
  @expand
  implicit def Hash_T_UpdateOpG[@expand.args(OpMulMatrix, OpSet, OpSub, OpAdd, OpMulScalar, OpDiv, OpMod, OpPow) Op <: OpType, T:Field:ClassTag]
  : Op.InPlaceImpl2[HashMatrix[T],T] = {
    updateFromPure_Hash_T( implicitly[Zero[T]]  , implicitly[Op.Impl2[HashMatrix[T], T, Matrix[T]]])
  }

  //used for HashMatrix := HashMatrix
  //not sure what the other expansions were used for in CSCMatrix, so kept the expand macro format for future implementation when I figure it out
  @expand
  implicit def Hash_Hash_UpdateOpG[@expand.args(OpSet) Op <: OpType, T:Field:ClassTag]
  : Op.InPlaceImpl2[HashMatrix[T],HashMatrix[T]] = {
    updateFromPure_Hash_T( implicitly[Zero[T]]  , implicitly[Op.Impl2[HashMatrix[T], HashMatrix[T], HashMatrix[T]]])
  }

  def test = {
    import breeze.linalg._
    val hmd = HashMatrix.zeros[Double](3,3)
    hmd(1,1) = 1.0
    val dmd = hmd.toDenseMatrix
    val dv = DenseVector.zeros[Double](3)
    val hmd2 = hmd.copy
    //new BinaryUpdateRegistry[Matrix[T], Matrix[T], Op.type]
   hmd.:=(hmd2)
  }

}

trait HashMatrixOpsLowPrio {
  //this: CSCMatrixOps =>
  implicit def canMulM_V_def[T, A, B <: Vector[T]](implicit bb: B <:< Vector[T], op: OpMulMatrix.Impl2[HashMatrix[T], Vector[T], Vector[T]]) ={
    implicitly[OpMulMatrix.Impl2[HashMatrix[T], Vector[T], Vector[T]]].asInstanceOf[breeze.linalg.operators.OpMulMatrix.Impl2[A, B, Vector[T]]]
}
  // ibid.

implicit def canMulM_M_def[T, B <: Matrix[T]](implicit bb: B <:< Matrix[T], op: OpMulMatrix.Impl2[HashMatrix[T], Matrix[T], HashMatrix[T]]) ={
    op.asInstanceOf[OpMulMatrix.Impl2[HashMatrix[T], B, HashMatrix[T]]]
	}

}

