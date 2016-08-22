package breeze.linalg
package operators

import breeze.math.Ring
import breeze.math.Semiring
import breeze.storage.Zero
import breeze.collection.mutable.OpenAddressHashArray
import breeze.math.Field

import scala.reflect.ClassTag


//OpAdd, OpSub, OpMulScalar, OpMulMatrix, OpDiv, OpSet, OpMod, OpPow)

//operam elemento a elemento numa matriz e são caracterizadas pelos seguintes símbolos

//{_ + _},  {_ - _}, {_ * _}, {_ * _}, {_ / _}, {(a,b) => b}, {_ % _}, {_ pow _}

trait HashMatrixOps_Ring extends HashMatrixOpsLowPrio /*with SerializableLogging*/ {
  //this: CSCMatrixOps =>
implicit def hash_OpNeg[T:Ring]: OpNeg.Impl[HashMatrix[T], HashMatrix[T]] = {
    new OpNeg.Impl[HashMatrix[T], HashMatrix[T]] {
	  println("usei "+
			" implicit def hash_OpNeg[T:Ring:ClassTag]: OpNeg.Impl[HashMatrix[T], HashMatrix[T]]"+
			"de HashMatrixOps_Ring")
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

  //This is used for creating some vector spaces and MutableFiniteCoordinateField.make
  //There is some stuff defined explicitly for each matrix type in math\src\main\scala\breeze\math\VectorSpace.scala
  //It is not clear what is this used for precisely
  /*
  implicit def hashScaleAdd[T: Semiring]: scaleAdd.InPlaceImpl3[HashMatrix[T], T, HashMatrix[T]] = {
    new scaleAdd.InPlaceImpl3[HashMatrix[T], T, HashMatrix[T]] {
      override def apply(a: HashMatrix[T], s: T, b: HashMatrix[T]): Unit = {
	  println("usei "+
			" implicit def hashScaleAdd[T: Semiring : ClassTag]: scaleAdd.InPlaceImpl3[HashMatrix[T], T, HashMatrix[T]]"+
			"de HashMatrixOps_Ring")
        val ring = implicitly[Semiring[T]]
        require(a.rows == b.rows, "Matrices must have same number of rows!")
        require(a.cols == b.cols, "Matrices must have same number of cols!")
        val rows = a.rows
        val cols = a.cols

        if (cols == 0 || rows == 0) return

		for(((f,t),v) <- b.activeIterator){
			a(f,t) = ring.+(a(f,t), ring.*(s, v))
		}
      }
    }
  }
  */

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
		//val re = dm + v2
		//re
		dm + v2
      }
    }


}
//I don't know what is this used for
//this trait is pure cargo cult programming
trait HashMatrixOpsLowPrio {
  //this: CSCMatrixOps =>
  implicit def canMulM_V_def[T, A, B <: Vector[T]](implicit bb: B <:< Vector[T], op: OpMulMatrix.Impl2[HashMatrix[T], Vector[T], Vector[T]]) ={
	println("usei "+
			"  implicit def canMulM_V_def[T, A, B <: Vector[T]](implicit bb: B <:< Vector[T], op: OpMulMatrix.Impl2[HashMatrix[T], Vector[T], Vector[T]]) ={"+
			"de HashMatrixOps")
    implicitly[OpMulMatrix.Impl2[HashMatrix[T], Vector[T], Vector[T]]].asInstanceOf[breeze.linalg.operators.OpMulMatrix.Impl2[A, B, Vector[T]]]
}
  // ibid.
  

implicit def canMulM_M_def[T, B <: Matrix[T]](implicit bb: B <:< Matrix[T], op: OpMulMatrix.Impl2[HashMatrix[T], Matrix[T], HashMatrix[T]]) ={
	println("usei "+
		"canMulM_M_def[T, B <: Matrix[T]](implicit bb: B <:< Matrix[T], op: OpMulMatrix.Impl2[HashMatrix[T], Matrix[T], HashMatrix[T]])"
		+"de HashMatrixOps")
    op.asInstanceOf[OpMulMatrix.Impl2[HashMatrix[T], B, HashMatrix[T]]]
	}

}

