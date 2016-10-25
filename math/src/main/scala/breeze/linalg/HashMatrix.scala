package breeze.linalg

//import breeze.linalg.HashVector
import scala.{specialized => spec}
import breeze.collection.mutable.OpenAddressHashArray
import breeze.storage.Zero

import scala.reflect.ClassTag
import breeze.linalg.operators.HashMatrixOpsLowPrio
import breeze.linalg.operators.MatrixOps
import breeze.linalg.operators.HashMatrixOps_Ring
import breeze.linalg.support.{CanMapValues, ScalarOf}
import breeze.math.{Field, Semiring}

import scalaxy.debug._

class HashMatrix[@spec(Double, Int, Float, Long) V: Zero](val harray: OpenAddressHashArray[V],
                                                         val rows: Int,
                                                         val cols: Int) 
  extends Matrix[V] with MatrixLike[V, HashMatrix[V]] with Serializable{

  def apply(row: Int, col: Int): V = {
    if(row < - rows || row >= rows) throw new IndexOutOfBoundsException((row,col) + " not in [-"+rows+","+rows+") x [-"+cols+"," + cols+")")
    if(col < - cols || col >= cols) throw new IndexOutOfBoundsException((row,col) + " not in [-"+rows+","+rows+") x [-"+cols+"," + cols+")")
    val trueRow = if(row<0) row + rows else row
    val trueCol = if(col<0) col + cols else col
    harray(linearIndex(trueRow, trueCol))
  }

  def linearIndex(row: Int, col: Int): Int = {
      row + col * rows
  }

 def cartesianIndex(lindex: Int): (Int,Int) = {
      val col = lindex / rows
	  val row = lindex - col * rows
     (row , col)
  }

  def update(row: Int, col: Int, e: V): Unit = {
    if(row < - rows || row >= rows) throw new IndexOutOfBoundsException((row,col) + " not in [-"+rows+","+rows+") x [-"+cols+"," + cols+")")
    if(col < - cols || col >= cols) throw new IndexOutOfBoundsException((row,col) + " not in [-"+rows+","+rows+") x [-"+cols+"," + cols+")")
    val trueRow = if(row<0) row + rows else row
    val trueCol = if(col<0) col + cols else col
    harray(linearIndex(trueRow, trueCol)) = e
}



  def flatten(view: View=View.Prefer): HashVector[V] = view match {
    case View.Require =>
      new HashVector(harray)
    case View.Copy =>
      new HashVector(harray.copy)
    case View.Prefer =>
      new HashVector(harray)
}


  def copy: HashMatrix[V] = new HashMatrix(harray.copy , rows , cols)

  def repr: HashMatrix[V] = this

  def activeIterator: Iterator[((Int, Int), V)] = {
	harray.activeIterator.map({x => (cartesianIndex(x._1),x._2)})
  }

  def activeKeysIterator: Iterator[(Int, Int)] = {
    harray.activeKeysIterator.map( cartesianIndex _)
  }

  def activeValuesIterator: Iterator[V] =  harray.activeValuesIterator

  def activeSize: Int = harray.activeSize

  override def toString(maxLines: Int, maxWidth: Int): String = {
    val buf = new StringBuilder()
    buf ++= ("%d x %d HashMatrix".format(rows, cols))
    activeIterator.take(maxLines - 1).foreach { case ((r,c),v) =>
      buf += '\n'
      buf ++= "(%d,%d) ".format(r,c)
      buf ++= v.toString
    }
    buf.toString()
  }

  def zero = implicitly[Zero[V]].zero

  //TODO break hashmatrix in two parts : a algebric and a pure data container
  //the methodo above could then be broke into a zip and a separate map
  def zipMapHeterogeneousVals[S , R: Field  : ClassTag : Semiring : Zero](that : HashMatrix[S],  fn: (V, S) => R ) : HashMatrix[R] = {
    require(rows == that.rows, "Matrices must have same number of rows!")
    require(cols == that.cols, "Matrices must have same number of cols!")

    val res =  HashMatrix.zeros[R](rows, cols)
    val rzero = implicitly[Zero[R]].zero

    val fn00 = fn(this.zero,that.zero)
    if(fn00!=rzero) { //a no-sparse operation, to be avoided
      res := fn00
    }
    val commonActive = (activeKeysIterator.toSet).union(that.activeKeysIterator.toSet)
    for((f,t) <- commonActive){
      res(f,t) = fn(this(f,t),that(f,t))
    }
    res
  }
}


object HashMatrix extends MatrixConstructors[HashMatrix] with HashMatrixOps_Ring/*with MatrixOpswith HashMatrixOpsLowPrio*/{


  def zeros[@spec(Double, Int, Float, Long) V:ClassTag:Zero](rows: Int, cols: Int): HashMatrix[V] = {
	new HashMatrix[V](new OpenAddressHashArray[V](rows * cols),rows,cols)
  }

	def create[@spec(Double, Int, Float, Long) V: Zero](rows: Int, cols: Int, data: Array[V]): HashMatrix[V] = {
		val z = implicitly[Zero[V]].zero
		implicit val man = ClassTag[V](data.getClass.getComponentType.asInstanceOf[Class[V]])
		val res = zeros(rows, cols)
		var i = 0
		for(c <- 0 until cols; r <- 0 until rows) {
		  val v = data(i)
		  i += 1
		  if ( v != z) {
			res(r, c) = v
		  }
		}
		res
	}

  implicit def canMapValues[V: Zero, R:Field:ClassTag:Zero:Semiring]:CanMapValues[HashMatrix[V], V, R, HashMatrix[R]] = {
    val z = implicitly[Zero[R]].zero
    new CanMapValues[HashMatrix[V], V, R, HashMatrix[R]] {
      override def apply(from: HashMatrix[V], fn: (V => R)) = {
        val fz = fn(implicitly[Zero[V]].zero)
        val fzIsNotZero = fz != z
        val res = HashMatrix.zeros[R](from.rows, from.cols)
        if(fzIsNotZero){
          res := fz
        }
        for(((f,t),v) <- from.activeIterator){
          res(f,t) = fn(v)
        }
        res
      }
    }
  }
  //TODO shouldn't this be in Matrix?
  implicit def scalarOf[T]: ScalarOf[HashMatrix[T], T] = ScalarOf.dummy

}
