package breeze.linalg


import breeze.linalg.operators.OpMulMatrix
import breeze.math.{Complex, Field}
import breeze.storage.Zero
import org.scalatest._
import org.scalatest.junit._
import org.scalatest.prop._
import org.junit.runner.RunWith

import scala.reflect.ClassTag

@RunWith(classOf[JUnitRunner])
class HashMatrixTest extends FunSuite with Checkers {
  test("Multiply") {
    val a = HashMatrix((1.0, 2.0, 3.0),(4.0, 5.0, 6.0))
    val ad = DenseMatrix((1.0, 2.0, 3.0),(4.0, 5.0, 6.0))
    val b = HashMatrix((7.0, -2.0, 8.0),(-3.0, -3.0, 1.0),(12.0, 0.0, 5.0))
    val bd = DenseMatrix((7.0, -2.0, 8.0),(-3.0, -3.0, 1.0),(12.0, 0.0, 5.0))
    val c = DenseVector(6.0,2.0,3.0)

    assert( (a * b: HashMatrix[Double]) === HashMatrix((37.0, -8.0, 25.0), (85.0, -23.0, 67.0)))
    assert((a * bd :HashMatrix[Double])=== HashMatrix((37.0, -8.0, 25.0), (85.0, -23.0, 67.0)))
    assert((ad * b :DenseMatrix[Double])=== DenseMatrix((37.0, -8.0, 25.0), (85.0, -23.0, 67.0)))
    val prod = implicitly[OpMulMatrix.Impl2[HashMatrix[Double], DenseVector[Double], Vector[Double]]]
    assert(a * c === DenseVector(19.0,52.0))
    assert(b * c === DenseVector(62.0, -21.0, 87.0))


    //given HashMatrix has a somewhat different behaviour then CSCMatrix
    //because it is expected to be used in some extreme sparse ccases
    //the following tests are no longer valid
    //TODO Should one change the behaviour of HashMatrix or rewrite those tests?
    //    assert(b.t * c === DenseVector(72.0, -18.0, 65.0))
    //    assert(a.t * DenseVector(4.0, 3.0) === DenseVector(16.0, 23.0, 30.0))

        // should be dense
    //    val x = a * a.t
    //    assert(x === DenseMatrix((14.0,32.0),(32.0,77.0)))

        // should be dense
    //    val y = a.t * a
    //    assert(y === DenseMatrix((17.0,22.0,27.0),(22.0,29.0,36.0),(27.0,36.0,45.0)))

    //    val z : DenseMatrix[Double] = b * (b + 1.0)
    //    assert(z === DenseMatrix((164.0,5.0,107.0),(-5.0,10.0,-27.0),(161.0,-7.0,138.0)))
  }

  test("Multiply Int") {
    val a = HashMatrix((1, 2, 3),(4, 5, 6))
    val b = HashMatrix((7, -2, 8),(-3, -3, 1),(12, 0, 5))
    val bd = DenseMatrix((7, -2, 8),(-3, -3, 1),(12, 0, 5))
    val c = DenseVector(6,2,3)
    val cs = SparseVector(3)( (1,2))
    assert(a * b === HashMatrix((37, -8, 25), (85, -23, 67)))
    assert(a * bd === DenseMatrix((37, -8, 25), (85, -23, 67)))
    assert(a * c === DenseVector(19,52))
    assert(b * c === DenseVector(62, -21, 87))
    assert(a * cs === SparseVector(4, 10))
    assert(b * cs === SparseVector(3)((0, -4), (1, -6)))

    //given HashMatrix has a somewhat different behaviour then CSCMatrix
    //because it is expected to be used in some extreme sparse ccases
    //the following tests are no longer valid
    //TODO Should one change the behaviour of HashMatrix or rewrite those tests?

//    assert(b.t * c === DenseVector(72, -18, 65))
//    assert(a.t * DenseVector(4, 3) === DenseVector(16, 23, 30))

    // should be dense
//    val x = a * a.t
//    assert(x === DenseMatrix((14,32),(32,77)))

    // should be dense
//    val y = a.t * a
//    assert(y === DenseMatrix((17,22,27),(22,29,36),(27,36,45)))

//    val z : DenseMatrix[Double] = b * (b + 1.0)
//    assert(z === DenseMatrix((164,5,107),(-5,10,-27),(161,-7,138)))
  }

  test("Multiply Complex") {

    val a = HashMatrix((Complex(1,1), Complex(2,2), Complex(3,3)),
      (Complex(4,4), Complex(5,5), Complex(6,6)))
    val b = HashMatrix((Complex(7,7), Complex(-2,-2), Complex(8,8)),
      (Complex(-3,-3), Complex(-3,-3), Complex(1,1)),
      (Complex(12,12), Complex(0,0), Complex(5,5)))
    val c = DenseVector(Complex(6,0), Complex(2,0), Complex(3,0))
    val cs = SparseVector(Complex(6,0), Complex(2,0), Complex(3,0))
    val value: HashMatrix[Complex] = a * b
    assert(value === HashMatrix((Complex(0,74), Complex(0,-16), Complex(0,50)),
      (Complex(0,170), Complex(0,-46), Complex(0,134))))
    assert(b * c === DenseVector(Complex(62,62), Complex(-21,-21), Complex(87,87)))//TODO error
    assert(b * cs === DenseVector(Complex(62,62), Complex(-21,-21), Complex(87,87)))
    assert(b.t * c === DenseVector(Complex(72,-72), Complex(-18,18), Complex(65,-65)))
  }
  
  test("Transpose") {
    val a = HashMatrix.zeros[Int](2,3)
    a(0,0) = 1
    a(1,2) = 2
    
    val expected = HashMatrix.zeros[Int](3,2)
    expected(0,0) = 1
    expected(2,1) = 2
    
    assert(a.t === expected)//TODO error
  }
  
  test("Transpose Complex") {
    val a = HashMatrix.zeros[Complex](2,3)
    a(0,0) = Complex(1,1)
    a(1,2) = Complex(-2,-2)
    
    val expected = HashMatrix.zeros[Complex](3,2)
    expected(0,0) = Complex(1,-1)
    expected(2,1) = Complex(-2,2)
    
    assert(a.t === expected)
  }

  test("Generic Hash ops") {
    // mostly for coverage
    val a = HashMatrix.create[String](1,1, Array("SSS"))
    intercept[IndexOutOfBoundsException] {
      a(3,3) = ":("
      assert(false, "Shouldn't be here!")
    }
    assert(a(0,0) === "SSS")
    intercept[IndexOutOfBoundsException] {
      a(3,3)
      assert(false, "Shouldn't be here!")
    }
    a(0,0) = ":("
    assert(a(0,0) === ":(")
  }


  test("MapValues") {
    val a : HashMatrix[Int] = HashMatrix((1,0,0),(2,3,-1))

    val b1 : HashMatrix[Int] = a.mapValues(_ + 1)
    assert(b1 === HashMatrix((2,1,1),(3,4,0)))

    val b2 : HashMatrix[Double] = a.mapValues(_ + 1.0)
    assert(b2 === HashMatrix((2.0,1.0,1.0),(3.0,4.0,0.0)))
  }

  test("addition/subtraction") {
    val a : HashMatrix[Int] = HashMatrix((1,0,0),(2,3,-1))
    val b : HashMatrix[Int] = HashMatrix((0,1,0),(2,3,-1))
    assert(a + b === HashMatrix((1, 1, 0), (4,6,-2)))
    assert(a - b === HashMatrix((1, -1, 0), (0,0,0)))
  }

  test("addition/subtraction hash/dm") {
    val a : HashMatrix[Int] = HashMatrix((1,0,0),(2,3,-1))
    val b : DenseMatrix[Int] = DenseMatrix((0,1,0),(2,3,-1))
    assert(a + b === DenseMatrix((1, 1, 0), (4,6,-2)))
    assert(a - b === DenseMatrix((1, -1, 0), (0,0,0)))
    assert(b - a === -DenseMatrix((1, -1, 0), (0,0,0)))
  }

  test("inplace addition/subtraction") {
    val a : HashMatrix[Int] = HashMatrix((1,0,0),(2,3,-1))
    val b : HashMatrix[Int] = HashMatrix((0,1,0),(2,3,-1))
    a += b
    assert(a === HashMatrix((1, 1, 0), (4,6,-2)))
    //zero cells are never "returned" anyway, could be an improvement
    //TODO reclaim zeroed cells?

    // assert(a.activeSize === 5)
    a -= b
    a -= b
    assert(a === HashMatrix((1, -1, 0), (0,0,0)))
   // assert(a.activeSize === 2)
  }

  test("inplace set dm/hash") {
    val a : HashMatrix[Int] = HashMatrix((1,0,0),(2,3,-1))
    val b : DenseMatrix[Int] = DenseMatrix((0,1,0),(2,3,-1))
    b := a
    assert(a == b)
  }

  test("InPlace Ops") {
    var a = HashMatrix((1.0, 2.0, 3.0),(4.0, 5.0, 6.0))
    val b = HashMatrix((7.0, -2.0, 8.0),(-3.0, -3.0, 1.0))

    a :*= b
    assert(a === HashMatrix((7.0,-4.0,24.0),(-12.0,-15.0,6.0)))

    a = HashMatrix((1.0, 2.0, 3.0),(4.0, 5.0, 6.0))
    a :/= b
    assert(a === HashMatrix((1.0/7.0,-1.0,3.0/8.0),(4.0/(-3.0),5.0/(-3.0),6.0)))
  }

  test("hash scalar \"bad\" ops") {
    val a : HashMatrix[Int] = HashMatrix((1,0,0),(2,3,-1))
    assert(a :/ 3 === HashMatrix((0, 0, 0), (0, 1, 0)))

    val b : HashMatrix[Complex] = HashMatrix((Complex(1,0), Complex(0, 0), Complex(0, 0)), (Complex(2, 0), Complex(3, 0), Complex(-1, 0)))
    assert(b :/ Complex(3, 0) === HashMatrix((Complex(1.0 / 3.0, 0), Complex(0, 0), Complex(0, 0)), (Complex(2.0 / 3.0, 0), Complex(1, 0), Complex(-1.0 / 3.0, 0))))
  }

  test("hash scalar \"bad\" pow ops") {
    val a = HashMatrix((1.0, 2.0, 3.0),(4.0, 5.0, 6.0))
    val b = HashMatrix((7.0, -2.0, 8.0),(-3.0, -3.0, 1.0))

    assert(a :^ b === a.toDenseMatrix :^ b.toDenseMatrix)
    val ac = convert(a, Complex)
    val bc = convert(b, Complex)
    assert(ac :^ bc === ac.toDenseMatrix :^ bc.toDenseMatrix)
  }


  test("flatten") {
    val a = HashMatrix((1.0, 2.0, 3.0),(4.0, 5.0, 6.0))
    val b = HashMatrix.zeros[Double](3,2)
    b(0,1) = 1.0; b(2,1) = 3.0
    val z = HashMatrix.zeros[Double](5,3)
    assert(a.flatten() === SparseVector(1.0,2.0,3.0,4.0,5.0,6.0))
    assert(z.flatten() === SparseVector.zeros[Double](15))
    assert(b.flatten() === SparseVector(6)((1,1.0),(5,3.0)))
  }

  test("HashxHash: OpAddInPlace2:Field") {

    val hashA = HashMatrix.zeros[Double](3,4)
    val hashB = HashMatrix.zeros[Double](3,4)
    hashB(1,1) = 1.3
    hashB(0,0) = 1.0
    hashB(2,3) = 1.8
    hashB(2,0) = 1.6
    def testAddInPlace[T:Field:Zero:ClassTag](a: HashMatrix[T],b: HashMatrix[T]) = {
      a += b
    }
    testAddInPlace[Double](hashA,hashB)
    assert(hashA === hashB)
    testAddInPlace[Double](hashA,hashB)
    assert(hashA === hashB * 2.0)
    testAddInPlace[Double](hashA,HashMatrix.zeros[Double](3,4))
    assert(hashA === hashB * 2.0)
  }

  test("HashxHash: OpSubInPlace2:Field") {
    def testSubInPlace[T:Field:Zero:ClassTag](a: HashMatrix[T],b: HashMatrix[T]) = {
      a -= b
    }
    val hashA = HashMatrix.zeros[Double](3,4)
    val hashB = HashMatrix.zeros[Double](3,4)
    hashB(1,1) = 1.3
    hashB(0,0) = 1.0
    hashB(2,3) = 1.8
    hashB(2,0) = 1.6
    testSubInPlace[Double](hashA,hashB)
    assert(hashA === hashB * -1.0)
    testSubInPlace[Double](hashA,hashB)
    assert(hashA === hashB * -2.0)
    testSubInPlace[Double](hashA,HashMatrix.zeros[Double](3,4))
    assert(hashA === hashB * -2.0)
  }
  test("HashxHash: OpMulScalarInPlace2:Field") {
    def testMulScalarInPlace[T:Field:Zero:ClassTag](a: HashMatrix[T],b: HashMatrix[T]) = {
      a *= b
    }
    val hashA = HashMatrix.zeros[Double](3,4)
    val hashB = HashMatrix.zeros[Double](3,4)
    hashB(1,1) = 1.3
    hashB(0,0) = 1.0
    hashB(2,3) = 1.8
    hashB(2,0) = 1.6
    testMulScalarInPlace[Double](hashA,hashB)
    assert(hashA === hashA)
    hashA(1,1) = 2.0
    hashA(0,0) = 2.0
    hashA(1,0) = 2.0
    testMulScalarInPlace[Double](hashA,hashB)
    val hashR = HashMatrix.zeros[Double](3,4)
    hashR(1,1) = 2.6
    hashR(0,0) = 2.0
    assert(hashA === hashR)
    testMulScalarInPlace[Double](hashA,HashMatrix.zeros[Double](3,4))
    assert(hashA === HashMatrix.zeros[Double](3,4))
  }

  test("HashxHash: OpSetInPlace:Scalar:Field") {
    def testSetInPlace[T:Field:Zero:ClassTag](a: HashMatrix[T], b: T) = {
      a := b
    }
    val hashA = HashMatrix.zeros[Double](3,4)
    val b = 4
    testSetInPlace[Double](hashA,b)
    assert(hashA === HashMatrix.fill(3, 4)(b))
    testSetInPlace[Double](hashA,0)
    assert(hashA === HashMatrix.zeros[Double](3, 4))
  }
  test("ZipMapVals Test") {
    def testZipMap[T:Field:Zero:ClassTag](a: HashMatrix[T],b: HashMatrix[T]): HashMatrix[T] = {
      val f = implicitly[Field[T]]
      val addMapFn = (t1: T, t2: T) => f.+(t1,t2)
      val zmv = HashMatrix.zipMapVals[T,T]
      zmv.map(a,b,addMapFn)
    }
    val hashA = HashMatrix.zeros[Double](3,4)
    val hashB = HashMatrix.zeros[Double](3,4)
    hashB(1,1) = 1.3
    hashB(0,0) = 1.0
    hashB(2,3) = 1.8
    hashB(2,0) = 1.6
    val hashR = testZipMap(hashA,hashB)
    assert(hashR === hashB)
    val hashR1 = testZipMap(hashR,hashB)
    assert(hashR1 === hashB * 2.0)
    hashR(1,0) = 1.1
    hashB(0,1) = 1.2
    val hashR2 = testZipMap(hashR,hashB)
    val hashR3 = hashR * 2.0
    hashR3(1,0) = 1.1
    hashR3(0,1) = 1.2
    assert(hashR2 === hashR3)
    val hashR4 = testZipMap(hashB,hashA)
    assert(hashR4 === hashB)
  }



  test("#348") {
    val a = DenseMatrix((2,2),(3,3))
    val b = HashMatrix((2,2),(3,3))
    assert(a + b === a + b.toDenseMatrix)
  }

  //TODO implement solve for HashMatrix
/*
  test("HashMatrix Solve") {
    val r2 : DenseVector[Double] = HashMatrix((1.0,3.0,4.0),(2.0,0.0,6.0)) \ DenseVector(1.0,3.0)
    import breeze.numerics.inf
    assert( norm(r2 - DenseVector(0.1813186813186811, -0.3131868131868131, 0.43956043956043944), inf) < 1E-5)
  }
  */
}

