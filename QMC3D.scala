package QMC
import scala.util.Random
object Step{
	def apply() = 0.001
}
class Walker(a:Double, b:Double, c:Double){
	def V(x:Double,y:Double,z:Double):Double=(-1.0/Math.sqrt(x*x+y*y+z*z))
	def forward=(_:Double)+Math.sqrt(Step())*(Random.nextGaussian())
	var x:Double=a
	var y:Double=b
	var z:Double=c
	def getV():Double= V(x,y,z)
	def upDate()={
		x=forward(x)
		y=forward(y)
		z=forward(z)
	}
}
object Main extends App{
	val w_Ini = 1000  //initial walkers number
 	val step_Ini = 30000 //step number to go
	var Walkers=List.fill(w_Ini)(new Walker(Random.nextDouble()-0.5,Random.nextDouble()-0.5,Random.nextDouble()-0.5))
	var ER = (Walkers.map(_.getV()).foldLeft(0.0)(_+_))/Walkers.length
	var watcher=List[Double]()
	def Proliferator(a:Double)={
		if (Random.nextDouble<(a-a.toInt))
		a.toInt+1
		else
		a.toInt
	}
	0 to step_Ini foreach{
		_=>{
			var shadowWalkers = List[Walker]()
			println(ER)
			watcher:+=ER
			Walkers.foreach(_.upDate)
			Walkers.foreach{a =>
				val P_factor=Math.exp(Step()*(ER-a.getV()))
				shadowWalkers:::=List.fill(Proliferator(P_factor))(a)	
			}	
			Walkers = shadowWalkers
			ER+=Math.log(w_Ini.toDouble/Walkers.length)/10
		}
	}
	val drop_rate=0.2 // drop unstable sample data from head
	val sample=watcher.drop(step_Ini*drop_rate.toInt)
	val mean=sample.sum/sample.length
	def variance(data:List[Double],avg:Double)=data.map(a=>Math.pow(a-avg, 2)).sum/data.length
	println("avg="+mean+"+/-"+Math.sqrt(variance(sample, mean)/sample.length))
	println("variance="+variance(sample, mean))
	val L=40 //block size
	println("resampling variance="+variance(sample.grouped(L).toList.map(a=>a.sum/a.length),mean))
}
