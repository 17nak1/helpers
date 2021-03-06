
var mathLib = {}
let exp = 2.718281828
let pi = 3.141592654
var erf = require('math-erf')

var libUnif = require('lib-r-math.js');
const {
    R: { numberPrecision },
    rng: { MersenneTwister, timeseed }
} = libUnif
var U = new MersenneTwister(0)
 // console.log(U.unif_rand());//0.8966972001362592

const libR = require('lib-r-math.js')
//* Set the seed for rbinom-In R: RNGkind("Knuth-TAOCP-2002");set.seed(1234)
const {
  Binomial,
  rng: { KnuthTAOCP2002 }
} = libR
const kn = new KnuthTAOCP2002(1234)
const { rnorm, rbinom } = Binomial(kn)
// console.log(rbinom(2,40,.5))

mathLib.pnorm = function (x, mu = 0, sd = 1, lower_tail = true, give_log = false) {
  if (sd < 0) {
    return NaN
  }
  let ans = 1 / 2 * (1 + erf((x - mu) / sd / Math.sqrt(2)))
  if (!lower_tail) {
    ans = 1 - ans
  }
  if (give_log) {
    ans = Math.log(ans)
  }
  return ans
}

// random generation for the normal distribution
mathLib.rnorm = function (mu = 0, sd = 1) {
  var val = Math.sqrt(-2 * Math.log(rng())) * Math.cos(2 * Math.PI * rng())
  return val * sd + mu
}


mathLib.numEulerSteps = function(t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  if (t1 >= t2) {
    dt = 0
    nstep = 0
  } else if (t1 + dt >= t2) {
    dt = t2 - t1
    nstep = 1
  } else {
    nstep = Math.ceil((t2 - t1) / dt /(1 + tol))
    dt = (t2 - t1) / nstep
  }
  return nstep
}

mathLib.numMapSteps = function (t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}

mathLib.nosortResamp = function (nw, w, np, p, offset) {
  for (j = 1; j < nw; j++) {
   w[j] += w[j-1]
 }
  if (w[nw - 1] <= 0) {
    throw "in 'systematic_resampling': non-positive sum of weight"
  }
  var du = w[nw - 1] / np
  var u = -du * U.unif_rand()//Math.random()

  for (j = 0, i = 0; j < np; j++) {
    u += du
    while ((u > w[i]) && (i < nw - 1)) i++;//looking for the low weight
    p[j] = i
  }
  if (offset){// add offset if needed
    for (j = 0; j < np; j++) p[j] += offset
  }
}


mathLib.reulermultinom = function (m = 1, size, rateAdd, dt, transAdd, rate, trans) {
  var p = 0
  var j, k
  if ((size < 0) || (dt < 0) || (Math.floor(size + 0.5) !== size)) {
    for (k = 0; k < m; k++) trans[k + transAdd] = NaN
    return 0
  }
  for (k = 0; k < m; k++) {
    if (rate[k + rateAdd] < 0.0) {
      for (j = 0; j < m; j++) trans[j + transAdd] = NaN
      return 0
    }
    p += rate[k + rateAdd]// total event rate
  }
  if (p > 0) {
    size = rbinom(1, size, 1 - Math.exp(-p * dt)) // total number of events
    if (!(isFinite(size)))
      throw 'result of binomial draw is not finite.'
    m -= 1
    for (k = 0; k < m; k++) {
      if (rate[k + rateAdd] > p) p = rate[k + rateAdd]
      trans[k + transAdd] = ((size > 0) && (p > 0)) ? rbinom(1, size, rate[k + rateAdd] / p) : 0
      if (!(isFinite(size) && isFinite(p) && isFinite(rate[k + rateAdd]) && isFinite(trans[k + transAdd]))) {
        throw 'result of binomial draw is not finite.'
      }
      size -= trans[k + transAdd]
      p -= rate[k + rateAdd]
    }
    trans[m + transAdd] = size
  } else {
    for (k = 0; k < m; k++) trans[k + transAdd] = 0
  }
}

mathLib.rpois = function (lambda = 1) {
  var k = 0; p = 1; l= Math.exp(-lambda)
  while (p > l) { 
    k += 1
    p = p * Math.random()
  }
  return k-1
}


mathLib.interpolator = function (points) {
  var first, n = points.length - 1,
    interpolated,
    leftExtrapolated,
    rightExtrapolated;

  if (points.length === 0) {
    return function () {
      return 0
    }
  }

  if (points.length === 1) {
    return function () {
      return points[0][1]
    }
  }

  points = points.sort(function (a, b) {
    return a[0] - b[0]
  })
  first = points[0]

  leftExtrapolated = function (x) {
    var a = points[0], b = points[1];
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  interpolated = function (x, a, b) {
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  rightExtrapolated = function (x) {
    var a = points[n - 1], b = points[n];
    return b[1] + (x - b[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  return function (x) {
    var i
    if (x <= first[0]) {
      return leftExtrapolated(x)
    }
    for (i = 0; i < n; i += 1) {
      if (x > points[i][0] && x <= points[i + 1][0]) {
        return interpolated(x, points[i], points[i + 1])
      }
    }
    return rightExtrapolated(x);
  }
}

/** Vector calculus
  * sum: adding; sp: scalar product; abs : Euclidean norm.
  */
sum = function (array) {
  var sum = []  
  for(i = 0; i < array[0].length; i++){
    var s= 0
    for (j = 0; j < array.length; j++) {
       s += array[j][i] 
    }
    sum.push(s)
  }
  return sum
}

sp = function (scalar, array) {
  var sum = []
  for(i = 0; i < array.length; i++){
   sum.push(scalar * array[i]);
  }
  return sum
}
abs = function (array) {
  var sum = 0
  for(i = 0; i < array.length; i++){
   sum += Math.pow(Math.abs(array[i]), 2)
  }
  return Math.sqrt(sum)
}

module.exports = mathLib;


